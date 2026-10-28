mod alignment;
mod depth;
mod jobs;
mod output;
mod samples;
mod scan;

use self::{alignment::*, depth::*, jobs::*, output::*, samples::*, scan::*};

use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

use anyhow::{Context, Result, anyhow, bail};
use noodles_bam as bam;
use noodles_core::{Position, Region};
use noodles_sam::{
    self as sam,
    alignment::{
        Record as SamRecord,
        record::data::field::{Tag, Value},
    },
};
use rayon::{ThreadPoolBuilder, prelude::*};
use serde::Serialize;
use tracing::{debug, info, warn};

use crate::{
    annotation::{GeneSpan, breakpoint_annotation, load_target_spans},
    cli::Args,
    fasta::FastaWriter,
    vcf::{Junction, StructuralVariant, cluster_consensus, write_vcf},
};

#[derive(Debug, Clone, Serialize)]
struct TsvRecord {
    query_gene: String,
    matched_partner_gene: Option<String>,
    query_transcript_id: String,
    partner_transcript_id: String,
    breakpoint_estimate: String,
    read_name: String,
    read_flags: u16,
    reference_name: String,
    alignment_start: usize,
    alignment_end: usize,
    cigar: String,
    mapping_quality: Option<u8>,
    mate_reference_name: Option<String>,
    mate_alignment_start: Option<usize>,
    inferred_partner_reference: String,
    inferred_partner_start: usize,
    inferred_partner_strand: String,
    sa_tag: String,
    sample: String,
}

pub fn run(args: Args) -> Result<()> {
    validate_inputs(&args)?;

    if args.threads > 0 {
        ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .map_err(|error| anyhow!("failed to configure rayon thread pool: {error}"))?;
    }

    if args.min_mapq == 0 {
        warn!(
            "--min-mapq is 0: taking every alignment regardless of mapping quality. \
             Output may include low-quality and multi-mapping reads, which are a common \
             source of spurious fusion candidates; raise --min-mapq to filter them."
        );
    }

    // Depth only exists for consensus calls, so the filter is a silent no-op
    // without a VCF to filter.
    if args.min_depth > 0 && args.output_vcf.is_none() {
        warn!("--min-depth only filters the consensus VCF and has no effect without --output-vcf");
    }

    // --partner-gene is a single global constraint; in batch mode each row
    // carries its own, so the flag has no meaning there.
    if args.loci.is_some() && args.partner_gene.is_some() {
        warn!("--partner-gene is ignored in --loci mode; each row carries its own partner");
    }

    let samples = open_bam_samples(&args.bam)?;
    info!("processing {} BAM sample(s)", samples.len());

    let jobs = build_jobs(&args)?;

    // Load only the referenced genes when every job names a partner; a job
    // without one annotates against any overlapping gene and so needs them all.
    let needs_all_spans = jobs.iter().any(|job| job.partner_gene.is_none());
    let annotation_genes = if needs_all_spans {
        Vec::new()
    } else {
        referenced_genes(&jobs)
    };
    let all_spans = load_target_spans(&args.annotation, &annotation_genes)?;

    let genes_token = output_genes_token(&args);
    let output_tsv = resolve_output_path(args.output_tsv.as_deref(), &samples, &genes_token, "tsv");
    let output_fasta = resolve_output_path(
        args.output_fasta.as_deref(),
        &samples,
        &genes_token,
        "fasta.gz",
    );
    let output_vcf = args
        .output_vcf
        .as_ref()
        .map(|explicit| resolve_output_path(explicit.as_deref(), &samples, &genes_token, "vcf"));

    let plans: Vec<JobPlan> = jobs
        .iter()
        .map(|job| build_job_plan(job, &all_spans, output_vcf.is_some()))
        .collect();
    let query_interval_count: usize = plans.iter().map(|plan| plan.query_spans.len()).sum();
    info!(
        "loaded {} query interval(s) across {} job(s)",
        query_interval_count,
        plans.len()
    );

    let tsv_writer = Arc::new(Mutex::new(create_tsv_writer(&output_tsv)?));
    let fasta_writer = Arc::new(Mutex::new(FastaWriter::create(&output_fasta)?));
    write_tsv_header(&tsv_writer)?;

    let scan_options = ScanOptions::from_args(&args);

    // Scan every (job, sample, query interval) in parallel. Scoped so the
    // borrow of `plans` ends before its junctions are consumed below.
    {
        let ctx = ScanContext {
            all_spans: &all_spans,
            tsv_writer: &tsv_writer,
            fasta_writer: &fasta_writer,
            options: scan_options,
        };
        let work: Vec<(&JobPlan, &BamSample, &GeneSpan)> = plans
            .iter()
            .flat_map(|plan| {
                samples.iter().flat_map(move |sample| {
                    plan.query_spans
                        .iter()
                        .map(move |span| (plan, sample, span))
                })
            })
            .collect();
        work.par_iter()
            .try_for_each(|&(plan, sample, span)| process_span(sample, span, plan, &ctx))?;
    }

    let fasta_writer = Arc::into_inner(fasta_writer)
        .ok_or_else(|| anyhow!("failed to reclaim FASTA writer"))?
        .into_inner()
        .map_err(|_| anyhow!("FASTA writer lock was poisoned"))?;
    fasta_writer.finish()?;

    let mut tsv_writer = Arc::into_inner(tsv_writer)
        .ok_or_else(|| anyhow!("failed to reclaim TSV writer"))?
        .into_inner()
        .map_err(|_| anyhow!("TSV writer lock was poisoned"))?;
    tsv_writer.flush()?;

    info!(
        "finished writing {} and {}",
        output_tsv.display(),
        output_fasta.display()
    );

    if let Some(vcf_path) = output_vcf.as_ref() {
        // Cluster each job's junctions with that job's own tolerance, then
        // combine the calls and order them for a tidy VCF.
        let mut variants: Vec<StructuralVariant> = Vec::new();
        for plan in plans {
            if let Some(collector) = plan.junctions {
                let junctions = collector
                    .into_inner()
                    .map_err(|_| anyhow!("junction collector lock was poisoned"))?;
                variants.extend(cluster_consensus(junctions, plan.slop));
            }
        }
        variants.sort_by(|a, b| {
            a.chrom1
                .cmp(&b.chrom1)
                .then(a.pos1.cmp(&b.pos1))
                .then(a.chrom2.cmp(&b.chrom2))
                .then(a.pos2.cmp(&b.pos2))
        });

        // Depth is only knowable once the consensus breakpoint is fixed, so it
        // is measured here rather than during the streaming scan. Each sample
        // opens its reader once and walks every breakpoint, and samples run in
        // parallel.
        let per_sample: Vec<(String, Vec<usize>)> = samples
            .par_iter()
            .map(|sample| {
                let depths = sample_depths(sample, &variants, scan_options)?;
                Ok((sample.name.clone(), depths))
            })
            .collect::<Result<Vec<_>>>()?;

        for (sample_name, depths) in per_sample {
            for (variant, depth) in variants.iter_mut().zip(depths) {
                variant.depth_by_sample.insert(sample_name.clone(), depth);
            }
        }

        // Drop calls whose breakpoint is too thinly covered to interpret: an
        // allele fraction over a denominator of one or two reads says little.
        // A call survives if any single sample is deep enough, so one shallow
        // sample cannot sink a cohort-wide call.
        if args.min_depth > 0 {
            let before = variants.len();
            variants.retain(|variant| {
                samples
                    .iter()
                    .any(|sample| variant.depth(&sample.name) >= args.min_depth)
            });

            let dropped = before - variants.len();
            if dropped > 0 {
                info!(
                    "dropped {dropped} consensus call(s) below --min-depth {}",
                    args.min_depth
                );
            }
        }

        let sample_names = sample_names(&samples);
        let contigs = union_contigs(&samples);
        write_vcf(vcf_path, &variants, &sample_names, &contigs)?;
        info!(
            "wrote {} consensus structural variant(s) to {}",
            variants.len(),
            vcf_path.display()
        );
    }

    Ok(())
}

pub fn validate_output_paths(args: &Args) -> Result<()> {
    output::validate_output_paths(args)
}

#[cfg(test)]
#[path = "pipeline/tests.rs"]
mod tests;
