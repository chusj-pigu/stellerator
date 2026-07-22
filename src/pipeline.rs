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

/// Count reads spanning a breakpoint in one sample, applying the same record
/// filters as the main scan.
///
/// This is the denominator behind allele fraction: without it a support count
/// cannot distinguish a clonal event from a handful of artefacts. Reads are
/// counted by name so a read with several alignments over the position counts
/// once, and supplementary alignments are skipped for the same reason.
fn sample_depths(
    sample: &BamSample,
    variants: &[StructuralVariant],
    options: ScanOptions,
) -> Result<Vec<usize>> {
    if variants.is_empty() {
        return Ok(Vec::new());
    }

    // Opening an indexed reader parses the whole BAM index, so it is done once
    // per sample and reused for every breakpoint rather than per call.
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&sample.path)
        .with_context(|| format!("failed to open indexed BAM {}", sample.path.display()))?;

    let known_references: std::collections::HashSet<String> = sample
        .header
        .reference_sequences()
        .keys()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .collect();

    let mut depths = Vec::with_capacity(variants.len());

    for variant in variants {
        let mut names: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();

        // A contig absent from this sample's header contributes no depth.
        if known_references.contains(&variant.chrom1) {
            let start = Position::try_from(variant.pos1)
                .map_err(|_| anyhow!("invalid breakpoint position {}", variant.pos1))?;
            let region = Region::new(variant.chrom1.clone(), start..=start);

            let query = reader.query(&sample.header, &region)?;
            for result in query.records() {
                let record = result?;
                let flags = record.flags();
                if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                    continue;
                }
                if flags.is_duplicate() && !options.include_duplicates {
                    continue;
                }
                if options.min_mapq > 0
                    && !matches!(record.mapping_quality(), Some(mapq) if u8::from(mapq) >= options.min_mapq)
                {
                    continue;
                }
                if let Some(name) = record.name() {
                    names.insert(String::from_utf8_lossy(name.as_ref()).into_owned());
                }
            }
        }

        // A supporting read spans the junction by definition, even when
        // clipping leaves its alignment ending short of the consensus
        // breakpoint. Without this, AD would not sum to DP.
        if let Some(support) = variant.support_reads(&sample.name) {
            names.extend(support.iter().cloned());
        }

        depths.push(names.len());
    }

    Ok(depths)
}

/// Deterministically ordered sample names for VCF genotype columns.
fn sample_names(samples: &[BamSample]) -> Vec<String> {
    let mut names: Vec<String> = samples.iter().map(|sample| sample.name.clone()).collect();
    names.sort();
    names
}

/// Union of reference sequences across all sample headers, in first-seen order,
/// used to emit `##contig` lines in the VCF.
fn union_contigs(samples: &[BamSample]) -> Vec<(String, Option<usize>)> {
    let mut seen = std::collections::HashSet::new();
    let mut contigs = Vec::new();

    for sample in samples {
        for (name, reference) in sample.header.reference_sequences() {
            let name = String::from_utf8_lossy(name.as_ref()).into_owned();
            if seen.insert(name.clone()) {
                contigs.push((name, Some(usize::from(reference.length()))));
            }
        }
    }

    contigs
}

/// Resolve an output path: use `explicit` when provided, otherwise build a
/// default from the BAM basename and requested genes (e.g. `sample.BCR_ABL1.tsv`).
fn resolve_output_path(
    explicit: Option<&Path>,
    samples: &[BamSample],
    genes_token: &str,
    extension: &str,
) -> PathBuf {
    if let Some(path) = explicit {
        return path.to_path_buf();
    }

    let paths: Vec<PathBuf> = samples.iter().map(|sample| sample.path.clone()).collect();
    let names: Vec<String> = samples.iter().map(|sample| sample.name.clone()).collect();
    let bam_token = bam_basename_token(&paths, &names);
    default_output_path(&bam_token, genes_token, extension)
}

/// The gene component of default output names: the requested genes in single
/// mode, or the loci file stem in batch mode.
fn output_genes_token(args: &Args) -> String {
    match &args.loci {
        Some(path) => path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .map(sanitize_token)
            .filter(|token| !token.is_empty())
            .unwrap_or_else(|| "loci".to_string()),
        None => genes_basename_token(&args.genes, args.partner_gene.as_deref()),
    }
}

fn default_output_path(bam_token: &str, genes_token: &str, extension: &str) -> PathBuf {
    PathBuf::from(format!("{bam_token}.{genes_token}.{extension}"))
}

/// Representative BAM basename for default output names: the sample stem for a
/// single BAM, the shared parent directory name for several BAMs in one place,
/// otherwise the first sample's name.
fn bam_basename_token(paths: &[PathBuf], names: &[String]) -> String {
    if names.len() == 1 {
        return sanitize_token(&names[0]);
    }

    if let Some(parent) = paths.first().and_then(|path| path.parent())
        && !parent.as_os_str().is_empty()
        && paths.iter().all(|path| path.parent() == Some(parent))
        && let Some(name) = parent.file_name().and_then(|name| name.to_str())
        && !name.is_empty()
    {
        return sanitize_token(name);
    }

    names
        .first()
        .map(|name| sanitize_token(name))
        .unwrap_or_else(|| "stellerator".to_string())
}

/// Join the requested genes (and partner gene, if distinct) into a filename
/// token, e.g. `BCR_ABL1`.
fn genes_basename_token(genes: &[String], partner_gene: Option<&str>) -> String {
    let mut parts: Vec<String> = Vec::new();
    for gene in genes {
        push_unique_token(&mut parts, gene);
    }
    if let Some(partner) = partner_gene {
        push_unique_token(&mut parts, partner);
    }

    if parts.is_empty() {
        "genes".to_string()
    } else {
        parts.join("_")
    }
}

fn push_unique_token(parts: &mut Vec<String>, value: &str) {
    let token = sanitize_token(value);
    if !token.is_empty()
        && !parts
            .iter()
            .any(|existing| existing.eq_ignore_ascii_case(&token))
    {
        parts.push(token);
    }
}

/// Replace filename-hostile characters so derived output paths stay valid.
fn sanitize_token(value: &str) -> String {
    value
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn validate_inputs(args: &Args) -> Result<()> {
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }

    if !args.annotation.exists() {
        bail!(
            "annotation input does not exist: {}",
            args.annotation.display()
        );
    }

    Ok(())
}

/// A single indexed BAM to scan, carrying the sample name used for output
/// provenance and its parsed header.
struct BamSample {
    path: PathBuf,
    name: String,
    header: sam::Header,
}

impl BamSample {
    fn open(path: &Path) -> Result<Self> {
        let header = load_sam_header(path)?;
        Ok(Self {
            path: path.to_path_buf(),
            name: sample_name(path),
            header,
        })
    }
}

/// Resolve the requested BAM inputs (files and/or directories) into a validated,
/// deduplicated set of indexed BAMs and open each as a [`BamSample`].
fn open_bam_samples(inputs: &[PathBuf]) -> Result<Vec<BamSample>> {
    let paths = resolve_bam_inputs(inputs)?;
    check_unique_sample_names(&paths)?;
    paths.iter().map(|path| BamSample::open(path)).collect()
}

/// Expand directory inputs to their `.bam` contents, confirm every BAM exists
/// and is indexed, then sort and deduplicate the result.
fn resolve_bam_inputs(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut resolved: Vec<PathBuf> = Vec::new();

    for input in inputs {
        if input.is_dir() {
            let mut found = directory_bam_files(input)?;
            if found.is_empty() {
                bail!("no .bam files found in directory {}", input.display());
            }
            resolved.append(&mut found);
        } else if input.exists() {
            resolved.push(input.clone());
        } else {
            bail!("BAM input does not exist: {}", input.display());
        }
    }

    resolved.sort();
    resolved.dedup();

    if resolved.is_empty() {
        bail!("no BAM inputs provided");
    }

    for path in &resolved {
        if !has_associated_index(path) {
            bail!(
                "indexed BAM required; expected {path}.bai or {path}.csi",
                path = path.display()
            );
        }
    }

    Ok(resolved)
}

/// List the `.bam` files directly inside `dir`, sorted by path.
fn directory_bam_files(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in
        fs::read_dir(dir).with_context(|| format!("failed to read directory {}", dir.display()))?
    {
        let path = entry?.path();
        if path.is_file() && has_bam_extension(&path) {
            files.push(path);
        }
    }

    files.sort();
    Ok(files)
}

fn has_bam_extension(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
}

/// Derive the sample name for a BAM from its file stem, falling back to the
/// full path when no stem is available.
fn sample_name(path: &Path) -> String {
    path.file_stem()
        .and_then(|stem| stem.to_str())
        .map(|stem| stem.to_string())
        .unwrap_or_else(|| path.to_string_lossy().into_owned())
}

/// Reject input sets where two BAMs would collapse to the same sample name,
/// which would make per-sample output provenance ambiguous.
fn check_unique_sample_names(paths: &[PathBuf]) -> Result<()> {
    let mut seen: HashMap<String, &Path> = HashMap::new();

    for path in paths {
        let name = sample_name(path);
        if let Some(existing) = seen.insert(name.clone(), path) {
            bail!(
                "duplicate sample name {name:?} derived from {} and {}; rename one of the BAM files to disambiguate",
                existing.display(),
                path.display()
            );
        }
    }

    Ok(())
}

fn load_sam_header(path: &Path) -> Result<sam::Header> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .with_context(|| format!("failed to open BAM {}", path.display()))?;
    reader.read_header().context("failed to read BAM header")
}

fn has_associated_index(bam_path: &Path) -> bool {
    let bam_bai = bam_path.with_extension(format!(
        "{}.bai",
        bam_path
            .extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or_default()
    ));
    let bam_csi = bam_path.with_extension(format!(
        "{}.csi",
        bam_path
            .extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or_default()
    ));

    bam_bai.exists() || bam_csi.exists()
}

/// Shared, read-only state for scanning gene spans: where results go and how
/// records are filtered. Held once per run and borrowed by every worker.
/// One unit of work: query genes, an optional partner constraint, and the
/// clustering tolerance for the calls it produces. Single-CLI mode is one job;
/// `--loci` mode is one per row.
struct Job {
    query_genes: Vec<String>,
    partner_gene: Option<String>,
    slop: usize,
}

/// A job with its annotation spans resolved and its own junction collector,
/// ready to scan.
struct JobPlan {
    query_spans: Vec<GeneSpan>,
    /// Spans the partner side is restricted to; `None` annotates against every
    /// span (the `--partner-gene`-omitted behaviour) and requires no match.
    partner_filter: Option<Vec<GeneSpan>>,
    slop: usize,
    junctions: Option<Mutex<Vec<Junction>>>,
}

/// Scan-wide state shared by every worker, independent of the job.
struct ScanContext<'a> {
    all_spans: &'a [GeneSpan],
    tsv_writer: &'a Arc<Mutex<BufWriter<File>>>,
    fasta_writer: &'a Arc<Mutex<FastaWriter>>,
    options: ScanOptions,
}

fn build_jobs(args: &Args) -> Result<Vec<Job>> {
    match &args.loci {
        Some(path) => {
            let requests = crate::loci::parse_loci_file(path)?;
            info!(
                "loaded {} locus job(s) from {}",
                requests.len(),
                path.display()
            );
            Ok(requests
                .into_iter()
                .map(|request| Job {
                    query_genes: vec![request.gene],
                    partner_gene: request.partner,
                    slop: request.tolerance.unwrap_or(args.sv_slop),
                })
                .collect())
        }
        None => Ok(vec![Job {
            query_genes: args.genes.clone(),
            partner_gene: args.partner_gene.clone(),
            slop: args.sv_slop,
        }]),
    }
}

/// Every gene named across the jobs, query and partner sides both.
fn referenced_genes(jobs: &[Job]) -> Vec<String> {
    let mut names = Vec::new();
    for job in jobs {
        names.extend(job.query_genes.iter().cloned());
        if let Some(partner) = &job.partner_gene {
            names.push(partner.clone());
        }
    }
    names
}

fn build_job_plan(job: &Job, all_spans: &[GeneSpan], collect_junctions: bool) -> JobPlan {
    let query_spans = filter_spans_by_gene_names(all_spans, &job.query_genes);
    for gene in &job.query_genes {
        if !query_spans
            .iter()
            .any(|span| span.gene.eq_ignore_ascii_case(gene))
        {
            warn!("gene {gene:?} has no intervals in the annotation; it will produce no output");
        }
    }

    let partner_filter = resolve_partner_filter(all_spans, job.partner_gene.as_deref());
    let junctions = collect_junctions.then(|| Mutex::new(Vec::<Junction>::new()));

    JobPlan {
        query_spans,
        partner_filter,
        slop: job.slop,
        junctions,
    }
}

/// Restrict partner annotation to a named gene, or `None` to annotate against
/// every span.
fn resolve_partner_filter(
    all_spans: &[GeneSpan],
    partner_gene: Option<&str>,
) -> Option<Vec<GeneSpan>> {
    partner_gene.map(|partner| {
        filter_spans_by_gene_names(all_spans, std::slice::from_ref(&partner.to_string()))
    })
}

fn process_span(
    sample: &BamSample,
    span: &GeneSpan,
    plan: &JobPlan,
    ctx: &ScanContext<'_>,
) -> Result<()> {
    debug!(
        sample = sample.name,
        gene = span.gene,
        reference = span.reference_name,
        start = span.start,
        end = span.end,
        "querying gene interval"
    );

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(&sample.path)
        .with_context(|| format!("failed to open indexed BAM {}", sample.path.display()))?;

    // A job without a partner filter annotates against every span.
    let partner_spans = plan.partner_filter.as_deref().unwrap_or(ctx.all_spans);
    let require_partner_match = plan.partner_filter.is_some();

    let region = build_region(span)?;
    let query = reader.query(&sample.header, &region)?;
    for result in query.records() {
        let record = result?;
        let RecordHits { sequence, hits } = classify_record(
            &sample.header,
            &sample.name,
            span,
            Some(partner_spans),
            require_partner_match,
            &record,
            ctx.options,
        );

        for (index, hit) in hits.into_iter().enumerate() {
            write_tsv_row(ctx.tsv_writer, &hit)?;

            // One FASTA record per read: a chimeric read supporting several
            // junctions gets one row per junction in the TSV, but its sequence
            // is written once. The header describes the first junction.
            if index == 0 {
                write_fasta_record(ctx.fasta_writer, &hit.fasta_header, &sequence)?;
            }

            if let Some(collector) = &plan.junctions {
                collector
                    .lock()
                    .map_err(|_| anyhow!("junction collector lock was poisoned"))?
                    .push(hit.junction);
            }
        }
    }

    Ok(())
}

fn build_region(span: &GeneSpan) -> Result<Region> {
    let start = Position::try_from(span.start as usize)
        .map_err(|_| anyhow!("invalid start coordinate for {}", span.gene))?;
    let end = Position::try_from(span.end as usize)
        .map_err(|_| anyhow!("invalid end coordinate for {}", span.gene))?;
    Ok(Region::new(span.reference_name.clone(), start..=end))
}

/// Per-read fields shared across every supplementary alignment of one record.
struct ReadContext {
    read_name: String,
    sequence: String,
    read_flags: u16,
    query_strand: char,
    reference_name: String,
    alignment_start: usize,
    alignment_end: usize,
    query_breakpoint_position: usize,
    cigar: String,
    mapping_quality: Option<u8>,
    mate_reference_name: Option<String>,
    mate_alignment_start: Option<usize>,
    sa_tag: String,
}

/// Read-filtering options applied to every scanned record.
#[derive(Debug, Clone, Copy)]
struct ScanOptions {
    include_duplicates: bool,
    min_mapq: u8,
}

impl ScanOptions {
    fn from_args(args: &Args) -> Self {
        Self {
            include_duplicates: args.include_duplicates,
            min_mapq: args.min_mapq,
        }
    }
}

fn classify_record(
    header: &sam::Header,
    sample: &str,
    span: &GeneSpan,
    partner_spans: Option<&[GeneSpan]>,
    require_partner_match: bool,
    record: &bam::Record,
    options: ScanOptions,
) -> RecordHits {
    let Some(context) = read_context(header, record, options) else {
        return RecordHits {
            sequence: String::new(),
            hits: Vec::new(),
        };
    };

    let hits = parse_sa_entries(&context.sa_tag)
        .into_iter()
        .filter_map(|partner| {
            build_hit(
                sample,
                span,
                partner_spans,
                require_partner_match,
                &context,
                partner,
            )
        })
        .collect();

    RecordHits {
        sequence: context.sequence,
        hits,
    }
}

/// Compute the per-read context, or `None` if the record cannot support a fusion
/// call (unmapped, secondary, duplicate, missing `SA` tag, or no sequence).
fn read_context(
    header: &sam::Header,
    record: &bam::Record,
    options: ScanOptions,
) -> Option<ReadContext> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() {
        return None;
    }

    // PCR/optical duplicates inflate apparent support, so drop them unless the
    // caller opts in.
    if flags.is_duplicate() && !options.include_duplicates {
        return None;
    }

    // A floor of 0 takes everything. Above that, drop alignments below the
    // threshold, including records with no reported MAPQ since their quality
    // cannot be verified.
    if options.min_mapq > 0
        && !matches!(record.mapping_quality(), Some(mapq) if u8::from(mapq) >= options.min_mapq)
    {
        return None;
    }

    let sa_tag = extract_sa(record)?;
    let read_name = String::from_utf8_lossy(record.name()?.as_ref()).into_owned();
    let sequence: String = record.sequence().iter().map(char::from).collect();
    if sequence.is_empty() {
        return None;
    }

    let reference_name = reference_name_for_id(header, record.reference_sequence_id()?.ok()?)?;
    let alignment_start = usize::from(record.alignment_start()?.ok()?);
    let alignment_end = usize::from(SamRecord::alignment_end(record)?.ok()?);
    let query_breakpoint_position =
        estimate_query_breakpoint_position(record, alignment_start, alignment_end)?;
    let cigar = cigar_to_string(&record.cigar()).ok()?;
    let mapping_quality = record.mapping_quality().map(u8::from);
    let mate_reference_name = record
        .mate_reference_sequence_id()
        .and_then(|id| id.ok())
        .and_then(|id| reference_name_for_id(header, id));
    let mate_alignment_start = record
        .mate_alignment_start()
        .and_then(|position| position.ok())
        .map(usize::from);
    let query_strand = if flags.is_reverse_complemented() {
        '-'
    } else {
        '+'
    };

    Some(ReadContext {
        read_name,
        sequence,
        read_flags: flags.bits(),
        query_strand,
        reference_name,
        alignment_start,
        alignment_end,
        query_breakpoint_position,
        cigar,
        mapping_quality,
        mate_reference_name,
        mate_alignment_start,
        sa_tag,
    })
}

/// Build a hit for one supplementary alignment, applying the same-gene and
/// partner-match filters. Returns `None` when the partner does not qualify.
fn build_hit(
    sample: &str,
    span: &GeneSpan,
    partner_spans: Option<&[GeneSpan]>,
    require_partner_match: bool,
    context: &ReadContext,
    partner: PartnerAlignment,
) -> Option<Hit> {
    if partner.reference_name == span.reference_name
        && partner.start >= span.start as usize
        && partner.start <= span.end as usize
    {
        return None;
    }

    let matched_partner_span = partner_spans.and_then(|spans| {
        find_overlapping_span(spans, &partner.reference_name, partner.breakpoint)
            .or_else(|| find_overlapping_span(spans, &partner.reference_name, partner.start))
            .filter(|partner_span| partner_span.gene != span.gene)
    });

    if require_partner_match && matched_partner_span.is_none() {
        return None;
    }

    let query_breakpoint = breakpoint_annotation(span, context.query_breakpoint_position);
    let partner_breakpoint = matched_partner_span
        .and_then(|partner_span| breakpoint_annotation(partner_span, partner.breakpoint));
    let query_transcript_id = query_breakpoint
        .as_ref()
        .map(|annotation| annotation.transcript_id.clone())
        .unwrap_or_else(|| "unknown".to_string());
    let partner_transcript_id = partner_breakpoint
        .as_ref()
        .map(|annotation| annotation.transcript_id.clone())
        .unwrap_or_else(|| "unknown".to_string());
    let query_breakpoint_region = query_breakpoint
        .as_ref()
        .map(|annotation| annotation.region.clone())
        .unwrap_or_else(|| "unknown".to_string());
    let partner_breakpoint_region = partner_breakpoint
        .as_ref()
        .map(|annotation| annotation.region.clone())
        .unwrap_or_else(|| "unknown".to_string());
    let matched_partner_gene = matched_partner_span.map(|partner_span| partner_span.gene.clone());
    let breakpoint_estimate = format!("{query_breakpoint_region}/{partner_breakpoint_region}");

    let junction = Junction {
        sample: sample.to_string(),
        read_name: context.read_name.clone(),
        query_gene: span.gene.clone(),
        partner_gene: matched_partner_gene.clone(),
        query_transcript: query_transcript_id.clone(),
        partner_transcript: partner_transcript_id.clone(),
        query_region: query_breakpoint_region,
        partner_region: partner_breakpoint_region,
        chrom1: context.reference_name.clone(),
        pos1: context.query_breakpoint_position,
        strand1: context.query_strand,
        chrom2: partner.reference_name.clone(),
        pos2: partner.breakpoint,
        strand2: partner.strand,
    };

    Some(Hit {
        tsv: TsvRecord {
            query_gene: span.gene.clone(),
            matched_partner_gene: matched_partner_gene.clone(),
            query_transcript_id: query_transcript_id.clone(),
            partner_transcript_id: partner_transcript_id.clone(),
            breakpoint_estimate: breakpoint_estimate.clone(),
            read_name: context.read_name.clone(),
            read_flags: context.read_flags,
            reference_name: context.reference_name.clone(),
            alignment_start: context.alignment_start,
            alignment_end: context.alignment_end,
            cigar: context.cigar.clone(),
            mapping_quality: context.mapping_quality,
            mate_reference_name: context.mate_reference_name.clone(),
            mate_alignment_start: context.mate_alignment_start,
            inferred_partner_reference: partner.reference_name.clone(),
            inferred_partner_start: partner.start,
            inferred_partner_strand: partner.strand.to_string(),
            sa_tag: context.sa_tag.clone(),
            sample: sample.to_string(),
        },
        fasta_header: format!(
            "{} gene={} matched_partner_gene={} query_transcript_id={} partner_transcript_id={} breakpoint_estimate={} partner={}:{} strand={} sample={}",
            context.read_name,
            span.gene,
            matched_partner_gene.unwrap_or_else(|| "NA".to_string()),
            query_transcript_id,
            partner_transcript_id,
            breakpoint_estimate,
            partner.reference_name,
            partner.start,
            partner.strand,
            sample
        ),
        junction,
    })
}

fn extract_sa(record: &bam::Record) -> Option<String> {
    let data = record.data();
    let value = data.get(&Tag::OTHER_ALIGNMENTS)?.ok()?;

    match value {
        Value::String(raw) => Some(String::from_utf8_lossy(raw.as_ref()).into_owned()),
        _ => None,
    }
}

fn reference_name_for_id(header: &sam::Header, id: usize) -> Option<String> {
    header
        .reference_sequences()
        .get_index(id)
        .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned())
}

fn filter_spans_by_gene_names(all_spans: &[GeneSpan], names: &[String]) -> Vec<GeneSpan> {
    let wanted: Vec<String> = names.iter().map(|name| name.to_ascii_lowercase()).collect();

    all_spans
        .iter()
        .filter(|span| wanted.contains(&span.gene.to_ascii_lowercase()))
        .cloned()
        .collect()
}

fn find_overlapping_span<'a>(
    spans: &'a [GeneSpan],
    reference_name: &str,
    position: usize,
) -> Option<&'a GeneSpan> {
    spans.iter().find(|span| {
        span.reference_name == reference_name
            && position >= span.start as usize
            && position <= span.end as usize
    })
}

fn estimate_query_breakpoint_position(
    record: &bam::Record,
    alignment_start: usize,
    alignment_end: usize,
) -> Option<usize> {
    use noodles_sam::alignment::record::cigar::op::Kind;

    let ops: Vec<_> = record.cigar().iter().collect::<Result<_, _>>().ok()?;
    let first = ops.first()?.kind();
    let last = ops.last()?.kind();
    let left_clipped = matches!(first, Kind::SoftClip | Kind::HardClip);
    let right_clipped = matches!(last, Kind::SoftClip | Kind::HardClip);

    Some(match (left_clipped, right_clipped) {
        (true, false) => alignment_start,
        (false, true) => alignment_end,
        _ => alignment_end,
    })
}

fn cigar_to_string(cigar: &bam::record::Cigar<'_>) -> Result<String, std::io::Error> {
    use noodles_sam::alignment::record::cigar::op::Kind;

    let mut rendered = String::new();

    for op in cigar.iter() {
        let op = op?;
        rendered.push_str(&op.len().to_string());
        rendered.push(match op.kind() {
            Kind::Match => 'M',
            Kind::Insertion => 'I',
            Kind::Deletion => 'D',
            Kind::Skip => 'N',
            Kind::SoftClip => 'S',
            Kind::HardClip => 'H',
            Kind::Pad => 'P',
            Kind::SequenceMatch => '=',
            Kind::SequenceMismatch => 'X',
        });
    }

    Ok(rendered)
}

/// Parse every entry in an `SA` tag into a partner alignment.
fn parse_sa_entries(raw: &str) -> Vec<PartnerAlignment> {
    raw.split(';')
        .filter(|entry| !entry.trim().is_empty())
        .filter_map(parse_sa_entry)
        .collect()
}

fn parse_sa_entry(entry: &str) -> Option<PartnerAlignment> {
    let mut fields = entry.split(',');
    let reference_name = fields.next()?.trim().to_string();
    let start: usize = fields.next()?.trim().parse().ok()?;
    let strand = fields.next()?.trim().chars().next()?;
    let cigar = fields.next().unwrap_or("").trim();
    let breakpoint = estimate_sa_breakpoint_position(start, cigar).unwrap_or(start);
    Some(PartnerAlignment {
        reference_name,
        start,
        strand,
        breakpoint,
    })
}

/// Estimate the partner-side breakpoint from an `SA` CIGAR: the aligned start
/// when the alignment is left-clipped, otherwise its aligned end. Mirrors the
/// query-side estimate so partner exon/intron labels use a comparable position.
fn estimate_sa_breakpoint_position(start: usize, cigar: &str) -> Option<usize> {
    let ops = parse_cigar_ops(cigar);
    if ops.is_empty() {
        return None;
    }

    let reference_span: usize = ops
        .iter()
        .filter(|(_, kind)| matches!(kind, 'M' | 'D' | 'N' | '=' | 'X'))
        .map(|(len, _)| *len)
        .sum();
    let end = start + reference_span.saturating_sub(1);

    let left_clipped = matches!(ops.first()?.1, 'S' | 'H');
    let right_clipped = matches!(ops.last()?.1, 'S' | 'H');

    Some(match (left_clipped, right_clipped) {
        (true, false) => start,
        _ => end,
    })
}

fn parse_cigar_ops(cigar: &str) -> Vec<(usize, char)> {
    let mut ops = Vec::new();
    let mut length = String::new();

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            length.push(ch);
        } else {
            if let Ok(len) = length.parse::<usize>() {
                ops.push((len, ch));
            }
            length.clear();
        }
    }

    ops
}

fn create_tsv_writer(path: &Path) -> Result<BufWriter<File>> {
    let file = File::create(path)
        .with_context(|| format!("failed to create TSV output {}", path.display()))?;
    Ok(BufWriter::new(file))
}

fn write_tsv_header(writer: &Arc<Mutex<BufWriter<File>>>) -> Result<()> {
    let mut writer = writer
        .lock()
        .map_err(|_| anyhow!("TSV writer lock was poisoned"))?;
    writeln!(
        writer,
        "query_gene\tmatched_partner_gene\tquery_transcript_id\tpartner_transcript_id\tbreakpoint_estimate\tread_name\tread_flags\treference_name\talignment_start\talignment_end\tcigar\tmapping_quality\tmate_reference_name\tmate_alignment_start\tinferred_partner_reference\tinferred_partner_start\tinferred_partner_strand\tsa_tag\tsample"
    )?;
    Ok(())
}

fn write_tsv_row(tsv_writer: &Arc<Mutex<BufWriter<File>>>, hit: &Hit) -> Result<()> {
    {
        let mut writer = tsv_writer
            .lock()
            .map_err(|_| anyhow!("TSV writer lock was poisoned"))?;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            hit.tsv.query_gene,
            hit.tsv.matched_partner_gene.clone().unwrap_or_default(),
            hit.tsv.query_transcript_id,
            hit.tsv.partner_transcript_id,
            hit.tsv.breakpoint_estimate,
            hit.tsv.read_name,
            hit.tsv.read_flags,
            hit.tsv.reference_name,
            hit.tsv.alignment_start,
            hit.tsv.alignment_end,
            hit.tsv.cigar,
            hit.tsv
                .mapping_quality
                .map(|value| value.to_string())
                .unwrap_or_default(),
            hit.tsv.mate_reference_name.clone().unwrap_or_default(),
            hit.tsv
                .mate_alignment_start
                .map(|value| value.to_string())
                .unwrap_or_default(),
            hit.tsv.inferred_partner_reference,
            hit.tsv.inferred_partner_start,
            hit.tsv.inferred_partner_strand,
            hit.tsv.sa_tag,
            hit.tsv.sample
        )?;
    }

    Ok(())
}

/// Write one FASTA record for a read. Called once per record even when the read
/// supports several junctions, so the sequence is not duplicated.
fn write_fasta_record(
    fasta_writer: &Arc<Mutex<FastaWriter>>,
    header: &str,
    sequence: &str,
) -> Result<()> {
    let mut writer = fasta_writer
        .lock()
        .map_err(|_| anyhow!("FASTA writer lock was poisoned"))?;
    writer.write_record(header, sequence)?;
    Ok(())
}

struct Hit {
    tsv: TsvRecord,
    fasta_header: String,
    junction: Junction,
}

/// Everything one alignment record contributes: a hit per supplementary
/// alignment, plus the read sequence held once.
///
/// The sequence lives here rather than on each [`Hit`] because a chimeric long
/// read can carry many `SA` entries, and copying a multi-kilobase sequence per
/// junction is expensive in both memory and FASTA volume.
struct RecordHits {
    sequence: String,
    hits: Vec<Hit>,
}

struct PartnerAlignment {
    reference_name: String,
    start: usize,
    strand: char,
    breakpoint: usize,
}

#[cfg(test)]
mod tests {
    use super::{
        bam_basename_token, check_unique_sample_names, default_output_path,
        estimate_sa_breakpoint_position, filter_spans_by_gene_names, find_overlapping_span,
        genes_basename_token, has_bam_extension, parse_sa_entries, resolve_bam_inputs,
        resolve_partner_filter, sample_name, sanitize_token,
    };
    use crate::annotation::{Exon, GeneSpan, Transcript};
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};

    fn unique_temp_dir() -> PathBuf {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let suffix = COUNTER.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("stellerator_unit_{nanos}_{suffix}"));
        std::fs::create_dir(&dir).unwrap();
        dir
    }

    #[test]
    fn sample_name_uses_file_stem() {
        assert_eq!(sample_name(Path::new("/data/sampleA.bam")), "sampleA");
        assert_eq!(
            sample_name(Path::new("/data/sampleA.sorted.bam")),
            "sampleA.sorted"
        );
    }

    #[test]
    fn has_bam_extension_is_case_insensitive() {
        assert!(has_bam_extension(Path::new("a.bam")));
        assert!(has_bam_extension(Path::new("a.BAM")));
        assert!(!has_bam_extension(Path::new("a.bai")));
        assert!(!has_bam_extension(Path::new("notes.txt")));
    }

    #[test]
    fn check_unique_sample_names_detects_collisions() {
        let colliding = [
            PathBuf::from("/a/sample.bam"),
            PathBuf::from("/b/sample.bam"),
        ];
        let error = check_unique_sample_names(&colliding).unwrap_err();
        assert!(error.to_string().contains("duplicate sample name"));

        let distinct = [PathBuf::from("/a/one.bam"), PathBuf::from("/b/two.bam")];
        assert!(check_unique_sample_names(&distinct).is_ok());
    }

    #[test]
    fn resolve_bam_inputs_requires_index() {
        let dir = unique_temp_dir();
        let bam = dir.join("a.bam");
        std::fs::write(&bam, b"").unwrap();

        let error = resolve_bam_inputs(std::slice::from_ref(&bam)).unwrap_err();
        assert!(error.to_string().contains("indexed BAM required"));

        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn resolve_bam_inputs_expands_directory_and_dedups() {
        let dir = unique_temp_dir();
        for name in ["b.bam", "a.bam"] {
            std::fs::write(dir.join(name), b"").unwrap();
            std::fs::write(dir.join(format!("{name}.bai")), b"").unwrap();
        }
        std::fs::write(dir.join("notes.txt"), b"").unwrap();

        // Directory plus an explicit duplicate of one of its BAMs.
        let resolved = resolve_bam_inputs(&[dir.clone(), dir.join("a.bam")]).unwrap();
        let names: Vec<_> = resolved
            .iter()
            .map(|path| path.file_name().unwrap().to_str().unwrap().to_string())
            .collect();
        assert_eq!(names, vec!["a.bam", "b.bam"]);

        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn resolve_bam_inputs_rejects_empty_directory() {
        let dir = unique_temp_dir();
        let error = resolve_bam_inputs(std::slice::from_ref(&dir)).unwrap_err();
        assert!(error.to_string().contains("no .bam files found"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn genes_token_joins_genes_and_appends_distinct_partner() {
        assert_eq!(genes_basename_token(&["BCR".to_string()], None), "BCR");
        assert_eq!(
            genes_basename_token(&["BCR".to_string(), "ABL1".to_string()], None),
            "BCR_ABL1"
        );
        assert_eq!(
            genes_basename_token(&["BCR".to_string()], Some("ABL1")),
            "BCR_ABL1"
        );
        // A partner already present as a query gene is not repeated.
        assert_eq!(
            genes_basename_token(&["BCR".to_string(), "ABL1".to_string()], Some("abl1")),
            "BCR_ABL1"
        );
    }

    #[test]
    fn bam_token_uses_stem_then_shared_parent_directory() {
        // A single BAM contributes its own stem.
        assert_eq!(
            bam_basename_token(
                &[PathBuf::from("/data/sampleA.bam")],
                &["sampleA".to_string()]
            ),
            "sampleA"
        );

        // Several BAMs in one directory collapse to that directory name.
        assert_eq!(
            bam_basename_token(
                &[
                    PathBuf::from("/data/cohort/a.bam"),
                    PathBuf::from("/data/cohort/b.bam"),
                ],
                &["a".to_string(), "b".to_string()]
            ),
            "cohort"
        );

        // With no shared parent, fall back to the first sample name.
        assert_eq!(
            bam_basename_token(
                &[PathBuf::from("/x/a.bam"), PathBuf::from("/y/b.bam")],
                &["a".to_string(), "b".to_string()]
            ),
            "a"
        );
    }

    #[test]
    fn default_output_path_combines_bam_and_gene_tokens() {
        assert_eq!(
            default_output_path("sampleA", "BCR_ABL1", "tsv"),
            PathBuf::from("sampleA.BCR_ABL1.tsv")
        );
        assert_eq!(
            default_output_path("cohort", "BCR", "fasta.gz"),
            PathBuf::from("cohort.BCR.fasta.gz")
        );
    }

    #[test]
    fn sanitize_token_replaces_path_hostile_characters() {
        assert_eq!(sanitize_token("BCR/ABL1"), "BCR_ABL1");
        assert_eq!(sanitize_token("gene with space"), "gene_with_space");
        // Hyphens and dots are already filename-safe and are preserved.
        assert_eq!(sanitize_token("HLA-DRB1.v2"), "HLA-DRB1.v2");
    }

    #[test]
    fn parses_every_sa_entry() {
        let partners = parse_sa_entries("chr9,420,-,20S80M,60,0;chr1,100,+,30M70S,55,1;");
        assert_eq!(partners.len(), 2);

        assert_eq!(partners[0].reference_name, "chr9");
        assert_eq!(partners[0].start, 420);
        assert_eq!(partners[0].strand, '-');
        // 20S80M is left-clipped, so the breakpoint sits at the aligned start.
        assert_eq!(partners[0].breakpoint, 420);

        assert_eq!(partners[1].reference_name, "chr1");
        // 30M70S is right-clipped, so the breakpoint sits at the aligned end.
        assert_eq!(partners[1].breakpoint, 129);
    }

    #[test]
    fn sa_breakpoint_follows_clip_side() {
        assert_eq!(estimate_sa_breakpoint_position(500, "40S60M"), Some(500));
        assert_eq!(estimate_sa_breakpoint_position(500, "60M40S"), Some(559));
        // Reference-consuming ops (D/N) extend the aligned end.
        assert_eq!(
            estimate_sa_breakpoint_position(500, "10M5D10M20S"),
            Some(524)
        );
    }

    #[test]
    fn finds_overlapping_partner_gene() {
        let spans = vec![
            GeneSpan {
                gene: "BCR".to_string(),
                reference_name: "chr22".to_string(),
                start: 100,
                end: 200,
                strand: Some('+'),
                transcripts: vec![Transcript {
                    id: "tx1".to_string(),
                    exons: vec![Exon {
                        start: 100,
                        end: 200,
                    }],
                }],
            },
            GeneSpan {
                gene: "ABL1".to_string(),
                reference_name: "chr9".to_string(),
                start: 300,
                end: 450,
                strand: Some('+'),
                transcripts: vec![Transcript {
                    id: "tx1".to_string(),
                    exons: vec![Exon {
                        start: 300,
                        end: 450,
                    }],
                }],
            },
        ];

        assert_eq!(
            find_overlapping_span(&spans, "chr9", 350).map(|span| span.gene.clone()),
            Some("ABL1".to_string())
        );
        assert_eq!(find_overlapping_span(&spans, "chr9", 900), None);
    }

    #[test]
    fn partner_filter_is_none_without_a_partner_gene() {
        let spans = vec![
            GeneSpan {
                gene: "BCR".to_string(),
                reference_name: "chr22".to_string(),
                start: 100,
                end: 200,
                strand: Some('+'),
                transcripts: vec![],
            },
            GeneSpan {
                gene: "ABL1".to_string(),
                reference_name: "chr9".to_string(),
                start: 300,
                end: 450,
                strand: Some('+'),
                transcripts: vec![],
            },
        ];

        // No partner gene => no filter, so the partner side is annotated against
        // every span.
        assert!(resolve_partner_filter(&spans, None).is_none());

        // A partner gene restricts the partner side to that gene's spans.
        let filtered = resolve_partner_filter(&spans, Some("ABL1")).expect("a filter");
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].gene, "ABL1");
    }

    #[test]
    fn restricts_query_spans_to_requested_genes() {
        let spans = vec![
            GeneSpan {
                gene: "BCR".to_string(),
                reference_name: "chr22".to_string(),
                start: 100,
                end: 200,
                strand: Some('+'),
                transcripts: vec![],
            },
            GeneSpan {
                gene: "ABL1".to_string(),
                reference_name: "chr9".to_string(),
                start: 300,
                end: 450,
                strand: Some('+'),
                transcripts: vec![],
            },
        ];

        let filtered = filter_spans_by_gene_names(&spans, &["BCR".to_string()]);
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].gene, "BCR");
    }
}
