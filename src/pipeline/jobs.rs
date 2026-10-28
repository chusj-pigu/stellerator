use super::*;

/// Shared, read-only state for scanning gene spans: where results go and how
/// records are filtered. Held once per run and borrowed by every worker.
/// One unit of work: query genes, an optional partner constraint, and the
/// clustering tolerance for the calls it produces. Single-CLI mode is one job;
/// `--loci` mode is one per row.
pub(super) struct Job {
    pub(super) query_genes: Vec<String>,
    pub(super) partner_gene: Option<String>,
    pub(super) slop: usize,
}

/// A job with its annotation spans resolved and its own junction collector,
/// ready to scan.
pub(super) struct JobPlan {
    pub(super) query_spans: Vec<GeneSpan>,
    /// Spans the partner side is restricted to; `None` annotates against every
    /// span (the `--partner-gene`-omitted behaviour) and requires no match.
    pub(super) partner_filter: Option<Vec<GeneSpan>>,
    pub(super) slop: usize,
    pub(super) junctions: Option<Mutex<Vec<Junction>>>,
}

/// Scan-wide state shared by every worker, independent of the job.
pub(super) struct ScanContext<'a> {
    pub(super) all_spans: &'a [GeneSpan],
    pub(super) tsv_writer: &'a Arc<Mutex<BufWriter<File>>>,
    pub(super) fasta_writer: &'a Arc<Mutex<FastaWriter>>,
    pub(super) options: ScanOptions,
}

pub(super) fn build_jobs(args: &Args) -> Result<Vec<Job>> {
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
pub(super) fn referenced_genes(jobs: &[Job]) -> Vec<String> {
    let mut names = Vec::new();
    for job in jobs {
        names.extend(job.query_genes.iter().cloned());
        if let Some(partner) = &job.partner_gene {
            names.push(partner.clone());
        }
    }
    names
}

pub(super) fn build_job_plan(
    job: &Job,
    all_spans: &[GeneSpan],
    collect_junctions: bool,
) -> JobPlan {
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
pub(super) fn resolve_partner_filter(
    all_spans: &[GeneSpan],
    partner_gene: Option<&str>,
) -> Option<Vec<GeneSpan>> {
    partner_gene.map(|partner| {
        filter_spans_by_gene_names(all_spans, std::slice::from_ref(&partner.to_string()))
    })
}
