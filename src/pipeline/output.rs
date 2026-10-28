use super::*;

/// Resolve an output path: use `explicit` when provided, otherwise build a
/// default from the BAM basename and requested genes (e.g. `sample.BCR_ABL1.tsv`).
pub(super) fn resolve_output_path(
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

/// Check every output before logging or result writers can truncate a file.
pub(super) fn validate_output_paths(args: &Args) -> Result<()> {
    let samples = open_bam_samples(&args.bam)?;
    let genes_token = output_genes_token(args);
    let mut outputs = vec![
        resolve_output_path(args.output_tsv.as_deref(), &samples, &genes_token, "tsv"),
        resolve_output_path(
            args.output_fasta.as_deref(),
            &samples,
            &genes_token,
            "fasta.gz",
        ),
    ];
    if let Some(vcf) = &args.output_vcf {
        outputs.push(resolve_output_path(
            vcf.as_deref(),
            &samples,
            &genes_token,
            "vcf",
        ));
    }
    if let Some(log) = &args.log_file {
        outputs.push(log.clone());
    }

    let mut inputs = vec![args.annotation.clone()];
    inputs.extend(args.loci.iter().cloned());
    for sample in &samples {
        inputs.push(sample.path.clone());
        let mut bai = sample.path.as_os_str().to_os_string();
        bai.push(".bai");
        let mut csi = sample.path.as_os_str().to_os_string();
        csi.push(".csi");
        for index in [
            PathBuf::from(bai),
            PathBuf::from(csi),
            sample.path.with_extension("bai"),
            sample.path.with_extension("csi"),
        ] {
            if index.exists() {
                inputs.push(index);
            }
        }
    }

    for (index, output) in outputs.iter().enumerate() {
        for other in inputs.iter().chain(outputs[..index].iter()) {
            if paths_alias(output, other)? {
                bail!(
                    "output path {} conflicts with {}",
                    output.display(),
                    other.display()
                );
            }
        }
    }
    Ok(())
}

pub(super) fn paths_alias(a: &Path, b: &Path) -> Result<bool> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(left), Ok(right)) = (fs::metadata(a), fs::metadata(b))
            && left.dev() == right.dev()
            && left.ino() == right.ino()
        {
            return Ok(true);
        }
    }
    Ok(canonical_output_path(a)? == canonical_output_path(b)?)
}

pub(super) fn canonical_output_path(path: &Path) -> Result<PathBuf> {
    if fs::symlink_metadata(path).is_ok() {
        return fs::canonicalize(path)
            .with_context(|| format!("failed to resolve path {}", path.display()));
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent = fs::canonicalize(parent)
        .with_context(|| format!("failed to resolve parent of {}", path.display()))?;
    let name = path
        .file_name()
        .ok_or_else(|| anyhow!("invalid path {}", path.display()))?;
    Ok(parent.join(name))
}

/// The gene component of default output names: the requested genes in single
/// mode, or the loci file stem in batch mode.
pub(super) fn output_genes_token(args: &Args) -> String {
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

pub(super) fn default_output_path(bam_token: &str, genes_token: &str, extension: &str) -> PathBuf {
    PathBuf::from(format!("{bam_token}.{genes_token}.{extension}"))
}

/// Representative BAM basename for default output names: the sample stem for a
/// single BAM, the shared parent directory name for several BAMs in one place,
/// otherwise the first sample's name.
pub(super) fn bam_basename_token(paths: &[PathBuf], names: &[String]) -> String {
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
pub(super) fn genes_basename_token(genes: &[String], partner_gene: Option<&str>) -> String {
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

pub(super) fn push_unique_token(parts: &mut Vec<String>, value: &str) {
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
pub(super) fn sanitize_token(value: &str) -> String {
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

pub(super) fn create_tsv_writer(path: &Path) -> Result<BufWriter<File>> {
    let file = File::create(path)
        .with_context(|| format!("failed to create TSV output {}", path.display()))?;
    Ok(BufWriter::new(file))
}

pub(super) fn write_tsv_header(writer: &Arc<Mutex<BufWriter<File>>>) -> Result<()> {
    let mut writer = writer
        .lock()
        .map_err(|_| anyhow!("TSV writer lock was poisoned"))?;
    writeln!(
        writer,
        "query_gene\tmatched_partner_gene\tquery_transcript_id\tpartner_transcript_id\tbreakpoint_estimate\tread_name\tread_flags\treference_name\talignment_start\talignment_end\tcigar\tmapping_quality\tmate_reference_name\tmate_alignment_start\tinferred_partner_reference\tinferred_partner_start\tinferred_partner_strand\tsa_tag\tsample"
    )?;
    Ok(())
}

pub(super) fn write_tsv_row(tsv_writer: &Arc<Mutex<BufWriter<File>>>, hit: &Hit) -> Result<()> {
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
pub(super) fn write_fasta_record(
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
