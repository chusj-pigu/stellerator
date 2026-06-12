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
use tracing::{debug, info};

use crate::{
    annotation::{GeneSpan, breakpoint_annotation, load_target_spans},
    cli::Args,
    fasta::FastaWriter,
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

    let samples = open_bam_samples(&args.bam)?;
    info!("processing {} BAM sample(s)", samples.len());

    let requested_genes = requested_gene_names(&args);
    let annotation_genes = if args.partner_gene.is_some() {
        requested_genes.clone()
    } else {
        Vec::new()
    };
    let all_spans = load_target_spans(&args.annotation, &annotation_genes)?;
    let query_spans = query_spans(&all_spans, &args);
    let partner_spans = partner_spans(&all_spans, &args);
    let require_partner_match = args.partner_gene.is_some();
    info!("loaded {} query intervals", query_spans.len());

    let tsv_writer = Arc::new(Mutex::new(create_tsv_writer(&args.output_tsv)?));
    let fasta_writer = Arc::new(Mutex::new(FastaWriter::create(&args.output_fasta)?));

    write_tsv_header(&tsv_writer)?;

    let work: Vec<(&BamSample, &GeneSpan)> = samples
        .iter()
        .flat_map(|sample| query_spans.iter().map(move |span| (sample, span)))
        .collect();

    work.par_iter().try_for_each(|&(sample, span)| {
        process_span(
            sample,
            span,
            partner_spans.as_deref(),
            require_partner_match,
            &tsv_writer,
            &fasta_writer,
        )
    })?;

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
        args.output_tsv.display(),
        args.output_fasta.display()
    );

    Ok(())
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

fn process_span(
    sample: &BamSample,
    span: &GeneSpan,
    partner_spans: Option<&[GeneSpan]>,
    require_partner_match: bool,
    tsv_writer: &Arc<Mutex<BufWriter<File>>>,
    fasta_writer: &Arc<Mutex<FastaWriter>>,
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

    let region = build_region(span)?;
    let query = reader.query(&sample.header, &region)?;
    for result in query.records() {
        let record = result?;
        if let Some(hit) = classify_record(
            &sample.header,
            &sample.name,
            span,
            partner_spans,
            require_partner_match,
            &record,
        ) {
            write_hit(tsv_writer, fasta_writer, &hit)?;
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

fn classify_record(
    header: &sam::Header,
    sample: &str,
    span: &GeneSpan,
    partner_spans: Option<&[GeneSpan]>,
    require_partner_match: bool,
    record: &bam::Record,
) -> Option<Hit> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() {
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

    let partner = parse_sa_entry(&sa_tag)?;
    if partner.reference_name == span.reference_name
        && partner.start >= span.start as usize
        && partner.start <= span.end as usize
    {
        return None;
    }

    let matched_partner_span = partner_spans.and_then(|spans| {
        find_overlapping_span(spans, &partner.reference_name, partner.start)
            .filter(|partner_span| partner_span.gene != span.gene)
    });

    if require_partner_match && matched_partner_span.is_none() {
        return None;
    }

    let query_breakpoint = breakpoint_annotation(span, query_breakpoint_position);
    let partner_breakpoint = matched_partner_span
        .and_then(|partner_span| breakpoint_annotation(partner_span, partner.start));
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

    Some(Hit {
        tsv: TsvRecord {
            query_gene: span.gene.clone(),
            matched_partner_gene: matched_partner_gene.clone(),
            query_transcript_id: query_transcript_id.clone(),
            partner_transcript_id: partner_transcript_id.clone(),
            breakpoint_estimate: breakpoint_estimate.clone(),
            read_name: read_name.clone(),
            read_flags: flags.bits(),
            reference_name,
            alignment_start,
            alignment_end,
            cigar,
            mapping_quality,
            mate_reference_name,
            mate_alignment_start,
            inferred_partner_reference: partner.reference_name.clone(),
            inferred_partner_start: partner.start,
            inferred_partner_strand: partner.strand.to_string(),
            sa_tag,
            sample: sample.to_string(),
        },
        fasta_header: format!(
            "{} gene={} matched_partner_gene={} query_transcript_id={} partner_transcript_id={} breakpoint_estimate={} partner={}:{} strand={} sample={}",
            read_name,
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
        fasta_sequence: sequence,
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

fn requested_gene_names(args: &Args) -> Vec<String> {
    let mut genes = args.genes.clone();

    if let Some(partner_gene) = &args.partner_gene {
        genes.push(partner_gene.clone());
    }

    genes
}

fn query_spans(all_spans: &[GeneSpan], args: &Args) -> Vec<GeneSpan> {
    filter_spans_by_gene_names(all_spans, &args.genes)
}

fn partner_spans(all_spans: &[GeneSpan], args: &Args) -> Option<Vec<GeneSpan>> {
    match args.partner_gene.as_ref() {
        Some(partner_gene) => Some(filter_spans_by_gene_names(
            all_spans,
            std::slice::from_ref(partner_gene),
        )),
        None => Some(all_spans.to_vec()),
    }
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

fn parse_sa_entry(raw: &str) -> Option<PartnerAlignment> {
    let first = raw.split(';').find(|entry| !entry.trim().is_empty())?;
    let mut fields = first.split(',');
    let reference_name = fields.next()?.to_string();
    let start = fields.next()?.parse().ok()?;
    let strand = fields.next()?.chars().next()?;
    Some(PartnerAlignment {
        reference_name,
        start,
        strand,
    })
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

fn write_hit(
    tsv_writer: &Arc<Mutex<BufWriter<File>>>,
    fasta_writer: &Arc<Mutex<FastaWriter>>,
    hit: &Hit,
) -> Result<()> {
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

    let mut writer = fasta_writer
        .lock()
        .map_err(|_| anyhow!("FASTA writer lock was poisoned"))?;
    writer.write_record(&hit.fasta_header, &hit.fasta_sequence)?;
    Ok(())
}

struct Hit {
    tsv: TsvRecord,
    fasta_header: String,
    fasta_sequence: String,
}

struct PartnerAlignment {
    reference_name: String,
    start: usize,
    strand: char,
}

#[cfg(test)]
mod tests {
    use super::{
        check_unique_sample_names, filter_spans_by_gene_names, find_overlapping_span,
        has_bam_extension, partner_spans, resolve_bam_inputs, sample_name,
    };
    use crate::annotation::{Exon, GeneSpan, Transcript};
    use crate::cli::Args;
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
    fn uses_all_annotation_spans_when_partner_gene_is_not_provided() {
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

        let args = Args {
            bam: vec![PathBuf::from("sample.bam")],
            annotation: PathBuf::from("genes.gtf"),
            genes: vec!["BCR".to_string()],
            partner_gene: None,
            output_tsv: PathBuf::from("out.tsv"),
            output_fasta: PathBuf::from("out.fa.gz"),
            threads: 1,
            verbose: false,
            log_file: None,
        };

        let partners = partner_spans(&spans, &args).unwrap();
        assert_eq!(partners.len(), 2);
        assert_eq!(
            find_overlapping_span(&partners, "chr9", 350).map(|span| span.gene.clone()),
            Some("ABL1".to_string())
        );
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

    #[test]
    fn does_not_add_partner_gene_to_query_spans() {
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

        let args = Args {
            bam: vec![PathBuf::from("sample.bam")],
            annotation: PathBuf::from("genes.gtf"),
            genes: vec!["BCR".to_string()],
            partner_gene: Some("ABL1".to_string()),
            output_tsv: PathBuf::from("out.tsv"),
            output_fasta: PathBuf::from("out.fa.gz"),
            threads: 1,
            verbose: false,
            log_file: None,
        };

        let queries = super::query_spans(&spans, &args);
        assert_eq!(queries.len(), 1);
        assert_eq!(queries[0].gene, "BCR");
    }
}
