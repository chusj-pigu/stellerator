use super::*;

pub(super) fn process_span(
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

pub(super) fn build_region(span: &GeneSpan) -> Result<Region> {
    let start = Position::try_from(span.start as usize)
        .map_err(|_| anyhow!("invalid start coordinate for {}", span.gene))?;
    let end = Position::try_from(span.end as usize)
        .map_err(|_| anyhow!("invalid end coordinate for {}", span.gene))?;
    Ok(Region::new(span.reference_name.clone(), start..=end))
}

/// Per-read fields shared across every supplementary alignment of one record.
pub(super) struct ReadContext {
    pub(super) read_name: String,
    pub(super) sequence: String,
    pub(super) read_flags: u16,
    pub(super) query_strand: char,
    pub(super) reference_name: String,
    pub(super) alignment_start: usize,
    pub(super) alignment_end: usize,
    pub(super) query_breakpoint_position: usize,
    pub(super) cigar: String,
    pub(super) mapping_quality: Option<u8>,
    pub(super) mate_reference_name: Option<String>,
    pub(super) mate_alignment_start: Option<usize>,
    pub(super) sa_tag: String,
}

/// Read-filtering options applied to every scanned record.
#[derive(Debug, Clone, Copy)]
pub(super) struct ScanOptions {
    pub(super) include_duplicates: bool,
    pub(super) min_mapq: u8,
}

impl ScanOptions {
    pub(super) fn from_args(args: &Args) -> Self {
        Self {
            include_duplicates: args.include_duplicates,
            min_mapq: args.min_mapq,
        }
    }
}

pub(super) fn classify_record(
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
pub(super) fn read_context(
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
pub(super) fn build_hit(
    sample: &str,
    span: &GeneSpan,
    partner_spans: Option<&[GeneSpan]>,
    require_partner_match: bool,
    context: &ReadContext,
    partner: PartnerAlignment,
) -> Option<Hit> {
    let matched_partner_span = partner_spans.and_then(|spans| {
        find_overlapping_partner(
            spans,
            &partner.reference_name,
            partner.breakpoint,
            &span.gene,
        )
        .or_else(|| {
            find_overlapping_partner(spans, &partner.reference_name, partner.start, &span.gene)
        })
    });

    if matched_partner_span.is_none()
        && partner.reference_name == span.reference_name
        && partner.start >= span.start as usize
        && partner.start <= span.end as usize
    {
        return None;
    }

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

pub(super) fn extract_sa(record: &bam::Record) -> Option<String> {
    let data = record.data();
    let value = data.get(&Tag::OTHER_ALIGNMENTS)?.ok()?;

    match value {
        Value::String(raw) => Some(String::from_utf8_lossy(raw.as_ref()).into_owned()),
        _ => None,
    }
}

pub(super) fn reference_name_for_id(header: &sam::Header, id: usize) -> Option<String> {
    header
        .reference_sequences()
        .get_index(id)
        .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned())
}

pub(super) fn filter_spans_by_gene_names(
    all_spans: &[GeneSpan],
    names: &[String],
) -> Vec<GeneSpan> {
    let wanted: Vec<String> = names.iter().map(|name| name.to_ascii_lowercase()).collect();

    all_spans
        .iter()
        .filter(|span| wanted.contains(&span.gene.to_ascii_lowercase()))
        .cloned()
        .collect()
}

pub(super) fn find_overlapping_partner<'a>(
    spans: &'a [GeneSpan],
    reference_name: &str,
    position: usize,
    query_gene: &str,
) -> Option<&'a GeneSpan> {
    spans.iter().find(|span| {
        !span.gene.eq_ignore_ascii_case(query_gene)
            && span.reference_name == reference_name
            && position >= span.start as usize
            && position <= span.end as usize
    })
}

pub(super) struct Hit {
    pub(super) tsv: TsvRecord,
    pub(super) fasta_header: String,
    pub(super) junction: Junction,
}

/// Everything one alignment record contributes: a hit per supplementary
/// alignment, plus the read sequence held once.
///
/// The sequence lives here rather than on each [`Hit`] because a chimeric long
/// read can carry many `SA` entries, and copying a multi-kilobase sequence per
/// junction is expensive in both memory and FASTA volume.
pub(super) struct RecordHits {
    pub(super) sequence: String,
    pub(super) hits: Vec<Hit>,
}

pub(super) struct PartnerAlignment {
    pub(super) reference_name: String,
    pub(super) start: usize,
    pub(super) strand: char,
    pub(super) breakpoint: usize,
}
