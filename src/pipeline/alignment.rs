use super::*;

pub(super) fn estimate_query_breakpoint_position(
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

pub(super) fn cigar_to_string(cigar: &bam::record::Cigar<'_>) -> Result<String, std::io::Error> {
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
pub(super) fn parse_sa_entries(raw: &str) -> Vec<PartnerAlignment> {
    raw.split(';')
        .filter(|entry| !entry.trim().is_empty())
        .filter_map(parse_sa_entry)
        .collect()
}

pub(super) fn parse_sa_entry(entry: &str) -> Option<PartnerAlignment> {
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
pub(super) fn estimate_sa_breakpoint_position(start: usize, cigar: &str) -> Option<usize> {
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

pub(super) fn parse_cigar_ops(cigar: &str) -> Vec<(usize, char)> {
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
