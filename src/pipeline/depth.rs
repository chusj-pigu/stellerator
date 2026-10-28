use super::*;

/// Count reads spanning a breakpoint in one sample, applying the same record
/// filters as the main scan.
///
/// This is the denominator behind allele fraction: without it a support count
/// cannot distinguish a clonal event from a handful of artefacts. Reads are
/// counted by name so a read with several alignments over the position counts
/// once, and supplementary alignments are skipped for the same reason.
pub(super) fn sample_depths(
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
pub(super) fn sample_names(samples: &[BamSample]) -> Vec<String> {
    let mut names: Vec<String> = samples.iter().map(|sample| sample.name.clone()).collect();
    names.sort();
    names
}

/// Union of reference sequences across all sample headers, in first-seen order,
/// used to emit `##contig` lines in the VCF.
pub(super) fn union_contigs(samples: &[BamSample]) -> Vec<(String, Option<usize>)> {
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
