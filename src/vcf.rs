use std::{
    collections::{BTreeMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::{Context, Result};

/// How a supporting read connects the two breakends.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Evidence {
    /// A supplementary (`SA`) alignment crosses the junction.
    Split,
    /// The read's mate maps to the partner locus.
    Discordant,
}

impl Evidence {
    pub fn as_str(self) -> &'static str {
        match self {
            Evidence::Split => "split",
            Evidence::Discordant => "discordant",
        }
    }
}

/// Per-sample supporting-read counts, split by evidence class.
#[derive(Debug, Clone, Copy, Default)]
pub struct SampleSupport {
    pub split: usize,
    pub discordant: usize,
}

/// A single supporting-read breakend observation produced by the pipeline.
///
/// Each candidate fusion-supporting read contributes one junction per
/// supplementary (`SA`) alignment that escapes the queried gene interval, plus
/// one more when its mate provides discordant-pair evidence.
#[derive(Clone)]
pub struct Junction {
    pub sample: String,
    pub read_name: String,
    pub evidence: Evidence,
    pub query_gene: String,
    pub partner_gene: Option<String>,
    pub query_transcript: String,
    pub partner_transcript: String,
    pub query_region: String,
    pub partner_region: String,
    pub chrom1: String,
    pub pos1: usize,
    pub strand1: char,
    pub chrom2: String,
    pub pos2: usize,
    pub strand2: char,
}

/// A consensus structural variant clustered from one or more [`Junction`]s.
pub struct StructuralVariant {
    pub chrom1: String,
    pub pos1: usize,
    pub strand1: char,
    pub chrom2: String,
    pub pos2: usize,
    pub strand2: char,
    pub gene1: String,
    pub gene2: String,
    pub transcript1: String,
    pub transcript2: String,
    pub region: String,
    pub support_split: usize,
    pub support_discordant: usize,
    pub support_by_sample: BTreeMap<String, SampleSupport>,
}

/// Discrete key that breakends must share before position-tolerance clustering.
type JunctionKey = (String, String, char, char);

fn junction_key(junction: &Junction) -> JunctionKey {
    (
        junction.chrom1.clone(),
        junction.chrom2.clone(),
        junction.strand1,
        junction.strand2,
    )
}

/// Cluster raw breakend observations into consensus structural variants.
///
/// Two breakends join the same call when they share both chromosomes and both
/// strands and their breakpoints each fall within `slop` bp of the cluster
/// anchor. Support is the count of distinct supporting reads (a read that hits
/// the same junction more than once is counted once).
pub fn cluster_consensus(mut junctions: Vec<Junction>, slop: usize) -> Vec<StructuralVariant> {
    junctions.sort_by(|a, b| {
        junction_key(a)
            .cmp(&junction_key(b))
            .then(a.pos1.cmp(&b.pos1))
            .then(a.pos2.cmp(&b.pos2))
    });

    let mut variants = Vec::new();
    let mut cluster: Vec<Junction> = Vec::new();

    for junction in junctions {
        let starts_new = match cluster.first() {
            None => false,
            Some(anchor) => {
                junction_key(anchor) != junction_key(&junction)
                    || junction.pos1.abs_diff(anchor.pos1) > slop
                    || junction.pos2.abs_diff(anchor.pos2) > slop
            }
        };

        if starts_new {
            variants.push(finalize_cluster(&cluster));
            cluster.clear();
        }
        cluster.push(junction);
    }

    if !cluster.is_empty() {
        variants.push(finalize_cluster(&cluster));
    }

    variants.sort_by(|a, b| {
        a.chrom1
            .cmp(&b.chrom1)
            .then(a.pos1.cmp(&b.pos1))
            .then(a.chrom2.cmp(&b.chrom2))
            .then(a.pos2.cmp(&b.pos2))
    });

    variants
}

fn finalize_cluster(cluster: &[Junction]) -> StructuralVariant {
    let anchor = &cluster[0];

    let mut support_by_sample: BTreeMap<String, SampleSupport> = BTreeMap::new();
    let mut seen: HashSet<(&str, &str, Evidence)> = HashSet::new();
    let mut support_split = 0;
    let mut support_discordant = 0;
    for junction in cluster {
        // A read counts once per evidence class, so a read that is both split
        // and discordant contributes to each tally exactly once.
        if seen.insert((
            junction.sample.as_str(),
            junction.read_name.as_str(),
            junction.evidence,
        )) {
            let entry = support_by_sample
                .entry(junction.sample.clone())
                .or_default();
            match junction.evidence {
                Evidence::Split => {
                    entry.split += 1;
                    support_split += 1;
                }
                Evidence::Discordant => {
                    entry.discordant += 1;
                    support_discordant += 1;
                }
            }
        }
    }

    StructuralVariant {
        chrom1: anchor.chrom1.clone(),
        pos1: median(cluster.iter().map(|junction| junction.pos1)),
        strand1: anchor.strand1,
        chrom2: anchor.chrom2.clone(),
        pos2: median(cluster.iter().map(|junction| junction.pos2)),
        strand2: anchor.strand2,
        gene1: anchor.query_gene.clone(),
        gene2: anchor
            .partner_gene
            .clone()
            .unwrap_or_else(|| "NA".to_string()),
        transcript1: anchor.query_transcript.clone(),
        transcript2: anchor.partner_transcript.clone(),
        region: format!("{}/{}", anchor.query_region, anchor.partner_region),
        support_split,
        support_discordant,
        support_by_sample,
    }
}

fn median(values: impl Iterator<Item = usize>) -> usize {
    let mut values: Vec<usize> = values.collect();
    values.sort_unstable();
    let mid = values.len() / 2;
    if values.len().is_multiple_of(2) {
        (values[mid - 1] + values[mid]) / 2
    } else {
        values[mid]
    }
}

/// Write the consensus structural variants as a multi-sample VCF. Per-sample
/// supporting-read counts populate the `SR` genotype field; `samples` defines
/// the (deterministic) column order.
pub fn write_vcf(
    path: &Path,
    variants: &[StructuralVariant],
    samples: &[String],
    contigs: &[(String, Option<usize>)],
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create VCF output {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##source=Stellerator-{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(
        writer,
        "##ALT=<ID=BND,Description=\"Breakend / fusion junction\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome of the mate breakend\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Position of the mate breakend\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strand orientation of the query and partner breakends\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=GENE1,Number=1,Type=String,Description=\"Gene at the query breakend\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=GENE2,Number=1,Type=String,Description=\"Gene at the partner breakend (NA if unannotated)\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=TRANSCRIPT1,Number=1,Type=String,Description=\"Transcript used for the query breakend label\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=TRANSCRIPT2,Number=1,Type=String,Description=\"Transcript used for the partner breakend label\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=REGION,Number=1,Type=String,Description=\"Breakpoint region labels (query/partner)\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split reads supporting the junction across all samples\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Discordant read pairs supporting the junction across all samples\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Split reads supporting the junction in the sample\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Discordant read pairs supporting the junction in the sample\">"
    )?;
    for (name, length) in contigs {
        match length {
            Some(length) => writeln!(writer, "##contig=<ID={name},length={length}>")?,
            None => writeln!(writer, "##contig=<ID={name}>")?,
        }
    }

    write!(
        writer,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    )?;
    for sample in samples {
        write!(writer, "\t{sample}")?;
    }
    writeln!(writer)?;

    for (index, variant) in variants.iter().enumerate() {
        let info = format!(
            "SVTYPE=BND;CHR2={};POS2={};STRANDS={}{};GENE1={};GENE2={};TRANSCRIPT1={};TRANSCRIPT2={};REGION={};SR={};PE={}",
            variant.chrom2,
            variant.pos2,
            variant.strand1,
            variant.strand2,
            variant.gene1,
            variant.gene2,
            variant.transcript1,
            variant.transcript2,
            variant.region,
            variant.support_split,
            variant.support_discordant,
        );
        write!(
            writer,
            "{}\t{}\tSTL_BND_{}\tN\t<BND>\t.\tPASS\t{}\tSR:PE",
            variant.chrom1,
            variant.pos1,
            index + 1,
            info,
        )?;
        for sample in samples {
            let support = variant
                .support_by_sample
                .get(sample)
                .copied()
                .unwrap_or_default();
            write!(writer, "\t{}:{}", support.split, support.discordant)?;
        }
        writeln!(writer)?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn junction(sample: &str, read: &str, pos1: usize, pos2: usize) -> Junction {
        evidence_junction(sample, read, pos1, pos2, Evidence::Split)
    }

    fn evidence_junction(
        sample: &str,
        read: &str,
        pos1: usize,
        pos2: usize,
        evidence: Evidence,
    ) -> Junction {
        Junction {
            sample: sample.to_string(),
            read_name: read.to_string(),
            evidence,
            query_gene: "BCR".to_string(),
            partner_gene: Some("ABL1".to_string()),
            query_transcript: "txBCR".to_string(),
            partner_transcript: "txABL1".to_string(),
            query_region: "exon1".to_string(),
            partner_region: "exon1".to_string(),
            chrom1: "chr22".to_string(),
            pos1,
            strand1: '+',
            chrom2: "chr9".to_string(),
            pos2,
            strand2: '-',
        }
    }

    #[test]
    fn clusters_nearby_breakends_and_counts_distinct_reads() {
        let junctions = vec![
            junction("s1", "r1", 100, 400),
            junction("s1", "r2", 103, 402),
            junction("s2", "r3", 101, 401),
        ];

        let variants = cluster_consensus(junctions, 10);
        assert_eq!(variants.len(), 1);

        let variant = &variants[0];
        assert_eq!(variant.support_split, 3);
        assert_eq!(variant.support_discordant, 0);
        assert_eq!(variant.support_by_sample["s1"].split, 2);
        assert_eq!(variant.support_by_sample["s2"].split, 1);
        assert_eq!(variant.gene1, "BCR");
        assert_eq!(variant.gene2, "ABL1");
        assert_eq!(variant.pos1, 101);
        assert_eq!(variant.pos2, 401);
    }

    #[test]
    fn separates_breakends_beyond_slop() {
        let junctions = vec![
            junction("s1", "r1", 100, 400),
            junction("s1", "r2", 200, 400),
        ];

        let variants = cluster_consensus(junctions, 10);
        assert_eq!(variants.len(), 2);
    }

    #[test]
    fn deduplicates_reads_within_a_cluster() {
        // A single read with two SA hits to the same locus counts once.
        let junctions = vec![
            junction("s1", "r1", 100, 400),
            junction("s1", "r1", 101, 401),
        ];

        let variants = cluster_consensus(junctions, 10);
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].support_split, 1);
        assert_eq!(variants[0].support_by_sample["s1"].split, 1);
    }

    #[test]
    fn tallies_split_and_discordant_support_separately() {
        // The same read supplies both evidence classes for one junction; each
        // class counts it once.
        let junctions = vec![
            evidence_junction("s1", "r1", 100, 400, Evidence::Split),
            evidence_junction("s1", "r1", 101, 401, Evidence::Discordant),
            evidence_junction("s1", "r2", 102, 402, Evidence::Discordant),
        ];

        let variants = cluster_consensus(junctions, 10);
        assert_eq!(variants.len(), 1);

        let variant = &variants[0];
        assert_eq!(variant.support_split, 1);
        assert_eq!(variant.support_discordant, 2);
        assert_eq!(variant.support_by_sample["s1"].split, 1);
        assert_eq!(variant.support_by_sample["s1"].discordant, 2);
    }

    #[test]
    fn separates_distinct_partner_chromosomes() {
        let mut other = junction("s1", "r2", 100, 400);
        other.chrom2 = "chr1".to_string();

        let variants = cluster_consensus(vec![junction("s1", "r1", 100, 400), other], 10);
        assert_eq!(variants.len(), 2);
    }

    #[test]
    fn median_handles_even_and_odd_counts() {
        assert_eq!(median([10usize, 30, 20].into_iter()), 20);
        assert_eq!(median([10usize, 20].into_iter()), 15);
    }
}
