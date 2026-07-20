use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::{Context, Result};

/// A single supporting-read breakend observation produced by the pipeline.
///
/// Each candidate fusion-supporting read contributes one junction per
/// supplementary (`SA`) alignment that escapes the queried gene interval.
#[derive(Clone)]
pub struct Junction {
    pub sample: String,
    pub read_name: String,
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
    /// Lowest and highest query breakpoint among the supporting reads, showing
    /// how far the cluster actually scattered.
    pub pos1_range: (usize, usize),
    /// Lowest and highest partner breakpoint among the supporting reads.
    pub pos2_range: (usize, usize),
    pub support_total: usize,
    /// Names of the reads supporting the junction, per sample. Kept rather than
    /// a bare count so depth can treat them as spanning by definition.
    pub support_reads_by_sample: BTreeMap<String, BTreeSet<String>>,
    /// Reads spanning the consensus breakpoint per sample, filled in after
    /// clustering by re-querying each BAM. Empty until then.
    pub depth_by_sample: BTreeMap<String, usize>,
}

impl StructuralVariant {
    /// Supporting reads for a sample (the ALT allele depth).
    pub fn support(&self, sample: &str) -> usize {
        self.support_reads_by_sample
            .get(sample)
            .map_or(0, |reads| reads.len())
    }

    /// Names of the reads supporting the junction in a sample.
    pub fn support_reads(&self, sample: &str) -> Option<&BTreeSet<String>> {
        self.support_reads_by_sample.get(sample)
    }

    /// Reads spanning the breakpoint in a sample (total depth).
    pub fn depth(&self, sample: &str) -> usize {
        self.depth_by_sample.get(sample).copied().unwrap_or(0)
    }

    /// Fraction of spanning reads that support the junction.
    pub fn allele_fraction(&self, sample: &str) -> f64 {
        let depth = self.depth(sample);
        if depth == 0 {
            return 0.0;
        }
        (self.support(sample) as f64 / depth as f64).min(1.0)
    }

    pub fn total_depth(&self) -> usize {
        self.depth_by_sample.values().sum()
    }

    pub fn total_allele_fraction(&self) -> f64 {
        let depth = self.total_depth();
        if depth == 0 {
            return 0.0;
        }
        (self.support_total as f64 / depth as f64).min(1.0)
    }
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

    // Collect names rather than counts: a read seen twice in the same cluster
    // still supports the junction once, and depth needs the names later.
    let mut support_reads_by_sample: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for junction in cluster {
        support_reads_by_sample
            .entry(junction.sample.clone())
            .or_default()
            .insert(junction.read_name.clone());
    }
    let support_total = support_reads_by_sample
        .values()
        .map(|reads| reads.len())
        .sum();

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
        pos1_range: position_range(cluster.iter().map(|junction| junction.pos1)),
        pos2_range: position_range(cluster.iter().map(|junction| junction.pos2)),
        support_total,
        support_reads_by_sample,
        depth_by_sample: BTreeMap::new(),
    }
}

/// Lowest and highest value in a cluster, used to expose breakpoint scatter.
fn position_range(values: impl Iterator<Item = usize>) -> (usize, usize) {
    let values: Vec<usize> = values.collect();
    let low = values.iter().copied().min().unwrap_or(0);
    let high = values.iter().copied().max().unwrap_or(0);
    (low, high)
}

/// Offsets from an anchor position to the ends of a range, as VCF confidence
/// intervals are expressed relative to the reported position.
fn range_offsets(anchor: usize, range: (usize, usize)) -> (i64, i64) {
    (
        range.0 as i64 - anchor as i64,
        range.1 as i64 - anchor as i64,
    )
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
        "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Total supporting reads across samples\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Offsets from POS to the lowest and highest supporting breakpoint in the cluster; the width shows the observed breakpoint scatter\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=CIPOS2,Number=2,Type=Integer,Description=\"Offsets from POS2 to the lowest and highest supporting partner breakpoint in the cluster\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Reads spanning the breakend across all samples\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Fraction of spanning reads supporting the junction across all samples\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype; nominal only, low-frequency fusions are not diploid states - use AF and AD\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Reads spanning the breakend in the sample\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Spanning reads not supporting, and supporting, the junction\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Fraction of spanning reads supporting the junction in the sample\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Supporting reads in the sample\">"
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
            "SVTYPE=BND;CHR2={};POS2={};STRANDS={}{};GENE1={};GENE2={};TRANSCRIPT1={};TRANSCRIPT2={};REGION={};SR={}",
            variant.chrom2,
            variant.pos2,
            variant.strand1,
            variant.strand2,
            variant.gene1,
            variant.gene2,
            variant.transcript1,
            variant.transcript2,
            variant.region,
            variant.support_total,
        );
        let (cipos_low, cipos_high) = range_offsets(variant.pos1, variant.pos1_range);
        let (cipos2_low, cipos2_high) = range_offsets(variant.pos2, variant.pos2_range);
        let info = format!(
            "{info};CIPOS={cipos_low},{cipos_high};CIPOS2={cipos2_low},{cipos2_high};DP={};AF={:.6}",
            variant.total_depth(),
            variant.total_allele_fraction()
        );
        write!(
            writer,
            "{}\t{}\tSTL_BND_{}\tN\t<BND>\t.\tPASS\t{}\tGT:DP:AD:AF:SR",
            variant.chrom1,
            variant.pos1,
            index + 1,
            info,
        )?;
        for sample in samples {
            let support = variant.support(sample);
            let depth = variant.depth(sample);
            // Depth is measured at the consensus breakpoint, which can sit a
            // base or two off an individual read's end, so clamp rather than
            // underflow.
            let reference = depth.saturating_sub(support);
            let genotype = if support > 0 { "0/1" } else { "0/0" };
            write!(
                writer,
                "\t{genotype}:{depth}:{reference},{support}:{:.6}:{support}",
                variant.allele_fraction(sample)
            )?;
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
        Junction {
            sample: sample.to_string(),
            read_name: read.to_string(),
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
        assert_eq!(variant.support_total, 3);
        assert_eq!(variant.support("s1"), 2);
        assert_eq!(variant.support("s2"), 1);
        assert_eq!(variant.gene1, "BCR");
        assert_eq!(variant.gene2, "ABL1");
        assert_eq!(variant.pos1, 101);
        assert_eq!(variant.pos2, 401);

        // The scatter of the member breakpoints is retained, so a generous
        // --sv-slop reveals how wide real clusters are.
        assert_eq!(variant.pos1_range, (100, 103));
        assert_eq!(variant.pos2_range, (400, 402));
        assert_eq!(range_offsets(variant.pos1, variant.pos1_range), (-1, 2));
        assert_eq!(range_offsets(variant.pos2, variant.pos2_range), (-1, 1));
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
        assert_eq!(variants[0].support_total, 1);
        assert_eq!(variants[0].support("s1"), 1);
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
