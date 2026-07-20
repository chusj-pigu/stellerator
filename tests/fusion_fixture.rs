use std::{
    fs::{self, File},
    io::Read,
    path::{Path, PathBuf},
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

use flate2::read::GzDecoder;
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::{Cigar, Data, QualityScores, Sequence, data::field::Value},
    },
};

#[test]
fn extracts_fusion_read_and_infers_partner_annotation() -> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("fusion.bam");
    let annotation_path = fixture_dir.join("genes.gff3");
    let output_tsv = fixture_dir.join("stellerator.tsv");
    let output_fasta = fixture_dir.join("stellerator.fasta.gz");

    write_annotation(&annotation_path)?;
    write_indexed_bam(&bam_path, "fusion-read-1")?;

    let output = Command::new(env!("CARGO_BIN_EXE_stellerator"))
        .arg("--bam")
        .arg(&bam_path)
        .arg("--annotation")
        .arg(&annotation_path)
        .arg("--gene")
        .arg("BCR")
        .arg("--output-tsv")
        .arg(&output_tsv)
        .arg("--output-fasta")
        .arg(&output_fasta)
        .arg("--threads")
        .arg("1")
        .output()?;

    assert!(
        output.status.success(),
        "stellerator failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tsv = fs::read_to_string(&output_tsv)?;
    let rows: Vec<_> = tsv.lines().collect();
    assert_eq!(rows.len(), 2, "{tsv}");

    let values: Vec<_> = rows[1].split('\t').collect();
    assert_eq!(values[0], "BCR");
    assert_eq!(values[1], "ABL1");
    assert_eq!(values[2], "txBCR");
    assert_eq!(values[3], "txABL1");
    assert_eq!(values[4], "exon1/exon1");
    assert_eq!(values[5], "fusion-read-1");
    assert_eq!(values[7], "chr22");
    assert_eq!(values[14], "chr9");
    assert_eq!(values[15], "420");
    assert_eq!(values[16], "-");
    assert_eq!(values[18], "fusion");

    let fasta = read_gzip_to_string(&output_fasta)?;
    assert!(fasta.contains(">fusion-read-1 gene=BCR matched_partner_gene=ABL1"));
    assert!(fasta.contains("partner_transcript_id=txABL1"));
    assert!(fasta.contains("breakpoint_estimate=exon1/exon1"));
    assert!(fasta.contains("partner=chr9:420 strand=-"));
    assert!(fasta.contains("sample=fusion"));
    assert!(fasta.contains("ACGTACGTACGT"));

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn aggregates_multiple_bams_from_directory_with_sample_provenance()
-> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_dir = fixture_dir.join("bams");
    fs::create_dir(&bam_dir)?;
    let annotation_path = fixture_dir.join("genes.gff3");
    let output_tsv = fixture_dir.join("stellerator.tsv");
    let output_fasta = fixture_dir.join("stellerator.fasta.gz");

    write_annotation(&annotation_path)?;
    write_indexed_bam(&bam_dir.join("sampleA.bam"), "read-A")?;
    write_indexed_bam(&bam_dir.join("sampleB.bam"), "read-B")?;

    let output = Command::new(env!("CARGO_BIN_EXE_stellerator"))
        .arg("--bam")
        .arg(&bam_dir)
        .arg("--annotation")
        .arg(&annotation_path)
        .arg("--gene")
        .arg("BCR")
        .arg("--output-tsv")
        .arg(&output_tsv)
        .arg("--output-fasta")
        .arg(&output_fasta)
        .arg("--threads")
        .arg("1")
        .output()?;

    assert!(
        output.status.success(),
        "stellerator failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let tsv = fs::read_to_string(&output_tsv)?;
    let rows: Vec<_> = tsv.lines().collect();
    assert_eq!(
        rows.len(),
        3,
        "expected header plus one row per sample\n{tsv}"
    );

    let mut samples: Vec<&str> = rows[1..]
        .iter()
        .map(|row| row.split('\t').next_back().unwrap())
        .collect();
    samples.sort_unstable();
    assert_eq!(samples, vec!["sampleA", "sampleB"]);

    let fasta = read_gzip_to_string(&output_fasta)?;
    assert!(fasta.contains(">read-A gene=BCR"));
    assert!(fasta.contains(">read-B gene=BCR"));
    assert!(fasta.contains("sample=sampleA"));
    assert!(fasta.contains("sample=sampleB"));

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn emits_consensus_sv_vcf_with_gene_annotation_and_support()
-> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("fusion.bam");
    let annotation_path = fixture_dir.join("genes.gff3");
    let output_tsv = fixture_dir.join("stellerator.tsv");
    let output_fasta = fixture_dir.join("stellerator.fasta.gz");
    let output_vcf = fixture_dir.join("stellerator.vcf");

    write_annotation(&annotation_path)?;
    write_indexed_bam(&bam_path, "fusion-read-1")?;

    let output = Command::new(env!("CARGO_BIN_EXE_stellerator"))
        .arg("--bam")
        .arg(&bam_path)
        .arg("--annotation")
        .arg(&annotation_path)
        .arg("--gene")
        .arg("BCR")
        .arg("--output-tsv")
        .arg(&output_tsv)
        .arg("--output-fasta")
        .arg(&output_fasta)
        .arg("--output-vcf")
        .arg(&output_vcf)
        .arg("--threads")
        .arg("1")
        .output()?;

    assert!(
        output.status.success(),
        "stellerator failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    let vcf = fs::read_to_string(&output_vcf)?;
    assert!(vcf.contains("##fileformat=VCFv4.2"));
    assert!(vcf.contains("##contig=<ID=chr22,length=1000>"));

    let header_line = vcf
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("VCF column header");
    assert!(header_line.ends_with("\tfusion"), "{header_line}");

    let record = vcf
        .lines()
        .find(|line| !line.starts_with('#'))
        .expect("a VCF record");
    let cols: Vec<_> = record.split('\t').collect();
    assert_eq!(cols[0], "chr22"); // query-side chromosome
    assert_eq!(cols[3], "N"); // REF
    assert_eq!(cols[4], "<BND>"); // ALT
    assert_eq!(cols[6], "PASS"); // FILTER
    assert!(cols[7].contains("SVTYPE=BND"));
    assert!(cols[7].contains("CHR2=chr9"));
    assert!(cols[7].contains("POS2=420"));
    assert!(cols[7].contains("GENE1=BCR"));
    assert!(cols[7].contains("GENE2=ABL1"));
    assert!(cols[7].contains("SR=1"));
    assert!(cols[7].contains("DP=1"));
    assert_eq!(cols[8], "GT:DP:AD:AF:SR"); // FORMAT
    // The lone spanning read is the supporting read, so AF is 1.
    assert_eq!(cols[9], "0/1:1:0,1:1.000000:1");

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn derives_default_output_paths_from_bam_and_genes() -> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("fusion.bam");
    let annotation_path = fixture_dir.join("genes.gff3");

    write_annotation(&annotation_path)?;
    write_indexed_bam(&bam_path, "fusion-read-1")?;

    // No --output-* flags at all; run inside the fixture dir so the relative
    // default paths land there.
    let output = Command::new(env!("CARGO_BIN_EXE_stellerator"))
        .current_dir(&fixture_dir)
        .arg("--bam")
        .arg(&bam_path)
        .arg("--annotation")
        .arg(&annotation_path)
        .arg("--gene")
        .arg("BCR")
        .arg("--threads")
        .arg("1")
        .output()?;

    assert!(
        output.status.success(),
        "stellerator failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    // Derived from the BAM stem (`fusion`) and the requested gene (`BCR`).
    assert!(
        fixture_dir.join("fusion.BCR.tsv").exists(),
        "expected derived TSV path"
    );
    assert!(
        fixture_dir.join("fusion.BCR.fasta.gz").exists(),
        "expected derived FASTA path"
    );
    // The VCF stays opt-in, so it must not be written without the flag.
    assert!(
        !fixture_dir.join("fusion.BCR.vcf").exists(),
        "VCF should not be written unless --output-vcf is given"
    );

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn bare_output_vcf_flag_uses_derived_default_path() -> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("fusion.bam");
    let annotation_path = fixture_dir.join("genes.gff3");

    write_annotation(&annotation_path)?;
    write_indexed_bam(&bam_path, "fusion-read-1")?;

    let output = Command::new(env!("CARGO_BIN_EXE_stellerator"))
        .current_dir(&fixture_dir)
        .arg("--bam")
        .arg(&bam_path)
        .arg("--annotation")
        .arg(&annotation_path)
        .arg("--gene")
        .arg("BCR")
        .arg("--partner-gene")
        .arg("ABL1")
        // Flag given with no value: the path is derived.
        .arg("--output-vcf")
        .arg("--threads")
        .arg("1")
        .output()?;

    assert!(
        output.status.success(),
        "stellerator failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    // Partner gene joins the gene token, so outputs are `fusion.BCR_ABL1.*`.
    assert!(
        fixture_dir.join("fusion.BCR_ABL1.vcf").exists(),
        "expected derived VCF path"
    );
    assert!(fixture_dir.join("fusion.BCR_ABL1.tsv").exists());
    assert!(fixture_dir.join("fusion.BCR_ABL1.fasta.gz").exists());

    let vcf = fs::read_to_string(fixture_dir.join("fusion.BCR_ABL1.vcf"))?;
    assert!(vcf.contains("##fileformat=VCFv4.2"));
    assert!(vcf.contains("GENE1=BCR"));
    assert!(vcf.contains("GENE2=ABL1"));

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn skips_duplicate_flagged_reads_unless_included() -> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("dups.bam");
    let annotation_path = fixture_dir.join("genes.gff3");

    write_annotation(&annotation_path)?;
    write_indexed_bam_reads(
        &bam_path,
        &[
            ("keep-read", Flags::empty(), 60),
            ("dup-read", Flags::DUPLICATE, 60),
        ],
    )?;

    // By default the duplicate-flagged read must not contribute support.
    let default_tsv = fixture_dir.join("default.tsv");
    run_extract(
        &bam_path,
        &annotation_path,
        &default_tsv,
        &fixture_dir.join("default.fasta.gz"),
        &[],
    )?;
    let tsv = fs::read_to_string(&default_tsv)?;
    assert_eq!(tsv.lines().count(), 2, "expected header + 1 row\n{tsv}");
    assert!(tsv.contains("keep-read"), "{tsv}");
    assert!(!tsv.contains("dup-read"), "{tsv}");

    // Opting in restores it.
    let included_tsv = fixture_dir.join("included.tsv");
    run_extract(
        &bam_path,
        &annotation_path,
        &included_tsv,
        &fixture_dir.join("included.fasta.gz"),
        &["--include-duplicates"],
    )?;
    let tsv = fs::read_to_string(&included_tsv)?;
    assert_eq!(tsv.lines().count(), 3, "expected header + 2 rows\n{tsv}");
    assert!(tsv.contains("keep-read"), "{tsv}");
    assert!(tsv.contains("dup-read"), "{tsv}");

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn min_mapq_defaults_to_taking_everything_and_warns() -> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("mapq.bam");
    let annotation_path = fixture_dir.join("genes.gff3");

    write_annotation(&annotation_path)?;
    write_indexed_bam_reads(
        &bam_path,
        &[
            ("high-mapq", Flags::empty(), 60),
            ("low-mapq", Flags::empty(), 5),
        ],
    )?;

    // Default takes everything, including the poorly mapped read.
    let default_tsv = fixture_dir.join("default.tsv");
    let output = Command::new(env!("CARGO_BIN_EXE_stellerator"))
        .arg("--bam")
        .arg(&bam_path)
        .arg("--annotation")
        .arg(&annotation_path)
        .arg("--gene")
        .arg("BCR")
        .arg("--output-tsv")
        .arg(&default_tsv)
        .arg("--output-fasta")
        .arg(fixture_dir.join("default.fasta.gz"))
        .arg("--threads")
        .arg("1")
        .output()?;
    assert!(output.status.success());

    let tsv = fs::read_to_string(&default_tsv)?;
    assert_eq!(tsv.lines().count(), 3, "expected header + 2 rows\n{tsv}");
    assert!(tsv.contains("high-mapq"), "{tsv}");
    assert!(tsv.contains("low-mapq"), "{tsv}");

    // ...and says so, so an unfiltered run is never silent.
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--min-mapq is 0"),
        "expected a warning about taking every alignment\n{stderr}"
    );

    // Raising the floor drops the low-quality alignment.
    let filtered_tsv = fixture_dir.join("filtered.tsv");
    run_extract(
        &bam_path,
        &annotation_path,
        &filtered_tsv,
        &fixture_dir.join("filtered.fasta.gz"),
        &["--min-mapq", "20"],
    )?;
    let tsv = fs::read_to_string(&filtered_tsv)?;
    assert_eq!(tsv.lines().count(), 2, "expected header + 1 row\n{tsv}");
    assert!(tsv.contains("high-mapq"), "{tsv}");
    assert!(!tsv.contains("low-mapq"), "{tsv}");

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn reports_low_allele_fraction_against_spanning_depth() -> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("lowfreq.bam");
    let annotation_path = fixture_dir.join("genes.gff3");
    let output_vcf = fixture_dir.join("lowfreq.vcf");

    write_annotation(&annotation_path)?;
    // One supporting read against three non-supporting spanning reads.
    write_indexed_bam_with_background(&bam_path, 3)?;

    run_extract(
        &bam_path,
        &annotation_path,
        &fixture_dir.join("lowfreq.tsv"),
        &fixture_dir.join("lowfreq.fasta.gz"),
        &["--output-vcf", output_vcf.to_str().unwrap()],
    )?;

    let vcf = fs::read_to_string(&output_vcf)?;
    let record = vcf
        .lines()
        .find(|line| !line.starts_with('#'))
        .expect("a VCF record");
    let cols: Vec<&str> = record.split('\t').collect();

    // Only the split read supports the junction, but four reads span it.
    assert!(cols[7].contains("SR=1"), "{record}");
    assert!(cols[7].contains("DP=4"), "{record}");
    assert!(cols[7].contains("AF=0.250000"), "{record}");
    assert_eq!(cols[8], "GT:DP:AD:AF:SR");
    assert_eq!(cols[9], "0/1:4:3,1:0.250000:1", "{record}");

    // The background reads must not be reported as candidates themselves.
    let tsv = fs::read_to_string(fixture_dir.join("lowfreq.tsv"))?;
    assert_eq!(tsv.lines().count(), 2, "expected header + 1 row\n{tsv}");
    assert!(!tsv.contains("background-"), "{tsv}");

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn scattered_cluster_reports_consistent_depth_and_scatter() -> Result<(), Box<dyn std::error::Error>>
{
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("scatter.bam");
    let annotation_path = fixture_dir.join("genes.gff3");
    let output_vcf = fixture_dir.join("scatter.vcf");

    write_annotation(&annotation_path)?;
    write_indexed_scattered_bam(&bam_path)?;

    // A generous tolerance keeps the three reads in one cluster.
    run_extract(
        &bam_path,
        &annotation_path,
        &fixture_dir.join("scatter.tsv"),
        &fixture_dir.join("scatter.fasta.gz"),
        &[
            "--output-vcf",
            output_vcf.to_str().unwrap(),
            "--sv-slop",
            "100",
        ],
    )?;

    let vcf = fs::read_to_string(&output_vcf)?;
    let records: Vec<&str> = vcf.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(records.len(), 1, "expected a single merged call\n{vcf}");

    let cols: Vec<&str> = records[0].split('\t').collect();
    let info = cols[7];

    // The observed scatter is reported, not hidden by the consensus position.
    assert!(info.contains("CIPOS=-5,6"), "{info}");
    assert!(info.contains("CIPOS2=-1,2"), "{info}");

    // read-1 aligns 120-139 and so does not reach the consensus breakpoint at
    // 144, but it supports the junction and must still count toward depth.
    assert!(info.contains("SR=3"), "{info}");
    assert!(info.contains("DP=3"), "{info}");
    assert_eq!(cols[8], "GT:DP:AD:AF:SR");
    assert_eq!(cols[9], "0/1:3:0,3:1.000000:3", "{records:?}");

    // AD must sum to DP.
    let sample: Vec<&str> = cols[9].split(':').collect();
    let depth: usize = sample[1].parse()?;
    let allele_depths: Vec<usize> = sample[2]
        .split(',')
        .map(|value| value.parse().unwrap())
        .collect();
    assert_eq!(
        allele_depths.iter().sum::<usize>(),
        depth,
        "AD must sum to DP"
    );

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn writes_one_fasta_record_per_read_across_multiple_sa_entries()
-> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("chimeric.bam");
    let annotation_path = fixture_dir.join("genes.gff3");
    let output_tsv = fixture_dir.join("chimeric.tsv");
    let output_fasta = fixture_dir.join("chimeric.fasta.gz");

    write_annotation(&annotation_path)?;
    write_indexed_multi_sa_bam(&bam_path, "chimeric-read")?;

    run_extract(&bam_path, &annotation_path, &output_tsv, &output_fasta, &[])?;

    // One row per supplementary alignment: the junction detail is preserved.
    let tsv = fs::read_to_string(&output_tsv)?;
    assert_eq!(tsv.lines().count(), 3, "expected header + 2 rows\n{tsv}");
    assert_eq!(tsv.matches("chimeric-read").count(), 2, "{tsv}");

    // ...but the sequence is written once, not once per junction.
    let fasta = read_gzip_to_string(&output_fasta)?;
    assert_eq!(
        fasta.matches('>').count(),
        1,
        "expected a single FASTA record\n{fasta}"
    );
    // Exactly one sequence line, so the read is not emitted once per junction.
    let sequence_lines = fasta
        .lines()
        .filter(|line| !line.starts_with('>') && !line.is_empty())
        .count();
    assert_eq!(
        sequence_lines, 1,
        "read sequence must not be duplicated\n{fasta}"
    );

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

#[test]
fn min_depth_drops_thinly_covered_calls() -> Result<(), Box<dyn std::error::Error>> {
    let fixture_dir = unique_fixture_dir()?;
    let bam_path = fixture_dir.join("lowdepth.bam");
    let annotation_path = fixture_dir.join("genes.gff3");

    write_annotation(&annotation_path)?;
    // One supporting read plus three spanning reads gives the call DP=4.
    write_indexed_bam_with_background(&bam_path, 3)?;

    let record_count =
        |vcf: &str| -> usize { vcf.lines().filter(|line| !line.starts_with('#')).count() };

    // At the depth of the call, it survives.
    let kept_vcf = fixture_dir.join("kept.vcf");
    run_extract(
        &bam_path,
        &annotation_path,
        &fixture_dir.join("kept.tsv"),
        &fixture_dir.join("kept.fasta.gz"),
        &[
            "--output-vcf",
            kept_vcf.to_str().unwrap(),
            "--min-depth",
            "4",
        ],
    )?;
    let vcf = fs::read_to_string(&kept_vcf)?;
    assert_eq!(record_count(&vcf), 1, "expected the call to survive\n{vcf}");

    // One above it, the call is dropped but the VCF is still well formed.
    let dropped_vcf = fixture_dir.join("dropped.vcf");
    run_extract(
        &bam_path,
        &annotation_path,
        &fixture_dir.join("dropped.tsv"),
        &fixture_dir.join("dropped.fasta.gz"),
        &[
            "--output-vcf",
            dropped_vcf.to_str().unwrap(),
            "--min-depth",
            "5",
        ],
    )?;
    let vcf = fs::read_to_string(&dropped_vcf)?;
    assert_eq!(
        record_count(&vcf),
        0,
        "expected the call to be dropped\n{vcf}"
    );
    assert!(vcf.contains("##fileformat=VCFv4.2"), "{vcf}");
    assert!(vcf.contains("#CHROM"), "{vcf}");

    // Filtering the VCF must not touch the per-read candidate outputs.
    let tsv = fs::read_to_string(fixture_dir.join("dropped.tsv"))?;
    assert_eq!(tsv.lines().count(), 2, "expected header + 1 row\n{tsv}");

    fs::remove_dir_all(fixture_dir)?;
    Ok(())
}

fn write_annotation(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    fs::write(
        path,
        "\
##gff-version 3
chr22\tstellerator\tgene\t100\t250\t.\t+\t.\tID=gene-BCR;Name=BCR
chr22\tstellerator\tmRNA\t100\t250\t.\t+\t.\tID=txBCR;Parent=gene-BCR
chr22\tstellerator\texon\t100\t150\t.\t+\t.\tID=exon-BCR-1;Parent=txBCR
chr22\tstellerator\texon\t200\t250\t.\t+\t.\tID=exon-BCR-2;Parent=txBCR
chr9\tstellerator\tgene\t300\t450\t.\t-\t.\tID=gene-ABL1;Name=ABL1
chr9\tstellerator\tmRNA\t300\t450\t.\t-\t.\tID=txABL1;Parent=gene-ABL1
chr9\tstellerator\texon\t300\t350\t.\t-\t.\tID=exon-ABL1-1;Parent=txABL1
chr9\tstellerator\texon\t400\t450\t.\t-\t.\tID=exon-ABL1-2;Parent=txABL1
",
    )?;
    Ok(())
}

fn write_indexed_bam(path: &Path, read_name: &str) -> Result<(), Box<dyn std::error::Error>> {
    write_indexed_bam_reads(path, &[(read_name, Flags::empty(), 60)])
}

/// Write an indexed BAM containing one split-read record per entry, each with
/// the given name and flags.
fn write_indexed_bam_reads(
    path: &Path,
    reads: &[(&str, Flags, u8)],
) -> Result<(), Box<dyn std::error::Error>> {
    let header: sam::Header = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr9\tLN:1000
@SQ\tSN:chr22\tLN:1000
"
    .parse()?;

    let sequence = "ACGT".repeat(25);
    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header)?;

    for (read_name, flags, mapq) in reads {
        let cigar: Cigar = [Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, 80)]
            .into_iter()
            .collect();
        let data: Data = [(
            Tag::OTHER_ALIGNMENTS,
            Value::from("chr9,420,-,20S80M,60,0;"),
        )]
        .into_iter()
        .collect();

        let record = RecordBuf::builder()
            .set_name(*read_name)
            .set_flags(*flags)
            .set_reference_sequence_id(1)
            .set_alignment_start(Position::try_from(120)?)
            .set_mapping_quality(MappingQuality::new(*mapq).expect("valid MAPQ"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(sequence.as_bytes()))
            .set_quality_scores(QualityScores::from(vec![30; sequence.len()]))
            .set_data(data)
            .build();

        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish()?;

    let index = bam::fs::index(path)?;
    bam::bai::fs::write(path.with_extension("bam.bai"), &index)?;

    Ok(())
}

/// Write an indexed BAM with one split read plus `background` plain reads that
/// span the same breakpoint without supporting any fusion, so the supporting
/// read sits at a known low allele fraction.
fn write_indexed_bam_with_background(
    path: &Path,
    background: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let header: sam::Header = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr9\tLN:1000
@SQ\tSN:chr22\tLN:1000
"
    .parse()?;

    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header)?;

    // Background first to keep the file coordinate sorted: 100M from 100 spans
    // the breakpoint at 139 without carrying an SA tag.
    let background_sequence = "ACGT".repeat(25);
    for index in 0..background {
        let cigar: Cigar = [Op::new(Kind::Match, 100)].into_iter().collect();
        let record = RecordBuf::builder()
            .set_name(format!("background-{index}"))
            .set_flags(Flags::empty())
            .set_reference_sequence_id(1)
            .set_alignment_start(Position::try_from(100)?)
            .set_mapping_quality(MappingQuality::new(60).expect("valid MAPQ"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(background_sequence.as_bytes()))
            .set_quality_scores(QualityScores::from(vec![30; background_sequence.len()]))
            .build();
        writer.write_alignment_record(&header, &record)?;
    }

    let sequence = "ACGT".repeat(25);
    let cigar: Cigar = [Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, 80)]
        .into_iter()
        .collect();
    let data: Data = [(
        Tag::OTHER_ALIGNMENTS,
        Value::from("chr9,420,-,20S80M,60,0;"),
    )]
    .into_iter()
    .collect();
    let record = RecordBuf::builder()
        .set_name("fusion-read-1")
        .set_flags(Flags::empty())
        .set_reference_sequence_id(1)
        .set_alignment_start(Position::try_from(120)?)
        .set_mapping_quality(MappingQuality::new(60).expect("valid MAPQ"))
        .set_cigar(cigar)
        .set_sequence(Sequence::from(sequence.as_bytes()))
        .set_quality_scores(QualityScores::from(vec![30; sequence.len()]))
        .set_data(data)
        .build();
    writer.write_alignment_record(&header, &record)?;

    writer.try_finish()?;
    let index = bam::fs::index(path)?;
    bam::bai::fs::write(path.with_extension("bam.bai"), &index)?;

    Ok(())
}

/// Write an indexed BAM where three reads support one junction but place it a
/// few bases apart, as aligner wobble does. The earliest read's alignment ends
/// before the consensus breakpoint, so it only counts toward depth if
/// supporting reads are treated as spanning.
fn write_indexed_scattered_bam(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let header: sam::Header = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr9\tLN:1000
@SQ\tSN:chr22\tLN:1000
"
    .parse()?;

    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header)?;

    let sequence = "ACGT".repeat(25);
    for (name, start, partner) in [
        ("read-1", 120usize, "chr9,420,-,20S80M,60,0;"),
        ("read-2", 125, "chr9,422,-,20S80M,60,0;"),
        ("read-3", 131, "chr9,419,-,20S80M,60,0;"),
    ] {
        let cigar: Cigar = [Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, 80)]
            .into_iter()
            .collect();
        let data: Data = [(Tag::OTHER_ALIGNMENTS, Value::from(partner))]
            .into_iter()
            .collect();
        let record = RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::empty())
            .set_reference_sequence_id(1)
            .set_alignment_start(Position::try_from(start)?)
            .set_mapping_quality(MappingQuality::new(60).expect("valid MAPQ"))
            .set_cigar(cigar)
            .set_sequence(Sequence::from(sequence.as_bytes()))
            .set_quality_scores(QualityScores::from(vec![30; sequence.len()]))
            .set_data(data)
            .build();
        writer.write_alignment_record(&header, &record)?;
    }

    writer.try_finish()?;
    let index = bam::fs::index(path)?;
    bam::bai::fs::write(path.with_extension("bam.bai"), &index)?;

    Ok(())
}

/// Write an indexed BAM holding one read whose `SA` tag lists two supplementary
/// alignments, as a chimeric long read would.
fn write_indexed_multi_sa_bam(
    path: &Path,
    read_name: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let header: sam::Header = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr9\tLN:1000
@SQ\tSN:chr22\tLN:1000
"
    .parse()?;

    let sequence = "ACGT".repeat(25);
    let cigar: Cigar = [Op::new(Kind::Match, 20), Op::new(Kind::SoftClip, 80)]
        .into_iter()
        .collect();
    let data: Data = [(
        Tag::OTHER_ALIGNMENTS,
        Value::from("chr9,420,-,20S80M,60,0;chr9,800,+,30S70M,60,0;"),
    )]
    .into_iter()
    .collect();

    let record = RecordBuf::builder()
        .set_name(read_name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(1)
        .set_alignment_start(Position::try_from(120)?)
        .set_mapping_quality(MappingQuality::new(60).expect("valid MAPQ"))
        .set_cigar(cigar)
        .set_sequence(Sequence::from(sequence.as_bytes()))
        .set_quality_scores(QualityScores::from(vec![30; sequence.len()]))
        .set_data(data)
        .build();

    let file = File::create(path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header)?;
    writer.write_alignment_record(&header, &record)?;
    writer.try_finish()?;

    let index = bam::fs::index(path)?;
    bam::bai::fs::write(path.with_extension("bam.bai"), &index)?;

    Ok(())
}

/// Run the CLI against a fixture, asserting it succeeded.
fn run_extract(
    bam: &Path,
    annotation: &Path,
    output_tsv: &Path,
    output_fasta: &Path,
    extra: &[&str],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_stellerator"));
    command
        .arg("--bam")
        .arg(bam)
        .arg("--annotation")
        .arg(annotation)
        .arg("--gene")
        .arg("BCR")
        .arg("--output-tsv")
        .arg(output_tsv)
        .arg("--output-fasta")
        .arg(output_fasta)
        .arg("--threads")
        .arg("1");
    for arg in extra {
        command.arg(arg);
    }

    let output = command.output()?;
    assert!(
        output.status.success(),
        "stellerator failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );

    Ok(())
}

fn read_gzip_to_string(path: &Path) -> Result<String, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut decoder = GzDecoder::new(file);
    let mut contents = String::new();
    decoder.read_to_string(&mut contents)?;
    Ok(contents)
}

fn unique_fixture_dir() -> Result<PathBuf, Box<dyn std::error::Error>> {
    let nanos = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let path = std::env::temp_dir().join(format!("stellerator_fixture_{nanos}"));
    fs::create_dir(&path)?;
    Ok(path)
}
