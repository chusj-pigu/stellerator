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
        Value::from("chr9,420,-,20S80M,60,0;"),
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
