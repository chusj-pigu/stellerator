use super::{
    PartnerAlignment, ReadContext, bam_basename_token, build_hit, check_unique_sample_names,
    default_output_path, estimate_sa_breakpoint_position, filter_spans_by_gene_names,
    find_overlapping_partner, genes_basename_token, has_bam_extension, parse_sa_entries,
    resolve_bam_inputs, resolve_partner_filter, sample_name, sanitize_token,
};
use crate::annotation::{Exon, GeneSpan, Transcript};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

#[test]
fn chooses_other_overlapping_gene_at_breakpoint_or_start() {
    let query = GeneSpan {
        gene: "BCR".into(),
        reference_name: "chr22".into(),
        start: 100,
        end: 200,
        strand: Some('+'),
        transcripts: vec![],
    };
    let same_gene = GeneSpan {
        gene: "BCR".into(),
        reference_name: "chr9".into(),
        start: 400,
        end: 500,
        strand: Some('+'),
        transcripts: vec![],
    };
    let other_gene = GeneSpan {
        gene: "ABL1".into(),
        reference_name: "chr9".into(),
        start: 400,
        end: 500,
        strand: Some('+'),
        transcripts: vec![],
    };
    let spans = vec![same_gene, other_gene];
    let context = ReadContext {
        read_name: "read".into(),
        sequence: "ACGT".into(),
        read_flags: 0,
        query_strand: '+',
        reference_name: "chr22".into(),
        alignment_start: 100,
        alignment_end: 150,
        query_breakpoint_position: 150,
        cigar: "50M".into(),
        mapping_quality: Some(60),
        mate_reference_name: None,
        mate_alignment_start: None,
        sa_tag: "chr9,420,+,50M".into(),
    };
    for breakpoint in [450, 600] {
        let partner = PartnerAlignment {
            reference_name: "chr9".into(),
            start: 420,
            strand: '+',
            breakpoint,
        };
        let hit = build_hit("sample", &query, Some(&spans), true, &context, partner)
            .expect("eligible overlapping partner");
        assert_eq!(hit.tsv.matched_partner_gene.as_deref(), Some("ABL1"));
    }

    let mut overlapping = query.clone();
    overlapping.gene = "BCR-AS1".into();
    let same_locus = vec![query.clone(), overlapping];
    let partner = PartnerAlignment {
        reference_name: "chr22".into(),
        start: 120,
        strand: '+',
        breakpoint: 150,
    };
    let hit = build_hit("sample", &query, Some(&same_locus), true, &context, partner)
        .expect("overlapping partner at query locus");
    assert_eq!(hit.tsv.matched_partner_gene.as_deref(), Some("BCR-AS1"));
}

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
fn genes_token_joins_genes_and_appends_distinct_partner() {
    assert_eq!(genes_basename_token(&["BCR".to_string()], None), "BCR");
    assert_eq!(
        genes_basename_token(&["BCR".to_string(), "ABL1".to_string()], None),
        "BCR_ABL1"
    );
    assert_eq!(
        genes_basename_token(&["BCR".to_string()], Some("ABL1")),
        "BCR_ABL1"
    );
    // A partner already present as a query gene is not repeated.
    assert_eq!(
        genes_basename_token(&["BCR".to_string(), "ABL1".to_string()], Some("abl1")),
        "BCR_ABL1"
    );
}

#[test]
fn bam_token_uses_stem_then_shared_parent_directory() {
    // A single BAM contributes its own stem.
    assert_eq!(
        bam_basename_token(
            &[PathBuf::from("/data/sampleA.bam")],
            &["sampleA".to_string()]
        ),
        "sampleA"
    );

    // Several BAMs in one directory collapse to that directory name.
    assert_eq!(
        bam_basename_token(
            &[
                PathBuf::from("/data/cohort/a.bam"),
                PathBuf::from("/data/cohort/b.bam"),
            ],
            &["a".to_string(), "b".to_string()]
        ),
        "cohort"
    );

    // With no shared parent, fall back to the first sample name.
    assert_eq!(
        bam_basename_token(
            &[PathBuf::from("/x/a.bam"), PathBuf::from("/y/b.bam")],
            &["a".to_string(), "b".to_string()]
        ),
        "a"
    );
}

#[test]
fn default_output_path_combines_bam_and_gene_tokens() {
    assert_eq!(
        default_output_path("sampleA", "BCR_ABL1", "tsv"),
        PathBuf::from("sampleA.BCR_ABL1.tsv")
    );
    assert_eq!(
        default_output_path("cohort", "BCR", "fasta.gz"),
        PathBuf::from("cohort.BCR.fasta.gz")
    );
}

#[test]
fn sanitize_token_replaces_path_hostile_characters() {
    assert_eq!(sanitize_token("BCR/ABL1"), "BCR_ABL1");
    assert_eq!(sanitize_token("gene with space"), "gene_with_space");
    // Hyphens and dots are already filename-safe and are preserved.
    assert_eq!(sanitize_token("HLA-DRB1.v2"), "HLA-DRB1.v2");
}

#[test]
fn parses_every_sa_entry() {
    let partners = parse_sa_entries("chr9,420,-,20S80M,60,0;chr1,100,+,30M70S,55,1;");
    assert_eq!(partners.len(), 2);

    assert_eq!(partners[0].reference_name, "chr9");
    assert_eq!(partners[0].start, 420);
    assert_eq!(partners[0].strand, '-');
    // 20S80M is left-clipped, so the breakpoint sits at the aligned start.
    assert_eq!(partners[0].breakpoint, 420);

    assert_eq!(partners[1].reference_name, "chr1");
    // 30M70S is right-clipped, so the breakpoint sits at the aligned end.
    assert_eq!(partners[1].breakpoint, 129);
}

#[test]
fn sa_breakpoint_follows_clip_side() {
    assert_eq!(estimate_sa_breakpoint_position(500, "40S60M"), Some(500));
    assert_eq!(estimate_sa_breakpoint_position(500, "60M40S"), Some(559));
    // Reference-consuming ops (D/N) extend the aligned end.
    assert_eq!(
        estimate_sa_breakpoint_position(500, "10M5D10M20S"),
        Some(524)
    );
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
        find_overlapping_partner(&spans, "chr9", 350, "BCR").map(|span| span.gene.clone()),
        Some("ABL1".to_string())
    );
    assert_eq!(find_overlapping_partner(&spans, "chr9", 900, "BCR"), None);
    let mut overlapping = spans[0].clone();
    overlapping.gene = "BCR-AS1".to_string();
    let overlaps = vec![spans[0].clone(), overlapping];
    assert_eq!(
        find_overlapping_partner(&overlaps, "chr22", 150, "bcr")
            .unwrap()
            .gene,
        "BCR-AS1"
    );
}

#[test]
fn partner_filter_is_none_without_a_partner_gene() {
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

    // No partner gene => no filter, so the partner side is annotated against
    // every span.
    assert!(resolve_partner_filter(&spans, None).is_none());

    // A partner gene restricts the partner side to that gene's spans.
    let filtered = resolve_partner_filter(&spans, Some("ABL1")).expect("a filter");
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].gene, "ABL1");
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
