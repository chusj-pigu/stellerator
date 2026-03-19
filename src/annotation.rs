use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use anyhow::{Context, Result, bail};
use serde::Serialize;

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct Exon {
    pub start: i32,
    pub end: i32,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    pub exons: Vec<Exon>,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct BreakpointAnnotation {
    pub transcript_id: String,
    pub region: String,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct GeneSpan {
    pub gene: String,
    pub reference_name: String,
    pub start: i32,
    pub end: i32,
    pub strand: Option<char>,
    pub transcripts: Vec<Transcript>,
}

pub fn load_target_spans(path: &Path, requested_genes: &[String]) -> Result<Vec<GeneSpan>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open annotation file {}", path.display()))?;
    let reader = BufReader::new(file);
    let wanted: BTreeSet<String> = requested_genes
        .iter()
        .map(|gene| gene.to_ascii_lowercase())
        .collect();
    let filter_requested = !wanted.is_empty();

    let mut spans: BTreeMap<(String, String), GeneSpan> = BTreeMap::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line =
            line.with_context(|| format!("failed reading annotation line {}", line_number + 1))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            continue;
        }

        let reference_name = fields[0].trim();
        let feature_type = fields[2].trim().to_ascii_lowercase();
        let start = match fields[3].parse::<i32>() {
            Ok(value) => value,
            Err(_) => continue,
        };
        let end = match fields[4].parse::<i32>() {
            Ok(value) => value,
            Err(_) => continue,
        };

        if end < start {
            continue;
        }

        let attributes = parse_attributes(fields[8]);
        let Some(matched_gene) = resolve_gene_name(&attributes, &wanted, filter_requested) else {
            continue;
        };

        let strand = match fields[6].trim() {
            "+" => Some('+'),
            "-" => Some('-'),
            _ => None,
        };

        let key = (matched_gene.clone(), reference_name.to_owned());
        let gene_span = spans.entry(key).or_insert_with(|| GeneSpan {
            gene: matched_gene,
            reference_name: reference_name.to_owned(),
            start,
            end,
            strand,
            transcripts: Vec::new(),
        });

        gene_span.start = gene_span.start.min(start);
        gene_span.end = gene_span.end.max(end);
        if gene_span.strand.is_none() {
            gene_span.strand = strand;
        }

        if feature_type == "exon" {
            let transcript_id = resolve_transcript_id(&attributes)
                .unwrap_or_else(|| format!("{}_{}_{}", gene_span.gene, start, end));
            let transcript = gene_span
                .transcripts
                .iter_mut()
                .find(|transcript| transcript.id == transcript_id);

            if let Some(transcript) = transcript {
                transcript.exons.push(Exon { start, end });
            } else {
                gene_span.transcripts.push(Transcript {
                    id: transcript_id,
                    exons: vec![Exon { start, end }],
                });
            }
        }
    }

    if spans.is_empty() && filter_requested {
        bail!(
            "none of the requested genes were found in {}",
            path.display()
        );
    }

    let mut models: Vec<GeneSpan> = spans
        .into_values()
        .map(|mut span| {
            for transcript in &mut span.transcripts {
                transcript.exons.sort_by_key(|exon| exon.start);
                transcript
                    .exons
                    .dedup_by(|left, right| left.start == right.start && left.end == right.end);
            }
            span.transcripts
                .sort_by(|left, right| left.id.cmp(&right.id));
            span
        })
        .collect();

    models.sort_by(|left, right| {
        left.gene
            .cmp(&right.gene)
            .then(left.reference_name.cmp(&right.reference_name))
            .then(left.start.cmp(&right.start))
    });

    Ok(models)
}

pub fn breakpoint_annotation(span: &GeneSpan, position: usize) -> Option<BreakpointAnnotation> {
    let transcript = select_longest_transcript(span)?;
    let region = label_within_transcript(transcript, span.strand, position as i32)?;
    Some(BreakpointAnnotation {
        transcript_id: transcript.id.clone(),
        region,
    })
}

fn label_within_transcript(
    transcript: &Transcript,
    strand: Option<char>,
    position: i32,
) -> Option<String> {
    if transcript.exons.is_empty() {
        return None;
    }

    let ordered_exons: Vec<&Exon> = match strand {
        Some('-') => transcript.exons.iter().rev().collect(),
        _ => transcript.exons.iter().collect(),
    };

    for (index, exon) in ordered_exons.iter().enumerate() {
        if position >= exon.start && position <= exon.end {
            return Some(format!("exon{}", index + 1));
        }

        if let Some(next_exon) = ordered_exons.get(index + 1) {
            let intron_start = exon.end.min(next_exon.end);
            let intron_end = exon.start.max(next_exon.start);

            if position > intron_start && position < intron_end {
                return Some(format!("intron{}", index + 1));
            }
        }
    }

    None
}

fn select_longest_transcript(span: &GeneSpan) -> Option<&Transcript> {
    span.transcripts.iter().max_by(|left, right| {
        transcript_span_bases(left)
            .cmp(&transcript_span_bases(right))
            .then_with(|| right.id.cmp(&left.id))
    })
}

fn transcript_span_bases(transcript: &Transcript) -> i32 {
    transcript
        .exons
        .iter()
        .map(|exon| exon.end - exon.start + 1)
        .sum()
}

fn resolve_gene_name(
    attributes: &BTreeMap<String, String>,
    wanted: &BTreeSet<String>,
    filter_requested: bool,
) -> Option<String> {
    const KEYS: &[&str] = &[
        "gene_name",
        "gene",
        "gene_id",
        "Name",
        "ID",
        "transcript_id",
    ];

    KEYS.iter().find_map(|key| {
        attributes.get(*key).and_then(|value| {
            if !filter_requested || wanted.contains(&value.to_ascii_lowercase()) {
                Some(value.clone())
            } else {
                None
            }
        })
    })
}

fn resolve_transcript_id(attributes: &BTreeMap<String, String>) -> Option<String> {
    const KEYS: &[&str] = &["transcript_id", "transcript", "Parent", "ID"];
    KEYS.iter().find_map(|key| attributes.get(*key).cloned())
}

fn parse_attributes(raw: &str) -> BTreeMap<String, String> {
    if raw.contains('=') {
        parse_gff_attributes(raw)
    } else {
        parse_gtf_attributes(raw)
    }
}

fn parse_gff_attributes(raw: &str) -> BTreeMap<String, String> {
    raw.split(';')
        .filter_map(|entry| {
            let entry = entry.trim();
            if entry.is_empty() {
                return None;
            }

            let (key, value) = entry.split_once('=')?;
            Some((
                key.trim().to_string(),
                value.trim().trim_matches('"').to_string(),
            ))
        })
        .collect()
}

fn parse_gtf_attributes(raw: &str) -> BTreeMap<String, String> {
    raw.split(';')
        .filter_map(|entry| {
            let entry = entry.trim();
            if entry.is_empty() {
                return None;
            }

            let mut parts = entry.splitn(2, char::is_whitespace);
            let key = parts.next()?.trim();
            let value = parts.next()?.trim().trim_matches('"');
            Some((key.to_string(), value.to_string()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{
        BreakpointAnnotation, Exon, GeneSpan, Transcript, breakpoint_annotation, load_target_spans,
        parse_attributes, resolve_gene_name,
    };
    use std::{
        collections::{BTreeMap, BTreeSet},
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    #[test]
    fn parses_gff_attributes() {
        let attrs = parse_attributes("ID=gene-1;Name=BCR;gene_name=BCR");
        assert_eq!(attrs.get("Name"), Some(&"BCR".to_string()));
        assert_eq!(attrs.get("gene_name"), Some(&"BCR".to_string()));
    }

    #[test]
    fn parses_gtf_attributes() {
        let attrs = parse_attributes("gene_id \"ABL1\"; gene_name \"ABL1\";");
        assert_eq!(attrs.get("gene_id"), Some(&"ABL1".to_string()));
        assert_eq!(attrs.get("gene_name"), Some(&"ABL1".to_string()));
    }

    #[test]
    fn resolves_requested_gene_case_insensitively() {
        let mut attrs = BTreeMap::new();
        attrs.insert("gene_name".to_string(), "BCR".to_string());
        let wanted = BTreeSet::from(["bcr".to_string()]);
        assert_eq!(
            resolve_gene_name(&attrs, &wanted, true),
            Some("BCR".to_string())
        );
    }

    #[test]
    fn resolves_first_available_gene_when_not_filtering() {
        let mut attrs = BTreeMap::new();
        attrs.insert("gene_name".to_string(), "BCR".to_string());
        let wanted = BTreeSet::new();
        assert_eq!(
            resolve_gene_name(&attrs, &wanted, false),
            Some("BCR".to_string())
        );
    }

    #[test]
    fn labels_exons_and_introns_on_plus_strand() {
        let span = GeneSpan {
            gene: "BCR".to_string(),
            reference_name: "chr22".to_string(),
            start: 100,
            end: 300,
            strand: Some('+'),
            transcripts: vec![Transcript {
                id: "tx1".to_string(),
                exons: vec![
                    Exon {
                        start: 100,
                        end: 150,
                    },
                    Exon {
                        start: 200,
                        end: 250,
                    },
                ],
            }],
        };

        assert_eq!(
            breakpoint_annotation(&span, 120),
            Some(BreakpointAnnotation {
                transcript_id: "tx1".to_string(),
                region: "exon1".to_string(),
            })
        );
        assert_eq!(
            breakpoint_annotation(&span, 175),
            Some(BreakpointAnnotation {
                transcript_id: "tx1".to_string(),
                region: "intron1".to_string(),
            })
        );
        assert_eq!(
            breakpoint_annotation(&span, 225),
            Some(BreakpointAnnotation {
                transcript_id: "tx1".to_string(),
                region: "exon2".to_string(),
            })
        );
    }

    #[test]
    fn labels_exons_in_transcript_order_on_minus_strand() {
        let span = GeneSpan {
            gene: "ETV6".to_string(),
            reference_name: "chr12".to_string(),
            start: 100,
            end: 300,
            strand: Some('-'),
            transcripts: vec![Transcript {
                id: "tx1".to_string(),
                exons: vec![
                    Exon {
                        start: 100,
                        end: 150,
                    },
                    Exon {
                        start: 200,
                        end: 250,
                    },
                ],
            }],
        };

        assert_eq!(
            breakpoint_annotation(&span, 225),
            Some(BreakpointAnnotation {
                transcript_id: "tx1".to_string(),
                region: "exon1".to_string(),
            })
        );
        assert_eq!(
            breakpoint_annotation(&span, 175),
            Some(BreakpointAnnotation {
                transcript_id: "tx1".to_string(),
                region: "intron1".to_string(),
            })
        );
        assert_eq!(
            breakpoint_annotation(&span, 120),
            Some(BreakpointAnnotation {
                transcript_id: "tx1".to_string(),
                region: "exon2".to_string(),
            })
        );
    }

    #[test]
    fn uses_longest_transcript_for_breakpoint_annotation() {
        let span = GeneSpan {
            gene: "BCR".to_string(),
            reference_name: "chr22".to_string(),
            start: 100,
            end: 500,
            strand: Some('+'),
            transcripts: vec![
                Transcript {
                    id: "tx_short".to_string(),
                    exons: vec![
                        Exon {
                            start: 100,
                            end: 150,
                        },
                        Exon {
                            start: 300,
                            end: 350,
                        },
                    ],
                },
                Transcript {
                    id: "tx_long".to_string(),
                    exons: vec![
                        Exon {
                            start: 100,
                            end: 180,
                        },
                        Exon {
                            start: 260,
                            end: 340,
                        },
                        Exon {
                            start: 420,
                            end: 500,
                        },
                    ],
                },
            ],
        };

        assert_eq!(
            breakpoint_annotation(&span, 170),
            Some(BreakpointAnnotation {
                transcript_id: "tx_long".to_string(),
                region: "exon1".to_string(),
            })
        );
    }

    #[test]
    fn breaks_ties_by_transcript_id() {
        let span = GeneSpan {
            gene: "BCR".to_string(),
            reference_name: "chr22".to_string(),
            start: 100,
            end: 300,
            strand: Some('+'),
            transcripts: vec![
                Transcript {
                    id: "tx_b".to_string(),
                    exons: vec![Exon {
                        start: 100,
                        end: 150,
                    }],
                },
                Transcript {
                    id: "tx_a".to_string(),
                    exons: vec![Exon {
                        start: 100,
                        end: 150,
                    }],
                },
            ],
        };

        assert_eq!(
            breakpoint_annotation(&span, 120),
            Some(BreakpointAnnotation {
                transcript_id: "tx_a".to_string(),
                region: "exon1".to_string(),
            })
        );
    }

    #[test]
    fn loads_all_gene_spans_when_requested_gene_filter_is_empty() {
        let path = unique_test_path("all_genes.gtf");
        fs::write(
            &path,
            "\
chr1\tsrc\texon\t100\t150\t.\t+\t.\tgene_name \"BCR\"; transcript_id \"tx1\";\n\
chr1\tsrc\texon\t200\t250\t.\t+\t.\tgene_name \"ABL1\"; transcript_id \"tx2\";\n",
        )
        .unwrap();

        let spans = load_target_spans(&path, &[]).unwrap();
        let genes: Vec<_> = spans.iter().map(|span| span.gene.as_str()).collect();
        assert_eq!(genes, vec!["ABL1", "BCR"]);

        fs::remove_file(path).unwrap();
    }

    fn unique_test_path(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("stellerator_{nanos}_{name}"))
    }
}
