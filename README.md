![Stellerator Logo]((https://docs.rs/crate/Stellerator/latest/source/images/logo.png))

# Stellerator

Stellerator is a Rust command-line tool for extracting candidate fusion-supporting reads from an indexed BAM file for one or more target genes. It scans the requested gene intervals from a GFF/GTF annotation, inspects supplementary alignments from the `SA` tag, and emits tabular and FASTA outputs for downstream review.

## What It Does

- Queries indexed BAM alignments across the requested gene interval(s)
- Uses supplementary alignments to identify candidate split-read evidence
- Writes a TSV summary of candidate reads and inferred partner loci
- Writes a gzipped FASTA file containing the supporting read sequences
- Annotates breakpoint regions as exon or intron labels using the longest transcript model per gene
- When `--partner-gene` is omitted, annotates supplementary loci against overlapping features from the annotation file when available

## Requirements

- Rust toolchain with Cargo
- Coordinate-sorted BAM with a sibling `.bai` or `.csi` index
- Gene annotation in GFF3 or GTF format

## Build

```bash
cargo build
```

## Run

Minimum example:

```bash
cargo run -- \
  --bam /path/to/sample.bam \
  --annotation /path/to/genes.gtf \
  --gene BCR
```

Example with multiple genes and a required partner gene:

```bash
cargo run -- \
  --bam /path/to/sample.bam \
  --annotation /path/to/genes.gtf \
  --gene BCR \
  --gene ABL1 \
  --partner-gene ABL1 \
  --output-tsv results/stellerator.tsv \
  --output-fasta results/stellerator.fasta.gz \
  --threads 4 \
  --verbose
```

## CLI Arguments

- `--bam`: input BAM file
- `--annotation`: input GFF3 or GTF file
- `--gene`: target gene to query; repeat for multiple genes
- `--partner-gene`: optional partner gene constraint
- `--output-tsv`: TSV output path
- `--output-fasta`: gzipped FASTA output path
- `--threads`: rayon worker count
- `--verbose`: enable debug logging
- `--log-file`: optional log file path

## Output

### TSV

The TSV includes:

- query gene name
- matched partner gene name, if annotated
- query transcript ID used for exon/intron labeling
- partner transcript ID used for exon/intron labeling
- breakpoint estimate in `query_region/partner_region` form
- read name, flags, coordinates, CIGAR, mapping quality, mate placement
- inferred partner reference, position, strand, and raw `SA` tag

### FASTA

The gzipped FASTA output contains the supporting read sequences. Each FASTA header includes the query gene, matched partner gene if available, transcript IDs used for labeling, breakpoint estimate, and inferred partner locus.

## Development

Run the standard checks from the repository root:

```bash
cargo test
cargo fmt -- --check
cargo clippy --all-targets --all-features -- -D warnings
```

## Tracked Work

This repository uses `bd` for task tracking. Per repository policy, open work is not maintained as a Markdown TODO list.

Use:

```bash
bd ready --json
bd list --json --status open
```

Current tracked themes include integration fixtures and expanded fusion heuristics / partner annotation behavior.
