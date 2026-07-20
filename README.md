![Stellerator Logo](https://docs.rs/crate/Stellerator/latest/source/images/logo.png)

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

Multiple samples in one run, via shell-expanded BAMs or a directory of BAMs:

```bash
# Several BAMs (the shell expands the glob into multiple --bam values)
cargo run -- \
  --bam cohort/*.bam \
  --annotation /path/to/genes.gtf \
  --gene BCR

# A directory; every *.bam inside (each with a .bai/.csi index) is scanned
cargo run -- \
  --bam cohort/ \
  --annotation /path/to/genes.gtf \
  --gene BCR
```

All samples are aggregated into the single `--output-tsv` and `--output-fasta`,
with a `sample` column (and `sample=` FASTA header field) recording the source
BAM. The sample name is the BAM file stem (`cohort/lib1.bam` becomes `lib1`); a
run aborts if two inputs would collapse to the same sample name.

Add `--output-vcf` to additionally emit consensus structural variants:

```bash
cargo run -- \
  --bam cohort/ \
  --annotation /path/to/genes.gtf \
  --gene BCR \
  --output-vcf results/stellerator.vcf \
  --sv-slop 10
```

## CLI Arguments

- `--bam`: one or more indexed BAM files, or directories of BAMs; repeat the flag or pass multiple paths (e.g. `--bam *.bam`)
- `--annotation`: input GFF3 or GTF file
- `--gene`: target gene to query; repeat for multiple genes
- `--partner-gene`: optional partner gene constraint
- `--output-tsv`: TSV output path (default: `<bam-basename>.<genes>.tsv`)
- `--output-fasta`: gzipped FASTA output path (default: `<bam-basename>.<genes>.fasta.gz`)
- `--output-vcf`: VCF output of consensus structural variants. Pass a path, or give the flag alone to use `<bam-basename>.<genes>.vcf`; omit the flag entirely to skip the VCF
- `--sv-slop`: breakpoint clustering tolerance in bp for consensus SV calling (default 10)
- `--threads`: rayon worker count
- `--verbose`: enable debug logging
- `--log-file`: optional log file path

### Default output names

When an output path is omitted, Stellerator builds one from the input BAM and
the requested genes, so parallel runs over different genes or samples do not
overwrite each other:

```
<bam-basename>.<genes>.<ext>
```

- `<bam-basename>` is the BAM file stem for a single input (`cohort/lib1.bam` gives
  `lib1`). For several BAMs it is their shared parent directory name (`--bam cohort/`
  or `--bam cohort/*.bam` gives `cohort`), falling back to the first sample name
  when the inputs have no common parent.
- `<genes>` joins the requested `--gene` values with `_`, appending `--partner-gene`
  when it is not already among them.

So `--bam sampleA.bam --gene BCR --partner-gene ABL1` writes
`sampleA.BCR_ABL1.tsv` and `sampleA.BCR_ABL1.fasta.gz`. Characters that are
awkward in filenames are replaced with `_`.

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
- `sample` name identifying the source BAM (final column)

### FASTA

The gzipped FASTA output contains the supporting read sequences. Each FASTA header includes the query gene, matched partner gene if available, transcript IDs used for labeling, breakpoint estimate, inferred partner locus, and the source `sample` name.

### VCF

When `--output-vcf` is given, supporting reads are clustered into consensus
structural variants and written as a multi-sample VCF (4.2). Every supplementary
(`SA`) alignment of each read is considered, and reads whose query and partner
breakpoints fall within `--sv-slop` bp (and share both chromosomes and strands)
are merged into one call. Each record is a `BND` breakend with:

- `CHROM`/`POS`: consensus query-side breakpoint (median of supporting reads)
- `INFO`: `SVTYPE=BND`, mate locus (`CHR2`/`POS2`), `STRANDS`, gene annotations
  (`GENE1`/`GENE2`), transcripts used for labeling, breakpoint `REGION` labels,
  and total supporting reads (`SR`)
- one genotype column per sample with the per-sample supporting-read count
  (`FORMAT/SR`)

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
