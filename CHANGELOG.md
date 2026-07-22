# Changelog

Notable changes to Stellerator. Versions follow [Semantic Versioning](https://semver.org/).

## 0.2.1

### Changed

- Refreshed semver-compatible dependencies (`anyhow`, `clap`, `rayon`, `serde`,
  `serde_json`, `thiserror`, `tracing-appender`, and transitives). No behaviour
  or API changes. The `noodles` crates were intentionally left at their pinned
  versions, since newer minor releases are breaking and warrant a dedicated
  migration.

## 0.2.0

Reorients the tool toward detecting low-frequency fusions in long-read
(ONT/PacBio) data. **Output formats changed in ways that break downstream
parsers** — see the migration notes below.

### Added

- Multiple BAMs in one run: `--bam` accepts several paths or a directory, each
  validated for a `.bai`/`.csi` index. All samples aggregate into one set of
  outputs, with a `sample` column and `sample=` FASTA field recording the source.
- Consensus structural-variant VCF via `--output-vcf`, clustering supporting
  reads into `BND` breakends annotated with genes, transcripts and breakpoint
  regions. `--sv-slop` sets the clustering tolerance (default 10 bp).
- Per-sample depth and allele fraction in the VCF: `FORMAT/GT:DP:AD:AF:SR` plus
  cohort-wide `INFO/DP` and `INFO/AF`, so support can be read against a real
  denominator instead of interpreted blind.
- `INFO/CIPOS` and `INFO/CIPOS2`, reporting how far supporting breakpoints
  scattered around the consensus position. Running with a generous `--sv-slop`
  makes the natural scatter measurable rather than guessed.
- `--min-mapq` (default `0`) to filter low-quality and multi-mapping
  alignments. At the default every alignment is kept and a warning is logged.
- `--min-depth` (default `0`) to drop consensus calls whose breakpoint is too
  thinly covered for the allele fraction to be meaningful. It thresholds on
  depth, never on support, so low-frequency events are not filtered away.
- `--include-duplicates` to opt back into reads flagged as PCR/optical
  duplicates, which are now skipped by default.

### Changed

- Every alignment in the `SA` tag is now used, not only the first, so a read
  split across several loci yields one candidate per supplementary alignment.
- Partner exon/intron and transcript labels use an SA-CIGAR-derived breakpoint
  position rather than the raw `SA` start.
- Reads flagged as PCR/optical duplicates no longer count as support.
- Default output paths derive from the input BAM and requested genes
  (`<bam-basename>.<genes>.<ext>`) instead of the fixed `stellerator.tsv` and
  `stellerator.fasta.gz`.
- The FASTA holds one record per read. A chimeric read supporting several
  junctions still produces one TSV row per junction, but its sequence is written
  once rather than duplicated per junction.

### Fixed

- Depth now counts supporting reads as spanning the breakpoint, which they do by
  definition. Previously a scattered cluster could report more support than
  depth, leaving `AD` not summing to `DP` and overstating `AF`.
- The depth pass opens one BAM reader per sample instead of one per call, and
  runs samples in parallel. Re-parsing the index for every call dominated
  runtime on long-read BAMs.

### Migration notes

- **TSV** gains a trailing `sample` column. Existing column positions are
  unchanged, so index-based parsers still work if they ignore trailing fields.
- **FASTA** headers gain `sample=`, and a chimeric read now yields one record
  rather than one per supplementary alignment.
- **Output filenames** change when `--output-tsv`/`--output-fasta` are omitted.
  Pass them explicitly to keep the previous fixed names.
- Duplicate-flagged reads are excluded by default; pass `--include-duplicates`
  to restore the previous behaviour.

## 0.1.2

Initial published release: single-BAM extraction of candidate fusion-supporting
reads with TSV and gzipped FASTA output, and exon/intron breakpoint labels from
the longest transcript per gene.
