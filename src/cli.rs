use std::path::PathBuf;

use clap::{ArgAction, Parser};

#[derive(Debug, Clone, Parser)]
#[command(
    name = "stellerator",
    version,
    about = "Extract candidate fusion-supporting reads for target genes from an indexed BAM"
)]
pub struct Args {
    #[arg(
        long,
        value_name = "BAM",
        required = true,
        num_args = 1..,
        help = "One or more indexed BAM files, or directories of BAMs; repeat the flag or pass multiple paths (e.g. --bam *.bam)"
    )]
    pub bam: Vec<PathBuf>,
    #[arg(long, value_name = "GFF_OR_GTF")]
    pub annotation: PathBuf,
    #[arg(
        long = "gene",
        value_name = "GENE",
        num_args = 1..,
        required_unless_present = "loci",
        conflicts_with = "loci",
        help = "Target gene to query; repeat for multiple genes. Mutually exclusive with --loci"
    )]
    pub genes: Vec<String>,
    #[arg(
        long,
        value_name = "GENE",
        help = "Restrict output to reads connecting the query gene(s) to this partner. Ignored in --loci mode, where each row carries its own partner"
    )]
    pub partner_gene: Option<String>,
    #[arg(
        long,
        value_name = "LOCI",
        help = "Batch file of `gene [partner] [tolerance]` rows, one job per line; each is processed with its own partner constraint and clustering tolerance. Mutually exclusive with --gene"
    )]
    pub loci: Option<PathBuf>,
    #[arg(
        long,
        value_name = "TSV",
        help = "TSV output path (default: <bam-basename>.<genes>.tsv)"
    )]
    pub output_tsv: Option<PathBuf>,
    #[arg(
        long,
        value_name = "FASTA_GZ",
        help = "Gzipped FASTA output path (default: <bam-basename>.<genes>.fasta.gz)"
    )]
    pub output_fasta: Option<PathBuf>,
    #[arg(
        long,
        value_name = "VCF",
        num_args = 0..=1,
        help = "VCF output of consensus structural variants. Pass a path, or give the flag alone for <bam-basename>.<genes>.vcf; omit to skip the VCF."
    )]
    pub output_vcf: Option<Option<PathBuf>>,
    #[arg(
        long,
        value_name = "BP",
        default_value_t = 200,
        help = "Breakpoint clustering tolerance in bp for consensus SV calling (also the per-row fallback in --loci mode)"
    )]
    pub sv_slop: usize,
    #[arg(
        long,
        action = ArgAction::SetTrue,
        help = "Include reads flagged as PCR/optical duplicates; they are skipped by default"
    )]
    pub include_duplicates: bool,
    #[arg(
        long,
        value_name = "MAPQ",
        default_value_t = 0,
        help = "Minimum mapping quality for a read to be considered; 0 (the default) takes every alignment"
    )]
    pub min_mapq: u8,
    #[arg(
        long,
        value_name = "N",
        default_value_t = 0,
        help = "Drop consensus VCF calls whose breakpoint has fewer than N spanning reads; 0 (the default) keeps every call"
    )]
    pub min_depth: usize,
    #[arg(long, value_name = "N", default_value_t = 1)]
    pub threads: usize,
    #[arg(long, action = ArgAction::SetTrue)]
    pub verbose: bool,
    #[arg(long, value_name = "LOG")]
    pub log_file: Option<PathBuf>,
}

impl Args {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
