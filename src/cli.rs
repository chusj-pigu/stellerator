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
    #[arg(long = "gene", value_name = "GENE", required = true, num_args = 1..)]
    pub genes: Vec<String>,
    #[arg(long, value_name = "GENE")]
    pub partner_gene: Option<String>,
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
        default_value_t = 10,
        help = "Breakpoint clustering tolerance in bp for consensus SV calling"
    )]
    pub sv_slop: usize,
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
