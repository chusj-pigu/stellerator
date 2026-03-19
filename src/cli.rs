use std::path::PathBuf;

use clap::{ArgAction, Parser};

#[derive(Debug, Clone, Parser)]
#[command(
    name = "stellerator",
    version,
    about = "Extract candidate fusion-supporting reads for target genes from an indexed BAM"
)]
pub struct Args {
    #[arg(long, value_name = "BAM")]
    pub bam: PathBuf,
    #[arg(long, value_name = "GFF_OR_GTF")]
    pub annotation: PathBuf,
    #[arg(long = "gene", value_name = "GENE", required = true, num_args = 1..)]
    pub genes: Vec<String>,
    #[arg(long, value_name = "GENE")]
    pub partner_gene: Option<String>,
    #[arg(long, value_name = "TSV", default_value = "stellerator.tsv")]
    pub output_tsv: PathBuf,
    #[arg(long, value_name = "FASTA_GZ", default_value = "stellerator.fasta.gz")]
    pub output_fasta: PathBuf,
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
