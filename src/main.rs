mod annotation;
mod cli;
mod fasta;
mod loci;
mod logging;
mod pipeline;
mod vcf;

use anyhow::Result;

fn main() -> Result<()> {
    let args = cli::Args::parse_args();
    pipeline::validate_output_paths(&args)?;
    logging::init(&args)?;
    pipeline::run(args)
}
