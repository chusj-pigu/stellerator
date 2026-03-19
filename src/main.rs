mod annotation;
mod cli;
mod fasta;
mod logging;
mod pipeline;

use anyhow::Result;

fn main() -> Result<()> {
    let args = cli::Args::parse_args();
    logging::init(&args)?;
    pipeline::run(args)
}
