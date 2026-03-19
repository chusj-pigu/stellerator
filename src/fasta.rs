use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::{Context, Result};
use flate2::{Compression, write::GzEncoder};

pub struct FastaWriter {
    writer: GzEncoder<BufWriter<File>>,
}

impl FastaWriter {
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("failed to create FASTA output {}", path.display()))?;
        let writer = GzEncoder::new(BufWriter::new(file), Compression::default());
        Ok(Self { writer })
    }

    pub fn write_record(&mut self, header: &str, sequence: &str) -> Result<()> {
        writeln!(self.writer, ">{header}")?;
        writeln!(self.writer, "{sequence}")?;
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.writer.try_finish()?;
        Ok(())
    }
}
