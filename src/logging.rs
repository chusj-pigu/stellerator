use std::{
    fs::File,
    io::{self, Stderr, Write},
    path::Path,
    sync::{Arc, Mutex},
};

use anyhow::{Context, Result};
use tracing::level_filters::LevelFilter;
use tracing_subscriber::{fmt, fmt::writer::MakeWriter};

use crate::cli::Args;

pub fn init(args: &Args) -> Result<()> {
    let level = if args.verbose {
        LevelFilter::DEBUG
    } else {
        LevelFilter::INFO
    };

    let file = args
        .log_file
        .as_ref()
        .map(|path| open_log_file(path))
        .transpose()?;

    let writer = BroadcastMakeWriter { file };

    fmt()
        .with_max_level(level)
        .with_target(false)
        .with_ansi(false)
        .with_writer(writer)
        .init();

    Ok(())
}

fn open_log_file(path: &Path) -> Result<Arc<Mutex<File>>> {
    let file = File::create(path)
        .with_context(|| format!("failed to create log file {}", path.display()))?;
    Ok(Arc::new(Mutex::new(file)))
}

#[derive(Clone, Default)]
struct BroadcastMakeWriter {
    file: Option<Arc<Mutex<File>>>,
}

impl<'a> MakeWriter<'a> for BroadcastMakeWriter {
    type Writer = BroadcastWriter;

    fn make_writer(&'a self) -> Self::Writer {
        BroadcastWriter {
            stderr: io::stderr(),
            file: self.file.clone(),
        }
    }
}

struct BroadcastWriter {
    stderr: Stderr,
    file: Option<Arc<Mutex<File>>>,
}

impl Write for BroadcastWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.stderr.write_all(buf)?;

        if let Some(file) = &self.file {
            let mut file = file
                .lock()
                .map_err(|_| io::Error::other("log file mutex was poisoned"))?;
            file.write_all(buf)?;
        }

        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.stderr.flush()?;

        if let Some(file) = &self.file {
            let mut file = file
                .lock()
                .map_err(|_| io::Error::other("log file mutex was poisoned"))?;
            file.flush()?;
        }

        Ok(())
    }
}
