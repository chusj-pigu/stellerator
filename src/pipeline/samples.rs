use super::*;

pub(super) fn validate_inputs(args: &Args) -> Result<()> {
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }

    if !args.annotation.exists() {
        bail!(
            "annotation input does not exist: {}",
            args.annotation.display()
        );
    }

    Ok(())
}

/// A single indexed BAM to scan, carrying the sample name used for output
/// provenance and its parsed header.
pub(super) struct BamSample {
    pub(super) path: PathBuf,
    pub(super) name: String,
    pub(super) header: sam::Header,
}

impl BamSample {
    fn open(path: &Path) -> Result<Self> {
        let header = load_sam_header(path)?;
        Ok(Self {
            path: path.to_path_buf(),
            name: sample_name(path),
            header,
        })
    }
}

/// Resolve the requested BAM inputs (files and/or directories) into a validated,
/// deduplicated set of indexed BAMs and open each as a [`BamSample`].
pub(super) fn open_bam_samples(inputs: &[PathBuf]) -> Result<Vec<BamSample>> {
    let paths = resolve_bam_inputs(inputs)?;
    check_unique_sample_names(&paths)?;
    paths.iter().map(|path| BamSample::open(path)).collect()
}

/// Expand directory inputs to their `.bam` contents, confirm every BAM exists
/// and is indexed, then sort and deduplicate the result.
pub(super) fn resolve_bam_inputs(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut resolved: Vec<PathBuf> = Vec::new();

    for input in inputs {
        if input.is_dir() {
            let mut found = directory_bam_files(input)?;
            if found.is_empty() {
                bail!("no .bam files found in directory {}", input.display());
            }
            resolved.append(&mut found);
        } else if input.exists() {
            resolved.push(input.clone());
        } else {
            bail!("BAM input does not exist: {}", input.display());
        }
    }

    resolved.sort();
    resolved.dedup();

    if resolved.is_empty() {
        bail!("no BAM inputs provided");
    }

    for path in &resolved {
        if !has_associated_index(path) {
            bail!(
                "indexed BAM required; expected {path}.bai or {path}.csi",
                path = path.display()
            );
        }
    }

    Ok(resolved)
}

/// List the `.bam` files directly inside `dir`, sorted by path.
pub(super) fn directory_bam_files(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in
        fs::read_dir(dir).with_context(|| format!("failed to read directory {}", dir.display()))?
    {
        let path = entry?.path();
        if path.is_file() && has_bam_extension(&path) {
            files.push(path);
        }
    }

    files.sort();
    Ok(files)
}

pub(super) fn has_bam_extension(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
}

/// Derive the sample name for a BAM from its file stem, falling back to the
/// full path when no stem is available.
pub(super) fn sample_name(path: &Path) -> String {
    path.file_stem()
        .and_then(|stem| stem.to_str())
        .map(|stem| stem.to_string())
        .unwrap_or_else(|| path.to_string_lossy().into_owned())
}

/// Reject input sets where two BAMs would collapse to the same sample name,
/// which would make per-sample output provenance ambiguous.
pub(super) fn check_unique_sample_names(paths: &[PathBuf]) -> Result<()> {
    let mut seen: HashMap<String, &Path> = HashMap::new();

    for path in paths {
        let name = sample_name(path);
        if let Some(existing) = seen.insert(name.clone(), path) {
            bail!(
                "duplicate sample name {name:?} derived from {} and {}; rename one of the BAM files to disambiguate",
                existing.display(),
                path.display()
            );
        }
    }

    Ok(())
}

pub(super) fn load_sam_header(path: &Path) -> Result<sam::Header> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .with_context(|| format!("failed to open BAM {}", path.display()))?;
    reader.read_header().context("failed to read BAM header")
}

pub(super) fn has_associated_index(bam_path: &Path) -> bool {
    let bam_bai = bam_path.with_extension(format!(
        "{}.bai",
        bam_path
            .extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or_default()
    ));
    let bam_csi = bam_path.with_extension(format!(
        "{}.csi",
        bam_path
            .extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or_default()
    ));

    bam_bai.exists() || bam_csi.exists()
}
