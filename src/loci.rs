use std::{fs, path::Path};

use anyhow::{Context, Result, bail};

/// One batch job read from a loci file: a query gene, an optional partner
/// constraint, and an optional clustering tolerance.
///
/// A `None` partner means "annotate the partner side against any overlapping
/// gene", matching the behaviour of omitting `--partner-gene`. A `None`
/// tolerance falls back to the global `--sv-slop`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocusRequest {
    pub gene: String,
    pub partner: Option<String>,
    pub tolerance: Option<usize>,
}

/// Read and parse a loci file into batch jobs.
pub fn parse_loci_file(path: &Path) -> Result<Vec<LocusRequest>> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("failed to read loci file {}", path.display()))?;
    parse_loci(&contents).with_context(|| format!("in loci file {}", path.display()))
}

/// Parse loci-file contents.
///
/// Each non-empty line is `gene [partner] [tolerance]`, split on any
/// whitespace. `#` starts a comment (inline or whole-line), blank lines are
/// skipped, and an optional header line (first field `gene`, case-insensitive)
/// is ignored. A partner field of `-`, `.`, `NA`/`na`, or an absent column
/// means "no partner constraint"; a tolerance that is absent or `-`/`.` falls
/// back to the global tolerance.
fn parse_loci(contents: &str) -> Result<Vec<LocusRequest>> {
    let mut requests = Vec::new();

    for (index, raw) in contents.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }

        let mut fields = line.split_whitespace();
        let gene = fields.next().expect("a non-empty line has a first field");

        // Skip a single optional header row before any data has been read.
        if requests.is_empty() && gene.eq_ignore_ascii_case("gene") {
            continue;
        }

        let partner = fields.next().and_then(value_or_none).map(str::to_string);
        let tolerance = match fields.next().and_then(value_or_none) {
            Some(value) => Some(
                value
                    .parse::<usize>()
                    .with_context(|| format!("line {}: invalid tolerance {value:?}", index + 1))?,
            ),
            None => None,
        };

        requests.push(LocusRequest {
            gene: gene.to_string(),
            partner,
            tolerance,
        });
    }

    if requests.is_empty() {
        bail!("no loci found (only comments, blank lines, or a header)");
    }

    Ok(requests)
}

/// Map an "absent value" sentinel to `None`, otherwise keep the field.
fn value_or_none(field: &str) -> Option<&str> {
    match field {
        "-" | "." | "NA" | "na" | "" => None,
        other => Some(other),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request(gene: &str, partner: Option<&str>, tolerance: Option<usize>) -> LocusRequest {
        LocusRequest {
            gene: gene.to_string(),
            partner: partner.map(str::to_string),
            tolerance,
        }
    }

    #[test]
    fn parses_full_rows_tab_or_space_separated() {
        let requests = parse_loci("BCR\tABL1\t10\nEWSR1 FLI1 200\n").unwrap();
        assert_eq!(
            requests,
            vec![
                request("BCR", Some("ABL1"), Some(10)),
                request("EWSR1", Some("FLI1"), Some(200)),
            ]
        );
    }

    #[test]
    fn skips_comments_blank_lines_and_header() {
        let text = "\
# a panel
gene\tpartner\ttolerance
BCR\tABL1\t10

MYC   # gene only
";
        let requests = parse_loci(text).unwrap();
        assert_eq!(
            requests,
            vec![
                request("BCR", Some("ABL1"), Some(10)),
                request("MYC", None, None),
            ]
        );
    }

    #[test]
    fn treats_sentinels_as_absent() {
        let requests = parse_loci("MYC - 50\nFOXO1 . -\nTP53 NA\n").unwrap();
        assert_eq!(
            requests,
            vec![
                request("MYC", None, Some(50)),
                request("FOXO1", None, None),
                request("TP53", None, None),
            ]
        );
    }

    #[test]
    fn a_partner_without_a_tolerance_falls_back() {
        let requests = parse_loci("BCR ABL1\n").unwrap();
        assert_eq!(requests, vec![request("BCR", Some("ABL1"), None)]);
    }

    #[test]
    fn rejects_a_non_numeric_tolerance() {
        let error = parse_loci("BCR ABL1 wide\n").unwrap_err();
        assert!(error.to_string().contains("invalid tolerance"), "{error}");
    }

    #[test]
    fn rejects_a_file_with_no_loci() {
        let error = parse_loci("# only comments\n\n").unwrap_err();
        assert!(error.to_string().contains("no loci"), "{error}");
    }
}
