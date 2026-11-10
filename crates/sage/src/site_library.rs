//! Site libraries: modification sites known on proteins, for example the
//! `sites.tsv` of an earlier search with `localization`, or a curated list.
//! A variable modification with `site_mode: library` is then searched only
//! at the sites the library lists for it (after the PTM site libraries of
//! Sage Plus, MIT):
//!
//! ```json
//! "database": {
//!   "site_library": {"path": "sites.tsv"},
//!   "variable_mods": {
//!     "S": [{"mass": 79.966331, "name": "UNIMOD:21", "site_mode": "library"}]
//!   }
//! }
//! ```
//!
//! The library is a tab-separated table with the columns `protein` (the
//! FASTA accession), `position` (1-based on the protein), `residue` and
//! `modification`: the modification's name (`UNIMOD:21`) or mass
//! (`+79.966331`), as `sites.tsv` writes it. Other columns are ignored, but
//! with a `site_q` column only sites at or below `max_q` (default 0.01) are
//! used. A site is found on a peptide in each protein it belongs to; decoys
//! are reversed after modification, so they keep their target's sites.

use crate::modification::{ModId, Modification};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

fn default_max_q() -> f32 {
    0.01
}

/// `database.site_library`
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct SiteLibrarySettings {
    /// The table (a local path or URL)
    pub path: String,
    /// With a `site_q` column, use only the sites at or below it
    #[serde(default = "default_max_q")]
    pub max_q: f32,
}

/// One row of a site library
#[derive(Clone, Debug, PartialEq)]
pub struct LibrarySite {
    pub protein: String,
    /// 1-based position on the protein
    pub position: usize,
    pub residue: u8,
    pub modification: String,
}

/// Read a site library table
pub fn parse(text: &str, max_q: f32) -> Result<Vec<LibrarySite>, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let header: Vec<String> = lines
        .next()
        .ok_or("the site library is empty")?
        .split('\t')
        .map(|h| h.trim().trim_start_matches('\u{feff}').to_ascii_lowercase())
        .collect();
    let column = |name: &str| header.iter().position(|h| h == name);
    let (Some(protein), Some(position), Some(residue), Some(modification)) =
        (column("protein"), column("position"), column("residue"), column("modification"))
    else {
        return Err("a site library needs the columns protein, position, residue and modification".into());
    };
    let site_q = column("site_q");
    let mut sites = Vec::new();
    for (n, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').map(str::trim).collect();
        let field = |i: usize| fields.get(i).copied().unwrap_or("");
        let row = n + 2;
        if let Some(q) = site_q {
            let q: f32 = field(q).parse().map_err(|_| format!("site library row {row}: site_q is not a number"))?;
            if q > max_q {
                continue;
            }
        }
        let at: usize = field(position)
            .parse()
            .ok()
            .filter(|&p| p > 0)
            .ok_or_else(|| format!("site library row {row}: position must be a 1-based number"))?;
        let r = field(residue).as_bytes();
        if r.len() != 1 || !crate::mass::VALID_AA.contains(&r[0].to_ascii_uppercase()) {
            return Err(format!("site library row {row}: residue must be one amino acid"));
        }
        if field(protein).is_empty() || field(modification).is_empty() {
            return Err(format!("site library row {row}: protein and modification may not be empty"));
        }
        sites.push(LibrarySite {
            protein: field(protein).to_string(),
            position: at,
            residue: r[0].to_ascii_uppercase(),
            modification: field(modification).to_string(),
        });
    }
    Ok(sites)
}

/// Does a library's modification label name this modification? Its name
/// (ignoring case), or its mass within 0.001 Da
pub fn names(modification: &Modification, label: &str) -> bool {
    let label = label.trim();
    if modification.name.as_ref().is_some_and(|n| n.as_str().eq_ignore_ascii_case(label)) {
        return true;
    }
    label
        .trim_start_matches('+')
        .parse::<f32>()
        .is_ok_and(|mass| (mass - modification.mass).abs() <= 1e-3)
}

/// A site library resolved to a search's modifications
#[derive(Clone, Debug, Default)]
pub struct ResolvedLibrary {
    /// By protein accession: (0-based position, residue, modification)
    pub sites: HashMap<String, Vec<(usize, u8, ModId)>>,
    /// Rows read, and rows used (naming a `site_mode: library` modification
    /// on their residue)
    pub rows: usize,
    pub used: usize,
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn reads_sites_tsv() {
        let text = "protein\tposition\tresidue\tmodification\tprobability\tpsms\tpeptides\tlocalization_q\tsite_q\n\
                    P1\t12\tS\tUNIMOD:21\t0.99\t3\t1\t0\t0.001\n\
                    P1\t15\tT\t+79.96633\t0.6\t0\t0\t1\t1\n\
                    P2\t3\ty\tUNIMOD:21\t0.9\t1\t1\t0\t0.01\n";
        let sites = parse(text, 0.01).unwrap();
        assert_eq!(sites.len(), 2, "the site above max_q is left out");
        assert_eq!(sites[1], LibrarySite { protein: "P2".into(), position: 3, residue: b'Y', modification: "UNIMOD:21".into() });
        assert_eq!(parse(text, 1.0).unwrap().len(), 3);
        // without site_q, every row
        let plain = "Protein\tPosition\tResidue\tModification\nP1\t1\tM\t+15.9949\n";
        assert_eq!(parse(plain, 0.01).unwrap()[0].position, 1);
        assert!(parse("protein\tposition\n", 0.01).is_err());
        assert!(parse("protein\tposition\tresidue\tmodification\nP1\t0\tS\tx\n", 0.01).is_err());
        assert!(parse("protein\tposition\tresidue\tmodification\nP1\t2\tSS\tx\n", 0.01).is_err());

        let phospho = Modification { mass: 79.966331, name: Some("UNIMOD:21".parse().unwrap()) };
        assert!(names(&phospho, "unimod:21") && names(&phospho, "+79.96633") && names(&phospho, "79.9663"));
        assert!(!names(&phospho, "UNIMOD:35") && !names(&phospho, "+15.9949"));
    }
}
