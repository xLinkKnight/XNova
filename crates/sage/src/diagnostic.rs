//! Diagnostic ions: low-mass fragments that mark a modification class, such
//! as glycan oxonium ions or modified-residue immonium ions.
//!
//! Each ion is looked up with one binary search per spectrum on the raw,
//! sorted fragment m/z array, before deisotoping or peak trimming.
//!
//! Adapted from Sage Plus (`crates/sage/src/diagnostic.rs`, MIT,
//! <https://github.com/pgarrett-scripps/sage-plus>).

use crate::mass::Tolerance;
use serde::{Deserialize, Serialize};

/// Default matching tolerance of a diagnostic ion
pub const DEFAULT_TOLERANCE: Tolerance = Tolerance::Ppm(-20.0, 20.0);

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiagnosticIon {
    /// Name reported in the output
    pub name: String,
    /// Singly charged m/z
    pub mz: f32,
    /// Matching tolerance; 20 ppm when omitted
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<Tolerance>,
}

impl DiagnosticIon {
    fn new(name: &str, mz: f32) -> Self {
        Self { name: name.into(), mz, tolerance: None }
    }

    pub fn tolerance(&self) -> Tolerance {
        self.tolerance.unwrap_or(DEFAULT_TOLERANCE)
    }
}

/// Built-in ions: glycan oxonium ions, and the acetyl-lysine and
/// phosphotyrosine immonium ions
pub fn default_ions() -> Vec<DiagnosticIon> {
    vec![
        DiagnosticIon::new("HexNAc", 204.0867),
        DiagnosticIon::new("HexNAc_fragment", 138.055),
        DiagnosticIon::new("Hex", 163.0601),
        DiagnosticIon::new("NeuAc", 292.1027),
        DiagnosticIon::new("acetyl_K_immonium", 126.0913),
        DiagnosticIon::new("phospho_Y_immonium", 216.0426),
    ]
}

/// `diagnostic_ions` configuration: `true` for the built-in ions, or a list
/// of ions that replaces them
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(untagged)]
pub enum DiagnosticIonsConfig {
    Enabled(bool),
    Ions(Vec<DiagnosticIon>),
}

impl DiagnosticIonsConfig {
    /// Ions to search for; `None` when disabled
    pub fn resolve(self) -> Option<Vec<DiagnosticIon>> {
        match self {
            Self::Enabled(false) => None,
            Self::Enabled(true) => Some(default_ions()),
            Self::Ions(ions) => Some(ions),
        }
    }
}

/// One diagnostic ion found in a spectrum
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiagnosticHit {
    /// Index into the ion list
    pub ion: usize,
    /// Observed m/z of the matched peak
    pub mz: f32,
    /// Matched peak intensity divided by the spectrum's summed intensity
    pub relative_intensity: f32,
}

/// Find `ions` in a spectrum. `mz` must be sorted ascending and parallel to
/// `intensity`. When several peaks fall in an ion's window, the most intense
/// one is reported.
pub fn find_ions(ions: &[DiagnosticIon], mz: &[f32], intensity: &[f32]) -> Vec<DiagnosticHit> {
    let total: f32 = intensity.iter().sum();
    if total <= 0.0 {
        return Vec::new();
    }
    ions.iter()
        .enumerate()
        .filter_map(|(index, ion)| {
            let (low, high) = ion.tolerance().bounds(ion.mz);
            let start = mz.partition_point(|&value| value < low);
            let best = (start..mz.len())
                .take_while(|&peak| mz[peak] <= high)
                .max_by(|&left, &right| intensity[left].total_cmp(&intensity[right]))?;
            Some(DiagnosticHit { ion: index, mz: mz[best], relative_intensity: intensity[best] / total })
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn configuration() {
        let parse = |v: serde_json::Value| serde_json::from_value::<DiagnosticIonsConfig>(v).unwrap().resolve();
        assert_eq!(parse(serde_json::json!(false)), None);
        assert_eq!(parse(serde_json::json!(true)).unwrap().len(), 6);
        let custom = parse(serde_json::json!([{"name": "TMT126", "mz": 126.127726, "tolerance": {"da": [-0.01, 0.01]}}])).unwrap();
        assert_eq!((custom[0].name.as_str(), custom[0].tolerance()), ("TMT126", Tolerance::Da(-0.01, 0.01)));
        assert_eq!(default_ions()[0].tolerance(), DEFAULT_TOLERANCE);
    }

    #[test]
    fn the_most_intense_peak_in_the_window() {
        let ions = default_ions();
        let mz = [126.0911, 126.0914, 204.0867, 204.5, 300.0];
        let intensity = [5.0, 10.0, 20.0, 15.0, 50.0];
        let hits = find_ions(&ions, &mz, &intensity);
        let names: Vec<&str> = hits.iter().map(|h| ions[h.ion].name.as_str()).collect();
        assert_eq!(names, ["HexNAc", "acetyl_K_immonium"]);
        let ack = hits.iter().find(|h| ions[h.ion].name == "acetyl_K_immonium").unwrap();
        assert_eq!(ack.mz, 126.0914);
        assert!((ack.relative_intensity - 0.1).abs() < 1e-6);
        assert!(find_ions(&ions, &mz, &[0.0; 5]).is_empty());
    }
}
