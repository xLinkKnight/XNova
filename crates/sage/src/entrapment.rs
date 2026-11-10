//! Entrapment: how many false discoveries a reported FDR really lets through,
//! estimated from proteins that cannot be in the sample (for example another
//! species' proteome) added to the FASTA. Any discovery that maps only to
//! such proteins is false. Following Wen et al. 2025 ("Assessment of false
//! discovery rate control in tandem mass spectrometry analysis using
//! entrapment", Nature Methods):
//!
//! * lower bound: `N_e / (N_t + N_e)`
//! * combined estimate: `N_e (1 + 1/r) / (N_t + N_e)`
//!
//! where `N_t` counts discoveries of the original proteins, `N_e` of the
//! entrapment proteins, and `r` is the size of the entrapment part of the
//! database relative to the original part. A discovery shared by an original
//! and an entrapment protein counts as original.

use serde::Serialize;
use std::sync::Arc;

/// `entrapment` configuration: proteins whose accession contains `tag`
#[derive(Clone, Debug, PartialEq, serde::Deserialize, Serialize)]
pub struct EntrapmentConfig {
    pub tag: String,
}

/// Does this discovery map only to entrapment proteins?
pub fn is_entrapment(proteins: &[Arc<str>], tag: &str) -> bool {
    !proteins.is_empty() && proteins.iter().all(|p| p.contains(tag))
}

/// One level (PSMs, peptides, proteins, crosslinks) at one reported FDR
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct EntrapmentEstimate {
    pub level: String,
    /// The reported FDR (q-value cutoff)
    pub fdr: f32,
    /// Discoveries of the original proteins
    pub targets: usize,
    /// Discoveries of the entrapment proteins only
    pub entrapments: usize,
    /// Entrapment-to-original database size ratio `r`
    pub ratio: f64,
    /// `N_e / (N_t + N_e)`: at least this share is false
    pub lower_bound: f64,
    /// `N_e (1 + 1/r) / (N_t + N_e)`: the estimated false discovery proportion
    pub combined: f64,
}

pub fn estimate(level: &str, fdr: f32, targets: usize, entrapments: usize, ratio: f64) -> EntrapmentEstimate {
    let all = (targets + entrapments) as f64;
    let (lower_bound, combined) = if all > 0.0 {
        let lower = entrapments as f64 / all;
        let combined = if ratio > 0.0 { entrapments as f64 * (1.0 + 1.0 / ratio) / all } else { f64::NAN };
        (lower, combined.min(1.0))
    } else {
        (0.0, 0.0)
    };
    EntrapmentEstimate { level: level.into(), fdr, targets, entrapments, ratio, lower_bound, combined }
}

/// The ratio for pairs of peptides (crosslinks) when each partner is drawn
/// at random: the share of pairs with an entrapment partner is
/// `1 - 1/(1 + r)^2`, so the pair ratio is `(1 + r)^2 - 1`. An approximation:
/// real crosslinks are not random pairs.
pub fn pair_ratio(ratio: f64) -> f64 {
    (1.0 + ratio).powi(2) - 1.0
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn estimators() {
        // 990 original and 10 entrapment discoveries, entrapment database half the size
        let e = estimate("peptides", 0.01, 990, 10, 0.5);
        assert!((e.lower_bound - 0.01).abs() < 1e-12);
        assert!((e.combined - 0.03).abs() < 1e-12);
        assert_eq!(estimate("psms", 0.01, 0, 0, 1.0).combined, 0.0);
        assert!(estimate("psms", 0.01, 1, 9, 0.01).combined <= 1.0);
        assert!((pair_ratio(1.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn shared_discoveries_are_original() {
        let p = |v: &[&str]| v.iter().map(|s| Arc::from(*s)).collect::<Vec<Arc<str>>>();
        assert!(is_entrapment(&p(&["ENTRAP_P1", "ENTRAP_P2"]), "ENTRAP_"));
        assert!(!is_entrapment(&p(&["ENTRAP_P1", "sp|P02769|ALBU_BOVIN"]), "ENTRAP_"));
        assert!(!is_entrapment(&p(&[]), "ENTRAP_"));
    }
}
