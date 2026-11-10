//! Polymer contamination in MS1 spectra.
//!
//! Detergents and plastics ionize as ladders of peaks spaced by one repeat
//! unit. A spectrum's polymer signal is the intensity of the peaks that belong
//! to a run of consecutive ladder members, for any charge state and adduct.
//! Masses are monoisotopic, from the elemental formulas of the repeat units.
//!
//! Adapted from Sage Plus (`crates/sage/src/polymer.rs`, MIT,
//! <https://github.com/pgarrett-scripps/sage-plus>).

use serde::Serialize;

/// Monoisotopic mass of water, the end groups of linear PEG and PPG
const WATER: f64 = 18.010_565;

/// A polymer series: `end_group + n * repeat` for `n >= 1`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Polymer {
    pub name: &'static str,
    /// Neutral mass of the non-repeating part (zero for cyclic polymers)
    pub end_group: f64,
    /// Neutral mass of one repeat unit
    pub repeat: f64,
}

/// Polyethylene glycol, H-(C2H4O)n-OH
pub const PEG: Polymer = Polymer { name: "PEG", end_group: WATER, repeat: 44.026_215 };

/// Polypropylene glycol, H-(C3H6O)n-OH
pub const PPG: Polymer = Polymer { name: "PPG", end_group: WATER, repeat: 58.041_865 };

/// Cyclic polydimethylsiloxane, (C2H6OSi)n, e.g. the 445.12 lock mass
pub const POLYSILOXANE: Polymer = Polymer { name: "polysiloxane", end_group: 0.0, repeat: 74.018_792 };

pub const POLYMERS: [Polymer; 3] = [PEG, PPG, POLYSILOXANE];

/// Charge carriers: proton, sodium and ammonium cations
const ADDUCTS: [f64; 3] = [1.007_276, 22.989_218, 18.033_823];

/// MS1 matching tolerance (ppm)
pub const TOLERANCE_PPM: f64 = 10.0;

/// Consecutive ladder members required before any of them count
pub const MIN_LADDER: usize = 4;

/// Highest charge state searched
pub const MAX_CHARGE: u8 = 3;

/// A polymer carrying more than this percent of a file's MS1 signal is warned about
pub const WARN_PCT: f64 = 5.0;

/// Index of the most intense peak within `TOLERANCE_PPM` of `target`;
/// `mz` sorted ascending
fn best_peak(mz: &[f32], intensity: &[f32], target: f64) -> Option<usize> {
    let delta = target * TOLERANCE_PPM * 1e-6;
    let (low, high) = (target - delta, target + delta);
    let start = mz.partition_point(|&value| (value as f64) < low);
    (start..mz.len())
        .take_while(|&index| mz[index] as f64 <= high)
        .max_by(|&left, &right| intensity[left].total_cmp(&intensity[right]))
}

/// Intensity of the peaks on each polymer's ladders, in the order of
/// [`POLYMERS`]. `mz` must be sorted ascending and parallel to `intensity`.
/// A peak counts once per polymer, even when several charge states or adducts
/// claim it.
pub fn scan(mz: &[f32], intensity: &[f32]) -> [f64; 3] {
    let mut out = [0.0; 3];
    let (Some(&first), Some(&last)) = (mz.first(), mz.last()) else {
        return out;
    };
    // widen the peak range by the tolerance so ladder ends are kept
    let slack = TOLERANCE_PPM * 1e-6;
    let (first, last) = (first as f64 * (1.0 - slack), last as f64 * (1.0 + slack));
    let mut claimed = vec![false; mz.len()];
    let mut run = Vec::new();
    for (k, polymer) in POLYMERS.iter().enumerate() {
        claimed.iter_mut().for_each(|flag| *flag = false);
        for charge in 1..=MAX_CHARGE {
            let z = charge as f64;
            for adduct in ADDUCTS {
                let offset = polymer.end_group + z * adduct;
                let start = (((first * z - offset) / polymer.repeat).ceil() as i64).max(1);
                let mut n = start as f64;
                run.clear();
                loop {
                    let target = (offset + n * polymer.repeat) / z;
                    let hit = if target > last { None } else { best_peak(mz, intensity, target) };
                    match hit {
                        Some(index) => run.push(index),
                        None => {
                            if run.len() >= MIN_LADDER {
                                run.iter().for_each(|&index| claimed[index] = true);
                            }
                            run.clear();
                        }
                    }
                    if target > last {
                        break;
                    }
                    n += 1.0;
                }
            }
        }
        out[k] = claimed.iter().zip(intensity).filter(|(flag, _)| **flag).map(|(_, &v)| v as f64).sum();
    }
    out
}

/// Polymer signal of one file's MS1 spectra
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct PolymerFileStats {
    /// Centroided MS1 spectra scanned
    pub ms1_spectra: usize,
    /// Profile MS1 spectra, which are skipped
    pub skipped_profile_spectra: usize,
    /// Summed intensity of all scanned MS1 peaks
    pub total_ion_current: f64,
    /// Each polymer's ladder intensity and its percent of `total_ion_current`
    pub polymers: Vec<PolymerShare>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct PolymerShare {
    pub name: String,
    pub intensity: f64,
    pub tic_pct: f64,
}

impl PolymerFileStats {
    pub fn new() -> Self {
        Self {
            polymers: POLYMERS.iter().map(|p| PolymerShare { name: p.name.into(), ..Default::default() }).collect(),
            ..Default::default()
        }
    }

    /// Add one centroided spectrum
    pub fn add(&mut self, mz: &[f32], intensity: &[f32]) {
        self.ms1_spectra += 1;
        self.total_ion_current += intensity.iter().map(|&v| v as f64).sum::<f64>();
        for (share, value) in self.polymers.iter_mut().zip(scan(mz, intensity)) {
            share.intensity += value;
        }
    }

    pub fn merge(&mut self, other: &Self) {
        self.ms1_spectra += other.ms1_spectra;
        self.skipped_profile_spectra += other.skipped_profile_spectra;
        self.total_ion_current += other.total_ion_current;
        for (share, other) in self.polymers.iter_mut().zip(&other.polymers) {
            share.intensity += other.intensity;
        }
    }

    /// Fill in `tic_pct` from the accumulated intensities
    pub fn finish(&mut self) {
        for share in &mut self.polymers {
            share.tic_pct =
                if self.total_ion_current > 0.0 { 100.0 * share.intensity / self.total_ion_current } else { 0.0 };
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// m/z of member `n` of a ladder
    fn member(p: Polymer, n: f64, z: f64, adduct: f64) -> f32 {
        ((p.end_group + n * p.repeat + z * adduct) / z) as f32
    }

    #[test]
    fn ladders_count_from_four_members() {
        // a PEG ladder [M+H]+ n = 6..=11, a lone PPG peak, a short siloxane run, background
        let mut peaks: Vec<(f32, f32)> = (6..=11).map(|n| (member(PEG, n as f64, 1.0, ADDUCTS[0]), 100.0)).collect();
        peaks.push((member(PPG, 7.0, 1.0, ADDUCTS[0]), 50.0));
        peaks.extend((5..=7).map(|n| (member(POLYSILOXANE, n as f64, 1.0, ADDUCTS[0]), 30.0)));
        peaks.extend([(500.2731, 1000.0), (612.3302, 400.0)]);
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (mz, intensity): (Vec<f32>, Vec<f32>) = peaks.into_iter().unzip();
        let found = scan(&mz, &intensity);
        assert_eq!(found[0], 600.0, "the six PEG peaks");
        assert_eq!(found[1], 0.0, "one PPG peak is not a ladder");
        assert_eq!(found[2], 0.0, "three siloxane peaks are too few");

        let mut stats = PolymerFileStats::new();
        stats.add(&mz, &intensity);
        stats.finish();
        assert_eq!(stats.ms1_spectra, 1);
        let total: f64 = intensity.iter().map(|&v| v as f64).sum();
        assert!((stats.polymers[0].tic_pct - 100.0 * 600.0 / total).abs() < 1e-9);
    }

    #[test]
    fn doubly_charged_sodiated_ladders_count() {
        let mut peaks: Vec<(f32, f32)> = (20..=25).map(|n| (member(PEG, n as f64, 2.0, ADDUCTS[1]), 10.0)).collect();
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (mz, intensity): (Vec<f32>, Vec<f32>) = peaks.into_iter().unzip();
        assert_eq!(scan(&mz, &intensity)[0], 60.0);
        assert_eq!(scan(&[], &[]), [0.0; 3]);
    }
}
