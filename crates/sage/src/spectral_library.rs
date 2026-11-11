//! Spectral libraries from search results, for DIA searches in DIA-NN,
//! Spectronaut and PEAKS.
//!
//! The linear library follows EasyPQP's method (`easypqp convertsage` +
//! `easypqp library`, defaults): target PSMs at 1% PSM, peptide and protein
//! FDR, proteotypic peptides only; for each precursor (modified peptide and
//! charge) the best PSM of every run, by posterior error; retention times on
//! one scale across runs (here XNova's own run alignment, scaled 0-100, as
//! EasyPQP scales its reference run 0-100) and ion mobilities calibrated to
//! the run with the most precursors; then a consensus: the median retention
//! time and ion mobility of the runs, and each fragment's median intensity
//! (every PSM's matched fragments scaled to a maximum of 10,000).
//!
//! Crosslinks (`xlibrary`) are written as Rappsilber Laboratory's
//! xiDIA-library writes them for Spectronaut: one precursor per crosslinked
//! pair, the two peptides in alphabetical order joined by `_`, the
//! crosslinker as a modification on each linked residue
//! (`PEPK[DSSO]TIDE_LIK[DSSO]ER`), with the crosslink's own precursor m/z and
//! fragment m/z, and the peptide each fragment comes from. No DIA tool knows
//! crosslinks: they search such a library like any other, by the masses it
//! gives. DIA-NN uses the decoys of a library that has some, instead of
//! making its own (which would be meaningless for a pair), so the OpenSWATH
//! form gets one decoy per crosslink, both peptides reversed; a decoy peptide
//! never has a target peptide's sequence (I and L alike), so that no true hit
//! is taken for a decoy. The linear library has no decoys: DIA tools make
//! their own.
//!
//! EasyPQP: https://github.com/grosenberger/easypqp (George Rosenberger,
//! BSD-3-Clause). xiDIA-library:
//! https://github.com/Rappsilber-Laboratory/xiDIA-library (Apache-2.0;
//! Muller et al., Mol Cell Proteomics 2019, doi:10.1074/mcp.TIR118.001276).
//! Their methods and formats are reimplemented here; no code of either is
//! included.
//!
//! Two layouts are written:
//! * OpenSWATH (`library.tsv`): EasyPQP's columns plus the fragment columns of
//!   its DIA-NN option. DIA-NN reads it as it reads FragPipe's libraries,
//!   PEAKS Studio's Library Viewer as an "OpenMS generated library", OpenSWATH
//!   through TargetedFileConverter.
//! * Spectronaut (`library.spectronaut.xls`): the columns of a Spectronaut
//!   library export, tab-separated. Spectronaut imports it, DIA-NN reads
//!   Spectronaut .xls libraries as they are, and PEAKS Studio's Library Viewer
//!   reads "Spectronaut libraries with the same format and column headers".

use crate::crosslink::{CrosslinkMatch, LinkKind, LinkedFragment, PairClass};
use crate::database::{isobaric_key, IndexedDatabase, PeptideIx};
use crate::ion_series::{IonSeries, Kind};
use crate::mass::PROTON;
use crate::modification::{ModId, Modification, Modifications};
use crate::peptide::Peptide;
use crate::scoring::Feature;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// `spectral_library`: on with `{}`
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LibrarySettings {
    /// Layouts to write (default both)
    #[serde(default = "default_formats")]
    pub formats: Vec<LibraryFormat>,
    /// Maximum `spectrum_q` of a PSM in the library
    #[serde(default = "default_fdr")]
    pub psm_fdr: f32,
    /// Maximum `peptide_q`
    #[serde(default = "default_fdr")]
    pub peptide_fdr: f32,
    /// Maximum `protein_q`
    #[serde(default = "default_fdr")]
    pub protein_fdr: f32,
    /// Crosslinks: maximum `spectrum_q` and `crosslink_q` of a CSM
    #[serde(default = "default_fdr")]
    pub crosslink_fdr: f32,
    /// Only peptides of a single protein (EasyPQP's default)
    #[serde(default = "default_true")]
    pub proteotypic: bool,
    /// Median over the runs' best PSMs (EasyPQP's default); otherwise the
    /// single best PSM of all runs
    #[serde(default = "default_true")]
    pub consensus: bool,
    /// Precursors with fewer matched fragments are left out
    #[serde(default = "default_min_fragments")]
    pub min_fragments: usize,
    /// The most intense fragments kept per precursor
    #[serde(default = "default_max_fragments")]
    pub max_fragments: usize,
    /// Retention times in iRT, calibrated with Biognosys' iRT kit peptides
    /// found by the search; otherwise (and when too few are found) the runs'
    /// aligned times scaled 0-100
    #[serde(default)]
    pub irt: bool,
    /// Crosslinks: one decoy per crosslink in the OpenSWATH layout (for DIA-NN)
    #[serde(default = "default_true")]
    pub crosslink_decoys: bool,
}

fn default_formats() -> Vec<LibraryFormat> {
    vec![LibraryFormat::Openswath, LibraryFormat::Spectronaut]
}
fn default_fdr() -> f32 {
    0.01
}
fn default_true() -> bool {
    true
}
fn default_min_fragments() -> usize {
    3
}
fn default_max_fragments() -> usize {
    20
}

impl Default for LibrarySettings {
    fn default() -> Self {
        Self {
            formats: default_formats(),
            psm_fdr: default_fdr(),
            peptide_fdr: default_fdr(),
            protein_fdr: default_fdr(),
            crosslink_fdr: default_fdr(),
            proteotypic: true,
            consensus: true,
            min_fragments: default_min_fragments(),
            max_fragments: default_max_fragments(),
            irt: false,
            crosslink_decoys: true,
        }
    }
}

impl LibrarySettings {
    pub fn validate(&self) -> Result<(), String> {
        for (name, q) in [
            ("psm_fdr", self.psm_fdr),
            ("peptide_fdr", self.peptide_fdr),
            ("protein_fdr", self.protein_fdr),
            ("crosslink_fdr", self.crosslink_fdr),
        ] {
            if !(q > 0.0 && q <= 1.0) {
                return Err(format!("spectral_library.{name} must be above 0 and at most 1, not {q}"));
            }
        }
        if self.min_fragments == 0 || self.max_fragments < self.min_fragments {
            return Err(format!(
                "spectral_library: min_fragments ({}) must be at least 1 and at most max_fragments ({})",
                self.min_fragments, self.max_fragments
            ));
        }
        if self.formats.is_empty() {
            return Err("spectral_library.formats is empty".into());
        }
        Ok(())
    }

    pub fn writes(&self, format: LibraryFormat) -> bool {
        self.formats.contains(&format)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum LibraryFormat {
    /// `library.tsv`: DIA-NN, PEAKS (OpenMS library), OpenSWATH
    Openswath,
    /// `library.spectronaut.xls`: Spectronaut, DIA-NN, PEAKS
    Spectronaut,
}

/// A fragment ion of a library precursor
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LibraryFragment {
    pub kind: Kind,
    /// Series number: residues the ion contains
    pub number: u16,
    pub charge: u8,
    /// Calculated m/z
    pub mz: f64,
    /// Relative intensity, at most 10,000
    pub intensity: f32,
    /// Crosslinks: the peptide it comes from (0: the first in the library's order)
    pub peptide: u8,
    /// Crosslinks: carries the linker (with the partner peptide or a stub)
    pub linked: bool,
    /// Neutral mass added to the plain backbone ion (linked ions)
    pub shift: f32,
}

/// A peptide of a library precursor
#[derive(Clone, Debug)]
pub struct LibraryPeptide {
    pub peptide: Peptide,
    /// 0-based linked residues (crosslinks)
    pub links: Vec<usize>,
}

/// A library precursor with its fragments
#[derive(Clone, Debug)]
pub struct LibraryPrecursor {
    /// One peptide (linear, loop-link) or two (crosslinked pair, alphabetical)
    pub peptides: Vec<LibraryPeptide>,
    pub charge: u8,
    pub mz: f64,
    /// Retention time: 0-100 across the aligned runs, or iRT
    pub rt: f64,
    /// 1/K0 (timsTOF), if the runs report one
    pub ion_mobility: Option<f32>,
    pub decoy: bool,
    pub fragments: Vec<LibraryFragment>,
    /// Runs whose best PSM went into the entry
    pub runs: usize,
}

impl LibraryPrecursor {
    pub fn is_crosslink(&self) -> bool {
        self.peptides.iter().any(|p| !p.links.is_empty())
    }
}

/// Protein accession and entry name from a FASTA identifier (`db|ACC|ENTRY`)
pub fn accession_and_entry(token: &str) -> (String, String) {
    let parts: Vec<&str> = token.split('|').collect();
    match parts.len() {
        n if n >= 3 => (parts[1].to_string(), parts[2].to_string()),
        2 => (parts[1].to_string(), String::new()),
        _ => (token.to_string(), String::new()),
    }
}

fn accessions(peptide: &Peptide) -> (Vec<String>, Vec<String>) {
    let mut accs: Vec<String> = Vec::new();
    let mut entries: Vec<String> = Vec::new();
    for protein in &peptide.proteins {
        let (acc, entry) = accession_and_entry(protein);
        if !accs.contains(&acc) {
            accs.push(acc);
        }
        if !entry.is_empty() && !entries.contains(&entry) {
            entries.push(entry);
        }
    }
    (accs, entries)
}

/// A peptide's precursor m/z at `charge`, in f64
pub fn precursor_mz(monoisotopic: f64, charge: u8) -> f64 {
    (monoisotopic + charge as f64 * PROTON as f64) / charge as f64
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    match n {
        0 => f64::NAN,
        _ if n % 2 == 1 => values[n / 2],
        _ => 0.5 * (values[n / 2 - 1] + values[n / 2]),
    }
}

/// Least-squares line `y = slope x + intercept` with one round of outlier
/// rejection (residuals beyond 3 median absolute deviations); `None` with
/// fewer than `min_points` points
pub fn robust_line(points: &[(f64, f64)], min_points: usize) -> Option<(f64, f64, f64)> {
    let fit = |pts: &[(f64, f64)]| -> Option<(f64, f64)> {
        let n = pts.len() as f64;
        if pts.len() < 2 {
            return None;
        }
        let (mx, my) = pts.iter().fold((0.0, 0.0), |(a, b), (x, y)| (a + x, b + y));
        let (mx, my) = (mx / n, my / n);
        let sxx: f64 = pts.iter().map(|(x, _)| (x - mx).powi(2)).sum();
        if sxx <= f64::EPSILON {
            return None;
        }
        let sxy: f64 = pts.iter().map(|(x, y)| (x - mx) * (y - my)).sum();
        let slope = sxy / sxx;
        Some((slope, my - slope * mx))
    };
    if points.len() < min_points {
        return None;
    }
    let (slope, intercept) = fit(points)?;
    let mut residuals: Vec<f64> = points.iter().map(|(x, y)| (y - (slope * x + intercept)).abs()).collect();
    let mad = median(&mut residuals.clone()).max(1e-9);
    let kept: Vec<(f64, f64)> = points
        .iter()
        .zip(residuals.iter_mut())
        .filter(|(_, r)| **r <= 3.0 * 1.4826 * mad)
        .map(|(p, _)| *p)
        .collect();
    let (slope, intercept) = if kept.len() >= min_points { fit(&kept)? } else { (slope, intercept) };
    // coefficient of determination on the kept points
    let pts = if kept.len() >= min_points { &kept[..] } else { points };
    let my = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
    let ss_tot: f64 = pts.iter().map(|(_, y)| (y - my).powi(2)).sum();
    let ss_res: f64 = pts.iter().map(|(x, y)| (y - (slope * x + intercept)).powi(2)).sum();
    let r2 = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 1.0 };
    Some((slope, intercept, r2))
}

/// Biognosys iRT kit peptides and their iRT values (Escher et al. 2012)
pub const IRT_KIT: [(&str, f64); 11] = [
    ("LGGNEQVTR", -24.92),
    ("GAGSSEPVTGLDAK", 0.0),
    ("VEATFGVDESNANK", 12.39),
    ("YILAGVENSK", 19.79),
    ("TPVISGGPYEYR", 28.71),
    ("TPVITGAPYEYR", 33.38),
    ("DGLDAASYYAPVR", 42.26),
    ("ADVTPADFSEWSK", 54.62),
    ("GTFIIDPGGVIR", 70.52),
    ("GTFIIDPAAVIR", 87.23),
    ("LFLQFGAQGSPFLK", 100.0),
];

/// Retention-time calibration of a library onto iRT: `iRT = slope rt + intercept`
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IrtCalibration {
    pub slope: f64,
    pub intercept: f64,
    pub r2: f64,
    /// Kit peptides used
    pub peptides: usize,
}

impl IrtCalibration {
    pub fn apply(&self, rt: f64) -> f64 {
        self.slope * rt + self.intercept
    }

    /// From the unmodified kit peptides among a library's precursors (their
    /// median time per sequence); `None` with fewer than 4 of them or a poor fit
    pub fn fit(precursors: &[LibraryPrecursor]) -> Result<Self, String> {
        let mut times: HashMap<&str, Vec<f64>> = HashMap::new();
        for p in precursors.iter().filter(|p| !p.decoy && p.peptides.len() == 1) {
            let peptide = &p.peptides[0].peptide;
            let unmodified = peptide.nterm.is_none()
                && peptide.cterm.is_none()
                && peptide.modifications.iter().all(|m| m.is_none());
            if !unmodified {
                continue;
            }
            if let Some((seq, _)) = IRT_KIT.iter().find(|(s, _)| s.as_bytes() == peptide.sequence.as_ref()) {
                times.entry(seq).or_default().push(p.rt);
            }
        }
        let points: Vec<(f64, f64)> = IRT_KIT
            .iter()
            .filter_map(|(seq, irt)| times.get_mut(seq).map(|t| (median(t), *irt)))
            .collect();
        let found = points.len();
        let (slope, intercept, r2) = robust_line(&points, 4)
            .ok_or_else(|| format!("{found} of the 11 iRT kit peptides were found (4 needed)"))?;
        if slope <= 0.0 || r2 < 0.9 {
            return Err(format!(
                "the {found} iRT kit peptides found do not fit a line (slope {slope:.3}, r2 {r2:.3})"
            ));
        }
        Ok(IrtCalibration { slope, intercept, r2, peptides: found })
    }
}

/// Per-run ion mobility calibration onto a reference run (linear)
#[derive(Clone, Debug, Default)]
pub struct ImCalibration {
    /// `(slope, intercept)` per file; identity where none was fitted
    pub fits: Vec<(f64, f64)>,
    pub reference: usize,
}

impl ImCalibration {
    pub fn apply(&self, file_id: usize, im: f32) -> f32 {
        match self.fits.get(file_id) {
            Some(&(slope, intercept)) => (slope * im as f64 + intercept) as f32,
            None => im,
        }
    }

    /// From `(precursor, file, im)` observations (each run's best PSM per
    /// precursor): the run with the most precursors is the reference, and every
    /// other run sharing at least 5 precursors with it is fitted to it
    fn fit<K: std::hash::Hash + Eq + Copy>(observations: &[(K, usize, f32)], n_files: usize) -> Self {
        let mut per_file: Vec<HashMap<K, f32>> = vec![HashMap::new(); n_files];
        for &(key, file, im) in observations {
            if let Some(map) = per_file.get_mut(file) {
                map.insert(key, im);
            }
        }
        let reference = (0..n_files).max_by_key(|&f| per_file[f].len()).unwrap_or(0);
        let fits = (0..n_files)
            .map(|file| {
                if file == reference {
                    return (1.0, 0.0);
                }
                let points: Vec<(f64, f64)> = per_file[file]
                    .iter()
                    .filter_map(|(key, &im)| per_file[reference].get(key).map(|&r| (im as f64, r as f64)))
                    .collect();
                match robust_line(&points, 5) {
                    Some((slope, intercept, _)) if slope > 0.5 && slope < 2.0 => (slope, intercept),
                    _ => (1.0, 0.0),
                }
            })
            .collect();
        ImCalibration { fits, reference }
    }
}

/// One run's best observation of a precursor
struct Observation {
    file_id: usize,
    rt: f64,
    im: f32,
    /// Lower is better
    rank_key: f64,
    /// Neutral mass of the precursor
    mass: f64,
    fragments: Vec<LibraryFragment>,
}

/// Fragments scaled to a maximum of 10,000, one per ion (the most intense)
fn normalized(mut fragments: Vec<LibraryFragment>) -> Vec<LibraryFragment> {
    fragments.sort_by(|a, b| {
        (a.peptide, a.kind as u8, a.number, a.charge)
            .cmp(&(b.peptide, b.kind as u8, b.number, b.charge))
            .then(b.intensity.total_cmp(&a.intensity))
    });
    fragments.dedup_by(|b, a| (a.peptide, a.kind, a.number, a.charge) == (b.peptide, b.kind, b.number, b.charge));
    let max = fragments.iter().map(|f| f.intensity).fold(0.0f32, f32::max);
    if max > 0.0 {
        for f in &mut fragments {
            f.intensity = f.intensity / max * 10_000.0;
        }
    }
    fragments.retain(|f| f.intensity > 0.0);
    fragments
}

/// The consensus (or best) entry of a precursor's observations
fn combine(mut observations: Vec<Observation>, settings: &LibrarySettings, im: &ImCalibration) -> Option<(f64, Option<f32>, Vec<LibraryFragment>, usize)> {
    observations.sort_by(|a, b| a.rank_key.total_cmp(&b.rank_key));
    let has_im = observations.iter().any(|o| o.im > 0.0);
    let calibrated = |o: &Observation| im.apply(o.file_id, o.im);
    let (rt, ion_mobility, fragments, runs) = if settings.consensus && observations.len() > 1 {
        let rt = median(&mut observations.iter().map(|o| o.rt).collect::<Vec<_>>());
        let ion_mobility = has_im.then(|| {
            median(&mut observations.iter().filter(|o| o.im > 0.0).map(|o| calibrated(o) as f64).collect::<Vec<_>>()) as f32
        });
        let mut per_ion: HashMap<(u8, u8, u16, u8), (LibraryFragment, Vec<f64>)> = HashMap::new();
        for o in &observations {
            for f in &o.fragments {
                per_ion
                    .entry((f.peptide, f.kind as u8, f.number, f.charge))
                    .or_insert_with(|| (*f, Vec::new()))
                    .1
                    .push(f.intensity as f64);
            }
        }
        let fragments: Vec<LibraryFragment> = per_ion
            .into_values()
            .map(|(mut f, mut intensities)| {
                f.intensity = median(&mut intensities) as f32;
                f
            })
            .collect();
        (rt, ion_mobility, normalized(fragments), observations.len())
    } else {
        let best = &observations[0];
        let ion_mobility = (best.im > 0.0).then(|| calibrated(best));
        (best.rt, ion_mobility, best.fragments.clone(), 1)
    };
    let mut fragments = fragments;
    fragments.sort_by(|a, b| b.intensity.total_cmp(&a.intensity));
    fragments.truncate(settings.max_fragments);
    if fragments.len() < settings.min_fragments {
        return None;
    }
    Some((rt, ion_mobility, fragments, runs))
}

/// Statistics of a library build
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LibraryStats {
    /// PSMs (CSMs) passing the filters
    pub matches: usize,
    pub precursors: usize,
    pub peptides: usize,
    pub fragments: usize,
    /// Precursors dropped for having fewer than `min_fragments`
    pub too_few_fragments: usize,
    /// No target protein reached `protein_fdr` (a small database cannot):
    /// the library is filtered at PSM and peptide level only
    pub protein_filter_skipped: bool,
}

/// The linear library: see the module documentation. `features` must carry
/// their q-values, aligned retention times (`aligned_rt`, fractions of the
/// run) and matched fragments.
pub fn linear_library(
    db: &IndexedDatabase,
    features: &[Feature],
    n_files: usize,
    settings: &LibrarySettings,
) -> (Vec<LibraryPrecursor>, ImCalibration, LibraryStats) {
    let mut stats = LibraryStats::default();
    // a database of a few proteins cannot reach a low protein FDR (picked
    // protein q-values start at 1 / proteins): the protein filter then goes
    let protein_filter = features
        .iter()
        .any(|f| f.rank == 1 && f.label == 1 && f.protein_q <= settings.protein_fdr);
    stats.protein_filter_skipped = !protein_filter;
    // each run's best PSM per precursor
    let mut best: HashMap<(PeptideIx, u8, usize), &Feature> = HashMap::new();
    for f in features {
        let passes = f.rank == 1
            && f.label == 1
            && f.spectrum_q <= settings.psm_fdr
            && f.peptide_q <= settings.peptide_fdr
            && (!protein_filter || f.protein_q <= settings.protein_fdr)
            && f.fragments.as_ref().map_or(false, |frags| !frags.kinds.is_empty());
        if !passes {
            continue;
        }
        let peptide = &db[f.peptide_idx];
        if peptide.decoy || (settings.proteotypic && accessions(peptide).0.len() != 1) {
            continue;
        }
        stats.matches += 1;
        let key = (f.peptide_idx, f.charge, f.file_id);
        let better = |old: &&Feature| {
            (f.posterior_error, -f.discriminant_score) < (old.posterior_error, -old.discriminant_score)
        };
        match best.get(&key) {
            Some(old) if !better(old) => {}
            _ => {
                best.insert(key, f);
            }
        }
    }

    let im_observations: Vec<((PeptideIx, u8), usize, f32)> = best
        .iter()
        .filter(|(_, f)| f.ims > 0.0)
        .map(|(&(ix, z, file), f)| ((ix, z), file, f.ims))
        .collect();
    let im = ImCalibration::fit(&im_observations, n_files.max(1));

    let mut per_precursor: HashMap<(PeptideIx, u8), Vec<Observation>> = HashMap::new();
    for (&(ix, z, file_id), f) in &best {
        let frags = f.fragments.as_ref().expect("filtered above");
        let fragments: Vec<LibraryFragment> = (0..frags.kinds.len())
            .map(|i| LibraryFragment {
                kind: frags.kinds[i],
                number: frags.fragment_ordinals[i].max(0) as u16,
                charge: frags.charges[i].clamp(1, u8::MAX as i32) as u8,
                mz: frags.mz_calculated[i] as f64,
                intensity: frags.intensities[i],
                peptide: 0,
                linked: false,
                shift: 0.0,
            })
            .collect();
        per_precursor.entry((ix, z)).or_default().push(Observation {
            file_id,
            rt: f.aligned_rt as f64 * 100.0,
            im: f.ims,
            rank_key: f.posterior_error as f64 - 1e-6 * f.discriminant_score as f64,
            mass: db[ix].monoisotopic as f64,
            fragments: normalized(fragments),
        });
    }

    let mut precursors: Vec<LibraryPrecursor> = Vec::with_capacity(per_precursor.len());
    let mut keys: Vec<(PeptideIx, u8)> = per_precursor.keys().copied().collect();
    keys.sort();
    for key in keys {
        let observations = per_precursor.remove(&key).unwrap_or_default();
        let peptide = db[key.0].clone();
        let Some((rt, ion_mobility, fragments, runs)) = combine(observations, settings, &im) else {
            stats.too_few_fragments += 1;
            continue;
        };
        precursors.push(LibraryPrecursor {
            mz: precursor_mz(peptide.monoisotopic as f64, key.1),
            peptides: vec![LibraryPeptide { peptide, links: Vec::new() }],
            charge: key.1,
            rt,
            ion_mobility,
            decoy: false,
            fragments,
            runs,
        });
    }
    finish_stats(&mut stats, &precursors);
    (precursors, im, stats)
}

fn finish_stats(stats: &mut LibraryStats, precursors: &[LibraryPrecursor]) {
    let targets: Vec<&LibraryPrecursor> = precursors.iter().filter(|p| !p.decoy).collect();
    stats.precursors = targets.len();
    stats.fragments = targets.iter().map(|p| p.fragments.len()).sum();
    let mut forms: Vec<String> = targets
        .iter()
        .map(|p| {
            p.peptides
                .iter()
                .map(|lp| format!("{:?}{:?}{:?}{:?}{:?}", lp.peptide.sequence, lp.peptide.modifications, lp.peptide.nterm, lp.peptide.cterm, lp.links))
                .collect::<Vec<_>>()
                .join("_")
        })
        .collect();
    forms.sort();
    forms.dedup();
    stats.peptides = forms.len();
}

/// Key of a crosslinked species for the library: peptides in library order
/// with their link sites, and the charge
type CrosslinkKey = (Vec<(PeptideIx, Vec<usize>)>, u8);

/// Library order of a crosslinked pair: alphabetical by sequence (then
/// modifications, then link site), as xiDIA-library orders it. `true` when
/// beta comes first.
fn swapped(db: &IndexedDatabase, m: &CrosslinkMatch) -> bool {
    match m.beta {
        None => false,
        Some(beta) => {
            let (a, b) = (&db[m.alpha], &db[beta]);
            match a.initial_sort(b) {
                std::cmp::Ordering::Greater => true,
                std::cmp::Ordering::Less => false,
                std::cmp::Ordering::Equal => m.beta_site < m.alpha_site,
            }
        }
    }
}

fn crosslink_key(db: &IndexedDatabase, m: &CrosslinkMatch) -> CrosslinkKey {
    let mut peptides = match m.beta {
        None => {
            let (s1, s2) = (m.alpha_site.min(m.beta_site) as usize, m.alpha_site.max(m.beta_site) as usize);
            vec![(m.alpha, vec![s1, s2])]
        }
        Some(beta) => vec![(m.alpha, vec![m.alpha_site as usize]), (beta, vec![m.beta_site as usize])],
    };
    if swapped(db, m) {
        peptides.swap(0, 1);
    }
    (peptides, m.charge)
}

/// The crosslink library (`xlibrary`): target pairs and loop-links at
/// `crosslink_fdr` (CSM and unique-crosslink q-values), each crosslink's
/// runs combined like the linear library's precursors. `rt` gives a match's
/// retention time on the library's 0-100 scale; `im` calibrates ion
/// mobilities (the linear library's, if there is one). Multivalent species
/// are left out.
pub fn crosslink_library(
    db: &IndexedDatabase,
    matches: &[CrosslinkMatch],
    rt: impl Fn(&CrosslinkMatch) -> f64,
    im: &ImCalibration,
    settings: &LibrarySettings,
) -> (Vec<LibraryPrecursor>, LibraryStats) {
    let mut stats = LibraryStats::default();
    let mut best: HashMap<(CrosslinkKey, usize), &CrosslinkMatch> = HashMap::new();
    for m in matches {
        let passes = m.rank == 1
            && m.class == PairClass::TT
            && m.kind != LinkKind::Multi
            && m.spectrum_q <= settings.crosslink_fdr
            && m.crosslink_q <= settings.crosslink_fdr
            && m.fragments.as_ref().map_or(false, |f| !f.is_empty());
        if !passes {
            continue;
        }
        stats.matches += 1;
        let key = (crosslink_key(db, m), m.file_id);
        match best.get(&key) {
            Some(old) if old.hyperscore >= m.hyperscore => {}
            _ => {
                best.insert(key, m);
            }
        }
    }
    let mut per_precursor: HashMap<CrosslinkKey, Vec<Observation>> = HashMap::new();
    for ((key, file_id), m) in best {
        let swap = swapped(db, m);
        let fragments: Vec<LibraryFragment> = m
            .fragments
            .as_ref()
            .expect("filtered above")
            .iter()
            .map(|f: &LinkedFragment| LibraryFragment {
                kind: f.kind,
                number: f.number,
                charge: f.charge,
                mz: f.mz as f64,
                intensity: f.intensity,
                peptide: if swap { 1 - f.peptide.min(1) } else { f.peptide },
                linked: f.linked,
                shift: f.shift,
            })
            .collect();
        per_precursor.entry(key).or_default().push(Observation {
            file_id,
            rt: rt(m),
            im: m.ims,
            rank_key: -m.hyperscore,
            mass: m.calcmass as f64,
            fragments: normalized(fragments),
        });
    }
    let mut keys: Vec<CrosslinkKey> = per_precursor.keys().cloned().collect();
    keys.sort();
    let mut precursors = Vec::with_capacity(keys.len());
    for key in keys {
        let observations = per_precursor.remove(&key).unwrap_or_default();
        let mass = observations.first().map(|o| o.mass).unwrap_or_default();
        let Some((rt, ion_mobility, fragments, runs)) = combine(observations, settings, im) else {
            stats.too_few_fragments += 1;
            continue;
        };
        let (peptides, charge) = key;
        let peptides: Vec<LibraryPeptide> = peptides
            .into_iter()
            .map(|(ix, links)| LibraryPeptide { peptide: db[ix].clone(), links })
            .collect();
        precursors.push(LibraryPrecursor {
            peptides,
            charge,
            mz: precursor_mz(mass, charge),
            rt,
            ion_mobility,
            decoy: false,
            fragments,
            runs,
        });
    }
    finish_stats(&mut stats, &precursors);
    (precursors, stats)
}

/// Where a residue moves when [`Peptide::reverse`] makes a decoy
fn reversed_position(site: usize, len: usize) -> usize {
    let n = len.saturating_sub(1);
    if n > 1 && site >= 1 && site < n {
        n - site
    } else {
        site
    }
}

/// A decoy of a peptide: the residues between the first and the last
/// reversed, as the search makes its decoys; where that gives the same
/// sequence (`FVKVR`) or a target's (`is_target`, I and L alike), rotated by
/// one instead; `None` when that is a target's too. Returns the decoy and
/// where each residue (and so each link site) moved.
fn decoy_of(peptide: &Peptide, is_target: &dyn Fn(&[u8]) -> bool) -> Option<(Peptide, Vec<usize>)> {
    let len = peptide.sequence.len();
    let moved = |map: &[usize]| -> Peptide {
        let mut sequence = vec![0u8; len];
        let mut modifications = vec![None; len];
        for (old, &new) in map.iter().enumerate() {
            sequence[new] = peptide.sequence[old];
            modifications[new] = peptide.modifications.get(old).copied().flatten();
        }
        Peptide {
            decoy: !peptide.decoy,
            sequence: std::sync::Arc::from(sequence.into_boxed_slice()),
            modifications: modifications.into_boxed_slice(),
            ..peptide.clone()
        }
    };
    let own = isobaric_key(&peptide.sequence);
    let usable = |p: &Peptide| {
        let key = isobaric_key(&p.sequence);
        key != own && !is_target(&key)
    };
    let reversed: Vec<usize> = (0..len).map(|i| reversed_position(i, len)).collect();
    let decoy = moved(&reversed);
    if usable(&decoy) {
        return Some((decoy, reversed));
    }
    if len < 4 {
        return None;
    }
    let n = len - 1;
    let rotated: Vec<usize> = (0..len)
        .map(|i| match i {
            0 => 0,
            _ if i == n => n,
            1 => n - 1,
            _ => i - 1,
        })
        .collect();
    let decoy = moved(&rotated);
    usable(&decoy).then_some((decoy, rotated))
}

/// Decoys for a crosslink library: both peptides reversed (as the search's
/// decoys are made; link sites move with their residues), with the targets'
/// fragment annotations, intensities and shifts and the masses of the
/// reversed peptides. The precursor m/z does not change. A decoy peptide
/// must not have the sequence of a target peptide (`is_target`, given the
/// [`isobaric_key`], I and L alike), or the decoy would take that peptide's
/// true hits: it is then rotated instead, and a crosslink whose peptides
/// cannot be made into decoys gets none. Returns the decoys and the number
/// of crosslinks without one.
pub fn crosslink_decoys(
    targets: &[LibraryPrecursor],
    mods: &Modifications,
    is_target: &dyn Fn(&[u8]) -> bool,
) -> (Vec<LibraryPrecursor>, usize) {
    let mut without = 0;
    let decoys = targets
        .iter()
        .filter(|p| !p.decoy)
        .filter_map(|target| {
            let peptides: Option<Vec<LibraryPeptide>> = target
                .peptides
                .iter()
                .map(|lp| {
                    let (peptide, moved) = decoy_of(&lp.peptide, is_target)?;
                    Some(LibraryPeptide { peptide, links: lp.links.iter().map(|&s| moved[s]).collect() })
                })
                .collect();
            let Some(peptides) = peptides else {
                without += 1;
                return None;
            };
            let fragments = target
                .fragments
                .iter()
                .filter_map(|f| {
                    let lp = peptides.get(f.peptide as usize)?;
                    let mass = ion_mass(&lp.peptide, mods, f.kind, f.number)?;
                    Some(LibraryFragment { mz: (mass + f.shift) as f64 / f.charge as f64 + PROTON as f64, ..*f })
                })
                .collect();
            Some(LibraryPrecursor { peptides, decoy: true, fragments, ..target.clone() })
        })
        .collect();
    (decoys, without)
}

/// Neutral mass of a backbone ion by series number
fn ion_mass(peptide: &Peptide, mods: &Modifications, kind: Kind, number: u16) -> Option<f32> {
    let len = peptide.sequence.len();
    let number = number as usize;
    if number == 0 || number >= len {
        return None;
    }
    let idx = match kind {
        Kind::A | Kind::B | Kind::C => number - 1,
        _ => len - 1 - number,
    };
    IonSeries::new(peptide, kind, mods).nth(idx).map(|ion| ion.monoisotopic_mass)
}

/// Retention times of a library onto iRT
pub fn apply_irt(precursors: &mut [LibraryPrecursor], irt: &IrtCalibration) {
    for p in precursors {
        p.rt = irt.apply(p.rt);
    }
}

// ---------------------------------------------------------------- modifications

/// Where a modification sits, for naming it
#[derive(Clone, Copy, Debug, PartialEq)]
enum ModSite {
    Residue(u8),
    Nterm,
    Cterm,
}

struct KnownMod {
    unimod: u32,
    mass: f32,
    /// Residues (or `^` for N-terminal, `$` for C-terminal)
    sites: &'static str,
    /// Spectronaut's name for it
    spectronaut: &'static str,
}

/// Common modifications, recognized by mass and site where a modification has
/// no `UNIMOD:` name
const KNOWN_MODS: &[KnownMod] = &[
    KnownMod { unimod: 4, mass: 57.021464, sites: "C", spectronaut: "Carbamidomethyl (C)" },
    KnownMod { unimod: 35, mass: 15.994915, sites: "MW", spectronaut: "Oxidation (M)" },
    KnownMod { unimod: 1, mass: 42.010565, sites: "^", spectronaut: "Acetyl (Protein N-term)" },
    KnownMod { unimod: 1, mass: 42.010565, sites: "KSTY", spectronaut: "Acetyl (K)" },
    KnownMod { unimod: 21, mass: 79.966331, sites: "STYH", spectronaut: "Phospho (STY)" },
    KnownMod { unimod: 7, mass: 0.984016, sites: "NQR", spectronaut: "Deamidation (NQ)" },
    KnownMod { unimod: 121, mass: 114.042927, sites: "KSTC", spectronaut: "GlyGly (K)" },
    KnownMod { unimod: 28, mass: -17.026549, sites: "^Q", spectronaut: "Gln->pyro-Glu" },
    KnownMod { unimod: 27, mass: -18.010565, sites: "^E", spectronaut: "Glu->pyro-Glu" },
    KnownMod { unimod: 26, mass: 39.994915, sites: "^C", spectronaut: "Pyro-carbamidomethyl" },
    KnownMod { unimod: 5, mass: 43.005814, sites: "^KR", spectronaut: "Carbamyl" },
    KnownMod { unimod: 34, mass: 14.01565, sites: "^KRDEHC", spectronaut: "Methyl" },
    KnownMod { unimod: 36, mass: 28.0313, sites: "^KR", spectronaut: "Dimethyl" },
    KnownMod { unimod: 199, mass: 32.056407, sites: "^K", spectronaut: "Dimethyl:2H(4)" },
    KnownMod { unimod: 510, mass: 34.063117, sites: "^K", spectronaut: "Dimethyl:2H(4)13C(2)" },
    KnownMod { unimod: 330, mass: 36.07567, sites: "^K", spectronaut: "Dimethyl:2H(6)13C(2)" },
    KnownMod { unimod: 259, mass: 8.014199, sites: "K", spectronaut: "Label:13C(6)15N(2)" },
    KnownMod { unimod: 267, mass: 10.008269, sites: "R", spectronaut: "Label:13C(6)15N(4)" },
    KnownMod { unimod: 188, mass: 6.020129, sites: "KR", spectronaut: "Label:13C(6)" },
    KnownMod { unimod: 481, mass: 4.025107, sites: "K", spectronaut: "Label:2H(4)" },
    KnownMod { unimod: 737, mass: 229.162932, sites: "^KST", spectronaut: "TMT6plex" },
    KnownMod { unimod: 2016, mass: 304.207146, sites: "^KST", spectronaut: "TMTpro" },
    KnownMod { unimod: 214, mass: 144.102063, sites: "^KY", spectronaut: "iTRAQ4plex" },
    KnownMod { unimod: 730, mass: 304.20536, sites: "^KY", spectronaut: "iTRAQ8plex" },
];

fn site_matches(sites: &str, site: ModSite) -> bool {
    match site {
        ModSite::Residue(r) => sites.as_bytes().contains(&r),
        ModSite::Nterm => sites.contains('^'),
        ModSite::Cterm => sites.contains('$'),
    }
}

fn known(m: &Modification, site: ModSite) -> Option<&'static KnownMod> {
    if let Some(id) = unimod_from_name(m) {
        return KNOWN_MODS
            .iter()
            .find(|k| k.unimod == id && site_matches(k.sites, site))
            .or_else(|| KNOWN_MODS.iter().find(|k| k.unimod == id));
    }
    KNOWN_MODS
        .iter()
        .find(|k| (k.mass - m.mass).abs() < 0.0015 && site_matches(k.sites, site))
}

fn unimod_from_name(m: &Modification) -> Option<u32> {
    let name = m.name.as_ref()?.to_string();
    let upper = name.to_ascii_uppercase();
    upper.strip_prefix("UNIMOD:").and_then(|n| n.parse().ok())
}

/// UniMod accession of a modification at a site: its `UNIMOD:` name, or a
/// common modification of that mass and site
fn unimod(m: &Modification, site: ModSite) -> Option<u32> {
    unimod_from_name(m).or_else(|| known(m, site).map(|k| k.unimod))
}

fn mass_tag(mass: f32) -> String {
    format!("{:+.6}", mass)
}

/// A modification as written in a layout's modified sequence
fn mod_tag(mods: &Modifications, id: ModId, site: ModSite, format: LibraryFormat) -> String {
    let m = mods.get(id);
    match format {
        LibraryFormat::Openswath => match unimod(m, site) {
            Some(n) => format!("(UniMod:{n})"),
            None => format!("[{}]", mass_tag(m.mass)),
        },
        LibraryFormat::Spectronaut => match known(m, site) {
            Some(k) => format!("[{}]", k.spectronaut),
            None => match unimod_from_name(m) {
                Some(n) => format!("[UniMod:{n}]"),
                None => format!("[{}]", mass_tag(m.mass)),
            },
        },
    }
}

/// Modified sequence of one peptide in a layout: `skip` leaves out label
/// modifications (Spectronaut's ModifiedPeptide), `linker` is written after
/// each linked residue
fn modified_sequence(
    lp: &LibraryPeptide,
    mods: &Modifications,
    format: LibraryFormat,
    skip: &[ModId],
    linker: Option<&str>,
) -> String {
    let p = &lp.peptide;
    let mut s = String::new();
    let keep = |id: &ModId| !skip.contains(id);
    if let Some(id) = p.nterm.filter(keep) {
        let tag = mod_tag(mods, id, ModSite::Nterm, format);
        match format {
            LibraryFormat::Openswath => {
                s.push('.');
                s.push_str(&tag);
            }
            LibraryFormat::Spectronaut => s.push_str(&tag),
        }
    }
    for (i, &residue) in p.sequence.iter().enumerate() {
        s.push(residue as char);
        if let Some(id) = p.modifications.get(i).copied().flatten().filter(keep) {
            s.push_str(&mod_tag(mods, id, ModSite::Residue(residue), format));
        }
        if let Some(linker) = linker {
            if lp.links.contains(&i) {
                match format {
                    LibraryFormat::Openswath => s.push_str(&format!("({linker})")),
                    LibraryFormat::Spectronaut => s.push_str(&format!("[{linker}]")),
                }
            }
        }
    }
    if let Some(id) = p.cterm.filter(keep) {
        let tag = mod_tag(mods, id, ModSite::Cterm, format);
        match format {
            LibraryFormat::Openswath => {
                s.push('.');
                s.push_str(&tag);
            }
            LibraryFormat::Spectronaut => s.push_str(&tag),
        }
    }
    s
}

fn stripped(p: &LibraryPrecursor) -> String {
    p.peptides
        .iter()
        .map(|lp| String::from_utf8_lossy(&lp.peptide.sequence).to_string())
        .collect::<Vec<_>>()
        .join("_")
}

fn joined_sequences(p: &LibraryPrecursor, mods: &Modifications, format: LibraryFormat, skip: &[ModId], linker: &str) -> String {
    let parts: Vec<String> = p
        .peptides
        .iter()
        .map(|lp| modified_sequence(lp, mods, format, skip, Some(linker)))
        .collect();
    parts.join("_")
}

fn kind_letter(kind: Kind) -> &'static str {
    match kind {
        Kind::A => "a",
        Kind::B => "b",
        Kind::C => "c",
        Kind::X => "x",
        Kind::Y => "y",
        Kind::Z | Kind::ZRadical => "z",
    }
}

/// DIA-NN's fragment types: `b` for a, b and c ions, `y` for x, y and z
fn diann_kind(kind: Kind) -> &'static str {
    match kind {
        Kind::A | Kind::B | Kind::C => "b",
        _ => "y",
    }
}

/// A number as written: to 1e-6, without trailing zeros
fn number(x: f64) -> String {
    format!("{}", (x * 1e6).round() / 1e6)
}

/// Proteins and entry names of a precursor: each peptide's accessions joined
/// with `;`, a crosslinked pair's two peptides with `_` (as xiDIA-library)
fn proteins(p: &LibraryPrecursor) -> (String, String, bool) {
    let per: Vec<(Vec<String>, Vec<String>)> = p.peptides.iter().map(|lp| accessions(&lp.peptide)).collect();
    let proteotypic = per.iter().all(|(a, _)| a.len() == 1);
    let accs = per.iter().map(|(a, _)| a.join(";")).collect::<Vec<_>>().join("_");
    let genes = per.iter().map(|(_, g)| g.join(";")).collect::<Vec<_>>().join("_");
    (accs, genes, proteotypic)
}

/// Write a library in the OpenSWATH layout
pub fn write_openswath(precursors: &[LibraryPrecursor], mods: &Modifications, linker: &str) -> Vec<u8> {
    let has_im = precursors.iter().any(|p| p.ion_mobility.is_some());
    let crosslinks = precursors.iter().any(|p| p.is_crosslink());
    let mut header = vec![
        "PrecursorMz",
        "ProductMz",
        "Annotation",
        "ProteinId",
        "GeneName",
        "PeptideSequence",
        "ModifiedPeptideSequence",
        "PrecursorCharge",
        "LibraryIntensity",
        "NormalizedRetentionTime",
    ];
    if has_im {
        header.push("PrecursorIonMobility");
    }
    header.extend(["FragmentType", "FragmentSeriesNumber", "FragmentCharge", "FragmentLossType", "Proteotypic", "Decoy"]);
    if crosslinks {
        header.extend(["FragmentPepId", "CLContainingFragment"]);
    }
    let mut out = String::new();
    out.push_str(&header.join("\t"));
    out.push('\n');
    for p in precursors {
        let modified = joined_sequences(p, mods, LibraryFormat::Openswath, &[], linker);
        let (accs, genes, proteotypic) = proteins(p);
        let seq = stripped(p);
        let precursor_mz = number(p.mz);
        let rt = number(p.rt);
        let im = p.ion_mobility.map(|im| number(im as f64)).unwrap_or_default();
        for f in &p.fragments {
            let mut row: Vec<String> = vec![
                precursor_mz.clone(),
                number(f.mz),
                format!("{}{}^{}", kind_letter(f.kind), f.number, f.charge),
                accs.clone(),
                genes.clone(),
                seq.clone(),
                modified.clone(),
                p.charge.to_string(),
                number(f.intensity as f64),
                rt.clone(),
            ];
            if has_im {
                row.push(im.clone());
            }
            row.extend([
                diann_kind(f.kind).to_string(),
                f.number.to_string(),
                f.charge.to_string(),
                "noloss".to_string(),
                u8::from(proteotypic).to_string(),
                u8::from(p.decoy).to_string(),
            ]);
            if crosslinks {
                row.extend([f.peptide.to_string(), u8::from(f.linked).to_string()]);
            }
            out.push_str(&row.join("\t"));
            out.push('\n');
        }
    }
    out.into_bytes()
}

/// Write a library in Spectronaut's layout (targets only: Spectronaut makes
/// its own decoys). `labels` are isotope-label modifications, left out of
/// `ModifiedPeptide` and kept in `LabeledPeptide`.
pub fn write_spectronaut(precursors: &[LibraryPrecursor], mods: &Modifications, labels: &[ModId], linker: &str) -> Vec<u8> {
    let has_im = precursors.iter().any(|p| p.ion_mobility.is_some());
    let crosslinks = precursors.iter().any(|p| p.is_crosslink());
    let mut header = vec!["PrecursorCharge", "ModifiedPeptide", "StrippedPeptide", "iRT"];
    if has_im {
        header.push("IonMobility");
    }
    header.extend([
        "IsProteotypic",
        "LabeledPeptide",
        "PrecursorMz",
        "FragmentLossType",
        "FragmentNumber",
        "FragmentType",
        "FragmentCharge",
        "FragmentMz",
        "RelativeIntensity",
        "ExcludeFromAssay",
        "ProteinGroups",
        "UniProtIds",
        "Genes",
    ]);
    if crosslinks {
        header.extend(["FragmentPepId", "CLContainingFragment"]);
    }
    let mut out = String::new();
    out.push_str(&header.join("\t"));
    out.push('\n');
    for p in precursors.iter().filter(|p| !p.decoy) {
        let (modified, labeled) = if p.is_crosslink() {
            // xiDIA-library's form, without Spectronaut's `_` ends
            (
                joined_sequences(p, mods, LibraryFormat::Spectronaut, labels, linker),
                joined_sequences(p, mods, LibraryFormat::Spectronaut, &[], linker),
            )
        } else {
            let lp = &p.peptides[0];
            (
                format!("_{}_", modified_sequence(lp, mods, LibraryFormat::Spectronaut, labels, None)),
                format!("_{}_", modified_sequence(lp, mods, LibraryFormat::Spectronaut, &[], None)),
            )
        };
        let (accs, genes, proteotypic) = proteins(p);
        let max = p.fragments.iter().map(|f| f.intensity).fold(0.0f32, f32::max).max(f32::MIN_POSITIVE);
        for f in &p.fragments {
            let mut row: Vec<String> = vec![p.charge.to_string(), modified.clone(), stripped(p), number(p.rt)];
            if has_im {
                row.push(p.ion_mobility.map(|im| number(im as f64)).unwrap_or_default());
            }
            row.extend([
                if proteotypic { "True" } else { "False" }.to_string(),
                labeled.clone(),
                number(p.mz),
                "noloss".to_string(),
                f.number.to_string(),
                kind_letter(f.kind).to_string(),
                f.charge.to_string(),
                number(f.mz),
                number((f.intensity / max * 100.0) as f64),
                "False".to_string(),
                accs.clone(),
                accs.clone(),
                genes.clone(),
            ]);
            if crosslinks {
                row.extend([f.peptide.to_string(), if f.linked { "True" } else { "False" }.to_string()]);
            }
            out.push_str(&row.join("\t"));
            out.push('\n');
        }
    }
    out.into_bytes()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::Builder;
    use crate::scoring::Fragments;

    fn database(variable: serde_json::Value, fasta: &str) -> IndexedDatabase {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 0, "min_len": 5, "cleave_at": "KR", "restrict": "P"},
            "static_mods": {"C": 57.021464},
            "variable_mods": variable,
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        params.build(crate::fasta::Fasta::parse(fasta.into(), "rev_", false))
    }

    fn find(db: &IndexedDatabase, seq: &str, modified: bool) -> PeptideIx {
        let i = db
            .peptides
            .iter()
            .position(|p| !p.decoy && p.sequence.as_ref() == seq.as_bytes() && p.modifications.iter().any(|m| m.is_some()) == modified)
            .unwrap_or_else(|| panic!("{seq} not in the database"));
        PeptideIx(i as u32)
    }

    fn psm(ix: PeptideIx, file_id: usize, rt: f32, pep: f32, intensities: &[f32]) -> Feature {
        let n = intensities.len();
        Feature {
            peptide_idx: ix,
            file_id,
            rank: 1,
            label: 1,
            charge: 2,
            aligned_rt: rt,
            posterior_error: pep,
            fragments: Some(Fragments {
                charges: vec![1; n],
                kinds: vec![Kind::Y; n],
                fragment_ordinals: (1..=n as i32).collect(),
                intensities: intensities.to_vec(),
                mz_calculated: (0..n).map(|i| 200.0 + i as f32).collect(),
                mz_experimental: vec![0.0; n],
            }),
            ..Default::default()
        }
    }

    /// UniMod accessions and Spectronaut names from names or masses
    #[test]
    fn modification_tags() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 0, "min_len": 5},
            "static_mods": {"C": 57.021464},
            "variable_mods": {"M": [15.994915], "S": [{"mass": 79.966331, "name": "UNIMOD:21"}], "^": [42.010565], "K": [176.01433]},
            "max_variable_mods": 4,
        }))
        .unwrap();
        let db = builder
            .make_parameters()
            .unwrap()
            .build(crate::fasta::Fasta::parse(">sp|P1|ONE_HUMAN\nMSACKDEFGR\n".into(), "rev_", false));
        let mods = &db.modifications;
        let peptide = db
            .peptides
            .iter()
            .find(|p| {
                !p.decoy
                    && p.sequence.as_ref() == b"MSACK"
                    && p.nterm.is_some()
                    && [0, 1, 3, 4].iter().all(|&i| p.modifications[i].is_some())
            })
            .expect("MSACK with all its modifications")
            .clone();
        let lp = LibraryPeptide { peptide, links: Vec::new() };
        let openswath = modified_sequence(&lp, mods, LibraryFormat::Openswath, &[], None);
        assert_eq!(openswath, ".(UniMod:1)M(UniMod:35)S(UniMod:21)AC(UniMod:4)K[+176.014328]");
        let spectronaut = modified_sequence(&lp, mods, LibraryFormat::Spectronaut, &[], None);
        assert_eq!(
            spectronaut,
            "[Acetyl (Protein N-term)]M[Oxidation (M)]S[Phospho (STY)]AC[Carbamidomethyl (C)]K[+176.014328]"
        );
        assert_eq!(accession_and_entry("sp|P02769|ALBU_BOVIN"), ("P02769".into(), "ALBU_BOVIN".into()));
        assert_eq!(accession_and_entry("human:sp|P05023|AT1A1_HUMAN"), ("P05023".into(), "AT1A1_HUMAN".into()));
        assert_eq!(accession_and_entry("ENTRAP_42"), ("ENTRAP_42".into(), String::new()));
    }

    /// EasyPQP's consensus: per-run best PSMs, median time and intensity,
    /// proteotypic peptides, fragment limits
    #[test]
    fn consensus_library() {
        let db = database(serde_json::json!({}), ">sp|P1|ONE_HUMAN\nLESLIEKAPEPTIDEKSAMPLER\n>sp|P2|TWO_HUMAN\nSAMPLER\n");
        let unique = find(&db, "APEPTIDEK", false);
        let shared = find(&db, "SAMPLER", false);
        let features = vec![
            // run 0: two PSMs, the better (lower PEP) counts
            psm(unique, 0, 0.30, 0.001, &[100.0, 50.0, 25.0, 5.0]),
            psm(unique, 0, 0.90, 0.2, &[1.0, 100.0, 1.0, 1.0]),
            // run 1 and run 2
            psm(unique, 1, 0.32, 0.01, &[100.0, 30.0, 25.0, 0.0]),
            psm(unique, 2, 0.34, 0.01, &[100.0, 70.0, 25.0, 5.0]),
            // in two proteins: not proteotypic
            psm(shared, 0, 0.5, 0.001, &[100.0, 50.0, 25.0, 5.0]),
        ];
        let settings = LibrarySettings::default();
        let (library, _, stats) = linear_library(&db, &features, 3, &settings);
        assert_eq!(library.len(), 1, "{library:?}");
        let p = &library[0];
        assert_eq!(p.runs, 3);
        assert!((p.rt - 32.0).abs() < 1e-4, "median of 30, 32, 34: {}", p.rt);
        let y2 = p.fragments.iter().find(|f| f.number == 2).unwrap();
        assert!((y2.intensity - 5000.0).abs() < 1e-2, "median of 5000, 3000, 7000: {}", y2.intensity);
        let y4 = p.fragments.iter().find(|f| f.number == 4).unwrap();
        assert!((y4.intensity - 500.0).abs() < 1e-2, "absent in run 1: median of the two");
        assert_eq!(p.fragments[0].intensity, 10_000.0, "most intense first");
        assert_eq!(stats.matches, 4);
        assert!((p.mz - precursor_mz(db[unique].monoisotopic as f64, 2)).abs() < 1e-9);

        // shared peptides kept without the proteotypic filter, and fragment limits
        let open = LibrarySettings { proteotypic: false, min_fragments: 1, max_fragments: 2, ..LibrarySettings::default() };
        let (library, _, _) = linear_library(&db, &features, 3, &open);
        assert_eq!(library.len(), 2);
        assert!(library.iter().all(|p| p.fragments.len() == 2));
        let strict = LibrarySettings { min_fragments: 5, ..LibrarySettings::default() };
        let (library, _, stats) = linear_library(&db, &features, 3, &strict);
        assert!(library.is_empty());
        assert_eq!(stats.too_few_fragments, 1);
        // best replicate instead of consensus
        let best = LibrarySettings { consensus: false, ..LibrarySettings::default() };
        let (library, _, _) = linear_library(&db, &features, 3, &best);
        assert!((library[0].rt - 30.0).abs() < 1e-4);
    }

    /// Ion mobility: other runs calibrated onto the run with most precursors
    #[test]
    fn ion_mobility_calibration() {
        // run 1 reads 2% high
        let observations: Vec<(usize, usize, f32)> = (0..20)
            .flat_map(|k| {
                let im = 0.7 + k as f32 * 0.03;
                [(k, 0, im), (k, 1, im * 1.02)]
            })
            .chain([(99, 0, 1.0)])
            .collect();
        let calibration = ImCalibration::fit(&observations, 2);
        assert_eq!(calibration.reference, 0);
        assert!((calibration.apply(1, 1.02) - 1.0).abs() < 1e-4);
        assert_eq!(calibration.apply(0, 0.9), 0.9);
    }

    /// iRT from the kit peptides: a line through their times
    #[test]
    fn irt_calibration() {
        let fasta: String = IRT_KIT.iter().enumerate().map(|(i, (s, _))| format!(">irt{i}\n{s}\n")).collect();
        let db = database(serde_json::json!({}), &fasta);
        let precursors: Vec<LibraryPrecursor> = IRT_KIT
            .iter()
            .take(6)
            .map(|(seq, irt)| LibraryPrecursor {
                peptides: vec![LibraryPeptide { peptide: db[find(&db, seq, false)].clone(), links: Vec::new() }],
                charge: 2,
                mz: 500.0,
                // the library's own 0-100 scale: rt = (irt + 30) / 2
                rt: (irt + 30.0) / 2.0,
                ion_mobility: None,
                decoy: false,
                fragments: Vec::new(),
                runs: 1,
            })
            .collect();
        let fit = IrtCalibration::fit(&precursors).unwrap();
        assert_eq!(fit.peptides, 6);
        assert!((fit.apply(50.0) - 70.0).abs() < 1e-6, "{fit:?}");
        assert!(IrtCalibration::fit(&precursors[..3]).is_err());
    }

    /// The layouts' columns, and crosslinks written as xiDIA-library writes them
    #[test]
    fn layouts() {
        let db = database(serde_json::json!({}), ">sp|P1|ONE_HUMAN\nPEPKPTIDEKLIKPER\n");
        let a = db[find(&db, "PEPKPTIDEK", false)].clone();
        let b = db[find(&db, "LIKPER", false)].clone();
        let fragment = |peptide: u8, number: u16, linked: bool| LibraryFragment {
            kind: Kind::Y,
            number,
            charge: 1,
            mz: 300.5,
            intensity: 5000.0,
            peptide,
            linked,
            shift: if linked { 100.0 } else { 0.0 },
        };
        let linear = LibraryPrecursor {
            peptides: vec![LibraryPeptide { peptide: a.clone(), links: Vec::new() }],
            charge: 2,
            mz: 512.25,
            rt: 33.5,
            ion_mobility: Some(0.95),
            decoy: false,
            fragments: vec![LibraryFragment { intensity: 10_000.0, ..fragment(0, 3, false) }, fragment(0, 2, false)],
            runs: 2,
        };
        let tsv = String::from_utf8(write_openswath(std::slice::from_ref(&linear), &db.modifications, "XL")).unwrap();
        let mut lines = tsv.lines();
        assert_eq!(
            lines.next().unwrap(),
            "PrecursorMz\tProductMz\tAnnotation\tProteinId\tGeneName\tPeptideSequence\tModifiedPeptideSequence\tPrecursorCharge\tLibraryIntensity\tNormalizedRetentionTime\tPrecursorIonMobility\tFragmentType\tFragmentSeriesNumber\tFragmentCharge\tFragmentLossType\tProteotypic\tDecoy"
        );
        assert_eq!(lines.next().unwrap(), "512.25\t300.5\ty3^1\tP1\tONE_HUMAN\tPEPKPTIDEK\tPEPKPTIDEK\t2\t10000\t33.5\t0.95\ty\t3\t1\tnoloss\t1\t0");
        let xls = String::from_utf8(write_spectronaut(std::slice::from_ref(&linear), &db.modifications, &[], "XL")).unwrap();
        let mut lines = xls.lines();
        assert_eq!(
            lines.next().unwrap(),
            "PrecursorCharge\tModifiedPeptide\tStrippedPeptide\tiRT\tIonMobility\tIsProteotypic\tLabeledPeptide\tPrecursorMz\tFragmentLossType\tFragmentNumber\tFragmentType\tFragmentCharge\tFragmentMz\tRelativeIntensity\tExcludeFromAssay\tProteinGroups\tUniProtIds\tGenes"
        );
        let second = lines.nth(1).unwrap();
        assert!(second.starts_with("2\t_PEPKPTIDEK_\tPEPKPTIDEK\t33.5\t"), "{second}");
        assert!(second.contains("\t50\tFalse\tP1\tP1\tONE_HUMAN"), "relative intensity in percent: {second}");

        // a crosslinked pair
        let pair = LibraryPrecursor {
            peptides: vec![
                LibraryPeptide { peptide: b.clone(), links: vec![2] },
                LibraryPeptide { peptide: a.clone(), links: vec![3] },
            ],
            ion_mobility: None,
            fragments: vec![fragment(0, 2, true), fragment(1, 2, false)],
            ..linear.clone()
        };
        let none = |_: &[u8]| false;
        let (decoys, without) = crosslink_decoys(std::slice::from_ref(&pair), &db.modifications, &none);
        assert_eq!((decoys.len(), without), (1, 0));
        let decoy = &decoys[0];
        assert!(decoy.decoy);
        assert_eq!(decoy.mz, pair.mz, "reversing keeps the precursor mass");
        assert_eq!(&*decoy.peptides[1].peptide.sequence, b"PEDITPKPEK".as_slice());
        assert_eq!(decoy.peptides[1].links, vec![6], "the linked K moves with its residue");
        assert_eq!(&*decoy.peptides[0].peptide.sequence, b"LEPKIR".as_slice());
        assert_eq!(decoy.peptides[0].links, vec![3]);
        // y2 of the reversed LIKPER, plus the shift
        let y2 = ion_mass(&decoy.peptides[0].peptide, &db.modifications, Kind::Y, 2).unwrap();
        assert!((decoy.fragments[0].mz - ((y2 + 100.0) as f64 + PROTON as f64)).abs() < 1e-3);
        let all = [pair.clone(), decoy.clone()];
        let tsv = String::from_utf8(write_openswath(&all, &db.modifications, "DSSO")).unwrap();
        let row = tsv.lines().nth(1).unwrap();
        assert!(row.contains("\tLIKPER_PEPKPTIDEK\tLIK(DSSO)PER_PEPK(DSSO)PTIDEK\t"), "{row}");
        assert!(row.ends_with("\t1\t0\t0\t1"), "proteotypic, target, first peptide, linked: {row}");
        assert_eq!(tsv.lines().filter(|l| l.contains("\t1\t1\t")).count(), 2, "the decoy's rows");
        let xls = String::from_utf8(write_spectronaut(&all, &db.modifications, &[], "DSSO")).unwrap();
        assert_eq!(xls.lines().count(), 3, "header and the target's two fragments, no decoy");

        // a peptide that reversing leaves as it is gets its inner residues rotated
        let palindrome = Peptide {
            sequence: std::sync::Arc::from(b"FVKVR".to_vec().into_boxed_slice()),
            modifications: vec![None; 5].into_boxed_slice(),
            ..b.clone()
        };
        let (decoy, moved) = decoy_of(&palindrome, &none).unwrap();
        assert_eq!(&*decoy.sequence, b"FKVVR".as_slice());
        assert_eq!(moved[2], 1, "the K");
        assert!(decoy.decoy);

        // a decoy peptide with a target's sequence (I and L alike) is not used:
        // rotated instead, and with no usable form the crosslink gets no decoy
        let leptkpek = isobaric_key(b"PEDITPKPEK");
        let target_reversed = |key: &[u8]| key == leptkpek.as_slice();
        let (decoys, without) = crosslink_decoys(std::slice::from_ref(&pair), &db.modifications, &target_reversed);
        assert_eq!((decoys.len(), without), (1, 0));
        assert_eq!(&*decoys[0].peptides[1].peptide.sequence, b"PPKPTIDEEK".as_slice(), "rotated");
        let any = |_: &[u8]| true;
        let (decoys, without) = crosslink_decoys(std::slice::from_ref(&pair), &db.modifications, &any);
        assert_eq!((decoys.len(), without), (0, 1));
        let row = xls.lines().nth(1).unwrap();
        assert!(row.starts_with("2\tLIK[DSSO]PER_PEPK[DSSO]PTIDEK\tLIKPER_PEPKPTIDEK\t"), "{row}");
        assert!(row.ends_with("\tP1_P1\tP1_P1\tONE_HUMAN_ONE_HUMAN\t0\tTrue"), "{row}");
    }

    #[test]
    fn settings() {
        let s: LibrarySettings = serde_json::from_value(serde_json::json!({})).unwrap();
        assert_eq!(s, LibrarySettings::default());
        assert!(s.validate().is_ok());
        let s: LibrarySettings = serde_json::from_value(serde_json::json!({"formats": ["spectronaut"], "irt": true})).unwrap();
        assert!(s.writes(LibraryFormat::Spectronaut) && !s.writes(LibraryFormat::Openswath) && s.irt);
        assert!(serde_json::from_value::<LibrarySettings>(serde_json::json!({"format": []})).is_err());
        let bad = LibrarySettings { min_fragments: 6, max_fragments: 4, ..LibrarySettings::default() };
        assert!(bad.validate().is_err());
        let bad = LibrarySettings { psm_fdr: 0.0, ..LibrarySettings::default() };
        assert!(bad.validate().is_err());
    }
}
