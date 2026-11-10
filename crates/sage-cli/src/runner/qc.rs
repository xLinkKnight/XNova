//! Quality control of a search: polymer contamination in MS1 and diagnostic
//! ions in MS2 (collected while the spectra are read), and the digestion
//! summary of the peptides identified at 1% FDR. Written to `digestion.tsv`,
//! `diagnostic_ions.tsv` (when diagnostic ions are searched) and the `qc`
//! section of `results.json`, and summarised in the log.

use super::Runner;
use log::{info, warn};
use rayon::prelude::*;
use sage_core::diagnostic::{find_ions, DiagnosticIon};
use sage_core::digestion::{summarize, Classifier, DigestionRow, Identified};
use sage_core::enzyme::EnzymeParameters;
use sage_core::polymer::{PolymerFileStats, WARN_PCT};
use sage_core::scoring::Feature;
use sage_core::spectrum::{RawSpectrum, Representation};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap, HashSet};

/// What the spectra showed, gathered while they are read
#[derive(Default)]
pub struct QcCollector {
    /// Files already scanned (a crosslink search reads them again)
    scanned: HashSet<usize>,
    polymers: BTreeMap<usize, PolymerFileStats>,
    /// (file, spectrum id, ion, m/z, relative intensity)
    hits: Vec<(usize, String, usize, f32, f32)>,
    /// MS2 spectra scanned for diagnostic ions, per file
    ms2: BTreeMap<usize, usize>,
}

/// The `qc` section of `results.json`
#[derive(Serialize, Default, Clone, Debug)]
pub struct QcReport {
    /// One row per file, then `total`
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub digestion: Vec<DigestionRow>,
    /// File name → MS1 polymer signal
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub polymers: BTreeMap<String, PolymerFileStats>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub diagnostic_ions: Vec<DiagnosticSummary>,
    /// False discoveries at the reported FDRs, from entrapment proteins
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub entrapment: Vec<sage_core::entrapment::EntrapmentEstimate>,
}

/// How often one diagnostic ion was seen in one file's MS2 spectra
#[derive(Serialize, Clone, Debug)]
pub struct DiagnosticSummary {
    pub file: String,
    pub ion: String,
    pub spectra: usize,
    pub pct: f64,
}

/// `(mz, intensity)` sorted by m/z (most readers' are already)
fn sorted(spectrum: &RawSpectrum) -> std::borrow::Cow<'_, RawSpectrum> {
    if spectrum.mz.windows(2).all(|w| w[0] <= w[1]) {
        return std::borrow::Cow::Borrowed(spectrum);
    }
    let mut order: Vec<usize> = (0..spectrum.mz.len()).collect();
    order.sort_by(|&a, &b| spectrum.mz[a].total_cmp(&spectrum.mz[b]));
    let mut copy = spectrum.clone();
    copy.mz = order.iter().map(|&i| spectrum.mz[i]).collect();
    copy.intensity = order.iter().map(|&i| spectrum.intensity[i]).collect();
    std::borrow::Cow::Owned(copy)
}

impl QcCollector {
    /// Scan the raw spectra of files not scanned yet: MS1 for polymer
    /// ladders, MS2 for `ions` (before deisotoping and peak trimming)
    pub fn collect(&mut self, ms1: &[RawSpectrum], msn: &[RawSpectrum], ions: Option<&[DiagnosticIon]>) {
        let new: HashSet<usize> =
            ms1.iter().chain(msn).map(|s| s.file_id).filter(|f| !self.scanned.contains(f)).collect();
        if new.is_empty() {
            return;
        }
        let per_file = ms1
            .par_iter()
            .filter(|s| new.contains(&s.file_id))
            .fold(BTreeMap::new, |mut acc: BTreeMap<usize, PolymerFileStats>, s| {
                let stats = acc.entry(s.file_id).or_insert_with(PolymerFileStats::new);
                match s.representation {
                    Representation::Centroid => {
                        let s = sorted(s);
                        stats.add(&s.mz, &s.intensity)
                    }
                    Representation::Profile => stats.skipped_profile_spectra += 1,
                }
                acc
            })
            .reduce(BTreeMap::new, |mut a, b| {
                for (file, stats) in b {
                    a.entry(file).or_insert_with(PolymerFileStats::new).merge(&stats);
                }
                a
            });
        for (file, stats) in per_file {
            self.polymers.entry(file).or_insert_with(PolymerFileStats::new).merge(&stats);
        }
        if let Some(ions) = ions {
            let hits: Vec<_> = msn
                .par_iter()
                .filter(|s| new.contains(&s.file_id))
                .flat_map_iter(|s| {
                    let peaks = sorted(s);
                    find_ions(ions, &peaks.mz, &peaks.intensity)
                        .into_iter()
                        .map(|h| (s.file_id, s.id.clone(), h.ion, h.mz, h.relative_intensity))
                        .collect::<Vec<_>>()
                })
                .collect();
            self.hits.extend(hits);
            for s in msn.iter().filter(|s| new.contains(&s.file_id)) {
                *self.ms2.entry(s.file_id).or_default() += 1;
            }
        }
        self.scanned.extend(new);
    }

    /// The data gathered so far, keeping the record of scanned files
    fn take(&mut self) -> (BTreeMap<usize, PolymerFileStats>, Vec<(usize, String, usize, f32, f32)>, BTreeMap<usize, usize>) {
        (std::mem::take(&mut self.polymers), std::mem::take(&mut self.hits), std::mem::take(&mut self.ms2))
    }
}

/// Distinct sequences, as `(sequence, proteins, decoy)` of rank-1 PSMs at 1%
/// spectrum and peptide FDR
fn accepted(features: &[Feature]) -> impl Iterator<Item = &Feature> {
    features.iter().filter(|f| f.rank == 1 && f.spectrum_q <= 0.01 && f.peptide_q <= 0.01)
}

impl Runner {
    /// Turn what was gathered into `results.json`'s `qc` section, the log
    /// summary and the QC tables. `features`: the linear search's PSMs (none
    /// for a crosslink-only search, which has no digestion summary).
    pub(super) fn finish_qc(&mut self, features: Option<&[Feature]>, filenames: &[String]) -> anyhow::Result<()> {
        let (polymers, mut hits, ms2) = match self.qc.lock() {
            Ok(mut qc) => qc.take(),
            Err(_) => return Ok(()),
        };
        let name = |file: usize| filenames.get(file).cloned().unwrap_or_else(|| format!("file {file}"));
        let mut report = QcReport::default();

        for (file, mut stats) in polymers {
            stats.finish();
            if stats.ms1_spectra == 0 {
                continue;
            }
            let shares: Vec<String> = stats.polymers.iter().map(|p| format!("{} {:.2}%", p.name, p.tic_pct)).collect();
            info!("- MS1 polymer signal in {}: {} of the MS1 ion current ({} spectra)", name(file), shares.join(", "), stats.ms1_spectra);
            for p in &stats.polymers {
                if p.tic_pct > WARN_PCT {
                    warn!(
                        "{}: {} ladders carry {:.1}% of the MS1 ion current: polymer contamination (detergent or plastic)",
                        name(file),
                        p.name,
                        p.tic_pct
                    );
                }
            }
            report.polymers.insert(name(file), stats);
        }

        if let Some(ions) = self.parameters.diagnostic_ions.clone() {
            hits.sort_by(|a, b| (a.0, &a.1, a.2).cmp(&(b.0, &b.1, b.2)));
            let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
            wtr.write_record(["file", "scannr", "ion", "mz", "relative_intensity"])?;
            for (file, id, ion, mz, rel) in &hits {
                wtr.write_record([name(*file), id.clone(), ions[*ion].name.clone(), mz.to_string(), rel.to_string()])?;
            }
            let path = self.make_path("diagnostic_ions.tsv");
            sage_cloudpath::write_bytes_sync(&path, wtr.into_inner()?)?;
            self.parameters.output_paths.push(path);
            for (&file, &spectra) in &ms2 {
                let mut line = Vec::new();
                for (i, ion) in ions.iter().enumerate() {
                    let n = hits.iter().filter(|h| h.0 == file && h.2 == i).count();
                    let pct = if spectra > 0 { 100.0 * n as f64 / spectra as f64 } else { 0.0 };
                    line.push(format!("{} {:.1}%", ion.name, pct));
                    report.diagnostic_ions.push(DiagnosticSummary { file: name(file), ion: ion.name.clone(), spectra: n, pct });
                }
                info!("- diagnostic ions in {} (share of {} MS2 spectra): {}", name(file), spectra, line.join(", "));
            }
        }

        if let Some(features) = features {
            report.entrapment = self.entrapment_linear(features);
            self.write_entrapment(&report.entrapment.clone())?;
            report.digestion = self.digestion_summary(features, filenames)?;
            if let Some(total) = report.digestion.last() {
                info!(
                    "- digestion of {} distinct peptide sequences at 1% FDR (decoy-corrected): {:.1}% with missed cleavages, {:.1}% semi-enzymatic (N {:.1}%, C {:.1}%), {:.1}% non-enzymatic",
                    total.peptides,
                    total.missed_cleavage_pct,
                    total.semi_n_pct + total.semi_c_pct,
                    total.semi_n_pct,
                    total.semi_c_pct,
                    total.non_enzymatic_pct
                );
            }
        }
        self.parameters.qc = Some(report);
        Ok(())
    }

    /// `digestion.tsv`: one row per file and a `total` row (each sequence once)
    fn digestion_summary(&mut self, features: &[Feature], filenames: &[String]) -> anyhow::Result<Vec<DigestionRow>> {
        let wanted: HashSet<&str> = accepted(features)
            .flat_map(|f| self.database[f.peptide_idx].proteins.iter().map(|p| p.as_ref()))
            .collect();
        let fasta = super::read_protein_databases(&mut self.parameters.database.clone())?;
        let proteins: HashMap<String, String> = fasta
            .targets
            .iter()
            .filter(|(acc, _)| wanted.contains(acc.as_ref()))
            .map(|(acc, seq)| (acc.to_string(), seq.clone()))
            .collect();
        let enzyme: EnzymeParameters = self.parameters.database.enzyme.clone().into();
        let classifier = Classifier::new(&enzyme, &proteins, self.parameters.database.clip_n_term_met);
        let db = &self.database;
        let identified = |f: &Feature| {
            let p = &db[f.peptide_idx];
            Identified { sequence: std::str::from_utf8(&p.sequence).unwrap_or(""), proteins: &p.proteins, decoy: p.decoy }
        };
        let mut rows: Vec<DigestionRow> = (0..filenames.len())
            .map(|file| summarize(&filenames[file], accepted(features).filter(|f| f.file_id == file).map(identified), &classifier))
            .collect();
        rows.push(summarize("total", accepted(features).map(identified), &classifier));

        let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
        wtr.write_record([
            "file", "target_peptides", "decoy_peptides", "peptides", "missed_cleavages_0", "missed_cleavages_1",
            "missed_cleavages_2_plus", "semi_n", "semi_c", "non_enzymatic", "missed_cleavage_pct", "semi_n_pct",
            "semi_c_pct", "non_enzymatic_pct",
        ])?;
        let pct = |x: f64| format!("{x:.2}");
        for r in &rows {
            wtr.write_record([
                r.file.clone(),
                r.target_peptides.to_string(),
                r.decoy_peptides.to_string(),
                r.peptides.to_string(),
                r.missed_cleavages_0.to_string(),
                r.missed_cleavages_1.to_string(),
                r.missed_cleavages_2_plus.to_string(),
                r.semi_n.to_string(),
                r.semi_c.to_string(),
                r.non_enzymatic.to_string(),
                pct(r.missed_cleavage_pct),
                pct(r.semi_n_pct),
                pct(r.semi_c_pct),
                pct(r.non_enzymatic_pct),
            ])?;
        }
        let path = self.make_path("digestion.tsv");
        sage_cloudpath::write_bytes_sync(&path, wtr.into_inner()?)?;
        self.parameters.output_paths.push(path);
        Ok(rows)
    }
}
