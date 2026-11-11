//! Spectral libraries of the search results (`spectral_library`, see
//! `sage_core::spectral_library`): `library.tsv` (OpenSWATH layout: DIA-NN,
//! PEAKS, OpenSWATH) and `library.spectronaut.xls` (Spectronaut layout:
//! Spectronaut, DIA-NN, PEAKS) from the linear search; `xlibrary.tsv` and
//! `xlibrary.spectronaut.xls` from a crosslink search.

use super::Runner;
use log::{info, warn};
use sage_core::crosslink::{CrosslinkMatch, CrosslinkSettings, LinkKind, PairClass};
use sage_core::ml::retention_alignment::Alignment;
use sage_core::scoring::Feature;
use sage_core::spectral_library::{
    apply_irt, crosslink_decoys, crosslink_library, linear_library, write_openswath, write_spectronaut, ImCalibration,
    IrtCalibration, LibraryFormat, LibrarySettings, LibraryStats,
};
use std::collections::HashMap;

/// What the linear library leaves for the crosslink library of the same
/// search: the runs' alignment, the iRT calibration and the ion mobility
/// calibration, so both libraries share one scale
#[derive(Default)]
pub struct LibraryContext {
    pub alignments: Option<Vec<Alignment>>,
    pub irt: Option<IrtCalibration>,
    pub im: ImCalibration,
}

fn summary(stats: &LibraryStats, what: &str, q: f32) -> String {
    format!(
        "{} precursors ({} {what}, {} fragments) from {} matches at {}% FDR",
        stats.precursors,
        stats.peptides,
        stats.fragments,
        stats.matches,
        q * 100.0
    )
}

/// The crosslinker's name in a library's modified sequences
fn linker_tag(settings: &CrosslinkSettings) -> String {
    match (settings.reagent, settings.mass) {
        (Some(reagent), _) => format!("{reagent:?}").to_uppercase(),
        (None, Some(_)) => "XL".to_string(),
        (None, None) => "DSS".to_string(),
    }
}

impl Runner {
    /// The linear library, from features with q-values, aligned retention
    /// times and matched fragments
    pub(super) fn write_linear_library(
        &mut self,
        features: &[Feature],
        settings: &LibrarySettings,
        alignments: Vec<Alignment>,
    ) -> anyhow::Result<()> {
        let n_files = self.parameters.mzml_paths.len();
        let (mut precursors, im, stats) = linear_library(&self.database, features, n_files, settings);
        let irt = match settings.irt {
            false => None,
            true => match IrtCalibration::fit(&precursors) {
                Ok(fit) => {
                    info!(
                        "spectral library: iRT from {} kit peptides, iRT = {:.4} x + {:.3} (r2 {:.3})",
                        fit.peptides, fit.slope, fit.intercept, fit.r2
                    );
                    apply_irt(&mut precursors, &fit);
                    Some(fit)
                }
                Err(why) => {
                    warn!("spectral library: no iRT calibration ({why}); retention times stay on the 0-100 scale of the aligned runs");
                    None
                }
            },
        };
        if stats.protein_filter_skipped {
            warn!(
                "spectral library: no protein reaches {}% protein FDR (a database of few proteins cannot), so the library is filtered at PSM and peptide level only",
                settings.protein_fdr * 100.0
            );
        }
        let mods = &self.database.modifications;
        let labels = self.database.isotope_mods();
        if settings.writes(LibraryFormat::Openswath) {
            let path = self.make_path("library.tsv");
            sage_cloudpath::write_bytes_sync(&path, write_openswath(&precursors, mods, ""))?;
            self.parameters.output_paths.push(path);
        }
        if settings.writes(LibraryFormat::Spectronaut) {
            let path = self.make_path("library.spectronaut.xls");
            sage_cloudpath::write_bytes_sync(&path, write_spectronaut(&precursors, mods, &labels, ""))?;
            self.parameters.output_paths.push(path);
        }
        let with_im = precursors.iter().filter(|p| p.ion_mobility.is_some()).count();
        info!(
            "spectral library: {}{}{}",
            summary(&stats, "peptides", settings.psm_fdr),
            if with_im > 0 { format!(", {with_im} with ion mobility") } else { String::new() },
            if stats.too_few_fragments > 0 {
                format!("; {} precursors with fewer than {} fragments left out", stats.too_few_fragments, settings.min_fragments)
            } else {
                String::new()
            }
        );
        self.library = Some(LibraryContext { alignments: Some(alignments), irt, im });
        Ok(())
    }

    /// The crosslink library, from CSMs with q-values and matched fragments
    pub(super) fn write_crosslink_library(
        &mut self,
        matches: &[CrosslinkMatch],
        crosslink: &CrosslinkSettings,
        settings: &LibrarySettings,
    ) -> anyhow::Result<()> {
        let n_files = self.parameters.mzml_paths.len();
        let context = self.library.take().unwrap_or_default();
        let alignments = match context.alignments {
            Some(alignments) => alignments,
            None => {
                // no linear search: align the runs on the crosslinks themselves
                let mut species: HashMap<(u32, u16, u32, u16), Vec<usize>> = HashMap::new();
                for (i, m) in matches.iter().enumerate() {
                    if m.rank == 1
                        && m.class == PairClass::TT
                        && m.kind != LinkKind::Multi
                        && m.spectrum_q <= settings.crosslink_fdr
                    {
                        let beta = m.beta.map_or(u32::MAX, |b| b.0);
                        species.entry((m.alpha.0, m.alpha_site, beta, m.beta_site)).or_default().push(i);
                    }
                }
                let mut species: Vec<Vec<usize>> = species.into_values().collect();
                for s in &mut species {
                    s.sort_by(|&a, &b| matches[b].hyperscore.total_cmp(&matches[a].hyperscore));
                }
                sage_core::crosslink_quant::align_runs(
                    matches,
                    &species,
                    &[],
                    n_files,
                    self.parameters.quant.lfq_settings.alignment,
                )
            }
        };
        let rt = |m: &CrosslinkMatch| alignments.get(m.file_id).map_or(0.0, |a| a.transform(m.rt) as f64 * 100.0);
        let (mut targets, stats) = crosslink_library(&self.database, matches, rt, &context.im, settings);
        match (&context.irt, settings.irt) {
            (Some(irt), _) => apply_irt(&mut targets, irt),
            (None, true) => warn!(
                "crosslink library: no iRT calibration (it comes from the iRT kit peptides of the linear search); retention times stay on the 0-100 scale of the aligned runs"
            ),
            (None, false) => {}
        }
        let mods = &self.database.modifications;
        let tag = linker_tag(crosslink);
        if settings.writes(LibraryFormat::Openswath) {
            let mut all = targets.clone();
            if settings.crosslink_decoys {
                // no decoy peptide may be a target peptide of the search (I and L alike)
                let target_keys: std::collections::HashSet<Vec<u8>> = self
                    .database
                    .peptides
                    .iter()
                    .filter(|p| !p.decoy)
                    .map(|p| sage_core::database::isobaric_key(&p.sequence))
                    .collect();
                let is_target = |key: &[u8]| target_keys.contains(key);
                let (decoys, without) = crosslink_decoys(&targets, mods, &is_target);
                if without > 0 {
                    info!("crosslink library: {without} crosslinks have no decoy (their peptides reversed or rotated are target peptides)");
                }
                all.extend(decoys);
            }
            let path = self.make_path("xlibrary.tsv");
            sage_cloudpath::write_bytes_sync(&path, write_openswath(&all, mods, &tag))?;
            self.parameters.output_paths.push(path);
        }
        if settings.writes(LibraryFormat::Spectronaut) {
            let labels = self.database.isotope_mods();
            let path = self.make_path("xlibrary.spectronaut.xls");
            sage_cloudpath::write_bytes_sync(&path, write_spectronaut(&targets, mods, &labels, &tag))?;
            self.parameters.output_paths.push(path);
        }
        let loops = targets.iter().filter(|p| p.peptides.len() == 1).count();
        let multi = matches
            .iter()
            .filter(|m| m.kind == LinkKind::Multi && m.class == PairClass::TT && m.spectrum_q <= settings.crosslink_fdr)
            .count();
        info!(
            "crosslink library: {} ({loops} of them loop-links){}{}",
            summary(&stats, "crosslinks", settings.crosslink_fdr),
            if settings.crosslink_decoys && settings.writes(LibraryFormat::Openswath) {
                ", with one decoy each in xlibrary.tsv"
            } else {
                ""
            },
            if multi > 0 { format!("; {multi} multivalent CSMs are not in it") } else { String::new() }
        );
        Ok(())
    }
}
