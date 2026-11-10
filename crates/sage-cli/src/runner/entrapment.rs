//! Entrapment estimates (see `sage_core::entrapment`): `entrapment.tsv`,
//! the log, and the `qc` section of `results.json`.

use super::Runner;
use log::info;
use sage_core::crosslink::{CrosslinkMatch, PairClass};
use sage_core::entrapment::{estimate, is_entrapment, pair_ratio, EntrapmentEstimate};
use sage_core::scoring::Feature;
use std::collections::HashSet;

/// Reported FDRs the estimates are made at
const LEVELS: [f32; 2] = [0.01, 0.05];

impl Runner {
    fn entrapment_tag(&self) -> Option<String> {
        self.parameters.entrapment.as_ref().map(|e| e.tag.clone()).filter(|t| !t.is_empty())
    }

    /// Entrapment-to-original ratios of the database: target peptide forms
    /// and the proteins they come from
    fn entrapment_ratios(&self, tag: &str) -> (f64, f64) {
        let (mut pep_e, mut pep_t) = (0usize, 0usize);
        let mut proteins: HashSet<&str> = HashSet::new();
        for p in self.database.peptides.iter().filter(|p| !p.decoy) {
            if is_entrapment(&p.proteins, tag) {
                pep_e += 1;
            } else {
                pep_t += 1;
            }
            proteins.extend(p.proteins.iter().map(|a| a.as_ref()));
        }
        let prot_e = proteins.iter().filter(|a| a.contains(tag)).count();
        let prot_t = proteins.len() - prot_e;
        let ratio = |e: usize, t: usize| if t > 0 { e as f64 / t as f64 } else { 0.0 };
        (ratio(pep_e, pep_t), ratio(prot_e, prot_t))
    }

    /// PSMs, peptides and proteins of the linear search
    pub(super) fn entrapment_linear(&self, features: &[Feature]) -> Vec<EntrapmentEstimate> {
        let Some(tag) = self.entrapment_tag() else { return Vec::new() };
        let (r_peptides, r_proteins) = self.entrapment_ratios(&tag);
        let db = &self.database;
        let targets: Vec<&Feature> = features.iter().filter(|f| f.rank == 1 && f.label == 1).collect();
        let mut out = Vec::new();
        for fdr in LEVELS {
            let count = |items: &mut dyn Iterator<Item = bool>| {
                items.fold((0, 0), |(t, e), entrapment| if entrapment { (t, e + 1) } else { (t + 1, e) })
            };
            let (t, e) = count(
                &mut targets.iter().filter(|f| f.spectrum_q <= fdr).map(|f| is_entrapment(&db[f.peptide_idx].proteins, &tag)),
            );
            out.push(estimate("PSMs", fdr, t, e, r_peptides));
            let peptides: HashSet<_> =
                targets.iter().filter(|f| f.peptide_q <= fdr).map(|f| f.peptide_idx).collect();
            let (t, e) = count(&mut peptides.iter().map(|&ix| is_entrapment(&db[ix].proteins, &tag)));
            out.push(estimate("peptides", fdr, t, e, r_peptides));
            // proteins with protein_q: those of peptides unique to one protein
            let proteins: HashSet<&str> = targets
                .iter()
                .filter(|f| f.protein_q <= fdr && db[f.peptide_idx].proteins.len() == 1)
                .map(|f| db[f.peptide_idx].proteins[0].as_ref())
                .collect();
            let (t, e) = count(&mut proteins.iter().map(|a| a.contains(tag.as_str())));
            out.push(estimate("proteins", fdr, t, e, r_proteins));
        }
        self.log_entrapment(&out, &tag);
        out
    }

    /// Crosslink-spectrum matches: a crosslink is entrapment when either
    /// peptide maps only to entrapment proteins. The combined estimate uses
    /// the random-pairing ratio (an approximation for crosslinks).
    pub(super) fn entrapment_crosslinks(&self, matches: &[CrosslinkMatch]) -> Vec<EntrapmentEstimate> {
        let Some(tag) = self.entrapment_tag() else { return Vec::new() };
        let (r_peptides, _) = self.entrapment_ratios(&tag);
        let db = &self.database;
        let mut out = Vec::new();
        for fdr in LEVELS {
            let (mut t, mut e) = (0, 0);
            for m in matches.iter().filter(|m| m.class == PairClass::TT && m.spectrum_q <= fdr) {
                let entrapment = std::iter::once(m.alpha).chain(m.beta).any(|ix| is_entrapment(&db[ix].proteins, &tag));
                if entrapment {
                    e += 1;
                } else {
                    t += 1;
                }
            }
            out.push(estimate("crosslink-spectrum matches", fdr, t, e, pair_ratio(r_peptides)));
        }
        self.log_entrapment(&out, &tag);
        out
    }

    fn log_entrapment(&self, estimates: &[EntrapmentEstimate], tag: &str) {
        for fdr in LEVELS {
            let parts: Vec<String> = estimates
                .iter()
                .filter(|e| e.fdr == fdr)
                .map(|e| {
                    format!(
                        "{} {:.2}% (at least {:.2}%; {} of {})",
                        e.level,
                        100.0 * e.combined,
                        100.0 * e.lower_bound,
                        e.entrapments,
                        e.targets + e.entrapments
                    )
                })
                .collect();
            if !parts.is_empty() {
                info!("- entrapment ({tag}) at a reported {:.0}% FDR: {}", 100.0 * fdr, parts.join(", "));
            }
        }
    }

    /// `entrapment.tsv`: every estimate so far (written again when crosslinks add theirs)
    pub(super) fn write_entrapment(&mut self, estimates: &[EntrapmentEstimate]) -> anyhow::Result<()> {
        if estimates.is_empty() {
            return Ok(());
        }
        let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
        wtr.write_record(["level", "reported_fdr", "targets", "entrapments", "ratio", "lower_bound", "combined_fdp"])?;
        for e in estimates {
            wtr.write_record([
                e.level.clone(),
                e.fdr.to_string(),
                e.targets.to_string(),
                e.entrapments.to_string(),
                format!("{:.4}", e.ratio),
                format!("{:.5}", e.lower_bound),
                format!("{:.5}", e.combined),
            ])?;
        }
        let path = self.make_path("entrapment.tsv");
        sage_cloudpath::write_bytes_sync(&path, wtr.into_inner()?)?;
        if !self.parameters.output_paths.contains(&path) {
            self.parameters.output_paths.push(path);
        }
        Ok(())
    }
}
