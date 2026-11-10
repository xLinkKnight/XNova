//! `sites.tsv`: modification sites on proteins from localized PSMs (see
//! `sage_core::localization`), with a site-level target-decoy q-value.
//!
//! A PSM at 1% FDR supports a site when it places the modification there
//! with probability 0.75 or more; a site's score is the best discriminant
//! score of its supporting PSMs, and sites of decoy PSMs give the q-value,
//! `(decoys + 1) / targets` ranked by score. Every target site that a PSM
//! places with probability 0.5 or more is listed.

use super::Runner;
use log::info;
use sage_cloudpath::Url;
use sage_core::database::PeptideIx;
use sage_core::scoring::Feature;
use std::collections::{HashMap, HashSet};

const PSM_Q: f32 = 0.01;
const SUPPORT: f32 = 0.75;
const LISTED: f32 = 0.5;

#[derive(Default)]
struct Site {
    probability: f32,
    psms: usize,
    peptides: HashSet<PeptideIx>,
    localization_q: f32,
    score: f32,
    q: f32,
}

impl Runner {
    pub(super) fn write_sites(&mut self, features: &[Feature]) -> anyhow::Result<Option<Url>> {
        let fasta = super::read_protein_databases(&mut self.parameters.database.clone())?;
        let proteins: HashMap<&str, &[u8]> = fasta.targets.iter().map(|(a, s)| (a.as_ref(), s.as_bytes())).collect();
        let db = &self.database;
        // (protein, 1-based position, residue, modification) -> site; decoys by peptide
        let mut targets: HashMap<(String, usize, char, String), Site> = HashMap::new();
        let mut decoys: HashMap<(String, usize, String), f32> = HashMap::new();
        let (mut localized, mut confident) = (0usize, 0usize);
        for f in features.iter().filter(|f| f.rank == 1 && f.spectrum_q <= PSM_Q) {
            let Some(loc) = &f.localization else { continue };
            let peptide = &db[f.peptide_idx];
            if f.label == 1 {
                localized += 1;
                confident += usize::from(loc.q <= 0.01);
            }
            for group in &loc.groups {
                let name = db.modifications.get(group.modification).to_string();
                for &(i, p) in &group.sites {
                    let i = i as usize;
                    if f.label != 1 {
                        if p >= SUPPORT {
                            let key = (String::from_utf8_lossy(&peptide.sequence).to_string(), i, name.clone());
                            let best = decoys.entry(key).or_insert(f32::NEG_INFINITY);
                            *best = best.max(f.discriminant_score);
                        }
                        continue;
                    }
                    if p < LISTED {
                        continue;
                    }
                    for accession in &peptide.proteins {
                        let Some(protein) = proteins.get(accession.as_ref()) else { continue };
                        let Some(start) = protein.windows(peptide.sequence.len()).position(|w| w == peptide.sequence.as_ref())
                        else {
                            continue;
                        };
                        let key = (accession.to_string(), start + i + 1, peptide.sequence[i] as char, name.clone());
                        let site = targets.entry(key).or_insert_with(|| Site {
                            localization_q: 1.0,
                            score: f32::NEG_INFINITY,
                            ..Default::default()
                        });
                        site.probability = site.probability.max(p);
                        site.localization_q = site.localization_q.min(group.q);
                        if p >= SUPPORT {
                            site.psms += 1;
                            site.peptides.insert(f.peptide_idx);
                            site.score = site.score.max(f.discriminant_score);
                        }
                    }
                }
            }
        }
        if targets.is_empty() {
            info!("- localization: no localized modification sites at {:.0}% FDR", 100.0 * PSM_Q);
            return Ok(None);
        }

        // site q-values: supported target sites against decoy sites, by score
        let mut ranked: Vec<(f32, bool, Option<(String, usize, char, String)>)> = targets
            .iter()
            .filter(|(_, s)| s.psms > 0)
            .map(|(k, s)| (s.score, false, Some(k.clone())))
            .chain(decoys.values().map(|&score| (score, true, None)))
            .collect();
        ranked.sort_by(|a, b| b.0.total_cmp(&a.0));
        let (mut t, mut d) = (0.0f32, 1.0f32);
        let mut qs: Vec<f32> = Vec::with_capacity(ranked.len());
        for r in &ranked {
            if r.1 {
                d += 1.0;
            } else {
                t += 1.0;
            }
            qs.push(if t > 0.0 { (d / t).min(1.0) } else { 1.0 });
        }
        let mut minimum = 1.0f32;
        for (r, q) in ranked.iter().zip(qs.iter()).rev() {
            minimum = minimum.min(*q);
            if let Some(key) = &r.2 {
                if let Some(site) = targets.get_mut(key) {
                    site.q = minimum;
                }
            }
        }
        for site in targets.values_mut().filter(|s| s.psms == 0) {
            site.q = 1.0;
        }

        let mut rows: Vec<_> = targets.into_iter().collect();
        rows.sort_by(|a, b| a.0 .0.cmp(&b.0 .0).then(a.0 .1.cmp(&b.0 .1)).then(a.0 .3.cmp(&b.0 .3)));
        let passing = rows.iter().filter(|(_, s)| s.q <= 0.01).count();
        let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
        wtr.write_record([
            "protein",
            "position",
            "residue",
            "modification",
            "probability",
            "psms",
            "peptides",
            "localization_q",
            "site_q",
        ])?;
        for ((protein, position, residue, modification), site) in &rows {
            wtr.write_record([
                protein.clone(),
                position.to_string(),
                residue.to_string(),
                modification.clone(),
                format!("{:.4}", site.probability),
                site.psms.to_string(),
                site.peptides.len().to_string(),
                ryu::Buffer::new().format(site.localization_q).to_string(),
                ryu::Buffer::new().format(site.q).to_string(),
            ])?;
        }
        let path = self.make_path("sites.tsv");
        sage_cloudpath::write_bytes_sync(&path, wtr.into_inner()?)?;
        info!(
            "- localization: {localized} target PSMs with localizable modifications at {:.0}% FDR, {confident} at 1% FLR; \
             {} sites in sites.tsv, {passing} at 1% site FDR",
            100.0 * PSM_Q,
            rows.len()
        );
        Ok(Some(path))
    }
}
