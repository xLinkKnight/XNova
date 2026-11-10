//! Modification site localization (adapted from Sage Plus, MIT).
//!
//! For the best PSM of each spectrum, every variable modification on a
//! residue is moved over the residues that could carry it (its residues,
//! within its position rule and motif) and over *decoy residues* that cannot
//! (by default alanine): each arrangement is rescored with the search's own
//! fragment matching, and the arrangements' hyperscores give probabilities
//! (`exp(h_i - h_best)`, normalized). A site's probability sums the
//! arrangements that place the modification there; reported probabilities
//! are over the real (target) arrangements, decoy residues only serve the
//! false localization rate.
//!
//! The false localization rate (FLR) comes from the decoy residues, for each
//! modification on its own (Ramsbottom et al. 2022's decoy amino acid
//! approach, per PTM): a PSM whose best arrangement of a modification puts it
//! on a decoy residue is a known false localization of that modification.
//! Ranked by the probability of their best arrangement among all (decoy
//! residues included), the FLR at each PSM is `D * (n_t / n_d) / T`
//! (decoy-residue and target winners above it, corrected by the numbers of
//! candidate target and decoy positions), made monotone: the modification's
//! localization q-value. A PSM's localization q-value is that of its least
//! certain modification. Pooled, a modification that is often misplaced on
//! a decoy residue (Met oxidation next to alanine, say) would raise the
//! q-values of all the others.
//!
//! `localized_peptide` is written in ProForma 2.0 with localization scores
//! (`T[#g1(0.09)]S[UNIMOD:21#g1(0.90)]`) where a site is not certain.

use crate::database::ResolvedModifications;
use crate::modification::{ModId, ModificationSpecificity, Modifications};
use crate::peptide::Peptide;
use crate::scoring::{Feature, Scorer};
use crate::site_rules::SiteRule;
use crate::spectrum::ProcessedSpectrum;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

/// `localization`: on with `{}`
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub struct LocalizationSettings {
    /// Residues no modification may sit on, used as decoy sites (default `A`)
    #[serde(default = "default_decoy_residues")]
    pub decoy_residues: String,
    /// Arrangements scored per modification at most (more: not localized)
    #[serde(default = "default_max_arrangements")]
    pub max_arrangements: usize,
}

fn default_decoy_residues() -> String {
    "A".into()
}

fn default_max_arrangements() -> usize {
    256
}

impl Default for LocalizationSettings {
    fn default() -> Self {
        Self { decoy_residues: default_decoy_residues(), max_arrangements: default_max_arrangements() }
    }
}

/// One localized modification type of a PSM
#[derive(Clone, Debug, PartialEq)]
pub struct SiteGroup {
    pub modification: ModId,
    /// Target residue positions (in the peptide) and their probabilities,
    /// most likely first
    pub sites: Vec<(u16, f32)>,
    /// Positions of the best target arrangement
    pub best: Vec<u16>,
    /// Probability of the best arrangement among all, decoy residues
    /// included: the FLR ranks by it
    pub score: f32,
    /// Hyperscore of the best target arrangement less the next best (any)
    pub delta: f32,
    /// Arrangements scored (1: nothing to localize)
    pub arrangements: u32,
    /// The best arrangement puts the modification on a decoy residue
    pub decoy_site: bool,
    /// Candidate target and decoy positions, for the FLR correction
    pub target_positions: u32,
    pub decoy_positions: u32,
    /// Localization q-value (FLR) of this modification, set by
    /// [`assign_localization_q`]
    pub q: f32,
}

/// Localization of a PSM's variable modifications
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Localization {
    pub groups: Vec<SiteGroup>,
    /// Probability of the best target arrangement among the target
    /// arrangements (product over groups)
    pub probability: f32,
    /// Arrangements scored (product over groups; 1: nothing to localize)
    pub arrangements: u32,
    /// Localization q-value: the highest of the groups' (its least certain
    /// modification), set by [`assign_localization_q`]
    pub q: f32,
}

/// A localizable modification: where it may sit, and its decoy residues
#[derive(Clone, Debug)]
struct Localizable {
    id: ModId,
    /// Residue and the rule on it
    residues: Vec<(u8, SiteRule)>,
    decoys: Vec<u8>,
}

/// Localizes the variable modifications of a database's peptides
#[derive(Clone, Debug)]
pub struct Localizer {
    mods: Vec<Localizable>,
    max_arrangements: usize,
}

impl Localizer {
    /// From the search's configured modifications (`Parameters::modifications`:
    /// their ids are the database's)
    pub fn new(mods: &ResolvedModifications, settings: &LocalizationSettings) -> Localizer {
        let statics: Vec<u8> = mods
            .static_mods
            .iter()
            .filter_map(|&(s, _)| match s {
                ModificationSpecificity::Residue(r) => Some(r),
                _ => None,
            })
            .collect();
        let mut out: Vec<Localizable> = Vec::new();
        // indexed variable modifications, then mass offsets; a site library
        // limits where the search places a modification, not where its
        // localization may move it
        let rules = mods.variable_mods.iter().enumerate().map(|(i, &(specificity, id, _))| {
            (specificity, id, mods.variable_rules.get(i).cloned().unwrap_or_default())
        });
        let offsets = mods.offset_mods.iter().map(|(specificity, id, rule)| (*specificity, *id, rule.clone()));
        for (specificity, id, rule) in rules.chain(offsets) {
            let ModificationSpecificity::Residue(r) = specificity else { continue };
            let rule = rule.without_library();
            match out.iter_mut().find(|l| l.id == id) {
                Some(l) => l.residues.push((r, rule)),
                None => out.push(Localizable { id, residues: vec![(r, rule)], decoys: Vec::new() }),
            }
        }
        for l in &mut out {
            l.decoys = settings
                .decoy_residues
                .bytes()
                .map(|b| b.to_ascii_uppercase())
                .filter(|b| !l.residues.iter().any(|(r, _)| r == b) && !statics.contains(b))
                .unique()
                .collect();
        }
        Localizer { mods: out, max_arrangements: settings.max_arrangements.max(2) }
    }

    /// Localize the best PSM of `query` (rank 1)
    pub fn localize(&self, scorer: &Scorer, query: &ProcessedSpectrum, features: &mut [Feature]) {
        let Some(feature) = features.iter_mut().find(|f| f.rank == 1) else { return };
        let peptide = scorer.db.resolve(feature);
        let charge = feature.charge;
        feature.localization = self.localize_peptide(&peptide, |p| scorer.rescore(query, p, charge));
    }

    /// Localize `peptide`'s modifications with `score` (hyperscore of an
    /// arrangement)
    pub fn localize_peptide(&self, peptide: &Peptide, mut score: impl FnMut(&Peptide) -> f64) -> Option<Localization> {
        let mut loc = Localization { probability: 1.0, arrangements: 1, ..Default::default() };
        let mut current = peptide.clone();
        let sequence = peptide.sequence.clone();
        for m in &self.mods {
            let carrying: Vec<usize> = (0..sequence.len()).filter(|&i| current.modifications[i] == Some(m.id)).collect();
            let k = carrying.len();
            if k == 0 {
                continue;
            }
            let free = |i: usize| current.modifications[i].is_none() || current.modifications[i] == Some(m.id);
            let targets: Vec<usize> = (0..sequence.len())
                .filter(|&i| free(i))
                .filter(|&i| m.residues.iter().any(|(r, rule)| *r == sequence[i] && rule.allows(i, &sequence, None)))
                .collect();
            // the configured position must stay a candidate (a motif may have
            // used protein residues the peptide alone lacks)
            let targets: Vec<usize> = targets.into_iter().chain(carrying.iter().copied()).sorted().dedup().collect();
            let decoys: Vec<usize> = (0..sequence.len())
                .filter(|&i| current.modifications[i].is_none() && m.decoys.contains(&sequence[i]))
                .collect();
            let candidates: Vec<usize> = targets.iter().chain(&decoys).copied().sorted().collect();
            let n = candidates.len();
            let combinations = binomial(n, k);
            if combinations > self.max_arrangements as u64 {
                return None;
            }
            let mut scored: Vec<(Vec<usize>, bool, f64)> = Vec::with_capacity(combinations as usize);
            for combination in candidates.iter().copied().combinations(k) {
                let mut p = current.clone();
                for &i in &carrying {
                    p.modifications[i] = None;
                }
                for &i in &combination {
                    p.modifications[i] = Some(m.id);
                }
                let decoy = combination.iter().any(|i| decoys.contains(i));
                let h = if combination == carrying { score(&current) } else { score(&p) };
                scored.push((combination, decoy, h));
            }
            loc.arrangements = loc.arrangements.saturating_mul(scored.len() as u32);
            let best_any = scored.iter().map(|s| s.2).fold(f64::NEG_INFINITY, f64::max);
            let weights: Vec<f64> = scored.iter().map(|s| (s.2 - best_any).exp()).collect();
            let total: f64 = weights.iter().sum();
            let target_total: f64 = scored.iter().zip(&weights).filter(|(s, _)| !s.1).map(|(_, w)| w).sum();
            let (best_target, best_h) = scored
                .iter()
                .enumerate()
                .filter(|(_, s)| !s.1)
                .max_by(|a, b| a.1 .2.total_cmp(&b.1 .2).then(b.0.cmp(&a.0)))
                .map(|(i, s)| (i, s.2))?;
            let next = scored
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != best_target)
                .map(|(_, s)| s.2)
                .fold(f64::NEG_INFINITY, f64::max);
            loc.probability *= (weights[best_target] / target_total) as f32;
            // the best arrangement has weight 1
            let best_score = (1.0 / total) as f32;
            let delta = if next.is_finite() { (best_h - next) as f32 } else { 0.0 };
            let decoy_site = scored.iter().any(|s| s.1 && s.2 > best_h);
            let mut sites: Vec<(u16, f32)> = targets
                .iter()
                .map(|&i| {
                    let p: f64 =
                        scored.iter().zip(&weights).filter(|(s, _)| !s.1 && s.0.contains(&i)).map(|(_, w)| w).sum();
                    (i as u16, (p / target_total) as f32)
                })
                .collect();
            sites.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
            let best: Vec<usize> = scored[best_target].0.clone();
            for &i in &carrying {
                current.modifications[i] = None;
            }
            for &i in &best {
                current.modifications[i] = Some(m.id);
            }
            loc.groups.push(SiteGroup {
                modification: m.id,
                sites,
                best: best.iter().map(|&i| i as u16).collect(),
                score: best_score,
                delta,
                arrangements: scored.len() as u32,
                decoy_site,
                target_positions: targets.len() as u32,
                decoy_positions: decoys.len() as u32,
                q: 1.0,
            });
        }
        if loc.groups.is_empty() {
            return None;
        }
        Some(loc)
    }
}

fn binomial(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    (0..k).fold(1u64, |acc, i| acc.saturating_mul((n - i) as u64) / (i as u64 + 1))
}

/// Localization q-values (FLR) of the target PSMs at `psm_q` (rank 1), for
/// each modification on its own: see the module documentation. A
/// modification with nothing to localize gets 0; one localized on a decoy
/// residue, and those of all other PSMs, 1. A PSM's q-value is the highest
/// of its modifications'.
pub fn assign_localization_q(features: &mut [Feature], psm_q: f32) {
    /// One localized modification of a PSM, as the FLR ranks it
    struct Ranked {
        feature: usize,
        group: usize,
        score: f32,
        delta: f32,
        decoy_site: bool,
    }
    let mut by_mod: std::collections::BTreeMap<ModId, (Vec<Ranked>, u64, u64)> = Default::default();
    for (i, f) in features.iter_mut().enumerate() {
        let passing = f.rank == 1 && f.label == 1 && f.spectrum_q <= psm_q;
        let Some(loc) = f.localization.as_mut() else { continue };
        for (g, group) in loc.groups.iter_mut().enumerate() {
            group.q = 1.0;
            if !passing {
                continue;
            }
            if group.arrangements <= 1 {
                group.q = 0.0;
                continue;
            }
            let (ranked, t_pos, d_pos) = by_mod.entry(group.modification).or_default();
            ranked.push(Ranked { feature: i, group: g, score: group.score, delta: group.delta, decoy_site: group.decoy_site });
            *t_pos += group.target_positions as u64;
            *d_pos += group.decoy_positions as u64;
        }
    }
    for (mut ranked, t_pos, d_pos) in by_mod.into_values() {
        if d_pos == 0 {
            // no decoy residue in any peptide with this modification: its
            // FLR cannot be estimated (q stays 1)
            continue;
        }
        ranked.sort_by(|a, b| b.score.total_cmp(&a.score).then(b.delta.total_cmp(&a.delta)));
        let correction = t_pos as f64 / d_pos as f64;
        let (mut targets, mut decoys) = (0.0f64, 0.0f64);
        let mut flr = Vec::with_capacity(ranked.len());
        for r in &ranked {
            match r.decoy_site {
                true => decoys += 1.0,
                false => targets += 1.0,
            }
            flr.push(if targets > 0.0 { (decoys * correction / targets).min(1.0) } else { 1.0 } as f32);
        }
        let mut minimum = 1.0f32;
        for (r, &f) in ranked.iter().zip(&flr).rev() {
            minimum = minimum.min(f);
            if let Some(loc) = features[r.feature].localization.as_mut() {
                loc.groups[r.group].q = if r.decoy_site { 1.0 } else { minimum };
            }
        }
    }
    for f in features.iter_mut() {
        if let Some(loc) = f.localization.as_mut() {
            loc.q = loc.groups.iter().map(|g| g.q).fold(0.0, f32::max);
        }
    }
}

/// ProForma 2.0 of `peptide` with its localization: modifications of a
/// localized group are written at the best target arrangement, with
/// `#g` localization scores where a site is not certain
pub fn localized_proforma(peptide: &Peptide, loc: &Localization, mods: &Modifications) -> String {
    use std::fmt::Write;
    let label = |id: ModId| mods.get(id).proforma_label();
    let n = peptide.sequence.len();
    // annotations per residue
    let mut at: Vec<Vec<String>> = vec![Vec::new(); n];
    let localized: Vec<ModId> = loc.groups.iter().map(|g| g.modification).collect();
    for (i, m) in peptide.modifications.iter().enumerate() {
        if let Some(m) = m.filter(|m| !localized.contains(m)) {
            at[i].push(label(m));
        }
    }
    for (g, group) in loc.groups.iter().enumerate() {
        let certain = group.sites.iter().filter(|s| group.best.contains(&s.0)).all(|s| s.1 >= 0.995);
        for &(i, p) in &group.sites {
            let i = i as usize;
            let best = group.best.contains(&(i as u16));
            match (certain, best) {
                (true, true) => at[i].push(label(group.modification)),
                (true, false) => {}
                (false, true) => at[i].push(format!("{}#g{}({:.2})", label(group.modification), g + 1, p)),
                (false, false) if p >= 0.005 => at[i].push(format!("#g{}({:.2})", g + 1, p)),
                (false, false) => {}
            }
        }
    }
    let mut s = String::new();
    if let Some(m) = peptide.nterm {
        let _ = write!(s, "[{}]-", label(m));
    }
    for (i, &r) in peptide.sequence.iter().enumerate() {
        s.push(r as char);
        for a in &at[i] {
            let _ = write!(s, "[{a}]");
        }
    }
    if let Some(m) = peptide.cterm {
        let _ = write!(s, "-[{}]", label(m));
    }
    s
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::Builder;

    fn localizer(
        variable: serde_json::Value,
        decoys: &str,
    ) -> (crate::database::Parameters, Localizer, crate::database::IndexedDatabase) {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 0, "min_len": 5},
            "static_mods": {"C": 57.021465},
            "variable_mods": variable,
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        let fasta = crate::fasta::Fasta::parse(">p\nMAGSPTSAEKLLLSAYR\n".into(), "rev_", false);
        let db = params.clone().build(fasta);
        let settings = LocalizationSettings { decoy_residues: decoys.into(), ..Default::default() };
        let localizer = Localizer::new(&params.modifications, &settings);
        (params, localizer, db)
    }

    /// A score that favors one site: the arrangements' probabilities, the
    /// best site, decoy residues, and the ProForma localization scores
    #[test]
    fn arrangements_and_probabilities() {
        let (params, localizer, db) =
            localizer(serde_json::json!({"S": [{"mass": 79.966331, "name": "UNIMOD:21"}], "T": [{"mass": 79.966331, "name": "UNIMOD:21"}]}), "A");
        let phospho = params.modifications.variable_mods[0].1;
        // MAGSPTSAEK with phospho on S3
        let peptide = db
            .peptides
            .iter()
            .find(|p| !p.decoy && p.sequence.as_ref() == b"MAGSPTSAEK" && p.modifications[3] == Some(phospho))
            .expect("phospho S3");
        // hyperscore: 30 at S3, 29 at T5, 20 elsewhere (A positions are decoys)
        let score = |p: &Peptide| match p.modifications.iter().position(|m| *m == Some(phospho)) {
            Some(3) => 30.0,
            Some(5) => 29.0,
            _ => 20.0,
        };
        let loc = localizer.localize_peptide(peptide, score).unwrap();
        assert_eq!(loc.groups.len(), 1);
        let g = &loc.groups[0];
        assert_eq!(g.best, vec![3]);
        assert_eq!(g.sites[0].0, 3);
        let e = std::f64::consts::E;
        // targets S3, T5, S6; decoys A1, A7: 5 arrangements
        assert_eq!((g.target_positions, g.decoy_positions, g.arrangements), (3, 2, 5));
        assert_eq!(loc.arrangements, 5);
        // reported over the targets (S3, T5, S6); the FLR score over all
        let targets = 1.0 + 1.0 / e + (-10.0f64).exp();
        assert!((loc.probability as f64 - 1.0 / targets).abs() < 1e-5, "{loc:?}");
        let all = targets + 2.0 * (-10.0f64).exp();
        assert!((g.score as f64 - 1.0 / all).abs() < 1e-5, "{loc:?}");
        assert!((g.delta - 1.0).abs() < 1e-6);
        assert!(!g.decoy_site);
        let proforma = localized_proforma(peptide, &loc, &db.modifications);
        assert!(proforma.starts_with("MAGS[UNIMOD:21#g1(0.73)]PT[#g1(0.27)]SAEK"), "{proforma}");
        // one real site (S3 of LLLSAYR; A4 a decoy): reported certain, but
        // still ranked for the FLR by its probability among all arrangements
        let only = db
            .peptides
            .iter()
            .find(|p| !p.decoy && p.sequence.as_ref() == b"LLLSAYR" && p.modifications[3] == Some(phospho))
            .expect("phospho S3 of LLLSAYR");
        let one = localizer.localize_peptide(only, |p| if p.modifications[3] == Some(phospho) { 30.0 } else { 25.0 }).unwrap();
        let g = &one.groups[0];
        assert_eq!((g.target_positions, g.decoy_positions, g.arrangements), (1, 1, 2));
        assert_eq!(one.probability, 1.0);
        assert!((g.score as f64 - 1.0 / (1.0 + (-5.0f64).exp())).abs() < 1e-5);
        assert_eq!(localized_proforma(only, &one, &db.modifications), "LLLS[UNIMOD:21]AYR");
        // a decoy residue wins
        let decoy = |p: &Peptide| match p.modifications.iter().position(|m| *m == Some(phospho)) {
            Some(1) => 40.0,
            _ => 20.0,
        };
        assert!(localizer.localize_peptide(peptide, decoy).unwrap().groups[0].decoy_site);
        // certain: plain ProForma
        let sure = |p: &Peptide| if p.modifications[3] == Some(phospho) { 60.0 } else { 20.0 };
        let loc = localizer.localize_peptide(peptide, sure).unwrap();
        assert_eq!(localized_proforma(peptide, &loc, &db.modifications), "MAGS[UNIMOD:21]PTSAEK");
    }

    /// One localized modification: 4 arrangements over 3 target and 1 decoy
    /// positions
    fn group(modification: ModId, score: f32, decoy_site: bool) -> SiteGroup {
        SiteGroup {
            modification,
            sites: Vec::new(),
            best: Vec::new(),
            score,
            delta: 0.0,
            arrangements: 4,
            decoy_site,
            target_positions: 3,
            decoy_positions: 1,
            q: 1.0,
        }
    }

    fn localized(spectrum_q: f32, groups: Vec<SiteGroup>) -> Feature {
        Feature {
            rank: 1,
            label: 1,
            spectrum_q,
            localization: Some(Localization { groups, probability: 1.0, arrangements: 4, q: 1.0 }),
            ..Default::default()
        }
    }

    /// The ids of phospho (S) and oxidation (M)
    fn two_mods() -> (ModId, ModId) {
        let (params, _, _) = localizer(
            serde_json::json!({"S": [{"mass": 79.966331, "name": "UNIMOD:21"}], "M": [{"mass": 15.994915, "name": "UNIMOD:35"}]}),
            "A",
        );
        let ids: Vec<ModId> = params.modifications.variable_mods.iter().map(|m| m.1).sorted().dedup().collect();
        assert_eq!(ids.len(), 2);
        (ids[0], ids[1])
    }

    #[test]
    fn flr_from_decoy_residues() {
        let (a, _) = two_mods();
        let feature = |q: f32, score: f32, decoy_site: bool| localized(q, vec![group(a, score, decoy_site)]);
        let mut features: Vec<Feature> = (0..30).map(|i| feature(0.001, 1.0 - i as f32 * 0.01, false)).collect();
        features.push(feature(0.001, 0.5, true));
        features.push(feature(0.001, 0.4, false));
        features.push(feature(0.5, 0.99, false));
        assign_localization_q(&mut features, 0.01);
        let q = |i: usize| features[i].localization.as_ref().unwrap().q;
        assert_eq!(q(0), 0.0);
        assert_eq!(q(30), 1.0, "decoy-site winner");
        // after the decoy: 1 * 3 / 31
        assert!((q(31) - 3.0 / 31.0).abs() < 1e-6, "{}", q(31));
        assert_eq!(q(32), 1.0, "PSM above 1% FDR");
    }

    /// Each modification has its own FLR: decoy-residue winners of one (say
    /// Met oxidation) do not raise the q-values of another
    #[test]
    fn flr_per_modification() {
        let (a, b) = two_mods();
        // a: 40 PSMs without a decoy-residue winner
        let mut features: Vec<Feature> =
            (0..40).map(|i| localized(0.001, vec![group(a, 0.99 - i as f32 * 0.001, false)])).collect();
        // b: 10 PSMs, the best 5 of them decoy-residue winners
        features.extend((0..10).map(|i| localized(0.001, vec![group(b, 1.0 - i as f32 * 0.001, i < 5)])));
        // both: a certain, b ranked after b's decoy-residue winners
        features.push(localized(0.001, vec![group(a, 0.995, false), group(b, 0.98, false)]));
        // nothing to localize
        let mut single = group(a, 1.0, false);
        single.arrangements = 1;
        features.push(localized(0.001, vec![single]));
        assign_localization_q(&mut features, 0.01);
        let loc = |i: usize| features[i].localization.as_ref().unwrap();
        assert!((0..40).all(|i| loc(i).q == 0.0), "b's decoy-residue winners count against b only");
        assert!((40..50).all(|i| loc(i).q == 1.0), "b: 5 decoys * 3 / at most 6 targets");
        assert_eq!((loc(50).groups[0].q, loc(50).groups[1].q, loc(50).q), (0.0, 1.0, 1.0), "the least certain group");
        assert_eq!(loc(51).q, 0.0);
    }
}
