//! False discovery rate control using double-competition (picked-peptide &
//! picked-protein) approaches
//!
//! Lin et al., https://pubmed.ncbi.nlm.nih.gov/36166314/
//! Savitski et al., https://pubmed.ncbi.nlm.nih.gov/25987413/

use crate::database::{IndexedDatabase, PeptideIx};
use crate::lfq::PrecursorId;
use crate::modification::ModId;
use crate::scoring::Feature;
use fnv::FnvHashMap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::BuildHasher;
use std::sync::Arc;

#[derive(Copy, Clone, Debug)]
pub struct Competition<Ix> {
    pub forward: f32,
    pub foward_ix: Option<Ix>,
    pub reverse: f32,
    pub reverse_ix: Option<Ix>,
}

struct Row<Ix> {
    ix: Ix,
    decoy: bool,
    score: f32,
    q: f32,
}

impl<Ix: Default + Send> Default for Competition<Ix> {
    fn default() -> Self {
        Self {
            forward: f32::MIN,
            reverse: f32::MIN,
            foward_ix: None,
            reverse_ix: None,
        }
    }
}

impl<Ix: Default + Send> Competition<Ix> {
    fn score(&self) -> f32 {
        self.forward.max(self.reverse)
    }

    fn is_decoy(&self) -> bool {
        self.reverse >= self.forward
    }

    /// Picked target-decoy q-values (Savitski et al. 2015; Lin et al. 2022):
    /// each target competes with its own decoy, only the winner of each pair
    /// is ranked (a tie goes to the decoy), and q is `(decoys + 1) / targets`
    /// after each complete group of tied scores, made monotone. The loser of
    /// each pair gets q = 1, so a target beaten by its own decoy never passes.
    ///
    /// Upstream Sage ranked both members of every pair and counted decoys as a
    /// sum of KDE posterior error probabilities; Sage Plus found that this lets
    /// a target beaten by its own decoy pass, and replaced it with this
    /// estimator (adapted from Sage Plus, MIT).
    fn assign_q_value<K, B>(
        scores: HashMap<K, Self, B>,
        threshold: f32,
    ) -> (HashMap<Ix, f32, B>, usize)
    where
        K: Eq + std::hash::Hash + Send,
        Ix: Eq + std::hash::Hash,
        B: BuildHasher + Default + Send,
    {
        let mut winners = Vec::with_capacity(scores.len());
        let mut losers = Vec::new();
        for (_, comp) in scores {
            let decoy = comp.is_decoy();
            let score = comp.score();
            let (winner, loser) = match decoy {
                true => (comp.reverse_ix, comp.foward_ix),
                false => (comp.foward_ix, comp.reverse_ix),
            };
            if let Some(ix) = winner {
                winners.push(Row { ix, decoy, score, q: 1.0 });
            }
            losers.extend(loser);
        }
        let passing = assign_count_q_values(&mut winners, threshold);
        let mut q = winners.into_iter().map(|row| (row.ix, row.q)).collect::<HashMap<_, _, B>>();
        q.extend(losers.into_iter().map(|ix| (ix, 1.0)));
        (q, passing)
    }
}

/// Count-based q-values at complete score thresholds, ties kept together;
/// the number of targets at or below `threshold`
fn assign_count_q_values<Ix>(scores: &mut [Row<Ix>], threshold: f32) -> usize {
    scores.sort_by(|a, b| b.score.total_cmp(&a.score));
    let mut decoys = 1.0_f64;
    let mut targets = 0.0_f64;
    let mut start = 0;
    while start < scores.len() {
        let mut end = start + 1;
        while end < scores.len() && scores[end].score == scores[start].score {
            end += 1;
        }
        for row in &scores[start..end] {
            if row.decoy {
                decoys += 1.0;
            } else {
                targets += 1.0;
            }
        }
        let q = (decoys / targets).min(1.0) as f32;
        for row in &mut scores[start..end] {
            row.q = q;
        }
        start = end;
    }
    let mut minimum = 1.0_f32;
    let mut passing = 0;
    for row in scores.iter_mut().rev() {
        minimum = minimum.min(row.q);
        row.q = minimum;
        passing += usize::from(!row.decoy && row.q <= threshold);
    }
    passing
}

/// Sequence and modifications of a peptide, used to pair targets with decoys
type PeptideKey = (
    Arc<[u8]>,
    Box<[Option<ModId>]>,
    Option<ModId>,
    Option<ModId>,
);

pub fn picked_peptide(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    let mut map: FnvHashMap<PeptideKey, Competition<PeptideIx>> = FnvHashMap::default();
    for feat in features.iter() {
        let peptide = &db[feat.peptide_idx];
        // Only reverse the peptide sequence if we generated decoys ourselves
        let key = match db.generate_decoys && peptide.decoy {
            true => {
                let reversed = peptide.reverse();
                (
                    reversed.sequence,
                    reversed.modifications,
                    reversed.nterm,
                    reversed.cterm,
                )
            }
            false => (
                peptide.sequence.clone(),
                peptide.modifications.clone(),
                peptide.nterm,
                peptide.cterm,
            ),
        };

        let entry = map.entry(key).or_default();
        match peptide.decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = Some(feat.peptide_idx);
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = Some(feat.peptide_idx);
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    features.par_iter_mut().for_each(|feat| {
        feat.peptide_q = scores[&feat.peptide_idx];
    });

    passing
}

pub fn picked_protein(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    // Critical: All non-proteotypic, non-unique, or shared peptides are discarded
    // else the assumptions of picked protein FDR are invalid. Shared peptides are
    // still reported, albeit with protein FDR = 1.0
    let mut map: FnvHashMap<_, Competition<String>> = FnvHashMap::default();
    for feat in features
        .iter()
        .filter(|x| db[x.peptide_idx].proteins.len() == 1)
    {
        let decoy = db[feat.peptide_idx].decoy;
        let entry = map.entry(&db[feat.peptide_idx].proteins).or_default();
        let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
        match decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = Some(proteins);
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = Some(proteins);
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    features
        .par_iter_mut()
        .filter(|x| db[x.peptide_idx].proteins.len() == 1)
        .for_each(|feat| {
            let proteins = db[feat.peptide_idx].proteins(&db.decoy_tag, db.generate_decoys);
            feat.protein_q = scores[&proteins];
        });

    passing
}

pub fn picked_protein_group(db: &IndexedDatabase, features: &mut [Feature]) -> usize {
    // Critical: All non-proteotypic, non-unique, or shared peptides are discarded
    // else the assumptions of picked group FDR are invalid. Shared peptides are
    // still reported, albeit with protein group FDR = 1.0
    let mut map: FnvHashMap<_, Competition<String>> = FnvHashMap::default();
    for feat in features
        .iter()
        .filter(|x| x.num_protein_groups == 1 && x.protein_groups.is_some())
    {
        let decoy = db[feat.peptide_idx].decoy;
        let entry = map.entry(feat.protein_groups.clone()).or_default();
        match decoy {
            true => {
                entry.reverse = entry.reverse.max(feat.discriminant_score);
                entry.reverse_ix = feat.protein_groups.clone();
            }
            false => {
                entry.forward = entry.forward.max(feat.discriminant_score);
                entry.foward_ix = feat.protein_groups.clone();
            }
        }
    }

    let (scores, passing) = Competition::assign_q_value(map, 0.01);

    features
        .par_iter_mut()
        .filter(|x| x.num_protein_groups == 1 && x.protein_groups.is_some())
        .for_each(|feat| {
            let protein_groups = feat.protein_groups.as_deref().unwrap().to_string();
            feat.protein_group_q = scores[&protein_groups];
        });

    passing
}

pub fn picked_precursor(
    peaks: &mut FnvHashMap<(PrecursorId, bool), (crate::lfq::Peak, Vec<f64>)>,
) -> usize {
    // let mut map: FnvHashMap<PeptideIx, Competition<(PeptideIx, bool)>> = FnvHashMap::default();
    // for (key, (peak, _)) in peaks.iter() {
    //     let entry = map.entry(key.0).or_default();
    //     match key.1 {
    //         true => {
    //             entry.reverse = entry.reverse.max(peak.score as f32);
    //             entry.reverse_ix = Some(*key);
    //         }
    //         false => {
    //             entry.forward = entry.forward.max(peak.score as f32);
    //             entry.foward_ix = Some(*key);
    //         }
    //     }
    // }
    let mut scores = peaks
        .par_iter()
        .map(|(&(ix, decoy), (peak, _))| Row {
            ix,
            decoy,
            score: peak.score as f32,
            q: 1.0,
        })
        .collect::<Vec<_>>();

    scores.par_sort_by(|a, b| b.score.total_cmp(&a.score));

    let mut decoy = 1.0;
    let mut target = 0.0;
    for score in scores.iter_mut() {
        match score.decoy {
            true => decoy += 1.0,
            false => target += 1.0,
        };
        score.q = decoy / target;
    }
    // Q-value is the minimum q-value at any given score threshold
    // `q = q[::-1].cummin()[::-1] in python`
    let mut q_min = 1.0f32;
    let mut passing = 0;
    for score in scores.iter_mut().rev() {
        q_min = q_min.min(score.q);
        score.q = q_min;
        if q_min <= 0.05 && !score.decoy {
            passing += 1;
        }
    }

    let scores = scores
        .into_par_iter()
        .map(|score| ((score.ix, score.decoy), score.q))
        .collect::<FnvHashMap<_, _>>();

    peaks.par_iter_mut().for_each(|(ix, (peak, _))| {
        peak.q_value = scores[ix];
    });
    passing
}

/// Row-level extraction q-values for `lfq_long.tsv` (adapted from Sage
/// Plus, MIT): in each run, a precursor's extraction competes with its
/// paired decoy's (the same precursor 11.06 Th up and earlier in time, traced
/// and integrated the same way) on [`crate::lfq::PeakFile::score`]. Winners
/// only, `(decoys + 1) / targets` at complete score thresholds, made
/// monotone; the loser of each pair, and an extraction without signal, get
/// q = 1. A diagnostic of how clean each run's extraction is, not a
/// calibrated match-between-runs FDR.
///
/// Returns the number of target rows at 1%.
pub fn picked_precursor_rows(
    peaks: &mut FnvHashMap<(PrecursorId, bool), (crate::lfq::Peak, Vec<f64>)>,
) -> usize {
    let files = peaks.values().map(|(peak, _)| peak.files.len()).max().unwrap_or(0);
    let ids: std::collections::BTreeSet<PrecursorId> = peaks.keys().map(|&(id, _)| id).collect();
    let score = |id: PrecursorId, decoy: bool, file: usize| {
        peaks
            .get(&(id, decoy))
            .and_then(|(peak, _)| peak.files.get(file))
            .map(|f| f.score as f32)
            // an extraction without signal does not compete (q = 1)
            .filter(|&s| s > 0.0)
    };
    let mut winners = Vec::new();
    for id in ids {
        for file in 0..files {
            let winner = match (score(id, false, file), score(id, true, file)) {
                (Some(t), Some(d)) if t > d => Some((false, t)),
                (Some(_), Some(d)) => Some((true, d)),
                (Some(t), None) => Some((false, t)),
                (None, Some(d)) => Some((true, d)),
                (None, None) => None,
            };
            if let Some((decoy, score)) = winner {
                winners.push(Row { ix: (id, decoy, file), decoy, score, q: 1.0 });
            }
        }
    }
    let passing = assign_count_q_values(&mut winners, 0.01);
    // everyone starts as a loser (q = 1); winners get theirs
    for (peak, _) in peaks.values_mut() {
        peak.files.iter_mut().for_each(|f| f.q_value = 1.0);
    }
    for row in winners {
        let (id, decoy, file) = row.ix;
        if let Some((peak, _)) = peaks.get_mut(&(id, decoy)) {
            peak.files[file].q_value = row.q;
        }
    }
    passing
}

#[cfg(test)]
mod test {
    use super::*;

    fn pair(forward: Option<(usize, f32)>, reverse: Option<(usize, f32)>) -> Competition<usize> {
        Competition {
            forward: forward.map_or(f32::MIN, |f| f.1),
            foward_ix: forward.map(|f| f.0),
            reverse: reverse.map_or(f32::MIN, |r| r.1),
            reverse_ix: reverse.map(|r| r.0),
        }
    }

    /// Winners only, (decoys + 1) / targets; losers get q = 1
    #[test]
    fn picked_q_values() {
        let mut map: HashMap<usize, Competition<usize>, fnv::FnvBuildHasher> = HashMap::default();
        // ten targets that beat their decoys, scores 100 down to 91
        for k in 0..10 {
            map.insert(k, pair(Some((k, 100.0 - k as f32)), Some((100 + k, 10.0))));
        }
        // a target beaten by its own decoy (decoy 95.5 wins)
        map.insert(10, pair(Some((10, 95.0)), Some((110, 95.5))));
        // a lone target and a lone decoy below everything
        map.insert(11, pair(Some((11, 50.0)), None));
        map.insert(12, pair(None, Some((112, 49.0))));
        let (q, passing) = Competition::assign_q_value(map, 0.2);
        // the target that lost to its decoy never passes, nor does any loser
        assert_eq!(q[&10], 1.0);
        assert!((100..110).all(|ix| q[&ix] == 1.0));
        // ranked winners: 5 targets, the decoy at 95.5, 5 targets, the target
        // at 50 (q = (1 + 1) / 11), the decoy at 49 ((2 + 1) / 11); made
        // monotone from the bottom, every target above 50 has 2/11 too
        for ix in [0, 4, 5, 9, 11] {
            assert!((q[&ix] - 2.0 / 11.0).abs() < 1e-6, "{ix}: {}", q[&ix]);
        }
        assert!((q[&112] - 3.0 / 11.0).abs() < 1e-6, "{}", q[&112]);
        assert_eq!(passing, 11);
    }

    #[test]
    fn ties_share_one_q_value() {
        let mut rows: Vec<Row<usize>> = vec![
            Row { ix: 0, decoy: false, score: 5.0, q: 1.0 },
            Row { ix: 1, decoy: true, score: 5.0, q: 1.0 },
            Row { ix: 2, decoy: false, score: 5.0, q: 1.0 },
            Row { ix: 3, decoy: false, score: 4.0, q: 1.0 },
        ];
        assign_count_q_values(&mut rows, 1.0);
        let q = |ix: usize| rows.iter().find(|r| r.ix == ix).unwrap().q;
        // the tie at 5.0 counts its decoy for everyone in it: (1 + 1) / 2 = 1,
        // then (1 + 1) / 3 below; monotone, so the tie gets 2/3 too, together
        for ix in 0..4 {
            assert!((q(ix) - 2.0 / 3.0).abs() < 1e-6, "{ix}: {}", q(ix));
        }
        // without monotonicity the tie would not be split: one value for all three
        assert_eq!((q(0), q(1)), (q(2), q(2)));
    }

    /// Each run's extraction competes with its own decoy's in that run
    #[test]
    fn picked_precursor_row_q_values() {
        use crate::lfq::{Peak, PeakFile};
        let file = |score: f64| PeakFile {
            apex_rt: 0.5,
            start_rt: 0.4,
            end_rt: 0.6,
            half_start_rt: 0.45,
            half_end_rt: 0.55,
            spectral_angle: 0.9,
            score,
            q_value: 0.0,
        };
        let peak = |scores: &[f64]| Peak { files: scores.iter().map(|&s| file(s)).collect(), ..Default::default() };
        let mut peaks: FnvHashMap<(PrecursorId, bool), (Peak, Vec<f64>)> = FnvHashMap::default();
        // 20 precursors, clean in run 0 (score 10 + k vs decoy 1); in run 1, precursor 0 loses to its decoy
        for k in 0..20u32 {
            let id = PrecursorId::Combined(PeptideIx(k));
            let run1 = if k == 0 { 2.0 } else { 5.0 };
            peaks.insert((id, false), (peak(&[10.0 + k as f64, run1]), vec![1.0, 1.0]));
            peaks.insert((id, true), (peak(&[1.0, 3.0]), vec![1.0, 1.0]));
        }
        // a run without signal for a precursor whose decoy has none either
        let empty = PrecursorId::Combined(PeptideIx(99));
        peaks.insert((empty, false), (peak(&[0.0, 0.0]), vec![0.0, 0.0]));
        let passing = picked_precursor_rows(&mut peaks);
        assert!(peaks[&(empty, false)].0.files.iter().all(|f| f.q_value == 1.0));
        let target = |k: u32| &peaks[&(PrecursorId::Combined(PeptideIx(k)), false)].0.files;
        let decoy = |k: u32| &peaks[&(PrecursorId::Combined(PeptideIx(k)), true)].0.files;
        // the target beaten by its decoy in run 1 never passes, nor does any losing decoy
        assert_eq!(target(0)[1].q_value, 1.0);
        assert!((0..20).all(|k| decoy(k)[0].q_value == 1.0));
        assert!((1..20).all(|k| decoy(k)[1].q_value == 1.0));
        // run-0 rows outrank everything: (0 + 1) / targets at their thresholds
        assert!((0..20).all(|k| target(k)[0].q_value <= 1.0 / 20.0 + 1e-6));
        // 39 winning targets, one winning decoy (run 1, precursor 0) at score 3, below all of them
        assert!((1..20).all(|k| (target(k)[1].q_value - 1.0 / 39.0).abs() < 1e-6));
        assert_eq!(passing, 0, "(0 + 1) / 39 > 1%");
    }
}
