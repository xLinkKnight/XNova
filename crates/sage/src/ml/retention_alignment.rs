//! Global retention time alignment, based on Chen AT, et al. DART-ID
//! https://pubmed.ncbi.nlm.nih.gov/31260443/ (Fig 1A/B)
//!
//! 0) Transform all RTs into unit-less fractions of each run (0.0 - 1.0)
//! 1) Estimate each peptide's consensus RT across runs
//! 2) For each run, map its observed RTs onto the consensus
//!
//! Two methods (`retention_time_alignment`):
//! * `nonlinear` (default): the consensus is the median over the runs that
//!   share the peptide; a robust affine fit (least median of squares,
//!   RANSAC-style) rejects inconsistent landmarks, then a monotone
//!   piecewise-linear warp is fitted to binned medians of the rest (isotonic
//!   regression), following local gradient distortion. Runs with fewer than 16
//!   shared landmarks, or landmarks spanning less than a quarter of the run,
//!   keep the robust affine fit alone.
//! * `linear`: the mean over all runs and an ordinary least-squares line, as
//!   in upstream Sage and XNova 1.2 and earlier.
//!
//! The nonlinear method is adapted from Sage Plus
//! (`crates/sage/src/ml/retention_alignment.rs`, MIT,
//! <https://github.com/pgarrett-scripps/sage-plus>).
//!
//! If LFQ is enabled, MS1 traces are placed on the same aligned scale.

use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::sync::atomic::AtomicU32;

use super::matrix::Matrix;
use crate::database::PeptideIx;
use crate::scoring::Feature;
use dashmap::DashMap;
use fnv::FnvHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

type FnvDashMap<K, V> = DashMap<K, V, BuildHasherDefault<FnvHasher>>;

const MIN_NONLINEAR_LANDMARKS: usize = 16;
const MIN_NONLINEAR_SPAN: f64 = 0.25;
const MAX_RANSAC_POINTS: usize = 2_048;
const RANSAC_TRIALS: usize = 256;
const MIN_X_SEPARATION: f64 = 0.05;

/// How runs are aligned (`retention_time_alignment`)
#[derive(Copy, Clone, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum AlignmentMethod {
    /// Mean consensus and an ordinary least-squares line (upstream Sage)
    Linear,
    /// Median consensus, robust outlier rejection and a monotone piecewise-linear warp
    #[default]
    Nonlinear,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Alignment {
    pub file_id: usize,
    pub max_rt: f32,
    pub slope: f32,
    pub intercept: f32,
    /// Monotone `(observed, consensus)` knots in fractions of the run; empty
    /// for an affine alignment
    pub knots: Vec<(f32, f32)>,
}

impl Alignment {
    /// No alignment: the run's own (fractional) time scale
    pub fn identity(file_id: usize, max_rt: f32) -> Self {
        Self { file_id, max_rt, slope: 1.0, intercept: 0.0, knots: Vec::new() }
    }

    /// A retention time of this run (its own units) on the aligned scale
    pub fn transform(&self, rt: f32) -> f32 {
        let x = rt / self.max_rt;
        if self.knots.len() < 2 {
            return x * self.slope + self.intercept;
        }
        let upper = self.knots.partition_point(|&(knot_x, _)| knot_x <= x);
        let (left, right) = match upper {
            0 => (self.knots[0], self.knots[1]),
            n if n == self.knots.len() => (self.knots[n - 2], self.knots[n - 1]),
            n => (self.knots[n - 1], self.knots[n]),
        };
        let width = right.0 - left.0;
        if width <= f32::EPSILON {
            left.1
        } else {
            left.1 + (x - left.0) * (right.1 - left.1) / width
        }
    }

    /// An aligned time back in this run's own units: the inverse of
    /// [`Alignment::transform`] (flat warp segments map to their left end)
    pub fn inverse(&self, aligned: f32) -> f32 {
        let x = if self.knots.len() < 2 {
            if self.slope.abs() > f32::EPSILON {
                (aligned - self.intercept) / self.slope
            } else {
                aligned
            }
        } else {
            let upper = self.knots.partition_point(|&(_, knot_y)| knot_y <= aligned);
            let (left, right) = match upper {
                0 => (self.knots[0], self.knots[1]),
                n if n == self.knots.len() => (self.knots[n - 2], self.knots[n - 1]),
                n => (self.knots[n - 1], self.knots[n]),
            };
            let height = right.1 - left.1;
            if height <= f32::EPSILON {
                left.0
            } else {
                left.0 + (aligned - left.1) * (right.0 - left.0) / height
            }
        };
        x * self.max_rt
    }
}

fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mid = values.len() / 2;
    values.select_nth_unstable_by(mid, f64::total_cmp);
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        let lower = values[..mid].iter().copied().max_by(f64::total_cmp).unwrap_or(values[mid]);
        (lower + values[mid]) / 2.0
    }
}

/// The consensus of one landmark's fractional times in several runs
pub fn consensus(values: &mut [f64], method: AlignmentMethod) -> f64 {
    match method {
        AlignmentMethod::Linear => values.iter().sum::<f64>() / values.len() as f64,
        AlignmentMethod::Nonlinear => median(values),
    }
}

fn max_rt_by_file(features: &[Feature], n_files: usize) -> Vec<f64> {
    let max_rt = (0..n_files).map(|_| AtomicU32::new(0)).collect::<Vec<_>>();
    features.par_iter().for_each(|feat| {
        max_rt[feat.file_id].fetch_max(feat.rt.ceil() as u32, std::sync::atomic::Ordering::SeqCst);
    });
    max_rt
        .into_iter()
        .map(|v| v.load(std::sync::atomic::Ordering::Acquire) as f64)
        .collect()
}

/// Each confidently identified peptide's earliest RT in each file
fn rt_by_file(features: &[Feature]) -> FnvDashMap<PeptideIx, HashMap<usize, f64>> {
    let rts: FnvDashMap<PeptideIx, HashMap<usize, f64>> = DashMap::default();
    features
        .par_iter()
        .filter(|feat| feat.label == 1 && feat.spectrum_q <= 0.01)
        .for_each(|feat| {
            rts.entry(feat.peptide_idx)
                .or_default()
                .entry(feat.file_id)
                .and_modify(|f| *f = f.min(feat.rt as f64))
                .or_insert(feat.rt as f64);
        });
    rts
}

/// Fractional RTs: one row per peptide, one column per file (NaN where
/// absent); the nonlinear method uses only peptides seen in two or more runs
fn rt_matrix(features: &[Feature], max_rt: &[f64], method: AlignmentMethod) -> Matrix {
    let rts = rt_by_file(features);
    let mat: Vec<Vec<f64>> = rts
        .par_iter()
        .filter(|entry| method == AlignmentMethod::Linear || entry.value().len() >= 2)
        .map(|entry| {
            let mut row = vec![f64::NAN; max_rt.len()];
            for (&file_id, &rt) in entry.value() {
                row[file_id] = rt / max_rt[file_id];
            }
            row
        })
        .filter(|row| row.iter().any(|v| v.is_normal()))
        .collect();
    let n = mat.len();
    let mat: Vec<f64> = mat.into_par_iter().flatten().collect();
    Matrix::new(mat, n, max_rt.len())
}

fn ordinary_least_squares(points: &[(f64, f64)]) -> (f64, f64) {
    if points.len() < 2 {
        return (1.0, 0.0);
    }
    let n = points.len() as f64;
    let (sum_x, sum_y) = points.iter().fold((0.0, 0.0), |(sx, sy), &(x, y)| (sx + x, sy + y));
    let (mean_x, mean_y) = (sum_x / n, sum_y / n);
    let (ss_xy, ss_xx) = points
        .iter()
        .fold((0.0, 1e-8), |(xy, xx), &(x, y)| (xy + (x - mean_x) * (y - mean_y), xx + (x - mean_x).powi(2)));
    let slope = ss_xy / ss_xx;
    let intercept = mean_y - slope * mean_x;
    if slope.is_finite() && slope > 0.0 && intercept.is_finite() {
        (slope, intercept)
    } else {
        (1.0, 0.0)
    }
}

/// An affine trend from a deterministic least-median-of-squares variant of
/// RANSAC, and the indices of the landmarks consistent with it
fn robust_affine_inliers(points: &[(f64, f64)]) -> ((f64, f64), Vec<usize>) {
    if points.len() < 3 {
        return (ordinary_least_squares(points), (0..points.len()).collect());
    }
    let sample: Vec<_> = if points.len() <= MAX_RANSAC_POINTS {
        points.to_vec()
    } else {
        (0..MAX_RANSAC_POINTS).map(|i| points[i * (points.len() - 1) / (MAX_RANSAC_POINTS - 1)]).collect()
    };
    let mut state = 0x9e37_79b9_7f4a_7c15_u64 ^ points.len() as u64;
    let mut next_index = || {
        state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (state as usize) % sample.len()
    };
    let residuals_of =
        |(slope, intercept): (f64, f64)| -> Vec<f64> { sample.iter().map(|&(x, y)| (y - (slope * x + intercept)).abs()).collect() };
    let mut best = ordinary_least_squares(&sample);
    let mut best_median = median(&mut residuals_of(best));
    for _ in 0..RANSAC_TRIALS {
        let a = sample[next_index()];
        let b = sample[next_index()];
        if (b.0 - a.0).abs() < MIN_X_SEPARATION {
            continue;
        }
        let slope = (b.1 - a.1) / (b.0 - a.0);
        if !slope.is_finite() || slope <= 0.0 {
            continue;
        }
        let candidate = (slope, a.1 - slope * a.0);
        let candidate_median = median(&mut residuals_of(candidate));
        if candidate_median < best_median {
            best = candidate;
            best_median = candidate_median;
        }
    }
    let residuals: Vec<f64> = points.iter().map(|&(x, y)| (y - (best.0 * x + best.1)).abs()).collect();
    let residual_median = median(&mut residuals.clone());
    let mad = median(&mut residuals.iter().map(|r| (r - residual_median).abs()).collect::<Vec<_>>());
    // a liberal adaptive cutoff keeps real curvature for the nonlinear stage
    let cutoff = 0.02_f64.max(residual_median * 4.0).max(residual_median + 4.0 * 1.4826 * mad);
    let inliers = residuals.iter().enumerate().filter_map(|(i, &r)| (r <= cutoff).then_some(i)).collect();
    (best, inliers)
}

/// Pool-adjacent-violators isotonic regression
fn isotonic(values: &[f64], weights: &[usize]) -> Vec<f64> {
    struct Block {
        start: usize,
        end: usize,
        weighted_sum: f64,
        weight: usize,
    }
    let mut blocks: Vec<Block> = Vec::with_capacity(values.len());
    for (i, (&value, &weight)) in values.iter().zip(weights).enumerate() {
        blocks.push(Block { start: i, end: i + 1, weighted_sum: value * weight as f64, weight });
        while blocks.len() >= 2 {
            let n = blocks.len();
            let left = blocks[n - 2].weighted_sum / blocks[n - 2].weight as f64;
            let right = blocks[n - 1].weighted_sum / blocks[n - 1].weight as f64;
            if left <= right {
                break;
            }
            let right = blocks.pop().unwrap();
            let left = blocks.last_mut().unwrap();
            left.end = right.end;
            left.weighted_sum += right.weighted_sum;
            left.weight += right.weight;
        }
    }
    let mut fitted = vec![0.0; values.len()];
    for block in blocks {
        fitted[block.start..block.end].fill(block.weighted_sum / block.weight as f64);
    }
    fitted
}

/// Monotone knots through the binned medians of the inlier landmarks; none
/// when there are too few, or they cover too little of the run
fn nonlinear_knots(points: &[(f64, f64)], inliers: &[usize]) -> Vec<(f32, f32)> {
    if inliers.len() < MIN_NONLINEAR_LANDMARKS {
        return Vec::new();
    }
    let mut clean: Vec<_> = inliers.iter().map(|&i| points[i]).collect();
    clean.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    if clean.last().unwrap().0 - clean.first().unwrap().0 < MIN_NONLINEAR_SPAN {
        return Vec::new();
    }
    let bins = (clean.len() / 8).clamp(4, 24);
    let (mut xs, mut ys, mut weights) = (Vec::with_capacity(bins), Vec::with_capacity(bins), Vec::with_capacity(bins));
    for bin in 0..bins {
        let (start, end) = (bin * clean.len() / bins, (bin + 1) * clean.len() / bins);
        xs.push(median(&mut clean[start..end].iter().map(|p| p.0).collect::<Vec<_>>()));
        ys.push(median(&mut clean[start..end].iter().map(|p| p.1).collect::<Vec<_>>()));
        weights.push(end - start);
    }
    let ys = isotonic(&ys, &weights);
    xs.into_iter()
        .zip(ys)
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .map(|(x, y)| (x as f32, y as f32))
        .collect()
}

/// Fit one run's alignment from `(observed, consensus)` landmarks in
/// fractions of the run. With fewer than `min_points` landmarks, or an affine
/// slope outside `slope_range` (a wild fit from a handful of points is worse
/// than none), the run is left unaligned.
pub fn fit_alignment(
    file_id: usize,
    max_rt: f32,
    mut points: Vec<(f64, f64)>,
    method: AlignmentMethod,
    min_points: usize,
    slope_range: Option<(f64, f64)>,
) -> Alignment {
    points.retain(|(x, y)| x.is_finite() && y.is_finite());
    if points.len() < min_points.max(2) {
        return Alignment::identity(file_id, max_rt);
    }
    let (affine, knots) = match method {
        AlignmentMethod::Linear => (ordinary_least_squares(&points), Vec::new()),
        AlignmentMethod::Nonlinear => {
            points.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
            let (robust, mut inliers) = robust_affine_inliers(&points);
            let mut knots = nonlinear_knots(&points, &inliers);
            if knots.len() >= 2 {
                // A strongly bent gradient leaves real landmarks far from any
                // straight line, and the affine stage rejects them. Take back
                // those that follow the warp (extended past its last knots),
                // and fit it again. (Not in Sage Plus.)
                let warp = Alignment { file_id, max_rt: 1.0, slope: 1.0, intercept: 0.0, knots: knots.clone() };
                let residuals: Vec<f64> = points.iter().map(|&(x, y)| (y - warp.transform(x as f32) as f64).abs()).collect();
                let residual_median = median(&mut residuals.clone());
                let mad = median(&mut residuals.iter().map(|r| (r - residual_median).abs()).collect::<Vec<_>>());
                let cutoff = 0.02_f64.max(residual_median * 4.0).max(residual_median + 4.0 * 1.4826 * mad);
                inliers = residuals.iter().enumerate().filter_map(|(i, &r)| (r <= cutoff).then_some(i)).collect();
                knots = nonlinear_knots(&points, &inliers);
            }
            let kept: Vec<_> = inliers.iter().map(|&i| points[i]).collect();
            let affine = if kept.len() >= 2 { ordinary_least_squares(&kept) } else { robust };
            (affine, knots)
        }
    };
    if let Some((lo, hi)) = slope_range {
        if !(lo..=hi).contains(&affine.0) {
            return Alignment::identity(file_id, max_rt);
        }
    }
    Alignment { file_id, max_rt, slope: affine.0 as f32, intercept: affine.1 as f32, knots }
}

/// Upstream Sage's least-squares line, kept bit for bit for the linear method
fn legacy_linear(file_id: usize, max_rt: f64, points: &[(f64, f64)]) -> Alignment {
    let len = points.len();
    let (dot, sum_x, sum_y) =
        points.iter().fold((0.0f64, 0.0f64, 0.0f64), |(dot, sx, sy), &(x, y)| (dot + x * y, sx + x, sy + y));
    let x_mean = sum_x / len as f64;
    let y_mean = sum_y / len as f64;
    let ssxy = dot - len as f64 * x_mean * y_mean;
    let sx2 = points.iter().fold(1E-8f64, |sum, &(x, _)| sum + (x - x_mean).powi(2));
    let mut slope = ssxy / sx2;
    let mut intercept = y_mean - slope * x_mean;
    if !slope.is_finite() {
        slope = 1.0;
    }
    if !intercept.is_finite() {
        intercept = 0.0;
    }
    Alignment { file_id, max_rt: max_rt as f32, slope: slope as f32, intercept: intercept as f32, knots: Vec::new() }
}

/// Align runs with the default (nonlinear) method
pub fn global_alignment(features: &mut [Feature], n_files: usize) -> Vec<Alignment> {
    global_alignment_with_method(features, n_files, AlignmentMethod::default())
}

pub fn global_alignment_with_method(features: &mut [Feature], n_files: usize, method: AlignmentMethod) -> Vec<Alignment> {
    let max_rt = max_rt_by_file(features, n_files);
    let rt = rt_matrix(features, &max_rt, method);

    let consensus_rts: Vec<f64> = (0..rt.rows)
        .into_par_iter()
        .map(|row| consensus(&mut rt.row(row).filter(|rt| rt.is_finite()).collect::<Vec<_>>(), method))
        .collect();

    let alignments = (0..n_files)
        .into_par_iter()
        .map(|file_id| {
            let points: Vec<(f64, f64)> = rt
                .col(file_id)
                .zip(consensus_rts.iter().copied())
                .filter(|(x, y)| x.is_finite() && y.is_finite())
                .collect();
            let alignment = match method {
                AlignmentMethod::Linear => legacy_linear(file_id, max_rt[file_id], &points),
                AlignmentMethod::Nonlinear => fit_alignment(file_id, max_rt[file_id] as f32, points.clone(), method, 2, None),
            };
            // A run with few identifications (e.g. a sparse fraction) has few
            // landmarks, and a larger share of them are false matches
            let consistent = points
                .iter()
                .filter(|&&(x, y)| ((alignment.transform((x * max_rt[file_id]) as f32) as f64) - y).abs() <= 0.02)
                .count();
            if n_files > 1 && (points.len() < MIN_NONLINEAR_LANDMARKS || 2 * consistent < points.len()) {
                log::warn!(
                    "aligning file #{file_id}: {consistent} of {} landmark peptides agree with its alignment;                      with this few, its aligned retention times (and MS1 quantification in it) may be unreliable",
                    points.len()
                );
            }
            match method {
                AlignmentMethod::Linear => log::info!(
                    "aligning file #{file}: y = {m:.4}x + {b:.4}",
                    file = file_id,
                    m = alignment.slope,
                    b = alignment.intercept
                ),
                AlignmentMethod::Nonlinear => log::info!(
                    "aligning file #{file}: robust affine y = {m:.4}x + {b:.4}, nonlinear knots = {knots}",
                    file = file_id,
                    m = alignment.slope,
                    b = alignment.intercept,
                    knots = alignment.knots.len()
                ),
            }
            alignment
        })
        .collect::<Vec<Alignment>>();

    log::info!("aligned retention times across {} files", n_files);

    features.par_iter_mut().for_each(|feature| {
        feature.aligned_rt = alignments[feature.file_id].transform(feature.rt);
    });

    alignments
}

#[cfg(test)]
mod test {
    use super::*;

    /// A run whose gradient bends (consensus = observed + 0.08 sin(pi observed)),
    /// with a few wild landmarks
    fn bent(n: usize) -> Vec<(f64, f64)> {
        let mut points: Vec<(f64, f64)> = (1..=n)
            .map(|i| {
                let x = i as f64 / (n + 1) as f64;
                (x, x + 0.08 * (std::f64::consts::PI * x).sin())
            })
            .collect();
        for k in 0..5 {
            let x = 0.1 + 0.15 * k as f64;
            points.push((x, 1.0 - x));
        }
        points
    }

    #[test]
    fn nonlinear_follows_a_bent_gradient() {
        let points = bent(200);
        let linear = fit_alignment(0, 100.0, points.clone(), AlignmentMethod::Linear, 2, None);
        let warp = fit_alignment(0, 100.0, points.clone(), AlignmentMethod::Nonlinear, 2, None);
        assert!(warp.knots.len() >= 4, "{:?}", warp.knots);
        // error on the clean landmarks, in fractions of the run
        let error = |a: &Alignment| {
            points[..200].iter().map(|&(x, y)| (a.transform(x as f32 * 100.0) as f64 - y).abs()).fold(0.0, f64::max)
        };
        assert!(error(&warp) < 0.02, "warp {}", error(&warp));
        assert!(error(&linear) > 2.0 * error(&warp), "linear {} vs warp {}", error(&linear), error(&warp));
        // monotone, and inverse undoes transform
        let t: Vec<f32> = (0..=100).map(|i| warp.transform(i as f32)).collect();
        assert!(t.windows(2).all(|w| w[1] >= w[0]));
        for rt in [5.0f32, 37.5, 80.0] {
            assert!((warp.inverse(warp.transform(rt)) - rt).abs() < 0.05, "{rt}");
        }
    }

    #[test]
    fn few_landmarks_keep_an_affine_fit() {
        let points: Vec<(f64, f64)> = (0..10).map(|i| (i as f64 / 10.0, 0.02 + i as f64 / 10.0)).collect();
        let a = fit_alignment(3, 60.0, points, AlignmentMethod::Nonlinear, 3, Some((0.5, 2.0)));
        assert!(a.knots.is_empty());
        assert!((a.slope - 1.0).abs() < 1e-3 && (a.intercept - 0.02).abs() < 1e-3);
        // too few points, or a wild slope: unaligned
        assert_eq!(fit_alignment(1, 60.0, vec![(0.1, 0.1), (0.2, 0.2)], AlignmentMethod::Nonlinear, 3, None), Alignment::identity(1, 60.0));
        let wild: Vec<(f64, f64)> = (0..5).map(|i| (i as f64 / 10.0, i as f64 / 2.0)).collect();
        assert_eq!(fit_alignment(2, 60.0, wild, AlignmentMethod::Linear, 3, Some((0.5, 2.0))), Alignment::identity(2, 60.0));
    }

    #[test]
    fn methods_and_consensus() {
        assert_eq!(serde_json::to_value(AlignmentMethod::default()).unwrap(), serde_json::json!("nonlinear"));
        assert_eq!(serde_json::from_value::<AlignmentMethod>(serde_json::json!("linear")).unwrap(), AlignmentMethod::Linear);
        assert_eq!(consensus(&mut [0.1, 0.2, 0.9], AlignmentMethod::Nonlinear), 0.2);
        assert!((consensus(&mut [0.1, 0.2, 0.9], AlignmentMethod::Linear) - 0.4).abs() < 1e-12);
        let a = Alignment::identity(0, 50.0);
        assert_eq!((a.transform(25.0), a.inverse(0.5)), (0.5, 25.0));
    }
}
