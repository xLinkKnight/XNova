//! Linear Algebra, Machine Learning & FDR refinement

pub mod gauss;
pub mod kde;
pub mod linear_discriminant;
pub mod matrix;
pub mod mobility_model;
pub mod qvalue;
pub mod regression;
pub mod retention_alignment;
pub mod retention_model;

#[allow(dead_code)]
fn all_close(lhs: &[f64], rhs: &[f64], eps: f64) -> bool {
    lhs.iter()
        .zip(rhs.iter())
        .all(|(l, r)| (l - r).abs() <= eps)
}

/// Items per block of [`ordered_sum`] and [`ordered_fold`]
pub const ORDERED_BLOCK: usize = 1024;

/// The sum of `f` over `items`, computed in parallel but always in the same
/// order: fixed blocks of [`ORDERED_BLOCK`] items, each summed in order, then
/// the blocks' sums in order. A parallel `sum()` adds in an order that depends
/// on how rayon split the work, so its last bits, and everything fitted from
/// it, changed from run to run.
pub fn ordered_sum<T: Sync>(items: &[T], f: impl Fn(&T) -> f64 + Sync) -> f64 {
    use rayon::prelude::*;
    let blocks: Vec<f64> = items
        .par_chunks(ORDERED_BLOCK)
        .map(|block| block.iter().map(&f).sum::<f64>())
        .collect();
    blocks.iter().sum()
}

/// Fold `items` in parallel into accumulators, always in the same order (see
/// [`ordered_sum`]): each fixed block folds in order from `zero()`, then the
/// blocks' accumulators merge in order
pub fn ordered_fold<T: Sync, A: Send>(
    items: &[T],
    zero: impl Fn() -> A + Sync,
    add: impl Fn(A, &T) -> A + Sync,
    merge: impl Fn(A, A) -> A,
) -> A {
    use rayon::prelude::*;
    let blocks: Vec<A> = items
        .par_chunks(ORDERED_BLOCK)
        .map(|block| block.iter().fold(zero(), &add))
        .collect();
    blocks.into_iter().fold(zero(), merge)
}

pub fn norm(slice: &[f64]) -> f64 {
    slice.iter().fold(0.0, |acc, x| acc + x.powi(2)).sqrt()
}

pub fn mean(slice: &[f64]) -> f64 {
    slice.iter().sum::<f64>() / slice.len() as f64
}

pub fn std(slice: &[f64]) -> f64 {
    let mean = mean(slice);
    let x = slice.iter().fold(0.0, |acc, x| acc + (x - mean).powi(2));
    (x / slice.len() as f64).sqrt()
}

#[cfg(test)]
mod ordered_test {
    use super::*;

    /// The same bits on every call, whatever rayon does with the work; and
    /// the same sum as adding the blocks one by one
    #[test]
    fn ordered_sums_repeat_exactly() {
        let items: Vec<f64> = (0..100_000).map(|i| 1.0 / (1.0 + i as f64) * if i % 3 == 0 { -1.0 } else { 1.0 }).collect();
        let first = ordered_sum(&items, |x| x * 1.000001);
        for _ in 0..20 {
            assert_eq!(ordered_sum(&items, |x| x * 1.000001).to_bits(), first.to_bits());
        }
        let by_blocks: f64 = items
            .chunks(ORDERED_BLOCK)
            .map(|b| b.iter().map(|x| x * 1.000001).sum::<f64>())
            .sum();
        assert_eq!(first.to_bits(), by_blocks.to_bits());
        let folded = ordered_fold(&items, || (0.0f64, 0usize), |(s, n), x| (s + x * 1.000001, n + 1), |a, b| (a.0 + b.0, a.1 + b.1));
        assert_eq!(folded.1, items.len());
        assert_eq!(folded.0.to_bits(), first.to_bits());
    }
}
