//! L^3-Reduction in floating point (L3FP) and the deep-insertion variant.
//!
//! Algorithm from:
//! "Lattice Basis Reduction: Improved Practical Algorithms and
//! Solving Subset Sum Problems" - C. P. Schnorr, M. Euchner.
//!
//! The basis itself is kept in exact representation (integer entries),
//! while the orthogonal Gram-Schmidt scratch (`bd`, `mu`, `c`) is held in
//! `f64`. When the floating-point view of a column drifts from the exact
//! value, we re-fetch from the basis (via `Basis::fdot` and `Basis::fget`).
//!
//! Layout note: `bd` and `mu` are stored as flat `Vec<f64>` rather than
//! `Vec<Vec<f64>>` to keep rows contiguous in memory (column-major is fine
//! since each "column" of the basis is a contiguous f64 row of length
//! `dimension`). This improves cache locality and lets `dot()` auto-vectorize.

use crate::basis::Basis;
use crate::{TWO_POW_HALF_TOR, TWO_POW_NEG_HALF_TOR};

/// Sum-of-products of two equal-length slices.
#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let n = a.len();
    let mut acc = 0.0f64;
    // Manual 4-way unroll. Helps the optimizer issue independent FMAs.
    let mut i = 0;
    let chunks = n / 4;
    let mut s0 = 0.0f64;
    let mut s1 = 0.0f64;
    let mut s2 = 0.0f64;
    let mut s3 = 0.0f64;
    for _ in 0..chunks {
        // Bounds-check elision via the loop bound.
        s0 += a[i] * b[i];
        s1 += a[i + 1] * b[i + 1];
        s2 += a[i + 2] * b[i + 2];
        s3 += a[i + 3] * b[i + 3];
        i += 4;
    }
    acc += s0 + s1 + s2 + s3;
    while i < n {
        acc += a[i] * b[i];
        i += 1;
    }
    acc
}

#[inline]
fn abs(v: f64) -> f64 {
    if v < 0.0 { -v } else { v }
}

/// Flat scratch space for L3FP. `bd` is a flattened `rank x dim` matrix and
/// `mu` is a flattened `rank x rank` matrix, both row-major (`row * stride`).
struct Scratch {
    bd: Vec<f64>,
    mu: Vec<f64>,
    c: Vec<f64>,
    rank: usize,
    dim: usize,
}

impl Scratch {
    fn new<B: Basis>(basis: &B) -> Self {
        let rank = basis.rank();
        let dim = basis.dimension();
        let mut bd = vec![0.0f64; rank * dim];
        let mu = vec![0.0f64; rank * rank];
        let c = vec![0.0f64; rank];
        for k in 0..rank {
            let row = &mut bd[k * dim..(k + 1) * dim];
            for i in 0..dim {
                row[i] = basis.fget(k, i);
            }
        }
        Self {
            bd,
            mu,
            c,
            rank,
            dim,
        }
    }

    #[inline]
    fn bd_row(&self, k: usize) -> &[f64] {
        &self.bd[k * self.dim..(k + 1) * self.dim]
    }

    #[inline]
    fn mu_row(&self, k: usize) -> &[f64] {
        &self.mu[k * self.rank..(k + 1) * self.rank]
    }

    #[inline]
    fn refresh<B: Basis>(&mut self, basis: &B, k: usize) {
        let row = &mut self.bd[k * self.dim..(k + 1) * self.dim];
        basis.fcopy_column_into(k, row);
    }

    #[inline]
    fn swap_bd(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        let dim = self.dim;
        let (lo_idx, hi_idx) = if a < b { (a, b) } else { (b, a) };
        let (lo_part, hi_part) = self.bd.split_at_mut(hi_idx * dim);
        let lo = &mut lo_part[lo_idx * dim..lo_idx * dim + dim];
        let hi = &mut hi_part[..dim];
        lo.swap_with_slice(hi);
    }
}

/// GSO update for column `k`. Mirrors step 2 of L3FP.
#[inline]
fn gso_update<B: Basis>(basis: &B, s: &mut Scratch, k: usize) {
    let dim = s.dim;
    let rank = s.rank;

    // Cache `dot(bd[k], bd[k])` once for the whole inner loop.
    let bd_k = &s.bd[k * dim..k * dim + dim];
    let dot_k_k = dot(bd_k, bd_k);
    s.c[k] = dot_k_k;
    if k == 1 {
        let bd_0 = &s.bd[0..dim];
        s.c[0] = dot(bd_0, bd_0);
    }

    for j in 0..k {
        let bd_k = &s.bd[k * dim..k * dim + dim];
        let bd_j = &s.bd[j * dim..j * dim + dim];

        let raw = dot(bd_k, bd_j);
        let dot_j_j = dot(bd_j, bd_j);
        let bound = TWO_POW_NEG_HALF_TOR * (dot_k_k * dot_j_j).sqrt();
        let mut v = if abs(raw) < bound {
            basis.fdot(k, j)
        } else {
            raw
        };

        // v -= sum_{i<j} mu[j][i] * mu[k][i] * c[i]
        let mu_j = &s.mu[j * rank..j * rank + j];
        let mu_k = &s.mu[k * rank..k * rank + j];
        let c = &s.c[..j];
        let mut sub = 0.0f64;
        for i in 0..j {
            sub += mu_j[i] * mu_k[i] * c[i];
        }
        v -= sub;

        v /= s.c[j];
        s.mu[k * rank + j] = v;
        s.c[k] -= v * v * s.c[j];
    }
}

/// Size-reduce column `k` against columns `0..k` (step 3). Returns whether
/// the driver should restart at `k - 1` (the original `Fc` flag).
///
/// We pick between two strategies for keeping `bd[k]` in sync with the
/// integer basis:
///
/// * `Basis::cheap_fget() == true` (e.g. `Int64Basis`): re-read `bd[k]`
///   from the basis after every integer reduction. The reread is just a
///   batch of `i64`-to-`f64` casts and stays exact.
///
/// * `Basis::cheap_fget() == false` (e.g. `BigBasis`): update `bd[k]
///   -= m * bd[j]` in-place. Avoids `O(dim)` `BigInt::to_f64` conversions
///   per step. `bd` may drift slightly from the exact basis but the algo
///   tolerates this through its `Fc` overflow check and `bound` check.
#[inline]
fn size_reduce<B: Basis>(basis: &mut B, s: &mut Scratch, k: usize) -> bool {
    let mut fc = false;
    let rank = s.rank;
    let dim = s.dim;
    let cheap = basis.cheap_fget();
    for j in (0..k).rev() {
        let mu_kj = s.mu[k * rank + j];
        if abs(mu_kj) > 0.5 {
            let m = (mu_kj + 0.5).floor();
            if abs(m) > TWO_POW_HALF_TOR {
                fc = true;
            }
            // mu[k][i] -= m * mu[j][i] for i in 0..j  (k != j by construction).
            {
                let (a, b) = if k < j {
                    let (lo, hi) = s.mu.split_at_mut(j * rank);
                    (&mut lo[k * rank..k * rank + j], &hi[..j])
                } else {
                    let (lo, hi) = s.mu.split_at_mut(k * rank);
                    (&mut hi[..j], &lo[j * rank..j * rank + j])
                };
                for i in 0..j {
                    a[i] -= m * b[i];
                }
            }
            s.mu[k * rank + j] -= m;

            basis.column_reduce_i64(k, j, m as i64);

            if cheap {
                s.refresh(basis, k);
            } else {
                // bd[k] -= m * bd[j] in-place.
                let (a, b) = if k < j {
                    let (lo, hi) = s.bd.split_at_mut(j * dim);
                    (&mut lo[k * dim..k * dim + dim], &hi[..dim])
                } else {
                    let (lo, hi) = s.bd.split_at_mut(k * dim);
                    (&mut hi[..dim], &lo[j * dim..j * dim + dim])
                };
                for i in 0..dim {
                    a[i] -= m * b[i];
                }
            }
        }
    }
    fc
}

/// L^3-Reduction (floating-point) of `basis` with parameter `delta`.
/// `delta` should satisfy `0.5 < delta < 1`.
pub fn l3fp<B: Basis>(basis: &B, delta: f64) -> B {
    let mut basis = basis.clone();
    let rank = basis.rank();
    if rank < 2 {
        return basis;
    }

    let mut s = Scratch::new(&basis);
    let mut k = 1usize;

    while k < rank {
        gso_update(&basis, &mut s, k);

        if size_reduce(&mut basis, &mut s, k) {
            k = if k > 1 { k - 1 } else { 1 };
            continue;
        }

        let mu_k_km1 = s.mu[k * rank + (k - 1)];
        if delta * s.c[k - 1] > s.c[k] + mu_k_km1 * mu_k_km1 * s.c[k - 1] {
            basis.column_swap(k, k - 1);
            s.swap_bd(k, k - 1);
            k = if k > 1 { k - 1 } else { 1 };
        } else {
            k += 1;
        }
    }

    basis
}

/// L^3-Reduction with deep insertion. `deep` controls how far back we will
/// reinsert a column when its projected length is shorter than an earlier
/// one.
pub fn l3fp_deep<B: Basis>(basis: &B, delta: f64, deep: usize) -> B {
    let mut basis = basis.clone();
    let rank = basis.rank();
    if rank < 2 {
        return basis;
    }

    let mut s = Scratch::new(&basis);
    let mut k = 1usize;

    'outer: while k < rank {
        gso_update(&basis, &mut s, k);

        if size_reduce(&mut basis, &mut s, k) {
            k = if k > 1 { k - 1 } else { 1 };
            continue;
        }

        if deep > 0 {
            let bd_k = s.bd_row(k);
            let mut cur = dot(bd_k, bd_k);
            let mut i = 0usize;
            let mu_k = s.mu_row(k);
            while i < k && delta * s.c[i] <= cur {
                cur -= mu_k[i] * mu_k[i] * s.c[i];
                i += 1;
            }

            // BUG (preserved from the original Go): the `k - i < deep` arm
            // of the condition causes infinite loops; mirror the Go and only
            // act on `i < deep`.
            if i < k && i < deep {
                // Cyclic right-shift of columns `(i..rank)`.
                for j in (i + 1..rank).rev() {
                    basis.column_swap(j, j - 1);
                    s.swap_bd(j, j - 1);
                }
                k = if i > 0 { i } else { 1 };
                continue 'outer;
            }
        }

        let mu_k_km1 = s.mu[k * rank + (k - 1)];
        if delta * s.c[k - 1] > s.c[k] + mu_k_km1 * mu_k_km1 * s.c[k - 1] {
            basis.column_swap(k, k - 1);
            s.swap_bd(k, k - 1);
            k = if k > 1 { k - 1 } else { 1 };
        } else {
            k += 1;
        }
    }

    basis
}
