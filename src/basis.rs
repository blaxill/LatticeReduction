//! Lattice basis representations.
//!
//! A `Basis` is a collection of column vectors. The trait abstracts over the
//! underlying integer type used to store the entries (currently `i64` and
//! arbitrary-precision `BigInt`).
//!
//! Note: As in the original Go implementation, the `column_reduce` operation
//! takes an `i64` multiplier. The driving algorithm guards against overflow
//! by tracking when `|mu|` exceeds `2^(tau/2)` and triggering a recovery step.

use num_bigint::BigInt;
use num_traits::ToPrimitive;
use std::fmt;

/// Trait shared by all lattice basis representations.
///
/// The "columns" terminology mirrors the original Go API: each `column` is a
/// basis vector, and `dimension` is the length of each vector.
pub trait Basis: Clone {
    /// Number of basis vectors.
    fn rank(&self) -> usize;

    /// Length of each basis vector.
    fn dimension(&self) -> usize;

    /// Dot product of two basis vectors as `f64` (may lose precision for `BigBasis`).
    fn fdot(&self, c1: usize, c2: usize) -> f64;

    /// Returns `(||v[c1] - v[c2]||^2, ||v[c1] + v[c2]||^2)` as `f64`.
    fn fpair_size(&self, c1: usize, c2: usize) -> (f64, f64);

    /// `column[c1] -= mu * column[c2]`.
    fn column_reduce_i64(&mut self, c1: usize, c2: usize, mu: i64);

    /// Swap two columns.
    fn column_swap(&mut self, c1: usize, c2: usize);

    /// Single entry as `f64`.
    fn fget(&self, column: usize, row: usize) -> f64;

    /// Copy the entire column into `dst` as `f64`. The default implementation
    /// loops over `fget`; `Int64Basis` overrides this with a vectorisable
    /// integer-to-float cast.
    #[inline]
    fn fcopy_column_into(&self, column: usize, dst: &mut [f64]) {
        debug_assert_eq!(dst.len(), self.dimension());
        for i in 0..self.dimension() {
            dst[i] = self.fget(column, i);
        }
    }

    /// Returns `true` if reading column entries via `fget` is cheap (e.g. an
    /// `i64`-to-`f64` cast). When `false`, callers may prefer to update an
    /// f64 mirror in-place rather than re-reading from the basis. The
    /// `BigBasis` impl returns `false` because reading involves a
    /// `BigInt::to_f64` conversion.
    #[inline]
    fn cheap_fget(&self) -> bool {
        true
    }
}

// =============================================================================
// Int64Basis
// =============================================================================

/// Lattice basis with `i64` entries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Int64Basis {
    /// Column-major storage. `data[col * dimension + row]`.
    data: Vec<i64>,
    rank: usize,
    dimension: usize,
}

impl Int64Basis {
    /// Construct from a vector of column vectors. All columns must have the
    /// same length.
    pub fn from_columns(cols: Vec<Vec<i64>>) -> Self {
        let rank = cols.len();
        let dimension = if rank == 0 { 0 } else { cols[0].len() };
        assert!(
            cols.iter().all(|c| c.len() == dimension),
            "Basis vectors have inconsistent dimension."
        );
        let mut data = Vec::with_capacity(rank * dimension);
        for c in cols {
            data.extend_from_slice(&c);
        }
        Self {
            data,
            rank,
            dimension,
        }
    }

    /// Construct an empty basis with the given shape; all entries are zero.
    pub fn zeros(rank: usize, dimension: usize) -> Self {
        Self {
            data: vec![0; rank * dimension],
            rank,
            dimension,
        }
    }

    #[inline]
    pub fn column(&self, c: usize) -> &[i64] {
        let start = c * self.dimension;
        &self.data[start..start + self.dimension]
    }

    #[inline]
    pub fn column_mut(&mut self, c: usize) -> &mut [i64] {
        let start = c * self.dimension;
        &mut self.data[start..start + self.dimension]
    }

    #[inline]
    pub fn get(&self, column: usize, row: usize) -> i64 {
        self.data[column * self.dimension + row]
    }

    #[inline]
    pub fn set(&mut self, column: usize, row: usize, value: i64) {
        self.data[column * self.dimension + row] = value;
    }

    /// Promote to a `BigBasis` (arbitrary-precision integer entries).
    pub fn promote_to_big(&self) -> BigBasis {
        let mut cols = Vec::with_capacity(self.rank);
        for c in 0..self.rank {
            let col = self
                .column(c)
                .iter()
                .map(|&v| BigInt::from(v))
                .collect::<Vec<_>>();
            cols.push(col);
        }
        BigBasis::from_columns(cols)
    }

    /// Mean squared entry value across the entire basis. Used by tests as a
    /// rough quality measure.
    pub fn modulus(&self) -> f64 {
        let mut acc = 0.0f64;
        for &v in &self.data {
            let f = v as f64;
            acc += f * f;
        }
        acc / (self.rank as f64 * self.dimension as f64)
    }
}

impl Basis for Int64Basis {
    #[inline]
    fn rank(&self) -> usize {
        self.rank
    }

    #[inline]
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn fdot(&self, c1: usize, c2: usize) -> f64 {
        // Match the Go `Int64Basis.FDot` semantics: convert to f64 first then
        // multiply, so that overflow in the integer multiply doesn't lose
        // sign information (the high-magnitude product is approximated in
        // f64).
        let a = self.column(c1);
        let b = self.column(c2);
        let mut r = 0.0f64;
        for i in 0..self.dimension {
            r += (a[i] as f64) * (b[i] as f64);
        }
        r
    }

    fn fpair_size(&self, c1: usize, c2: usize) -> (f64, f64) {
        let a = self.column(c1);
        let b = self.column(c2);
        let mut sub = 0.0f64;
        let mut add = 0.0f64;
        for i in 0..self.dimension {
            // Use wrapping arithmetic to match Go's signed-overflow behavior.
            let s = a[i].wrapping_sub(b[i]) as f64;
            let p = a[i].wrapping_add(b[i]) as f64;
            sub += s * s;
            add += p * p;
        }
        (sub, add)
    }

    fn column_reduce_i64(&mut self, c1: usize, c2: usize, mu: i64) {
        let dim = self.dimension;
        let (a_start, b_start) = (c1 * dim, c2 * dim);
        // column[c1] -= mu * column[c2]. Wrapping to match Go.
        if c1 == c2 {
            for i in 0..dim {
                let v = self.data[a_start + i];
                self.data[a_start + i] = v.wrapping_sub(mu.wrapping_mul(v));
            }
            return;
        }
        // Split into disjoint mutable slices.
        let (lower, upper) = if a_start < b_start {
            let (l, u) = self.data.split_at_mut(b_start);
            (&mut l[a_start..a_start + dim], &mut u[..dim])
        } else {
            let (l, u) = self.data.split_at_mut(a_start);
            (&mut u[..dim], &mut l[b_start..b_start + dim])
        };
        for i in 0..dim {
            lower[i] = lower[i].wrapping_sub(mu.wrapping_mul(upper[i]));
        }
    }

    fn column_swap(&mut self, c1: usize, c2: usize) {
        if c1 == c2 {
            return;
        }
        let dim = self.dimension;
        let (a_start, b_start) = (c1 * dim, c2 * dim);
        // Pick split point and use disjoint mutable slices to swap.
        let (lo, hi) = if a_start < b_start {
            let (l, u) = self.data.split_at_mut(b_start);
            (&mut l[a_start..a_start + dim], &mut u[..dim])
        } else {
            let (l, u) = self.data.split_at_mut(a_start);
            (&mut u[..dim], &mut l[b_start..b_start + dim])
        };
        lo.swap_with_slice(hi);
    }

    fn fget(&self, column: usize, row: usize) -> f64 {
        self.get(column, row) as f64
    }

    #[inline]
    fn fcopy_column_into(&self, column: usize, dst: &mut [f64]) {
        let src = self.column(column);
        debug_assert_eq!(src.len(), dst.len());
        for i in 0..src.len() {
            dst[i] = src[i] as f64;
        }
    }

    #[inline]
    fn cheap_fget(&self) -> bool {
        true
    }
}

impl fmt::Display for Int64Basis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for c in 0..self.rank {
            writeln!(f, "{:?}", self.column(c))?;
        }
        Ok(())
    }
}

// =============================================================================
// BigBasis
// =============================================================================

/// Lattice basis with arbitrary-precision integer entries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BigBasis {
    columns: Vec<Vec<BigInt>>,
    rank: usize,
    dimension: usize,
}

impl BigBasis {
    pub fn from_columns(cols: Vec<Vec<BigInt>>) -> Self {
        let rank = cols.len();
        let dimension = if rank == 0 { 0 } else { cols[0].len() };
        assert!(
            cols.iter().all(|c| c.len() == dimension),
            "Basis vectors have inconsistent dimension."
        );
        Self {
            columns: cols,
            rank,
            dimension,
        }
    }

    #[inline]
    pub fn column(&self, c: usize) -> &[BigInt] {
        &self.columns[c]
    }

    #[inline]
    pub fn get(&self, column: usize, row: usize) -> &BigInt {
        &self.columns[column][row]
    }
}

#[inline]
fn bigint_to_f64(v: &BigInt) -> f64 {
    // BigInt -> f64 conversion. `to_f64` returns `Some(infty)` for overflow
    // which matches Go's `big.Rat.Float64` behavior closely enough for the
    // floating-point bookkeeping in the L3FP algorithm.
    v.to_f64().unwrap_or(f64::INFINITY)
}

impl Basis for BigBasis {
    #[inline]
    fn rank(&self) -> usize {
        self.rank
    }

    #[inline]
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn fdot(&self, c1: usize, c2: usize) -> f64 {
        let a = &self.columns[c1];
        let b = &self.columns[c2];
        let mut r = 0.0f64;
        for i in 0..self.dimension {
            r += bigint_to_f64(&(&a[i] * &b[i]));
        }
        r
    }

    fn fpair_size(&self, c1: usize, c2: usize) -> (f64, f64) {
        let a = &self.columns[c1];
        let b = &self.columns[c2];
        let mut sub = 0.0f64;
        let mut add = 0.0f64;
        for i in 0..self.dimension {
            let s = &a[i] - &b[i];
            let p = &a[i] + &b[i];
            sub += bigint_to_f64(&(&s * &s));
            add += bigint_to_f64(&(&p * &p));
        }
        (sub, add)
    }

    fn column_reduce_i64(&mut self, c1: usize, c2: usize, mu: i64) {
        // Out-of-place to dodge borrow constraints; the reduction is
        // already O(dimension * BigInt-mul) so the extra alloc is minor.
        let mu_big = BigInt::from(mu);
        let dim = self.dimension;
        // SAFETY: c1 != c2 in practice (checked algorithmically), but we don't
        // assume it. Read column[c2] entries through a clone-free reference by
        // taking the columns out and re-inserting.
        if c1 == c2 {
            for i in 0..dim {
                let prod = &mu_big * &self.columns[c1][i];
                self.columns[c1][i] = &self.columns[c1][i] - &prod;
            }
            return;
        }
        let (a, b) = if c1 < c2 {
            let (lo, hi) = self.columns.split_at_mut(c2);
            (&mut lo[c1], &hi[0])
        } else {
            let (lo, hi) = self.columns.split_at_mut(c1);
            (&mut hi[0], &lo[c2])
        };
        for i in 0..dim {
            let prod = &mu_big * &b[i];
            a[i] -= prod;
        }
    }

    fn column_swap(&mut self, c1: usize, c2: usize) {
        self.columns.swap(c1, c2);
    }

    fn fget(&self, column: usize, row: usize) -> f64 {
        bigint_to_f64(&self.columns[column][row])
    }

    #[inline]
    fn cheap_fget(&self) -> bool {
        false
    }
}

impl BigBasis {
    /// Demote to `Int64Basis` if all entries fit in `i64`.
    pub fn try_demote_to_i64(&self) -> Option<Int64Basis> {
        let mut data = Vec::with_capacity(self.rank * self.dimension);
        for c in 0..self.rank {
            for v in &self.columns[c] {
                data.push(v.to_i64()?);
            }
        }
        Some(Int64Basis {
            data,
            rank: self.rank,
            dimension: self.dimension,
        })
    }
}

impl fmt::Display for BigBasis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for c in 0..self.rank {
            let row: Vec<String> = self.columns[c].iter().map(|v| v.to_string()).collect();
            writeln!(f, "[{}]", row.join(" "))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small() -> Int64Basis {
        Int64Basis::from_columns(vec![vec![1, 1, 1], vec![-1, 0, 2], vec![3, 5, 6]])
    }

    #[test]
    fn rank_and_dimension() {
        let b = small();
        assert_eq!(b.rank(), 3);
        assert_eq!(b.dimension(), 3);
    }

    #[test]
    fn fdot_matches_manual() {
        let b = small();
        // [1,1,1] . [-1,0,2] = -1 + 0 + 2 = 1
        assert_eq!(b.fdot(0, 1), 1.0);
        // [1,1,1] . [3,5,6] = 3 + 5 + 6 = 14
        assert_eq!(b.fdot(0, 2), 14.0);
        // [-1,0,2] . [3,5,6] = -3 + 0 + 12 = 9
        assert_eq!(b.fdot(1, 2), 9.0);
    }

    #[test]
    fn fpair_size_matches_manual() {
        let b = small();
        // col0=[1,1,1], col1=[-1,0,2]
        // sub = [2,1,-1] -> 4+1+1 = 6
        // add = [0,1,3]  -> 0+1+9 = 10
        assert_eq!(b.fpair_size(0, 1), (6.0, 10.0));
    }

    #[test]
    fn column_swap() {
        let mut b = small();
        b.column_swap(0, 2);
        assert_eq!(b.column(0), &[3, 5, 6]);
        assert_eq!(b.column(2), &[1, 1, 1]);
    }

    #[test]
    fn column_reduce_subtracts_multiple() {
        let mut b = small();
        // col0 -= 2 * col1  ->  [1,1,1] - 2*[-1,0,2] = [3,1,-3]
        b.column_reduce_i64(0, 1, 2);
        assert_eq!(b.column(0), &[3, 1, -3]);
    }

    #[test]
    fn promote_to_big_round_trip() {
        let b = small();
        let big = b.promote_to_big();
        let back = big.try_demote_to_i64().unwrap();
        assert_eq!(back, b);
    }

    #[test]
    fn big_basis_fdot_agrees_with_int64() {
        let b = small();
        let big = b.promote_to_big();
        for i in 0..b.rank() {
            for j in 0..b.rank() {
                assert_eq!(big.fdot(i, j), b.fdot(i, j), "fdot({i},{j}) mismatch");
            }
        }
    }

    #[test]
    fn big_basis_column_reduce_agrees_with_int64() {
        let mut a = small();
        let mut b = small().promote_to_big();
        a.column_reduce_i64(0, 1, 2);
        b.column_reduce_i64(0, 1, 2);
        for r in 0..3 {
            assert_eq!(a.get(0, r), b.get(0, r).try_into().unwrap());
        }
    }
}

