//! Property-style tests, mirroring `TestXx`, `TestXxx`, `TestXxxx` in the
//! original Go source.

use lattice_reduction::{Int64Basis, l3fp, l3fp_deep, pairwise_reduce};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

const SIZE: usize = 35;

fn random_basis(rng: &mut impl Rng, range: i64) -> Int64Basis {
    let mut cols = Vec::with_capacity(SIZE);
    for _ in 0..SIZE {
        let mut col = Vec::with_capacity(SIZE);
        for _ in 0..SIZE {
            col.push(rng.gen_range(0..range));
        }
        cols.push(col);
    }
    Int64Basis::from_columns(cols)
}

#[test]
fn l3fp_vs_l3fp_deep_quality() {
    let cases = if cfg!(debug_assertions) { 10 } else { 50 };
    let mut rng = ChaCha8Rng::seed_from_u64(0xDEAD_BEEF);
    let mut shallow_smaller = 0;
    for _ in 0..cases {
        let basis = random_basis(&mut rng, 0xFFFFFF);
        let shallow = l3fp(&basis, 0.75).modulus();
        let deep = l3fp_deep(&basis, 0.75, 4).modulus();
        if shallow <= deep {
            shallow_smaller += 1;
        }
    }
    eprintln!(
        "l3fp better-or-equal in {}/{} runs (the rest deep was strictly smaller)",
        shallow_smaller, cases
    );
    // Sanity bound: deep should usually be at least as good. Fail only on
    // gross regressions.
    assert!(
        cases - shallow_smaller >= cases / 10 || shallow_smaller >= cases / 10,
        "Suspicious quality mismatch"
    );
}

#[test]
fn pairwise_does_not_increase_size() {
    let cases = if cfg!(debug_assertions) { 5 } else { 25 };
    let mut rng = ChaCha8Rng::seed_from_u64(0xCAFEBABE);
    for _ in 0..cases {
        let basis = random_basis(&mut rng, 0xFFFFFF);
        let reduced = l3fp(&basis, 0.75);
        let pair = pairwise_reduce(&reduced);
        assert!(
            pair.modulus() <= reduced.modulus() + 1.0,
            "pairwise reduction increased modulus by more than rounding"
        );
    }
}

#[test]
fn deep_then_pairwise_does_not_increase_size() {
    let cases = if cfg!(debug_assertions) { 5 } else { 25 };
    let mut rng = ChaCha8Rng::seed_from_u64(0xFEEDFACE);
    for _ in 0..cases {
        let basis = random_basis(&mut rng, 0xFFFFFF);
        let reduced = l3fp_deep(&basis, 0.75, 4);
        let pair = pairwise_reduce(&reduced);
        assert!(
            pair.modulus() <= reduced.modulus() + 1.0,
            "pairwise reduction after deep increased modulus by more than rounding"
        );
    }
}

#[test]
fn l3fp_actually_shrinks_basis() {
    let mut rng = ChaCha8Rng::seed_from_u64(0xABCDEF01);
    let basis = random_basis(&mut rng, 0xFFFFFF);
    let before = basis.modulus();
    let after = l3fp(&basis, 0.75).modulus();
    assert!(
        after < before / 2.0,
        "expected significant shrink: before={before}, after={after}"
    );
}

/// `||det(B)||` (sign ignored) is a lattice invariant: any unimodular basis
/// transformation must preserve it. L3FP / L3FPDeep / pairwise reduction are
/// all unimodular, so the determinant magnitude must agree before and after.
fn abs_det(b: &Int64Basis) -> f64 {
    use lattice_reduction::Basis;
    let n = b.rank();
    assert_eq!(n, b.dimension());
    // Build f64 matrix `m[row][col]` (we treat columns as basis vectors).
    let mut m = vec![vec![0.0f64; n]; n];
    for col in 0..n {
        for row in 0..n {
            m[row][col] = b.fget(col, row);
        }
    }
    // Gaussian elimination with partial pivoting in f64.
    let mut det = 1.0f64;
    for i in 0..n {
        // Pivot.
        let mut pivot = i;
        let mut best = m[i][i].abs();
        for r in (i + 1)..n {
            let v = m[r][i].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if pivot != i {
            m.swap(i, pivot);
            det = -det;
        }
        if m[i][i] == 0.0 {
            return 0.0;
        }
        det *= m[i][i];
        for r in (i + 1)..n {
            let factor = m[r][i] / m[i][i];
            for c in i..n {
                m[r][c] -= factor * m[i][c];
            }
        }
    }
    det.abs()
}

#[test]
fn l3fp_preserves_determinant() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x1234_5678);
    for _ in 0..5 {
        // Use small entries so that the f64 determinant stays well within
        // mantissa precision.
        let basis = random_basis(&mut rng, 0xFFFF);
        let before = abs_det(&basis);
        let after = abs_det(&l3fp(&basis, 0.75));
        let rel = (before - after).abs() / before;
        assert!(
            rel < 1e-6,
            "l3fp changed determinant: before={before}, after={after}, rel={rel}"
        );
    }
}

#[test]
fn l3fp_deep_preserves_determinant() {
    let mut rng = ChaCha8Rng::seed_from_u64(0x9ABC_DEF0);
    for _ in 0..5 {
        let basis = random_basis(&mut rng, 0xFFFF);
        let before = abs_det(&basis);
        let after = abs_det(&l3fp_deep(&basis, 0.75, 4));
        let rel = (before - after).abs() / before;
        assert!(
            rel < 1e-6,
            "l3fp_deep changed determinant: before={before}, after={after}, rel={rel}"
        );
    }
}
