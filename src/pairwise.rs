//! Pairwise reduction.
//!
//! For each column `i`, look for any other column `k` whose addition or
//! subtraction reduces `||v[i]||^2`. If found, apply the reduction and
//! restart. This is a simple post-processing pass that often shaves a few
//! more units off after L3FP.

use crate::basis::Basis;

/// Pairwise reduce a basis. Returns a new basis (the input is left
/// untouched).
pub fn pairwise_reduce<B: Basis>(basis: &B) -> B {
    let mut basis = basis.clone();
    let rank = basis.rank();
    if rank < 2 {
        return basis;
    }

    'restart: loop {
        for i in 0..rank {
            let v = basis.fdot(i, i);
            for k in 0..rank {
                if k == i {
                    continue;
                }
                let (mut s, a) = basis.fpair_size(i, k);
                let mut m: i64 = 1;
                if a < s {
                    m = -1;
                    s = a;
                }
                if s < v {
                    basis.column_reduce_i64(i, k, m);
                    continue 'restart;
                }
            }
        }
        break;
    }

    basis
}
