//! Lattice basis reduction tools.
//!
//! Currently implements L^3-Reduction in floating point (L3FP) and
//! L3FP with deep insertion, plus pairwise reduction.
//!
//! See: "Lattice Basis Reduction: Improved Practical Algorithms and
//! Solving Subset Sum Problems" - C. P. Schnorr, M. Euchner.
//!
//! Lattice reduction is a technique for finding a small basis for a
//! lattice. Many problems in cryptography can be posed as problems on
//! a lattice with the solution being a small (smallest) vector in the
//! lattice (e.g. the knapsack problem).

pub mod basis;
pub mod l3fp;
pub mod pairwise;

pub use basis::{BigBasis, Basis, Int64Basis};
pub use l3fp::{l3fp, l3fp_deep};
pub use pairwise::pairwise_reduce;

/// `2^(-tau/2)` where `tau = 53` is the f64 mantissa precision in bits.
pub(crate) const TWO_POW_NEG_HALF_TOR: f64 = 1.110223e-16;
/// `2^(tau/2)` where `tau = 53` is the f64 mantissa precision in bits.
pub(crate) const TWO_POW_HALF_TOR: f64 = 2.2517998e+15;
