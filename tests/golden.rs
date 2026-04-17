//! Golden tests mirroring the `Example*` tests in the original Go source.
//!
//! Reduces the same `[[1,1,1],[-1,0,2],[3,5,6]]` example and checks the
//! result matches the Go implementation's output:
//!
//!   [0 1 0]
//!   [1 0 1]
//!   [-1 0 2]

use lattice_reduction::{Int64Basis, l3fp, l3fp_deep, pairwise_reduce};

fn small_basis() -> Int64Basis {
    Int64Basis::from_columns(vec![
        vec![1, 1, 1],
        vec![-1, 0, 2],
        vec![3, 5, 6],
    ])
}

fn expected_reduced() -> Int64Basis {
    Int64Basis::from_columns(vec![
        vec![0, 1, 0],
        vec![1, 0, 1],
        vec![-1, 0, 2],
    ])
}

#[test]
fn small_basis_l3fp_matches_go_golden() {
    let reduced = l3fp(&small_basis(), 0.75);
    assert_eq!(reduced, expected_reduced(), "got:\n{}", reduced);
}

#[test]
fn small_basis_l3fp_deep_matches_go_golden() {
    let reduced = l3fp_deep(&small_basis(), 0.75, 4);
    assert_eq!(reduced, expected_reduced(), "got:\n{}", reduced);
}

#[test]
fn small_basis_l3fp_big_matches_go_golden() {
    let big = small_basis().promote_to_big();
    let reduced = l3fp(&big, 0.75);
    let demoted = reduced
        .try_demote_to_i64()
        .expect("entries should fit in i64");
    assert_eq!(demoted, expected_reduced(), "got:\n{}", reduced);
}

#[test]
fn small_basis_l3fp_deep_big_matches_go_golden() {
    let big = small_basis().promote_to_big();
    let reduced = l3fp_deep(&big, 0.75, 4);
    let demoted = reduced
        .try_demote_to_i64()
        .expect("entries should fit in i64");
    assert_eq!(demoted, expected_reduced(), "got:\n{}", reduced);
}

#[test]
fn pairwise_reduce_is_idempotent_on_reduced_basis() {
    let reduced = l3fp(&small_basis(), 0.75);
    let twice = pairwise_reduce(&reduced);
    assert_eq!(twice.modulus(), reduced.modulus());
}
