//! Criterion benchmarks mirroring the Go `BenchmarkSmallBasis*` /
//! `BenchmarkLargeBasis*` suite.

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use lattice_reduction::{Int64Basis, l3fp, l3fp_deep};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

const SIZE: usize = 35;
// Match the Go benchmarks: random entries in [0, 0x7FFFFFFFFFFFFFF).
const RANGE: i64 = 0x7FFFFFFFFFFFFFF;

fn random_basis(seed: u64) -> Int64Basis {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut cols = Vec::with_capacity(SIZE);
    for _ in 0..SIZE {
        let mut col = Vec::with_capacity(SIZE);
        for _ in 0..SIZE {
            col.push(rng.gen_range(0..RANGE));
        }
        cols.push(col);
    }
    Int64Basis::from_columns(cols)
}

fn bench_l3fp_int64(c: &mut Criterion) {
    let mut counter = 0u64;
    c.bench_function("l3fp int64 35x35", |b| {
        b.iter_batched(
            || {
                counter = counter.wrapping_add(1);
                random_basis(counter)
            },
            |basis| {
                let _ = black_box(l3fp(&basis, 0.75));
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_l3fp_deep_int64(c: &mut Criterion) {
    let mut counter = 0u64;
    c.bench_function("l3fp_deep int64 35x35", |b| {
        b.iter_batched(
            || {
                counter = counter.wrapping_add(1);
                random_basis(counter)
            },
            |basis| {
                let _ = black_box(l3fp_deep(&basis, 0.75, 4));
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_l3fp_big(c: &mut Criterion) {
    let mut counter = 0u64;
    c.bench_function("l3fp big 35x35", |b| {
        b.iter_batched(
            || {
                counter = counter.wrapping_add(1);
                random_basis(counter).promote_to_big()
            },
            |basis| {
                let _ = black_box(l3fp(&basis, 0.75));
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

fn bench_l3fp_deep_big(c: &mut Criterion) {
    let mut counter = 0u64;
    c.bench_function("l3fp_deep big 35x35", |b| {
        b.iter_batched(
            || {
                counter = counter.wrapping_add(1);
                random_basis(counter).promote_to_big()
            },
            |basis| {
                let _ = black_box(l3fp_deep(&basis, 0.75, 4));
            },
            criterion::BatchSize::SmallInput,
        )
    });
}

criterion_group!(
    benches,
    bench_l3fp_int64,
    bench_l3fp_deep_int64,
    bench_l3fp_big,
    bench_l3fp_deep_big
);
criterion_main!(benches);
