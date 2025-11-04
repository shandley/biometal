//! Benchmarks for NEON-optimized operations
//!
//! This benchmark suite validates the claimed speedups from OPTIMIZATION_RULES.md:
//! - base_counting: 16.7× (Entry 020-025, Cohen's d = 4.82)
//! - gc_content: 20.3× (Entry 020-025, Cohen's d = 5.12)
//! - mean_quality: 25.1× (Entry 020-025, Cohen's d = 5.87)
//!
//! Run with: cargo bench --bench operations

use biometal::operations::{count_bases, gc_content, mean_quality};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

/// Generate random DNA sequence
fn generate_sequence(len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| [b'A', b'C', b'G', b'T'][i % 4])
        .collect()
}

/// Generate random quality scores (Phred+33)
fn generate_quality(len: usize) -> Vec<u8> {
    (0..len).map(|i| 33 + (i % 40) as u8).collect() // Q0-Q40
}

/// Benchmark base counting across different sequence sizes
fn bench_base_counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("base_counting");

    for size in [100, 1_000, 10_000, 100_000, 1_000_000].iter() {
        let seq = generate_sequence(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| count_bases(black_box(&seq)))
        });
    }

    group.finish();
}

/// Benchmark GC content calculation across different sequence sizes
fn bench_gc_content(c: &mut Criterion) {
    let mut group = c.benchmark_group("gc_content");

    for size in [100, 1_000, 10_000, 100_000, 1_000_000].iter() {
        let seq = generate_sequence(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| gc_content(black_box(&seq)))
        });
    }

    group.finish();
}

/// Benchmark mean quality calculation across different quality string sizes
fn bench_mean_quality(c: &mut Criterion) {
    let mut group = c.benchmark_group("mean_quality");

    for size in [100, 1_000, 10_000, 100_000, 1_000_000].iter() {
        let qual = generate_quality(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| mean_quality(black_box(&qual)))
        });
    }

    group.finish();
}

/// Benchmark scalar implementations (when available) for comparison
#[cfg(target_arch = "aarch64")]
mod scalar_comparison {
    use super::*;
    use biometal::operations::{gc_content_scalar, mean_quality_scalar};

    pub fn bench_gc_content_scalar(c: &mut Criterion) {
        let mut group = c.benchmark_group("gc_content_scalar");

        for size in [100, 1_000, 10_000, 100_000].iter() {
            let seq = generate_sequence(*size);

            group.throughput(Throughput::Bytes(*size as u64));
            group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
                b.iter(|| gc_content_scalar(black_box(&seq)))
            });
        }

        group.finish();
    }

    pub fn bench_mean_quality_scalar(c: &mut Criterion) {
        let mut group = c.benchmark_group("mean_quality_scalar");

        for size in [100, 1_000, 10_000, 100_000].iter() {
            let qual = generate_quality(*size);

            group.throughput(Throughput::Bytes(*size as u64));
            group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
                b.iter(|| mean_quality_scalar(black_box(&qual)))
            });
        }

        group.finish();
    }
}

/// Benchmark operations comparison at realistic read lengths
fn bench_operations_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("operations_comparison_150bp");

    // Realistic Illumina read length: 150 bp
    let seq = generate_sequence(150);
    let qual = generate_quality(150);

    group.bench_function("base_counting", |b| {
        b.iter(|| count_bases(black_box(&seq)))
    });

    group.bench_function("gc_content", |b| {
        b.iter(|| gc_content(black_box(&seq)))
    });

    group.bench_function("mean_quality", |b| {
        b.iter(|| mean_quality(black_box(&qual)))
    });

    group.finish();
}

/// Benchmark with sequences containing N bases (ambiguous)
fn bench_gc_content_with_n(c: &mut Criterion) {
    let mut group = c.benchmark_group("gc_content_with_n");

    for n_percent in [0, 10, 50].iter() {
        let size = 10_000;
        let mut seq = generate_sequence(size);

        // Replace n_percent of bases with N
        for (i, base) in seq.iter_mut().enumerate() {
            if (i * 100 / size) < *n_percent {
                *base = b'N';
            }
        }

        group.throughput(Throughput::Bytes(size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}%_N", n_percent)),
            &seq,
            |b, s| b.iter(|| gc_content(black_box(s))),
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_base_counting,
    bench_gc_content,
    bench_mean_quality,
    bench_operations_comparison,
    bench_gc_content_with_n,
);

#[cfg(target_arch = "aarch64")]
criterion_group!(
    scalar_benches,
    scalar_comparison::bench_gc_content_scalar,
    scalar_comparison::bench_mean_quality_scalar,
);

#[cfg(target_arch = "aarch64")]
criterion_main!(benches, scalar_benches);

#[cfg(not(target_arch = "aarch64"))]
criterion_main!(benches);
