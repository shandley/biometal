//! Performance benchmarks for [Experiment Name]
//!
//! # Benchmarking Protocol
//!
//! - Sample size: N=30 (statistical validity)
//! - Warmup: 3 seconds
//! - Measurement: 5 seconds
//! - Statistics: Mean, 95% CI, Cohen's d
//!
//! # Success Criteria
//!
//! - Target: [X]× speedup (NEON vs scalar)
//! - Minimum: [Y]× speedup for GO decision
//!
//! # Run with:
//!
//! ```bash
//! cargo bench
//! ```

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use experiment_name::{operation, operation_scalar};

#[cfg(target_arch = "aarch64")]
use experiment_name::neon::operation_neon;

/// Generate test data for benchmarks
fn generate_test_data(size: usize) -> Vec<u8> {
    // TODO: Generate appropriate test data for your experiment
    vec![0u8; size]
}

/// Benchmark: Compare scalar vs NEON implementations
fn benchmark_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("operation_comparison");

    // Sample size: N=30 for statistical validity
    group.sample_size(30);

    // Test with different data sizes
    for size in [1_000, 10_000, 100_000, 1_000_000].iter() {
        let data = generate_test_data(*size);

        // Set throughput for MB/s calculation
        group.throughput(Throughput::Bytes(*size as u64));

        // Benchmark scalar implementation
        group.bench_with_input(BenchmarkId::new("scalar", size), size, |b, _| {
            b.iter(|| operation_scalar(black_box(&data)))
        });

        // Benchmark NEON implementation (ARM only)
        #[cfg(target_arch = "aarch64")]
        group.bench_with_input(BenchmarkId::new("neon", size), size, |b, _| {
            b.iter(|| unsafe { operation_neon(black_box(&data)) })
        });

        // Benchmark dispatched version (uses NEON on ARM, scalar elsewhere)
        group.bench_with_input(BenchmarkId::new("optimized", size), size, |b, _| {
            b.iter(|| operation(black_box(&data)))
        });
    }

    group.finish();
}

/// Benchmark: Scalar implementation across different input sizes
///
/// This establishes the baseline performance.
fn benchmark_scalar_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("scalar_scaling");

    group.sample_size(30);

    for size in [100, 1_000, 10_000, 100_000, 1_000_000].iter() {
        let data = generate_test_data(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| operation_scalar(black_box(&data)))
        });
    }

    group.finish();
}

/// Benchmark: NEON implementation across different input sizes
///
/// This measures if NEON speedup is consistent across sizes.
#[cfg(target_arch = "aarch64")]
fn benchmark_neon_scaling(c: &mut Criterion) {
    let mut group = c.benchmark_group("neon_scaling");

    group.sample_size(30);

    for size in [100, 1_000, 10_000, 100_000, 1_000_000].iter() {
        let data = generate_test_data(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| unsafe { operation_neon(black_box(&data)) })
        });
    }

    group.finish();
}

/// Benchmark: Edge cases and special inputs
fn benchmark_edge_cases(c: &mut Criterion) {
    let mut group = c.benchmark_group("edge_cases");

    group.sample_size(30);

    // Empty input
    let empty: Vec<u8> = vec![];
    group.bench_function("empty", |b| {
        b.iter(|| operation_scalar(black_box(&empty)))
    });

    // Small input (< SIMD width)
    let small = generate_test_data(8);
    group.bench_function("small_8", |b| {
        b.iter(|| operation_scalar(black_box(&small)))
    });

    // Exactly SIMD width
    let exact = generate_test_data(16);
    group.bench_function("exact_16", |b| {
        b.iter(|| operation_scalar(black_box(&exact)))
    });

    // Misaligned (17 bytes = 16 + 1)
    let misaligned = generate_test_data(17);
    group.bench_function("misaligned_17", |b| {
        b.iter(|| operation_scalar(black_box(&misaligned)))
    });

    group.finish();
}

// Register benchmarks
criterion_group!(
    benches,
    benchmark_comparison,
    benchmark_scalar_scaling,
    benchmark_edge_cases,
);

// On ARM, also run NEON scaling benchmarks
#[cfg(target_arch = "aarch64")]
criterion_group!(neon_benches, benchmark_neon_scaling,);

#[cfg(target_arch = "aarch64")]
criterion_main!(benches, neon_benches);

#[cfg(not(target_arch = "aarch64"))]
criterion_main!(benches);
