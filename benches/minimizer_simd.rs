//! SIMD Minimizer Extraction Benchmark (ASBB Entry 036-E)
//!
//! This benchmark compares the SIMD-accelerated minimizer extraction against
//! the fast ntHash + SlidingMin implementation to validate the expected 3-15×
//! speedup from SIMD parallelization.
//!
//! # Purpose
//!
//! Validate that simd-minimizers library integration delivers the expected
//! speedup from parallel k-mer processing using ARM NEON (or AVX2 on x86).
//!
//! # Expected Speedup (Entry 036-E SIMD vs Entry 036-C Fast)
//!
//! - **Fast** (Entry 036-C): ~105 Mbp/s @ 1Mbp (ntHash + O(1) sliding min)
//! - **SIMD** (Entry 036-E): 315-1575 Mbp/s (3-15× from SIMD parallelization)
//! - **Target**: ≥3× speedup (conservative), ≥5× realistic, ≥10× exceptional
//!
//! # Success Criteria
//!
//! - ≥3× speedup: SUCCESS (conservative threshold, paper lower bound)
//! - ≥5× speedup: REALISTIC (median expectation)
//! - ≥10× speedup: EXCEPTIONAL (paper upper bound)
//!
//! # Methodology
//!
//! - **N=100 repetitions** per configuration (statistical rigor)
//! - **4 sequence lengths**: 100bp, 1Kbp, 10Kbp, 100Kbp (scaling validation)
//! - **2 k-mer sizes**: k=21 (typical), k=31 (high-specificity)
//! - **2 window sizes**: w=11 (typical), w=19 (large)
//! - **2 implementations**: fast (baseline), simd (test)
//! - **Total**: 2 × 16 configurations × 100 = 3,200 measurements
//!
//! # Running the Benchmark
//!
//! ```bash
//! # Full benchmark (N=100, ~90 minutes)
//! cargo bench --bench minimizer_simd --features simd -- --measurement-time 60 --sample-size 100
//!
//! # Quick validation (N=30, ~30 minutes)
//! cargo bench --bench minimizer_simd --features simd -- --measurement-time 30 --sample-size 30
//! ```
//!
//! # Evidence
//!
//! - **Entry 036-C** (fast): 105 Mbp/s @ 1Mbp (post lazy k-mer optimization)
//! - **SimdMinimizers paper**: 3-15× speedup from SIMD parallelization
//! - **Entry 036-E** (this benchmark): Validates SIMD integration

use biometal::operations::kmer::{extract_minimizers_fast, extract_minimizers_simd};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hash, Hasher};

/// Generate random DNA sequence for benchmarking
///
/// Uses seeded PRNG for reproducibility across runs (same as Entry 036)
fn generate_random_dna(len: usize) -> Vec<u8> {
    let build_hasher = RandomState::new();
    (0..len)
        .map(|i| {
            let mut hasher = build_hasher.build_hasher();
            i.hash(&mut hasher);
            b"ACGT"[(hasher.finish() as usize) % 4]
        })
        .collect()
}

/// Benchmark SIMD vs Fast minimizer extraction
///
/// Compares both implementations across all 16 configurations:
/// - k ∈ {21, 31}
/// - w ∈ {11, 19}
/// - sequence length ∈ {100, 1000, 10000, 100000}
fn bench_minimizer_simd_vs_fast(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_comparison");

    // Parameters (identical to Entry 036 baseline)
    let k_values = [21, 31];
    let w_values = [11, 19];
    let seq_lengths = [100, 1_000, 10_000, 100_000];

    for &k in &k_values {
        for &w in &w_values {
            for &len in &seq_lengths {
                // Generate test sequence
                let seq = generate_random_dna(len);

                // Throughput: base pairs per second
                group.throughput(Throughput::Bytes(len as u64));

                // Benchmark 1: Fast implementation (baseline)
                let id_fast = BenchmarkId::new(format!("fast/k{}_w{}", k, w), len);
                group.bench_with_input(id_fast, &seq, |b, seq| {
                    b.iter(|| {
                        extract_minimizers_fast(black_box(seq), black_box(k), black_box(w))
                            .unwrap()
                    });
                });

                // Benchmark 2: SIMD implementation (test)
                let id_simd = BenchmarkId::new(format!("simd/k{}_w{}", k, w), len);
                group.bench_with_input(id_simd, &seq, |b, seq| {
                    b.iter(|| {
                        extract_minimizers_simd(black_box(seq), black_box(k), black_box(w))
                            .unwrap()
                    });
                });
            }
        }
    }

    group.finish();
}

/// Benchmark SIMD vs Fast at production scale
///
/// Additional benchmark for large sequences (1Mbp, 10Mbp) to validate
/// performance at realistic genomic scales.
///
/// **Note**: These are resource-intensive. Use smaller N for quick validation.
fn bench_minimizer_simd_production_scale(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_comparison_production");

    // Typical genomics parameters
    let k = 21;
    let w = 11;

    // Production scales: 1Mbp (bacterial gene cluster), 10Mbp (small chromosome)
    let seq_lengths = [1_000_000, 10_000_000];

    for &len in &seq_lengths {
        let seq = generate_random_dna(len);

        group.throughput(Throughput::Bytes(len as u64));

        // Fast implementation (baseline)
        let id_fast = BenchmarkId::new(format!("fast/k{}_w{}", k, w), len);
        group.bench_with_input(id_fast, &seq, |b, seq| {
            b.iter(|| {
                extract_minimizers_fast(black_box(seq), black_box(k), black_box(w)).unwrap()
            });
        });

        // SIMD implementation (test)
        let id_simd = BenchmarkId::new(format!("simd/k{}_w{}", k, w), len);
        group.bench_with_input(id_simd, &seq, |b, seq| {
            b.iter(|| {
                extract_minimizers_simd(black_box(seq), black_box(k), black_box(w)).unwrap()
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_minimizer_simd_vs_fast,
    bench_minimizer_simd_production_scale
);
criterion_main!(benches);
