//! Benchmarks for k-mer operations
//!
//! This benchmark suite validates Entry 034 findings from OPTIMIZATION_RULES.md:
//! - extract_kmers: Scalar-only (2.19-2.38× parallel for large datasets)
//! - extract_minimizers: Scalar-only (1.26× max parallel, below threshold)
//! - kmer_spectrum: Scalar-only (0.95-1.88× parallel, sometimes SLOWER)
//!
//! Key finding: K-mer operations are data-structure-bound (hash+HashMap), not
//! compute-bound. NEON provides no benefit.
//!
//! Run with: cargo bench --bench kmer_operations

use biometal::operations::kmer::{
    extract_kmers, extract_minimizers, kmer_spectrum, KmerExtractor,
};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

/// Generate random DNA sequence
fn generate_sequence(len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| [b'A', b'C', b'G', b'T'][i % 4])
        .collect()
}

/// Benchmark k-mer extraction across different k values
fn bench_extract_kmers(c: &mut Criterion) {
    let mut group = c.benchmark_group("extract_kmers");

    // Realistic Illumina read length: 150 bp
    let seq = generate_sequence(150);

    for k in [3, 6, 9, 12].iter() {
        group.throughput(Throughput::Bytes(150));
        group.bench_with_input(BenchmarkId::from_parameter(format!("k={}", k)), k, |b, &k| {
            b.iter(|| extract_kmers(black_box(&seq), k))
        });
    }

    group.finish();
}

/// Benchmark k-mer extraction across different sequence sizes
fn bench_extract_kmers_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("extract_kmers_sizes");

    let k = 6; // Common k-value for BERT models

    for size in [100, 1_000, 10_000, 100_000].iter() {
        let seq = generate_sequence(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| extract_kmers(black_box(&seq), k))
        });
    }

    group.finish();
}

/// Benchmark minimizer extraction (validates Entry 034: scalar optimal)
fn bench_extract_minimizers(c: &mut Criterion) {
    let mut group = c.benchmark_group("extract_minimizers");

    // Longer sequence for meaningful minimizer density
    let seq = generate_sequence(1_000);

    for (k, w) in [(3, 5), (6, 10), (9, 15)].iter() {
        group.throughput(Throughput::Bytes(1_000));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("k={}_w={}", k, w)),
            &(*k, *w),
            |b, &(k, w)| b.iter(|| extract_minimizers(black_box(&seq), k, w)),
        );
    }

    group.finish();
}

/// Benchmark k-mer spectrum (validates Entry 034: parallel makes it SLOWER)
fn bench_kmer_spectrum(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_spectrum");

    let k = 6;

    for num_seqs in [10, 100, 1_000].iter() {
        let sequences: Vec<Vec<u8>> = (0..*num_seqs).map(|_| generate_sequence(150)).collect();
        let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();

        let total_bytes = sequences.len() * 150;

        group.throughput(Throughput::Bytes(total_bytes as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{}_seqs", num_seqs)),
            &seq_refs,
            |b, seqs| b.iter(|| kmer_spectrum(black_box(seqs), k)),
        );
    }

    group.finish();
}

/// Benchmark parallel vs scalar k-mer extraction (validates Entry 034: 2.2× speedup)
fn bench_parallel_vs_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_vs_scalar");

    let k = 6;

    // Test with different dataset sizes to show threshold behavior
    for num_seqs in [100, 1_000, 10_000].iter() {
        let sequences: Vec<Vec<u8>> = (0..*num_seqs).map(|_| generate_sequence(150)).collect();
        let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();

        let total_bytes = sequences.len() * 150;

        // Scalar
        let scalar = KmerExtractor::new();
        group.throughput(Throughput::Bytes(total_bytes as u64));
        group.bench_with_input(
            BenchmarkId::new("scalar", format!("{}_seqs", num_seqs)),
            &seq_refs,
            |b, seqs| b.iter(|| scalar.extract(black_box(seqs), k)),
        );

        // Parallel (only test for ≥1000 sequences, where it actually activates)
        if *num_seqs >= 1000 {
            let parallel = KmerExtractor::with_parallel(4);
            group.throughput(Throughput::Bytes(total_bytes as u64));
            group.bench_with_input(
                BenchmarkId::new("parallel_4t", format!("{}_seqs", num_seqs)),
                &seq_refs,
                |b, seqs| b.iter(|| parallel.extract(black_box(seqs), k)),
            );
        }
    }

    group.finish();
}

/// Benchmark operations comparison at realistic read lengths
fn bench_kmer_operations_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_operations_comparison_150bp");

    // Realistic Illumina read length: 150 bp
    let seq = generate_sequence(150);
    let k = 6;
    let w = 10;

    group.bench_function("extract_kmers", |b| {
        b.iter(|| extract_kmers(black_box(&seq), k))
    });

    group.bench_function("extract_minimizers", |b| {
        b.iter(|| extract_minimizers(black_box(&seq), k, w))
    });

    // Spectrum with small dataset (1 sequence)
    let seq_refs = vec![seq.as_slice()];
    group.bench_function("kmer_spectrum_1seq", |b| {
        b.iter(|| kmer_spectrum(black_box(&seq_refs), k))
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_extract_kmers,
    bench_extract_kmers_sizes,
    bench_extract_minimizers,
    bench_kmer_spectrum,
    bench_parallel_vs_scalar,
    bench_kmer_operations_comparison,
);

criterion_main!(benches);
