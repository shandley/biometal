//! Minimizer Optimization Micro-Benchmark
//!
//! Tests different optimization strategies to determine which provides
//! the best performance improvement over baseline.
//!
//! # Strategies Tested
//!
//! 1. **Current**: Full Minimizer with kmer Vec allocation
//! 2. **Minimal**: Position + hash only (no k-mer storage)
//! 3. **LazyKmer**: Store position/hash, extract k-mer on-demand via method
//!
//! # Purpose
//!
//! Experimentally determine which optimization approach provides the best
//! speedup for closing the gap with SimdMinimizers (820 Mbp/s).

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hash, Hasher};

// Import the fast implementation components
use biometal::operations::nthash::NtHashIterator;
use biometal::operations::sliding_min::SlidingMin;

/// Generate random DNA sequence for benchmarking
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

// ===== Strategy 1: Current (with k-mer Vec) =====

#[derive(Debug, Clone)]
struct MinimizerFull {
    position: usize,
    hash: u64,
    kmer: Vec<u8>,
}

fn extract_minimizers_full(sequence: &[u8], k: usize, w: usize) -> Vec<MinimizerFull> {
    if k == 0 || w == 0 || k > sequence.len() {
        return Vec::new();
    }

    let hash_iter = NtHashIterator::new(sequence, k).unwrap();
    let mut sliding_min = SlidingMin::new(w).unwrap();
    let mut minimizers = Vec::new();
    let mut last_pos = None;

    for (pos, hash) in hash_iter.enumerate() {
        if let Some(min_elem) = sliding_min.push(hash, pos) {
            if last_pos != Some(min_elem.pos) {
                let kmer = sequence[min_elem.pos..min_elem.pos + k].to_vec(); // ALLOCATION
                minimizers.push(MinimizerFull {
                    position: min_elem.pos,
                    hash: min_elem.val,
                    kmer,
                });
                last_pos = Some(min_elem.pos);
            }
        }
    }

    minimizers
}

// ===== Strategy 2: Minimal (position + hash only) =====

#[derive(Debug, Clone, Copy)]
struct MinimizerMinimal {
    position: usize,
    hash: u64,
}

fn extract_minimizers_minimal(sequence: &[u8], k: usize, w: usize) -> Vec<MinimizerMinimal> {
    if k == 0 || w == 0 || k > sequence.len() {
        return Vec::new();
    }

    let hash_iter = NtHashIterator::new(sequence, k).unwrap();
    let mut sliding_min = SlidingMin::new(w).unwrap();
    let mut minimizers = Vec::new();
    let mut last_pos = None;

    for (pos, hash) in hash_iter.enumerate() {
        if let Some(min_elem) = sliding_min.push(hash, pos) {
            if last_pos != Some(min_elem.pos) {
                minimizers.push(MinimizerMinimal {
                    position: min_elem.pos,
                    hash: min_elem.val,
                }); // NO ALLOCATION
                last_pos = Some(min_elem.pos);
            }
        }
    }

    minimizers
}

// ===== Strategy 3: Lazy K-mer (extract on-demand) =====

#[derive(Debug, Clone)]
struct MinimizerLazy<'a> {
    position: usize,
    hash: u64,
    sequence: &'a [u8], // Reference to original sequence
    k: usize,
}

impl<'a> MinimizerLazy<'a> {
    /// Extract k-mer on-demand (zero-copy)
    fn kmer(&self) -> &'a [u8] {
        &self.sequence[self.position..self.position + self.k]
    }

    /// Extract k-mer as owned Vec (when needed)
    fn kmer_owned(&self) -> Vec<u8> {
        self.kmer().to_vec()
    }
}

fn extract_minimizers_lazy(sequence: &[u8], k: usize, w: usize) -> Vec<MinimizerLazy> {
    if k == 0 || w == 0 || k > sequence.len() {
        return Vec::new();
    }

    let hash_iter = NtHashIterator::new(sequence, k).unwrap();
    let mut sliding_min = SlidingMin::new(w).unwrap();
    let mut minimizers = Vec::new();
    let mut last_pos = None;

    for (pos, hash) in hash_iter.enumerate() {
        if let Some(min_elem) = sliding_min.push(hash, pos) {
            if last_pos != Some(min_elem.pos) {
                minimizers.push(MinimizerLazy {
                    position: min_elem.pos,
                    hash: min_elem.val,
                    sequence, // Store reference
                    k,
                }); // NO ALLOCATION
                last_pos = Some(min_elem.pos);
            }
        }
    }

    minimizers
}

/// Benchmark all three strategies
fn bench_minimizer_strategies(c: &mut Criterion) {
    let mut group = c.benchmark_group("minimizer_strategies");

    // Test configurations (focus on 100Kbp where we saw 28-35× vs target 100×)
    let configs = [
        (21, 11, 100_000),  // Current: 28×, Target: 100×
        (21, 19, 100_000),  // Current: 33×, Target: 100×
        (31, 11, 100_000),  // Current: 32×, Target: 100×
    ];

    for &(k, w, len) in &configs {
        let seq = generate_random_dna(len);
        let id = format!("k{}_w{}", k, w);

        group.throughput(Throughput::Bytes(len as u64));

        // Strategy 1: Full (current implementation)
        group.bench_with_input(
            BenchmarkId::new("full", &id),
            &seq,
            |b, seq| {
                b.iter(|| extract_minimizers_full(black_box(seq), black_box(k), black_box(w)))
            },
        );

        // Strategy 2: Minimal (no k-mer storage)
        group.bench_with_input(
            BenchmarkId::new("minimal", &id),
            &seq,
            |b, seq| {
                b.iter(|| extract_minimizers_minimal(black_box(seq), black_box(k), black_box(w)))
            },
        );

        // Strategy 3: Lazy (zero-copy k-mer)
        group.bench_with_input(
            BenchmarkId::new("lazy", &id),
            &seq,
            |b, seq| {
                b.iter(|| extract_minimizers_lazy(black_box(seq), black_box(k), black_box(w)))
            },
        );
    }

    group.finish();
}

/// Benchmark k-mer extraction overhead specifically
fn bench_kmer_extraction_overhead(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_extraction");

    let seq = generate_random_dna(100_000);
    let k = 21;
    let w = 11;

    // Extract once to get minimizers
    let minimizers_lazy = extract_minimizers_lazy(&seq, k, w);

    group.throughput(Throughput::Elements(minimizers_lazy.len() as u64));

    // Benchmark lazy access (zero-copy)
    group.bench_function("lazy_access", |b| {
        b.iter(|| {
            for m in &minimizers_lazy {
                black_box(m.kmer()); // Zero-copy slice
            }
        })
    });

    // Benchmark owned extraction (when needed)
    group.bench_function("owned_extraction", |b| {
        b.iter(|| {
            for m in &minimizers_lazy {
                black_box(m.kmer_owned()); // Allocates Vec
            }
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_minimizer_strategies,
    bench_kmer_extraction_overhead
);
criterion_main!(benches);
