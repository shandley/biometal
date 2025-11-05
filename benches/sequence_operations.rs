//! Benchmarks for sequence manipulation operations (Phase 4)
//!
//! This benchmark suite evaluates performance of sequence manipulation primitives:
//! - Sequence transformations: reverse_complement, complement, reverse
//! - Record operations: extract_region, reverse_complement_record
//! - Trimming operations: quality-based trimming
//! - Masking operations: quality-based masking
//!
//! Purpose: Determine if NEON optimization is warranted (≥5× speedup threshold)
//!
//! Run with: cargo bench --bench sequence_operations
//! Run specific: cargo bench --bench sequence_operations -- reverse_complement

use biometal::operations::{
    complement, extract_region, mask_low_quality_copy, reverse, reverse_complement,
    reverse_complement_record, trim_quality_both, trim_quality_end,
};
use biometal::FastqRecord;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

/// Generate random DNA sequence
fn generate_sequence(len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| [b'A', b'C', b'G', b'T'][i % 4])
        .collect()
}

/// Generate random quality scores (Phred+33)
fn generate_quality(len: usize) -> Vec<u8> {
    // Mix of high and low quality
    (0..len)
        .map(|i| {
            if i % 5 == 0 {
                33 // Q0 (low quality, will be masked/trimmed)
            } else {
                73 // Q40 (high quality)
            }
        })
        .collect()
}

/// Generate FastqRecord for testing
fn generate_record(len: usize) -> FastqRecord {
    FastqRecord::new(
        format!("read_{}", len),
        generate_sequence(len),
        generate_quality(len),
    )
}

// ============================================================================
// Section 1: Sequence Transformations
// ============================================================================

/// Benchmark reverse_complement across different sequence sizes
fn bench_reverse_complement(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse_complement");

    // Realistic read lengths: 100bp (short), 150bp (Illumina), 300bp (long), 10K (PacBio)
    for size in [100, 150, 300, 1_000, 10_000, 100_000].iter() {
        let seq = generate_sequence(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| reverse_complement(black_box(&seq)))
        });
    }

    group.finish();
}

/// Benchmark complement across different sequence sizes
fn bench_complement(c: &mut Criterion) {
    let mut group = c.benchmark_group("complement");

    for size in [100, 150, 300, 1_000, 10_000].iter() {
        let seq = generate_sequence(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| complement(black_box(&seq)))
        });
    }

    group.finish();
}

/// Benchmark reverse across different sequence sizes
fn bench_reverse(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse");

    for size in [100, 150, 300, 1_000, 10_000].iter() {
        let seq = generate_sequence(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| reverse(black_box(&seq)))
        });
    }

    group.finish();
}

// ============================================================================
// Section 2: Record Operations
// ============================================================================

/// Benchmark reverse_complement_record
fn bench_reverse_complement_record(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse_complement_record");

    for size in [100, 150, 300, 1_000].iter() {
        let record = generate_record(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| reverse_complement_record(black_box(&record)))
        });
    }

    group.finish();
}

/// Benchmark extract_region
fn bench_extract_region(c: &mut Criterion) {
    let mut group = c.benchmark_group("extract_region");

    let record = generate_record(10_000);

    // Test different extraction sizes
    for extract_size in [10, 50, 100, 500, 1_000].iter() {
        group.throughput(Throughput::Bytes(*extract_size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(extract_size),
            extract_size,
            |b, &size| {
                let end = 100 + size; // Start at 100, extract 'size' bases
                b.iter(|| extract_region(black_box(&record), 100, end))
            },
        );
    }

    group.finish();
}

// ============================================================================
// Section 3: Quality-Based Trimming
// ============================================================================

/// Benchmark trim_quality_end
fn bench_trim_quality_end(c: &mut Criterion) {
    let mut group = c.benchmark_group("trim_quality_end");

    for size in [100, 150, 300, 1_000].iter() {
        let record = generate_record(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| trim_quality_end(black_box(&record), 20))
        });
    }

    group.finish();
}

/// Benchmark trim_quality_both
fn bench_trim_quality_both(c: &mut Criterion) {
    let mut group = c.benchmark_group("trim_quality_both");

    for size in [100, 150, 300, 1_000].iter() {
        let record = generate_record(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| trim_quality_both(black_box(&record), 20))
        });
    }

    group.finish();
}

// ============================================================================
// Section 4: Quality-Based Masking
// ============================================================================

/// Benchmark mask_low_quality
fn bench_mask_low_quality(c: &mut Criterion) {
    let mut group = c.benchmark_group("mask_low_quality");

    for size in [100, 150, 300, 1_000, 10_000].iter() {
        let record = generate_record(*size);

        group.throughput(Throughput::Bytes(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| mask_low_quality_copy(black_box(&record), 20))
        });
    }

    group.finish();
}

// ============================================================================
// Section 5: Realistic Workload Comparison
// ============================================================================

/// Compare all operations at realistic read length (150bp Illumina)
fn bench_operations_comparison_150bp(c: &mut Criterion) {
    let mut group = c.benchmark_group("operations_comparison_150bp");

    let seq = generate_sequence(150);
    let record = generate_record(150);

    group.bench_function("reverse_complement", |b| {
        b.iter(|| reverse_complement(black_box(&seq)))
    });

    group.bench_function("complement", |b| b.iter(|| complement(black_box(&seq))));

    group.bench_function("reverse", |b| b.iter(|| reverse(black_box(&seq))));

    group.bench_function("reverse_complement_record", |b| {
        b.iter(|| reverse_complement_record(black_box(&record)))
    });

    group.bench_function("trim_quality_end", |b| {
        b.iter(|| trim_quality_end(black_box(&record), 20))
    });

    group.bench_function("mask_low_quality", |b| {
        b.iter(|| mask_low_quality_copy(black_box(&record), 20))
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_reverse_complement,
    bench_complement,
    bench_reverse,
    bench_reverse_complement_record,
    bench_extract_region,
    bench_trim_quality_end,
    bench_trim_quality_both,
    bench_mask_low_quality,
    bench_operations_comparison_150bp,
);

criterion_main!(benches);
