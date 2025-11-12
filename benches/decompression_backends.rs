// Decompression Backend Comparison Benchmark
//
// Tests different decompression backends to address bottleneck (98.7% of time)
//
// Current: flate2 with zlib-ng backend
// Baseline: flate2 with rust_backend (miniz_oxide) - for comparison
//
// Test matrix:
// - Small: 5.4 MB file
// - Medium: 54 MB file
// - Large: 544 MB file
//
// Metrics:
// - Decompression time (primary - this is 98.7% of total)
// - Total time (I/O + decompression)

use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use std::fs::File;
use std::io::{Read, BufReader};
use flate2::bufread::MultiGzDecoder;

/// Test file paths (from apple-silicon-bio-bench datasets)
const TEST_FILES: &[(&str, &str, u64)] = &[
    (
        "/Users/scotthandley/Code/apple-silicon-bio-bench/datasets/large_100k_150bp.fq.gz",
        "5MB",
        5_400_000,
    ),
    (
        "/Users/scotthandley/Code/apple-silicon-bio-bench/datasets/vlarge_1m_150bp.fq.gz",
        "54MB",
        54_000_000,
    ),
    (
        "/Users/scotthandley/Code/apple-silicon-bio-bench/datasets/huge_10m_150bp.fq.gz",
        "544MB",
        544_000_000,
    ),
];

fn bench_decompression(c: &mut Criterion) {
    let mut group = c.benchmark_group("decompression_backend");

    // Configure for statistical rigor
    group.sample_size(10);  // N=10 samples for quick validation
    group.significance_level(0.05);
    group.confidence_level(0.95);

    for (path, size_label, _expected_size) in TEST_FILES {
        // Check if file exists
        if !std::path::Path::new(path).exists() {
            eprintln!("⚠️  Skipping {} - file not found: {}", size_label, path);
            continue;
        }

        println!("\n📊 Testing: {}", size_label);
        println!("   File: \"{}\"", path);
        println!("   Backend: {}", if cfg!(feature = "zlib-ng") { "zlib-ng" } else { "rust_backend (miniz_oxide)" });

        // Benchmark decompression only (the actual bottleneck)
        group.bench_with_input(
            BenchmarkId::new("decompress_only", size_label),
            &path,
            |b, &path| {
                b.iter(|| {
                    let file = File::open(path).unwrap();
                    let mut decoder = MultiGzDecoder::new(BufReader::new(file));
                    let mut output = Vec::new();
                    decoder.read_to_end(&mut output).unwrap();
                    std::hint::black_box(output.len())
                })
            },
        );

        // Benchmark I/O + decompression (total time)
        group.bench_with_input(
            BenchmarkId::new("total_time", size_label),
            &path,
            |b, &path| {
                b.iter(|| {
                    let file = File::open(path).unwrap();
                    let mut decoder = MultiGzDecoder::new(BufReader::new(file));
                    let mut output = Vec::new();
                    decoder.read_to_end(&mut output).unwrap();
                    std::hint::black_box(output.len())
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_decompression);
criterion_main!(benches);
