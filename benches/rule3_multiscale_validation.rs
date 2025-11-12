//! Multi-scale validation of Rule 3 (Parallel BGZF)
//!
//! Tests hypothesis: Bounded streaming (8 blocks) shows benefit at larger scales
//!
//! Test scales:
//! - Small (5.4 MB, 100K seqs): Overhead may dominate
//! - Medium (54 MB, 1M seqs): Threshold boundary
//! - Large (544 MB, 10M seqs): Should show parallel benefit

use biometal::io::compression::{CompressedReader, DataSource};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

/// Get path to test datasets
fn get_dataset_path(name: &str) -> Option<PathBuf> {
    let home = std::env::var("HOME").ok()?;
    let path = PathBuf::from(&home)
        .join("Code/apple-silicon-bio-bench/datasets")
        .join(name);

    if path.exists() {
        Some(path)
    } else {
        None
    }
}

/// Multi-scale benchmark: Test parallel BGZF across file sizes
fn bench_multiscale_parallel_bgzf(c: &mut Criterion) {
    let mut group = c.benchmark_group("multiscale_bgzf");

    // Set longer warmup and measurement time for large files
    group.sample_size(10); // N=10 for large files (faster iteration)

    let test_files = vec![
        ("large_100k_150bp.fq.gz", "Small: 5.4 MB, 100K sequences", 5_400_000),
        ("vlarge_1m_150bp.fq.gz", "Medium: 54 MB, 1M sequences", 54_000_000),
        ("huge_10m_150bp.fq.gz", "Large: 544 MB, 10M sequences", 544_000_000),
    ];

    for (filename, description, _expected_size) in test_files {
        let path = match get_dataset_path(filename) {
            Some(p) => p,
            None => {
                eprintln!("⚠️  Skipping {}: File not found", filename);
                continue;
            }
        };

        let file_size = std::fs::metadata(&path).unwrap().len();
        let size_mb = file_size as f64 / 1_048_576.0;

        // Estimate block count (typical: ~12 KB per block compressed)
        let estimated_blocks = file_size / 12_000;

        println!("\n📊 Testing: {}", description);
        println!("   File: {:?}", path);
        println!("   Actual size: {:.1} MB ({} bytes)", size_mb, file_size);
        println!("   Estimated blocks: ~{}", estimated_blocks);
        println!("   Parallel chunks: {} blocks at a time (biometal bounded streaming)", 8);

        // Sequential decompression (baseline)
        group.bench_with_input(
            BenchmarkId::new("sequential", format!("{:.0}MB", size_mb)),
            &path,
            |b, path| {
                b.iter(|| {
                    let file = File::open(path).unwrap();
                    let mut reader = MultiGzDecoder::new(file);
                    let mut output = Vec::new();
                    reader.read_to_end(&mut output).unwrap();
                    black_box(output.len())
                });
            },
        );

        // Parallel decompression (biometal bounded streaming)
        group.bench_with_input(
            BenchmarkId::new("parallel_bounded", format!("{:.0}MB", size_mb)),
            &path,
            |b, path| {
                b.iter(|| {
                    let source = DataSource::from_path(path);
                    let mut reader = CompressedReader::new(source).unwrap();
                    let mut output = Vec::new();
                    reader.read_to_end(&mut output).unwrap();
                    black_box(output.len())
                });
            },
        );
    }

    group.finish();
}

/// Test block count scaling hypothesis
fn bench_block_count_effect(c: &mut Criterion) {
    let mut group = c.benchmark_group("block_count_scaling");
    group.sample_size(10);

    let test_files = vec![
        ("large_100k_150bp.fq.gz", "~450 blocks"),
        ("vlarge_1m_150bp.fq.gz", "~4,500 blocks"),
        ("huge_10m_150bp.fq.gz", "~45,000 blocks"),
    ];

    for (filename, block_description) in test_files {
        let path = match get_dataset_path(filename) {
            Some(p) => p,
            None => continue,
        };

        println!("\n📊 Block scaling test: {}", block_description);

        group.bench_with_input(
            BenchmarkId::new("parallel_bounded", block_description),
            &path,
            |b, path| {
                b.iter(|| {
                    let source = DataSource::from_path(path);
                    let mut reader = CompressedReader::new(source).unwrap();
                    let mut output = Vec::new();
                    reader.read_to_end(&mut output).unwrap();
                    black_box(output.len())
                });
            },
        );
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = bench_multiscale_parallel_bgzf, bench_block_count_effect
}

criterion_main!(benches);
