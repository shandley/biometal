//! Benchmark to validate Rule 3 (Parallel BGZF) 6.5× speedup
//!
//! # Evidence Base
//!
//! **Entry 029** (apple-silicon-bio-bench):
//! - Small files (51 blocks): 5.48× speedup
//! - Large files (485 blocks): 6.50× speedup
//! - Method: Rayon-based parallel block decompression
//!
//! This benchmark replicates Entry 029 methodology with biometal's implementation.

use biometal::io::compression::{CompressedReader, DataSource};
use biometal::FastqStream;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

/// Get path to test datasets (apple-silicon-bio-bench datasets)
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

/// Benchmark parallel vs sequential BGZF decompression
///
/// Replicates Entry 029 methodology:
/// - Sequential: Standard gzip decompression (MultiGzDecoder)
/// - Parallel: Biometal's BoundedParallelBgzipReader (8 blocks in parallel)
fn bench_parallel_vs_sequential(c: &mut Criterion) {
    let mut group = c.benchmark_group("bgzf_decompression");

    // Test files matching Entry 029 methodology
    let test_files = vec![
        ("medium_10k_150bp.fq.gz", "10K sequences (Entry 029: ~51 blocks, 5.48× speedup)"),
        ("large_100k_150bp.fq.gz", "100K sequences (Entry 029: ~485 blocks, 6.50× speedup)"),
    ];

    for (filename, description) in test_files {
        let path = match get_dataset_path(filename) {
            Some(p) => p,
            None => {
                eprintln!(
                    "⚠️  Skipping {}: File not found. Run from apple-silicon-bio-bench repo.",
                    filename
                );
                continue;
            }
        };

        println!("\n📊 Benchmarking: {}", description);
        println!("   File: {:?}", path);

        let file_size = std::fs::metadata(&path).unwrap().len();
        println!("   Size: {} MB", file_size / 1_048_576);

        // Sequential decompression (baseline)
        group.bench_with_input(
            BenchmarkId::new("sequential", filename),
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

        // Parallel decompression (Rule 3)
        group.bench_with_input(
            BenchmarkId::new("parallel", filename),
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

/// Benchmark FASTQ parsing with parallel BGZF vs sequential
///
/// Real-world scenario: Parsing FASTQ records with decompression overhead
fn bench_fastq_parsing_with_decompression(c: &mut Criterion) {
    let mut group = c.benchmark_group("fastq_parsing_with_bgzf");

    let test_files = vec![
        ("medium_10k_150bp.fq.gz", "10K sequences"),
        ("large_100k_150bp.fq.gz", "100K sequences"),
    ];

    for (filename, description) in test_files {
        let path = match get_dataset_path(filename) {
            Some(p) => p,
            None => {
                eprintln!("⚠️  Skipping {}: File not found", filename);
                continue;
            }
        };

        println!("\n📊 FASTQ parsing: {}", description);

        // With parallel BGZF (biometal default)
        group.bench_with_input(
            BenchmarkId::new("parallel_bgzf", filename),
            &path,
            |b, path| {
                b.iter(|| {
                    let source = DataSource::from_path(path);
                    let mut stream = FastqStream::new(source).unwrap();
                    let mut count = 0;
                    for record in &mut stream {
                        let _record = record.unwrap();
                        count += 1;
                    }
                    black_box(count)
                });
            },
        );

        // With sequential decompression (for comparison)
        group.bench_with_input(
            BenchmarkId::new("sequential_gzip", filename),
            &path,
            |b, path| {
                b.iter(|| {
                    let file = File::open(path).unwrap();
                    let reader = BufReader::new(MultiGzDecoder::new(file));
                    let mut stream = FastqStream::from_reader(reader);
                    let mut count = 0;
                    for record in &mut stream {
                        let _record = record.unwrap();
                        count += 1;
                    }
                    black_box(count)
                });
            },
        );
    }

    group.finish();
}

/// Benchmark impact of block count on parallelism
///
/// Entry 029 showed speedup scales with number of blocks:
/// - 51 blocks → 5.48× speedup
/// - 485 blocks → 6.50× speedup
fn bench_scaling_with_block_count(c: &mut Criterion) {
    let mut group = c.benchmark_group("bgzf_scaling");

    let test_files = vec![
        ("tiny_100_150bp.fq.gz", "100 sequences (few blocks)"),
        ("small_1k_150bp.fq.gz", "1K sequences"),
        ("medium_10k_150bp.fq.gz", "10K sequences"),
        ("large_100k_150bp.fq.gz", "100K sequences (many blocks)"),
    ];

    for (filename, description) in test_files {
        let path = match get_dataset_path(filename) {
            Some(p) => p,
            None => {
                eprintln!("⚠️  Skipping {}: File not found", filename);
                continue;
            }
        };

        println!("\n📊 Scaling test: {}", description);

        group.bench_with_input(
            BenchmarkId::new("parallel", filename),
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
    config = Criterion::default().sample_size(30);  // N=30 for statistical rigor (Entry 029)
    targets = bench_parallel_vs_sequential, bench_fastq_parsing_with_decompression, bench_scaling_with_block_count
}

criterion_main!(benches);
