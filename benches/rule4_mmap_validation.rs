//! Rule 4 (Smart mmap) validation benchmark
//!
//! Validates Entry 032's 2.5× speedup claim for files ≥50 MB
//!
//! Expected results (from Entry 032):
//! - Small files (<50 MB): 0.99× (standard I/O faster, no mmap overhead)
//! - Large files (≥50 MB): 2.30-2.55× speedup (APFS prefetching + madvise)

use biometal::io::compression::{CompressedReader, DataSource};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

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

/// Benchmark Rule 4 (mmap) across file sizes
fn bench_rule4_mmap_threshold(c: &mut Criterion) {
    let mut group = c.benchmark_group("rule4_mmap");
    group.sample_size(10);

    let test_files = vec![
        ("large_100k_150bp.fq.gz", "Below threshold: 5.4 MB"),
        ("vlarge_1m_150bp.fq.gz", "At threshold: 54 MB"),
        ("huge_10m_150bp.fq.gz", "Above threshold: 544 MB"),
    ];

    for (filename, description) in test_files {
        let path = match get_dataset_path(filename) {
            Some(p) => p,
            None => {
                eprintln!("⚠️  Skipping {}: File not found", filename);
                continue;
            }
        };

        let file_size = std::fs::metadata(&path).unwrap().len();
        let size_mb = file_size as f64 / 1_048_576.0;

        println!("\n📊 Testing: {}", description);
        println!("   File: {:?}", path);
        println!("   Size: {:.1} MB", size_mb);
        println!("   Expected: {} mmap", if size_mb >= 50.0 { "USE" } else { "SKIP" });

        // Baseline: Standard I/O (no mmap)
        group.bench_with_input(
            BenchmarkId::new("standard_io", format!("{:.0}MB", size_mb)),
            &path,
            |b, path| {
                b.iter(|| {
                    let file = File::open(path).unwrap();
                    let mut reader = BufReader::new(file);
                    let mut output = Vec::new();
                    reader.read_to_end(&mut output).unwrap();
                    black_box(output.len())
                });
            },
        );

        // Rule 4: Smart mmap (threshold-based)
        group.bench_with_input(
            BenchmarkId::new("smart_mmap", format!("{:.0}MB", size_mb)),
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

/// Test decompression with mmap vs standard I/O
fn bench_decompression_with_mmap(c: &mut Criterion) {
    let mut group = c.benchmark_group("mmap_decompression");
    group.sample_size(10);

    // Focus on large file where mmap should shine
    let path = match get_dataset_path("huge_10m_150bp.fq.gz") {
        Some(p) => p,
        None => {
            eprintln!("⚠️  Skipping: huge_10m file not found");
            return;
        }
    };

    println!("\n📊 Decompression test: 544 MB file");
    println!("   Expected: 2.5× speedup from mmap on macOS");

    group.bench_function("with_mmap", |b| {
        b.iter(|| {
            let source = DataSource::from_path(&path);
            let mut reader = CompressedReader::new(source).unwrap();
            let mut output = Vec::new();
            reader.read_to_end(&mut output).unwrap();
            black_box(output.len())
        });
    });

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default();
    targets = bench_rule4_mmap_threshold, bench_decompression_with_mmap
}

criterion_main!(benches);
