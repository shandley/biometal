use biometal::Result;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufReader, Read, Write};

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

fn decompress_file(path: &str) -> Result<Vec<u8>> {
    let file = File::open(path)?;
    let mut decoder = flate2::read::MultiGzDecoder::new(BufReader::new(file));
    let mut uncompressed = Vec::new();
    decoder.read_to_end(&mut uncompressed)?;
    Ok(uncompressed)
}

fn bench_compression(c: &mut Criterion) {
    let mut group = c.benchmark_group("compression_backend");
    group.sample_size(10); // N=10 samples for consistency with decompression benchmark

    // Print backend info
    #[cfg(feature = "zlib-ng")]
    println!("   Backend: zlib-ng (ARM NEON optimized)");
    #[cfg(not(feature = "zlib-ng"))]
    println!("   Backend: rust_backend (miniz_oxide)");

    for (path, size_label, _expected_size) in TEST_FILES {
        println!("\n📊 Testing: {}", size_label);
        println!("   File: \"{}\"", path);

        // Decompress once to get uncompressed data for benchmarking
        let uncompressed_data = match decompress_file(path) {
            Ok(data) => data,
            Err(e) => {
                eprintln!("Failed to decompress {}: {}", path, e);
                continue;
            }
        };

        println!(
            "   Uncompressed size: {:.1} MB",
            uncompressed_data.len() as f64 / 1_000_000.0
        );

        // Benchmark compression with default level
        group.bench_with_input(
            BenchmarkId::new("compress_default", size_label),
            &uncompressed_data,
            |b, data| {
                b.iter(|| {
                    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
                    encoder.write_all(data).unwrap();
                    let compressed = encoder.finish().unwrap();
                    std::hint::black_box(compressed.len())
                })
            },
        );

        // Benchmark compression with fast level
        group.bench_with_input(
            BenchmarkId::new("compress_fast", size_label),
            &uncompressed_data,
            |b, data| {
                b.iter(|| {
                    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
                    encoder.write_all(data).unwrap();
                    let compressed = encoder.finish().unwrap();
                    std::hint::black_box(compressed.len())
                })
            },
        );

        // Benchmark compression with best level
        group.bench_with_input(
            BenchmarkId::new("compress_best", size_label),
            &uncompressed_data,
            |b, data| {
                b.iter(|| {
                    let mut encoder = GzEncoder::new(Vec::new(), Compression::best());
                    encoder.write_all(data).unwrap();
                    let compressed = encoder.finish().unwrap();
                    std::hint::black_box(compressed.len())
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_compression);
criterion_main!(benches);
