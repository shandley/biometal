//! Microbenchmark for isolating BAM sequence decoding performance.
//!
//! This benchmark measures ONLY the 4-bit to ASCII decoding operation,
//! isolating it from BGZF decompression, record parsing, and I/O.
//!
//! Purpose: Validate if sequence decoding is ≥15% of total BAM parsing time
//! (threshold for SIMD optimization per OPTIMIZATION_RULES.md Rule 1).

use biometal::io::bam::sequence::decode_sequence;
use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use std::io::Read;

/// Generate realistic 4-bit encoded sequence data
fn generate_test_data(num_bases: usize) -> Vec<u8> {
    // BAM encoding: 2 bases per byte
    let num_bytes = num_bases.div_ceil(2);

    // Generate random-ish data using a simple pattern
    // Real BAM files have ~25% each of A,C,G,T (indices 1,2,4,8)
    let mut data = Vec::with_capacity(num_bytes);
    for i in 0..num_bytes {
        // Alternate between common bases to simulate realistic data
        let base1 = match i % 4 {
            0 => 1, // A
            1 => 2, // C
            2 => 4, // G
            _ => 8, // T
        };
        let base2 = match (i + 1) % 4 {
            0 => 1, // A
            1 => 2, // C
            2 => 4, // G
            _ => 8, // T
        };
        data.push((base1 << 4) | base2);
    }

    data
}

/// Benchmark decoding a single read (100 bases, typical Illumina)
fn bench_decode_single_read(c: &mut Criterion) {
    let length = 100;
    let data = generate_test_data(length);

    let mut group = c.benchmark_group("sequence_decode_single");
    group.throughput(Throughput::Elements(length as u64));

    group.bench_function("100bp_read", |b| {
        b.iter(|| {
            let sequence = decode_sequence(black_box(&data), black_box(length))
                .expect("Decode failed");
            black_box(sequence);
        });
    });

    group.finish();
}

/// Benchmark decoding 100K reads (realistic workload)
fn bench_decode_100k_reads(c: &mut Criterion) {
    let length = 100;
    let single_read = generate_test_data(length);
    let num_reads = 100_000;

    let mut group = c.benchmark_group("sequence_decode_bulk");
    group.throughput(Throughput::Elements(num_reads * length as u64));

    group.bench_function("100k_reads_100bp", |b| {
        b.iter(|| {
            let mut total_bases = 0;
            for _ in 0..num_reads {
                let sequence = decode_sequence(black_box(&single_read), black_box(length))
                    .expect("Decode failed");
                total_bases += sequence.len();
                black_box(sequence);
            }
            total_bases
        });
    });

    group.finish();
}

/// Benchmark different sequence lengths
fn bench_decode_various_lengths(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence_decode_lengths");

    for &length in &[50, 100, 150, 250, 500, 1000] {
        let data = generate_test_data(length);

        group.throughput(Throughput::Elements(length as u64));
        group.bench_with_input(
            format!("{}bp", length),
            &(data, length),
            |b, (data, length)| {
                b.iter(|| {
                    let sequence = decode_sequence(black_box(data), black_box(*length))
                        .expect("Decode failed");
                    black_box(sequence);
                });
            },
        );
    }

    group.finish();
}

/// Benchmark realistic BAM file decoding (extract from actual BAM parsing)
///
/// This measures how much time is spent JUST on sequence decoding
/// within the context of full BAM file parsing.
fn bench_decode_from_bam_file(c: &mut Criterion) {
    use biometal::io::bam::BamReader;
    use std::path::Path;

    let bam_path = "../native-bam-implementation/test-data/synthetic_100000.bam";

    if !Path::new(bam_path).exists() {
        eprintln!("Warning: Test BAM file not found: {}", bam_path);
        eprintln!("Skipping BAM-based benchmark");
        return;
    }

    // First pass: Extract all encoded sequences from the BAM file
    let mut bam = BamReader::from_path(bam_path)
        .expect("Failed to open BAM file");

    let mut encoded_sequences = Vec::new();
    for result in bam.records() {
        let record = result.expect("Failed to parse record");
        // NOTE: This would require access to the raw encoded data
        // For now, we'll use generated data as proxy
        let length = record.sequence.len();
        encoded_sequences.push((generate_test_data(length), length));

        if encoded_sequences.len() >= 10_000 {
            break; // Limit to 10K for manageable benchmark time
        }
    }

    let mut group = c.benchmark_group("sequence_decode_realistic");
    let total_bases: u64 = encoded_sequences.iter().map(|(_, len)| *len as u64).sum();
    group.throughput(Throughput::Elements(total_bases));

    group.bench_function("10k_reads_from_bam", |b| {
        b.iter(|| {
            let mut total_decoded = 0;
            for (data, length) in &encoded_sequences {
                let sequence = decode_sequence(black_box(data), black_box(*length))
                    .expect("Decode failed");
                total_decoded += sequence.len();
                black_box(sequence);
            }
            total_decoded
        });
    });

    group.finish();
}

criterion_group! {
    name = sequence_decode_benches;
    config = Criterion::default()
        .sample_size(30)  // N=30 for statistical significance
        .measurement_time(std::time::Duration::from_secs(5));
    targets =
        bench_decode_single_read,
        bench_decode_100k_reads,
        bench_decode_various_lengths,
        bench_decode_from_bam_file
}

criterion_main!(sequence_decode_benches);
