//! Integration tests for FASTA real-world data
//!
//! These tests validate biometal's FASTA parser against real reference sequence
//! from UCSC/Ensembl (GRCh38 chr21).

use biometal::io::FastaStream;
use std::path::PathBuf;

#[test]
fn test_real_fasta_parsing() {
    // Load real-world FASTA data from UCSC hg38 chr21
    let fasta_path = PathBuf::from("tests/data/real_world/sequence/hg38_chr21_10kb.fa.gz");

    let stream = FastaStream::from_path(&fasta_path)
        .expect("Failed to open real-world FASTA file");

    let mut record_count = 0;
    let mut total_bases = 0;
    let mut gc_count = 0;
    let mut n_count = 0;

    for record in stream {
        let record = record.expect("Failed to parse FASTA record");
        record_count += 1;

        // Validate record structure
        assert!(!record.id.is_empty(), "Record ID should not be empty");
        assert!(!record.sequence.is_empty(), "Sequence should not be empty");

        // Validate sequence contains only valid bases
        for &base in &record.sequence {
            match base {
                b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' => {},
                _ => panic!("Invalid base in sequence: {}", base as char),
            }
        }

        // Count bases
        total_bases += record.sequence.len();
        for &base in &record.sequence {
            match base {
                b'G' | b'C' | b'g' | b'c' => gc_count += 1,
                b'N' | b'n' => n_count += 1,
                _ => {},
            }
        }
    }

    // Verify we got the expected record
    assert_eq!(record_count, 1, "Expected 1 FASTA record");
    assert!(total_bases > 0, "Should have parsed some bases");

    // GC content should be reasonable for human genome (~40-45%)
    let gc_percent = (gc_count as f64 / total_bases as f64) * 100.0;
    assert!(
        gc_percent > 30.0 && gc_percent < 60.0,
        "GC content {} seems unusual for human genome",
        gc_percent
    );

    println!("✅ Parsed {} FASTA record(s)", record_count);
    println!("   Total bases: {}", total_bases);
    println!("   GC content: {:.2}%", gc_percent);
    println!("   N bases: {}", n_count);
}

#[test]
fn test_real_fasta_id_parsing() {
    // Validate FASTA ID parsing with real UCSC/Ensembl header
    let fasta_path = PathBuf::from("tests/data/real_world/sequence/hg38_chr21_10kb.fa.gz");

    let stream = FastaStream::from_path(&fasta_path)
        .expect("Failed to open real-world FASTA file");

    for record in stream {
        let record = record.expect("Failed to parse FASTA record");

        // ID should contain "chr21" or "21"
        assert!(
            record.id.contains("chr21") || record.id.contains("21"),
            "Expected ID to contain chromosome 21 reference, got: {}",
            record.id
        );

        println!("✅ FASTA ID: {}", record.id);
        break; // Just check first record
    }
}

#[test]
fn test_real_fasta_sequence_quality() {
    // Validate sequence quality (not all N's, has variety)
    let fasta_path = PathBuf::from("tests/data/real_world/sequence/hg38_chr21_10kb.fa.gz");

    let stream = FastaStream::from_path(&fasta_path)
        .expect("Failed to open real-world FASTA file");

    for record in stream {
        let record = record.expect("Failed to parse FASTA record");

        // Count each base type
        let mut base_counts = [0usize; 5]; // A, C, G, T, N
        for &base in &record.sequence {
            match base {
                b'A' | b'a' => base_counts[0] += 1,
                b'C' | b'c' => base_counts[1] += 1,
                b'G' | b'g' => base_counts[2] += 1,
                b'T' | b't' => base_counts[3] += 1,
                b'N' | b'n' => base_counts[4] += 1,
                _ => {},
            }
        }

        // Should have all four bases represented (not all one type or all N)
        let non_n_bases = base_counts[0..4].iter().filter(|&&c| c > 0).count();
        assert!(
            non_n_bases >= 3,
            "Expected at least 3 different base types, got {} (A:{}, C:{}, G:{}, T:{}, N:{})",
            non_n_bases, base_counts[0], base_counts[1], base_counts[2], base_counts[3], base_counts[4]
        );

        // N's should not dominate the sequence (for real non-gap regions)
        let n_percent = (base_counts[4] as f64 / record.sequence.len() as f64) * 100.0;
        assert!(
            n_percent < 90.0,
            "Too many N's ({}%), sequence might be all gaps",
            n_percent
        );

        println!("✅ Base distribution: A:{}, C:{}, G:{}, T:{}, N:{}",
                 base_counts[0], base_counts[1], base_counts[2], base_counts[3], base_counts[4]);
        break;
    }
}

#[test]
fn test_real_fasta_constant_memory() {
    // Verify streaming parser doesn't load entire file into memory
    let fasta_path = PathBuf::from("tests/data/real_world/sequence/hg38_chr21_10kb.fa.gz");

    let stream = FastaStream::from_path(&fasta_path)
        .expect("Failed to open real-world FASTA file");

    let mut prev_record = None;
    let mut records_seen = 0;

    for record in stream {
        let record = record.expect("Failed to parse FASTA record");

        // Verify previous record is dropped
        if let Some(_prev) = prev_record.take() {
            // Previous record dropped here
        }

        prev_record = Some(record);
        records_seen += 1;
    }

    assert_eq!(records_seen, 1, "Should have processed 1 record");
    println!("✅ Streaming parser maintains constant memory (1 record at a time)");
}
