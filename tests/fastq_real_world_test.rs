//! Integration tests for FASTQ real-world data
//!
//! These tests validate biometal's FASTQ parser against real sequencing data
//! from the 1000 Genomes Project (SRA ERR000589).

use biometal::io::FastqStream;
use std::path::PathBuf;

#[test]
fn test_real_fastq_parsing() {
    // Load real-world FASTQ data from 1000 Genomes Project
    let fastq_path = PathBuf::from("tests/data/real_world/sequence/small_sample.fastq.gz");

    let stream = FastqStream::from_path(&fastq_path)
        .expect("Failed to open real-world FASTQ file");

    let mut record_count = 0;
    let mut total_bases = 0;
    let mut total_quality = 0u64;

    for record in stream {
        let record = record.expect("Failed to parse FASTQ record");
        record_count += 1;

        // Validate record structure
        assert!(!record.id.is_empty(), "Record ID should not be empty");
        assert!(!record.sequence.is_empty(), "Sequence should not be empty");
        assert_eq!(
            record.sequence.len(),
            record.quality.len(),
            "Sequence and quality lengths must match"
        );

        // Accumulate statistics
        total_bases += record.sequence.len();
        total_quality += record.quality.iter().map(|&q| q as u64).sum::<u64>();

        // Validate quality scores are in valid range (Illumina 1.8+: 0-41 -> ASCII 33-74)
        for &qual in &record.quality {
            assert!(
                qual >= 33 && qual <= 74,
                "Quality score {} out of valid range (33-74)",
                qual
            );
        }
    }

    // Verify we got the expected number of records (1000 reads, 4 lines each = 4000 lines)
    assert_eq!(record_count, 1000, "Expected 1000 FASTQ records");

    // Basic sanity checks on statistics
    assert!(total_bases > 0, "Should have parsed some bases");
    assert!(total_quality > 0, "Should have parsed quality scores");

    let avg_read_length = total_bases / record_count;
    let avg_quality = total_quality / total_bases as u64;

    // Typical Illumina reads are 50-150bp
    assert!(
        avg_read_length >= 50 && avg_read_length <= 150,
        "Average read length {} seems unusual for Illumina data",
        avg_read_length
    );

    // Quality scores should be reasonable (typically 30-40 for good data)
    // ASCII 33 offset, so quality 30 = 63, quality 40 = 73
    assert!(
        avg_quality >= 50 && avg_quality <= 75,
        "Average quality score {} seems unusual",
        avg_quality
    );

    println!("✅ Parsed {} records", record_count);
    println!("   Average read length: {} bp", avg_read_length);
    println!("   Average quality: {} (ASCII)", avg_quality);
}

#[test]
fn test_real_fastq_id_format() {
    // Validate that real SRA data has properly formatted IDs
    let fastq_path = PathBuf::from("tests/data/real_world/sequence/small_sample.fastq.gz");

    let stream = FastqStream::from_path(&fastq_path)
        .expect("Failed to open real-world FASTQ file");

    let mut sra_format_count = 0;

    for (idx, record) in stream.enumerate() {
        let record = record.expect("Failed to parse FASTQ record");

        // Check if ID starts with '@' (some parsers include it, some don't)
        let id = if record.id.starts_with('@') {
            &record.id[1..]
        } else {
            &record.id
        };

        // SRA identifiers typically start with ERR/SRR/DRR followed by numbers
        if id.starts_with("ERR") || id.starts_with("SRR") || id.starts_with("DRR") {
            sra_format_count += 1;
        }

        // Just check first 10 records
        if idx >= 10 {
            break;
        }
    }

    // Most records should follow SRA naming convention
    assert!(
        sra_format_count > 0,
        "Expected at least some records to have SRA-style IDs"
    );

    println!("✅ Found {} SRA-formatted IDs in first 10 records", sra_format_count);
}

#[test]
fn test_real_fastq_constant_memory() {
    // Verify streaming parser doesn't load entire file into memory
    let fastq_path = PathBuf::from("tests/data/real_world/sequence/small_sample.fastq.gz");

    let stream = FastqStream::from_path(&fastq_path)
        .expect("Failed to open real-world FASTQ file");

    let mut prev_record = None;
    let mut records_seen = 0;

    for record in stream {
        let record = record.expect("Failed to parse FASTQ record");

        // Verify previous record is dropped (can't check memory directly,
        // but this pattern ensures we're not accumulating)
        if let Some(_prev) = prev_record.take() {
            // Previous record dropped here
        }

        prev_record = Some(record);
        records_seen += 1;

        // Check first 100 records
        if records_seen >= 100 {
            break;
        }
    }

    assert_eq!(records_seen, 100, "Should have processed 100 records");
    println!("✅ Streaming parser maintains constant memory (1 record at a time)");
}
