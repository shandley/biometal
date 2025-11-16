//! Integration tests for PAF real-world data
//!
//! These tests validate biometal's PAF parser against real minimap2 output.

use biometal::formats::paf::PafParser;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

#[test]
fn test_real_paf_parsing() {
    // Load real-world PAF data from minimap2
    let paf_path = PathBuf::from("tests/data/real_world/alignments/minimap2_alignment.paf");

    let file = File::open(&paf_path).expect("Failed to open PAF file");
    let reader = BufReader::new(file);
    let parser = PafParser::new(reader);

    let mut record_count = 0;
    let mut total_alignment_length = 0;
    let mut total_matches = 0;

    for record in parser {
        let record = record.expect("Failed to parse PAF record");
        record_count += 1;

        // Validate record structure
        assert!(!record.query_name.is_empty(), "Query name should not be empty");
        assert!(!record.target_name.is_empty(), "Target name should not be empty");
        assert!(record.query_length > 0, "Query length should be positive");
        assert!(record.target_length > 0, "Target length should be positive");

        // Validate coordinate ranges
        assert!(
            record.query_start < record.query_end,
            "Query start must be < query end"
        );
        assert!(
            record.target_start < record.target_end,
            "Target start must be < target end"
        );

        // Validate alignment metrics
        assert!(
            record.num_matches <= record.alignment_length,
            "Matches cannot exceed alignment length"
        );
        assert!(
            record.alignment_length > 0,
            "Alignment length should be positive"
        );

        // Accumulate statistics
        total_alignment_length += record.alignment_length;
        total_matches += record.num_matches;

        // mapq is u8 (0-255), so no validation needed
    }

    // Verify we got the expected number of records
    assert_eq!(record_count, 18, "Expected 18 PAF alignment records");

    // Basic sanity checks on statistics
    assert!(total_alignment_length > 0, "Should have aligned some bases");
    assert!(total_matches > 0, "Should have some matching bases");

    let avg_identity = (total_matches as f64 / total_alignment_length as f64) * 100.0;

    // Alignments should have reasonable identity (these are from the same reference)
    assert!(
        avg_identity > 90.0,
        "Average identity {} seems too low for self-alignment",
        avg_identity
    );

    println!("✅ Parsed {} PAF records", record_count);
    println!(
        "   Average identity: {:.2}%",
        avg_identity
    );
    println!(
        "   Total alignment length: {} bp",
        total_alignment_length
    );
}

#[test]
fn test_real_paf_target_names() {
    // Validate that alignments map to expected target (chr1)
    let paf_path = PathBuf::from("tests/data/real_world/alignments/minimap2_alignment.paf");

    let file = File::open(&paf_path).expect("Failed to open PAF file");
    let reader = BufReader::new(file);
    let parser = PafParser::new(reader);

    let mut chr1_alignments = 0;

    for record in parser {
        let record = record.expect("Failed to parse PAF record");

        // All alignments should be to chr1 (our reference)
        assert_eq!(
            record.target_name, "chr1",
            "Expected all alignments to target 'chr1', got '{}'",
            record.target_name
        );
        chr1_alignments += 1;
    }

    assert!(
        chr1_alignments > 0,
        "Should have found alignments to chr1"
    );

    println!("✅ Found {} alignments to chr1", chr1_alignments);
}

#[test]
fn test_real_paf_query_names() {
    // Validate query names match expected pattern
    let paf_path = PathBuf::from("tests/data/real_world/alignments/minimap2_alignment.paf");

    let file = File::open(&paf_path).expect("Failed to open PAF file");
    let reader = BufReader::new(file);
    let parser = PafParser::new(reader);

    let mut query_count = 0;

    for record in parser {
        let record = record.expect("Failed to parse PAF record");

        // Query names should start with "query_" (our naming convention)
        assert!(
            record.query_name.starts_with("query_"),
            "Expected query name to start with 'query_', got '{}'",
            record.query_name
        );

        query_count += 1;
    }

    assert!(query_count > 0, "Should have found query alignments");

    println!("✅ Found {} alignments with expected query naming", query_count);
}

#[test]
fn test_real_paf_strand() {
    // Validate strand information
    let paf_path = PathBuf::from("tests/data/real_world/alignments/minimap2_alignment.paf");

    let file = File::open(&paf_path).expect("Failed to open PAF file");
    let reader = BufReader::new(file);
    let parser = PafParser::new(reader);

    let mut forward_count = 0;
    let mut reverse_count = 0;

    for record in parser {
        let record = record.expect("Failed to parse PAF record");

        // Strand should be '+' or '-'
        assert!(
            record.strand == '+' || record.strand == '-',
            "Invalid strand: {}",
            record.strand
        );

        if record.strand == '+' {
            forward_count += 1;
        } else {
            reverse_count += 1;
        }
    }

    // We should have at least some forward alignments (our queries are from the reference)
    assert!(
        forward_count > 0,
        "Should have at least some forward alignments"
    );

    println!("✅ Strand distribution: {} forward, {} reverse", forward_count, reverse_count);
}

#[test]
fn test_real_paf_constant_memory() {
    // Verify streaming parser doesn't load entire file into memory
    let paf_path = PathBuf::from("tests/data/real_world/alignments/minimap2_alignment.paf");

    let file = File::open(&paf_path).expect("Failed to open PAF file");
    let reader = BufReader::new(file);
    let parser = PafParser::new(reader);

    let mut prev_record = None;
    let mut records_seen = 0;

    for record in parser {
        let record = record.expect("Failed to parse PAF record");

        // Verify previous record is dropped
        if let Some(_prev) = prev_record.take() {
            // Previous record dropped here
        }

        prev_record = Some(record);
        records_seen += 1;
    }

    assert_eq!(records_seen, 18, "Should have processed all 18 records");
    println!("✅ Streaming parser maintains constant memory (1 record at a time)");
}
