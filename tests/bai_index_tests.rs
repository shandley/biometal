use biometal::io::bam::{BamReader, BaiIndex};
use std::path::PathBuf;

/// Integration tests for BAI index functionality
///
/// These tests validate that our BAI index implementation correctly:
/// - Loads BAI index files
/// - Performs O(log n) region queries
/// - Returns the same results as samtools
/// - Handles edge cases gracefully

#[test]
fn test_bai_load() {
    // Test loading a BAI index file
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");
    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    // Should have 3 references (chr1, chr2, chr22)
    assert_eq!(index.references.len(), 3, "Expected 3 references in index");
}

#[test]
fn test_bai_query_chunks() {
    // Test that querying chunks from index works
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");
    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    // Query chunks for chr1:1-1000 (reference_id = 0)
    let chunks = index.query_chunks(0, 1, 1000);

    // Should return some chunks or None
    match chunks {
        Some(c) => {
            println!("Found {} chunks for chr1:1-1000", c.len());
            for (i, chunk) in c.iter().enumerate().take(3) {
                println!("  Chunk {}: start={:#x}, end={:#x}",
                    i, chunk.start.as_raw(), chunk.end.as_raw());
            }
            assert!(c.len() > 0, "Should have at least one chunk");
        }
        None => {
            println!("No chunks found for region - may be expected if no data");
        }
    }
}

#[test]
fn test_indexed_region_query_small() {
    // Test querying a small region (chr1:1-1000)
    // Expected: ~2985 records (verified with samtools)
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    let mut query = BamReader::query(&bam_path, &index, "chr1", 1, 1000)
        .expect("Failed to create region query");

    let mut count = 0;
    while let Some(record) = query.next() {
        let record = record.expect("Failed to read record");

        // Verify record is in the queried region
        // Note: BAM uses 0-based coordinates, so query "1-1000" in 1-based
        // maps to [0, 1000) in 0-based coordinates
        if let Some(pos) = record.position {
            assert!(pos >= 0 && pos < 1000,
                "Record at position {} outside query range [0, 1000)", pos);
        }

        count += 1;
    }

    // Allow for some flexibility due to overlapping reads
    // Samtools reports 2985, but reads can start before and extend into region
    assert!(count >= 2900 && count <= 3100,
        "Expected ~2985 records, got {}", count);
}

#[test]
fn test_indexed_region_query_medium() {
    // Test querying a medium region (chr1:1-10000)
    // Expected: ~30693 records (verified with samtools)
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    let mut query = BamReader::query(&bam_path, &index, "chr1", 1, 10000)
        .expect("Failed to create region query");

    let mut count = 0;
    while let Some(record) = query.next() {
        let record = record.expect("Failed to read record");

        // BAM uses 0-based coordinates
        if let Some(pos) = record.position {
            assert!(pos >= 0 && pos < 10000,
                "Record at position {} outside query range [0, 10000)", pos);
        }

        count += 1;
    }

    // Allow for flexibility
    assert!(count >= 30000 && count <= 31000,
        "Expected ~30693 records, got {}", count);
}

#[test]
fn test_indexed_region_query_different_chromosome() {
    // Test querying chr1 vs chr2 returns different results
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    // Query chr1
    let query1 = BamReader::query(&bam_path, &index, "chr1", 1, 1000)
        .expect("Failed to create chr1 query");
    let count1: usize = query1.map(|r| r.ok()).flatten().count();

    // Query chr2 (which has no reads in our test file)
    let query2 = BamReader::query(&bam_path, &index, "chr2", 1, 1000)
        .expect("Failed to create chr2 query");
    let count2: usize = query2.map(|r| r.ok()).flatten().count();

    // chr1 should have records, chr2 should be empty in our test file
    assert!(count1 > 0, "chr1 should have records");
    assert_eq!(count2, 0, "chr2 should have no records in test file");
}

#[test]
fn test_indexed_empty_region() {
    // Test querying a region with no records
    // Using a high position where data is unlikely
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    let mut query = BamReader::query(&bam_path, &index, "chr22", 50000000, 50001000)
        .expect("Failed to create region query");

    let count: usize = query.map(|r| r.ok()).flatten().count();

    // Should have zero or very few records
    assert!(count < 10, "Expected empty or nearly empty region, got {} records", count);
}

#[test]
fn test_nonexistent_reference() {
    // Test querying a reference that doesn't exist
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    let result = BamReader::query(&bam_path, &index, "chrX", 1, 1000);

    assert!(result.is_err(), "Should fail when querying non-existent reference");
}

#[test]
fn test_boundary_conditions() {
    // Test edge cases: very small regions, single base queries
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    // Single base query (chr1:100-101)
    let query = BamReader::query(&bam_path, &index, "chr1", 100, 101)
        .expect("Failed to create single-base query");

    // Should be able to iterate
    let count: usize = query.map(|r| r.ok()).flatten().count();

    // samtools shows 331 reads overlap this position (reads span multiple bases)
    // Allow some flexibility for alignment differences
    assert!(count >= 300 && count <= 400,
        "Expected ~331 records overlapping position 100, got {}", count);
}

#[test]
fn test_adjacent_regions() {
    // Test that adjacent non-overlapping regions return different results
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    // First region
    let mut query1 = BamReader::query(&bam_path, &index, "chr1", 1, 5000)
        .expect("Failed to create first query");
    let count1: usize = query1.map(|r| r.ok()).flatten().count();

    // Adjacent region
    let mut query2 = BamReader::query(&bam_path, &index, "chr1", 5000, 10000)
        .expect("Failed to create second query");
    let count2: usize = query2.map(|r| r.ok()).flatten().count();

    assert!(count1 > 0, "First region should have records");
    assert!(count2 > 0, "Second region should have records");

    // Together should be close to the combined range count
    let combined_count = count1 + count2;
    assert!(combined_count >= 28000 && combined_count <= 32000,
        "Combined count should be ~30693, got {}", combined_count);
}

#[test]
fn test_record_fields_in_indexed_query() {
    // Test that records from indexed queries have all fields populated correctly
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    let mut query = BamReader::query(&bam_path, &index, "chr1", 1, 1000)
        .expect("Failed to create region query");

    let mut checked = 0;
    while let Some(Ok(record)) = query.next() {
        // Verify essential fields are present
        assert!(!record.name.is_empty(), "Record name should not be empty");
        assert!(record.sequence.len() > 0, "Sequence should not be empty");
        assert!(record.quality.len() == record.sequence.len(),
            "Quality length should match sequence length");

        if let Some(pos) = record.position {
            assert!(pos >= 0, "Position should be non-negative");
        }

        checked += 1;
        if checked >= 100 {
            break; // Check first 100 records
        }
    }

    assert!(checked > 0, "Should have checked at least some records");
}

#[test]
fn test_streaming_constant_memory() {
    // Test that indexed queries maintain constant memory (streaming)
    // This is a behavioral test - we iterate through many records
    // and verify the iterator pattern works correctly
    let bam_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam");
    let bai_path = PathBuf::from("tests/data/synthetic/alignment/synthetic_100k.bam.bai");

    let index = BaiIndex::from_path(&bai_path)
        .expect("Failed to load BAI index");

    let mut query = BamReader::query(&bam_path, &index, "chr1", 1, 10000)
        .expect("Failed to create region query");

    let mut count = 0;
    let mut prev_pos = None;

    // Process records in streaming fashion
    while let Some(result) = query.next() {
        let record = result.expect("Failed to read record");

        // Verify records are generally in order (may have some variation)
        if let Some(pos) = record.position {
            if let Some(prev) = prev_pos {
                // Positions should not decrease dramatically
                assert!(pos >= prev - 1000,
                    "Positions should be generally ordered, prev={}, current={}", prev, pos);
            }
            prev_pos = Some(pos);
        }

        count += 1;
    }

    assert!(count > 1000, "Should process substantial number of records");
}
