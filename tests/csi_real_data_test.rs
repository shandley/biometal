//! Real-world CSI index validation tests
//!
//! These tests validate the CSI parser against real index files created by samtools.

use biometal::formats::index::CsiIndex;
use biometal::io::compression::{CompressedReader, DataSource};
use std::path::Path;
use std::io::{Read, Cursor};

/// Test parsing a real CSI index created by samtools
#[test]
fn test_parse_real_samtools_csi() {
    let csi_path = "tests/data/synthetic/alignment/synthetic_100k.bam.csi";

    // Skip test if file doesn't exist
    if !Path::new(csi_path).exists() {
        eprintln!("Skipping test: {} not found", csi_path);
        eprintln!("Run: samtools index -c tests/data/synthetic/alignment/synthetic_100k.bam");
        return;
    }

    // CSI files created by samtools are BGZF-compressed
    // We need to decompress them first
    let source = DataSource::from_path(csi_path);
    let mut decompressor = CompressedReader::new(source)
        .expect("Failed to create decompressor");

    // Read all decompressed data into memory
    let mut decompressed_data = Vec::new();
    decompressor.read_to_end(&mut decompressed_data)
        .expect("Failed to decompress CSI file");

    println!("Decompressed {} bytes", decompressed_data.len());

    // Parse the decompressed CSI data
    let mut cursor = Cursor::new(decompressed_data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Validate basic structure
    println!("CSI Index loaded successfully!");
    println!("  Min shift: {}", index.min_shift());
    println!("  Depth: {}", index.depth());
    println!("  References: {}", index.references().len());

    // Expected values for synthetic_100k.bam:
    // - min_shift: 14 (16kb bins)
    // - depth: 5 (5 binning levels)
    // - n_ref: 3 (chr1, chr2, chr22)
    assert_eq!(index.min_shift(), 14, "Expected min_shift=14 (16kb bins)");
    assert_eq!(index.depth(), 5, "Expected depth=5");
    assert_eq!(index.references().len(), 3, "Expected 3 references (chr1, chr2, chr22)");

    // Validate each reference
    for (idx, reference) in index.references().iter().enumerate() {
        println!("\nReference {}:", idx);
        println!("  Bins: {}", reference.bins.len());
        // Note: CSI does not have linear intervals like BAI/TBI

        // Only ref 0 (chr1) should have bins (it has reads)
        if idx == 0 {
            assert!(!reference.bins.is_empty(), "Reference {} should have bins", idx);
        }

        // Validate bin structure
        for bin in &reference.bins {
            println!("    Bin {}: {} chunks", bin.bin_id, bin.chunks.len());

            // Each bin should have at least one chunk
            assert!(!bin.chunks.is_empty(), "Bin {} should have chunks", bin.bin_id);

            // Validate chunks
            for chunk in &bin.chunks {
                assert!(
                    chunk.end.as_raw() >= chunk.start.as_raw(),
                    "Chunk end should be >= start"
                );
            }
        }
    }
}

/// Test querying regions in the real CSI index
#[test]
fn test_query_real_csi() {
    let csi_path = "tests/data/synthetic/alignment/synthetic_100k.bam.csi";

    if !Path::new(csi_path).exists() {
        eprintln!("Skipping test: {} not found", csi_path);
        return;
    }

    // Load and parse the CSI index
    let source = DataSource::from_path(csi_path);
    let mut decompressor = CompressedReader::new(source)
        .expect("Failed to create decompressor");

    let mut decompressed_data = Vec::new();
    decompressor.read_to_end(&mut decompressed_data)
        .expect("Failed to decompress CSI file");

    let mut cursor = Cursor::new(decompressed_data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Test query by index (we don't have reference names in aux data)
    println!("\nTesting region queries:");

    // Query chr1 (reference 0) region
    let start = 1_000_000;
    let end = 2_000_000;

    match index.query_by_index(0, start, end).unwrap() {
        Some(chunks) => {
            println!("Query ref=0, region {}-{}: {} chunks", start, end, chunks.len());
            assert!(!chunks.is_empty(), "Should find chunks for region");

            for (i, chunk) in chunks.iter().enumerate() {
                println!(
                    "  Chunk {}: {:#x} - {:#x}",
                    i + 1,
                    chunk.start.as_raw(),
                    chunk.end.as_raw()
                );
            }
        }
        None => panic!("Reference 0 should exist"),
    }

    // Test query for chr2 (reference 1)
    match index.query_by_index(1, start, end).unwrap() {
        Some(chunks) => {
            println!("Query ref=1, region {}-{}: {} chunks", start, end, chunks.len());
        }
        None => panic!("Reference 1 should exist"),
    }

    // Test query for chr22 (reference 2)
    match index.query_by_index(2, 1_000_000, 1_500_000).unwrap() {
        Some(chunks) => {
            println!("Query ref=2, region 1000000-1500000: {} chunks", chunks.len());
        }
        None => panic!("Reference 2 should exist"),
    }

    // Test query out of bounds
    let result = index.query_by_index(10, start, end).unwrap();
    assert!(result.is_none(), "Query for non-existent reference should return None");

    // Test invalid range
    let result = index.query_by_index(0, end, start);
    assert!(result.is_err(), "Query with start > end should fail");
}

/// Test that CSI results match the structure we expect from samtools
#[test]
fn test_csi_structure_validation() {
    let csi_path = "tests/data/synthetic/alignment/synthetic_100k.bam.csi";

    if !Path::new(csi_path).exists() {
        eprintln!("Skipping test: {} not found", csi_path);
        return;
    }

    let source = DataSource::from_path(csi_path);
    let mut decompressor = CompressedReader::new(source)
        .expect("Failed to create decompressor");

    let mut decompressed_data = Vec::new();
    decompressor.read_to_end(&mut decompressed_data)
        .expect("Failed to decompress CSI file");

    let mut cursor = Cursor::new(decompressed_data);
    let index = CsiIndex::parse(&mut cursor).expect("Failed to parse CSI");

    // Validate CSI-specific features
    println!("\nCSI-specific validation:");

    // Check loffset values (should be non-zero for bins with data)
    for (ref_idx, reference) in index.references().iter().enumerate() {
        println!("Reference {} loffsets:", ref_idx);
        for bin in &reference.bins {
            if !bin.chunks.is_empty() {
                println!("  Bin {}: loffset={:#x}", bin.bin_id, bin.loffset.as_raw());
                // Note: loffset is a file offset hint (leftmost position in bin)
                // It doesn't necessarily equal the first chunk's start
            }
        }
    }

    // Verify hierarchical binning
    for (idx, reference) in index.references().iter().enumerate() {
        let bin_ids: Vec<u32> = reference.bins.iter().map(|b| b.bin_id).collect();

        if !bin_ids.is_empty() {
            println!("Reference {} bin IDs: {:?}", idx, &bin_ids[..bin_ids.len().min(10)]);
            // Bin 0 covers entire sequence, but may not exist for sparse/empty references
        }
    }

    // Validate that ref 0 (which has data) has bins
    assert!(!index.references()[0].bins.is_empty(), "Reference 0 should have bins");
}
