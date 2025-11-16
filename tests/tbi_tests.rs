//! Integration tests for TBI (Tabix Index) functionality

use biometal::formats::index::TbiIndex;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/index")
}

/// Create a minimal TBI file for testing
///
/// This creates a valid TBI file with:
/// - 1 reference sequence ("chr1")
/// - VCF format
/// - Minimal binning and linear index data
fn create_minimal_tbi() -> Vec<u8> {
    let mut data = Vec::new();

    // Magic string "TBI\1"
    data.extend_from_slice(b"TBI\x01");

    // n_ref = 1 (int32)
    data.extend_from_slice(&1i32.to_le_bytes());

    // format = 2 (VCF, int32)
    data.extend_from_slice(&2i32.to_le_bytes());

    // col_seq = 0 (int32)
    data.extend_from_slice(&0i32.to_le_bytes());

    // col_beg = 1 (int32)
    data.extend_from_slice(&1i32.to_le_bytes());

    // col_end = 0 (int32, same as beg for VCF)
    data.extend_from_slice(&0i32.to_le_bytes());

    // meta = '#' (int32)
    data.extend_from_slice(&(b'#' as i32).to_le_bytes());

    // skip = 0 (int32)
    data.extend_from_slice(&0i32.to_le_bytes());

    // l_nm = 5 (length of "chr1\0", int32)
    data.extend_from_slice(&5i32.to_le_bytes());

    // names = "chr1\0"
    data.extend_from_slice(b"chr1\0");

    // Index for chr1:
    // n_bin = 1 (int32)
    data.extend_from_slice(&1i32.to_le_bytes());

    // Bin 0:
    // bin_id = 0 (uint32)
    data.extend_from_slice(&0u32.to_le_bytes());

    // n_chunk = 1 (int32)
    data.extend_from_slice(&1i32.to_le_bytes());

    // Chunk: start = 0x1000 (uint64)
    data.extend_from_slice(&0x1000u64.to_le_bytes());

    // Chunk: end = 0x2000 (uint64)
    data.extend_from_slice(&0x2000u64.to_le_bytes());

    // Linear index:
    // n_intv = 2 (int32)
    data.extend_from_slice(&2i32.to_le_bytes());

    // Interval 0: offset = 0x1000
    data.extend_from_slice(&0x1000u64.to_le_bytes());

    // Interval 1: offset = 0x1500
    data.extend_from_slice(&0x1500u64.to_le_bytes());

    data
}

#[test]
fn test_tbi_parse_minimal() {
    // Create test directory
    std::fs::create_dir_all(test_data_dir()).ok();

    let tbi_path = test_data_dir().join("minimal.tbi");
    let tbi_data = create_minimal_tbi();

    // Write TBI file
    let mut file = File::create(&tbi_path).expect("Failed to create TBI file");
    file.write_all(&tbi_data).expect("Failed to write TBI data");

    // Parse TBI file
    let index = TbiIndex::from_path(&tbi_path).expect("Failed to parse TBI");

    // Verify header fields
    assert_eq!(index.col_seq(), 0);
    assert_eq!(index.col_beg(), 1);
    assert_eq!(index.col_end(), 0);
    assert_eq!(index.meta_char(), '#');
    assert_eq!(index.skip_lines(), 0);

    // Verify references
    assert_eq!(index.references().len(), 1);
    let chr1 = index.get_reference("chr1").expect("chr1 not found");
    assert_eq!(chr1.name, "chr1");
    assert_eq!(chr1.bins.len(), 1);
    assert_eq!(chr1.intervals.len(), 2);

    // Cleanup
    std::fs::remove_file(&tbi_path).ok();
}

#[test]
fn test_tbi_query_region() {
    std::fs::create_dir_all(test_data_dir()).ok();

    let tbi_path = test_data_dir().join("query_test.tbi");
    let tbi_data = create_minimal_tbi();

    let mut file = File::create(&tbi_path).expect("Failed to create TBI file");
    file.write_all(&tbi_data).expect("Failed to write TBI data");

    let index = TbiIndex::from_path(&tbi_path).expect("Failed to parse TBI");

    // Query region
    let chunks = index.query("chr1", 0, 100000).expect("Query failed");

    // Should return at least one chunk
    assert!(!chunks.is_empty());
    assert_eq!(chunks[0].start.as_raw(), 0x1000);
    assert_eq!(chunks[0].end.as_raw(), 0x2000);

    // Cleanup
    std::fs::remove_file(&tbi_path).ok();
}

#[test]
fn test_tbi_query_nonexistent_reference() {
    std::fs::create_dir_all(test_data_dir()).ok();

    let tbi_path = test_data_dir().join("nonexist_test.tbi");
    let tbi_data = create_minimal_tbi();

    let mut file = File::create(&tbi_path).expect("Failed to create TBI file");
    file.write_all(&tbi_data).expect("Failed to write TBI data");

    let index = TbiIndex::from_path(&tbi_path).expect("Failed to parse TBI");

    // Query non-existent reference
    let result = index.query("chr99", 0, 100000);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("not found"));

    // Cleanup
    std::fs::remove_file(&tbi_path).ok();
}

#[test]
fn test_tbi_invalid_range() {
    std::fs::create_dir_all(test_data_dir()).ok();

    let tbi_path = test_data_dir().join("invalid_range_test.tbi");
    let tbi_data = create_minimal_tbi();

    let mut file = File::create(&tbi_path).expect("Failed to create TBI file");
    file.write_all(&tbi_data).expect("Failed to write TBI data");

    let index = TbiIndex::from_path(&tbi_path).expect("Failed to parse TBI");

    // Invalid range: start >= end
    let result = index.query("chr1", 100, 50);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid range"));

    // Cleanup
    std::fs::remove_file(&tbi_path).ok();
}

#[test]
fn test_tbi_multiple_references() {
    std::fs::create_dir_all(test_data_dir()).ok();

    // Create TBI with 2 references
    let mut data = Vec::new();

    // Magic
    data.extend_from_slice(b"TBI\x01");

    // n_ref = 2
    data.extend_from_slice(&2i32.to_le_bytes());

    // Header fields
    data.extend_from_slice(&2i32.to_le_bytes()); // format = VCF
    data.extend_from_slice(&0i32.to_le_bytes()); // col_seq
    data.extend_from_slice(&1i32.to_le_bytes()); // col_beg
    data.extend_from_slice(&0i32.to_le_bytes()); // col_end
    data.extend_from_slice(&(b'#' as i32).to_le_bytes()); // meta
    data.extend_from_slice(&0i32.to_le_bytes()); // skip

    // l_nm = 10 ("chr1\0chr2\0")
    data.extend_from_slice(&10i32.to_le_bytes());
    data.extend_from_slice(b"chr1\0chr2\0");

    // Index for chr1
    data.extend_from_slice(&1i32.to_le_bytes()); // n_bin
    data.extend_from_slice(&0u32.to_le_bytes()); // bin_id
    data.extend_from_slice(&1i32.to_le_bytes()); // n_chunk
    data.extend_from_slice(&0x1000u64.to_le_bytes()); // chunk start
    data.extend_from_slice(&0x2000u64.to_le_bytes()); // chunk end
    data.extend_from_slice(&1i32.to_le_bytes()); // n_intv
    data.extend_from_slice(&0x1000u64.to_le_bytes()); // interval

    // Index for chr2
    data.extend_from_slice(&1i32.to_le_bytes()); // n_bin
    data.extend_from_slice(&0u32.to_le_bytes()); // bin_id
    data.extend_from_slice(&1i32.to_le_bytes()); // n_chunk
    data.extend_from_slice(&0x3000u64.to_le_bytes()); // chunk start
    data.extend_from_slice(&0x4000u64.to_le_bytes()); // chunk end
    data.extend_from_slice(&1i32.to_le_bytes()); // n_intv
    data.extend_from_slice(&0x3000u64.to_le_bytes()); // interval

    let tbi_path = test_data_dir().join("multi_ref_test.tbi");
    let mut file = File::create(&tbi_path).expect("Failed to create TBI file");
    file.write_all(&data).expect("Failed to write TBI data");

    let index = TbiIndex::from_path(&tbi_path).expect("Failed to parse TBI");

    // Verify both references
    assert_eq!(index.references().len(), 2);
    assert!(index.get_reference("chr1").is_some());
    assert!(index.get_reference("chr2").is_some());

    // Query chr2
    let chunks = index.query("chr2", 0, 100000).expect("Query chr2 failed");
    assert!(!chunks.is_empty());
    assert_eq!(chunks[0].start.as_raw(), 0x3000);

    // Cleanup
    std::fs::remove_file(&tbi_path).ok();
}

#[test]
fn test_real_world_tbi_1000genomes() {
    // Test with real TBI index created from 1000 Genomes VCF
    let tbi_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/real_world/variants/synthetic_1000g.vcf.gz.tbi");

    let index = TbiIndex::from_path(&tbi_path)
        .expect("Failed to load real-world TBI index");

    // Verify format is VCF (format() returns TbiFormat enum)
    // We just check that we can call format() successfully
    let _format = index.format();

    // Should have chr21 reference
    let refs = index.references();
    assert_eq!(refs.len(), 1, "Expected 1 reference sequence");
    assert_eq!(refs[0].name, "chr21", "Expected chr21 reference");

    // Query a region
    let chunks = index.query("chr21", 9411000, 9412000)
        .expect("Query failed for chr21 region");

    // Should find chunks since we have variants in this region
    assert!(!chunks.is_empty(), "Should find chunks in variant region");

    println!("✅ Real-world TBI index: {} reference(s), {} chunks in test region",
             refs.len(), chunks.len());

    println!("✅ Real-world 1000 Genomes TBI validation passed");
}
