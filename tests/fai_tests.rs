//! Integration tests for FAI (FASTA Index) functionality

use biometal::io::fasta::{FaiIndex, FastaStream};
use std::fs;
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/synthetic/sequence")
}

#[test]
fn test_fai_build_and_write() {
    let fasta_path = test_data_dir().join("test.fa");
    let fai_path = test_data_dir().join("test.fa.fai");

    // Build index from FASTA file
    let index = FaiIndex::build(&fasta_path).expect("Failed to build FAI index");

    // Verify we found all sequences
    assert_eq!(index.len(), 3);
    assert!(index.get("chr1").is_some());
    assert!(index.get("chr2").is_some());
    assert!(index.get("chr3").is_some());

    // Check chr1 details
    let chr1 = index.get("chr1").unwrap();
    assert_eq!(chr1.name, "chr1");
    assert_eq!(chr1.length, 44); // 28 + 16 bases
    assert_eq!(chr1.line_bases, 28);
    assert_eq!(chr1.line_width, 29); // 28 + newline

    // Check chr2 details
    let chr2 = index.get("chr2").unwrap();
    assert_eq!(chr2.name, "chr2");
    assert_eq!(chr2.length, 16);
    assert_eq!(chr2.line_bases, 16);
    assert_eq!(chr2.line_width, 17);

    // Check chr3 details
    let chr3 = index.get("chr3").unwrap();
    assert_eq!(chr3.name, "chr3");
    assert_eq!(chr3.length, 64); // 28 + 28 + 8 bases
    assert_eq!(chr3.line_bases, 28);
    assert_eq!(chr3.line_width, 29);

    // Write index to file
    index.write(&fai_path).expect("Failed to write FAI index");

    // Read it back
    let loaded_index = FaiIndex::from_path(&fai_path).expect("Failed to load FAI index");

    // Verify loaded index matches
    assert_eq!(loaded_index.len(), 3);

    let loaded_chr1 = loaded_index.get("chr1").unwrap();
    assert_eq!(loaded_chr1.name, chr1.name);
    assert_eq!(loaded_chr1.length, chr1.length);
    assert_eq!(loaded_chr1.offset, chr1.offset);
    assert_eq!(loaded_chr1.line_bases, chr1.line_bases);
    assert_eq!(loaded_chr1.line_width, chr1.line_width);

    // Cleanup
    fs::remove_file(&fai_path).ok();
}

#[test]
fn test_fai_fetch_entire_sequence() {
    let fasta_path = test_data_dir().join("test.fa");

    // Build index
    let index = FaiIndex::build(&fasta_path).expect("Failed to build index");

    // Fetch chr1 (entire sequence)
    let chr1_seq = index.fetch("chr1", &fasta_path).expect("Failed to fetch chr1");
    assert_eq!(chr1_seq.len(), 44);
    assert!(chr1_seq.starts_with("ACGTACGTACGTACGTACGTACGTACGT"));
    assert!(chr1_seq.ends_with("TGCATGCATGCATGCA"));

    // Fetch chr2
    let chr2_seq = index.fetch("chr2", &fasta_path).expect("Failed to fetch chr2");
    assert_eq!(chr2_seq, "GGGGCCCCAAAATTTT");

    // Fetch chr3
    let chr3_seq = index.fetch("chr3", &fasta_path).expect("Failed to fetch chr3");
    assert_eq!(chr3_seq.len(), 64);
    assert!(chr3_seq.starts_with("ATCGATCGATCGATCGATCGATCGATCG"));
    assert!(chr3_seq.ends_with("ATCGATCG"));
}

#[test]
fn test_fai_fetch_region() {
    let fasta_path = test_data_dir().join("test.fa");
    let index = FaiIndex::build(&fasta_path).expect("Failed to build index");

    // Fetch first 10 bases of chr1
    let region = index
        .fetch_region("chr1", 0, 10, &fasta_path)
        .expect("Failed to fetch region");
    assert_eq!(region, "ACGTACGTAC");

    // Fetch bases 28-38 (crosses line boundary)
    let region = index
        .fetch_region("chr1", 28, 38, &fasta_path)
        .expect("Failed to fetch region");
    assert_eq!(region, "TGCATGCATG");

    // Fetch last 5 bases of chr2 (positions 11-15, half-open [11, 16))
    let region = index
        .fetch_region("chr2", 11, 16, &fasta_path)
        .expect("Failed to fetch region");
    assert_eq!(region, "ATTTT");

    // Fetch region that exceeds length (should truncate)
    let region = index
        .fetch_region("chr2", 10, 100, &fasta_path)
        .expect("Failed to fetch region");
    assert_eq!(region, "AATTTT");
}

#[test]
fn test_fai_fetch_nonexistent_sequence() {
    let fasta_path = test_data_dir().join("test.fa");
    let index = FaiIndex::build(&fasta_path).expect("Failed to build index");

    // Try to fetch non-existent sequence
    let result = index.fetch("chr99", &fasta_path);
    assert!(result.is_err());
    assert!(result
        .unwrap_err()
        .to_string()
        .contains("not found in index"));
}

#[test]
fn test_fai_fetch_invalid_range() {
    let fasta_path = test_data_dir().join("test.fa");
    let index = FaiIndex::build(&fasta_path).expect("Failed to build index");

    // Invalid range: start >= end
    let result = index.fetch_region("chr1", 10, 5, &fasta_path);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid range"));

    // Start position beyond sequence length
    let result = index.fetch_region("chr1", 1000, 2000, &fasta_path);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("exceeds"));
}

#[test]
fn test_fai_sequence_order_preserved() {
    let fasta_path = test_data_dir().join("test.fa");
    let index = FaiIndex::build(&fasta_path).expect("Failed to build index");

    // Check that sequence order matches FASTA file
    assert_eq!(index.sequence_names[0], "chr1");
    assert_eq!(index.sequence_names[1], "chr2");
    assert_eq!(index.sequence_names[2], "chr3");
}

#[test]
fn test_fai_validates_against_streaming_parser() {
    let fasta_path = test_data_dir().join("test.fa");

    // Get sequences using streaming parser
    let stream = FastaStream::from_path(&fasta_path).expect("Failed to create stream");
    let mut streaming_seqs: Vec<(String, String)> = Vec::new();
    for record in stream {
        let record = record.expect("Failed to read record");
        let seq_string = String::from_utf8(record.sequence).expect("Invalid UTF-8");
        streaming_seqs.push((record.id, seq_string));
    }

    // Get sequences using FAI
    let index = FaiIndex::build(&fasta_path).expect("Failed to build index");

    // Compare results
    assert_eq!(streaming_seqs.len(), index.len());

    for (name, seq) in &streaming_seqs {
        let indexed_seq = index
            .fetch(name, &fasta_path)
            .expect(&format!("Failed to fetch {}", name));
        assert_eq!(
            seq, &indexed_seq,
            "Sequence mismatch for {} (streaming vs indexed)",
            name
        );
    }
}

#[test]
fn test_fai_round_trip() {
    let fasta_path = test_data_dir().join("test.fa");
    let fai_path = test_data_dir().join("roundtrip.fai");

    // Build, write, and reload
    let index1 = FaiIndex::build(&fasta_path).expect("Failed to build index");
    index1.write(&fai_path).expect("Failed to write index");
    let index2 = FaiIndex::from_path(&fai_path).expect("Failed to reload index");

    // Verify indexes are identical
    assert_eq!(index1.len(), index2.len());
    assert_eq!(index1.sequence_names, index2.sequence_names);

    for name in &index1.sequence_names {
        let entry1 = index1.get(name).unwrap();
        let entry2 = index2.get(name).unwrap();

        assert_eq!(entry1.name, entry2.name);
        assert_eq!(entry1.length, entry2.length);
        assert_eq!(entry1.offset, entry2.offset);
        assert_eq!(entry1.line_bases, entry2.line_bases);
        assert_eq!(entry1.line_width, entry2.line_width);
    }

    // Cleanup
    fs::remove_file(&fai_path).ok();
}
