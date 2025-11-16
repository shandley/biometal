//! Integration tests for BAM writer with real-world BAM files.
//!
//! These tests verify that BamWriter correctly handles real BAM data:
//! - Round-trip accuracy (read → write → read, verify identical)
//! - Filtering workflows (quality filtering, region subsetting)
//! - Large file handling (streaming 100K+ records)
//! - Complex CIGAR operations and tags

use biometal::io::bam::{BamReader, BamWriter};
use biometal::Result;
use std::io;
use tempfile::NamedTempFile;

/// Test round-trip: read real BAM, write to new file, verify identical
#[test]
fn test_round_trip_small_bam() -> Result<()> {
    let input_path = "tests/data/synthetic/alignment/synthetic_100k.bam";
    let temp_file = NamedTempFile::new()?;
    let output_path = temp_file.path();

    // Read input BAM and write to output
    {
        let mut reader = BamReader::from_path(input_path)?;
        let header = reader.header().clone();
        let mut writer = BamWriter::create(output_path, header)?;

        for record in reader.records() {
            let record = record?;
            writer.write_record(&record)?;
        }

        writer.finish()?;
    }

    // Read both files and compare
    let mut original_reader = BamReader::from_path(input_path)?;
    let mut roundtrip_reader = BamReader::from_path(output_path)?;

    // Compare headers
    assert_eq!(
        original_reader.header().references.len(),
        roundtrip_reader.header().references.len(),
        "Header reference count mismatch"
    );

    // Compare records
    let original_records: Vec<_> = original_reader
        .records()
        .collect::<io::Result<Vec<_>>>()?;
    let roundtrip_records: Vec<_> = roundtrip_reader
        .records()
        .collect::<io::Result<Vec<_>>>()?;

    assert_eq!(
        original_records.len(),
        roundtrip_records.len(),
        "Record count mismatch"
    );

    for (i, (orig, rt)) in original_records
        .iter()
        .zip(roundtrip_records.iter())
        .enumerate()
    {
        assert_eq!(orig.name, rt.name, "Record {}: name mismatch", i);
        assert_eq!(
            orig.reference_id, rt.reference_id,
            "Record {}: reference_id mismatch",
            i
        );
        assert_eq!(orig.position, rt.position, "Record {}: position mismatch", i);
        assert_eq!(orig.mapq, rt.mapq, "Record {}: MAPQ mismatch", i);
        assert_eq!(orig.flags, rt.flags, "Record {}: flags mismatch", i);
        assert_eq!(orig.sequence, rt.sequence, "Record {}: sequence mismatch", i);
        assert_eq!(orig.quality, rt.quality, "Record {}: quality mismatch", i);
        assert_eq!(orig.cigar, rt.cigar, "Record {}: CIGAR mismatch", i);
    }

    Ok(())
}

/// Test filtering workflow: read BAM, filter by MAPQ, write filtered BAM
#[test]
fn test_quality_filtering_workflow() -> Result<()> {
    let input_path = "tests/data/synthetic/alignment/synthetic_100k.bam";
    let temp_file = NamedTempFile::new()?;
    let filtered_path = temp_file.path();

    const MIN_MAPQ: u8 = 30;

    // Filter and write
    let filtered_count = {
        let mut reader = BamReader::from_path(input_path)?;
        let header = reader.header().clone();
        let mut writer = BamWriter::create(filtered_path, header)?;

        for record in reader.records() {
            let record = record?;
            if record.mapq.unwrap_or(0) >= MIN_MAPQ {
                writer.write_record(&record)?;
            }
        }

        let count = writer.records_written();
        writer.finish()?;
        count
    };

    // Verify all written records meet quality threshold
    let mut reader = BamReader::from_path(filtered_path)?;
    let mut verified_count = 0;

    for record in reader.records() {
        let record = record?;
        assert!(
            record.mapq.unwrap_or(0) >= MIN_MAPQ,
            "Found record with MAPQ < {} in filtered file",
            MIN_MAPQ
        );
        verified_count += 1;
    }

    assert_eq!(
        filtered_count, verified_count,
        "Written count doesn't match verified count"
    );

    // Verify we actually filtered some records
    let mut original_reader = BamReader::from_path(input_path)?;
    let original_count = original_reader.records().count();

    assert!(
        filtered_count < original_count,
        "Filtering should remove some records (filtered: {}, original: {})",
        filtered_count,
        original_count
    );

    Ok(())
}

/// Test unmapped reads filtering
#[test]
fn test_filter_mapped_reads() -> Result<()> {
    let input_path = "tests/data/synthetic/alignment/synthetic_100k.bam";
    let temp_file = NamedTempFile::new()?;
    let mapped_path = temp_file.path();

    // Filter for mapped reads only
    let mapped_count = {
        let mut reader = BamReader::from_path(input_path)?;
        let header = reader.header().clone();
        let mut writer = BamWriter::create(mapped_path, header)?;

        for record in reader.records() {
            let record = record?;
            if !record.is_unmapped() {
                writer.write_record(&record)?;
            }
        }

        let count = writer.records_written();
        writer.finish()?;
        count
    };

    // Verify all records are mapped
    let mut reader = BamReader::from_path(mapped_path)?;

    for record in reader.records() {
        let record = record?;
        assert!(
            !record.is_unmapped(),
            "Found unmapped record in mapped-only file: {}",
            record.name
        );
        assert!(record.position.is_some(), "Mapped record missing position");
    }

    println!("Filtered to {} mapped reads", mapped_count);

    Ok(())
}

/// Test region subsetting (write only records in specific position range)
#[test]
fn test_region_subsetting() -> Result<()> {
    let input_path = "tests/data/synthetic/alignment/synthetic_100k.bam";
    let temp_file = NamedTempFile::new()?;
    let subset_path = temp_file.path();

    const MIN_POS: i32 = 10000;
    const MAX_POS: i32 = 20000;

    // Write only records in position range
    let subset_count = {
        let mut reader = BamReader::from_path(input_path)?;
        let header = reader.header().clone();
        let mut writer = BamWriter::create(subset_path, header)?;

        for record in reader.records() {
            let record = record?;
            if let Some(pos) = record.position {
                if pos >= MIN_POS && pos <= MAX_POS {
                    writer.write_record(&record)?;
                }
            }
        }

        let count = writer.records_written();
        writer.finish()?;
        count
    };

    // Verify all records are in range
    let mut reader = BamReader::from_path(subset_path)?;

    for record in reader.records() {
        let record = record?;
        if let Some(pos) = record.position {
            assert!(
                pos >= MIN_POS && pos <= MAX_POS,
                "Record at position {} outside range [{}, {}]",
                pos,
                MIN_POS,
                MAX_POS
            );
        }
    }

    println!("Subset contains {} records in range", subset_count);

    Ok(())
}

/// Test streaming large BAM file (100K records)
#[test]
fn test_large_file_streaming() -> Result<()> {
    let input_path = "tests/data/synthetic/alignment/synthetic_100k.bam";
    let temp_file = NamedTempFile::new()?;
    let output_path = temp_file.path();

    // Stream large file: read and write
    let record_count = {
        let mut reader = BamReader::from_path(input_path)?;
        let header = reader.header().clone();
        let mut writer = BamWriter::create(output_path, header)?;

        for record in reader.records() {
            let record = record?;
            writer.write_record(&record)?;
        }

        let count = writer.records_written();
        writer.finish()?;
        count
    };

    // Verify count
    assert!(
        record_count >= 100_000,
        "Expected at least 100K records, got {}",
        record_count
    );

    // Verify output file is readable and has same count
    let mut reader = BamReader::from_path(output_path)?;
    let output_count = reader.records().count();

    assert_eq!(
        record_count, output_count,
        "Output record count mismatch"
    );

    println!("Successfully streamed {} records", record_count);

    Ok(())
}

/// Test preserving complex CIGAR operations
#[test]
fn test_complex_cigar_preservation() -> Result<()> {
    let input_path = "tests/data/synthetic/alignment/synthetic_100k.bam";
    let temp_file = NamedTempFile::new()?;
    let output_path = temp_file.path();

    // Find record with complex CIGAR (multiple operations)
    let mut complex_cigar_found = false;

    {
        let mut reader = BamReader::from_path(input_path)?;
        let header = reader.header().clone();
        let mut writer = BamWriter::create(output_path, header)?;

        for record in reader.records() {
            let record = record?;

            // Check for complex CIGAR (> 1 operation)
            if record.cigar.len() > 1 {
                complex_cigar_found = true;
            }

            writer.write_record(&record)?;
        }

        writer.finish()?;
    }

    if complex_cigar_found {
        // Verify CIGAR preserved correctly
        let mut original_reader = BamReader::from_path(input_path)?;
        let mut roundtrip_reader = BamReader::from_path(output_path)?;

        for (orig, rt) in original_reader
            .records()
            .zip(roundtrip_reader.records())
        {
            let orig = orig?;
            let rt = rt?;

            assert_eq!(
                orig.cigar.len(),
                rt.cigar.len(),
                "CIGAR operation count mismatch for read {}",
                orig.name
            );

            for (i, (orig_op, rt_op)) in orig.cigar.iter().zip(rt.cigar.iter()).enumerate() {
                assert_eq!(
                    orig_op, rt_op,
                    "CIGAR operation {} mismatch for read {}: {:?} vs {:?}",
                    i, orig.name, orig_op, rt_op
                );
            }
        }

        println!("Complex CIGAR operations preserved correctly");
    } else {
        println!("No complex CIGAR operations found in test file");
    }

    Ok(())
}

/// Test header preservation
#[test]
fn test_header_preservation() -> Result<()> {
    let input_path = "tests/data/synthetic/alignment/synthetic_100k.bam";
    let temp_file = NamedTempFile::new()?;
    let output_path = temp_file.path();

    // Copy BAM with header
    {
        let mut reader = BamReader::from_path(input_path)?;
        let header = reader.header().clone();
        let mut writer = BamWriter::create(output_path, header)?;

        for record in reader.records() {
            writer.write_record(&record?)?;
        }

        writer.finish()?;
    }

    // Compare headers
    let original_reader = BamReader::from_path(input_path)?;
    let roundtrip_reader = BamReader::from_path(output_path)?;

    let orig_header = original_reader.header();
    let rt_header = roundtrip_reader.header();

    // Compare reference count
    assert_eq!(
        orig_header.references.len(),
        rt_header.references.len(),
        "Reference count mismatch"
    );

    // Compare each reference
    for (orig_ref, rt_ref) in orig_header.references.iter().zip(rt_header.references.iter()) {
        assert_eq!(orig_ref.name, rt_ref.name, "Reference name mismatch");
        assert_eq!(orig_ref.length, rt_ref.length, "Reference length mismatch");
    }

    println!(
        "Header with {} references preserved correctly",
        orig_header.references.len()
    );

    Ok(())
}

/// Test empty BAM (header only, no records)
#[test]
fn test_empty_bam() -> Result<()> {
    use biometal::io::bam::{Header, Reference};

    let temp_file = NamedTempFile::new()?;
    let empty_path = temp_file.path();

    // Create empty BAM (header only)
    {
        let header = Header::new(
            "@HD\tVN:1.6\n".to_string(),
            vec![Reference::new("chr1".to_string(), 1000000)],
        );
        let writer = BamWriter::create(empty_path, header)?;
        writer.finish()?;
    }

    // Verify it can be read
    let mut reader = BamReader::from_path(empty_path)?;

    assert_eq!(reader.header().references.len(), 1);
    assert_eq!(reader.header().references[0].name, "chr1");

    let record_count = reader.records().count();
    assert_eq!(record_count, 0, "Empty BAM should have no records");

    Ok(())
}
