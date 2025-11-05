//! Integration tests for writing operations
//!
//! These tests demonstrate complete read → transform → write pipelines,
//! which are the foundation for CLI tool development.

use biometal::{FastqStream, FastqWriter, FastqRecord};
use biometal::io::DataSink;
use biometal::operations::{gc_content, mean_quality};
use std::io::Write;
use tempfile::TempDir;

/// Test complete pipeline: Read → Filter by GC content → Write
#[test]
fn test_pipeline_gc_filter() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = temp_dir.path().join("input.fq");
    let output_path = temp_dir.path().join("output.fq");

    // Create test data with varying GC content
    {
        let mut file = std::fs::File::create(&input_path).unwrap();
        // GC=0.00 (all AT) - should be filtered out
        writeln!(file, "@read1\nAAAATTTT\n+\nIIIIIIII").unwrap();
        // GC=0.50 (4G+4C out of 8) - should pass
        writeln!(file, "@read2\nATGCATGC\n+\nIIIIIIII").unwrap();
        // GC=1.00 (all GC) - should be filtered out
        writeln!(file, "@read3\nGCGCGCGC\n+\nIIIIIIII").unwrap();
        // GC=0.50 (2G+2C out of 8) - should pass
        writeln!(file, "@read4\nATGCGCAT\n+\nIIIIIIII").unwrap();
    }

    // Pipeline: Read → Filter → Write
    {
        let input = FastqStream::from_path(&input_path).unwrap();
        let sink = DataSink::from_path(&output_path);
        let mut writer = FastqWriter::new(sink).unwrap();

        let mut total = 0;
        let mut passed = 0;

        for record in input {
            let record = record.unwrap();
            total += 1;

            // Filter: Keep records with GC content between 40-60%
            let gc = gc_content(&record.sequence);
            if gc >= 0.4 && gc <= 0.6 {
                writer.write_record(&record).unwrap();
                passed += 1;
            }
        }

        writer.finish().unwrap();

        assert_eq!(total, 4, "Should process 4 input records");
        assert_eq!(passed, 2, "Should pass 2 records through filter");
    }

    // Verify output
    {
        let output = FastqStream::from_path(&output_path).unwrap();
        let records: Vec<_> = output.collect::<biometal::Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "read2");
        assert_eq!(records[1].id, "read4");
    }
}

/// Test complete pipeline: Read → Filter by quality → Write
#[test]
fn test_pipeline_quality_filter() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = temp_dir.path().join("input.fq");
    let output_path = temp_dir.path().join("output.fq");

    // Create test data with varying quality scores
    {
        let mut file = std::fs::File::create(&input_path).unwrap();
        // Q=40 (high quality) - should pass
        writeln!(file, "@read1\nATGC\n+\nIIII").unwrap();
        // Q=0 (low quality) - should be filtered
        writeln!(file, "@read2\nGCTA\n+\n!!!!").unwrap();
        // Q=30 (good quality) - should pass
        writeln!(file, "@read3\nTGCA\n+\n????").unwrap();
        // Q=10 (poor quality) - should be filtered
        writeln!(file, "@read4\nACGT\n+\n++++").unwrap();
    }

    // Pipeline: Read → Filter → Write
    {
        let input = FastqStream::from_path(&input_path).unwrap();
        let sink = DataSink::from_path(&output_path);
        let mut writer = FastqWriter::new(sink).unwrap();

        for record in input {
            let record = record.unwrap();

            // Filter: Keep records with mean quality ≥ 25
            let mean_q = mean_quality(&record.quality);
            if mean_q >= 25.0 {
                writer.write_record(&record).unwrap();
            }
        }

        writer.finish().unwrap();
    }

    // Verify output
    {
        let output = FastqStream::from_path(&output_path).unwrap();
        let records: Vec<_> = output.collect::<biometal::Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "read1");
        assert_eq!(records[1].id, "read3");
    }
}

/// Test complete pipeline: Read → Transform sequences → Write
#[test]
fn test_pipeline_transform() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = temp_dir.path().join("input.fq");
    let output_path = temp_dir.path().join("output.fq");

    // Create test data
    {
        let mut file = std::fs::File::create(&input_path).unwrap();
        writeln!(file, "@read1\nATGC\n+\nIIII").unwrap();
        writeln!(file, "@read2\nGCTA\n+\nHHHH").unwrap();
    }

    // Pipeline: Read → Transform → Write
    {
        let input = FastqStream::from_path(&input_path).unwrap();
        let sink = DataSink::from_path(&output_path);
        let mut writer = FastqWriter::new(sink).unwrap();

        for record in input {
            let mut record = record.unwrap();

            // Transform: Convert to uppercase (identity operation for test)
            record.sequence = record.sequence.iter().map(|&b| b).collect();

            // Add suffix to ID
            record.id = format!("{}_processed", record.id);

            writer.write_record(&record).unwrap();
        }

        writer.finish().unwrap();
    }

    // Verify transformations
    {
        let output = FastqStream::from_path(&output_path).unwrap();
        let records: Vec<_> = output.collect::<biometal::Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "read1_processed");
        assert_eq!(records[1].id, "read2_processed");
    }
}

/// Test large-scale pipeline: 10K records
#[test]
fn test_pipeline_large_scale() {
    let temp_dir = TempDir::new().unwrap();
    let input_path = temp_dir.path().join("input.fq");
    let output_path = temp_dir.path().join("output.fq");

    // Create 10K test records
    {
        let mut file = std::fs::File::create(&input_path).unwrap();
        for i in 0..10_000 {
            writeln!(
                file,
                "@read_{}\nATGCATGCATGC\n+\nIIIIIIIIIIII",
                i
            ).unwrap();
        }
    }

    // Pipeline: Read → Filter (every 10th record) → Write
    {
        let input = FastqStream::from_path(&input_path).unwrap();
        let sink = DataSink::from_path(&output_path);
        let mut writer = FastqWriter::new(sink).unwrap();

        for (i, record) in input.enumerate() {
            let record = record.unwrap();

            // Keep every 10th record
            if i % 10 == 0 {
                writer.write_record(&record).unwrap();
            }
        }

        writer.finish().unwrap();
    }

    // Verify count
    {
        let output = FastqStream::from_path(&output_path).unwrap();
        let count = output.count();
        assert_eq!(count, 1_000, "Should have 1000 records (every 10th of 10K)");
    }
}

/// Test that writer properly validates records
#[test]
fn test_pipeline_with_validation_errors() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.fq");

    let sink = DataSink::from_path(&output_path);
    let mut writer = FastqWriter::new(sink).unwrap();

    // Valid record - should succeed
    let good_record = FastqRecord::new(
        "read1".to_string(),
        b"ATGC".to_vec(),
        b"IIII".to_vec(),
    );
    assert!(writer.write_record(&good_record).is_ok());

    // Invalid record - should fail
    let bad_record = FastqRecord::new(
        "read2".to_string(),
        b"ATGC".to_vec(),
        b"III".to_vec(), // Too short!
    );
    assert!(writer.write_record(&bad_record).is_err());
}
