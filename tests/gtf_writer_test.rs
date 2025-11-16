//! Integration tests for GTF writer

use biometal::formats::gtf::GtfRecord;
use biometal::formats::gtf_writer::GtfWriter;
use biometal::formats::primitives::{Strand, TabDelimitedParser};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use tempfile::NamedTempFile;

#[test]
fn test_gtf_writer_basic_roundtrip() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write
    {
        let mut writer = GtfWriter::create(path).unwrap();

        let mut attributes = HashMap::new();
        attributes.insert("gene_id".to_string(), "ENSG00000223972".to_string());
        attributes.insert("transcript_id".to_string(), "ENST00000456328".to_string());
        attributes.insert("gene_name".to_string(), "DDX11L1".to_string());

        let record = GtfRecord {
            seqname: "chr1".to_string(),
            source: "HAVANA".to_string(),
            feature: "exon".to_string(),
            start: 11869,
            end: 12227,
            score: Some(100.0),
            strand: Strand::Forward,
            frame: None,
            attributes,
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, GtfRecord>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].seqname, "chr1");
    assert_eq!(records[0].source, "HAVANA");
    assert_eq!(records[0].feature, "exon");
    assert_eq!(records[0].start, 11869);
    assert_eq!(records[0].end, 12227);
    assert!((records[0].score.unwrap() - 100.0).abs() < 0.01);
    assert_eq!(records[0].strand, Strand::Forward);
    assert_eq!(records[0].frame, None);
    assert_eq!(records[0].gene_id(), "ENSG00000223972");
    assert_eq!(records[0].transcript_id(), Some("ENST00000456328"));
    assert_eq!(records[0].gene_name(), Some("DDX11L1"));
}

#[test]
fn test_gtf_writer_attribute_format() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write
    {
        let mut writer = GtfWriter::create(path).unwrap();

        let mut attributes = HashMap::new();
        attributes.insert("gene_id".to_string(), "ENSG001".to_string());
        attributes.insert("transcript_id".to_string(), "ENST001".to_string());

        let record = GtfRecord {
            seqname: "chr1".to_string(),
            source: "test".to_string(),
            feature: "exon".to_string(),
            start: 100,
            end: 200,
            score: None,
            strand: Strand::Forward,
            frame: None,
            attributes,
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Verify GTF attribute format: key "value";
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    // Should have GTF attribute syntax with quotes and semicolons
    assert!(lines[0].contains("gene_id \"ENSG001\""));
    assert!(lines[0].contains("transcript_id \"ENST001\""));
    assert!(lines[0].ends_with(";"));
}

#[test]
fn test_gtf_writer_gene_feature() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write gene feature (transcript_id is optional for gene features)
    {
        let mut writer = GtfWriter::create(path).unwrap();

        let mut attributes = HashMap::new();
        attributes.insert("gene_id".to_string(), "ENSG001".to_string());
        attributes.insert("gene_name".to_string(), "ABC1".to_string());

        let record = GtfRecord {
            seqname: "chr1".to_string(),
            source: "GENCODE".to_string(),
            feature: "gene".to_string(),
            start: 1000,
            end: 5000,
            score: None,
            strand: Strand::Forward,
            frame: None,
            attributes,
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, GtfRecord>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].feature, "gene");
    assert_eq!(records[0].gene_id(), "ENSG001");
    assert_eq!(records[0].transcript_id(), None);  // Optional for gene features
    assert_eq!(records[0].gene_name(), Some("ABC1"));
}

#[test]
fn test_gtf_writer_optional_fields() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write with all optional fields missing
    {
        let mut writer = GtfWriter::create(path).unwrap();

        let mut attributes = HashMap::new();
        attributes.insert("gene_id".to_string(), "ENSG001".to_string());

        let record = GtfRecord {
            seqname: "chr1".to_string(),
            source: "test".to_string(),
            feature: "gene".to_string(),
            start: 100,
            end: 200,
            score: None,
            strand: Strand::Unknown,
            frame: None,
            attributes,
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, GtfRecord>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].score, None);
    assert_eq!(records[0].strand, Strand::Unknown);
    assert_eq!(records[0].frame, None);
}

#[test]
fn test_gtf_writer_with_frame() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write CDS with frame
    {
        let mut writer = GtfWriter::create(path).unwrap();

        for frame in 0..=2 {
            let mut attributes = HashMap::new();
            attributes.insert("gene_id".to_string(), "ENSG001".to_string());
            attributes.insert("transcript_id".to_string(), "ENST001".to_string());

            let record = GtfRecord {
                seqname: "chr1".to_string(),
                source: "test".to_string(),
                feature: "CDS".to_string(),
                start: 100 + frame as u64 * 100,
                end: 200 + frame as u64 * 100,
                score: None,
                strand: Strand::Forward,
                frame: Some(frame),
                attributes,
            };
            writer.write_record(&record).unwrap();
        }

        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, GtfRecord>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 3);
    assert_eq!(records[0].frame, Some(0));
    assert_eq!(records[1].frame, Some(1));
    assert_eq!(records[2].frame, Some(2));
}

#[test]
fn test_gtf_writer_validation_empty_seqname() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = GtfWriter::create(path).unwrap();

    let mut attributes = HashMap::new();
    attributes.insert("gene_id".to_string(), "ENSG001".to_string());

    let record = GtfRecord {
        seqname: "".to_string(),  // Empty seqname
        source: "test".to_string(),
        feature: "gene".to_string(),
        start: 100,
        end: 200,
        score: None,
        strand: Strand::Forward,
        frame: None,
        attributes,
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("seqname cannot be empty"));
}

#[test]
fn test_gtf_writer_validation_invalid_interval() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = GtfWriter::create(path).unwrap();

    let mut attributes = HashMap::new();
    attributes.insert("gene_id".to_string(), "ENSG001".to_string());

    let record = GtfRecord {
        seqname: "chr1".to_string(),
        source: "test".to_string(),
        feature: "gene".to_string(),
        start: 200,  // start >= end
        end: 100,
        score: None,
        strand: Strand::Forward,
        frame: None,
        attributes,
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid interval"));
}

#[test]
fn test_gtf_writer_validation_invalid_frame() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = GtfWriter::create(path).unwrap();

    let mut attributes = HashMap::new();
    attributes.insert("gene_id".to_string(), "ENSG001".to_string());
    attributes.insert("transcript_id".to_string(), "ENST001".to_string());

    let record = GtfRecord {
        seqname: "chr1".to_string(),
        source: "test".to_string(),
        feature: "CDS".to_string(),
        start: 100,
        end: 200,
        score: None,
        strand: Strand::Forward,
        frame: Some(3),  // Invalid frame (must be 0, 1, or 2)
        attributes,
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Invalid frame"));
}

#[test]
fn test_gtf_writer_validation_missing_gene_id() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = GtfWriter::create(path).unwrap();

    let attributes = HashMap::new();  // No gene_id

    let record = GtfRecord {
        seqname: "chr1".to_string(),
        source: "test".to_string(),
        feature: "gene".to_string(),
        start: 100,
        end: 200,
        score: None,
        strand: Strand::Forward,
        frame: None,
        attributes,
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("gene_id attribute is required"));
}

#[test]
fn test_gtf_writer_validation_missing_transcript_id() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = GtfWriter::create(path).unwrap();

    let mut attributes = HashMap::new();
    attributes.insert("gene_id".to_string(), "ENSG001".to_string());
    // Missing transcript_id for exon feature

    let record = GtfRecord {
        seqname: "chr1".to_string(),
        source: "test".to_string(),
        feature: "exon".to_string(),  // Requires transcript_id
        start: 100,
        end: 200,
        score: None,
        strand: Strand::Forward,
        frame: None,
        attributes,
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("transcript_id attribute is required"));
}

#[test]
fn test_gtf_writer_multiple_records() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write multiple records
    {
        let mut writer = GtfWriter::create(path).unwrap();

        for i in 0..10 {
            let mut attributes = HashMap::new();
            attributes.insert("gene_id".to_string(), format!("ENSG{:05}", i));
            attributes.insert("gene_name".to_string(), format!("GENE{}", i));

            let record = GtfRecord {
                seqname: "chr1".to_string(),
                source: "test".to_string(),
                feature: "gene".to_string(),
                start: i * 1000,
                end: (i + 1) * 1000,
                score: Some(i as f64 * 10.0),
                strand: Strand::Forward,
                frame: None,
                attributes,
            };

            writer.write_record(&record).unwrap();
        }

        assert_eq!(writer.records_written(), 10);
        writer.finish().unwrap();
    }

    // Read back
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let parser = TabDelimitedParser::<_, GtfRecord>::new(reader);
    let records: Vec<_> = parser.collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 10);
    for i in 0..10 {
        assert_eq!(records[i].start, i as u64 * 1000);
        assert_eq!(records[i].end, (i as u64 + 1) * 1000);
        assert_eq!(records[i].gene_id(), format!("ENSG{:05}", i));
        assert_eq!(records[i].gene_name(), Some(format!("GENE{}", i).as_str()));
    }
}

#[test]
fn test_gtf_writer_records_written_counter() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = GtfWriter::create(path).unwrap();

    assert_eq!(writer.records_written(), 0);

    for i in 0..50 {
        let mut attributes = HashMap::new();
        attributes.insert("gene_id".to_string(), format!("ENSG{:05}", i));

        let record = GtfRecord {
            seqname: "chr1".to_string(),
            source: "test".to_string(),
            feature: "gene".to_string(),
            start: i * 1000,
            end: (i + 1) * 1000,
            score: None,
            strand: Strand::Forward,
            frame: None,
            attributes,
        };
        writer.write_record(&record).unwrap();
    }

    assert_eq!(writer.records_written(), 50);
    writer.finish().unwrap();
}
