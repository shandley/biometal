//! Integration tests for VCF writer

use biometal::formats::primitives::TabDelimitedParser;
use biometal::formats::vcf::{VcfHeader, VcfRecord};
use biometal::formats::vcf_writer::VcfWriter;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use tempfile::NamedTempFile;

#[test]
fn test_vcf_writer_basic_roundtrip() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    // Write
    {
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.info_fields.insert("DP".to_string(), "Total Depth".to_string());

        let mut writer = VcfWriter::create(path).unwrap();
        writer.write_header(&header).unwrap();

        let mut info = HashMap::new();
        info.insert("DP".to_string(), "100".to_string());

        let record = VcfRecord {
            chrom: "chr1".to_string(),
            pos: 12345,
            id: Some("rs123".to_string()),
            reference: "A".to_string(),
            alternate: vec!["T".to_string()],
            quality: Some(30.0),
            filter: Some("PASS".to_string()),
            info,
            format: None,
            samples: vec![],
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Read back (parse lines manually to verify format)
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    // Verify header
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert!(lines[1].contains("##INFO=<ID=DP"));
    assert!(lines[2].starts_with("#CHROM"));

    // Verify record
    assert!(lines[3].contains("chr1"));
    assert!(lines[3].contains("12345"));
    assert!(lines[3].contains("rs123"));
    assert!(lines[3].contains("A\tT"));
    assert!(lines[3].contains("30"));
    assert!(lines[3].contains("PASS"));
    assert!(lines[3].contains("DP=100"));
}

#[test]
fn test_vcf_writer_header_with_samples() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    {
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.samples = vec!["sample1".to_string(), "sample2".to_string()];
        header.format_fields.insert("GT".to_string(), "Genotype".to_string());

        let mut writer = VcfWriter::create(path).unwrap();
        writer.write_header(&header).unwrap();
        writer.finish().unwrap();
    }

    // Verify header
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert!(lines[1].contains("##FORMAT=<ID=GT"));
    assert!(lines[2].ends_with("FORMAT\tsample1\tsample2"));
}

#[test]
fn test_vcf_writer_multiple_alternate_alleles() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    {
        let header = VcfHeader::new("VCFv4.2".to_string());
        let mut writer = VcfWriter::create(path).unwrap();
        writer.write_header(&header).unwrap();

        let record = VcfRecord {
            chrom: "chr1".to_string(),
            pos: 100,
            id: None,
            reference: "A".to_string(),
            alternate: vec!["T".to_string(), "G".to_string(), "C".to_string()],
            quality: None,
            filter: None,
            info: HashMap::new(),
            format: None,
            samples: vec![],
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Verify alternate alleles are comma-separated
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    let data_line = &lines[2]; // Skip header lines
    assert!(data_line.contains("T,G,C"));
}

#[test]
fn test_vcf_writer_with_samples() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    {
        let mut header = VcfHeader::new("VCFv4.2".to_string());
        header.samples = vec!["sample1".to_string()];

        let mut writer = VcfWriter::create(path).unwrap();
        writer.write_header(&header).unwrap();

        let record = VcfRecord {
            chrom: "chr1".to_string(),
            pos: 100,
            id: None,
            reference: "A".to_string(),
            alternate: vec!["T".to_string()],
            quality: None,
            filter: None,
            info: HashMap::new(),
            format: Some("GT".to_string()),
            samples: vec!["0/1".to_string()],
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Verify FORMAT and sample columns
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    let data_line = &lines[2];
    assert!(data_line.ends_with("GT\t0/1"));
}

#[test]
fn test_vcf_writer_info_field_formatting() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    {
        let header = VcfHeader::new("VCFv4.2".to_string());
        let mut writer = VcfWriter::create(path).unwrap();
        writer.write_header(&header).unwrap();

        let mut info = HashMap::new();
        info.insert("DP".to_string(), "50".to_string());
        info.insert("AF".to_string(), "0.5".to_string());
        info.insert("DB".to_string(), "".to_string()); // Flag

        let record = VcfRecord {
            chrom: "chr1".to_string(),
            pos: 100,
            id: None,
            reference: "A".to_string(),
            alternate: vec!["T".to_string()],
            quality: None,
            filter: None,
            info,
            format: None,
            samples: vec![],
        };

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
    }

    // Verify INFO formatting
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    let data_line = &lines[2];
    // INFO fields should be sorted and semicolon-separated
    assert!(data_line.contains("AF=0.5"));
    assert!(data_line.contains("DB;")); // Flag with no value
    assert!(data_line.contains("DP=50"));
}

#[test]
fn test_vcf_writer_validation_empty_chrom() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let header = VcfHeader::new("VCFv4.2".to_string());
    let mut writer = VcfWriter::create(path).unwrap();
    writer.write_header(&header).unwrap();

    let record = VcfRecord {
        chrom: "".to_string(), // Empty chrom
        pos: 100,
        id: None,
        reference: "A".to_string(),
        alternate: vec!["T".to_string()],
        quality: None,
        filter: None,
        info: HashMap::new(),
        format: None,
        samples: vec![],
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("chrom cannot be empty"));
}

#[test]
fn test_vcf_writer_validation_zero_pos() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let header = VcfHeader::new("VCFv4.2".to_string());
    let mut writer = VcfWriter::create(path).unwrap();
    writer.write_header(&header).unwrap();

    let record = VcfRecord {
        chrom: "chr1".to_string(),
        pos: 0, // Invalid: VCF is 1-based
        id: None,
        reference: "A".to_string(),
        alternate: vec!["T".to_string()],
        quality: None,
        filter: None,
        info: HashMap::new(),
        format: None,
        samples: vec![],
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("pos must be >= 1"));
}

#[test]
fn test_vcf_writer_validation_empty_reference() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let header = VcfHeader::new("VCFv4.2".to_string());
    let mut writer = VcfWriter::create(path).unwrap();
    writer.write_header(&header).unwrap();

    let record = VcfRecord {
        chrom: "chr1".to_string(),
        pos: 100,
        id: None,
        reference: "".to_string(), // Empty reference
        alternate: vec!["T".to_string()],
        quality: None,
        filter: None,
        info: HashMap::new(),
        format: None,
        samples: vec![],
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("reference allele cannot be empty"));
}

#[test]
fn test_vcf_writer_validation_empty_alternate() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let header = VcfHeader::new("VCFv4.2".to_string());
    let mut writer = VcfWriter::create(path).unwrap();
    writer.write_header(&header).unwrap();

    let record = VcfRecord {
        chrom: "chr1".to_string(),
        pos: 100,
        id: None,
        reference: "A".to_string(),
        alternate: vec![], // Empty alternate
        quality: None,
        filter: None,
        info: HashMap::new(),
        format: None,
        samples: vec![],
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("alternate alleles cannot be empty"));
}

#[test]
fn test_vcf_writer_validation_header_not_written() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let mut writer = VcfWriter::create(path).unwrap();

    let record = VcfRecord {
        chrom: "chr1".to_string(),
        pos: 100,
        id: None,
        reference: "A".to_string(),
        alternate: vec!["T".to_string()],
        quality: None,
        filter: None,
        info: HashMap::new(),
        format: None,
        samples: vec![],
    };

    let result = writer.write_record(&record);
    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("Header must be written before records"));
}

#[test]
fn test_vcf_writer_multiple_records() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();

    {
        let header = VcfHeader::new("VCFv4.2".to_string());
        let mut writer = VcfWriter::create(path).unwrap();
        writer.write_header(&header).unwrap();

        for i in 0..10 {
            let record = VcfRecord {
                chrom: "chr1".to_string(),
                pos: (i + 1) * 100,
                id: Some(format!("var{}", i)),
                reference: "A".to_string(),
                alternate: vec!["T".to_string()],
                quality: Some(i as f64 * 10.0),
                filter: Some("PASS".to_string()),
                info: HashMap::new(),
                format: None,
                samples: vec![],
            };
            writer.write_record(&record).unwrap();
        }

        assert_eq!(writer.records_written(), 10);
        writer.finish().unwrap();
    }

    // Verify all records written
    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>().unwrap();

    // Header (2 lines) + 10 records
    assert_eq!(lines.len(), 12);

    // Verify positions
    for i in 0..10 {
        let pos = (i + 1) * 100;
        let line = &lines[i + 2]; // Skip header
        assert!(line.contains(&pos.to_string()));
        assert!(line.contains(&format!("var{}", i)));
    }
}

#[test]
fn test_vcf_writer_records_written_counter() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let header = VcfHeader::new("VCFv4.2".to_string());
    let mut writer = VcfWriter::create(path).unwrap();
    writer.write_header(&header).unwrap();

    assert_eq!(writer.records_written(), 0);

    for i in 0..50 {
        let record = VcfRecord {
            chrom: "chr1".to_string(),
            pos: (i + 1) * 100,
            id: None,
            reference: "A".to_string(),
            alternate: vec!["T".to_string()],
            quality: None,
            filter: None,
            info: HashMap::new(),
            format: None,
            samples: vec![],
        };
        writer.write_record(&record).unwrap();
    }

    assert_eq!(writer.records_written(), 50);
    writer.finish().unwrap();
}

#[test]
fn test_vcf_writer_header_written_check() {
    let temp_file = NamedTempFile::new().unwrap();
    let path = temp_file.path();
    let header = VcfHeader::new("VCFv4.2".to_string());
    let mut writer = VcfWriter::create(path).unwrap();

    assert!(!writer.header_written());
    writer.write_header(&header).unwrap();
    assert!(writer.header_written());
}
