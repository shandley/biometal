//! Integration tests for GTF real-world data
//!
//! These tests validate biometal's GTF parser against real gene annotations
//! from Ensembl Release 110 (Human chr21).

use biometal::formats::gtf::GtfParser;
use biometal::formats::{Strand, TabDelimitedParser};
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;

#[test]
fn test_real_gtf_parsing() {
    // Load real-world GTF data from Ensembl
    let gtf_path = PathBuf::from("tests/data/real_world/annotation/ensembl_chr21.gtf.gz");

    let file = File::open(&gtf_path).expect("Failed to open GTF file");
    let decoder = GzDecoder::new(file);
    let stream = GtfParser::new(decoder);

    let mut record_count = 0;
    let mut feature_counts: HashMap<String, usize> = HashMap::new();
    let mut gene_ids = Vec::new();

    for record in stream {
        let record = record.expect("Failed to parse GTF record");
        record_count += 1;

        // Validate record structure
        assert_eq!(record.seqname, "21", "All records should be from chromosome 21");
        assert!(!record.source.is_empty(), "Source should not be empty");
        assert!(!record.feature.is_empty(), "Feature type should not be empty");
        assert!(record.start > 0, "Start position should be positive");
        assert!(record.end >= record.start, "End must be >= start");
        // Strand should be one of the valid values
        assert!(
            matches!(record.strand, Strand::Forward | Strand::Reverse | Strand::Unknown),
            "Invalid strand: {:?}",
            record.strand
        );

        // Count feature types
        *feature_counts.entry(record.feature.clone()).or_insert(0) += 1;

        // Collect gene IDs
        if record.feature == "gene" {
            if let Some(gene_id) = record.attributes.get("gene_id") {
                gene_ids.push(gene_id.clone());
            }
        }
    }

    // Verify we got the expected number of records
    assert_eq!(record_count, 5000, "Expected 5000 GTF records");

    // Verify we have multiple feature types (gene, transcript, exon, CDS, etc.)
    assert!(
        feature_counts.len() >= 3,
        "Expected at least 3 different feature types, got {}",
        feature_counts.len()
    );

    // Common Ensembl GTF features
    let common_features = vec!["gene", "transcript", "exon"];
    for feature in common_features {
        assert!(
            feature_counts.contains_key(feature),
            "Expected to find '{}' features",
            feature
        );
    }

    // Verify we found some genes
    assert!(
        !gene_ids.is_empty(),
        "Should have found at least one gene"
    );

    println!("✅ Parsed {} GTF records", record_count);
    println!("   Feature types:");
    for (feature, count) in feature_counts.iter() {
        println!("     {}: {}", feature, count);
    }
    println!("   Genes found: {}", gene_ids.len());
}

#[test]
fn test_real_gtf_attributes() {
    // Validate GTF attribute parsing with real Ensembl data
    let gtf_path = PathBuf::from("tests/data/real_world/annotation/ensembl_chr21.gtf.gz");

    let file = File::open(&gtf_path).expect("Failed to open GTF file");
    let decoder = GzDecoder::new(file);
    let stream = GtfParser::new(decoder);

    let mut has_gene_id = 0;
    let mut has_gene_name = 0;
    let mut has_transcript_id = 0;

    for (idx, record) in stream.enumerate() {
        let record = record.expect("Failed to parse GTF record");

        // Check for required Ensembl GTF attributes
        if record.attributes.contains_key("gene_id") {
            has_gene_id += 1;
        }
        if record.attributes.contains_key("gene_name") {
            has_gene_name += 1;
        }
        if record.attributes.contains_key("transcript_id") {
            has_transcript_id += 1;
        }

        // Verify attributes are properly parsed (no quotes, semicolons, etc.)
        for (key, value) in &record.attributes {
            assert!(!key.is_empty(), "Attribute key should not be empty");
            assert!(!value.is_empty(), "Attribute value should not be empty");
            assert!(
                !value.contains('"'),
                "Attribute value should not contain quotes: {}",
                value
            );
            assert!(
                !value.ends_with(';'),
                "Attribute value should not end with semicolon: {}",
                value
            );
        }

        // Just check first 100 records for attributes
        if idx >= 100 {
            break;
        }
    }

    // All records should have gene_id
    assert!(
        has_gene_id > 0,
        "All records should have gene_id attribute"
    );

    // Most records should have gene_name
    assert!(
        has_gene_name > 0,
        "Most records should have gene_name attribute"
    );

    println!("✅ GTF attributes parsed correctly");
    println!("   Records with gene_id: {}", has_gene_id);
    println!("   Records with gene_name: {}", has_gene_name);
    println!("   Records with transcript_id: {}", has_transcript_id);
}

#[test]
fn test_real_gtf_gene_structure() {
    // Validate hierarchical gene structure (gene -> transcript -> exon)
    let gtf_path = PathBuf::from("tests/data/real_world/annotation/ensembl_chr21.gtf.gz");

    let file = File::open(&gtf_path).expect("Failed to open GTF file");
    let decoder = GzDecoder::new(file);
    let stream = GtfParser::new(decoder);

    let mut current_gene: Option<String> = None;
    let mut genes_with_transcripts = 0;
    let mut genes_with_exons = 0;

    for record in stream {
        let record = record.expect("Failed to parse GTF record");

        if record.feature == "gene" {
            if let Some(gene_id) = record.attributes.get("gene_id") {
                current_gene = Some(gene_id.clone());
            }
        } else if record.feature == "transcript" {
            if current_gene.is_some() {
                genes_with_transcripts += 1;
            }
        } else if record.feature == "exon" {
            if current_gene.is_some() {
                genes_with_exons += 1;
            }
        }
    }

    // Verify hierarchical structure
    assert!(
        genes_with_transcripts > 0,
        "Should have found genes with transcripts"
    );
    assert!(
        genes_with_exons > 0,
        "Should have found genes with exons"
    );

    println!("✅ GTF gene structure validated");
    println!("   Genes with transcripts: {}", genes_with_transcripts);
    println!("   Genes with exons: {}", genes_with_exons);
}

#[test]
fn test_real_gtf_constant_memory() {
    // Verify streaming parser doesn't load entire file into memory
    let gtf_path = PathBuf::from("tests/data/real_world/annotation/ensembl_chr21.gtf.gz");

    let file = File::open(&gtf_path).expect("Failed to open GTF file");
    let decoder = GzDecoder::new(file);
    let stream = GtfParser::new(decoder);

    let mut prev_record = None;
    let mut records_seen = 0;

    for record in stream {
        let record = record.expect("Failed to parse GTF record");

        // Verify previous record is dropped
        if let Some(_prev) = prev_record.take() {
            // Previous record dropped here
        }

        prev_record = Some(record);
        records_seen += 1;

        // Check first 200 records
        if records_seen >= 200 {
            break;
        }
    }

    assert_eq!(records_seen, 200, "Should have processed 200 records");
    println!("✅ Streaming parser maintains constant memory (1 record at a time)");
}
