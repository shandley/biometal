//! BLAST Tabular Format Parser
//!
//! BLAST (Basic Local Alignment Search Tool) tabular output provides alignment results
//! in a simple tab-delimited format, commonly used in bioinformatics pipelines.
//!
//! # Format Structure
//!
//! BLAST supports two tabular formats:
//! - **outfmt 6**: Tab-delimited, no headers (12 standard columns)
//! - **outfmt 7**: Tab-delimited with comment lines (starts with #)
//!
//! ## Standard 12 Columns
//!
//! 1. **qseqid**: Query sequence ID
//! 2. **sseqid**: Subject (reference) sequence ID
//! 3. **pident**: Percentage of identical matches (0-100)
//! 4. **length**: Alignment length
//! 5. **mismatch**: Number of mismatches
//! 6. **gapopen**: Number of gap openings
//! 7. **qstart**: Query sequence start position
//! 8. **qend**: Query sequence end position
//! 9. **sstart**: Subject sequence start position
//! 10. **send**: Subject sequence end position
//! 11. **evalue**: Expect value (statistical significance)
//! 12. **bitscore**: Bit score (alignment quality)
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::blast::BlastTabularParser;
//!
//! # fn main() -> biometal::Result<()> {
//! let parser = BlastTabularParser::from_path("alignments.blast")?;
//!
//! for record in parser {
//!     let record = record?;
//!     if record.evalue < 1e-10 && record.pident > 95.0 {
//!         println!("{} -> {}: {:.1}% identity, E={:.2e}",
//!             record.qseqid, record.sseqid, record.pident, record.evalue);
//!     }
//! }
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{
    fields::{parse_required, split_fields},
    TabDelimitedParser, TabDelimitedRecord, Result,
};

/// BLAST tabular alignment record (outfmt 6/7)
///
/// Represents a single alignment between a query and subject sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct BlastRecord {
    /// Query sequence ID
    pub qseqid: String,
    /// Subject (reference) sequence ID
    pub sseqid: String,
    /// Percentage of identical matches (0-100)
    pub pident: f64,
    /// Alignment length (number of aligned positions)
    pub length: u32,
    /// Number of mismatches
    pub mismatch: u32,
    /// Number of gap openings
    pub gapopen: u32,
    /// Query sequence start position (1-based)
    pub qstart: u32,
    /// Query sequence end position (1-based, inclusive)
    pub qend: u32,
    /// Subject sequence start position (1-based)
    pub sstart: u32,
    /// Subject sequence end position (1-based, inclusive)
    pub send: u32,
    /// Expect value (statistical significance, lower is better)
    pub evalue: f64,
    /// Bit score (alignment quality, higher is better)
    pub bitscore: f64,
}

impl TabDelimitedRecord for BlastRecord {
    fn from_line(line: &str) -> Result<Self> {
        // Skip comment lines in outfmt 7 (lines starting with #)
        if line.trim().starts_with('#') {
            return Err(crate::formats::primitives::FormatError::InvalidField {
                field: "line".to_string(),
                line: 0,
                reason: "Comment line should be skipped".to_string(),
            });
        }

        // BLAST tabular has 12 required fields
        let fields = split_fields(line, Some(12), 0)?;

        Ok(BlastRecord {
            qseqid: fields[0].to_string(),
            sseqid: fields[1].to_string(),
            pident: parse_required(fields[2], "pident", 0)?,
            length: parse_required(fields[3], "length", 0)?,
            mismatch: parse_required(fields[4], "mismatch", 0)?,
            gapopen: parse_required(fields[5], "gapopen", 0)?,
            qstart: parse_required(fields[6], "qstart", 0)?,
            qend: parse_required(fields[7], "qend", 0)?,
            sstart: parse_required(fields[8], "sstart", 0)?,
            send: parse_required(fields[9], "send", 0)?,
            evalue: parse_required(fields[10], "evalue", 0)?,
            bitscore: parse_required(fields[11], "bitscore", 0)?,
        })
    }

    fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qseqid,
            self.sseqid,
            self.pident,
            self.length,
            self.mismatch,
            self.gapopen,
            self.qstart,
            self.qend,
            self.sstart,
            self.send,
            self.evalue,
            self.bitscore
        )
    }
}

/// Streaming BLAST tabular parser (outfmt 6/7)
///
/// Reads BLAST tabular output with constant memory usage.
/// Supports both outfmt 6 (no headers) and outfmt 7 (with comment lines).
///
/// # Example
///
/// ```no_run
/// use biometal::formats::blast::BlastTabularParser;
///
/// # fn main() -> biometal::Result<()> {
/// let parser = BlastTabularParser::from_path("alignments.blast")?;
/// for result in parser {
///     let record = result?;
///     println!("{} -> {}: {:.1}% identity",
///         record.qseqid, record.sseqid, record.pident);
/// }
/// # Ok(())
/// # }
/// ```
pub type BlastTabularParser<R> = TabDelimitedParser<R, BlastRecord>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_blast_outfmt6() {
        let parser = BlastTabularParser::from_path("tests/data/blast/test_blastn.outfmt6")
            .expect("Failed to open test file");

        let records: Vec<_> = parser.collect();
        assert_eq!(records.len(), 5);

        // Check first record
        let record = records[0].as_ref().expect("Failed to parse record");
        assert_eq!(record.qseqid, "Query_1");
        assert_eq!(record.sseqid, "gi|123456|ref|NC_000001.1|");
        assert_eq!(record.pident, 100.00);
        assert_eq!(record.length, 100);
        assert_eq!(record.mismatch, 0);
        assert_eq!(record.gapopen, 0);
        assert_eq!(record.qstart, 1);
        assert_eq!(record.qend, 100);
        assert_eq!(record.sstart, 50000);
        assert_eq!(record.send, 50099);
        assert_eq!(record.evalue, 1e-50);
        assert_eq!(record.bitscore, 200.0);
    }

    #[test]
    fn test_parse_blast_outfmt7() {
        let parser = BlastTabularParser::from_path("tests/data/blast/test_blastn.outfmt7")
            .expect("Failed to open test file");

        let records: Vec<_> = parser.collect();
        assert_eq!(records.len(), 5); // Comment lines should be skipped

        // Check first record
        let record = records[0].as_ref().expect("Failed to parse record");
        assert_eq!(record.qseqid, "Query_1");
        assert_eq!(record.pident, 100.00);
    }

    #[test]
    fn test_blast_record_filtering() {
        let parser = BlastTabularParser::from_path("tests/data/blast/test_blastn.outfmt6")
            .expect("Failed to open test file");

        // Filter for high-quality alignments (evalue < 1e-30 AND pident > 95.0)
        let high_quality: Vec<_> = parser
            .filter_map(|r| r.ok())
            .filter(|r| r.evalue < 1e-30 && r.pident > 95.0)
            .collect();

        // Only 2 records pass: Query_1 first hit (1e-50, 100%) and Query_2 first hit (1e-95, 98.5%)
        assert_eq!(high_quality.len(), 2);

        // Check that filtered records meet criteria
        for record in &high_quality {
            assert!(record.evalue < 1e-30);
            assert!(record.pident > 95.0);
        }
    }

    #[test]
    fn test_blast_record_query_grouping() {
        let parser = BlastTabularParser::from_path("tests/data/blast/test_blastn.outfmt6")
            .expect("Failed to open test file");

        // Count hits per query
        let mut query_counts = std::collections::HashMap::new();
        for record in parser.filter_map(|r| r.ok()) {
            *query_counts.entry(record.qseqid).or_insert(0) += 1;
        }

        assert_eq!(query_counts.get("Query_1"), Some(&2));
        assert_eq!(query_counts.get("Query_2"), Some(&2));
        assert_eq!(query_counts.get("Query_3"), Some(&1));
    }

    #[test]
    fn test_parse_real_world_ecoli_blast() {
        // Test outfmt6 (no comments)
        let parser = BlastTabularParser::from_path("tests/data/real_world/alignments/blastn_ecoli_sample.outfmt6")
            .expect("Failed to open real E. coli BLAST file");

        let records: Vec<_> = parser.filter_map(|r| r.ok()).collect();
        assert_eq!(records.len(), 17, "Should have 17 BLAST hits");

        // Check first record (16S rRNA perfect match)
        let first = &records[0];
        assert_eq!(first.qseqid, "Query_16S_rRNA");
        assert_eq!(first.sseqid, "NR_102804.1");
        assert_eq!(first.pident, 99.82);
        assert_eq!(first.length, 564);
        assert_eq!(first.evalue, 0.0);
        assert_eq!(first.bitscore, 1037.0);

        // Check perfect match (gyrB)
        let perfect_match = records.iter().find(|r| r.pident == 100.0).expect("Should have perfect match");
        assert_eq!(perfect_match.qseqid, "Query_gyrB");
        assert_eq!(perfect_match.mismatch, 0);

        // Check weak hit (high E-value)
        let weak_hit = records.iter().filter(|r| r.evalue > 1.0).count();
        assert_eq!(weak_hit, 1, "Should have 1 weak hit with E>1");

        println!("✅ Successfully parsed real E. coli BLAST output: {} hits", records.len());
    }

    #[test]
    fn test_parse_real_world_ecoli_blast_outfmt7() {
        // Test outfmt7 (with comments)
        let parser = BlastTabularParser::from_path("tests/data/real_world/alignments/blastn_ecoli_sample.outfmt7")
            .expect("Failed to open real E. coli BLAST outfmt7 file");

        let records: Vec<_> = parser.filter_map(|r| r.ok()).collect();
        assert_eq!(records.len(), 17, "Should have 17 BLAST hits (comments skipped)");

        // Verify same data as outfmt6
        let first = &records[0];
        assert_eq!(first.qseqid, "Query_16S_rRNA");
        assert_eq!(first.pident, 99.82);

        println!("✅ Successfully parsed real E. coli BLAST outfmt7: {} hits", records.len());
    }
}
