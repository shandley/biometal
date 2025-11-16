//! GTF (Gene Transfer Format) parser.
//!
//! GTF is a widely-used format for gene annotations, especially in RNA-seq workflows:
//! - **Gene models**: Gene and transcript structures
//! - **Exons and CDS**: Coding and non-coding regions
//! - **UTRs**: Untranslated regions (5' and 3')
//!
//! # Format Specification
//!
//! GTF uses 9 tab-delimited columns (same structure as GFF):
//! 1. **seqname**: Chromosome/contig name
//! 2. **source**: Annotation source (e.g., GENCODE, Ensembl)
//! 3. **feature**: Feature type (gene, transcript, exon, CDS, etc.)
//! 4. **start**: Start position (1-based, inclusive)
//! 5. **end**: End position (1-based, inclusive)
//! 6. **score**: Confidence score (or `.`)
//! 7. **strand**: `+` or `-` (`.` not typically used)
//! 8. **frame**: CDS phase (0, 1, 2, or `.`)
//! 9. **attributes**: Semicolon-separated `key "value"` pairs
//!
//! # GTF vs GFF3 Differences
//!
//! - **Attribute syntax**: GTF uses `gene_id "value";` vs GFF3's `ID=value`
//! - **Required attributes**: `gene_id` mandatory for all; `transcript_id` mandatory for transcript-level features (exon, CDS, etc.), optional for gene features
//! - **Feature types**: Limited set (CDS, exon, start_codon, stop_codon, UTR, etc.)
//! - **Hierarchy**: Implicit via `gene_id`/`transcript_id` (not explicit Parent/ID)
//!
//! # Examples
//!
//! ## Basic GTF record
//!
//! ```
//! use biometal::formats::gtf::GtfRecord;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = r#"chr1    HAVANA    exon    11869    12227    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";"#;
//! let record = GtfRecord::from_line(line)?;
//!
//! assert_eq!(record.seqname, "chr1");
//! assert_eq!(record.feature, "exon");
//! assert_eq!(record.start, 11869);
//! assert_eq!(record.end, 12227);
//! assert_eq!(record.gene_id(), "ENSG00000223972");
//! assert_eq!(record.transcript_id(), Some("ENST00000456328"));
//! # Ok(())
//! # }
//! ```
//!
//! ## CDS with frame
//!
//! ```
//! use biometal::formats::gtf::GtfRecord;
//! use biometal::formats::TabDelimitedRecord;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let line = r#"chr1    HAVANA    CDS    12227    12612    .    +    0    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";"#;
//! let record = GtfRecord::from_line(line)?;
//!
//! assert_eq!(record.feature, "CDS");
//! assert_eq!(record.frame, Some(0));
//! # Ok(())
//! # }
//! ```
//!
//! ## Streaming parser
//!
//! ```no_run
//! use biometal::formats::gtf::GtfParser;
//! use std::fs::File;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let file = File::open("annotations.gtf")?;
//! let parser = GtfParser::new(file);
//!
//! for result in parser {
//!     let record = result?;
//!     if record.feature == "gene" {
//!         println!("Gene: {} at {}:{}-{}",
//!                  record.gene_id(),
//!                  record.seqname, record.start, record.end);
//!     }
//! }
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{
    fields::{parse_optional, parse_required, split_fields},
    FormatError, GenomicInterval, Result, Strand, TabDelimitedRecord,
};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::str::FromStr;

/// GTF feature record.
///
/// Represents a single genomic feature (gene, transcript, exon, CDS, etc.).
///
/// # Examples
///
/// ```
/// use biometal::formats::gtf::GtfRecord;
/// use biometal::formats::TabDelimitedRecord;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let line = r#"chr1    HAVANA    exon    11869    12227    .    +    .    gene_id "ENSG00000223972"; transcript_id "ENST00000456328";"#;
/// let record = GtfRecord::from_line(line)?;
/// assert_eq!(record.feature, "exon");
/// assert_eq!(record.gene_id(), "ENSG00000223972");
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct GtfRecord {
    /// Chromosome/contig name
    pub seqname: String,
    /// Annotation source (e.g., GENCODE, Ensembl, HAVANA)
    pub source: String,
    /// Feature type (gene, transcript, exon, CDS, UTR, etc.)
    pub feature: String,
    /// Start position (1-based, inclusive)
    pub start: u64,
    /// End position (1-based, inclusive)
    pub end: u64,
    /// Confidence score (None if `.`)
    pub score: Option<f64>,
    /// Strand (+, -)
    pub strand: Strand,
    /// CDS frame/phase (0, 1, 2, or None if `.`)
    pub frame: Option<u8>,
    /// Attributes (key-value pairs from GTF attribute column)
    pub attributes: HashMap<String, String>,
}

impl GtfRecord {
    /// Get gene_id (required attribute in GTF)
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let line = r#"chr1    .    exon    100    200    .    +    .    gene_id "GENE1"; transcript_id "TX1";"#;
    /// let record = GtfRecord::from_line(line)?;
    /// assert_eq!(record.gene_id(), "GENE1");
    /// # Ok(())
    /// # }
    /// ```
    pub fn gene_id(&self) -> &str {
        // Safe: validated in from_line()
        &self.attributes["gene_id"]
    }

    /// Get transcript_id (required for transcript-level features, optional for gene features)
    ///
    /// Returns `None` for gene features, `Some(&str)` for transcript/exon/CDS features.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let line = r#"chr1    .    exon    100    200    .    +    .    gene_id "GENE1"; transcript_id "TX1";"#;
    /// let record = GtfRecord::from_line(line)?;
    /// assert_eq!(record.transcript_id(), Some("TX1"));
    ///
    /// let gene_line = r#"chr1    .    gene    100    200    .    +    .    gene_id "GENE1";"#;
    /// let gene_record = GtfRecord::from_line(gene_line)?;
    /// assert_eq!(gene_record.transcript_id(), None);
    /// # Ok(())
    /// # }
    /// ```
    pub fn transcript_id(&self) -> Option<&str> {
        self.attributes.get("transcript_id").map(|s| s.as_str())
    }

    /// Get optional gene_name attribute
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let line = r#"chr1    .    gene    100    200    .    +    .    gene_id "ENSG001"; transcript_id "ENST001"; gene_name "ABC1";"#;
    /// let record = GtfRecord::from_line(line)?;
    /// assert_eq!(record.gene_name(), Some("ABC1"));
    /// # Ok(())
    /// # }
    /// ```
    pub fn gene_name(&self) -> Option<&str> {
        self.attributes.get("gene_name").map(|s| s.as_str())
    }

    /// Get optional gene_biotype attribute
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let line = r#"chr1    .    gene    100    200    .    +    .    gene_id "ENSG001"; transcript_id "ENST001"; gene_biotype "protein_coding";"#;
    /// let record = GtfRecord::from_line(line)?;
    /// assert_eq!(record.gene_biotype(), Some("protein_coding"));
    /// # Ok(())
    /// # }
    /// ```
    pub fn gene_biotype(&self) -> Option<&str> {
        self.attributes.get("gene_biotype").map(|s| s.as_str())
    }

    /// Get genomic interval in 0-based half-open coordinates [start, end)
    ///
    /// Converts from GTF's 1-based inclusive coordinates to standard 0-based.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let line = r#"chr1    .    exon    100    200    .    +    .    gene_id "G1"; transcript_id "T1";"#;
    /// let record = GtfRecord::from_line(line)?;
    /// let interval = record.interval()?;
    ///
    /// assert_eq!(interval.chrom, "chr1");
    /// assert_eq!(interval.start, 99);  // 0-based
    /// assert_eq!(interval.end, 200);   // half-open
    /// # Ok(())
    /// # }
    /// ```
    pub fn interval(&self) -> Result<GenomicInterval> {
        // GTF uses 1-based inclusive coordinates
        // Convert to 0-based half-open: [start-1, end)
        GenomicInterval::new(self.seqname.clone(), self.start - 1, self.end)
    }

    /// Calculate feature length in base pairs
    ///
    /// GTF uses 1-based inclusive coordinates, so length = end - start + 1
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::TabDelimitedRecord;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let line = r#"chr1    .    exon    1000    2000    .    +    .    gene_id "G1"; transcript_id "T1";"#;
    /// let record = GtfRecord::from_line(line)?;
    /// assert_eq!(record.length(), 1001); // 1000-2000 inclusive
    /// # Ok(())
    /// # }
    /// ```
    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }
}

impl TabDelimitedRecord for GtfRecord {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(9), 0)?;

        let seqname = fields[0].to_string();
        let source = fields[1].to_string();
        let feature = fields[2].to_string();
        let start: u64 = parse_required(fields[3], "start", 0)?;
        let end: u64 = parse_required(fields[4], "end", 0)?;
        let score: Option<f64> = parse_optional(fields[5], "score", 0)?;
        let strand = Strand::from_str(fields[6])?;

        // Parse frame (0, 1, 2, or .)
        let frame: Option<u8> = if fields[7] == "." {
            None
        } else {
            Some(parse_required(fields[7], "frame", 0)?)
        };

        // Parse GTF attributes (different format from GFF3)
        let attributes = parse_gtf_attributes(fields[8])?;

        // Validate required attributes
        if !attributes.contains_key("gene_id") {
            return Err(FormatError::InvalidField {
                field: "attributes".to_string(),
                line: 0,
                reason: "Missing required attribute 'gene_id'".to_string(),
            });
        }

        // transcript_id is only required for transcript-level features (not gene features)
        // Per GTF spec: gene features don't have transcript_id
        if feature != "gene" && !attributes.contains_key("transcript_id") {
            return Err(FormatError::InvalidField {
                field: "attributes".to_string(),
                line: 0,
                reason: format!("Missing required attribute 'transcript_id' for feature type '{}'", feature),
            });
        }

        Ok(GtfRecord {
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            frame,
            attributes,
        })
    }

    fn to_line(&self) -> String {
        let score_str = self
            .score
            .map(|s| s.to_string())
            .unwrap_or_else(|| ".".to_string());
        let frame_str = self
            .frame
            .map(|f| f.to_string())
            .unwrap_or_else(|| ".".to_string());
        let strand_str = match self.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        };

        // Format attributes in GTF style
        let mut attr_vec: Vec<_> = self.attributes.iter().collect();
        attr_vec.sort_by_key(|(k, _)| *k); // Deterministic order
        let attrs = attr_vec
            .iter()
            .map(|(k, v)| format!("{} \"{}\";", k, v))
            .collect::<Vec<_>>()
            .join(" ");

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seqname,
            self.source,
            self.feature,
            self.start,
            self.end,
            score_str,
            strand_str,
            frame_str,
            attrs
        )
    }
}

/// Parse GTF attribute string.
///
/// GTF attributes use the format: `key "value"; key2 "value2";`
/// Each attribute ends with a semicolon followed by a space.
/// Values are enclosed in double quotes.
///
/// # Format
///
/// ```text
/// gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1";
/// ```
fn parse_gtf_attributes(attr_str: &str) -> Result<HashMap<String, String>> {
    let mut attributes = HashMap::new();

    if attr_str == "." || attr_str.trim().is_empty() {
        return Ok(attributes);
    }

    // Split by semicolon
    for part in attr_str.split(';') {
        let trimmed = part.trim();
        if trimmed.is_empty() {
            continue;
        }

        // Find the first space (separates key from quoted value)
        if let Some(space_pos) = trimmed.find(' ') {
            let key = &trimmed[..space_pos];
            let value_part = trimmed[space_pos + 1..].trim();

            // Remove quotes from value
            let value = if value_part.starts_with('"') && value_part.ends_with('"') {
                &value_part[1..value_part.len() - 1]
            } else {
                // Handle case without quotes (non-standard but sometimes seen)
                value_part
            };

            attributes.insert(key.to_string(), value.to_string());
        } else {
            // Malformed attribute (no space separator)
            return Err(FormatError::InvalidField {
                field: "attributes".to_string(),
                line: 0,
                reason: format!("Malformed GTF attribute: {}", trimmed),
            });
        }
    }

    Ok(attributes)
}

/// GTF parser with streaming support.
///
/// Iterates through GTF records one at a time for constant memory usage.
///
/// # Examples
///
/// ```no_run
/// use biometal::formats::gtf::GtfParser;
/// use std::fs::File;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let file = File::open("annotations.gtf")?;
/// let parser = GtfParser::new(file);
///
/// let mut gene_count = 0;
/// for result in parser {
///     let record = result?;
///     if record.feature == "gene" {
///         gene_count += 1;
///     }
///}
/// println!("Found {} genes", gene_count);
/// # Ok(())
/// # }
/// ```
pub struct GtfParser<R: Read> {
    reader: BufReader<R>,
    line_buffer: String,
}

impl<R: Read> GtfParser<R> {
    /// Create a new GTF parser
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gtf::GtfParser;
    /// use std::fs::File;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let file = File::open("annotations.gtf")?;
    /// let parser = GtfParser::new(file);
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(reader: R) -> Self {
        GtfParser {
            reader: BufReader::new(reader),
            line_buffer: String::new(),
        }
    }
}

impl<R: Read> Iterator for GtfParser<R> {
    type Item = Result<GtfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    let line = self.line_buffer.trim();

                    // Skip empty lines and comments
                    if line.is_empty() || line.starts_with('#') {
                        continue;
                    }

                    return Some(GtfRecord::from_line(line));
                }
                Err(e) => return Some(Err(e.into()))
            }
        }
    }
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_parse_gtf_attributes() {
		let attr_str = r#"gene_id "ENSG001"; transcript_id "ENST001"; gene_name "ABC1";"#;
		let attrs = parse_gtf_attributes(attr_str).unwrap();

		assert_eq!(attrs.get("gene_id"), Some(&"ENSG001".to_string()));
		assert_eq!(attrs.get("transcript_id"), Some(&"ENST001".to_string()));
		assert_eq!(attrs.get("gene_name"), Some(&"ABC1".to_string()));
	}

	#[test]
	fn test_gtf_record_parse() {
		let line = r#"chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328";"#;
		let record = GtfRecord::from_line(line).unwrap();

		assert_eq!(record.seqname, "chr1");
		assert_eq!(record.source, "HAVANA");
		assert_eq!(record.feature, "exon");
		assert_eq!(record.start, 11869);
		assert_eq!(record.end, 12227);
		assert_eq!(record.strand, Strand::Forward);
		assert_eq!(record.gene_id(), "ENSG00000223972");
		assert_eq!(record.transcript_id(), Some("ENST00000456328"));
	}

	#[test]
	fn test_gtf_cds_with_frame() {
		let line = r#"chr1	HAVANA	CDS	12227	12612	.	+	0	gene_id "ENSG00000223972"; transcript_id "ENST00000456328";"#;
		let record = GtfRecord::from_line(line).unwrap();

		assert_eq!(record.feature, "CDS");
		assert_eq!(record.frame, Some(0));
	}

	#[test]
	fn test_gtf_round_trip() {
		let original = r#"chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG001"; transcript_id "ENST001";"#;
		let record = GtfRecord::from_line(original).unwrap();
		let serialized = record.to_line();

		// Parse it again
		let record2 = GtfRecord::from_line(&serialized).unwrap();
		assert_eq!(record, record2);
	}

	#[test]
	fn test_gtf_interval_conversion() {
		let line = r#"chr1	.	exon	100	200	.	+	.	gene_id "G1"; transcript_id "T1";"#;
		let record = GtfRecord::from_line(line).unwrap();
		let interval = record.interval().unwrap();

		assert_eq!(interval.start, 99); // 1-based → 0-based
		assert_eq!(interval.end, 200); // inclusive → half-open
	}

	#[test]
	fn test_gtf_missing_values() {
		let line = r#"chr1	.	gene	100	200	.	.	.	gene_id "G1"; transcript_id "T1";"#;
		let record = GtfRecord::from_line(line).unwrap();

		assert_eq!(record.source, ".");
		assert_eq!(record.score, None);
		assert_eq!(record.strand, Strand::Unknown);
		assert_eq!(record.frame, None);
	}

	#[test]
	fn test_gtf_optional_attributes() {
		let line = r#"chr1	.	gene	100	200	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "ABC1"; gene_biotype "protein_coding";"#;
		let record = GtfRecord::from_line(line).unwrap();

		assert_eq!(record.gene_name(), Some("ABC1"));
		assert_eq!(record.gene_biotype(), Some("protein_coding"));
	}

	#[test]
	fn test_gtf_different_features() {
		let features = vec!["gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon"];

		for feature in features {
			let line = format!(
				r#"chr1	.	{}	100	200	.	+	.	gene_id "G1"; transcript_id "T1";"#,
				feature
			);
			let record = GtfRecord::from_line(&line).unwrap();
			assert_eq!(record.feature, feature);
		}
	}

	#[test]
	fn test_gtf_score_and_frame() {
		let line = r#"chr1	.	CDS	100	200	95.5	+	2	gene_id "G1"; transcript_id "T1";"#;
		let record = GtfRecord::from_line(line).unwrap();

		assert_eq!(record.score, Some(95.5));
		assert_eq!(record.frame, Some(2));
	}

	#[test]
	fn test_gtf_length() {
		let line = r#"chr1	.	exon	1000	2000	.	+	.	gene_id "G1"; transcript_id "T1";"#;
		let record = GtfRecord::from_line(line).unwrap();

		// 1000-2000 inclusive = 1001 bp
		assert_eq!(record.length(), 1001);
	}

	#[test]
	fn test_gtf_parser_iteration() {
		let gtf_data = r#"chr1	.	gene	100	200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	.	transcript	100	200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	.	exon	100	150	.	+	.	gene_id "G1"; transcript_id "T1";
"#;

		let cursor = std::io::Cursor::new(gtf_data);
		let parser = GtfParser::new(cursor);

		let records: Vec<_> = parser.collect::<Result<Vec<_>>>().unwrap();
		assert_eq!(records.len(), 3);
		assert_eq!(records[0].feature, "gene");
		assert_eq!(records[1].feature, "transcript");
		assert_eq!(records[2].feature, "exon");
	}

	#[test]
	fn test_gtf_skip_comments() {
		let gtf_data = r#"# This is a comment
chr1	.	gene	100	200	.	+	.	gene_id "G1"; transcript_id "T1";
## Another comment
chr1	.	exon	100	150	.	+	.	gene_id "G1"; transcript_id "T1";
"#;

		let cursor = std::io::Cursor::new(gtf_data);
		let parser = GtfParser::new(cursor);

		let records: Vec<_> = parser.collect::<Result<Vec<_>>>().unwrap();
		assert_eq!(records.len(), 2);
	}

	#[test]
	fn test_gtf_attributes_without_quotes() {
		// Some GTF files have unquoted values (non-standard but seen in practice)
		let attr_str = r#"gene_id ENSG001; transcript_id ENST001;"#;
		let attrs = parse_gtf_attributes(attr_str).unwrap();

		assert_eq!(attrs.get("gene_id"), Some(&"ENSG001".to_string()));
		assert_eq!(attrs.get("transcript_id"), Some(&"ENST001".to_string()));
	}

	#[test]
	fn test_gtf_empty_attributes() {
		let attr_str = ".";
		let attrs = parse_gtf_attributes(attr_str).unwrap();
		assert!(attrs.is_empty());
	}

	#[test]
	fn test_gtf_attributes_trailing_semicolon() {
		let attr_str = r#"gene_id "G1"; transcript_id "T1";"#;
		let attrs = parse_gtf_attributes(attr_str).unwrap();

		assert_eq!(attrs.get("gene_id"), Some(&"G1".to_string()));
		assert_eq!(attrs.get("transcript_id"), Some(&"T1".to_string()));
	}

	#[test]
	fn test_gtf_strand_variations() {
		let strands = vec![
			("+", Strand::Forward),
			("-", Strand::Reverse),
			(".", Strand::Unknown),
		];

		for (strand_str, expected_strand) in strands {
			let line = format!(
				r#"chr1	.	gene	100	200	.	{}	.	gene_id "G1"; transcript_id "T1";"#,
				strand_str
			);
			let record = GtfRecord::from_line(&line).unwrap();
			assert_eq!(record.strand, expected_strand);
		}
	}

	#[test]
	fn test_gtf_required_attributes_present() {
		let line = r#"chr1	.	exon	100	200	.	+	.	gene_id "G1"; transcript_id "T1";"#;
		let record = GtfRecord::from_line(line).unwrap();

		// Required attributes should always be accessible
		assert_eq!(record.gene_id(), "G1");
		assert_eq!(record.transcript_id(), Some("T1"));
	}

	#[test]
	fn test_gtf_missing_gene_id() {
		// GTF record without gene_id should fail validation
		let line = "chr1\t.\texon\t100\t200\t.\t+\t.\ttranscript_id \"T1\";";
		let result = GtfRecord::from_line(line);
		assert!(result.is_err());

		let err = result.unwrap_err();
		assert!(err.to_string().contains("gene_id"));
	}

	#[test]
	fn test_gtf_missing_transcript_id() {
		// GTF record without transcript_id should fail validation
		let line = "chr1\t.\texon\t100\t200\t.\t+\t.\tgene_id \"G1\";";
		let result = GtfRecord::from_line(line);
		assert!(result.is_err());

		let err = result.unwrap_err();
		assert!(err.to_string().contains("transcript_id"));
	}

	#[test]
	fn test_gtf_missing_both_required_attributes() {
		// GTF record without both required attributes should fail
		let line = "chr1\t.\texon\t100\t200\t.\t+\t.\tgene_name \"ABC\";";
		let result = GtfRecord::from_line(line);
		assert!(result.is_err());

		// Should fail on gene_id first
		let err = result.unwrap_err();
		assert!(err.to_string().contains("gene_id"));
	}

	#[test]
	fn test_gtf_gene_without_transcript_id() {
		// Gene features should be valid without transcript_id
		let line = r#"chr1	havana	gene	100	500	.	+	.	gene_id "ENSG00000001"; gene_name "ABC1";"#;
		let record = GtfRecord::from_line(line).unwrap();

		assert_eq!(record.feature, "gene");
		assert_eq!(record.gene_id(), "ENSG00000001");
		assert_eq!(record.transcript_id(), None);  // Gene features don't have transcript_id
		assert_eq!(record.gene_name(), Some("ABC1"));
	}
}
