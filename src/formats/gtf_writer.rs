//! GTF format writer with compression support
//!
//! # Format Specification
//!
//! GTF (Gene Transfer Format) is a 9-column tab-delimited format for gene annotations:
//! 1. seqname - Chromosome/contig name
//! 2. source - Annotation source (e.g., GENCODE, Ensembl)
//! 3. feature - Feature type (gene, transcript, exon, CDS, etc.)
//! 4. start - Start position (1-based, inclusive)
//! 5. end - End position (1-based, inclusive)
//! 6. score - Confidence score or "." if missing
//! 7. strand - +, -, or .
//! 8. frame - CDS frame (0, 1, 2) or "." if missing
//! 9. attributes - Semicolon-separated key "value" pairs (gene_id "ENSG001"; transcript_id "ENST001";)
//!
//! # GTF vs GFF3 Differences
//!
//! - **Attribute syntax**: GTF uses `gene_id "value";` vs GFF3's `ID=value`
//! - **Required attributes**: `gene_id` mandatory for all; `transcript_id` mandatory for transcript-level features
//! - **Feature types**: Limited set (gene, transcript, exon, CDS, UTR, etc.)
//!
//! # Architecture
//!
//! This writer follows biometal patterns:
//! - Automatic compression based on file extension (.gz, .bgz)
//! - cloudflare_zlib backend (1.67× faster decompression)
//! - Validation (coordinates, required fields, required attributes)
//! - Streaming write (constant memory)
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::gtf::GtfRecord;
//! use biometal::formats::gtf_writer::GtfWriter;
//! use biometal::formats::primitives::Strand;
//! use std::collections::HashMap;
//!
//! # fn main() -> biometal::Result<()> {
//! let mut writer = GtfWriter::create("annotations.gtf.gz")?;
//!
//! let mut attributes = HashMap::new();
//! attributes.insert("gene_id".to_string(), "ENSG00000223972".to_string());
//! attributes.insert("transcript_id".to_string(), "ENST00000456328".to_string());
//! attributes.insert("gene_name".to_string(), "DDX11L1".to_string());
//!
//! let record = GtfRecord {
//!     seqname: "chr1".to_string(),
//!     source: "HAVANA".to_string(),
//!     feature: "exon".to_string(),
//!     start: 11869,
//!     end: 12227,
//!     score: None,
//!     strand: Strand::Forward,
//!     frame: None,
//!     attributes,
//! };
//!
//! writer.write_record(&record)?;
//! writer.finish()?;
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::formats::gtf::GtfRecord;
use crate::formats::primitives::Strand;
use crate::io::compression::CompressedWriter;
use crate::io::sink::DataSink;
use std::io::Write;
use std::path::Path;

/// GTF format writer with compression support
///
/// # Features
///
/// - Automatic compression (gzip, bgzip) based on file extension
/// - Validation (coordinates, required fields, required attributes)
/// - Streaming write (constant memory)
/// - Proper GTF attribute formatting (key "value";)
///
/// # Example
///
/// ```no_run
/// use biometal::formats::gtf::GtfRecord;
/// use biometal::formats::gtf_writer::GtfWriter;
/// use biometal::formats::primitives::Strand;
/// use std::collections::HashMap;
///
/// # fn main() -> biometal::Result<()> {
/// let mut writer = GtfWriter::create("annotations.gtf.gz")?;
///
/// let mut attributes = HashMap::new();
/// attributes.insert("gene_id".to_string(), "ENSG001".to_string());
/// attributes.insert("transcript_id".to_string(), "ENST001".to_string());
///
/// let record = GtfRecord {
///     seqname: "chr1".to_string(),
///     source: "GENCODE".to_string(),
///     feature: "exon".to_string(),
///     start: 1000,
///     end: 2000,
///     score: None,
///     strand: Strand::Forward,
///     frame: None,
///     attributes,
/// };
///
/// writer.write_record(&record)?;
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct GtfWriter {
    writer: CompressedWriter,
    records_written: usize,
}

impl GtfWriter {
    /// Create a new GTF writer from a data sink
    ///
    /// Automatically detects compression from file extension:
    /// - `.gz` → gzip compression
    /// - `.bgz` → bgzip compression
    /// - other → uncompressed
    pub fn new(sink: DataSink) -> Result<Self> {
        let writer = CompressedWriter::new(sink)
            .map_err(|e| BiometalError::Io(e))?;
        Ok(Self {
            writer,
            records_written: 0,
        })
    }

    /// Create a GTF writer from a file path
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gtf_writer::GtfWriter;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = GtfWriter::create("output.gtf.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(DataSink::from_path(path))
    }

    /// Create a GTF writer to stdout
    ///
    /// Useful for streaming pipelines:
    /// ```bash
    /// biometal filter input.gtf.gz | bedtools intersect -a stdin -b regions.bed
    /// ```
    pub fn stdout() -> Result<Self> {
        Self::new(DataSink::stdout())
    }

    /// Write a single GTF record
    ///
    /// # Arguments
    ///
    /// * `record` - GTF record to write
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - seqname is empty
    /// - source is empty
    /// - feature is empty
    /// - start >= end
    /// - frame is invalid (not 0, 1, or 2)
    /// - gene_id attribute is missing (required for all GTF records)
    /// - transcript_id is missing for transcript-level features (exon, CDS, etc.)
    /// - An I/O error occurs
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::gtf_writer::GtfWriter;
    /// use biometal::formats::primitives::Strand;
    /// use std::collections::HashMap;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = GtfWriter::create("output.gtf")?;
    ///
    /// let mut attributes = HashMap::new();
    /// attributes.insert("gene_id".to_string(), "ENSG001".to_string());
    /// attributes.insert("gene_name".to_string(), "ABC1".to_string());
    ///
    /// let record = GtfRecord {
    ///     seqname: "chr1".to_string(),
    ///     source: "GENCODE".to_string(),
    ///     feature: "gene".to_string(),
    ///     start: 1000,
    ///     end: 2000,
    ///     score: Some(50.0),
    ///     strand: Strand::Forward,
    ///     frame: None,
    ///     attributes,
    /// };
    ///
    /// writer.write_record(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_record(&mut self, record: &GtfRecord) -> Result<()> {
        // Validate required fields
        if record.seqname.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "GTF: seqname cannot be empty".to_string(),
            });
        }

        if record.source.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "GTF: source cannot be empty".to_string(),
            });
        }

        if record.feature.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "GTF: feature cannot be empty".to_string(),
            });
        }

        if record.start >= record.end {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "GTF: Invalid interval: start ({}) >= end ({})",
                    record.start, record.end
                ),
            });
        }

        // Validate frame if present
        if let Some(frame) = record.frame {
            if frame > 2 {
                return Err(BiometalError::InvalidInput {
                    msg: format!("GTF: Invalid frame: {} (must be 0, 1, or 2)", frame),
                });
            }
        }

        // Validate required GTF attributes
        if !record.attributes.contains_key("gene_id") {
            return Err(BiometalError::InvalidInput {
                msg: "GTF: gene_id attribute is required for all records".to_string(),
            });
        }

        // transcript_id is required for transcript-level features (exon, CDS, etc.)
        // but optional for gene features
        let transcript_level_features = ["transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon"];
        if transcript_level_features.contains(&record.feature.as_str()) {
            if !record.attributes.contains_key("transcript_id") {
                return Err(BiometalError::InvalidInput {
                    msg: format!(
                        "GTF: transcript_id attribute is required for {} features",
                        record.feature
                    ),
                });
            }
        }

        // Format fields
        let score_str = record.score.map_or(".".to_string(), |s| format!("{:.2}", s));

        let strand_str = match record.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        };

        let frame_str = record.frame.map_or(".".to_string(), |f| f.to_string());

        // Format attributes as semicolon-separated key "value" pairs (GTF syntax)
        let attributes_str = if record.attributes.is_empty() {
            ".".to_string()
        } else {
            let mut attrs: Vec<String> = record.attributes
                .iter()
                .map(|(k, v)| format!("{} \"{}\"", k, v))
                .collect();
            // Sort for consistent output
            attrs.sort();
            // Join with "; " and add trailing semicolon
            format!("{};", attrs.join("; "))
        };

        // Write GTF line (9 columns)
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.seqname,
            record.source,
            record.feature,
            record.start,
            record.end,
            score_str,
            strand_str,
            frame_str,
            attributes_str
        )
        .map_err(|e| BiometalError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Write multiple GTF records from an iterator
    ///
    /// Convenience method for writing many records. The iterator can be
    /// any type that yields `Result<GtfRecord>`.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gtf_writer::GtfWriter;
    /// use biometal::formats::primitives::TabDelimitedParser;
    /// use biometal::formats::gtf::GtfRecord;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let file = File::open("input.gtf")?;
    /// let reader = BufReader::new(file);
    /// let parser = TabDelimitedParser::<_, GtfRecord>::new(reader);
    ///
    /// let mut writer = GtfWriter::create("output.gtf.gz")?;
    /// writer.write_all(parser)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_all<I>(&mut self, records: I) -> Result<()>
    where
        I: IntoIterator<Item = Result<GtfRecord>>,
    {
        for record in records {
            self.write_record(&record?)?;
        }
        Ok(())
    }

    /// Get the number of records written so far
    pub fn records_written(&self) -> usize {
        self.records_written
    }

    /// Flush buffered data to disk
    ///
    /// You typically don't need to call this explicitly as `finish()`
    /// will flush automatically. However, it can be useful for
    /// long-running processes to ensure data is persisted.
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()
            .map_err(|e| BiometalError::Io(e))
    }

    /// Finish writing and flush all data
    ///
    /// This method MUST be called to ensure all data is written to disk.
    /// It flushes the internal buffers and closes the compression stream.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::gtf::GtfRecord;
    /// use biometal::formats::gtf_writer::GtfWriter;
    /// use biometal::formats::primitives::Strand;
    /// use std::collections::HashMap;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = GtfWriter::create("output.gtf.gz")?;
    ///
    /// let mut attributes = HashMap::new();
    /// attributes.insert("gene_id".to_string(), "ENSG001".to_string());
    ///
    /// let record = GtfRecord {
    ///     seqname: "chr1".to_string(),
    ///     source: "GENCODE".to_string(),
    ///     feature: "gene".to_string(),
    ///     start: 1000,
    ///     end: 2000,
    ///     score: None,
    ///     strand: Strand::Forward,
    ///     frame: None,
    ///     attributes,
    /// };
    ///
    /// writer.write_record(&record)?;
    /// writer.finish()?;  // IMPORTANT: Flush and close
    /// # Ok(())
    /// # }
    /// ```
    pub fn finish(mut self) -> Result<()> {
        self.writer.flush()
            .map_err(|e| BiometalError::Io(e))?;
        Ok(())
    }
}
