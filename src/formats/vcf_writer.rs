//! VCF format writer with compression support
//!
//! # Format Specification
//!
//! VCF (Variant Call Format) is a text format for storing genetic variants:
//! - **Header lines** (##): Metadata (fileformat, INFO, FORMAT, FILTER, contig)
//! - **Column header** (#): Field names and sample IDs
//! - **Data records**: Tab-delimited variant calls (8-N columns)
//!
//! # VCF Structure
//!
//! ```text
//! ##fileformat=VCFv4.2
//! ##INFO=<ID=DP,Description="Total Depth">
//! ##FORMAT=<ID=GT,Description="Genotype">
//! #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1
//! chr1    100  .   A    T    30    PASS    DP=50 GT      0/1
//! ```
//!
//! # Architecture
//!
//! This writer follows biometal patterns:
//! - Automatic compression based on file extension (.gz, .bgz)
//! - cloudflare_zlib backend (1.67× faster decompression)
//! - Header management (automatic or explicit)
//! - Validation (coordinates, required fields)
//! - Streaming write (constant memory)
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::vcf::{VcfHeader, VcfRecord};
//! use biometal::formats::vcf_writer::VcfWriter;
//! use std::collections::HashMap;
//!
//! # fn main() -> biometal::Result<()> {
//! // Create header
//! let mut header = VcfHeader::new("VCFv4.2".to_string());
//! header.info_fields.insert("DP".to_string(), "Total Depth".to_string());
//! header.samples = vec!["sample1".to_string()];
//!
//! // Create writer
//! let mut writer = VcfWriter::create("variants.vcf.gz")?;
//! writer.write_header(&header)?;
//!
//! // Write variant
//! let mut info = HashMap::new();
//! info.insert("DP".to_string(), "50".to_string());
//!
//! let record = VcfRecord {
//!     chrom: "chr1".to_string(),
//!     pos: 100,
//!     id: None,
//!     reference: "A".to_string(),
//!     alternate: vec!["T".to_string()],
//!     quality: Some(30.0),
//!     filter: Some("PASS".to_string()),
//!     info,
//!     format: Some("GT".to_string()),
//!     samples: vec!["0/1".to_string()],
//! };
//!
//! writer.write_record(&record)?;
//! writer.finish()?;
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::formats::primitives::TabDelimitedRecord;
use crate::formats::vcf::{VcfHeader, VcfRecord};
use crate::io::compression::CompressedWriter;
use crate::io::sink::DataSink;
use std::io::Write;
use std::path::Path;

/// VCF format writer with compression support
///
/// # Features
///
/// - Automatic compression (gzip, bgzip) based on file extension
/// - Header management (automatic or explicit)
/// - Validation (coordinates, required fields)
/// - Streaming write (constant memory)
/// - Proper VCF 4.2 formatting
///
/// # Example
///
/// ```no_run
/// use biometal::formats::vcf::{VcfHeader, VcfRecord};
/// use biometal::formats::vcf_writer::VcfWriter;
/// use std::collections::HashMap;
///
/// # fn main() -> biometal::Result<()> {
/// let mut header = VcfHeader::new("VCFv4.2".to_string());
/// let mut writer = VcfWriter::create("variants.vcf.gz")?;
/// writer.write_header(&header)?;
///
/// let mut info = HashMap::new();
/// info.insert("DP".to_string(), "100".to_string());
///
/// let record = VcfRecord {
///     chrom: "chr1".to_string(),
///     pos: 12345,
///     id: Some("rs123".to_string()),
///     reference: "A".to_string(),
///     alternate: vec!["T".to_string()],
///     quality: Some(30.0),
///     filter: Some("PASS".to_string()),
///     info,
///     format: None,
///     samples: vec![],
/// };
///
/// writer.write_record(&record)?;
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct VcfWriter {
    writer: CompressedWriter,
    records_written: usize,
    header_written: bool,
}

impl VcfWriter {
    /// Create a new VCF writer from a data sink
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
            header_written: false,
        })
    }

    /// Create a VCF writer from a file path
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::vcf_writer::VcfWriter;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut writer = VcfWriter::create("output.vcf.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(DataSink::from_path(path))
    }

    /// Create a VCF writer to stdout
    ///
    /// Useful for streaming pipelines:
    /// ```bash
    /// biometal filter input.vcf.gz | bcftools view -
    /// ```
    pub fn stdout() -> Result<Self> {
        Self::new(DataSink::stdout())
    }

    /// Write VCF header
    ///
    /// This writes all header lines (##) and the column header (#CHROM...).
    /// Must be called before writing records.
    ///
    /// # Arguments
    ///
    /// * `header` - VCF header with metadata and samples
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::vcf::VcfHeader;
    /// use biometal::formats::vcf_writer::VcfWriter;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut header = VcfHeader::new("VCFv4.2".to_string());
    /// header.samples = vec!["sample1".to_string()];
    ///
    /// let mut writer = VcfWriter::create("output.vcf")?;
    /// writer.write_header(&header)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_header(&mut self, header: &VcfHeader) -> Result<()> {
        if self.header_written {
            return Err(BiometalError::InvalidInput {
                msg: "VCF header already written".to_string(),
            });
        }

        // Get all header lines (includes column header)
        let header_lines = header.to_header_lines();

        // Write each line
        for line in header_lines {
            writeln!(self.writer, "{}", line)
                .map_err(|e| BiometalError::Io(e))?;
        }

        self.header_written = true;
        Ok(())
    }

    /// Write a single VCF record
    ///
    /// # Arguments
    ///
    /// * `record` - VCF record to write
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Header has not been written
    /// - chrom is empty
    /// - reference allele is empty
    /// - pos is 0 (VCF is 1-based)
    /// - alternate alleles are empty
    /// - An I/O error occurs
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::vcf::{VcfHeader, VcfRecord};
    /// use biometal::formats::vcf_writer::VcfWriter;
    /// use std::collections::HashMap;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let header = VcfHeader::new("VCFv4.2".to_string());
    /// let mut writer = VcfWriter::create("output.vcf")?;
    /// writer.write_header(&header)?;
    ///
    /// let mut info = HashMap::new();
    /// info.insert("DP".to_string(), "50".to_string());
    ///
    /// let record = VcfRecord {
    ///     chrom: "chr1".to_string(),
    ///     pos: 100,
    ///     id: None,
    ///     reference: "A".to_string(),
    ///     alternate: vec!["T".to_string()],
    ///     quality: None,
    ///     filter: Some("PASS".to_string()),
    ///     info,
    ///     format: None,
    ///     samples: vec![],
    /// };
    ///
    /// writer.write_record(&record)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_record(&mut self, record: &VcfRecord) -> Result<()> {
        // Ensure header was written
        if !self.header_written {
            return Err(BiometalError::InvalidInput {
                msg: "VCF: Header must be written before records".to_string(),
            });
        }

        // Validate required fields
        if record.chrom.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "VCF: chrom cannot be empty".to_string(),
            });
        }

        if record.pos == 0 {
            return Err(BiometalError::InvalidInput {
                msg: "VCF: pos must be >= 1 (1-based coordinates)".to_string(),
            });
        }

        if record.reference.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "VCF: reference allele cannot be empty".to_string(),
            });
        }

        if record.alternate.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "VCF: alternate alleles cannot be empty (use '.' for reference-only sites)".to_string(),
            });
        }

        // Use VcfRecord's to_line() method for serialization
        let line = record.to_line();

        writeln!(self.writer, "{}", line)
            .map_err(|e| BiometalError::Io(e))?;

        self.records_written += 1;
        Ok(())
    }

    /// Write multiple VCF records from an iterator
    ///
    /// Convenience method for writing many records. The iterator can be
    /// any type that yields `Result<VcfRecord>`.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::vcf_writer::VcfWriter;
    /// use biometal::formats::vcf::{VcfParser, VcfHeader};
    /// use std::fs::File;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let file = File::open("input.vcf")?;
    /// let mut parser = VcfParser::new(file);
    /// let header = parser.parse_header()?;
    ///
    /// let mut writer = VcfWriter::create("output.vcf.gz")?;
    /// writer.write_header(&header)?;
    /// writer.write_all(parser)?;
    /// writer.finish()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write_all<I>(&mut self, records: I) -> Result<()>
    where
        I: IntoIterator<Item = Result<VcfRecord>>,
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

    /// Check if header has been written
    pub fn header_written(&self) -> bool {
        self.header_written
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
    /// use biometal::formats::vcf::{VcfHeader, VcfRecord};
    /// use biometal::formats::vcf_writer::VcfWriter;
    /// use std::collections::HashMap;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let header = VcfHeader::new("VCFv4.2".to_string());
    /// let mut writer = VcfWriter::create("output.vcf.gz")?;
    /// writer.write_header(&header)?;
    ///
    /// let record = VcfRecord {
    ///     chrom: "chr1".to_string(),
    ///     pos: 100,
    ///     id: None,
    ///     reference: "A".to_string(),
    ///     alternate: vec!["T".to_string()],
    ///     quality: None,
    ///     filter: None,
    ///     info: HashMap::new(),
    ///     format: None,
    ///     samples: vec![],
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
