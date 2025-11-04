//! FASTQ streaming parser implementing Rules 2 and 5
//!
//! # Evidence Base
//!
//! - **Rule 2**: Block-based processing (10K records, Entry 027)
//!   - Record-by-record NEON = 82-86% overhead
//!   - Block-based (10K) = 4-8% overhead
//!   - Critical for preserving NEON speedup
//!
//! - **Rule 5**: Constant-memory streaming (~5 MB, Entry 026)
//!   - 99.5% memory reduction (1,344 MB → 5 MB @ 1M sequences)
//!   - Memory is constant regardless of dataset size
//!   - Enables 5TB analysis on consumer laptops

use crate::error::{BiometalError, Result};
use crate::io::compression::{CompressedReader, DataSource};
use crate::types::FastqRecord;
use std::io::BufRead;
use std::path::Path;

/// Block size for batch processing (10,000 records)
///
/// # Evidence
///
/// Entry 027 (1,440 measurements):
/// - Record-by-record: 82-86% NEON overhead
/// - Block size 10K: 4-8% overhead
/// - Sweet spot: ~1.5 MB for 150bp reads
pub const BLOCK_SIZE: usize = 10_000;

/// FASTQ streaming parser with block-based processing
///
/// # Memory Footprint
///
/// - Block buffer: ~1.5 MB (10K × 150bp)
/// - Line buffers: ~512 bytes
/// - **Total: ~5 MB constant** (Entry 026 validated)
///
/// # Example
///
/// ```no_run
/// use biometal::FastqStream;
/// use biometal::io::DataSource;
///
/// # fn main() -> biometal::Result<()> {
/// let source = DataSource::from_path("large.fq.gz");
/// let stream = FastqStream::new(source)?;
///
/// for record in stream {
///     let record = record?;
///     // Process one record at a time (constant memory)
/// }
/// # Ok(())
/// # }
/// ```
pub struct FastqStream<R: BufRead> {
    reader: R,
    block_buffer: Vec<FastqRecord>,
    block_position: usize,
    line1: String,
    line2: String,
    line3: String,
    line4: String,
    line_number: usize,
    finished: bool,
}

impl<R: BufRead> FastqStream<R> {
    /// Create a new FASTQ stream from a buffered reader
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader,
            block_buffer: Vec::with_capacity(BLOCK_SIZE),
            block_position: 0,
            line1: String::with_capacity(256),
            line2: String::with_capacity(256),
            line3: String::with_capacity(256),
            line4: String::with_capacity(256),
            line_number: 0,
            finished: false,
        }
    }

}

impl FastqStream<CompressedReader> {
    /// Create a FASTQ stream from a data source (with compression support)
    ///
    /// # Optimization Stack
    ///
    /// 1. DataSource abstraction (Rule 6)
    /// 2. Threshold-based mmap (Rule 4, if ≥50 MB)
    /// 3. Parallel bgzip decompression (Rule 3, 6.5× speedup)
    /// 4. Block-based streaming (Rule 2+5, preserves NEON gains)
    pub fn new(source: DataSource) -> Result<Self> {
        let compressed_reader = CompressedReader::new(source)?;
        Ok(Self::from_reader(compressed_reader))
    }

    /// Create a FASTQ stream from a file path
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let source = DataSource::from_path(path);
        Self::new(source)
    }
}

impl<R: BufRead> FastqStream<R> {

    /// Read one FASTQ record from the reader
    fn read_record(&mut self) -> Result<Option<FastqRecord>> {
        // Clear reusable buffers (Rule 5: buffer reuse)
        self.line1.clear();
        self.line2.clear();
        self.line3.clear();
        self.line4.clear();

        // Read 4 lines
        let n1 = self.reader.read_line(&mut self.line1)?;
        if n1 == 0 {
            return Ok(None);
        }
        self.line_number += 1;

        let n2 = self.reader.read_line(&mut self.line2)?;
        if n2 == 0 {
            return Err(BiometalError::InvalidFastqFormat {
                line: self.line_number,
                msg: "Unexpected end of file after header".to_string(),
            });
        }
        self.line_number += 1;

        let n3 = self.reader.read_line(&mut self.line3)?;
        if n3 == 0 {
            return Err(BiometalError::InvalidFastqFormat {
                line: self.line_number,
                msg: "Unexpected end of file after sequence".to_string(),
            });
        }
        self.line_number += 1;

        let n4 = self.reader.read_line(&mut self.line4)?;
        if n4 == 0 {
            return Err(BiometalError::InvalidFastqFormat {
                line: self.line_number,
                msg: "Unexpected end of file after separator".to_string(),
            });
        }
        self.line_number += 1;

        // Validate format
        if !self.line1.starts_with('@') {
            return Err(BiometalError::InvalidFastqFormat {
                line: self.line_number - 3,
                msg: format!("Expected '@' at start of header, got: {}", &self.line1[..1.min(self.line1.len())]),
            });
        }

        if !self.line3.starts_with('+') {
            return Err(BiometalError::InvalidFastqFormat {
                line: self.line_number - 1,
                msg: format!("Expected '+' at start of separator, got: {}", &self.line3[..1.min(self.line3.len())]),
            });
        }

        // Extract components
        let id = self.line1[1..].trim_end().to_string();
        let sequence = self.line2.trim_end().as_bytes().to_vec();
        let quality = self.line4.trim_end().as_bytes().to_vec();

        // Validate lengths match
        if sequence.len() != quality.len() {
            return Err(BiometalError::InvalidFastqFormat {
                line: self.line_number,
                msg: format!(
                    "Sequence length ({}) != quality length ({})",
                    sequence.len(),
                    quality.len()
                ),
            });
        }

        Ok(Some(FastqRecord {
            id,
            sequence,
            quality,
        }))
    }

    /// Fill the block buffer with up to BLOCK_SIZE records
    fn fill_block(&mut self) -> Result<usize> {
        self.block_buffer.clear();
        self.block_position = 0;

        while self.block_buffer.len() < BLOCK_SIZE {
            match self.read_record()? {
                Some(record) => self.block_buffer.push(record),
                None => {
                    self.finished = true;
                    break;
                }
            }
        }

        Ok(self.block_buffer.len())
    }
}

impl<R: BufRead> Iterator for FastqStream<R> {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.block_position >= self.block_buffer.len() {
            if self.finished {
                return None;
            }

            match self.fill_block() {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => return Some(Err(e)),
            }
        }

        if self.block_position < self.block_buffer.len() {
            let record = self.block_buffer[self.block_position].clone();
            self.block_position += 1;
            Some(Ok(record))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    #[test]
    fn test_block_size_constant() {
        assert_eq!(BLOCK_SIZE, 10_000);
    }

    #[test]
    fn test_parse_valid_fastq() {
        let data = b"@SEQ_ID\nGATTACA\n+\n!!!!!!!\n";
        let cursor = Cursor::new(data);
        let mut stream = FastqStream::from_reader(BufReader::new(cursor));

        let record = stream.next().unwrap().unwrap();
        assert_eq!(record.id, "SEQ_ID");
        assert_eq!(record.sequence, b"GATTACA");
        assert_eq!(record.quality, b"!!!!!!!");
    }

    #[test]
    fn test_parse_multiple_records() {
        let data = b"@SEQ1\nGAT\n+\n!!!\n@SEQ2\nTACA\n+\n!!!!\n";
        let cursor = Cursor::new(data);
        let stream = FastqStream::from_reader(BufReader::new(cursor));

        let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "SEQ1");
        assert_eq!(records[1].id, "SEQ2");
    }

    #[test]
    fn test_invalid_header() {
        let data = b"SEQ_ID\nGATTACA\n+\n!!!!!!!!\n";
        let cursor = Cursor::new(data);
        let mut stream = FastqStream::from_reader(BufReader::new(cursor));

        let result = stream.next().unwrap();
        assert!(matches!(result, Err(BiometalError::InvalidFastqFormat { .. })));
    }

    // Property-based tests
    use proptest::prelude::*;

    proptest! {
        /// Test that valid FASTQ records can be parsed correctly
        #[test]
        fn test_fastq_roundtrip(
            id in "[A-Za-z0-9_]{1,50}",
            seq in "[ACGTN]{1,500}",
        ) {
            let qual = "I".repeat(seq.len());
            let fastq = format!("@{}\n{}\n+\n{}\n", id, seq, qual);

            let stream = FastqStream::from_reader(BufReader::new(Cursor::new(fastq.as_bytes())));
            let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

            prop_assert_eq!(records.len(), 1);
            prop_assert_eq!(&records[0].id, &id);
            prop_assert_eq!(&records[0].sequence, seq.as_bytes());
            prop_assert_eq!(&records[0].quality, qual.as_bytes());
        }

        /// Test that FASTQ rejects mismatched sequence and quality lengths
        #[test]
        fn test_fastq_rejects_length_mismatch(
            id in "[A-Za-z0-9_]{1,50}",
            seq in "[ACGT]{10,20}",
            qual_len in 21..30usize, // Different length than seq
        ) {
            let qual = "I".repeat(qual_len);
            let fastq = format!("@{}\n{}\n+\n{}\n", id, seq, qual);

            let stream = FastqStream::from_reader(BufReader::new(Cursor::new(fastq.as_bytes())));
            let result: Result<Vec<_>> = stream.collect();

            // Should detect mismatch
            prop_assert!(result.is_err());
        }

        /// Test that multiple valid records can be parsed
        #[test]
        fn test_fastq_multiple_records(
            records_count in 1..10usize,
        ) {
            let mut fastq = String::new();
            for i in 0..records_count {
                let seq = "ACGT".repeat(10);
                let qual = "I".repeat(40);
                fastq.push_str(&format!("@read_{}\n{}\n+\n{}\n", i, seq, qual));
            }

            let stream = FastqStream::from_reader(BufReader::new(Cursor::new(fastq.as_bytes())));
            let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

            prop_assert_eq!(records.len(), records_count);
            for (i, record) in records.iter().enumerate() {
                prop_assert_eq!(&record.id, &format!("read_{}", i));
            }
        }

        /// Test that invalid headers are rejected
        #[test]
        fn test_fastq_invalid_header(
            id in "[A-Za-z0-9_]{1,50}",
            seq in "[ACGT]{10,20}",
        ) {
            let qual = "I".repeat(seq.len());
            // Missing '@' prefix
            let fastq = format!("{}\n{}\n+\n{}\n", id, seq, qual);

            let stream = FastqStream::from_reader(BufReader::new(Cursor::new(fastq.as_bytes())));
            let result: Result<Vec<_>> = stream.collect();

            prop_assert!(result.is_err());
        }
    }
}
