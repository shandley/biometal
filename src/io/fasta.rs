//! FASTA streaming parser with constant memory (Rule 5)
//!
//! # Format
//!
//! FASTA format consists of:
//! - Header line starting with '>' followed by sequence identifier
//! - One or more sequence lines (can be wrapped)
//!
//! Example:
//! ```text
//! >sequence1 description
//! GATTACAGATTACA
//! TGCATGCA
//! >sequence2
//! ACGTACGT
//! ```
//!
//! # Architecture
//!
//! This parser follows the same streaming principles as FastqStream:
//! - Constant memory: ~5 MB regardless of file size (Rule 5)
//! - Iterator-based: Process one record at a time
//! - Block-based buffering: Works with NEON operations (Rule 2)
//! - Compression support: Automatic bgzip decompression (Rule 3)

use crate::error::{BiometalError, Result};
use crate::io::compression::{CompressedReader, DataSource};
use crate::types::FastaRecord;
use std::io::BufRead;
use std::path::Path;

/// FASTA streaming parser with constant memory footprint
///
/// # Architecture
///
/// Unlike FASTQ (fixed 4 lines per record), FASTA records can span multiple
/// lines. This parser handles multi-line sequences while maintaining constant
/// memory usage.
///
/// # Memory Footprint (Rule 5)
///
/// - Line buffers: ~1 KB (reused across records)
/// - Current record: ~200 bytes average
/// - Total: ~5 MB including decompression buffer
///
/// # Example
///
/// ```no_run
/// use biometal::FastaStream;
///
/// let stream = FastaStream::from_path("genome.fa.gz")?;
/// for record in stream {
///     let record = record?;
///     println!("{}: {} bp", record.id, record.sequence.len());
/// }
/// # Ok::<(), biometal::error::BiometalError>(())
/// ```
pub struct FastaStream<R: BufRead> {
    reader: R,
    line_buffer: String,
    line_number: usize,
    finished: bool,
    /// Peek buffer for look-ahead (to detect next record start)
    next_line: Option<String>,
}

impl FastaStream<CompressedReader> {
    /// Create a FASTA stream from a data source
    ///
    /// Supports:
    /// - Local files (compressed or uncompressed)
    /// - Future: HTTP/HTTPS URLs (Week 3-4)
    /// - Future: SRA accessions (Week 3-4)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::FastaStream;
    /// use biometal::io::compression::DataSource;
    ///
    /// let source = DataSource::from_path("genome.fa.gz");
    /// let stream = FastaStream::new(source)?;
    /// # Ok::<(), biometal::error::BiometalError>(())
    /// ```
    pub fn new(source: DataSource) -> Result<Self> {
        let compressed_reader = CompressedReader::new(source)?;
        Ok(Self::from_reader(compressed_reader))
    }

    /// Create a FASTA stream from a local file path
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::FastaStream;
    ///
    /// let stream = FastaStream::from_path("genome.fa.gz")?;
    /// # Ok::<(), biometal::error::BiometalError>(())
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let source = DataSource::from_path(path);
        Self::new(source)
    }
}

impl<R: BufRead> FastaStream<R> {
    /// Create a FASTA stream from any buffered reader
    ///
    /// This is useful for testing or reading from in-memory sources.
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader,
            line_buffer: String::with_capacity(256),
            line_number: 0,
            finished: false,
            next_line: None,
        }
    }

    /// Read a single FASTA record
    fn read_record(&mut self) -> Result<Option<FastaRecord>> {
        if self.finished {
            return Ok(None);
        }

        // Get header line (either from peek buffer or read new)
        let header = if let Some(peeked) = self.next_line.take() {
            peeked
        } else {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => {
                    self.finished = true;
                    return Ok(None);
                }
                Ok(_) => {
                    self.line_number += 1;
                    self.line_buffer.trim_end().to_string()
                }
                Err(e) => return Err(BiometalError::Io(e)),
            }
        };

        // Check for empty file or leading whitespace
        if header.is_empty() {
            return self.read_record(); // Skip empty lines
        }

        // Validate header starts with '>'
        if !header.starts_with('>') {
            return Err(BiometalError::InvalidFastaFormat {
                line: self.line_number,
                msg: format!("Expected '>' at start of header, got: {}", header),
            });
        }

        // Extract ID (everything after '>', up to first whitespace)
        let id = header[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        // Read sequence lines until next header or EOF
        let mut sequence = Vec::new();

        loop {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => {
                    // EOF reached
                    self.finished = true;
                    break;
                }
                Ok(_) => {
                    self.line_number += 1;
                    let line = self.line_buffer.trim();

                    if line.is_empty() {
                        continue; // Skip empty lines
                    }

                    if line.starts_with('>') {
                        // Start of next record - save for next iteration
                        self.next_line = Some(line.to_string());
                        break;
                    }

                    // Sequence line - append to current record
                    sequence.extend_from_slice(line.as_bytes());
                }
                Err(e) => return Err(BiometalError::Io(e)),
            }
        }

        if sequence.is_empty() {
            return Err(BiometalError::InvalidFastaFormat {
                line: self.line_number,
                msg: "Record has no sequence".to_string(),
            });
        }

        Ok(Some(FastaRecord::new(id, sequence)))
    }
}

impl<R: BufRead> Iterator for FastaStream<R> {
    type Item = Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    #[test]
    fn test_parse_single_record() {
        let fasta = b">seq1\nGATTACA\n";
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let mut stream = FastaStream::from_reader(reader);

        let record = stream.next().unwrap().unwrap();
        assert_eq!(record.id, "seq1");
        assert_eq!(record.sequence, b"GATTACA");

        assert!(stream.next().is_none());
    }

    #[test]
    fn test_parse_multiple_records() {
        let fasta = b">seq1\nGATTACA\n>seq2\nACGT\n";
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let stream = FastaStream::from_reader(reader);

        let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, b"GATTACA");

        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, b"ACGT");
    }

    #[test]
    fn test_parse_multiline_sequence() {
        let fasta = b">seq1\nGATT\nACA\n>seq2\nACGT\n";
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let stream = FastaStream::from_reader(reader);

        let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, b"GATTACA"); // Multi-line concatenated

        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, b"ACGT");
    }

    #[test]
    fn test_parse_with_description() {
        let fasta = b">seq1 this is a description\nGATTACA\n";
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let mut stream = FastaStream::from_reader(reader);

        let record = stream.next().unwrap().unwrap();
        assert_eq!(record.id, "seq1"); // ID only, description stripped
        assert_eq!(record.sequence, b"GATTACA");
    }

    #[test]
    fn test_parse_with_empty_lines() {
        let fasta = b">seq1\n\nGATTACA\n\n>seq2\nACGT\n\n";
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let stream = FastaStream::from_reader(reader);

        let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, b"GATTACA");
    }

    #[test]
    fn test_invalid_no_header() {
        let fasta = b"GATTACA\n"; // No header
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let mut stream = FastaStream::from_reader(reader);

        let result = stream.next().unwrap();
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            BiometalError::InvalidFastaFormat { .. }
        ));
    }

    #[test]
    fn test_empty_sequence() {
        let fasta = b">seq1\n>seq2\nACGT\n"; // seq1 has no sequence
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let mut stream = FastaStream::from_reader(reader);

        let result = stream.next().unwrap();
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            BiometalError::InvalidFastaFormat { .. }
        ));
    }

    #[test]
    fn test_empty_file() {
        let fasta = b"";
        let cursor = Cursor::new(fasta);
        let reader = BufReader::new(cursor);
        let mut stream = FastaStream::from_reader(reader);

        assert!(stream.next().is_none());
    }

    // Property-based tests
    use proptest::prelude::*;

    proptest! {
        /// Test that valid FASTA records can be parsed correctly
        #[test]
        fn test_fasta_roundtrip(
            id in "[A-Za-z0-9_]{1,50}",
            seq in "[ACGTN]{1,500}",
        ) {
            let fasta = format!(">{}\n{}\n", id, seq);

            let stream = FastaStream::from_reader(BufReader::new(Cursor::new(fasta.as_bytes())));
            let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

            prop_assert_eq!(records.len(), 1);
            prop_assert_eq!(&records[0].id, &id);
            prop_assert_eq!(&records[0].sequence, seq.as_bytes());
        }

        /// Test that multi-line FASTA sequences are correctly joined
        #[test]
        fn test_fasta_multiline(
            id in "[A-Za-z0-9_]{1,50}",
            line_count in 2..10usize,
        ) {
            let mut fasta = format!(">{}\n", id);
            let line_seq = "ACGT".repeat(20); // 80 bp per line
            let full_seq = line_seq.repeat(line_count);

            for _ in 0..line_count {
                fasta.push_str(&line_seq);
                fasta.push('\n');
            }

            let stream = FastaStream::from_reader(BufReader::new(Cursor::new(fasta.as_bytes())));
            let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

            prop_assert_eq!(records.len(), 1);
            prop_assert_eq!(&records[0].sequence, full_seq.as_bytes());
        }

        /// Test that multiple FASTA records can be parsed
        #[test]
        fn test_fasta_multiple_records(
            records_count in 1..10usize,
        ) {
            let mut fasta = String::new();
            for i in 0..records_count {
                let seq = "ACGT".repeat(25);
                fasta.push_str(&format!(">seq_{}\n{}\n", i, seq));
            }

            let stream = FastaStream::from_reader(BufReader::new(Cursor::new(fasta.as_bytes())));
            let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

            prop_assert_eq!(records.len(), records_count);
            for (i, record) in records.iter().enumerate() {
                prop_assert_eq!(&record.id, &format!("seq_{}", i));
            }
        }

        /// Test that FASTA descriptions are stripped (only ID kept)
        #[test]
        fn test_fasta_description_stripped(
            id in "[A-Za-z0-9_]{1,50}",
            description in "[A-Za-z0-9 ]{1,100}",
            seq in "[ACGT]{10,100}",
        ) {
            let fasta = format!(">{} {}\n{}\n", id, description, seq);

            let stream = FastaStream::from_reader(BufReader::new(Cursor::new(fasta.as_bytes())));
            let records: Vec<_> = stream.collect::<Result<Vec<_>>>().unwrap();

            prop_assert_eq!(records.len(), 1);
            // Only the ID (before first space) should be kept
            prop_assert_eq!(&records[0].id, &id);
        }

        /// Test that invalid headers (missing '>') are rejected
        #[test]
        fn test_fasta_invalid_header(
            id in "[A-Za-z0-9_]{1,50}",
            seq in "[ACGT]{10,100}",
        ) {
            // Missing '>' prefix
            let fasta = format!("{}\n{}\n", id, seq);

            let stream = FastaStream::from_reader(BufReader::new(Cursor::new(fasta.as_bytes())));
            let result: Result<Vec<_>> = stream.collect();

            prop_assert!(result.is_err());
        }
    }
}
