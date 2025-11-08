//! BAM streaming reader.
//!
//! Provides a streaming interface for reading BAM files with constant memory.
//! Follows biometal's Rule 5 (streaming architecture, ~5 MB memory).
//!
//! # Design
//!
//! - Iterator-based interface (constant memory)
//! - Header read separately from records
//! - Records are not accumulated (streaming)
//! - Compatible with BGZF decompression (Phase 2 will add parallel BGZF)
//!
//! # Usage
//!
//! ```no_run
//! use biometal::io::bam::BamReader;
//! use std::fs::File;
//! use std::io::BufReader;
//!
//! # fn main() -> std::io::Result<()> {
//! let file = File::open("alignments.bam")?;
//! let reader = BufReader::new(file);
//! let mut bam = BamReader::new(reader)?;
//!
//! println!("Header: {} references", bam.header().reference_count());
//!
//! for result in bam.records() {
//!     let record = result?;
//!     println!("{} at {}", record.name, record.position.unwrap_or(-1));
//! }
//! # Ok(())
//! # }
//! ```

use super::header::{read_header, Header};
use super::record::{parse_record, Record};
use std::io::{self, BufRead};
use std::path::Path;
use std::fs::File;

/// BAM file reader with streaming interface.
///
/// Reads BAM files with constant memory (Rule 5: streaming architecture).
/// The header is read once during construction, then records are streamed.
///
/// # Buffer Reuse
///
/// Maintains an internal buffer that's reused across record reads to avoid
/// repeated allocations. This buffer grows to accommodate the largest record
/// seen, then stays at that size for subsequent reads.
pub struct BamReader<R> {
    /// Underlying reader (typically BufReader<File>)
    reader: R,
    /// BAM header (read during construction)
    header: Header,
    /// Reusable buffer for reading record data (avoids repeated allocations)
    buffer: Vec<u8>,
}

impl<R: BufRead> BamReader<R> {
    /// Create a new BAM reader.
    ///
    /// Reads and validates the BAM header immediately.
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Cannot read header
    /// - Invalid magic bytes
    /// - Header is malformed
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("alignments.bam")?;
    /// let reader = BufReader::new(file);
    /// let bam = BamReader::new(reader)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(mut reader: R) -> io::Result<Self> {
        let header = read_header(&mut reader)?;
        // Initialize buffer with reasonable capacity (typical record ~500 bytes)
        let buffer = Vec::with_capacity(512);
        Ok(Self {
            reader,
            header,
            buffer,
        })
    }

    /// Get a reference to the BAM header.
    ///
    /// The header contains SAM header text and reference sequence information.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Create an iterator over BAM records.
    ///
    /// Records are streamed with constant memory (not accumulated).
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("alignments.bam")?;
    /// let reader = BufReader::new(file);
    /// let mut bam = BamReader::new(reader)?;
    ///
    /// for result in bam.records() {
    ///     let record = result?;
    ///     // Process one record at a time (constant memory)
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records { reader: self }
    }

    /// Read a single record.
    ///
    /// Returns `Ok(None)` when EOF is reached.
    ///
    /// # Buffer Reuse
    ///
    /// This method reuses an internal buffer to avoid repeated allocations.
    /// The buffer grows to accommodate larger records but is not shrunk,
    /// making subsequent reads more efficient.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let file = File::open("alignments.bam")?;
    /// let reader = BufReader::new(file);
    /// let mut bam = BamReader::new(reader)?;
    ///
    /// while let Some(record) = bam.read_record()? {
    ///     println!("Read: {}", record.name);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn read_record(&mut self) -> io::Result<Option<Record>> {
        // Read block size (4 bytes, little-endian)
        let mut size_buf = [0u8; 4];
        match self.reader.read_exact(&mut size_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }

        let block_size = i32::from_le_bytes(size_buf);
        if block_size < 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid block size: {}", block_size),
            ));
        }

        let block_size = block_size as usize;

        // Reuse buffer, resizing if needed to fit block_size + 4 bytes for the size itself
        self.buffer.clear();
        self.buffer.reserve(block_size + 4);

        // Put block size at start (parse_record expects it)
        self.buffer.extend_from_slice(&size_buf);

        unsafe {
            // SAFETY: We're immediately reading into this memory
            let start_len = self.buffer.len();
            self.buffer.set_len(start_len + block_size);
        }

        // Read record data into buffer (after the block size)
        self.reader.read_exact(&mut self.buffer[4..])?;

        // Parse record from buffer (includes block size at start)
        let record = parse_record(&self.buffer)?;
        Ok(Some(record))
    }
}

impl BamReader<std::io::BufReader<File>> {
    /// Open a BAM file from a path.
    ///
    /// Convenience method that creates a BufReader internally.
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - File cannot be opened
    /// - Header is invalid
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::bam::BamReader;
    ///
    /// # fn main() -> std::io::Result<()> {
    /// let bam = BamReader::from_path("alignments.bam")?;
    /// println!("Opened BAM with {} references", bam.header().reference_count());
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::new(reader)
    }
}

/// Iterator over BAM records.
///
/// Created by [`BamReader::records()`]. Streams records with constant memory.
pub struct Records<'a, R> {
    reader: &'a mut BamReader<R>,
}

impl<'a, R: BufRead> Iterator for Records<'a, R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn create_minimal_bam() -> Vec<u8> {
        let mut data = Vec::new();

        // Magic
        data.extend_from_slice(b"BAM\x01");

        // SAM header text (empty)
        data.extend_from_slice(&0i32.to_le_bytes());

        // References (0)
        data.extend_from_slice(&0i32.to_le_bytes());

        // One minimal record
        let mut record_data = Vec::new();

        // Block size placeholder
        let block_size_pos = record_data.len();
        record_data.extend_from_slice(&0i32.to_le_bytes());

        // Reference ID (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // Position (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // l_read_name (5)
        record_data.push(5);

        // MAPQ (255)
        record_data.push(255);

        // BAI bin (0)
        record_data.extend_from_slice(&0u16.to_le_bytes());

        // n_cigar_op (0)
        record_data.extend_from_slice(&0u16.to_le_bytes());

        // FLAGS (4 = unmapped)
        record_data.extend_from_slice(&4u16.to_le_bytes());

        // l_seq (0)
        record_data.extend_from_slice(&0i32.to_le_bytes());

        // next_refID (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // next_pos (-1)
        record_data.extend_from_slice(&(-1i32).to_le_bytes());

        // tlen (0)
        record_data.extend_from_slice(&0i32.to_le_bytes());

        // read_name ("read\0")
        record_data.extend_from_slice(b"read\0");

        // Update block size
        let block_size = (record_data.len() - 4) as i32;
        record_data[block_size_pos..block_size_pos + 4]
            .copy_from_slice(&block_size.to_le_bytes());

        data.extend_from_slice(&record_data);

        data
    }

    #[test]
    fn test_bam_reader_new() {
        let bam_data = create_minimal_bam();
        let cursor = Cursor::new(bam_data);
        let bam = BamReader::new(cursor).unwrap();

        assert_eq!(bam.header().reference_count(), 0);
        assert_eq!(bam.header().text, "");
    }

    #[test]
    fn test_bam_reader_read_record() {
        let bam_data = create_minimal_bam();
        let cursor = Cursor::new(bam_data);
        let mut bam = BamReader::new(cursor).unwrap();

        // Read first record
        let record = bam.read_record().unwrap();
        assert!(record.is_some());
        let record = record.unwrap();
        assert_eq!(record.name, "read");

        // EOF
        let record = bam.read_record().unwrap();
        assert!(record.is_none());
    }

    #[test]
    fn test_bam_reader_records_iterator() {
        let bam_data = create_minimal_bam();
        let cursor = Cursor::new(bam_data);
        let mut bam = BamReader::new(cursor).unwrap();

        let records: Vec<_> = bam.records().collect::<io::Result<Vec<_>>>().unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "read");
    }

    #[test]
    fn test_invalid_magic() {
        let data = b"INVALID";
        let cursor = Cursor::new(data);
        let result = BamReader::new(cursor);
        assert!(result.is_err());
    }
}
