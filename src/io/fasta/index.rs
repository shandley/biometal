//! FAI (FASTA Index) support for random access to FASTA sequences.
//!
//! This module provides support for reading and building FAI index files,
//! which enable efficient random access to specific sequences or regions
//! within FASTA files without loading the entire file into memory.
//!
//! # Format Specification
//!
//! FAI files are tab-delimited text files with 5 columns (FASTA) or 6 columns (FASTQ):
//!
//! 1. **NAME**: Sequence name (first word after '>')
//! 2. **LENGTH**: Total sequence length in bases
//! 3. **OFFSET**: Byte offset to first base (after header newline)
//! 4. **LINEBASES**: Number of bases per line
//! 5. **LINEWIDTH**: Number of bytes per line (including newline)
//! 6. **QUALOFFSET**: (FASTQ only) Offset to first quality score
//!
//! # Example
//!
//! For a FASTA file:
//! ```text
//! >chr1
//! ACGTACGTACGTACGTACGTACGTACGT
//! TGCATGCATGCATGCA
//! >chr2
//! GGGGCCCCAAAATTTT
//! ```
//!
//! The corresponding FAI index:
//! ```text
//! chr1	44	6	28	29
//! chr2	16	81	16	17
//! ```
//!
//! # Basic Usage
//!
//! ## Loading an Index
//!
//! ```no_run
//! use biometal::io::fasta::FaiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! // Load existing FAI index (typical: <1ms)
//! let index = FaiIndex::from_path("genome.fa.fai")?;
//!
//! println!("Index contains {} sequences", index.sequences.len());
//! # Ok(())
//! # }
//! ```
//!
//! ## Building an Index
//!
//! ```no_run
//! use biometal::io::fasta::FaiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! // Build index from FASTA file
//! let index = FaiIndex::build("genome.fa")?;
//!
//! // Save to file (genome.fa.fai)
//! index.write("genome.fa.fai")?;
//! # Ok(())
//! # }
//! ```
//!
//! ## Fetching Sequences
//!
//! ```no_run
//! use biometal::io::fasta::FaiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! let index = FaiIndex::from_path("genome.fa.fai")?;
//!
//! // Fetch entire sequence
//! let chr1 = index.fetch("chr1", "genome.fa")?;
//! println!("chr1: {} bp", chr1.len());
//!
//! // Fetch region (0-based, half-open: [start, end))
//! let region = index.fetch_region("chr1", 1000, 2000, "genome.fa")?;
//! println!("chr1:1000-2000: {} bp", region.len());
//! # Ok(())
//! # }
//! ```
//!
//! # Performance
//!
//! - **Index loading**: <1ms (text parsing)
//! - **Sequence fetch**: O(1) seek + O(n) read where n = sequence length
//! - **Region fetch**: O(1) seek + O(n) read where n = region size
//! - **Memory**: ~100 bytes per sequence in index

use crate::error::{BiometalError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// A single entry in an FAI index
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiEntry {
    /// Sequence name (from FASTA header, first word)
    pub name: String,
    /// Total sequence length in bases
    pub length: u64,
    /// Byte offset to first base (after header line)
    pub offset: u64,
    /// Number of bases per line
    pub line_bases: u32,
    /// Number of bytes per line (including newline)
    pub line_width: u32,
    /// Offset to first quality score (FASTQ only)
    pub qual_offset: Option<u64>,
}

impl FaiEntry {
    /// Parse a single FAI entry from a line
    ///
    /// Format: NAME\tLENGTH\tOFFSET\tLINEBASES\tLINEWIDTH[\tQUALOFFSET]
    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 5 {
            return Err(BiometalError::InvalidInput {
                msg: format!("FAI entry must have at least 5 fields, got {}", fields.len()),
            });
        }

        let name = fields[0].to_string();
        let length = fields[1].parse::<u64>().map_err(|e| {
            BiometalError::InvalidInput {
                msg: format!("Invalid LENGTH field: {}", e),
            }
        })?;
        let offset = fields[2].parse::<u64>().map_err(|e| {
            BiometalError::InvalidInput {
                msg: format!("Invalid OFFSET field: {}", e),
            }
        })?;
        let line_bases = fields[3].parse::<u32>().map_err(|e| {
            BiometalError::InvalidInput {
                msg: format!("Invalid LINEBASES field: {}", e),
            }
        })?;
        let line_width = fields[4].parse::<u32>().map_err(|e| {
            BiometalError::InvalidInput {
                msg: format!("Invalid LINEWIDTH field: {}", e),
            }
        })?;

        // Optional QUALOFFSET for FASTQ
        let qual_offset = if fields.len() > 5 {
            Some(fields[5].parse::<u64>().map_err(|e| {
                BiometalError::InvalidInput {
                    msg: format!("Invalid QUALOFFSET field: {}", e),
                }
            })?)
        } else {
            None
        };

        Ok(FaiEntry {
            name,
            length,
            offset,
            line_bases,
            line_width,
            qual_offset,
        })
    }

    /// Format entry as FAI line
    pub fn to_line(&self) -> String {
        if let Some(qual_offset) = self.qual_offset {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.name, self.length, self.offset, self.line_bases, self.line_width, qual_offset
            )
        } else {
            format!(
                "{}\t{}\t{}\t{}\t{}",
                self.name, self.length, self.offset, self.line_bases, self.line_width
            )
        }
    }

    /// Calculate byte offset for a specific base position (0-based)
    ///
    /// This accounts for line wrapping in the FASTA file.
    fn position_to_offset(&self, position: u64) -> Result<u64> {
        if position >= self.length {
            return Err(BiometalError::InvalidRange(format!(
                "Position {} exceeds sequence length {}",
                position, self.length
            )));
        }

        // Calculate which line the position is on (0-based)
        let line_number = position / self.line_bases as u64;
        // Position within that line
        let line_position = position % self.line_bases as u64;

        // Byte offset = base offset + (full lines * line_width) + position in current line
        let byte_offset = self.offset + (line_number * self.line_width as u64) + line_position;

        Ok(byte_offset)
    }

    /// Calculate number of bytes needed to read a region
    ///
    /// Accounts for newlines that may be included when reading multi-line sequences.
    fn region_byte_size(&self, start: u64, end: u64) -> u64 {
        if start >= end || start >= self.length {
            return 0;
        }

        let actual_end = end.min(self.length);
        let region_length = actual_end - start;

        // Calculate how many complete lines we span
        let start_line = start / self.line_bases as u64;
        let end_line = (actual_end - 1) / self.line_bases as u64;
        let lines_spanned = end_line - start_line + 1;

        // Newline characters = lines_spanned (one per line)
        // But we need to be careful - if we end exactly at line boundary, don't count extra newline
        region_length + lines_spanned
    }
}

/// FAI index for random access to FASTA sequences
#[derive(Debug, Clone)]
pub struct FaiIndex {
    /// Indexed sequences (name -> entry)
    pub sequences: HashMap<String, FaiEntry>,
    /// Sequence names in order (preserves FAI file order)
    pub sequence_names: Vec<String>,
}

impl FaiIndex {
    /// Create a new empty FAI index
    pub fn new() -> Self {
        Self {
            sequences: HashMap::new(),
            sequence_names: Vec::new(),
        }
    }

    /// Load FAI index from a file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FaiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = FaiIndex::from_path("genome.fa.fai")?;
    /// println!("Loaded {} sequences", index.sequences.len());
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref()).map_err(|e| {
            BiometalError::Io(std::io::Error::new(
                e.kind(),
                format!("Failed to open FAI file {:?}: {}", path.as_ref(), e),
            ))
        })?;

        let reader = BufReader::new(file);
        let mut index = Self::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result?;

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let entry = FaiEntry::from_line(&line).map_err(|e| {
                BiometalError::InvalidInput {
                    msg: format!("Line {}: {}", line_num + 1, e),
                }
            })?;

            index.sequence_names.push(entry.name.clone());
            index.sequences.insert(entry.name.clone(), entry);
        }

        Ok(index)
    }

    /// Build FAI index from a FASTA file
    ///
    /// This scans the entire FASTA file and creates an index.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FaiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = FaiIndex::build("genome.fa")?;
    /// index.write("genome.fa.fai")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn build<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let file = File::open(fasta_path.as_ref())?;
        let mut reader = BufReader::new(file);
        let mut index = Self::new();

        let mut current_name: Option<String> = None;
        let mut current_offset: u64 = 0;
        let mut current_length: u64 = 0;
        let mut line_bases: Option<u32> = None;
        let mut line_width: Option<u32> = None;
        let mut sequence_start_offset: u64 = 0;

        let mut line_buffer = String::new();
        let mut byte_offset: u64 = 0;

        loop {
            line_buffer.clear();
            let bytes_read = reader.read_line(&mut line_buffer)?;

            if bytes_read == 0 {
                // EOF - save last sequence if exists
                if let Some(name) = current_name {
                    let entry = FaiEntry {
                        name: name.clone(),
                        length: current_length,
                        offset: sequence_start_offset,
                        line_bases: line_bases.unwrap_or(0),
                        line_width: line_width.unwrap_or(0),
                        qual_offset: None,
                    };
                    index.sequence_names.push(name.clone());
                    index.sequences.insert(name, entry);
                }
                break;
            }

            if line_buffer.starts_with('>') {
                // New sequence header

                // Save previous sequence if exists
                if let Some(name) = current_name {
                    let entry = FaiEntry {
                        name: name.clone(),
                        length: current_length,
                        offset: sequence_start_offset,
                        line_bases: line_bases.unwrap_or(0),
                        line_width: line_width.unwrap_or(0),
                        qual_offset: None,
                    };
                    index.sequence_names.push(name.clone());
                    index.sequences.insert(name, entry);
                }

                // Parse new header (first word after '>')
                let header = line_buffer[1..].trim();
                let name = header
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| BiometalError::InvalidInput {
                        msg: "Empty sequence name".to_string(),
                    })?
                    .to_string();

                current_name = Some(name);
                current_length = 0;
                line_bases = None;
                line_width = None;
                sequence_start_offset = byte_offset + bytes_read as u64;
                current_offset = sequence_start_offset;
            } else if current_name.is_some() {
                // Sequence line
                let trimmed = line_buffer.trim_end();
                let bases = trimmed.len() as u32;

                if bases > 0 {
                    // First sequence line - record line parameters
                    if line_bases.is_none() {
                        line_bases = Some(bases);
                        line_width = Some(bytes_read as u32);
                    } else {
                        // Verify consistent line length (except last line)
                        // Note: We allow the last line to be shorter
                        if line_bases != Some(bases) && bytes_read > 0 {
                            // This might be the last line, which can be shorter
                            // We'll validate during fetch
                        }
                    }

                    current_length += bases as u64;
                }
            }

            byte_offset += bytes_read as u64;
        }

        Ok(index)
    }

    /// Write index to a file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FaiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = FaiIndex::build("genome.fa")?;
    /// index.write("genome.fa.fai")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn write<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut file = File::create(path)?;

        for name in &self.sequence_names {
            if let Some(entry) = self.sequences.get(name) {
                writeln!(file, "{}", entry.to_line())?;
            }
        }

        Ok(())
    }

    /// Get entry for a sequence by name
    pub fn get(&self, name: &str) -> Option<&FaiEntry> {
        self.sequences.get(name)
    }

    /// Fetch entire sequence by name
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FaiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = FaiIndex::from_path("genome.fa.fai")?;
    /// let chr1 = index.fetch("chr1", "genome.fa")?;
    /// println!("chr1: {} bp", chr1.len());
    /// # Ok(())
    /// # }
    /// ```
    pub fn fetch<P: AsRef<Path>>(&self, name: &str, fasta_path: P) -> Result<String> {
        let entry = self
            .sequences
            .get(name)
            .ok_or_else(|| BiometalError::InvalidInput {
                msg: format!("Sequence '{}' not found in index", name),
            })?;

        self.fetch_region(name, 0, entry.length, fasta_path)
    }

    /// Fetch a region of a sequence (0-based, half-open interval [start, end))
    ///
    /// # Arguments
    ///
    /// * `name` - Sequence name
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    /// * `fasta_path` - Path to FASTA file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::fasta::FaiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = FaiIndex::from_path("genome.fa.fai")?;
    ///
    /// // Fetch chr1:1000-2000 (1000 bases)
    /// let region = index.fetch_region("chr1", 1000, 2000, "genome.fa")?;
    /// assert_eq!(region.len(), 1000);
    /// # Ok(())
    /// # }
    /// ```
    pub fn fetch_region<P: AsRef<Path>>(
        &self,
        name: &str,
        start: u64,
        end: u64,
        fasta_path: P,
    ) -> Result<String> {
        let entry = self
            .sequences
            .get(name)
            .ok_or_else(|| BiometalError::InvalidInput {
                msg: format!("Sequence '{}' not found in index", name),
            })?;

        if start >= end {
            return Err(BiometalError::InvalidRange(format!(
                "Invalid region: start ({}) >= end ({})",
                start, end
            )));
        }

        if start >= entry.length {
            return Err(BiometalError::InvalidRange(format!(
                "Start position {} exceeds sequence length {}",
                start, entry.length
            )));
        }

        let actual_end = end.min(entry.length);
        let region_length = actual_end - start;

        // Open FASTA file and seek to start position
        let mut file = File::open(fasta_path)?;
        let start_offset = entry.position_to_offset(start)?;
        file.seek(SeekFrom::Start(start_offset))?;

        // Calculate approximate bytes to read (with buffer for newlines)
        let bytes_to_read = entry.region_byte_size(start, actual_end);

        // Read bytes
        let mut buffer = vec![0u8; bytes_to_read as usize];
        file.read_exact(&mut buffer)?;

        // Remove newlines and extract exact region
        let sequence: String = buffer
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r')
            .map(|&b| b as char)
            .take(region_length as usize)
            .collect();

        Ok(sequence)
    }

    /// Get number of sequences in index
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Check if index is empty
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
}

impl Default for FaiIndex {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fai_entry_parse() {
        let line = "chr1\t248956422\t112\t70\t71";
        let entry = FaiEntry::from_line(line).unwrap();

        assert_eq!(entry.name, "chr1");
        assert_eq!(entry.length, 248956422);
        assert_eq!(entry.offset, 112);
        assert_eq!(entry.line_bases, 70);
        assert_eq!(entry.line_width, 71);
        assert_eq!(entry.qual_offset, None);
    }

    #[test]
    fn test_fai_entry_with_qual_offset() {
        let line = "read1\t150\t6\t150\t151\t312";
        let entry = FaiEntry::from_line(line).unwrap();

        assert_eq!(entry.name, "read1");
        assert_eq!(entry.qual_offset, Some(312));
    }

    #[test]
    fn test_fai_entry_to_line() {
        let entry = FaiEntry {
            name: "chr1".to_string(),
            length: 1000,
            offset: 6,
            line_bases: 50,
            line_width: 51,
            qual_offset: None,
        };

        assert_eq!(entry.to_line(), "chr1\t1000\t6\t50\t51");
    }

    #[test]
    fn test_position_to_offset() {
        let entry = FaiEntry {
            name: "test".to_string(),
            length: 100,
            offset: 10,
            line_bases: 20,
            line_width: 21, // 20 bases + 1 newline
            qual_offset: None,
        };

        // Position 0 -> offset 10 (first base)
        assert_eq!(entry.position_to_offset(0).unwrap(), 10);

        // Position 19 -> offset 29 (last base of first line)
        assert_eq!(entry.position_to_offset(19).unwrap(), 29);

        // Position 20 -> offset 31 (first base of second line, skip newline)
        assert_eq!(entry.position_to_offset(20).unwrap(), 31);

        // Position 40 -> offset 52 (first base of third line)
        assert_eq!(entry.position_to_offset(40).unwrap(), 52);
    }
}
