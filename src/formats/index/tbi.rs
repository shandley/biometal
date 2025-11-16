//! TBI (Tabix) index format support
//!
//! This module implements parsing and querying of Tabix index files (.tbi),
//! which enable fast random access to sorted, tab-delimited, BGZF-compressed
//! genomic files.
//!
//! # Format Specification
//!
//! TBI files are binary indexes with the following structure:
//!
//! ## Header
//! - Magic: "TBI\1" (4 bytes)
//! - n_ref: Number of reference sequences (int32)
//! - format: File format (int32: 0=generic, 1=SAM, 2=VCF, 3=BED3, 4=BED)
//! - col_seq: Column for sequence name (int32)
//! - col_beg: Column for start position (int32)
//! - col_end: Column for end position (int32, 0 if absent)
//! - meta: Comment character for header lines (int32)
//! - skip: Number of lines to skip (int32)
//! - l_nm: Length of concatenated sequence names (int32)
//! - names: Sequence names (null-terminated strings)
//!
//! ## Index Data (per reference)
//! - Binning index: Hierarchical bins with chunks
//! - Linear index: 16kb interval offsets
//!
//! # Binning Scheme
//!
//! Same as BAI: 37,450 bins covering 512 Mbp:
//! - Level 0: 1 bin (512 Mbp)
//! - Level 1: 8 bins (64 Mbp each)
//! - Level 2: 64 bins (8 Mbp each)
//! - Level 3: 512 bins (1 Mbp each)
//! - Level 4: 4096 bins (128 Kbp each)
//! - Level 5: 32768 bins (16 Kbp each)
//!
//! # Virtual File Offsets
//!
//! BGZF virtual offsets (64-bit):
//! - High 48 bits: Compressed file offset
//! - Low 16 bits: Uncompressed offset within block
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::index::TbiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! // Load index
//! let index = TbiIndex::from_path("data.vcf.gz.tbi")?;
//!
//! // Get metadata
//! println!("Format: {:?}", index.format());
//! println!("References: {}", index.references().len());
//!
//! // Query region
//! let chunks = index.query("chr1", 1000000, 2000000)?;
//! for chunk in chunks {
//!     println!("Chunk: {:016x} - {:016x}", chunk.start.as_raw(), chunk.end.as_raw());
//! }
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::io::bam::index::{Chunk, VirtualOffset};
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// TBI file format magic string
const TBI_MAGIC: &[u8; 4] = b"TBI\x01";

/// File format types recognized by tabix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TbiFormat {
    /// Generic tab-delimited file
    Generic = 0,
    /// SAM format
    Sam = 1,
    /// VCF format
    Vcf = 2,
}

impl TbiFormat {
    /// Parse format from integer
    fn from_i32(value: i32) -> Result<Self> {
        match value {
            0 => Ok(TbiFormat::Generic),
            1 => Ok(TbiFormat::Sam),
            2 => Ok(TbiFormat::Vcf),
            _ => Err(BiometalError::InvalidInput {
                msg: format!("Unknown TBI format code: {}", value),
            }),
        }
    }
}

/// A bin in the hierarchical binning index
#[derive(Debug, Clone)]
pub struct TbiBin {
    /// Bin number (0-37449)
    pub bin_id: u32,
    /// Chunks of data in this bin
    pub chunks: Vec<Chunk>,
}

impl TbiBin {
    /// Create a new bin
    pub fn new(bin_id: u32) -> Self {
        TbiBin {
            bin_id,
            chunks: Vec::new(),
        }
    }
}

/// Reference sequence index data
#[derive(Debug, Clone)]
pub struct TbiReference {
    /// Reference sequence name
    pub name: String,
    /// Bins for this reference (hierarchical spatial index)
    pub bins: Vec<TbiBin>,
    /// Linear index: virtual file offsets for 16kb intervals
    pub intervals: Vec<VirtualOffset>,
}

impl TbiReference {
    /// Create a new reference
    pub fn new(name: String) -> Self {
        TbiReference {
            name,
            bins: Vec::new(),
            intervals: Vec::new(),
        }
    }
}

/// TBI (Tabix) index
///
/// Provides fast random access to sorted, tab-delimited, BGZF-compressed files.
#[derive(Debug, Clone)]
pub struct TbiIndex {
    /// File format type
    format: TbiFormat,
    /// Column for sequence name (0-based)
    col_seq: i32,
    /// Column for start position (0-based)
    col_beg: i32,
    /// Column for end position (0-based, 0 if same as col_beg)
    col_end: i32,
    /// Comment character for header lines
    meta_char: char,
    /// Number of header lines to skip
    skip_lines: i32,
    /// Reference sequences
    references: Vec<TbiReference>,
    /// Reference name to index mapping
    ref_map: HashMap<String, usize>,
}

impl TbiIndex {
    /// Load TBI index from a file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::index::TbiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = TbiIndex::from_path("data.vcf.gz.tbi")?;
    /// println!("Loaded {} references", index.references().len());
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        let mut reader = BufReader::new(file);

        // Check if file is gzip-compressed by reading first 2 bytes
        let mut magic = [0u8; 2];
        reader.read_exact(&mut magic)?;

        // Gzip magic: 0x1f 0x8b
        if magic == [0x1f, 0x8b] {
            // File is gzip-compressed, decompress it
            // Need to reopen file since we can't rewind BufReader easily
            let file = File::open(path.as_ref())?;
            let gz_reader = GzDecoder::new(file);
            let mut buf_reader = BufReader::new(gz_reader);
            Self::parse(&mut buf_reader)
        } else {
            // File is raw TBI format
            // Need to reopen since we already consumed 2 bytes
            let file = File::open(path.as_ref())?;
            let mut reader = BufReader::new(file);
            Self::parse(&mut reader)
        }
    }

    /// Parse TBI index from a reader
    fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        // Read and verify magic string
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != TBI_MAGIC {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "Invalid TBI magic: expected {:?}, got {:?}",
                    TBI_MAGIC, magic
                ),
            });
        }

        // Read header fields
        let n_ref = read_i32(reader)?;
        let format = TbiFormat::from_i32(read_i32(reader)?)?;
        let col_seq = read_i32(reader)?;
        let col_beg = read_i32(reader)?;
        let col_end = read_i32(reader)?;
        let meta = read_i32(reader)?;
        let skip = read_i32(reader)?;
        let l_nm = read_i32(reader)?;

        // Read sequence names
        let mut names_buf = vec![0u8; l_nm as usize];
        reader.read_exact(&mut names_buf)?;
        let names = parse_sequence_names(&names_buf)?;

        if names.len() != n_ref as usize {
            return Err(BiometalError::InvalidInput {
                msg: format!(
                    "TBI header claims {} references but got {} names",
                    n_ref,
                    names.len()
                ),
            });
        }

        // Parse index data for each reference
        let mut references = Vec::new();
        let mut ref_map = HashMap::new();

        for (idx, name) in names.iter().enumerate() {
            let mut reference = TbiReference::new(name.clone());

            // Read binning index
            let n_bin = read_i32(reader)?;
            for _ in 0..n_bin {
                let bin_id = read_u32(reader)?;
                let n_chunk = read_i32(reader)?;

                let mut bin = TbiBin::new(bin_id);
                for _ in 0..n_chunk {
                    let chunk_beg = read_u64(reader)?;
                    let chunk_end = read_u64(reader)?;
                    bin.chunks.push(Chunk::new(
                        VirtualOffset::from_raw(chunk_beg),
                        VirtualOffset::from_raw(chunk_end),
                    ));
                }
                reference.bins.push(bin);
            }

            // Read linear index
            let n_intv = read_i32(reader)?;
            for _ in 0..n_intv {
                let ioff = read_u64(reader)?;
                reference.intervals.push(VirtualOffset::from_raw(ioff));
            }

            ref_map.insert(name.clone(), idx);
            references.push(reference);
        }

        Ok(TbiIndex {
            format,
            col_seq,
            col_beg,
            col_end,
            meta_char: meta as u8 as char,
            skip_lines: skip,
            references,
            ref_map,
        })
    }

    /// Get file format type
    pub fn format(&self) -> TbiFormat {
        self.format
    }

    /// Get column for sequence name (0-based)
    pub fn col_seq(&self) -> i32 {
        self.col_seq
    }

    /// Get column for start position (0-based)
    pub fn col_beg(&self) -> i32 {
        self.col_beg
    }

    /// Get column for end position (0-based)
    pub fn col_end(&self) -> i32 {
        self.col_end
    }

    /// Get comment character
    pub fn meta_char(&self) -> char {
        self.meta_char
    }

    /// Get number of header lines to skip
    pub fn skip_lines(&self) -> i32 {
        self.skip_lines
    }

    /// Get all references
    pub fn references(&self) -> &[TbiReference] {
        &self.references
    }

    /// Get reference by name
    pub fn get_reference(&self, name: &str) -> Option<&TbiReference> {
        self.ref_map.get(name).map(|&idx| &self.references[idx])
    }

    /// Query region and get chunks to read
    ///
    /// Returns list of file chunks that overlap the query region.
    ///
    /// # Arguments
    ///
    /// * `ref_name` - Reference sequence name
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::formats::index::TbiIndex;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let index = TbiIndex::from_path("data.vcf.gz.tbi")?;
    /// let chunks = index.query("chr1", 1000000, 2000000)?;
    /// println!("Need to read {} chunks", chunks.len());
    /// # Ok(())
    /// # }
    /// ```
    pub fn query(&self, ref_name: &str, start: u32, end: u32) -> Result<Vec<Chunk>> {
        let reference = self.get_reference(ref_name).ok_or_else(|| {
            BiometalError::InvalidInput {
                msg: format!("Reference '{}' not found in index", ref_name),
            }
        })?;

        if start >= end {
            return Err(BiometalError::InvalidRange(format!(
                "Invalid range: start ({}) >= end ({})",
                start, end
            )));
        }

        // Get candidate bins that overlap [start, end)
        let bins = reg2bins(start, end);

        // Collect chunks from overlapping bins
        let mut chunks = Vec::new();
        for bin_id in bins {
            if let Some(bin) = reference.bins.iter().find(|b| b.bin_id == bin_id) {
                chunks.extend_from_slice(&bin.chunks);
            }
        }

        // Use linear index to find minimum offset
        let min_offset = get_min_offset(&reference.intervals, start);

        // Filter chunks that are >= min_offset
        chunks.retain(|chunk| chunk.end.as_raw() > min_offset.as_raw());

        // Sort and merge overlapping chunks
        chunks.sort_by_key(|c| c.start.as_raw());
        let merged = merge_chunks(&chunks);

        Ok(merged)
    }
}

/// Parse null-terminated sequence names from buffer
fn parse_sequence_names(buf: &[u8]) -> Result<Vec<String>> {
    let mut names = Vec::new();
    let mut start = 0;

    for (i, &byte) in buf.iter().enumerate() {
        if byte == 0 {
            if i > start {
                let name = std::str::from_utf8(&buf[start..i])
                    .map_err(|e| BiometalError::InvalidInput {
                        msg: format!("Invalid UTF-8 in sequence name: {}", e),
                    })?
                    .to_string();
                names.push(name);
            }
            start = i + 1;
        }
    }

    Ok(names)
}

/// Calculate bin IDs that overlap a region [beg, end)
///
/// Uses the same binning scheme as BAI (37,450 bins)
fn reg2bins(beg: u32, end: u32) -> Vec<u32> {
    let mut bins = Vec::new();
    let end = end - 1; // Make end inclusive for calculation

    // Add bins at each level
    bins.push(0); // Level 0: entire sequence

    for level in 1..=5 {
        let offset = ((1 << (3 * level)) - 1) / 7;
        let shift = 29 - 3 * level;
        let beg_bin = offset + (beg >> shift);
        let end_bin = offset + (end >> shift);

        for bin in beg_bin..=end_bin {
            bins.push(bin);
        }
    }

    bins
}

/// Get minimum virtual offset from linear index
fn get_min_offset(intervals: &[VirtualOffset], beg: u32) -> VirtualOffset {
    let window = (beg >> 14) as usize; // 16kb windows
    if window < intervals.len() {
        intervals[window]
    } else {
        VirtualOffset::from_raw(0)
    }
}

/// Merge overlapping chunks
fn merge_chunks(chunks: &[Chunk]) -> Vec<Chunk> {
    if chunks.is_empty() {
        return Vec::new();
    }

    let mut merged = Vec::new();
    let mut current = chunks[0].clone();

    for chunk in &chunks[1..] {
        if chunk.start.as_raw() <= current.end.as_raw() {
            // Overlapping or adjacent - merge
            if chunk.end.as_raw() > current.end.as_raw() {
                current.end = chunk.end;
            }
        } else {
            // Non-overlapping - save current and start new
            merged.push(current.clone());
            current = chunk.clone();
        }
    }
    merged.push(current);

    merged
}

// Helper functions for reading binary data (little-endian)

fn read_i32<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reg2bins() {
        // Test a small region
        let bins = reg2bins(1000, 2000);
        assert!(bins.contains(&0)); // Level 0
        assert!(!bins.is_empty());

        // Bins should be sorted
        for i in 1..bins.len() {
            assert!(bins[i] >= bins[i - 1]);
        }
    }

    #[test]
    fn test_parse_sequence_names() {
        let buf = b"chr1\0chr2\0chr3\0";
        let names = parse_sequence_names(buf).unwrap();
        assert_eq!(names, vec!["chr1", "chr2", "chr3"]);
    }

    #[test]
    fn test_merge_chunks() {
        let chunks = vec![
            Chunk::new(VirtualOffset::from_raw(100), VirtualOffset::from_raw(200)),
            Chunk::new(VirtualOffset::from_raw(150), VirtualOffset::from_raw(250)),
            Chunk::new(VirtualOffset::from_raw(300), VirtualOffset::from_raw(400)),
        ];

        let merged = merge_chunks(&chunks);
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].start.as_raw(), 100);
        assert_eq!(merged[0].end.as_raw(), 250);
        assert_eq!(merged[1].start.as_raw(), 300);
        assert_eq!(merged[1].end.as_raw(), 400);
    }

    #[test]
    fn test_tbi_format_conversion() {
        assert_eq!(TbiFormat::from_i32(0).unwrap(), TbiFormat::Generic);
        assert_eq!(TbiFormat::from_i32(1).unwrap(), TbiFormat::Sam);
        assert_eq!(TbiFormat::from_i32(2).unwrap(), TbiFormat::Vcf);
        assert!(TbiFormat::from_i32(99).is_err());
    }
}
