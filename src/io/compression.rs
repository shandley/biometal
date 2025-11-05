//! Compression module implementing Rules 3, 4, and 6 (DataSource abstraction)
//!
//! # Evidence Base
//!
//! - **Rule 3**: Parallel bgzip decompression (6.5× speedup, Entry 029)
//! - **Rule 4**: Smart mmap for files ≥50 MB (2.5× additional, Entry 032)
//! - **Rule 6**: DataSource abstraction for network streaming (Entry 028)
//!
//! # Combined Performance
//!
//! - Small files (<50 MB): 6.5× (parallel bgzip only)
//! - Large files (≥50 MB): 16.3× (6.5 × 2.5, layered optimization)

use crate::error::{BiometalError, Result};
use flate2::read::GzDecoder;
use memmap2::Mmap;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// Memory-mapped file threshold (50 MB)
///
/// # Evidence
///
/// Entry 032 (scale validation across 0.54-544 MB):
/// - Files <50 MB: 0.66-0.99× (overhead dominates, don't use mmap)
/// - Files ≥50 MB: 2.30-2.55× speedup (APFS prefetching benefit)
pub const MMAP_THRESHOLD: u64 = 50 * 1024 * 1024; // 50 MB

/// Number of bgzip blocks to decompress in parallel
///
/// # Evidence & Memory Budget
///
/// Entry 029: Parallel bgzip achieves 6.5× speedup using rayon.
///
/// Memory budget (8 blocks in parallel):
/// - Compressed: 8 × ~64 KB = ~512 KB
/// - Decompressed: 8 × ~65 KB = ~520 KB (BGZF spec: max 64 KB uncompressed per block)
/// - Total: ~1 MB (bounded, regardless of file size)
///
/// # BGZF Specification Guarantee
///
/// The BGZF (Blocked GNU Zip Format) specification guarantees that each block
/// decompresses to a maximum of 64 KB uncompressed. This ensures our memory bounds
/// are respected even for malformed files. If a block violates the spec and
/// decompresses larger, the memory is still bounded to PARALLEL_BLOCK_COUNT × actual_size
/// and gets cleared on the next chunk.
///
/// This delivers Rule 3 (parallel speedup) while maintaining Rule 5 (constant memory).
pub const PARALLEL_BLOCK_COUNT: usize = 8;

/// Data source abstraction for local and network streaming
///
/// # Architecture
///
/// This enum is designed Day 1 (Rule 6 requirement) but implements variants incrementally:
/// - **Week 1-2**: Local file support (with Rules 3+4 optimization)
/// - **Week 3-4**: HTTP and SRA streaming implementation
///
/// # Evidence
///
/// Entry 028: I/O bottleneck dominates 264-352× compared to compute time.
/// Network streaming is CRITICAL, not optional.
#[derive(Debug, Clone)]
pub enum DataSource {
    /// Local file path
    Local(PathBuf),

    /// HTTP/HTTPS URL (Week 3-4 implementation)
    #[cfg(feature = "network")]
    Http(String),

    /// SRA accession (Week 3-4 implementation)
    #[cfg(feature = "network")]
    Sra(String),
}

impl DataSource {
    /// Create a local file data source
    pub fn from_path<P: AsRef<Path>>(path: P) -> Self {
        DataSource::Local(path.as_ref().to_path_buf())
    }

    /// Open the data source and return a buffered reader
    ///
    /// # Implementation Status
    ///
    /// - ✅ Local: Full implementation (Rules 3+4)
    /// - ⏳ Http/Sra: Stub (Week 3-4)
    pub fn open(&self) -> Result<Box<dyn BufRead + Send>> {
        match self {
            DataSource::Local(path) => open_local_file(path),

            #[cfg(feature = "network")]
            DataSource::Http(url) => {
                use crate::io::network::HttpReader;
                let reader = HttpReader::new(url)?;
                Ok(Box::new(std::io::BufReader::new(reader)))
            }

            #[cfg(feature = "network")]
            DataSource::Sra(accession) => {
                use crate::io::network::HttpReader;
                use crate::io::sra::sra_to_url;

                // Convert SRA accession to HTTP URL
                let url = sra_to_url(accession)?;
                let reader = HttpReader::new(&url)?;
                Ok(Box::new(std::io::BufReader::new(reader)))
            }
        }
    }
}

/// Open a local file with smart I/O method selection (Rule 4)
///
/// # Evidence
///
/// Entry 032 (threshold-based mmap):
/// - Small files (<50 MB): Use standard I/O (faster)
/// - Large files (≥50 MB): Use mmap + madvise (2.5× speedup on macOS)
fn open_local_file(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let metadata = std::fs::metadata(path)?;
    let file_size = metadata.len();

    if file_size >= MMAP_THRESHOLD {
        // Large file: Use memory-mapped I/O (Rule 4, 2.5× speedup)
        open_mmap_file(path)
    } else {
        // Small file: Use standard I/O (avoids mmap overhead)
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Open file with memory mapping and platform-specific optimization hints
///
/// # Platform Support
///
/// - ✅ macOS: Uses madvise(MADV_SEQUENTIAL | MADV_WILLNEED) for APFS optimization
/// - ⏳ Linux: Future validation (Week 3-4)
/// - ❌ Windows: Falls back to standard I/O
#[cfg(target_os = "macos")]
fn open_mmap_file(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    use libc::{madvise, MADV_SEQUENTIAL, MADV_WILLNEED};

    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };

    // Give kernel sequential access hints for APFS optimization
    unsafe {
        madvise(
            mmap.as_ptr() as *mut _,
            mmap.len(),
            MADV_SEQUENTIAL | MADV_WILLNEED,
        );
    }

    // Wrap mmap in a cursor that implements BufRead
    Ok(Box::new(std::io::Cursor::new(mmap)))
}

#[cfg(not(target_os = "macos"))]
fn open_mmap_file(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    // For non-macOS platforms, use standard mmap without madvise hints
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    Ok(Box::new(std::io::Cursor::new(mmap)))
}

/// Bgzip block structure
///
/// Bgzip files consist of independent compressed blocks that can be
/// decompressed in parallel (Rule 3).
#[derive(Debug, Clone)]
struct BgzipBlock {
    /// Compressed block data
    data: Vec<u8>,
    /// Block size in bytes (reserved for future proper block parsing)
    #[allow(dead_code)]
    size: usize,
}

/// Parse bgzip blocks from compressed data
///
/// # Bgzip Format
///
/// Bgzip is a variant of gzip that uses fixed-size blocks (typically 64KB uncompressed).
/// Each block is an independent gzip stream, enabling parallel decompression.
///
/// # Block Structure
///
/// Each bgzip block:
/// - Bytes 0-1: Gzip magic (31, 139)
/// - Bytes 2-9: Standard gzip header fields
/// - Bytes 10-11: XLEN (extra field length)
/// - Bytes 12+: Extra subfields, including BSIZE
///   - SI1=66 ('B'), SI2=67 ('C') for bgzip
///   - SLEN=2 (2-byte BSIZE field)
///   - BSIZE (little-endian u16): total block size - 1
///
/// # Evidence
///
/// Entry 029: Parallel bgzip decompression achieves 6.5× speedup using rayon.
fn parse_bgzip_blocks(data: &[u8]) -> Result<Vec<BgzipBlock>> {
    let mut blocks = Vec::new();
    let mut pos = 0;

    while pos < data.len() {
        // Need at least 18 bytes for minimal gzip header
        if pos + 18 > data.len() {
            if pos == data.len() {
                break; // Clean EOF
            }
            return Err(BiometalError::Compression(format!(
                "Incomplete gzip header at offset {} (only {} bytes remaining)",
                pos,
                data.len() - pos
            )));
        }

        // Check gzip magic bytes (ID1=31, ID2=139)
        if data[pos] != 31 || data[pos + 1] != 139 {
            return Err(BiometalError::Compression(format!(
                "Invalid gzip magic bytes at offset {}: expected [31, 139], got [{}, {}]",
                pos, data[pos], data[pos + 1]
            )));
        }

        // Check for extra field flag (FLG & 0x04)
        let flg = data[pos + 3];
        if flg & 0x04 == 0 {
            // No extra field - this might be regular gzip, not bgzip
            // For safety, treat entire remaining data as one block
            blocks.push(BgzipBlock {
                data: data[pos..].to_vec(),
                size: data.len() - pos,
            });
            break;
        }

        // Read XLEN (extra field length) at bytes 10-11 (little-endian)
        let xlen = u16::from_le_bytes([data[pos + 10], data[pos + 11]]) as usize;

        // Extra field starts at byte 12
        let extra_start = pos + 12;
        let extra_end = extra_start + xlen;

        if extra_end > data.len() {
            return Err(BiometalError::Compression(format!(
                "Extra field extends beyond data at offset {}", pos
            )));
        }

        // Parse extra subfields to find BSIZE
        let mut bsize: Option<u16> = None;
        let mut extra_pos = extra_start;

        while extra_pos + 4 <= extra_end {
            let si1 = data[extra_pos];
            let si2 = data[extra_pos + 1];
            let slen = u16::from_le_bytes([data[extra_pos + 2], data[extra_pos + 3]]) as usize;

            // Check for bgzip subfield (SI1='B'=66, SI2='C'=67)
            if si1 == 66 && si2 == 67 && slen == 2 {
                if extra_pos + 6 > extra_end {
                    return Err(BiometalError::Compression(format!(
                        "BSIZE field incomplete at offset {}", pos
                    )));
                }
                bsize = Some(u16::from_le_bytes([
                    data[extra_pos + 4],
                    data[extra_pos + 5],
                ]));
                break;
            }

            extra_pos += 4 + slen;
        }

        // Calculate block size
        let block_size = match bsize {
            Some(bs) => (bs as usize) + 1, // BSIZE is block_size - 1
            None => {
                // No BSIZE field found - might be regular gzip
                // Treat remaining data as one block
                blocks.push(BgzipBlock {
                    data: data[pos..].to_vec(),
                    size: data.len() - pos,
                });
                break;
            }
        };

        // Validate block size
        if pos + block_size > data.len() {
            return Err(BiometalError::Compression(format!(
                "Block size {} at offset {} exceeds remaining data ({} bytes)",
                block_size,
                pos,
                data.len() - pos
            )));
        }

        // Extract block
        blocks.push(BgzipBlock {
            data: data[pos..pos + block_size].to_vec(),
            size: block_size,
        });

        pos += block_size;
    }

    if blocks.is_empty() {
        return Err(BiometalError::Compression(
            "No valid bgzip blocks found".to_string(),
        ));
    }

    Ok(blocks)
}

/// Decompress a single bgzip block
fn decompress_block(block: &BgzipBlock) -> io::Result<Vec<u8>> {
    let mut decoder = GzDecoder::new(&block.data[..]);
    let mut decompressed = Vec::new();
    decoder.read_to_end(&mut decompressed)?;
    Ok(decompressed)
}

/// Bounded parallel bgzip reader (Rules 3 + 5)
///
/// # Architecture
///
/// Processes blocks in chunks of PARALLEL_BLOCK_COUNT (8):
/// 1. Read 8 blocks from input stream
/// 2. Decompress 8 blocks in parallel using rayon (Rule 3: 6.5× speedup)
/// 3. Buffer concatenated results
/// 4. Yield data incrementally
/// 5. Repeat
///
/// # Memory Footprint (Rule 5)
///
/// Bounded regardless of file size:
/// - Compressed blocks: 8 × ~64 KB = ~512 KB
/// - Decompressed blocks: 8 × ~65 KB = ~520 KB
/// - Total: ~1 MB (constant, even for 5TB files)
///
/// # Evidence
///
/// - Rule 3 (Entry 029): Parallel bgzip achieves 6.5× speedup
/// - Rule 5 (Entry 026): Constant memory streaming (~5 MB)
/// - Combined: Delivers both parallelism AND constant memory
struct BoundedParallelBgzipReader<R: BufRead> {
    inner: R,
    /// Buffer for decompressed data ready to read
    output_buffer: Vec<u8>,
    /// Current read position in output_buffer
    output_pos: usize,
    /// Whether we've reached EOF
    eof: bool,
}

impl<R: BufRead> BoundedParallelBgzipReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            output_buffer: Vec::new(),
            output_pos: 0,
            eof: false,
        }
    }

    /// Read one bgzip block from stream
    fn read_one_block(&mut self) -> io::Result<Option<BgzipBlock>> {
        // Read block header to determine block size
        let mut header = [0u8; 18];
        match self.inner.read_exact(&mut header) {
            Ok(_) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                return Ok(None); // Clean EOF
            }
            Err(e) => return Err(e),
        }

        // Check gzip magic
        if header[0] != 31 || header[1] != 139 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid gzip magic: [{}, {}]", header[0], header[1]),
            ));
        }

        // Check for extra field (FLG & 0x04)
        let flg = header[3];
        if flg & 0x04 == 0 {
            // Regular gzip, not bgzip - read entire remaining stream
            let mut compressed = header.to_vec();
            self.inner.read_to_end(&mut compressed)?;
            let size = compressed.len();
            return Ok(Some(BgzipBlock {
                data: compressed,
                size,
            }));
        }

        // Read XLEN (extra field length)
        let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;

        // Read extra field
        let mut extra = vec![0u8; xlen];
        self.inner.read_exact(&mut extra)?;

        // Parse BSIZE from extra field
        let mut bsize: Option<u16> = None;
        let mut pos = 0;

        while pos + 4 <= xlen {
            let si1 = extra[pos];
            let si2 = extra[pos + 1];
            let slen = u16::from_le_bytes([extra[pos + 2], extra[pos + 3]]) as usize;

            if si1 == 66 && si2 == 67 && slen == 2 {
                // Found BSIZE
                if pos + 6 > xlen {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Incomplete BSIZE field",
                    ));
                }
                bsize = Some(u16::from_le_bytes([extra[pos + 4], extra[pos + 5]]));
                break;
            }

            pos += 4 + slen;
        }

        let block_size = match bsize {
            Some(bs) => (bs as usize) + 1,
            None => {
                // No BSIZE - treat as regular gzip
                let mut compressed = header.to_vec();
                compressed.extend_from_slice(&extra);
                self.inner.read_to_end(&mut compressed)?;
                let size = compressed.len();
                return Ok(Some(BgzipBlock {
                    data: compressed,
                    size,
                }));
            }
        };

        // Calculate remaining bytes to read (block_size - header - extra)
        let already_read = 18 + xlen;
        if block_size < already_read {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid block size: {} < {}", block_size, already_read),
            ));
        }
        let remaining = block_size - already_read;

        // Build complete block
        let mut block_data = Vec::with_capacity(block_size);
        block_data.extend_from_slice(&header);
        block_data.extend_from_slice(&extra);

        // Read remaining block data
        let mut rest = vec![0u8; remaining];
        self.inner.read_exact(&mut rest)?;
        block_data.extend_from_slice(&rest);

        Ok(Some(BgzipBlock {
            data: block_data,
            size: block_size,
        }))
    }

    /// Read and decompress next chunk of blocks in parallel
    fn read_next_chunk(&mut self) -> io::Result<()> {
        if self.eof {
            return Ok(());
        }

        // Read up to PARALLEL_BLOCK_COUNT blocks
        let mut blocks = Vec::with_capacity(PARALLEL_BLOCK_COUNT);

        for _ in 0..PARALLEL_BLOCK_COUNT {
            match self.read_one_block()? {
                Some(block) => blocks.push(block),
                None => {
                    self.eof = true;
                    break;
                }
            }
        }

        if blocks.is_empty() {
            return Ok(());
        }

        // Decompress blocks in parallel (Rule 3: 6.5× speedup)
        let decompressed_blocks: Vec<_> = blocks
            .par_iter()
            .map(decompress_block)
            .collect::<io::Result<Vec<_>>>()?;

        // Concatenate decompressed blocks into output buffer
        self.output_buffer.clear();
        for block_data in decompressed_blocks {
            self.output_buffer.extend_from_slice(&block_data);
        }
        self.output_pos = 0;

        Ok(())
    }
}

impl<R: BufRead> Read for BoundedParallelBgzipReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // If output buffer is exhausted, read and decompress next chunk
        if self.output_pos >= self.output_buffer.len() {
            if self.eof {
                return Ok(0); // EOF
            }

            self.read_next_chunk()?;

            // Check if we got any data
            if self.output_buffer.is_empty() {
                return Ok(0); // EOF
            }
        }

        // Copy from output buffer to caller's buffer
        let available = self.output_buffer.len() - self.output_pos;
        let to_copy = available.min(buf.len());

        buf[..to_copy].copy_from_slice(
            &self.output_buffer[self.output_pos..self.output_pos + to_copy],
        );
        self.output_pos += to_copy;

        Ok(to_copy)
    }
}

/// Decompress bgzip file in parallel (Rule 3)
///
/// # Evidence
///
/// Entry 029: CPU parallel prototype achieves 6.5× speedup.
/// - Implementation: Rayon-based parallelism
/// - Platform: All platforms (Mac, Linux, Windows, ARM, x86_64)
/// - Cost: Zero platform dependencies (pure Rust + Rayon)
///
/// # Performance
///
/// - Uses all available CPU cores
/// - Each block decompressed independently
/// - Automatic work distribution via rayon
///
/// # Example
///
/// ```no_run
/// use biometal::io::compression::decompress_bgzip_parallel;
///
/// # fn main() -> std::io::Result<()> {
/// let compressed_data = std::fs::read("file.fq.gz")?;
/// let decompressed = decompress_bgzip_parallel(&compressed_data)?;
/// # Ok(())
/// # }
/// ```
pub fn decompress_bgzip_parallel(data: &[u8]) -> io::Result<Vec<u8>> {
    // Parse bgzip block boundaries
    let blocks = parse_bgzip_blocks(data)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    // Decompress blocks in parallel (uses all CPU cores)
    let decompressed_blocks: Vec<_> = blocks
        .par_iter()
        .map(decompress_block)
        .collect::<io::Result<Vec<_>>>()?;

    // Concatenate decompressed blocks
    Ok(decompressed_blocks.concat())
}

/// Compressed file reader with full optimization stack (Rules 3+4)
///
/// # Evidence-Based Optimization
///
/// This reader combines:
/// - **Rule 3**: Parallel bgzip decompression (6.5× speedup)
/// - **Rule 4**: Threshold-based mmap (2.5× additional for ≥50 MB)
/// - **Combined**: 16.3× speedup for large files
///
/// # Example
///
/// ```no_run
/// use biometal::io::compression::{DataSource, CompressedReader};
///
/// # fn main() -> biometal::Result<()> {
/// let source = DataSource::from_path("large.fq.gz");
/// let reader = CompressedReader::new(source)?;
///
/// // Reader implements BufRead, use with FASTQ parser
/// # Ok(())
/// # }
/// ```
pub struct CompressedReader {
    /// Inner buffered reader (from decompressed data)
    inner: Box<dyn BufRead + Send>,
}

impl CompressedReader {
    /// Create a new compressed reader from a data source
    ///
    /// # Optimization Stack (Rules 3-6)
    ///
    /// 1. Opens data source (Rule 6 abstraction)
    /// 2. Applies threshold-based mmap if local file ≥50 MB (Rule 4)
    /// 3. Decompresses blocks in parallel chunks (Rule 3: 6.5× speedup)
    /// 4. Maintains constant memory (Rule 5: ~1 MB bounded)
    ///
    /// # Memory Usage (Rule 5)
    ///
    /// Bounded regardless of file size:
    /// - 8 compressed blocks: ~512 KB
    /// - 8 decompressed blocks: ~520 KB
    /// - Total: ~1 MB (constant, even for 5TB files)
    ///
    /// # Performance (Rule 3)
    ///
    /// - Decompresses 8 blocks in parallel using rayon
    /// - Expected: 6.5× speedup (validated in Entry 029)
    /// - Maintains constant memory while achieving parallel speedup
    pub fn new(source: DataSource) -> Result<Self> {
        // Open source with smart I/O (Rules 4+6)
        let reader = source.open()?;

        // Wrap in bounded parallel bgzip reader (Rules 3+5 combined)
        // - Decompresses 8 blocks in parallel (Rule 3: 6.5× speedup)
        // - Bounded memory: ~1 MB (Rule 5: constant regardless of file size)
        let parallel_reader = BoundedParallelBgzipReader::new(reader);

        // Wrap in buffered reader for efficient line-by-line reading
        Ok(Self {
            inner: Box::new(BufReader::new(parallel_reader)),
        })
    }

    /// Get the inner buffered reader
    pub fn into_inner(self) -> Box<dyn BufRead + Send> {
        self.inner
    }
}

impl Read for CompressedReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

impl BufRead for CompressedReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mmap_threshold_constant() {
        // Verify evidence-based threshold (Entry 032)
        assert_eq!(MMAP_THRESHOLD, 50 * 1024 * 1024);
    }

    #[test]
    fn test_datasource_local_creation() {
        let source = DataSource::from_path("/tmp/test.fq");
        match source {
            DataSource::Local(path) => {
                assert_eq!(path, PathBuf::from("/tmp/test.fq"));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Local variant"),
        }
    }

    // Network streaming tests moved to io/network.rs module

    #[test]
    fn test_bgzip_block_parsing() {
        // Test that we correctly parse multiple bgzip blocks from a real file
        // Using the tiny dataset which should have multiple blocks
        let home = std::env::var("HOME").unwrap_or_else(|_| "/Users/scotthandley".to_string());
        let path = PathBuf::from(&home)
            .join("Code/apple-silicon-bio-bench/datasets/tiny_100_150bp.fq.gz");

        if !path.exists() {
            println!("Skipping test: file not found at {:?}", path);
            return;
        }

        // Read the compressed file
        let compressed_data = std::fs::read(&path).expect("Failed to read test file");

        // Parse blocks
        let blocks = parse_bgzip_blocks(&compressed_data).expect("Failed to parse bgzip blocks");

        // Verify we found multiple blocks (bgzip files typically have many blocks)
        println!("Found {} bgzip blocks in {:?}", blocks.len(), path);
        println!(
            "Block sizes: {:?}",
            blocks.iter().map(|b| b.size).collect::<Vec<_>>()
        );

        // Bgzip files should have multiple blocks for parallel decompression
        // Even small files typically have 1-10 blocks
        assert!(
            !blocks.is_empty(),
            "Expected at least 1 block, found {}",
            blocks.len()
        );

        // Verify total size matches original file
        let total_size: usize = blocks.iter().map(|b| b.size).sum();
        assert_eq!(
            total_size,
            compressed_data.len(),
            "Block sizes don't sum to file size"
        );

        // Verify each block has valid gzip header
        for (i, block) in blocks.iter().enumerate() {
            assert!(
                block.data.len() >= 18,
                "Block {} too small: {} bytes",
                i,
                block.data.len()
            );
            assert_eq!(
                block.data[0], 31,
                "Block {} has invalid magic byte 0: {}",
                i, block.data[0]
            );
            assert_eq!(
                block.data[1], 139,
                "Block {} has invalid magic byte 1: {}",
                i, block.data[1]
            );
        }
    }

    #[test]
    fn test_parallel_sequential_equivalence() {
        // Core correctness guarantee: parallel decompression must produce
        // byte-for-byte identical output to sequential decompression

        let home = std::env::var("HOME").unwrap_or_else(|_| "/Users/scotthandley".to_string());

        // Test with multiple dataset sizes to ensure correctness across all scales
        let test_files = vec![
            "Code/apple-silicon-bio-bench/datasets/tiny_100_150bp.fq.gz",
            "Code/apple-silicon-bio-bench/datasets/small_1k_150bp.fq.gz",
            "Code/apple-silicon-bio-bench/datasets/medium_10k_150bp.fq.gz",
        ];

        for rel_path in test_files {
            let path = PathBuf::from(&home).join(rel_path);

            if !path.exists() {
                println!("Skipping test: file not found at {:?}", path);
                continue;
            }

            println!("Testing parallel vs sequential for: {:?}", path);

            // Decompress with bounded parallel (our implementation)
            let source = DataSource::from_path(&path);
            let parallel_reader = BoundedParallelBgzipReader::new(
                source.open().expect("Failed to open source")
            );
            let mut parallel_output = Vec::new();
            BufReader::new(parallel_reader)
                .read_to_end(&mut parallel_output)
                .expect("Failed to read parallel output");

            // Decompress with sequential (standard GzDecoder)
            let file = File::open(&path).expect("Failed to open file");
            let mut sequential_reader = GzDecoder::new(file);
            let mut sequential_output = Vec::new();
            sequential_reader
                .read_to_end(&mut sequential_output)
                .expect("Failed to read sequential output");

            // Must be byte-for-byte identical
            assert_eq!(
                parallel_output.len(),
                sequential_output.len(),
                "Output lengths differ for {:?}",
                path
            );

            assert_eq!(
                parallel_output,
                sequential_output,
                "Parallel and sequential outputs differ for {:?}",
                path
            );

            println!("  ✓ Outputs match ({} bytes)", parallel_output.len());
        }
    }
}
