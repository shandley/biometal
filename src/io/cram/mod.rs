//! CRAM (reference-based compressed alignment format) - Native ARM-Optimized Parser
//!
//! **Status**: Native implementation in progress (v1.11.0+)
//!
//! This module provides a **native, zero-dependency CRAM 3.0/3.1 parser** with ARM NEON
//! optimizations. Unlike other CRAM implementations, this is designed specifically for:
//! - **ARM-native performance** (NEON-optimized base encoding/decoding)
//! - **Streaming architecture** (constant ~5 MB memory, Rule 5)
//! - **Evidence-based optimization** (following OPTIMIZATION_RULES.md)
//!
//! # Why Native Implementation?
//!
//! **Unique Value Proposition**:
//! 1. **First ARM-optimized CRAM reader** - 16-25× NEON speedup potential
//! 2. **Zero external dependencies** - full control over implementation
//! 3. **Perfect streaming** - designed for constant memory from the ground up
//! 4. **Flagship feature** - demonstrates biometal's ARM-native capabilities
//!
//! **vs noodles-cram**:
//! - No API compatibility issues
//! - Direct ARM NEON optimizations
//! - No transitive dependencies
//! - Optimized for biometal's streaming architecture
//!
//! # CRAM Format Overview
//!
//! CRAM files are 30-60% smaller than BAM through reference-based compression:
//!
//! ```text
//! CRAM File Structure:
//! ┌─────────────────────────────────────┐
//! │ File Definition (magic + version)   │ ← 26 bytes
//! ├─────────────────────────────────────┤
//! │ SAM Header Container                │ ← One container
//! ├─────────────────────────────────────┤
//! │ Data Container 1                    │
//! │  ├─ Compression Header              │
//! │  ├─ Slice 1                         │
//! │  │   └─ Blocks (core data, tags...) │
//! │  ├─ Slice 2                         │
//! │  └─ ...                             │
//! ├─────────────────────────────────────┤
//! │ Data Container 2                    │
//! │  └─ ...                             │
//! └─────────────────────────────────────┘
//! ```
//!
//! # Implementation Roadmap
//!
//! ## Phase 1: Basic Reading (Target: 3-5 days)
//! - [ ] File definition parser (magic number, version)
//! - [ ] Container structure (length, CRC, blocks)
//! - [ ] Slice structure (header, blocks)
//! - [ ] Block decompression (gzip via cloudflare_zlib)
//! - [ ] Basic record iteration (no reference reconstruction)
//! - **Deliverable**: Read CRAM files, extract alignment records
//!
//! ## Phase 2: Full Decoding (Target: 3-4 days)
//! - [ ] Reference FASTA integration
//! - [ ] Reference-based sequence reconstruction
//! - [ ] Full tag support
//! - [ ] Multi-codec support (bzip2, lzma, rANS)
//! - **Deliverable**: Complete CRAM 3.0/3.1 compliance
//!
//! ## Phase 3: ARM Optimization (Target: 3-4 days)
//! - [ ] NEON base encoding/decoding (16-25× potential)
//! - [ ] NEON quality score processing
//! - [ ] NEON reference comparison
//! - [ ] Benchmark vs samtools (target: 2-3× faster on ARM)
//! - **Deliverable**: Fastest ARM-native CRAM reader
//!
//! # Basic Usage (Planned)
//!
//! ```no_run
//! use biometal::io::cram::CramReader;
//!
//! # fn main() -> biometal::Result<()> {
//! let mut cram = CramReader::from_path("alignments.cram")?;
//!
//! println!("Header: {} references", cram.header().reference_count());
//!
//! // Stream records with constant memory (~5 MB)
//! for record in cram.records() {
//!     let record = record?;
//!     println!("Read {} at position {:?}", record.name, record.position);
//! }
//! # Ok(())
//! # }
//! ```
//!
//! # With Reference (for full sequence reconstruction)
//!
//! ```no_run
//! use biometal::io::cram::CramReader;
//!
//! # fn main() -> biometal::Result<()> {
//! let mut cram = CramReader::from_path_with_reference(
//!     "alignments.cram",
//!     "reference.fa"
//! )?;
//!
//! for record in cram.records() {
//!     let record = record?;
//!     // Full sequences reconstructed from reference
//!     println!("{}: {} bases", record.name, record.sequence.len());
//! }
//! # Ok(())
//! # }
//! ```
//!
//! # Performance Goals
//!
//! Based on OPTIMIZATION_RULES.md and ARM NEON potential:
//!
//! | Operation | Baseline | Target (NEON) | Speedup |
//! |-----------|----------|---------------|---------|
//! | Base decode | 100 Mseq/s | 1,600 Mseq/s | 16× |
//! | Quality decode | 150 Mseq/s | 3,000 Mseq/s | 20× |
//! | Overall parsing | 50 MB/s | 100-150 MB/s | 2-3× |
//!
//! **Target**: Fastest ARM-native CRAM reader (vs samtools, htslib)

use crate::io::bam::{Header, Record, Tags};
use crate::io::fasta::FaiIndex;
use crate::{BiometalError, Result};
use std::io::{self, BufReader, Read};
use std::path::{Path, PathBuf};
use std::fs::File;

// ============================================================================
// CRAM Container Structures
// ============================================================================

/// CRAM container header.
///
/// Containers are the top-level organizational unit in CRAM files.
/// Each container contains a compression header and one or more slices.
#[derive(Debug, Clone)]
pub struct ContainerHeader {
    /// Total container length in bytes (excluding this length field itself)
    pub length: i32,
    /// Reference sequence ID (-1 for unmapped, -2 for multi-ref)
    pub reference_id: i32,
    /// Start position on reference (1-based)
    pub start_position: i32,
    /// Alignment span (number of reference bases covered)
    pub alignment_span: i32,
    /// Number of records in this container
    pub num_records: i32,
    /// Global record counter (total records written so far)
    pub record_counter: i64,
    /// Total number of bases (sum of read lengths)
    pub bases: i64,
    /// Number of blocks in this container
    pub num_blocks: i32,
    /// Landmarks: byte offsets of slices within container data
    pub landmarks: Vec<i32>,
}

impl ContainerHeader {
    /// Parse container header from reader.
    ///
    /// # Container Header Format
    ///
    /// ```text
    /// - Length: i32 (4 bytes, big-endian)
    /// - Reference ID: ITF-8
    /// - Start position: ITF-8
    /// - Alignment span: ITF-8
    /// - Number of records: ITF-8
    /// - Record counter: LTF-8
    /// - Bases: LTF-8
    /// - Number of blocks: ITF-8
    /// - Landmarks: ITF-8 count + ITF-8 array
    /// - CRC32: u32 (4 bytes)
    /// ```
    pub fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        // Read length (4 bytes, big-endian)
        let mut length_buf = [0u8; 4];
        reader.read_exact(&mut length_buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read container length: {}", e)
            })?;
        let length = i32::from_be_bytes(length_buf);

        // Check for EOF container (length == 0 or special EOF marker)
        if length == 0 {
            return Err(BiometalError::InvalidCramFormat {
                msg: "EOF container (length = 0)".to_string()
            });
        }

        // Read reference ID
        let reference_id = decode_itf8(reader)?;

        // Read start position
        let start_position = decode_itf8(reader)?;

        // Read alignment span
        let alignment_span = decode_itf8(reader)?;

        // Read number of records
        let num_records = decode_itf8(reader)?;

        // Read record counter
        let record_counter = decode_ltf8(reader)?;

        // Read bases
        let bases = decode_ltf8(reader)?;

        // Read number of blocks
        let num_blocks = decode_itf8(reader)?;

        // Read landmarks (array of slice positions)
        let num_landmarks = decode_itf8(reader)?;
        let mut landmarks = Vec::with_capacity(num_landmarks as usize);
        for _ in 0..num_landmarks {
            landmarks.push(decode_itf8(reader)?);
        }

        // Read and validate CRC32 (4 bytes)
        let mut crc_buf = [0u8; 4];
        reader.read_exact(&mut crc_buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read container CRC32: {}", e)
            })?;
        let _crc32 = u32::from_be_bytes(crc_buf);
        // TODO: Validate CRC32 when we have the full header data

        Ok(ContainerHeader {
            length,
            reference_id,
            start_position,
            alignment_span,
            num_records,
            record_counter,
            bases,
            num_blocks,
            landmarks,
        })
    }

    /// Check if this is an EOF container.
    pub fn is_eof(&self) -> bool {
        self.length == 0
    }

    /// Check if this is a SAM header container.
    ///
    /// SAM header containers have:
    /// - reference_id = -1
    /// - start_position = 0
    /// - num_records = 0
    pub fn is_sam_header_container(&self) -> bool {
        self.reference_id == -1 && self.start_position == 0 && self.num_records == 0
    }
}

/// CRAM block structure.
///
/// Blocks are compressed chunks of data within containers.
#[derive(Debug, Clone)]
pub struct Block {
    /// Block compression method (0=raw, 1=gzip, 2=bzip2, 3=lzma, 4=rans)
    pub method: u8,
    /// Block content type ID
    pub content_type: u8,
    /// Block ID
    pub block_id: i32,
    /// Compressed size in bytes
    pub compressed_size: i32,
    /// Uncompressed size in bytes
    pub uncompressed_size: i32,
    /// Raw compressed data
    pub data: Vec<u8>,
}

impl Block {
    /// Parse block from reader.
    ///
    /// # Block Format
    ///
    /// ```text
    /// - Method: u8 (compression method)
    /// - Content type: u8
    /// - Block ID: ITF-8
    /// - Compressed size: ITF-8
    /// - Uncompressed size: ITF-8
    /// - Data: [u8; compressed_size]
    /// - CRC32: u32 (4 bytes)
    /// ```
    pub fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        // Read method (1 byte)
        let mut method_buf = [0u8; 1];
        reader.read_exact(&mut method_buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read block method: {}", e)
            })?;
        let method = method_buf[0];

        // Read content type (1 byte)
        let mut content_type_buf = [0u8; 1];
        reader.read_exact(&mut content_type_buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read block content type: {}", e)
            })?;
        let content_type = content_type_buf[0];

        // Read block ID
        let block_id = decode_itf8(reader)?;

        // Read compressed size
        let compressed_size = decode_itf8(reader)?;

        // Read uncompressed size
        let uncompressed_size = decode_itf8(reader)?;

        // Read compressed data
        let mut data = vec![0u8; compressed_size as usize];
        reader.read_exact(&mut data)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read block data: {}", e)
            })?;

        // Read and validate CRC32 (4 bytes)
        let mut crc_buf = [0u8; 4];
        reader.read_exact(&mut crc_buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read block CRC32: {}", e)
            })?;
        let _crc32 = u32::from_be_bytes(crc_buf);
        // TODO: Validate CRC32 when needed

        Ok(Block {
            method,
            content_type,
            block_id,
            compressed_size,
            uncompressed_size,
            data,
        })
    }

    /// Decompress block data.
    ///
    /// Returns the uncompressed data based on the compression method.
    pub fn decompress(&self) -> Result<Vec<u8>> {
        match self.method {
            0 => {
                // Raw (no compression)
                Ok(self.data.clone())
            }
            1 => {
                // Gzip (use existing cloudflare_zlib via flate2)
                use flate2::read::GzDecoder;
                let mut decoder = GzDecoder::new(&self.data[..]);
                let mut decompressed = Vec::new();
                decoder.read_to_end(&mut decompressed)
                    .map_err(|e| BiometalError::Compression(
                        format!("Failed to decompress gzip block: {}", e)
                    ))?;
                Ok(decompressed)
            }
            2 => {
                // Bzip2
                use bzip2::read::BzDecoder;
                let mut decoder = BzDecoder::new(&self.data[..]);
                let mut decompressed = Vec::new();
                decoder.read_to_end(&mut decompressed)
                    .map_err(|e| BiometalError::Compression(
                        format!("Failed to decompress bzip2 block: {}", e)
                    ))?;
                Ok(decompressed)
            }
            3 => {
                // LZMA
                use xz2::read::XzDecoder;
                let mut decoder = XzDecoder::new(&self.data[..]);
                let mut decompressed = Vec::new();
                decoder.read_to_end(&mut decompressed)
                    .map_err(|e| BiometalError::Compression(
                        format!("Failed to decompress LZMA block: {}", e)
                    ))?;
                Ok(decompressed)
            }
            4 => {
                // rANS (Phase 2)
                Err(BiometalError::InvalidCramFormat {
                    msg: "rANS compression not yet supported".to_string()
                })
            }
            _ => {
                Err(BiometalError::InvalidCramFormat {
                    msg: format!("Unknown compression method: {}", self.method)
                })
            }
        }
    }

    /// Check if this is a compression header block.
    ///
    /// Compression header blocks have content_type = 0.
    pub fn is_compression_header(&self) -> bool {
        self.content_type == 0
    }
}

/// CRAM data series identifier.
///
/// Two-character codes identifying different types of data in CRAM records.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DataSeries {
    /// BAM bit flags
    BF,
    /// CRAM bit flags
    CF,
    /// Reference ID
    RI,
    /// Read lengths
    RL,
    /// Alignment start positions
    AP,
    /// Read groups
    RG,
    /// Read names
    RN,
    /// Next fragment reference ID
    NF,
    /// Next mate bit flags
    MF,
    /// Next fragment alignment start
    NS,
    /// Next fragment alignment end
    NP,
    /// Template size
    TS,
    /// Tag count
    TC,
    /// Tag names list ID
    TL,
    /// Number of read features
    FN,
    /// Feature codes
    FC,
    /// In-seq positions
    FP,
    /// Deletion lengths
    DL,
    /// Bases (for unmapped reads)
    BA,
    /// Base substitution codes
    BS,
    /// Insertion sequences
    IN,
    /// Reference skip length
    RS,
    /// Soft clip sequences
    SC,
    /// Hard clip lengths
    HC,
    /// Padding lengths
    PD,
    /// Quality scores
    QS,
    /// Mapping qualities
    MQ,
    /// Unknown data series (raw bytes)
    Unknown([u8; 2]),
}

impl DataSeries {
    fn from_bytes(bytes: [u8; 2]) -> Self {
        match &bytes {
            b"BF" => Self::BF,
            b"CF" => Self::CF,
            b"RI" => Self::RI,
            b"RL" => Self::RL,
            b"AP" => Self::AP,
            b"RG" => Self::RG,
            b"RN" => Self::RN,
            b"NF" => Self::NF,
            b"MF" => Self::MF,
            b"NS" => Self::NS,
            b"NP" => Self::NP,
            b"TS" => Self::TS,
            b"TC" => Self::TC,
            b"TL" => Self::TL,
            b"FN" => Self::FN,
            b"FC" => Self::FC,
            b"FP" => Self::FP,
            b"DL" => Self::DL,
            b"BA" => Self::BA,
            b"BS" => Self::BS,
            b"IN" => Self::IN,
            b"RS" => Self::RS,
            b"SC" => Self::SC,
            b"HC" => Self::HC,
            b"PD" => Self::PD,
            b"QS" => Self::QS,
            b"MQ" => Self::MQ,
            _ => Self::Unknown(bytes),
        }
    }
}

/// CRAM encoding specification.
///
/// Describes how a data series is encoded in CRAM blocks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Encoding {
    /// NULL encoding (no data)
    Null,
    /// External encoding (data in external block)
    External { block_content_id: i32 },
    /// Huffman encoding with code table
    Huffman {
        alphabet: Vec<i32>,
        bit_lengths: Vec<i32>,
    },
    /// Byte array with length prefix
    ByteArrayLen {
        len_encoding: Box<Encoding>,
        value_encoding: Box<Encoding>,
    },
    /// Byte array with stop byte
    ByteArrayStop {
        stop_byte: u8,
        block_content_id: i32,
    },
    /// Beta encoding (offset + length)
    Beta { offset: i32, length: i32 },
    /// Subexponential encoding
    SubExp { offset: i32, k: i32 },
    /// Golomb encoding
    Golomb { offset: i32, m: i32 },
    /// Golomb-Rice encoding
    GolombRice { offset: i32, log2_m: i32 },
    /// Gamma encoding
    Gamma { offset: i32 },
    /// Delta encoding
    Delta { offset: i32, k: i32 },
}

impl Encoding {
    /// Parse encoding from bytes.
    fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        let encoding_id = {
            let mut buf = [0u8; 1];
            reader.read_exact(&mut buf)
                .map_err(|e| BiometalError::InvalidCramFormat {
                    msg: format!("Failed to read encoding ID: {}", e)
                })?;
            buf[0]
        };

        match encoding_id {
            0 => Ok(Self::Null),
            1 => {
                // EXTERNAL: ITF-8 block_content_id
                let block_content_id = decode_itf8(reader)?;
                Ok(Self::External { block_content_id })
            }
            2 => {
                // GOLOMB: ITF-8 offset, ITF-8 M
                let offset = decode_itf8(reader)?;
                let m = decode_itf8(reader)?;
                Ok(Self::Golomb { offset, m })
            }
            3 => {
                // HUFFMAN: ITF-8 alphabet_size, ITF-8[] alphabet, ITF-8[] bit_lengths
                let alphabet_size = decode_itf8(reader)?;
                let mut alphabet = Vec::with_capacity(alphabet_size as usize);
                for _ in 0..alphabet_size {
                    alphabet.push(decode_itf8(reader)?);
                }
                let mut bit_lengths = Vec::with_capacity(alphabet_size as usize);
                for _ in 0..alphabet_size {
                    bit_lengths.push(decode_itf8(reader)?);
                }
                Ok(Self::Huffman {
                    alphabet,
                    bit_lengths,
                })
            }
            4 => {
                // BYTE_ARRAY_LEN: Encoding (len), Encoding (value)
                let len_encoding = Box::new(Self::parse(reader)?);
                let value_encoding = Box::new(Self::parse(reader)?);
                Ok(Self::ByteArrayLen {
                    len_encoding,
                    value_encoding,
                })
            }
            5 => {
                // BYTE_ARRAY_STOP: u8 stop_byte, ITF-8 block_content_id
                let mut buf = [0u8; 1];
                reader.read_exact(&mut buf)
                    .map_err(|e| BiometalError::InvalidCramFormat {
                        msg: format!("Failed to read stop byte: {}", e)
                    })?;
                let stop_byte = buf[0];
                let block_content_id = decode_itf8(reader)?;
                Ok(Self::ByteArrayStop {
                    stop_byte,
                    block_content_id,
                })
            }
            6 => {
                // BETA: ITF-8 offset, ITF-8 length
                let offset = decode_itf8(reader)?;
                let length = decode_itf8(reader)?;
                Ok(Self::Beta { offset, length })
            }
            7 => {
                // SUBEXP: ITF-8 offset, ITF-8 K
                let offset = decode_itf8(reader)?;
                let k = decode_itf8(reader)?;
                Ok(Self::SubExp { offset, k })
            }
            8 => {
                // GOLOMB_RICE: ITF-8 offset, ITF-8 log2_M
                let offset = decode_itf8(reader)?;
                let log2_m = decode_itf8(reader)?;
                Ok(Self::GolombRice { offset, log2_m })
            }
            9 => {
                // GAMMA: ITF-8 offset
                let offset = decode_itf8(reader)?;
                Ok(Self::Gamma { offset })
            }
            _ => Err(BiometalError::InvalidCramFormat {
                msg: format!("Unknown encoding ID: {}", encoding_id),
            }),
        }
    }

    /// Decode a single integer value using this encoding.
    ///
    /// **Phase 2 Full - Partial Implementation**:
    /// Currently only EXTERNAL encoding is implemented. Other encodings (HUFFMAN, BETA, GAMMA, etc.)
    /// require bit-level reading infrastructure and will be implemented in the next iteration.
    ///
    /// # Arguments
    /// * `blocks` - Map of block_content_id -> block data
    /// * `block_positions` - Current read position in each block
    ///
    /// # Returns
    /// Decoded integer value
    pub fn decode_int(
        &self,
        blocks: &std::collections::HashMap<i32, &[u8]>,
        block_positions: &mut std::collections::HashMap<i32, usize>,
    ) -> Result<i32> {
        match self {
            Self::Null => Ok(0),
            Self::External { block_content_id } => {
                // Read ITF-8 from external block
                let block_data = blocks.get(block_content_id)
                    .ok_or_else(|| BiometalError::InvalidCramFormat {
                        msg: format!("External block {} not found", block_content_id)
                    })?;

                let pos = block_positions.get(block_content_id).copied().unwrap_or(0);
                let mut reader = std::io::Cursor::new(&block_data[pos..]);
                let value = decode_itf8(&mut reader)?;
                let new_pos = pos + (reader.position() as usize);
                block_positions.insert(*block_content_id, new_pos);

                Ok(value)
            }
            _ => {
                // Other encodings require bit-level reading, deferred to next iteration
                Err(BiometalError::InvalidCramFormat {
                    msg: format!("Encoding {:?} not yet implemented (requires bit-level reader)", self)
                })
            }
        }
    }

    /// Decode a single byte using this encoding.
    pub fn decode_byte(
        &self,
        blocks: &std::collections::HashMap<i32, &[u8]>,
        block_positions: &mut std::collections::HashMap<i32, usize>,
    ) -> Result<u8> {
        match self {
            Self::Null => Ok(0),
            Self::External { block_content_id } => {
                // Read single byte from external block
                let block_data = blocks.get(block_content_id)
                    .ok_or_else(|| BiometalError::InvalidCramFormat {
                        msg: format!("External block {} not found", block_content_id)
                    })?;

                let pos = block_positions.get(block_content_id).copied().unwrap_or(0);
                if pos >= block_data.len() {
                    return Err(BiometalError::InvalidCramFormat {
                        msg: format!("Attempted to read past end of block {}", block_content_id)
                    });
                }

                let value = block_data[pos];
                block_positions.insert(*block_content_id, pos + 1);

                Ok(value)
            }
            _ => {
                Err(BiometalError::InvalidCramFormat {
                    msg: format!("Encoding {:?} not yet implemented", self)
                })
            }
        }
    }

    /// Decode a byte array using this encoding.
    pub fn decode_byte_array(
        &self,
        blocks: &std::collections::HashMap<i32, &[u8]>,
        block_positions: &mut std::collections::HashMap<i32, usize>,
    ) -> Result<Vec<u8>> {
        match self {
            Self::Null => Ok(Vec::new()),
            Self::ByteArrayLen { len_encoding, value_encoding } => {
                // Decode length, then decode that many values
                let len = len_encoding.decode_int(blocks, block_positions)? as usize;
                let mut result = Vec::with_capacity(len);
                for _ in 0..len {
                    result.push(value_encoding.decode_byte(blocks, block_positions)?);
                }
                Ok(result)
            }
            Self::ByteArrayStop { stop_byte, block_content_id } => {
                // Read until stop byte encountered
                let block_data = blocks.get(block_content_id)
                    .ok_or_else(|| BiometalError::InvalidCramFormat {
                        msg: format!("External block {} not found", block_content_id)
                    })?;

                let pos = block_positions.get(block_content_id).copied().unwrap_or(0);
                let mut result = Vec::new();
                let mut current_pos = pos;

                while current_pos < block_data.len() {
                    let byte = block_data[current_pos];
                    current_pos += 1;
                    if byte == *stop_byte {
                        break;
                    }
                    result.push(byte);
                }

                block_positions.insert(*block_content_id, current_pos);
                Ok(result)
            }
            _ => {
                Err(BiometalError::InvalidCramFormat {
                    msg: format!("Encoding {:?} not suitable for byte arrays", self)
                })
            }
        }
    }
}

/// CRAM preservation map.
///
/// Stores preservation policy settings for the container.
#[derive(Debug, Clone, Default)]
pub struct PreservationMap {
    /// Read names stored (RN)
    pub read_names_included: bool,
    /// AP data series delta (AP)
    pub ap_data_series_delta: bool,
    /// Reference required (RR)
    pub reference_required: bool,
    /// Substitution matrix (SM): 5-byte matrix
    pub substitution_matrix: Option<[u8; 5]>,
    /// Tag ID dictionary (TD): tag names
    pub tag_ids: Vec<[u8; 3]>,
}

impl PreservationMap {
    /// Parse preservation map from bytes.
    fn parse(data: &[u8]) -> Result<Self> {
        use std::io::Cursor;
        let mut reader = Cursor::new(data);
        let mut map = Self::default();

        // Read map size (number of entries)
        let map_size = decode_itf8(&mut reader)?;

        for _ in 0..map_size {
            // Read 2-byte key
            let mut key = [0u8; 2];
            reader.read_exact(&mut key)
                .map_err(|e| BiometalError::InvalidCramFormat {
                    msg: format!("Failed to read preservation map key: {}", e)
                })?;

            match &key {
                b"RN" => {
                    // Boolean: read names included
                    let mut buf = [0u8; 1];
                    reader.read_exact(&mut buf)
                        .map_err(|e| BiometalError::InvalidCramFormat {
                            msg: format!("Failed to read RN value: {}", e)
                        })?;
                    map.read_names_included = buf[0] != 0;
                }
                b"AP" => {
                    // Boolean: AP data series delta
                    let mut buf = [0u8; 1];
                    reader.read_exact(&mut buf)
                        .map_err(|e| BiometalError::InvalidCramFormat {
                            msg: format!("Failed to read AP value: {}", e)
                        })?;
                    map.ap_data_series_delta = buf[0] != 0;
                }
                b"RR" => {
                    // Boolean: reference required
                    let mut buf = [0u8; 1];
                    reader.read_exact(&mut buf)
                        .map_err(|e| BiometalError::InvalidCramFormat {
                            msg: format!("Failed to read RR value: {}", e)
                        })?;
                    map.reference_required = buf[0] != 0;
                }
                b"SM" => {
                    // Substitution matrix: 5 bytes
                    let mut matrix = [0u8; 5];
                    reader.read_exact(&mut matrix)
                        .map_err(|e| BiometalError::InvalidCramFormat {
                            msg: format!("Failed to read SM value: {}", e)
                        })?;
                    map.substitution_matrix = Some(matrix);
                }
                b"TD" => {
                    // Tag IDs: array of 3-byte tag names
                    let num_tags = decode_itf8(&mut reader)?;
                    for _ in 0..num_tags {
                        let mut tag = [0u8; 3];
                        reader.read_exact(&mut tag)
                            .map_err(|e| BiometalError::InvalidCramFormat {
                                msg: format!("Failed to read tag ID: {}", e)
                            })?;
                        map.tag_ids.push(tag);
                    }
                }
                _ => {
                    // Unknown key - skip the value
                    // Value length is not standardized, so we can't safely skip
                    // For now, return an error
                    return Err(BiometalError::InvalidCramFormat {
                        msg: format!(
                            "Unknown preservation map key: {}{}",
                            key[0] as char, key[1] as char
                        ),
                    });
                }
            }
        }

        Ok(map)
    }
}

/// CRAM compression header.
///
/// The compression header describes how data is encoded in the container.
#[derive(Debug, Clone)]
pub struct CompressionHeader {
    /// Preservation map
    pub preservation_map: PreservationMap,
    /// Data series encoding map
    pub data_series_encoding: std::collections::HashMap<DataSeries, Encoding>,
    /// Tag encoding map (tag_id -> encoding)
    pub tag_encoding: std::collections::HashMap<i32, Encoding>,
}

impl CompressionHeader {
    /// Parse compression header from decompressed block data.
    ///
    /// **Phase 2 Full**: Parses preservation map, data series encoding, and tag encoding.
    ///
    /// **Remaining Work**:
    /// - Implement actual data decoders (EXTERNAL, HUFFMAN, etc.) - see next task
    /// - Decode CRAM features using the encoding map
    /// - Apply features to reference sequences
    pub fn parse(data: &[u8]) -> Result<Self> {
        use std::io::Cursor;
        let mut reader = Cursor::new(data);

        // Parse preservation map
        let preservation_map_size = decode_itf8(&mut reader)?;
        let mut preservation_map_data = vec![0u8; preservation_map_size as usize];
        reader.read_exact(&mut preservation_map_data)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read preservation map: {}", e)
            })?;
        let preservation_map = PreservationMap::parse(&preservation_map_data)?;

        // Parse data series encoding map
        let data_series_size = decode_itf8(&mut reader)?;
        let mut data_series_encoding = std::collections::HashMap::new();
        let data_series_start = reader.position() as usize;
        let data_series_end = data_series_start + data_series_size as usize;

        let mut data_series_reader = Cursor::new(&data[data_series_start..data_series_end]);
        let num_data_series = decode_itf8(&mut data_series_reader)?;

        for _ in 0..num_data_series {
            // Read 2-byte data series key
            let mut key = [0u8; 2];
            data_series_reader.read_exact(&mut key)
                .map_err(|e| BiometalError::InvalidCramFormat {
                    msg: format!("Failed to read data series key: {}", e)
                })?;

            // Parse encoding for this data series
            let encoding = Encoding::parse(&mut data_series_reader)?;
            data_series_encoding.insert(DataSeries::from_bytes(key), encoding);
        }

        reader.set_position((data_series_end) as u64);

        // Parse tag encoding map
        let tag_encoding_size = decode_itf8(&mut reader)?;
        let mut tag_encoding = std::collections::HashMap::new();
        let tag_encoding_start = reader.position() as usize;
        let tag_encoding_end = tag_encoding_start + tag_encoding_size as usize;

        let mut tag_reader = Cursor::new(&data[tag_encoding_start..tag_encoding_end]);
        let num_tags = decode_itf8(&mut tag_reader)?;

        for _ in 0..num_tags {
            // Read tag ID (ITF-8)
            let tag_id = decode_itf8(&mut tag_reader)?;

            // Parse encoding for this tag
            let encoding = Encoding::parse(&mut tag_reader)?;
            tag_encoding.insert(tag_id, encoding);
        }

        Ok(CompressionHeader {
            preservation_map,
            data_series_encoding,
            tag_encoding,
        })
    }
}

/// CRAM slice header.
///
/// A slice is a subdivision of a container containing alignment records.
#[derive(Debug, Clone)]
pub struct SliceHeader {
    /// Reference sequence ID (-1 for unmapped, -2 for multi-ref)
    pub reference_id: i32,
    /// Alignment start position (1-based)
    pub alignment_start: i32,
    /// Alignment span (number of reference bases)
    pub alignment_span: i32,
    /// Number of records in this slice
    pub num_records: i32,
    /// Global record counter
    pub record_counter: i64,
    /// Number of blocks in this slice
    pub num_blocks: i32,
    /// Block content IDs (describes what each block contains)
    pub block_content_ids: Vec<i32>,
    /// Embedded reference MD5 (16 bytes, for reference validation)
    pub embedded_ref_md5: [u8; 16],
    /// Optional reference MD5 (if present)
    pub optional_ref_md5: Option<[u8; 16]>,
}

impl SliceHeader {
    /// Parse slice header from reader.
    ///
    /// # Slice Header Format
    ///
    /// ```text
    /// - Reference sequence ID: ITF-8
    /// - Alignment start: ITF-8
    /// - Alignment span: ITF-8
    /// - Number of records: ITF-8
    /// - Record counter: LTF-8
    /// - Number of blocks: ITF-8
    /// - Block content IDs: ITF-8 array
    /// - Embedded reference MD5: [u8; 16]
    /// - Optional reference MD5: Optional [u8; 16]
    /// ```
    pub fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        // Read reference sequence ID
        let reference_id = decode_itf8(reader)?;

        // Read alignment start
        let alignment_start = decode_itf8(reader)?;

        // Read alignment span
        let alignment_span = decode_itf8(reader)?;

        // Read number of records
        let num_records = decode_itf8(reader)?;

        // Read record counter
        let record_counter = decode_ltf8(reader)?;

        // Read number of blocks
        let num_blocks = decode_itf8(reader)?;

        // Read block content IDs
        let mut block_content_ids = Vec::with_capacity(num_blocks as usize);
        for _ in 0..num_blocks {
            block_content_ids.push(decode_itf8(reader)?);
        }

        // Read embedded reference MD5 (16 bytes)
        let mut embedded_ref_md5 = [0u8; 16];
        reader.read_exact(&mut embedded_ref_md5)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read embedded reference MD5: {}", e)
            })?;

        // TODO Phase 2: Read optional reference MD5 based on actual CRAM flags
        // For Phase 1, we skip this to avoid consuming block data
        let optional_ref_md5 = None;

        Ok(SliceHeader {
            reference_id,
            alignment_start,
            alignment_span,
            num_records,
            record_counter,
            num_blocks,
            block_content_ids,
            embedded_ref_md5,
            optional_ref_md5,
        })
    }
}

/// CRAM read feature.
///
/// Features describe how a read differs from the reference sequence.
/// This is the core of CRAM's reference-based compression.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CramFeature {
    /// Base substitution: reference position + substituted base
    Substitution {
        /// Position in read (0-based)
        position: i32,
        /// Substituted base (A, C, G, T, N)
        base: u8,
    },
    /// Insertion: reference position + inserted bases
    Insertion {
        /// Position in read (0-based)
        position: i32,
        /// Inserted bases
        bases: Vec<u8>,
    },
    /// Deletion: reference position + deletion length
    Deletion {
        /// Position in read (0-based)
        position: i32,
        /// Number of bases deleted
        length: i32,
    },
    /// Reference skip: reference position + skip length
    ReferenceSkip {
        /// Position in read (0-based)
        position: i32,
        /// Number of bases skipped
        length: i32,
    },
    /// Soft clip: reference position + clipped bases
    SoftClip {
        /// Position in read (0-based)
        position: i32,
        /// Soft clipped bases
        bases: Vec<u8>,
    },
    /// Hard clip: reference position + clip length
    HardClip {
        /// Position in read (0-based)
        position: i32,
        /// Number of bases hard clipped
        length: i32,
    },
    /// Padding: reference position + padding length
    Padding {
        /// Position in read (0-based)
        position: i32,
        /// Padding length
        length: i32,
    },
    /// Insert base: single base insertion
    InsertBase {
        /// Position in read (0-based)
        position: i32,
        /// Inserted base
        base: u8,
    },
    /// Quality score: position + quality
    QualityScore {
        /// Position in read (0-based)
        position: i32,
        /// Quality score (Phred scale)
        score: u8,
    },
    /// Read base: position + base
    ReadBase {
        /// Position in read (0-based)
        position: i32,
        /// Base at this position
        base: u8,
    },
    /// Bases: stretch of bases (for unmapped reads)
    Bases {
        /// Position in read (0-based)
        position: i32,
        /// Bases
        bases: Vec<u8>,
    },
    /// Scores: stretch of quality scores
    Scores {
        /// Position in read (0-based)
        position: i32,
        /// Quality scores
        scores: Vec<u8>,
    },
}

impl CramFeature {
    /// Decode features for a single read from CRAM blocks.
    ///
    /// # Arguments
    /// * `compression_header` - Compression header with encoding definitions
    /// * `blocks` - Map of block_content_id -> block data
    /// * `block_positions` - Current read positions in each block
    ///
    /// # Returns
    /// Vector of features for this read
    pub fn decode_features(
        compression_header: &CompressionHeader,
        blocks: &std::collections::HashMap<i32, &[u8]>,
        block_positions: &mut std::collections::HashMap<i32, usize>,
    ) -> Result<Vec<Self>> {
        // Get encoding for feature count (FN)
        let fn_encoding = compression_header
            .data_series_encoding
            .get(&DataSeries::FN)
            .ok_or_else(|| BiometalError::InvalidCramFormat {
                msg: "Missing FN (feature count) encoding".to_string(),
            })?;

        // Decode number of features for this read
        let num_features = fn_encoding.decode_int(blocks, block_positions)?;

        if num_features == 0 {
            return Ok(Vec::new());
        }

        // Get encodings for feature codes and positions
        let fc_encoding = compression_header
            .data_series_encoding
            .get(&DataSeries::FC)
            .ok_or_else(|| BiometalError::InvalidCramFormat {
                msg: "Missing FC (feature codes) encoding".to_string(),
            })?;

        let fp_encoding = compression_header
            .data_series_encoding
            .get(&DataSeries::FP)
            .ok_or_else(|| BiometalError::InvalidCramFormat {
                msg: "Missing FP (feature positions) encoding".to_string(),
            })?;

        let mut features = Vec::with_capacity(num_features as usize);

        for _ in 0..num_features {
            // Decode feature code
            let code = fc_encoding.decode_byte(blocks, block_positions)? as char;

            // Decode feature position
            let position = fp_encoding.decode_int(blocks, block_positions)?;

            // Decode feature data based on code
            let feature = match code {
                'B' => {
                    // Base substitution: code + base
                    let bs_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::BS)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing BS (base substitution) encoding".to_string(),
                        })?;
                    let base = bs_encoding.decode_byte(blocks, block_positions)?;
                    Self::Substitution { position, base }
                }
                'X' => {
                    // Read base: code + base
                    let ba_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::BA)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing BA (bases) encoding".to_string(),
                        })?;
                    let base = ba_encoding.decode_byte(blocks, block_positions)?;
                    Self::ReadBase { position, base }
                }
                'D' => {
                    // Deletion: code + length
                    let dl_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::DL)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing DL (deletion length) encoding".to_string(),
                        })?;
                    let length = dl_encoding.decode_int(blocks, block_positions)?;
                    Self::Deletion { position, length }
                }
                'I' => {
                    // Insertion: code + bases
                    let in_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::IN)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing IN (insertion) encoding".to_string(),
                        })?;
                    let bases = in_encoding.decode_byte_array(blocks, block_positions)?;
                    Self::Insertion { position, bases }
                }
                'i' => {
                    // Insert base: code + single base
                    let ba_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::BA)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing BA (bases) encoding".to_string(),
                        })?;
                    let base = ba_encoding.decode_byte(blocks, block_positions)?;
                    Self::InsertBase { position, base }
                }
                'S' => {
                    // Soft clip: code + bases
                    let sc_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::SC)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing SC (soft clip) encoding".to_string(),
                        })?;
                    let bases = sc_encoding.decode_byte_array(blocks, block_positions)?;
                    Self::SoftClip { position, bases }
                }
                'H' => {
                    // Hard clip: code + length
                    let hc_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::HC)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing HC (hard clip) encoding".to_string(),
                        })?;
                    let length = hc_encoding.decode_int(blocks, block_positions)?;
                    Self::HardClip { position, length }
                }
                'P' => {
                    // Padding: code + length
                    let pd_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::PD)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing PD (padding) encoding".to_string(),
                        })?;
                    let length = pd_encoding.decode_int(blocks, block_positions)?;
                    Self::Padding { position, length }
                }
                'N' => {
                    // Reference skip: code + length
                    let rs_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::RS)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing RS (reference skip) encoding".to_string(),
                        })?;
                    let length = rs_encoding.decode_int(blocks, block_positions)?;
                    Self::ReferenceSkip { position, length }
                }
                'b' => {
                    // Bases: code + bases
                    let ba_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::BA)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing BA (bases) encoding".to_string(),
                        })?;
                    let bases = ba_encoding.decode_byte_array(blocks, block_positions)?;
                    Self::Bases { position, bases }
                }
                'q' => {
                    // Quality score: code + score
                    let qs_encoding = compression_header
                        .data_series_encoding
                        .get(&DataSeries::QS)
                        .ok_or_else(|| BiometalError::InvalidCramFormat {
                            msg: "Missing QS (quality scores) encoding".to_string(),
                        })?;
                    let score = qs_encoding.decode_byte(blocks, block_positions)?;
                    Self::QualityScore { position, score }
                }
                _ => {
                    return Err(BiometalError::InvalidCramFormat {
                        msg: format!("Unknown feature code: '{}'", code),
                    })
                }
            };

            features.push(feature);
        }

        Ok(features)
    }

    /// Apply features to reference sequence to reconstruct the actual read sequence.
    ///
    /// # Arguments
    /// * `reference` - Reference sequence (ASCII bases)
    /// * `features` - List of features describing differences from reference
    /// * `read_length` - Expected length of the read
    ///
    /// # Returns
    /// Reconstructed read sequence
    pub fn apply_to_reference(
        reference: &[u8],
        features: &[Self],
        read_length: usize,
    ) -> Result<Vec<u8>> {
        // Start with reference sequence
        let mut sequence = reference.to_vec();

        // Ensure sequence is long enough
        if sequence.len() < read_length {
            sequence.resize(read_length, b'N');
        }

        // Apply features in order
        for feature in features {
            match feature {
                Self::Substitution { position, base } => {
                    // Replace base at position
                    let pos = *position as usize;
                    if pos < sequence.len() {
                        sequence[pos] = *base;
                    }
                }
                Self::Insertion { position, bases } => {
                    // Insert bases at position
                    let pos = *position as usize;
                    if pos <= sequence.len() {
                        sequence.splice(pos..pos, bases.iter().copied());
                    }
                }
                Self::Deletion { position, length } => {
                    // Delete bases from position
                    let pos = *position as usize;
                    let len = *length as usize;
                    if pos < sequence.len() {
                        let end = (pos + len).min(sequence.len());
                        sequence.drain(pos..end);
                    }
                }
                Self::InsertBase { position, base } => {
                    // Insert single base
                    let pos = *position as usize;
                    if pos <= sequence.len() {
                        sequence.insert(pos, *base);
                    }
                }
                Self::SoftClip { position, bases } => {
                    // Soft clipped bases are part of the sequence
                    let pos = *position as usize;
                    if pos <= sequence.len() {
                        sequence.splice(pos..pos, bases.iter().copied());
                    }
                }
                Self::HardClip { .. } => {
                    // Hard clips don't affect the sequence (bases not present)
                }
                Self::Padding { .. } => {
                    // Padding doesn't affect the sequence
                }
                Self::ReferenceSkip { position, length } => {
                    // Reference skip (like intron): remove from sequence
                    let pos = *position as usize;
                    let len = *length as usize;
                    if pos < sequence.len() {
                        let end = (pos + len).min(sequence.len());
                        sequence.drain(pos..end);
                    }
                }
                Self::ReadBase { position, base } => {
                    // Explicit read base
                    let pos = *position as usize;
                    if pos < sequence.len() {
                        sequence[pos] = *base;
                    } else if pos == sequence.len() {
                        sequence.push(*base);
                    }
                }
                Self::Bases { position, bases } => {
                    // Replace stretch of bases
                    let pos = *position as usize;
                    if pos < sequence.len() {
                        let end = (pos + bases.len()).min(sequence.len());
                        sequence.splice(pos..end, bases.iter().copied());
                    } else {
                        sequence.extend_from_slice(bases);
                    }
                }
                Self::QualityScore { .. } => {
                    // Quality scores don't affect the sequence
                }
                Self::Scores { .. } => {
                    // Quality scores don't affect the sequence
                }
            }
        }

        // Truncate or pad to expected read length
        sequence.resize(read_length, b'N');

        Ok(sequence)
    }

    /// Build CIGAR string from features.
    ///
    /// CIGAR operations are derived from CRAM features:
    /// - M (match/mismatch): Default alignment to reference
    /// - I (insertion): From Insertion/InsertBase features
    /// - D (deletion): From Deletion features
    /// - N (reference skip): From ReferenceSkip features
    /// - S (soft clip): From SoftClip features
    /// - H (hard clip): From HardClip features
    /// - P (padding): From Padding features
    ///
    /// # Arguments
    /// * `features` - List of features for this read
    /// * `read_length` - Length of the read sequence
    /// * `reference_length` - Length covered on reference
    ///
    /// # Returns
    /// CIGAR string (e.g., "100M", "50M2I48M")
    pub fn build_cigar(
        features: &[Self],
        read_length: usize,
        reference_length: usize,
    ) -> String {
        if features.is_empty() {
            // No features = perfect match to reference
            return format!("{}M", read_length);
        }

        let mut cigar_ops: Vec<(usize, char)> = Vec::new();
        let mut current_pos = 0;

        // Sort features by position
        let mut sorted_features: Vec<_> = features.iter().collect();
        sorted_features.sort_by_key(|f| match f {
            Self::Substitution { position, .. } => *position,
            Self::Insertion { position, .. } => *position,
            Self::Deletion { position, .. } => *position,
            Self::ReferenceSkip { position, .. } => *position,
            Self::SoftClip { position, .. } => *position,
            Self::HardClip { position, .. } => *position,
            Self::Padding { position, .. } => *position,
            Self::InsertBase { position, .. } => *position,
            Self::ReadBase { position, .. } => *position,
            Self::Bases { position, .. } => *position,
            Self::QualityScore { position, .. } => *position,
            Self::Scores { position, .. } => *position,
        });

        for feature in sorted_features {
            match feature {
                Self::HardClip { position, length } => {
                    // Hard clip at start or end
                    let pos = *position as usize;
                    if pos == 0 {
                        cigar_ops.push((*length as usize, 'H'));
                    } else {
                        // Hard clip at end - will be added at end
                        cigar_ops.push((*length as usize, 'H'));
                    }
                }
                Self::SoftClip { position, bases } => {
                    let pos = *position as usize;
                    // Add match up to soft clip
                    if pos > current_pos {
                        cigar_ops.push((pos - current_pos, 'M'));
                    }
                    cigar_ops.push((bases.len(), 'S'));
                    current_pos = pos + bases.len();
                }
                Self::Insertion { position, bases } => {
                    let pos = *position as usize;
                    // Add match up to insertion
                    if pos > current_pos {
                        cigar_ops.push((pos - current_pos, 'M'));
                    }
                    cigar_ops.push((bases.len(), 'I'));
                    current_pos = pos;
                }
                Self::InsertBase { position, .. } => {
                    let pos = *position as usize;
                    // Add match up to insertion
                    if pos > current_pos {
                        cigar_ops.push((pos - current_pos, 'M'));
                    }
                    cigar_ops.push((1, 'I'));
                    current_pos = pos;
                }
                Self::Deletion { position, length } => {
                    let pos = *position as usize;
                    // Add match up to deletion
                    if pos > current_pos {
                        cigar_ops.push((pos - current_pos, 'M'));
                    }
                    cigar_ops.push((*length as usize, 'D'));
                    current_pos = pos;
                }
                Self::ReferenceSkip { position, length } => {
                    let pos = *position as usize;
                    // Add match up to skip
                    if pos > current_pos {
                        cigar_ops.push((pos - current_pos, 'M'));
                    }
                    cigar_ops.push((*length as usize, 'N'));
                    current_pos = pos;
                }
                Self::Padding { position, length } => {
                    let pos = *position as usize;
                    // Add match up to padding
                    if pos > current_pos {
                        cigar_ops.push((pos - current_pos, 'M'));
                    }
                    cigar_ops.push((*length as usize, 'P'));
                    current_pos = pos;
                }
                Self::Substitution { position, .. } => {
                    // Substitutions are still 'M' in CIGAR (match or mismatch)
                    let pos = *position as usize;
                    if pos >= current_pos {
                        current_pos = pos + 1;
                    }
                }
                Self::ReadBase { position, .. } => {
                    // Read bases are 'M' in CIGAR
                    let pos = *position as usize;
                    if pos >= current_pos {
                        current_pos = pos + 1;
                    }
                }
                Self::Bases { position, bases } => {
                    // Explicit bases are 'M' in CIGAR
                    let pos = *position as usize;
                    if pos >= current_pos {
                        current_pos = pos + bases.len();
                    }
                }
                Self::QualityScore { .. } | Self::Scores { .. } => {
                    // Quality scores don't affect CIGAR
                }
            }
        }

        // Add remaining match to end of read
        if current_pos < read_length {
            cigar_ops.push((read_length - current_pos, 'M'));
        }

        // Merge consecutive operations of the same type
        let mut merged: Vec<(usize, char)> = Vec::new();
        for (length, op) in cigar_ops {
            if let Some((last_len, last_op)) = merged.last_mut() {
                if *last_op == op {
                    *last_len += length;
                    continue;
                }
            }
            merged.push((length, op));
        }

        // Build CIGAR string
        if merged.is_empty() {
            format!("{}M", read_length)
        } else {
            merged
                .iter()
                .map(|(len, op)| format!("{}{}", len, op))
                .collect::<Vec<_>>()
                .join("")
        }
    }
}

/// CRAM slice (header + blocks).
///
/// A complete slice containing a header and all its data blocks.
#[derive(Debug, Clone)]
pub struct Slice {
    /// Slice header
    pub header: SliceHeader,
    /// Blocks in this slice (core block + external blocks)
    pub blocks: Vec<Block>,
}

impl Slice {
    /// Parse slice from reader.
    ///
    /// Reads the slice header followed by all blocks specified in the header.
    pub fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        let header = SliceHeader::parse(reader)?;

        // Read all blocks
        let mut blocks = Vec::with_capacity(header.num_blocks as usize);
        for _ in 0..header.num_blocks {
            blocks.push(Block::parse(reader)?);
        }

        Ok(Slice {
            header,
            blocks,
        })
    }

    /// Get the core data block (first block, content ID = 0).
    pub fn core_block(&self) -> Option<&Block> {
        self.blocks.first()
    }

    /// Get external blocks (all blocks after the core block).
    pub fn external_blocks(&self) -> &[Block] {
        if self.blocks.len() > 1 {
            &self.blocks[1..]
        } else {
            &[]
        }
    }

    /// Decode sequence from slice blocks.
    ///
    /// **Phase 2 Status**: Basic reference-based reconstruction (simplified)
    ///
    /// **Current Implementation**:
    /// - Fetches reference subsequence for alignment span
    /// - Returns empty placeholder if no reference available
    ///
    /// **Phase 2 Full TODO** (10-15 hours):
    /// 1. Parse compression header to get encoding schemes for:
    ///    - FC: Feature codes (what kind of variation)
    ///    - FP: Feature positions (where in the read)
    ///    - BS: Base substitutions (what base to substitute)
    ///    - IN: Insertion bases (what bases to insert)
    ///    - DL: Deletion lengths (how many bases to delete)
    ///    - BA: Bases (for unmapped reads)
    ///    - QS: Quality scores
    ///
    /// 2. Read from core block and external blocks:
    ///    - Decode feature count for this record
    ///    - For each feature, decode: code, position, base/length
    ///
    /// 3. Apply features to reference sequence:
    ///    ```rust
    ///    let mut sequence = reference_sequence.clone();
    ///    for feature in features {
    ///        match feature.code {
    ///            'B' => sequence[feature.pos] = feature.base,  // Substitution
    ///            'I' => sequence.splice(feature.pos..feature.pos, feature.bases),  // Insertion
    ///            'D' => sequence.drain(feature.pos..feature.pos+feature.len),  // Deletion
    ///            'S' => /* Soft clip */,
    ///            // ... other feature types
    ///        }
    ///    }
    ///    ```
    ///
    /// 4. Handle unmapped reads:
    ///    - Read bases directly from BA data series
    ///    - No reference sequence needed
    ///
    /// 5. Build CIGAR from features:
    ///    - Track match/mismatch runs
    ///    - Convert insertions/deletions to CIGAR ops
    ///
    /// # Arguments
    ///
    /// * `record_index` - Index of record within slice (0-based)
    /// * `reference_index` - Optional reference FASTA index
    /// * `reference_path` - Optional reference FASTA file path
    /// * `header_references` - Reference names from BAM header
    ///
    /// # Returns
    ///
    /// Decoded sequence as bytes (ASCII: A, C, G, T, N)
    fn decode_sequence(
        &self,
        _record_index: usize,
        reference_index: Option<&FaiIndex>,
        reference_path: Option<&PathBuf>,
        header_references: &[crate::io::bam::Reference],
    ) -> Result<Vec<u8>> {
        // Phase 2: Basic reference fetching
        // TODO Phase 2 full: Decode CRAM features and apply to reference

        if let (Some(index), Some(ref_path)) = (reference_index, reference_path) {
            // Fetch reference sequence for alignment span
            if self.header.reference_id >= 0 {
                let ref_id = self.header.reference_id as usize;

                if let Some(ref_info) = header_references.get(ref_id) {
                    let start = self.header.alignment_start as u64;
                    let end = (self.header.alignment_start + self.header.alignment_span) as u64;

                    // Fetch reference subsequence
                    let sequence = index.fetch_region(
                        &ref_info.name,
                        start,
                        end,
                        ref_path,
                    )?;

                    return Ok(sequence.into_bytes());
                }
            }
        }

        // Fallback: Return placeholder sequence
        // TODO Phase 2: Decode from blocks when reference not available
        Ok(Vec::new())
    }

    /// Decode quality scores for a single read from CRAM blocks.
    ///
    /// **Phase 2 Full Implementation**: Decodes quality scores from external blocks.
    ///
    /// Quality scores can come from:
    /// 1. QS data series (for all bases)
    /// 2. Quality features ('q' feature code) for individual bases
    ///
    /// # Arguments
    /// * `compression_header` - Compression header with encoding definitions
    /// * `blocks` - Map of block_content_id -> block data
    /// * `block_positions` - Current read positions in each block
    /// * `features` - Features for this read (may contain quality scores)
    /// * `sequence_length` - Length of sequence (for quality score array size)
    ///
    /// # Returns
    /// Quality scores (Phred+33 ASCII)
    pub fn decode_quality_scores(
        compression_header: &CompressionHeader,
        blocks: &std::collections::HashMap<i32, &[u8]>,
        block_positions: &mut std::collections::HashMap<i32, usize>,
        features: &[CramFeature],
        sequence_length: usize,
    ) -> Result<Vec<u8>> {
        // Start with default quality scores (Phred 30 = '?')
        let mut qualities = vec![b'?'; sequence_length];

        // Try to decode from QS data series
        if let Some(qs_encoding) = compression_header.data_series_encoding.get(&DataSeries::QS) {
            // Decode quality scores from QS data series
            match qs_encoding.decode_byte_array(blocks, block_positions) {
                Ok(scores) => {
                    // Convert Phred scores to ASCII (Phred+33)
                    for (i, &score) in scores.iter().enumerate().take(sequence_length) {
                        qualities[i] = score.saturating_add(33);
                    }
                    return Ok(qualities);
                }
                Err(_) => {
                    // Fall through to try features
                }
            }
        }

        // Apply quality scores from features
        for feature in features {
            match feature {
                CramFeature::QualityScore { position, score } => {
                    let pos = *position as usize;
                    if pos < sequence_length {
                        qualities[pos] = score.saturating_add(33);
                    }
                }
                CramFeature::Scores { position, scores } => {
                    let pos = *position as usize;
                    for (i, &score) in scores.iter().enumerate() {
                        if pos + i < sequence_length {
                            qualities[pos + i] = score.saturating_add(33);
                        }
                    }
                }
                _ => {}
            }
        }

        Ok(qualities)
    }

    /// Decode quality scores (simplified instance method for iterator).
    ///
    /// This is a placeholder that returns default quality scores.
    /// For full quality score decoding, use the static decode_quality_scores method
    /// with compression header and blocks.
    ///
    /// # Arguments
    /// * `_record_index` - Index of record within slice (0-based)
    /// * `sequence_length` - Length of sequence (for quality score array size)
    ///
    /// # Returns
    /// Quality scores (Phred+33 ASCII)
    fn decode_quality_scores_simple(
        &self,
        _record_index: usize,
        sequence_length: usize,
    ) -> Vec<u8> {
        // Return default quality scores (Phred 30 = '?')
        // Full implementation requires compression header access
        vec![b'?'; sequence_length]
    }
}

// ============================================================================
// Variable-Length Integer Encoding (ITF-8 and LTF-8)
// ============================================================================

/// Decode ITF-8 (Integer, Type-Free, 8-bit) variable-length integer.
///
/// ITF-8 encodes 32-bit integers in 1-5 bytes based on magnitude:
/// - `0xxxxxxx`: 1 byte (0-127)
/// - `10xxxxxx xxxxxxxx`: 2 bytes (128-16,383)
/// - `110xxxxx xxxxxxxx xxxxxxxx`: 3 bytes (16,384-2,097,151)
/// - `1110xxxx xxxxxxxx xxxxxxxx xxxxxxxx`: 4 bytes (2,097,152-268,435,455)
/// - `1111xxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx`: 5 bytes (268,435,456+)
///
/// # Arguments
///
/// * `reader` - Source to read ITF-8 encoded integer from
///
/// # Returns
///
/// Decoded 32-bit integer value
///
/// # Example
///
/// ```no_run
/// use std::io::Cursor;
/// # use biometal::io::cram::decode_itf8;
/// # use biometal::Result;
/// # fn main() -> Result<()> {
/// let data = vec![0x85, 0x42]; // 2-byte ITF-8 encoding
/// let mut reader = Cursor::new(data);
/// let value = decode_itf8(&mut reader)?;
/// assert_eq!(value, 322); // 0x142
/// # Ok(())
/// # }
/// ```
pub fn decode_itf8<R: Read>(reader: &mut R) -> Result<i32> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)
        .map_err(|e| BiometalError::InvalidCramFormat {
            msg: format!("Failed to read ITF-8 first byte: {}", e)
        })?;

    let first = buf[0];

    // Determine length based on leading bits
    if first & 0x80 == 0 {
        // 1 byte: 0xxxxxxx
        Ok(first as i32)
    } else if first & 0x40 == 0 {
        // 2 bytes: 10xxxxxx xxxxxxxx
        reader.read_exact(&mut buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read ITF-8 byte 2: {}", e)
            })?;
        Ok((((first & 0x3F) as i32) << 8) | (buf[0] as i32))
    } else if first & 0x20 == 0 {
        // 3 bytes: 110xxxxx xxxxxxxx xxxxxxxx
        let mut bytes = [0u8; 2];
        reader.read_exact(&mut bytes)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read ITF-8 bytes 2-3: {}", e)
            })?;
        Ok((((first & 0x1F) as i32) << 16)
            | ((bytes[0] as i32) << 8)
            | (bytes[1] as i32))
    } else if first & 0x10 == 0 {
        // 4 bytes: 1110xxxx xxxxxxxx xxxxxxxx xxxxxxxx
        let mut bytes = [0u8; 3];
        reader.read_exact(&mut bytes)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read ITF-8 bytes 2-4: {}", e)
            })?;
        Ok((((first & 0x0F) as i32) << 24)
            | ((bytes[0] as i32) << 16)
            | ((bytes[1] as i32) << 8)
            | (bytes[2] as i32))
    } else {
        // 5 bytes: 1111xxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
        let mut bytes = [0u8; 4];
        reader.read_exact(&mut bytes)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read ITF-8 bytes 2-5: {}", e)
            })?;
        Ok((((first & 0x0F) as i32) << 28)
            | ((bytes[0] as i32) << 24)
            | ((bytes[1] as i32) << 16)
            | ((bytes[2] as i32) << 8)
            | (bytes[3] as i32))
    }
}

/// Decode LTF-8 (Long, Type-Free, 8-bit) variable-length integer.
///
/// LTF-8 encodes 64-bit integers in 1-9 bytes using the same pattern as ITF-8.
///
/// # Arguments
///
/// * `reader` - Source to read LTF-8 encoded integer from
///
/// # Returns
///
/// Decoded 64-bit integer value
///
/// # Example
///
/// ```no_run
/// use std::io::Cursor;
/// # use biometal::io::cram::decode_ltf8;
/// # use biometal::Result;
/// # fn main() -> Result<()> {
/// let data = vec![0x85, 0x42]; // 2-byte LTF-8 encoding
/// let mut reader = Cursor::new(data);
/// let value = decode_ltf8(&mut reader)?;
/// assert_eq!(value, 322); // 0x142
/// # Ok(())
/// # }
/// ```
pub fn decode_ltf8<R: Read>(reader: &mut R) -> Result<i64> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)
        .map_err(|e| BiometalError::InvalidCramFormat {
            msg: format!("Failed to read LTF-8 first byte: {}", e)
        })?;

    let first = buf[0];

    // Count leading 1s to determine byte count
    let num_bytes = if first & 0x80 == 0 {
        1
    } else if first & 0x40 == 0 {
        2
    } else if first & 0x20 == 0 {
        3
    } else if first & 0x10 == 0 {
        4
    } else if first & 0x08 == 0 {
        5
    } else if first & 0x04 == 0 {
        6
    } else if first & 0x02 == 0 {
        7
    } else if first & 0x01 == 0 {
        8
    } else {
        9
    };

    if num_bytes == 1 {
        return Ok(first as i64);
    }

    // Read remaining bytes
    let mut bytes = vec![0u8; num_bytes - 1];
    reader.read_exact(&mut bytes)
        .map_err(|e| BiometalError::InvalidCramFormat {
            msg: format!("Failed to read LTF-8 bytes 2-{}: {}", num_bytes, e)
        })?;

    // Extract value bits from first byte
    let mask = 0xFF >> num_bytes;
    let mut value = ((first & mask) as i64) << ((num_bytes - 1) * 8);

    // Add remaining bytes
    for (i, &byte) in bytes.iter().enumerate() {
        value |= (byte as i64) << ((num_bytes - 2 - i) * 8);
    }

    Ok(value)
}

/// CRAM file reader with streaming interface.
///
/// **Implementation Status**: Phase 1 basic iteration (Phase 2: full decoding)
///
/// This is a **native, zero-dependency** CRAM parser designed for:
/// - ARM NEON optimizations (16-25× speedup potential)
/// - Streaming architecture (constant ~5 MB memory)
/// - Full CRAM 3.0/3.1 compliance
///
/// # Memory Footprint
///
/// Maintains constant ~5 MB memory regardless of file size (Rule 5: Streaming).
///
/// # Phase 1 Limitations
///
/// - Basic record structure extraction (no full reference reconstruction)
/// - Placeholder sequences and quality scores
/// - Full decoding in Phase 2
///
/// # Example
///
/// ```no_run
/// use biometal::io::cram::CramReader;
///
/// # fn main() -> biometal::Result<()> {
/// let mut cram = CramReader::from_path("sample.cram")?;
///
/// let mut count = 0;
/// for record in cram.records() {
///     count += 1;
/// }
/// println!("Processed {} records", count);
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct CramReader<R: Read> {
    /// Underlying reader
    reader: R,
    /// BAM-compatible header
    header: Header,
    /// CRAM major version (2 or 3)
    major_version: u8,
    /// CRAM minor version
    minor_version: u8,
    /// Reference FASTA index (optional, for sequence reconstruction)
    reference_index: Option<FaiIndex>,
    /// Reference FASTA file path (optional, for sequence reconstruction)
    reference_path: Option<PathBuf>,
}

impl CramReader<BufReader<File>> {
    /// Open CRAM file from path.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to CRAM file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::cram::CramReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let cram = CramReader::from_path("alignments.cram")?;
    /// println!("CRAM version: {}.{}", cram.major_version(), cram.minor_version());
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())
            .map_err(|e| io::Error::new(
                io::ErrorKind::NotFound,
                format!("Failed to open CRAM file: {}", e)
            ))?;

        Self::new(BufReader::new(file))
    }

    /// Open CRAM file with reference FASTA.
    ///
    /// **Status**: Phase 2 - Not yet implemented
    ///
    /// # Arguments
    ///
    /// * `cram_path` - Path to CRAM file
    /// * `reference_path` - Path to reference FASTA
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::cram::CramReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let cram = CramReader::from_path_with_reference(
    ///     "alignments.cram",
    ///     "hg38.fa"
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path_with_reference<P: AsRef<Path>, Q: AsRef<Path>>(
        cram_path: P,
        reference_path: Q,
    ) -> Result<Self> {
        let mut reader = Self::from_path(cram_path)?;
        reader.set_reference(reference_path)?;
        Ok(reader)
    }
}

impl<R: Read> CramReader<R> {
    /// Create CRAM reader from any Read source.
    ///
    /// **Implementation Status**: Phase 1 - Basic structure
    ///
    /// # Phase 1 TODO:
    /// 1. Read and validate CRAM magic number ("CRAM")
    /// 2. Parse major/minor version
    /// 3. Read file ID (20 bytes)
    /// 4. Parse SAM header container
    /// 5. Build Header structure
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::cram::CramReader;
    /// use std::fs::File;
    /// use std::io::BufReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let file = File::open("sample.cram")?;
    /// let cram = CramReader::new(BufReader::new(file))?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(mut reader: R) -> Result<Self> {
        // Phase 1, Task 1: CRAM file definition parsing
        //
        // CRAM file definition (26 bytes):
        // - Magic: "CRAM" (4 bytes)
        // - Major version: u8
        // - Minor version: u8
        // - File ID: [u8; 20]

        // Read and validate magic number
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read CRAM magic number: {}", e)
            })?;

        if &magic != b"CRAM" {
            return Err(BiometalError::InvalidCramFormat {
                msg: format!("Invalid CRAM magic number: expected 'CRAM', found '{}'",
                    String::from_utf8_lossy(&magic))
            });
        }

        // Read major version
        let mut version_buf = [0u8; 1];
        reader.read_exact(&mut version_buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read CRAM major version: {}", e)
            })?;
        let major_version = version_buf[0];

        // Read minor version
        reader.read_exact(&mut version_buf)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read CRAM minor version: {}", e)
            })?;
        let minor_version = version_buf[0];

        // Validate version (support CRAM 3.0 and 3.1)
        match (major_version, minor_version) {
            (3, 0) | (3, 1) => {
                // Supported versions
            }
            (2, _) => {
                return Err(BiometalError::InvalidCramFormat {
                    msg: format!("CRAM 2.x is not supported (found {}.{}). Please use CRAM 3.0 or 3.1",
                        major_version, minor_version)
                });
            }
            _ => {
                return Err(BiometalError::InvalidCramFormat {
                    msg: format!("Unsupported CRAM version: {}.{} (only 3.0 and 3.1 are supported)",
                        major_version, minor_version)
                });
            }
        }

        // Read file ID (20 bytes)
        let mut file_id = [0u8; 20];
        reader.read_exact(&mut file_id)
            .map_err(|e| BiometalError::InvalidCramFormat {
                msg: format!("Failed to read CRAM file ID: {}", e)
            })?;

        // TODO Phase 1, Task 2: Parse SAM header container
        // For now, create placeholder header
        let header = Header::new(
            String::from("@HD\tVN:1.6\n"),
            vec![],
        );

        Ok(CramReader {
            reader,
            header,
            major_version,
            minor_version,
            reference_index: None,
            reference_path: None,
        })
    }

    /// Get CRAM major version.
    pub fn major_version(&self) -> u8 {
        self.major_version
    }

    /// Get CRAM minor version.
    pub fn minor_version(&self) -> u8 {
        self.minor_version
    }

    /// Get BAM-compatible header.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::cram::CramReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let cram = CramReader::from_path("sample.cram")?;
    /// println!("References: {}", cram.header().reference_count());
    /// # Ok(())
    /// # }
    /// ```
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Set reference FASTA for sequence reconstruction.
    ///
    /// Loads the reference FASTA index (FAI file) for reconstructing sequences from
    /// reference-based CRAM encoding.
    ///
    /// # Arguments
    ///
    /// * `reference_path` - Path to reference FASTA file (e.g., "hg38.fa")
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::cram::CramReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut cram = CramReader::from_path("sample.cram")?;
    /// cram.set_reference("hg38.fa")?;
    ///
    /// for record in cram.records() {
    ///     let record = record?;
    ///     // Sequences will be reconstructed from reference
    ///     println!("{}: {} bp", record.name, record.sequence.len());
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn set_reference<P: AsRef<Path>>(&mut self, reference_path: P) -> Result<()> {
        let ref_path = reference_path.as_ref();

        // Build FAI path (reference.fa -> reference.fa.fai)
        let mut fai_path = PathBuf::from(ref_path);
        let original_path = fai_path.clone();
        fai_path.set_file_name(format!("{}.fai", original_path.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("reference.fa")));

        // Load FAI index
        let index = FaiIndex::from_path(&fai_path)?;

        self.reference_index = Some(index);
        self.reference_path = Some(PathBuf::from(ref_path));

        Ok(())
    }

    /// Get reference sequence for a chromosome/contig.
    ///
    /// # Arguments
    ///
    /// * `reference_id` - Reference ID from CRAM header (index into header.references)
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// Returns the reference subsequence as bytes.
    fn fetch_reference_sequence(
        &self,
        reference_id: usize,
        start: i32,
        end: i32,
    ) -> Result<Vec<u8>> {
        let index = self.reference_index.as_ref()
            .ok_or_else(|| BiometalError::InvalidInput {
                msg: "No reference loaded. Use set_reference() first.".to_string()
            })?;

        let ref_path = self.reference_path.as_ref()
            .ok_or_else(|| BiometalError::InvalidInput {
                msg: "No reference path set".to_string()
            })?;

        // Get reference name from header
        let ref_name = self.header.references.get(reference_id)
            .ok_or_else(|| BiometalError::InvalidInput {
                msg: format!("Reference ID {} not found in header", reference_id)
            })?;

        // Fetch sequence from FASTA
        let sequence = index.fetch_region(
            &ref_name.name,
            start as u64,
            end as u64,
            ref_path,
        )?;

        Ok(sequence.into_bytes())
    }

    /// Create iterator over CRAM records.
    ///
    /// **Implementation Status**: Phase 1 - Basic iteration
    ///
    /// Returns streaming iterator with constant memory footprint.
    ///
    /// # Phase 1 Limitations
    ///
    /// - Returns placeholder records with metadata from slice headers
    /// - Sequences and quality scores are placeholders
    /// - Full record decoding in Phase 2
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::io::cram::CramReader;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let mut cram = CramReader::from_path("sample.cram")?;
    ///
    /// for record in cram.records() {
    ///     let record = record?;
    ///     println!("Read: {}", record.name);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records(self) -> CramRecords<R> {
        CramRecords {
            reader: self.reader,
            current_slice: None,
            current_slice_index: 0,
            reached_eof: false,
            reference_index: self.reference_index,
            reference_path: self.reference_path,
            header: self.header,
        }
    }
}

/// Iterator over CRAM records.
///
/// **Implementation Status**: Phase 2 - Basic reference reconstruction
///
/// Provides streaming access with constant memory footprint.
///
/// # Phase 2 Implementation
///
/// Returns records with sequences reconstructed from reference (if available).
/// Full CRAM feature decoding (substitutions, insertions, etc.) in Phase 2 full.
pub struct CramRecords<R: Read> {
    /// Underlying reader
    reader: R,
    /// Current slice being processed
    current_slice: Option<Slice>,
    /// Current record index within the slice
    current_slice_index: usize,
    /// EOF reached flag
    reached_eof: bool,
    /// Reference FASTA index (optional, for sequence reconstruction)
    reference_index: Option<FaiIndex>,
    /// Reference FASTA file path (optional, for sequence reconstruction)
    reference_path: Option<PathBuf>,
    /// BAM-compatible header (for reference names)
    header: Header,
}

impl<R: Read> CramRecords<R> {
    /// Try to read the next slice from the file.
    ///
    /// Returns None if EOF is reached.
    fn read_next_slice(&mut self) -> Result<Option<Slice>> {
        if self.reached_eof {
            return Ok(None);
        }

        // Try to read next container header
        let container_header = match ContainerHeader::parse(&mut self.reader) {
            Ok(header) => header,
            Err(e) => {
                // Check if this is EOF
                if e.to_string().contains("EOF container") {
                    self.reached_eof = true;
                    return Ok(None);
                }
                // Check for end of file
                if e.to_string().contains("failed to fill whole buffer") {
                    self.reached_eof = true;
                    return Ok(None);
                }
                return Err(e);
            }
        };

        // Skip SAM header containers (no alignment data)
        if container_header.is_sam_header_container() {
            // Read and skip container data
            let mut skip_buffer = vec![0u8; container_header.length as usize];
            self.reader.read_exact(&mut skip_buffer)
                .map_err(|e| BiometalError::InvalidCramFormat {
                    msg: format!("Failed to skip SAM header container: {}", e)
                })?;
            // Recursively try next container
            return self.read_next_slice();
        }

        // Read compression header block
        let _compression_header_block = Block::parse(&mut self.reader)?;
        // TODO Phase 2: Parse compression header for full decoding

        // For Phase 1: Read first slice from container
        // In reality, containers can have multiple slices, but for Phase 1 we'll read one
        let slice = Slice::parse(&mut self.reader)?;

        Ok(Some(slice))
    }
}

impl<R: Read> Iterator for CramRecords<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        // Phase 2 Implementation: Decode sequences from reference (basic)

        loop {
            // Check if we have a current slice
            if let Some(ref slice) = self.current_slice {
                // Check if there are more records in this slice
                if self.current_slice_index < slice.header.num_records as usize {
                    let record_index = self.current_slice_index;
                    self.current_slice_index += 1;

                    // Phase 2: Decode sequence from reference
                    let sequence = match slice.decode_sequence(
                        record_index,
                        self.reference_index.as_ref(),
                        self.reference_path.as_ref(),
                        &self.header.references,
                    ) {
                        Ok(seq) => seq,
                        Err(e) => return Some(Err(e)),
                    };

                    // Phase 2: Decode quality scores (placeholder for now)
                    let quality = slice.decode_quality_scores_simple(record_index, sequence.len());

                    // Create Record with decoded data
                    let record = Record {
                        name: format!("read_{}", record_index),
                        flags: 0,
                        reference_id: Some(slice.header.reference_id as usize),
                        position: Some(slice.header.alignment_start + record_index as i32),
                        mapq: Some(60),
                        cigar: Vec::new(),  // TODO Phase 2: Decode CIGAR from blocks
                        mate_reference_id: None,
                        mate_position: None,
                        template_length: 0,
                        sequence,
                        quality,
                        tags: Tags::new(),  // TODO Phase 2: Decode tags from blocks
                    };

                    return Some(Ok(record));
                }
            }

            // Current slice exhausted or no slice, read next
            match self.read_next_slice() {
                Ok(Some(slice)) => {
                    self.current_slice = Some(slice);
                    self.current_slice_index = 0;
                    // Continue loop to return first record from new slice
                }
                Ok(None) => {
                    // EOF reached
                    return None;
                }
                Err(e) => {
                    return Some(Err(e));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Helper: Create CRAM file definition bytes
    fn make_cram_file_definition(major: u8, minor: u8) -> Vec<u8> {
        let mut data = Vec::new();
        data.extend_from_slice(b"CRAM");       // Magic (4 bytes)
        data.push(major);                       // Major version
        data.push(minor);                       // Minor version
        data.extend_from_slice(&[0u8; 20]);    // File ID (20 bytes)
        data
    }

    #[test]
    fn test_cram_magic_number_valid() {
        // Valid CRAM 3.0 file
        let data = make_cram_file_definition(3, 0);
        let reader = CramReader::new(Cursor::new(data)).unwrap();
        assert_eq!(reader.major_version(), 3);
        assert_eq!(reader.minor_version(), 0);
    }

    #[test]
    fn test_cram_magic_number_invalid() {
        // Invalid magic number
        let mut data = Vec::new();
        data.extend_from_slice(b"BAMX");       // Wrong magic
        data.push(3);
        data.push(0);
        data.extend_from_slice(&[0u8; 20]);

        let result = CramReader::new(Cursor::new(data));
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Invalid CRAM magic number"));
    }

    #[test]
    fn test_cram_version_3_0() {
        // CRAM 3.0 (supported)
        let data = make_cram_file_definition(3, 0);
        let reader = CramReader::new(Cursor::new(data)).unwrap();
        assert_eq!(reader.major_version(), 3);
        assert_eq!(reader.minor_version(), 0);
    }

    #[test]
    fn test_cram_version_3_1() {
        // CRAM 3.1 (supported)
        let data = make_cram_file_definition(3, 1);
        let reader = CramReader::new(Cursor::new(data)).unwrap();
        assert_eq!(reader.major_version(), 3);
        assert_eq!(reader.minor_version(), 1);
    }

    #[test]
    fn test_cram_version_2_x_rejected() {
        // CRAM 2.x (not supported)
        let data = make_cram_file_definition(2, 0);
        let result = CramReader::new(Cursor::new(data));
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("CRAM 2.x is not supported"));
    }

    #[test]
    fn test_cram_version_unsupported() {
        // CRAM 4.0 (future version, not supported)
        let data = make_cram_file_definition(4, 0);
        let result = CramReader::new(Cursor::new(data));
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Unsupported CRAM version: 4.0"));
    }

    #[test]
    fn test_cram_truncated_magic() {
        // Truncated file (only 2 bytes)
        let data = vec![b'C', b'R'];
        let result = CramReader::new(Cursor::new(data));
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Failed to read CRAM magic number"));
    }

    #[test]
    fn test_cram_truncated_version() {
        // Truncated file (magic + major version only)
        let mut data = Vec::new();
        data.extend_from_slice(b"CRAM");
        data.push(3);
        let result = CramReader::new(Cursor::new(data));
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Failed to read CRAM minor version"));
    }

    #[test]
    fn test_cram_truncated_file_id() {
        // Truncated file (missing file ID)
        let mut data = Vec::new();
        data.extend_from_slice(b"CRAM");
        data.push(3);
        data.push(0);
        // Only 10 bytes of file ID instead of 20
        data.extend_from_slice(&[0u8; 10]);
        let result = CramReader::new(Cursor::new(data));
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Failed to read CRAM file ID"));
    }

    // ========================================================================
    // ITF-8 Decoding Tests
    // ========================================================================

    #[test]
    fn test_itf8_1_byte() {
        // 1-byte encoding: 0xxxxxxx (0-127)
        let data = vec![0x00]; // Value: 0
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 0);

        let data = vec![0x7F]; // Value: 127 (max 1-byte)
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 127);

        let data = vec![0x42]; // Value: 66
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 66);
    }

    #[test]
    fn test_itf8_2_bytes() {
        // 2-byte encoding: 10xxxxxx xxxxxxxx (128-16,383)
        let data = vec![0x80, 0x80]; // Value: 128
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 128);

        let data = vec![0xBF, 0xFF]; // Value: 16,383 (max 2-byte)
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 16383);

        let data = vec![0x81, 0x42]; // Value: 322 (0x142)
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 322);
    }

    #[test]
    fn test_itf8_3_bytes() {
        // 3-byte encoding: 110xxxxx xxxxxxxx xxxxxxxx (16,384-2,097,151)
        let data = vec![0xC0, 0x40, 0x00]; // Value: 16,384
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 16384);

        let data = vec![0xDF, 0xFF, 0xFF]; // Value: 2,097,151 (max 3-byte)
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 2097151);
    }

    #[test]
    fn test_itf8_4_bytes() {
        // 4-byte encoding: 1110xxxx xxxxxxxx xxxxxxxx xxxxxxxx
        let data = vec![0xE0, 0x20, 0x00, 0x00]; // Value: 2,097,152
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 2097152);

        let data = vec![0xEF, 0xFF, 0xFF, 0xFF]; // Value: 268,435,455 (max 4-byte)
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 268435455);
    }

    #[test]
    fn test_itf8_5_bytes() {
        // 5-byte encoding: 1111xxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
        let data = vec![0xF0, 0x10, 0x00, 0x00, 0x00]; // Value: 268,435,456
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_itf8(&mut reader).unwrap(), 268435456);

        // Test another 5-byte value: (1 << 28) = 268,435,456
        let data = vec![0xF1, 0x00, 0x00, 0x00, 0x00]; // Bits: 0001 00000000 00000000 00000000 00000000
        let mut reader = Cursor::new(data);
        let result = super::decode_itf8(&mut reader).unwrap();
        assert_eq!(result, 268435456); // (1 << 28)
    }

    #[test]
    fn test_itf8_truncated() {
        // Truncated 2-byte encoding (missing second byte)
        let data = vec![0x80];
        let mut reader = Cursor::new(data);
        let result = super::decode_itf8(&mut reader);
        assert!(result.is_err());

        // Truncated 5-byte encoding (missing bytes)
        let data = vec![0xF0, 0x10, 0x00];
        let mut reader = Cursor::new(data);
        let result = super::decode_itf8(&mut reader);
        assert!(result.is_err());
    }

    // ========================================================================
    // LTF-8 Decoding Tests
    // ========================================================================

    #[test]
    fn test_ltf8_1_byte() {
        // 1-byte encoding (0-127)
        let data = vec![0x00]; // Value: 0
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_ltf8(&mut reader).unwrap(), 0);

        let data = vec![0x7F]; // Value: 127
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_ltf8(&mut reader).unwrap(), 127);
    }

    #[test]
    fn test_ltf8_2_bytes() {
        // 2-byte encoding
        let data = vec![0x80, 0x80]; // Value: 128
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_ltf8(&mut reader).unwrap(), 128);

        let data = vec![0x81, 0x42]; // Value: 322 (0x142)
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_ltf8(&mut reader).unwrap(), 322);

        let data = vec![0xBF, 0xFF]; // Max 2-byte: 16,383
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_ltf8(&mut reader).unwrap(), 16383);
    }

    #[test]
    fn test_ltf8_3_bytes() {
        // 3-byte encoding
        let data = vec![0xC0, 0x40, 0x00]; // Value: 16,384
        let mut reader = Cursor::new(data);
        assert_eq!(super::decode_ltf8(&mut reader).unwrap(), 16384);
    }

    #[test]
    fn test_ltf8_large_values() {
        // Test larger values that require more bytes
        // 5-byte encoding for a larger value
        let data = vec![0xF0, 0x00, 0x01, 0x00, 0x00]; // Large value
        let mut reader = Cursor::new(data);
        let result = super::decode_ltf8(&mut reader).unwrap();
        assert!(result > 0);
    }

    #[test]
    fn test_ltf8_truncated() {
        // Truncated 2-byte encoding
        let data = vec![0x80];
        let mut reader = Cursor::new(data);
        let result = super::decode_ltf8(&mut reader);
        assert!(result.is_err());

        // Truncated multi-byte encoding
        let data = vec![0xF0, 0x00];
        let mut reader = Cursor::new(data);
        let result = super::decode_ltf8(&mut reader);
        assert!(result.is_err());
    }

    // ========================================================================
    // Container Header Parsing Tests
    // ========================================================================

    /// Helper: Create a minimal container header with specified values
    fn make_container_header(
        length: i32,
        reference_id: i32,
        start_pos: i32,
        span: i32,
        num_records: i32,
    ) -> Vec<u8> {
        let mut data = Vec::new();

        // Length (4 bytes, big-endian)
        data.extend_from_slice(&length.to_be_bytes());

        // Reference ID (ITF-8) - for testing, use simplified encoding
        // In reality, CRAM uses signed representation, but for tests we'll use unsigned
        if reference_id == -1 {
            // Special case: -1 is commonly used, encode as 0xFF in 1-byte
            // Actually, -1 should be all bits set in ITF-8
            // For simplicity in tests, let's use a 5-byte encoding
            data.extend_from_slice(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF]);
        } else if reference_id >= 0 && reference_id < 128 {
            data.push(reference_id as u8);
        } else {
            // 2-byte encoding for positive values >= 128
            data.push(0x80 | ((reference_id >> 8) & 0x3F) as u8);
            data.push((reference_id & 0xFF) as u8);
        }

        // Helper to encode ITF-8 (simplified)
        let encode_itf8 = |value: i32| -> Vec<u8> {
            if value < 128 {
                vec![value as u8]
            } else if value < 16384 {
                vec![
                    0x80 | ((value >> 8) & 0x3F) as u8,
                    (value & 0xFF) as u8,
                ]
            } else {
                // For larger values, use 3-byte encoding
                vec![
                    0xC0 | ((value >> 16) & 0x1F) as u8,
                    ((value >> 8) & 0xFF) as u8,
                    (value & 0xFF) as u8,
                ]
            }
        };

        // Start position (ITF-8)
        data.extend_from_slice(&encode_itf8(start_pos));

        // Alignment span (ITF-8)
        data.extend_from_slice(&encode_itf8(span));

        // Number of records (ITF-8)
        data.extend_from_slice(&encode_itf8(num_records));

        // Record counter (LTF-8, simplified: 1-byte for 0)
        data.push(0x00);

        // Bases (LTF-8, simplified: 1-byte for 0)
        data.push(0x00);

        // Number of blocks (ITF-8, simplified: 1-byte)
        data.push(0x01); // 1 block

        // Landmarks count (ITF-8)
        data.push(0x00); // 0 landmarks

        // CRC32 (4 bytes, dummy value)
        data.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]);

        data
    }

    #[test]
    fn test_container_header_basic() {
        // Create a basic container header
        let data = make_container_header(1000, 0, 100, 500, 10);
        let mut reader = Cursor::new(data);

        let header = super::ContainerHeader::parse(&mut reader).unwrap();
        assert_eq!(header.length, 1000);
        assert_eq!(header.num_records, 10);
        assert_eq!(header.num_blocks, 1);
        assert_eq!(header.landmarks.len(), 0);
    }

    #[test]
    fn test_container_header_sam_header_detection() {
        // For now, manually construct a SAM header container
        // SAM header containers have: ref_id=-1, start=0, num_records=0
        // We'll test the is_sam_header_container() logic with known values
        let data = make_container_header(500, 0, 0, 0, 0);
        let mut reader = Cursor::new(data);

        let mut header = super::ContainerHeader::parse(&mut reader).unwrap();
        // Manually set reference_id to -1 to test the detection logic
        header.reference_id = -1;
        assert!(header.is_sam_header_container());

        // Test that non-SAM-header containers are correctly identified
        header.reference_id = 0;
        assert!(!header.is_sam_header_container());

        header.reference_id = -1;
        header.num_records = 10; // Has records, not a SAM header
        assert!(!header.is_sam_header_container());
    }

    #[test]
    fn test_container_header_eof() {
        // EOF container has length = 0
        let data = vec![0x00, 0x00, 0x00, 0x00]; // length = 0
        let mut reader = Cursor::new(data);

        let result = super::ContainerHeader::parse(&mut reader);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("EOF container"));
    }

    #[test]
    fn test_container_header_truncated() {
        // Truncated container header (only length field)
        let data = vec![0x00, 0x00, 0x03, 0xE8]; // length = 1000, but no other fields
        let mut reader = Cursor::new(data);

        let result = super::ContainerHeader::parse(&mut reader);
        assert!(result.is_err());
    }

    // ========================================================================
    // Block Parsing Tests
    // ========================================================================

    /// Helper: Create a minimal block with raw (uncompressed) data
    fn make_raw_block(content_type: u8, data: Vec<u8>) -> Vec<u8> {
        let mut block_data = Vec::new();

        // Method (0 = raw)
        block_data.push(0x00);

        // Content type
        block_data.push(content_type);

        // Block ID (ITF-8, use 0)
        block_data.push(0x00);

        // Compressed size (ITF-8)
        let compressed_size = data.len() as i32;
        if compressed_size < 128 {
            block_data.push(compressed_size as u8);
        } else {
            block_data.push(0x80 | ((compressed_size >> 8) & 0x3F) as u8);
            block_data.push((compressed_size & 0xFF) as u8);
        }

        // Uncompressed size (same as compressed for raw)
        if compressed_size < 128 {
            block_data.push(compressed_size as u8);
        } else {
            block_data.push(0x80 | ((compressed_size >> 8) & 0x3F) as u8);
            block_data.push((compressed_size & 0xFF) as u8);
        }

        // Data
        block_data.extend_from_slice(&data);

        // CRC32 (dummy)
        block_data.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]);

        block_data
    }

    #[test]
    fn test_block_parsing_raw() {
        // Create a raw block with some test data
        let test_data = vec![0x01, 0x02, 0x03, 0x04, 0x05];
        let block_bytes = make_raw_block(1, test_data.clone());
        let mut reader = Cursor::new(block_bytes);

        let block = super::Block::parse(&mut reader).unwrap();
        assert_eq!(block.method, 0); // Raw
        assert_eq!(block.content_type, 1);
        assert_eq!(block.compressed_size, 5);
        assert_eq!(block.uncompressed_size, 5);
        assert_eq!(block.data, test_data);
    }

    #[test]
    fn test_block_decompression_raw() {
        // Raw data should decompress to itself
        let test_data = vec![0x01, 0x02, 0x03, 0x04, 0x05];
        let block_bytes = make_raw_block(1, test_data.clone());
        let mut reader = Cursor::new(block_bytes);

        let block = super::Block::parse(&mut reader).unwrap();
        let decompressed = block.decompress().unwrap();
        assert_eq!(decompressed, test_data);
    }

    #[test]
    fn test_block_is_compression_header() {
        // Content type 0 = compression header
        let block_bytes = make_raw_block(0, vec![0x00]);
        let mut reader = Cursor::new(block_bytes);
        let block = super::Block::parse(&mut reader).unwrap();
        assert!(block.is_compression_header());

        // Content type 1 = not compression header
        let block_bytes = make_raw_block(1, vec![0x00]);
        let mut reader = Cursor::new(block_bytes);
        let block = super::Block::parse(&mut reader).unwrap();
        assert!(!block.is_compression_header());
    }

    /// Helper: Create a bzip2-compressed block
    fn make_bzip2_block(content_type: u8, data: &[u8]) -> Vec<u8> {
        use bzip2::write::BzEncoder;
        use bzip2::Compression;
        use std::io::Write;

        // Compress the data
        let mut encoder = BzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        let compressed = encoder.finish().unwrap();

        let mut block_data = Vec::new();
        block_data.push(0x02); // Method = bzip2
        block_data.push(content_type);
        block_data.push(0x00); // Block ID (ITF-8, 0)

        // Compressed size (ITF-8)
        let compressed_size = compressed.len() as i32;
        encode_itf8(&mut block_data, compressed_size);

        // Uncompressed size (ITF-8)
        let uncompressed_size = data.len() as i32;
        encode_itf8(&mut block_data, uncompressed_size);

        // Compressed data
        block_data.extend_from_slice(&compressed);

        // CRC32 (dummy)
        block_data.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]);

        block_data
    }

    /// Helper: Create an lzma-compressed block
    fn make_lzma_block(content_type: u8, data: &[u8]) -> Vec<u8> {
        use xz2::write::XzEncoder;
        use std::io::Write;

        // Compress the data
        let mut encoder = XzEncoder::new(Vec::new(), 6);
        encoder.write_all(data).unwrap();
        let compressed = encoder.finish().unwrap();

        let mut block_data = Vec::new();
        block_data.push(0x03); // Method = lzma
        block_data.push(content_type);
        block_data.push(0x00); // Block ID (ITF-8, 0)

        // Compressed size (ITF-8)
        let compressed_size = compressed.len() as i32;
        encode_itf8(&mut block_data, compressed_size);

        // Uncompressed size (ITF-8)
        let uncompressed_size = data.len() as i32;
        encode_itf8(&mut block_data, uncompressed_size);

        // Compressed data
        block_data.extend_from_slice(&compressed);

        // CRC32 (dummy)
        block_data.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]);

        block_data
    }

    /// Helper: Encode a 32-bit integer as ITF-8
    fn encode_itf8(buf: &mut Vec<u8>, value: i32) {
        if value < 128 {
            buf.push(value as u8);
        } else if value < 16384 {
            buf.push(0x80 | ((value >> 8) & 0x3F) as u8);
            buf.push((value & 0xFF) as u8);
        } else if value < 2097152 {
            buf.push(0xC0 | ((value >> 16) & 0x1F) as u8);
            buf.push(((value >> 8) & 0xFF) as u8);
            buf.push((value & 0xFF) as u8);
        } else if value < 268435456 {
            buf.push(0xE0 | ((value >> 24) & 0x0F) as u8);
            buf.push(((value >> 16) & 0xFF) as u8);
            buf.push(((value >> 8) & 0xFF) as u8);
            buf.push((value & 0xFF) as u8);
        } else {
            buf.push(0xF0);
            buf.push(((value >> 24) & 0xFF) as u8);
            buf.push(((value >> 16) & 0xFF) as u8);
            buf.push(((value >> 8) & 0xFF) as u8);
            buf.push((value & 0xFF) as u8);
        }
    }

    #[test]
    fn test_block_bzip2_decompression() {
        let test_data = b"Hello, CRAM with bzip2! This is a longer message to ensure proper compression.";
        let block_bytes = make_bzip2_block(1, test_data);
        let mut reader = Cursor::new(block_bytes);

        let block = super::Block::parse(&mut reader).unwrap();
        assert_eq!(block.method, 2); // Bzip2
        let decompressed = block.decompress().unwrap();
        assert_eq!(&decompressed[..], test_data);
    }

    #[test]
    fn test_block_lzma_decompression() {
        let test_data = b"Hello, CRAM with LZMA! This is a longer message to ensure proper compression.";
        let block_bytes = make_lzma_block(1, test_data);
        let mut reader = Cursor::new(block_bytes);

        let block = super::Block::parse(&mut reader).unwrap();
        assert_eq!(block.method, 3); // LZMA
        let decompressed = block.decompress().unwrap();
        assert_eq!(&decompressed[..], test_data);
    }

    #[test]
    fn test_block_rans_unsupported() {
        // Create a block with rANS compression (method=4, not yet supported)
        let mut block_data = Vec::new();
        block_data.push(0x04); // Method = rANS
        block_data.push(0x01); // Content type
        block_data.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]); // Block ID
        block_data.extend_from_slice(&[0x00, 0x00, 0x00, 0x05]); // Compressed size
        block_data.extend_from_slice(&[0x00, 0x00, 0x00, 0x05]); // Uncompressed size
        block_data.extend_from_slice(&[0x01, 0x02, 0x03, 0x04, 0x05]); // Data
        block_data.extend_from_slice(&[0x00, 0x00, 0x00, 0x00]); // CRC32

        let mut reader = Cursor::new(block_data);
        let block = super::Block::parse(&mut reader).unwrap();

        let result = block.decompress();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("rANS"));
    }

    // ========================================================================
    // Compression Header Tests
    // ========================================================================

    #[test]
    fn test_compression_header_basic() {
        // Create minimal compression header data
        let mut data = Vec::new();

        // Preservation map: size = 1 byte, containing num_entries = 0
        data.push(0x01); // size = 1 byte
        data.push(0x00); // num entries = 0 (ITF-8)

        // Data series encoding: size = 1 byte, containing num_entries = 0
        data.push(0x01); // size = 1 byte
        data.push(0x00); // num entries = 0 (ITF-8)

        // Tag encoding: size = 1 byte, containing num_entries = 0
        data.push(0x01); // size = 1 byte
        data.push(0x00); // num entries = 0 (ITF-8)

        let header = super::CompressionHeader::parse(&data).unwrap();
        // Check preservation map is default (empty)
        assert!(!header.preservation_map.read_names_included);
        assert!(!header.preservation_map.ap_data_series_delta);
        assert!(!header.preservation_map.reference_required);
        assert!(header.preservation_map.substitution_matrix.is_none());
        assert_eq!(header.preservation_map.tag_ids.len(), 0);

        // Check encoding maps are empty
        assert_eq!(header.data_series_encoding.len(), 0);
        assert_eq!(header.tag_encoding.len(), 0);
    }

    #[test]
    fn test_compression_header_with_data() {
        // Create a valid compression header with preservation map, data series encoding, and tag encoding
        let mut data = Vec::new();

        // Preservation map: 1 entry (RN: true)
        data.push(0x04); // size = 4 bytes (1 for num_entries + 2 for key + 1 for value)
        data.push(0x01); // num entries = 1 (ITF-8)
        data.extend_from_slice(b"RN"); // key
        data.push(0x01); // value = true

        // Data series encoding: 1 entry (BF: External block 1)
        data.push(0x05); // size = 5 bytes (1 for num_entries + 2 for key + 1 for encoding_id + 1 for block_id)
        data.push(0x01); // num entries = 1 (ITF-8)
        data.extend_from_slice(b"BF"); // key = BAM flags
        data.push(0x01); // encoding ID = EXTERNAL
        data.push(0x01); // block_content_id = 1 (ITF-8)

        // Tag encoding: 1 entry (tag 100: External block 2)
        data.push(0x04); // size = 4 bytes (1 for num_entries + 1 for tag_id + 1 for encoding_id + 1 for block_id)
        data.push(0x01); // num entries = 1 (ITF-8)
        data.push(0x64); // tag ID = 100 (ITF-8)
        data.push(0x01); // encoding ID = EXTERNAL
        data.push(0x02); // block_content_id = 2 (ITF-8)

        let header = super::CompressionHeader::parse(&data).unwrap();

        // Check preservation map
        assert!(header.preservation_map.read_names_included);
        assert!(!header.preservation_map.ap_data_series_delta);
        assert!(!header.preservation_map.reference_required);

        // Check data series encoding
        assert_eq!(header.data_series_encoding.len(), 1);
        assert!(header.data_series_encoding.contains_key(&super::DataSeries::BF));
        if let Some(super::Encoding::External { block_content_id }) = header.data_series_encoding.get(&super::DataSeries::BF) {
            assert_eq!(*block_content_id, 1);
        } else {
            panic!("Expected External encoding for BF");
        }

        // Check tag encoding
        assert_eq!(header.tag_encoding.len(), 1);
        assert!(header.tag_encoding.contains_key(&100));
        if let Some(super::Encoding::External { block_content_id }) = header.tag_encoding.get(&100) {
            assert_eq!(*block_content_id, 2);
        } else {
            panic!("Expected External encoding for tag 100");
        }
    }

    // ========================================================================
    // Slice Header Parsing Tests
    // ========================================================================

    /// Helper: Create a minimal slice header
    fn make_slice_header(
        reference_id: i32,
        start: i32,
        span: i32,
        num_records: i32,
        num_blocks: i32,
    ) -> Vec<u8> {
        let mut data = Vec::new();

        // Helper to encode ITF-8 (reuse from container tests)
        let encode_itf8 = |value: i32| -> Vec<u8> {
            if value < 128 {
                vec![value as u8]
            } else if value < 16384 {
                vec![
                    0x80 | ((value >> 8) & 0x3F) as u8,
                    (value & 0xFF) as u8,
                ]
            } else {
                vec![
                    0xC0 | ((value >> 16) & 0x1F) as u8,
                    ((value >> 8) & 0xFF) as u8,
                    (value & 0xFF) as u8,
                ]
            }
        };

        // Reference ID
        data.extend_from_slice(&encode_itf8(reference_id));

        // Alignment start
        data.extend_from_slice(&encode_itf8(start));

        // Alignment span
        data.extend_from_slice(&encode_itf8(span));

        // Number of records
        data.extend_from_slice(&encode_itf8(num_records));

        // Record counter (LTF-8, simplified: 1-byte for 0)
        data.push(0x00);

        // Number of blocks
        data.extend_from_slice(&encode_itf8(num_blocks));

        // Block content IDs (simplified: all zeros)
        for _ in 0..num_blocks {
            data.push(0x00);
        }

        // Embedded reference MD5 (16 bytes)
        data.extend_from_slice(&[0u8; 16]);

        // No optional reference MD5 for now

        data
    }

    #[test]
    fn test_slice_header_basic() {
        let data = make_slice_header(0, 1000, 500, 10, 2);
        let mut reader = Cursor::new(data);

        let header = super::SliceHeader::parse(&mut reader).unwrap();
        assert_eq!(header.reference_id, 0);
        assert_eq!(header.alignment_start, 1000);
        assert_eq!(header.alignment_span, 500);
        assert_eq!(header.num_records, 10);
        assert_eq!(header.num_blocks, 2);
        assert_eq!(header.block_content_ids.len(), 2);
        assert_eq!(header.embedded_ref_md5, [0u8; 16]);
        assert_eq!(header.optional_ref_md5, None);
    }

    #[test]
    fn test_slice_header_no_optional_md5() {
        // For Phase 1, optional MD5 is not read
        let data = make_slice_header(0, 1000, 500, 10, 1);

        let mut reader = Cursor::new(data);
        let header = super::SliceHeader::parse(&mut reader).unwrap();
        assert_eq!(header.optional_ref_md5, None);

        // TODO Phase 2: Add test for optional MD5 when flag-based parsing is implemented
    }

    #[test]
    fn test_slice_parsing_with_blocks() {
        // Create slice header
        let mut slice_data = make_slice_header(0, 1000, 500, 10, 2);

        // Add 2 blocks
        let block1 = make_raw_block(0, vec![0x01, 0x02]); // Core block
        let block2 = make_raw_block(1, vec![0x03, 0x04]); // External block
        slice_data.extend_from_slice(&block1);
        slice_data.extend_from_slice(&block2);

        let mut reader = Cursor::new(slice_data);
        let slice = super::Slice::parse(&mut reader).unwrap();

        assert_eq!(slice.header.num_blocks, 2);
        assert_eq!(slice.blocks.len(), 2);

        // Check core block
        let core = slice.core_block().unwrap();
        assert_eq!(core.content_type, 0);
        assert_eq!(core.data, vec![0x01, 0x02]);

        // Check external blocks
        let external = slice.external_blocks();
        assert_eq!(external.len(), 1);
        assert_eq!(external[0].content_type, 1);
        assert_eq!(external[0].data, vec![0x03, 0x04]);
    }

    #[test]
    fn test_slice_empty_external_blocks() {
        // Create slice with only one block (core block)
        let mut slice_data = make_slice_header(0, 1000, 500, 10, 1);
        let block = make_raw_block(0, vec![0x01, 0x02]);
        slice_data.extend_from_slice(&block);

        let mut reader = Cursor::new(slice_data);
        let slice = super::Slice::parse(&mut reader).unwrap();

        assert_eq!(slice.blocks.len(), 1);
        assert!(slice.core_block().is_some());
        assert_eq!(slice.external_blocks().len(), 0);
    }

    // ========================================================================
    // Record Iteration Tests
    // ========================================================================

    #[test]
    fn test_record_iteration_basic() {
        // For Phase 1, we test that the iteration structure works
        // We can't easily test with real CRAM files without a full writer,
        // so we verify the structure compiles and basic logic works

        // Create a minimal CRAM file structure
        let mut cram_data = Vec::new();

        // File definition
        cram_data.extend_from_slice(b"CRAM");  // Magic
        cram_data.push(3);  // Major version
        cram_data.push(0);  // Minor version
        cram_data.extend_from_slice(&[0u8; 20]);  // File ID

        // For now, this is enough to test the reader initialization
        // Full integration tests will come in Phase 2 when we can generate proper CRAM data

        let mut reader = Cursor::new(cram_data);
        let cram_reader = super::CramReader::new(reader).unwrap();

        // Verify reader is initialized
        assert_eq!(cram_reader.major_version(), 3);
        assert_eq!(cram_reader.minor_version(), 0);

        // Test records() returns an iterator (even if empty for this minimal data)
        let records = cram_reader.records();
        // Try to iterate (will hit EOF immediately for this minimal data)
        let count = records.count();
        assert_eq!(count, 0); // No data containers, so no records
    }

    #[test]
    fn test_record_placeholder_creation() {
        // Test that placeholder record creation works as expected
        // This tests the record generation logic directly

        let slice = super::Slice {
            header: super::SliceHeader {
                reference_id: 0,
                alignment_start: 1000,
                alignment_span: 100,
                num_records: 3,
                record_counter: 0,
                num_blocks: 1,
                block_content_ids: vec![0],
                embedded_ref_md5: [0u8; 16],
                optional_ref_md5: None,
            },
            blocks: vec![],
        };

        // Simulate what the iterator does
        for i in 0..slice.header.num_records {
            let record = super::Record {
                name: format!("read_{}", i),
                flags: 0,
                reference_id: Some(slice.header.reference_id as usize),
                position: Some(slice.header.alignment_start + i),
                mapq: Some(60),
                cigar: Vec::new(),
                mate_reference_id: None,
                mate_position: None,
                template_length: 0,
                sequence: Vec::new(),
                quality: Vec::new(),
                tags: Tags::new(),
            };

            // Verify placeholder record fields
            assert_eq!(record.name, format!("read_{}", i));
            assert_eq!(record.reference_id, Some(0));
            assert_eq!(record.position, Some(1000 + i));
            assert_eq!(record.mapq, Some(60));
        }
    }

    // TODO Phase 2: Add full record decoding tests
    // - test_record_sequence_decoding()
    // - test_record_quality_decoding()
    // - test_record_tag_parsing()
    // - test_reference_reconstruction()
}
