//! BAM record structure and parsing.
//!
//! A BAM record represents a single alignment (read mapped to reference).
//! Each record contains alignment information, sequence data, quality scores,
//! and optional tags.
//!
//! # Binary Format
//!
//! ```text
//! BAM Record (binary, little-endian):
//! - block_size (int32): Total record size in bytes (excluding this field)
//! - refID (int32): Reference sequence ID (-1 for unmapped)
//! - pos (int32): 0-based leftmost position (-1 for unmapped)
//! - l_read_name (uint8): Length of read name (includes null terminator)
//! - mapq (uint8): Mapping quality (255 = unavailable)
//! - bin (uint16): BAI index bin (for indexing)
//! - n_cigar_op (uint16): Number of CIGAR operations
//! - flag (uint16): Bitwise FLAGS
//! - l_seq (int32): Sequence length
//! - next_refID (int32): Reference ID of mate/next read
//! - next_pos (int32): Position of mate/next read
//! - tlen (int32): Template length
//! - read_name (char[l_read_name]): Null-terminated read name
//! - cigar (uint32[n_cigar_op]): CIGAR operations
//! - seq (uint8[(l_seq+1)/2]): 4-bit encoded sequence
//! - qual (char[l_seq]): Phred quality scores
//! - tags: Optional tags (variable length)
//! ```

use super::cigar::{parse_cigar, CigarOp};
use super::error::BamDecodeError;
use super::sequence::decode_sequence;
use super::tags::{parse_tags, Tags, TagValue, ArrayValue};
use std::io::{self, Read};

/// Safe helper to read i32 from little-endian bytes with bounds checking.
fn read_i32_le(data: &[u8], cursor: &mut usize) -> io::Result<i32> {
    if *cursor + 4 > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("Insufficient data at offset {}: need 4 bytes for i32, got {}", *cursor, data.len() - *cursor),
        ));
    }
    let value = i32::from_le_bytes([
        data[*cursor],
        data[*cursor + 1],
        data[*cursor + 2],
        data[*cursor + 3],
    ]);
    *cursor += 4;
    Ok(value)
}

/// Safe helper to read u16 from little-endian bytes with bounds checking.
fn read_u16_le(data: &[u8], cursor: &mut usize) -> io::Result<u16> {
    if *cursor + 2 > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("Insufficient data at offset {}: need 2 bytes for u16, got {}", *cursor, data.len() - *cursor),
        ));
    }
    let value = u16::from_le_bytes([data[*cursor], data[*cursor + 1]]);
    *cursor += 2;
    Ok(value)
}

/// Safe helper to read u8 with bounds checking.
fn read_u8(data: &[u8], cursor: &mut usize) -> io::Result<u8> {
    if *cursor >= data.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("Insufficient data at offset {}: need 1 byte", *cursor),
        ));
    }
    let value = data[*cursor];
    *cursor += 1;
    Ok(value)
}

/// Validate reference ID per BAM spec.
///
/// BAM spec allows only:
/// - -1 (unmapped)
/// - >= 0 (valid reference)
///
/// Values like -2, -3, etc. are invalid and must be rejected.
fn parse_reference_id(ref_id: i32, field_name: &str) -> io::Result<Option<usize>> {
    const UNMAPPED: i32 = -1;

    match ref_id {
        UNMAPPED => Ok(None),
        n if n >= 0 => Ok(Some(n as usize)),
        invalid => Err(BamDecodeError::InvalidReferenceId {
            value: invalid,
            field: field_name.to_string(),
        }
        .into()),
    }
}

/// BAM alignment record.
///
/// Represents a single read alignment with all associated information.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Read name/query name
    pub name: String,

    /// Reference sequence ID (index into header.references)
    /// None if unmapped (refID = -1)
    pub reference_id: Option<usize>,

    /// 0-based leftmost mapping position
    /// None if unmapped (pos = -1)
    pub position: Option<i32>,

    /// Mapping quality (0-254, 255 = unavailable)
    /// None if unavailable (mapq = 255)
    pub mapq: Option<u8>,

    /// Bitwise FLAGS (see SAM spec for flag meanings)
    pub flags: u16,

    /// Mate/next segment reference ID
    /// None if unavailable (next_refID = -1)
    pub mate_reference_id: Option<usize>,

    /// Mate/next segment position
    /// None if unavailable (next_pos = -1)
    pub mate_position: Option<i32>,

    /// Template length (TLEN)
    pub template_length: i32,

    /// Sequence bases (ASCII: A, C, G, T, N, etc.)
    pub sequence: Vec<u8>,

    /// Quality scores (Phred+33 ASCII)
    pub quality: Vec<u8>,

    /// CIGAR operations
    pub cigar: Vec<CigarOp>,

    /// Optional tags
    pub tags: Tags,
}

impl Record {
    /// Create a new empty record.
    pub fn new() -> Self {
        Self {
            name: String::new(),
            reference_id: None,
            position: None,
            mapq: None,
            flags: 0,
            mate_reference_id: None,
            mate_position: None,
            template_length: 0,
            sequence: Vec::new(),
            quality: Vec::new(),
            cigar: Vec::new(),
            tags: Tags::new(),
        }
    }

    /// Check if the read is unmapped.
    pub fn is_unmapped(&self) -> bool {
        self.flags & 0x4 != 0
    }

    /// Check if the read is paired.
    pub fn is_paired(&self) -> bool {
        self.flags & 0x1 != 0
    }

    /// Check if the read is a reverse complement.
    pub fn is_reverse_complement(&self) -> bool {
        self.flags & 0x10 != 0
    }

    /// Get sequence length.
    pub fn sequence_length(&self) -> usize {
        self.sequence.len()
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::new()
    }
}

/// Check for oversized CIGAR and extract from CG:B,I tag if present.
///
/// BAM format has a 16-bit field for n_cigar_op, limiting CIGAR to 65,535 operations.
/// Long reads (nanopore, PacBio) can exceed this limit. The workaround:
/// - Store placeholder CIGAR with 2 operations: kS*N (k = sequence length, S = soft clip, N = ref skip)
/// - Store actual CIGAR in CG:B,I tag (array of Int32)
///
/// This function detects the pattern and extracts the real CIGAR from the tag.
///
/// # Arguments
///
/// * `cigar` - Parsed CIGAR operations from n_cigar_op field
/// * `sequence_length` - Length of the sequence
/// * `tags` - Parsed auxiliary tags
///
/// # Returns
///
/// Returns the real CIGAR if oversized pattern detected and CG tag found, otherwise returns original CIGAR.
fn check_oversized_cigar(
    cigar: Vec<CigarOp>,
    sequence_length: usize,
    tags: &Tags,
) -> io::Result<Vec<CigarOp>> {
    // Check for oversized CIGAR pattern: 2 operations where first is kS (k = seq length)
    if cigar.len() != 2 {
        return Ok(cigar);
    }

    // First operation should be SoftClip with length = sequence_length
    if let CigarOp::SoftClip(len) = cigar[0] {
        if len as usize != sequence_length {
            return Ok(cigar);
        }
    } else {
        return Ok(cigar);
    }

    // Second operation should be RefSkip (any length)
    if !matches!(cigar[1], CigarOp::RefSkip(_)) {
        return Ok(cigar);
    }

    // Pattern matched! Look for CG:B,I tag
    match tags.get(b"CG")? {
        Some(tag) => match tag.value {
            TagValue::Array(ArrayValue::Int32(raw_cigar)) => {
                // Convert CG tag (Int32 array) to CigarOp vector
                // Each element is a packed 32-bit value: low 4 bits = op, high 28 bits = length
                let mut real_cigar = Vec::with_capacity(raw_cigar.len());

                for &packed in &raw_cigar {
                    let op_code = (packed & 0xF) as u8;
                    let length = (packed >> 4) as u32;

                    let cigar_op = match op_code {
                        0 => CigarOp::Match(length),
                        1 => CigarOp::Insertion(length),
                        2 => CigarOp::Deletion(length),
                        3 => CigarOp::RefSkip(length),
                        4 => CigarOp::SoftClip(length),
                        5 => CigarOp::HardClip(length),
                        6 => CigarOp::Padding(length),
                        7 => CigarOp::SeqMatch(length),
                        8 => CigarOp::SeqMismatch(length),
                        _ => {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Invalid CIGAR operation in CG tag: {}", op_code),
                            ));
                        }
                    };

                    real_cigar.push(cigar_op);
                }

                Ok(real_cigar)
            }
            _ => {
                // CG tag exists but wrong type - ignore and use placeholder CIGAR
                Ok(cigar)
            }
        },
        None => {
            // No CG tag found - this shouldn't happen if pattern matched, but use placeholder
            Ok(cigar)
        }
    }
}

/// Parse a BAM record from binary data.
///
/// # Arguments
///
/// * `data` - Complete record data (including block_size field)
///
/// # Returns
///
/// Parsed BAM record.
///
/// # Errors
///
/// Returns error if:
/// - Data is too short for fixed fields
/// - Sequence decoding fails
/// - CIGAR parsing fails
/// - Tag parsing fails
/// - Name is not valid UTF-8
///
/// # Format Details
///
/// See module documentation for complete binary format specification.
pub fn parse_record(data: &[u8]) -> io::Result<Record> {
    if data.len() < 32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Record too short: {} bytes (minimum 32)", data.len()),
        ));
    }

    let mut cursor = 0;

    // Block size (4 bytes) - we already know this from reading
    let _block_size = read_i32_le(data, &mut cursor)?;

    // Reference ID (4 bytes)
    let ref_id = read_i32_le(data, &mut cursor)?;

    // Position (4 bytes)
    let pos = read_i32_le(data, &mut cursor)?;

    // Read name length (1 byte)
    let l_read_name_u8 = read_u8(data, &mut cursor)?;
    let l_read_name = l_read_name_u8 as usize;

    // BAM spec requires l_read_name >= 1 (minimum "*\0" for unmapped reads)
    if l_read_name == 0 {
        return Err(BamDecodeError::InvalidReadNameLength {
            length: l_read_name_u8,
            offset: cursor - 1,
        }
        .into());
    }

    // MAPQ (1 byte)
    let mapq = read_u8(data, &mut cursor)?;

    // BAI bin (2 bytes) - skip for now
    let _bin = read_u16_le(data, &mut cursor)?;

    // Number of CIGAR operations (2 bytes)
    let n_cigar_op = read_u16_le(data, &mut cursor)? as usize;

    // FLAGS (2 bytes)
    let flags = read_u16_le(data, &mut cursor)?;

    // Sequence length (4 bytes)
    let l_seq = read_i32_le(data, &mut cursor)?;

    if l_seq < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid sequence length at offset {}: {}", cursor - 4, l_seq),
        ));
    }
    let l_seq = l_seq as usize;

    // Next reference ID (4 bytes)
    let next_ref_id = read_i32_le(data, &mut cursor)?;

    // Next position (4 bytes)
    let next_pos = read_i32_le(data, &mut cursor)?;

    // Template length (4 bytes)
    let tlen = read_i32_le(data, &mut cursor)?;

    // Read name (null-terminated)
    if cursor + l_read_name > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Insufficient data for read name at offset {}: need {} bytes, got {}",
                cursor,
                l_read_name,
                data.len() - cursor
            ),
        ));
    }

    let name_bytes = &data[cursor..cursor + l_read_name];
    // Remove null terminator
    let name = if name_bytes.last() == Some(&0) {
        String::from_utf8(name_bytes[..name_bytes.len() - 1].to_vec()).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid UTF-8 in read name at offset {}: {}", cursor, e),
            )
        })?
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Read name not null-terminated at offset {}", cursor),
        ));
    };
    cursor += l_read_name;

    // CIGAR (4 bytes per operation)
    // Use checked multiplication to prevent overflow
    let cigar_bytes_len = n_cigar_op
        .checked_mul(4)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("CIGAR operation count too large at offset {}: {}", cursor, n_cigar_op),
            )
        })?;

    if cursor + cigar_bytes_len > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Insufficient data for CIGAR at offset {}: need {} bytes for {} operations, got {}",
                cursor,
                cigar_bytes_len,
                n_cigar_op,
                data.len() - cursor
            ),
        ));
    }

    let cigar = parse_cigar(&data[cursor..cursor + cigar_bytes_len], n_cigar_op)?;
    cursor += cigar_bytes_len;

    // Sequence (4-bit encoded, 2 bases per byte)
    let seq_bytes_len = l_seq.div_ceil(2);
    if cursor + seq_bytes_len > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Insufficient data for sequence at offset {}: need {} bytes for {} bases, got {}",
                cursor,
                seq_bytes_len,
                l_seq,
                data.len() - cursor
            ),
        ));
    }

    let sequence = decode_sequence(&data[cursor..cursor + seq_bytes_len], l_seq)?;
    cursor += seq_bytes_len;

    // Quality scores
    if cursor + l_seq > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Insufficient data for quality scores at offset {}: need {} bytes, got {}",
                cursor,
                l_seq,
                data.len() - cursor
            ),
        ));
    }

    // Quality scores (Phred+0, 0xFF = missing)
    // BAM spec: When quality is omitted but sequence is not, qual is filled with 0xFF bytes
    let quality_bytes = &data[cursor..cursor + l_seq];
    let quality = if quality_bytes.iter().all(|&b| b == 0xFF) {
        Vec::new() // Missing quality scores
    } else {
        quality_bytes.to_vec()
    };
    cursor += l_seq;

    // Optional tags (remaining bytes)
    let tags = if cursor < data.len() {
        parse_tags(&data[cursor..])?
    } else {
        Tags::new()
    };

    // Check for oversized CIGAR (>65535 operations) in CG tag
    // Long reads store placeholder 2-op CIGAR (kS*N) and real CIGAR in CG:B,I tag
    let cigar = check_oversized_cigar(cigar, l_seq, &tags)?;

    Ok(Record {
        name,
        reference_id: parse_reference_id(ref_id, "read")?,
        position: if pos >= 0 { Some(pos) } else { None },
        mapq: if mapq != 255 { Some(mapq) } else { None },
        flags,
        mate_reference_id: parse_reference_id(next_ref_id, "mate")?,
        mate_position: if next_pos >= 0 {
            Some(next_pos)
        } else {
            None
        },
        template_length: tlen,
        sequence,
        quality,
        cigar,
        tags,
    })
}

/// Read a single BAM record from a reader.
///
/// Reads the block size, then the record data, then parses.
///
/// # Returns
///
/// - `Ok(Some(record))` - Successfully read a record
/// - `Ok(None)` - EOF (no more records)
/// - `Err(_)` - Parse error
pub fn read_record<R: Read>(reader: &mut R) -> io::Result<Option<Record>> {
    // Read block size
    let mut block_size_bytes = [0u8; 4];
    match reader.read_exact(&mut block_size_bytes) {
        Ok(_) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }

    let block_size = i32::from_le_bytes(block_size_bytes);

    if block_size < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid block size: {}", block_size),
        ));
    }

    let block_size = block_size as usize;

    // Read record data (including the block_size we just read)
    let mut record_data = vec![0u8; block_size + 4];
    record_data[0..4].copy_from_slice(&block_size_bytes);
    reader.read_exact(&mut record_data[4..])?;

    // Parse record
    parse_record(&record_data).map(Some)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_record_flags() {
        let mut record = Record::new();

        // Unmapped
        record.flags = 0x4;
        assert!(record.is_unmapped());
        assert!(!record.is_paired());

        // Paired
        record.flags = 0x1;
        assert!(record.is_paired());
        assert!(!record.is_unmapped());

        // Reverse complement
        record.flags = 0x10;
        assert!(record.is_reverse_complement());
    }

    #[test]
    fn test_parse_minimal_record() {
        let mut data = Vec::new();

        // Block size (32 bytes minimum + variable = let's say 40)
        data.extend_from_slice(&40i32.to_le_bytes());

        // Reference ID (-1 = unmapped)
        data.extend_from_slice(&(-1i32).to_le_bytes());

        // Position (-1 = unmapped)
        data.extend_from_slice(&(-1i32).to_le_bytes());

        // l_read_name (5 = "read" + null)
        data.push(5);

        // MAPQ (255 = unavailable)
        data.push(255);

        // BAI bin (0)
        data.extend_from_slice(&0u16.to_le_bytes());

        // n_cigar_op (0)
        data.extend_from_slice(&0u16.to_le_bytes());

        // FLAGS (4 = unmapped)
        data.extend_from_slice(&4u16.to_le_bytes());

        // l_seq (0)
        data.extend_from_slice(&0i32.to_le_bytes());

        // next_refID (-1)
        data.extend_from_slice(&(-1i32).to_le_bytes());

        // next_pos (-1)
        data.extend_from_slice(&(-1i32).to_le_bytes());

        // tlen (0)
        data.extend_from_slice(&0i32.to_le_bytes());

        // read_name ("read\0")
        data.extend_from_slice(b"read\0");

        // No CIGAR (n_cigar_op = 0)
        // No sequence (l_seq = 0)
        // No quality (l_seq = 0)
        // No tags

        let record = parse_record(&data).unwrap();
        assert_eq!(record.name, "read");
        assert_eq!(record.reference_id, None);
        assert_eq!(record.position, None);
        assert_eq!(record.mapq, None);
        assert_eq!(record.flags, 4);
        assert!(record.is_unmapped());
        assert_eq!(record.sequence.len(), 0);
        assert_eq!(record.quality.len(), 0);
        assert_eq!(record.cigar.len(), 0);
    }

    #[test]
    fn test_parse_record_with_sequence() {
        let mut data = Vec::new();

        // Block size (will calculate)
        let block_size_pos = data.len();
        data.extend_from_slice(&0i32.to_le_bytes()); // Placeholder

        // Reference ID (0)
        data.extend_from_slice(&0i32.to_le_bytes());

        // Position (100)
        data.extend_from_slice(&100i32.to_le_bytes());

        // l_read_name (5)
        data.push(5);

        // MAPQ (60)
        data.push(60);

        // BAI bin (0)
        data.extend_from_slice(&0u16.to_le_bytes());

        // n_cigar_op (1)
        data.extend_from_slice(&1u16.to_le_bytes());

        // FLAGS (0)
        data.extend_from_slice(&0u16.to_le_bytes());

        // l_seq (4)
        data.extend_from_slice(&4i32.to_le_bytes());

        // next_refID (-1)
        data.extend_from_slice(&(-1i32).to_le_bytes());

        // next_pos (-1)
        data.extend_from_slice(&(-1i32).to_le_bytes());

        // tlen (0)
        data.extend_from_slice(&0i32.to_le_bytes());

        // read_name ("test\0")
        data.extend_from_slice(b"test\0");

        // CIGAR (4M = 4 << 4 | 0 = 64)
        data.extend_from_slice(&64u32.to_le_bytes());

        // Sequence ("ACGT" = 0x12, 0x48)
        data.extend_from_slice(&[0x12, 0x48]);

        // Quality (IIII)
        data.extend_from_slice(b"IIII");

        // Update block size
        let block_size = (data.len() - 4) as i32;
        data[block_size_pos..block_size_pos + 4].copy_from_slice(&block_size.to_le_bytes());

        let record = parse_record(&data).unwrap();
        assert_eq!(record.name, "test");
        assert_eq!(record.reference_id, Some(0));
        assert_eq!(record.position, Some(100));
        assert_eq!(record.mapq, Some(60));
        assert_eq!(record.sequence, b"ACGT");
        assert_eq!(record.quality, b"IIII");
        assert_eq!(record.cigar.len(), 1);
    }

    // Malformed input tests

    #[test]
    fn test_integer_overflow_cigar_ops() {
        // Attempt integer overflow: n_cigar_op * 4 would overflow
        let mut data = Vec::new();

        // Block size
        data.extend_from_slice(&100i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(0);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op = u16::MAX (would overflow when * 4)
        data.extend_from_slice(&u16::MAX.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");

        let result = parse_record(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        eprintln!("Actual error: {}", err);
        assert!(err.to_string().contains("CIGAR operation count too large") || err.to_string().contains("Insufficient data"));
    }

    #[test]
    fn test_truncated_record() {
        // Record that claims to need more data than provided
        let mut data = Vec::new();

        // Block size
        data.extend_from_slice(&100i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(0);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op = 10 (but we won't provide the CIGAR data)
        data.extend_from_slice(&10u16.to_le_bytes());

        // Stop here - insufficient data

        let result = parse_record(&data);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_utf8_in_read_name() {
        let mut data = Vec::new();

        // Block size
        data.extend_from_slice(&50i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // l_read_name
        data.push(4);
        // mapq
        data.push(0);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op
        data.extend_from_slice(&0u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name with invalid UTF-8
        data.extend_from_slice(&[0xFF, 0xFE, 0xFD, 0x00]); // Invalid UTF-8 + null terminator

        let result = parse_record(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Invalid UTF-8"));
    }

    #[test]
    fn test_missing_null_terminator_in_name() {
        let mut data = Vec::new();

        // Block size
        data.extend_from_slice(&50i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // l_read_name
        data.push(4);
        // mapq
        data.push(0);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op
        data.extend_from_slice(&0u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name without null terminator
        data.extend_from_slice(b"read");

        let result = parse_record(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("null-terminated") || err.to_string().contains("null terminator"));
    }

    #[test]
    fn test_negative_sequence_length() {
        let mut data = Vec::new();

        // Block size
        data.extend_from_slice(&50i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(0);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op
        data.extend_from_slice(&0u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq = -1 (invalid negative length)
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // next_refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");

        let result = parse_record(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("sequence length") || err.to_string().contains("Invalid"));
    }

    #[test]
    fn test_extremely_large_sequence_length() {
        let mut data = Vec::new();

        // Block size
        data.extend_from_slice(&100i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(0);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op
        data.extend_from_slice(&0u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq = i32::MAX (unreasonably large)
        data.extend_from_slice(&i32::MAX.to_le_bytes());
        // next_refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // next_pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");

        let result = parse_record(&data);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Insufficient data"));
    }

    #[test]
    fn test_zero_read_name_length() {
        let mut data = Vec::new();

        // Block size
        data.extend_from_slice(&50i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&0i32.to_le_bytes());
        // l_read_name = 0 (invalid)
        data.push(0);

        let result = parse_record(&data);
        assert!(result.is_err());
    }

    // Noodles-inspired robustness tests

    #[test]
    fn test_invalid_reference_ids() {
        // Test various invalid reference IDs (only -1 or >= 0 allowed per spec)
        for invalid_id in [-2i32, -3, -100, i32::MIN] {
            let mut data = Vec::new();

            // Block size
            data.extend_from_slice(&50i32.to_le_bytes());
            // refID (invalid)
            data.extend_from_slice(&invalid_id.to_le_bytes());
            // pos
            data.extend_from_slice(&0i32.to_le_bytes());
            // l_read_name
            data.push(5);
            // mapq
            data.push(0);
            // bin
            data.extend_from_slice(&0u16.to_le_bytes());
            // n_cigar_op
            data.extend_from_slice(&0u16.to_le_bytes());
            // flag
            data.extend_from_slice(&0u16.to_le_bytes());
            // l_seq
            data.extend_from_slice(&0i32.to_le_bytes());
            // next_refID
            data.extend_from_slice(&(-1i32).to_le_bytes());
            // next_pos
            data.extend_from_slice(&0i32.to_le_bytes());
            // tlen
            data.extend_from_slice(&0i32.to_le_bytes());
            // read_name
            data.extend_from_slice(b"read\0");

            let result = parse_record(&data);
            assert!(result.is_err(), "Should reject ref_id={}", invalid_id);
            assert!(result.unwrap_err().to_string().contains("Invalid read reference ID"));
        }

        // Also test invalid mate reference IDs
        for invalid_id in [-2i32, -3, -100, i32::MIN] {
            let mut data = Vec::new();

            // Block size
            data.extend_from_slice(&50i32.to_le_bytes());
            // refID (valid)
            data.extend_from_slice(&0i32.to_le_bytes());
            // pos
            data.extend_from_slice(&0i32.to_le_bytes());
            // l_read_name
            data.push(5);
            // mapq
            data.push(0);
            // bin
            data.extend_from_slice(&0u16.to_le_bytes());
            // n_cigar_op
            data.extend_from_slice(&0u16.to_le_bytes());
            // flag
            data.extend_from_slice(&0u16.to_le_bytes());
            // l_seq
            data.extend_from_slice(&0i32.to_le_bytes());
            // next_refID (invalid)
            data.extend_from_slice(&invalid_id.to_le_bytes());
            // next_pos
            data.extend_from_slice(&0i32.to_le_bytes());
            // tlen
            data.extend_from_slice(&0i32.to_le_bytes());
            // read_name
            data.extend_from_slice(b"read\0");

            let result = parse_record(&data);
            assert!(result.is_err(), "Should reject next_ref_id={}", invalid_id);
            assert!(result.unwrap_err().to_string().contains("Invalid mate reference ID"));
        }
    }

    #[test]
    fn test_missing_quality_scores() {
        // Test that all-0xFF quality bytes are detected as missing quality scores
        let mut data = Vec::new();

        // Block size (will be updated)
        data.extend_from_slice(&0i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&100i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(60);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op
        data.extend_from_slice(&0u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&10i32.to_le_bytes());
        // next_refID
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // next_pos
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");
        // sequence (10 bases = 5 bytes)
        data.extend_from_slice(&[0x12, 0x48, 0x12, 0x48, 0x12]);
        // quality (all 0xFF = missing)
        data.extend_from_slice(&[0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF]);

        // Update block size
        let block_size = (data.len() - 4) as i32;
        data[0..4].copy_from_slice(&block_size.to_le_bytes());

        let result = parse_record(&data).expect("Should parse record with missing quality");

        // Quality should be empty when all-0xFF
        assert_eq!(result.quality.len(), 0, "All-0xFF quality should be detected as missing");
        assert_eq!(result.sequence.len(), 10, "Sequence should still be present");
    }

    #[test]
    fn test_present_quality_scores() {
        // Test that non-0xFF quality bytes are preserved
        let mut data = Vec::new();

        // Block size (will be updated)
        data.extend_from_slice(&0i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&100i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(60);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op
        data.extend_from_slice(&0u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&5i32.to_le_bytes());
        // next_refID
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // next_pos
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");
        // sequence (5 bases = 3 bytes)
        data.extend_from_slice(&[0x12, 0x48, 0x10]);
        // quality (real Phred scores)
        data.extend_from_slice(&[30, 35, 40, 38, 32]);

        // Update block size
        let block_size = (data.len() - 4) as i32;
        data[0..4].copy_from_slice(&block_size.to_le_bytes());

        let result = parse_record(&data).expect("Should parse record with quality scores");

        // Quality should be preserved
        assert_eq!(result.quality, vec![30, 35, 40, 38, 32], "Quality scores should be preserved");
        assert_eq!(result.sequence.len(), 5, "Sequence should be present");
    }

    #[test]
    fn test_oversized_cigar_from_cg_tag() {
        // Test that oversized CIGAR (>65535 ops) is extracted from CG:B,I tag
        // Pattern: n_cigar_op=2, first op is kS (k=seq_len), second is *N
        // Real CIGAR stored in CG:B,I tag as Int32 array

        let seq_len = 10;
        let mut data = Vec::new();

        // Block size (will be updated)
        data.extend_from_slice(&0i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&100i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(60);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op = 2 (placeholder for oversized)
        data.extend_from_slice(&2u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&(seq_len as i32).to_le_bytes());
        // next_refID
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // next_pos
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");

        // Placeholder CIGAR: kS*N (k=seq_len=10, S=soft clip op 4, *=100, N=ref skip op 3)
        // Each CIGAR op: low 4 bits = op, high 28 bits = length
        let op1: i32 = (seq_len << 4) | 4; // 10S
        let op2: i32 = (100 << 4) | 3;      // 100N
        data.extend_from_slice(&op1.to_le_bytes());
        data.extend_from_slice(&op2.to_le_bytes());

        // Sequence (10 bases = 5 bytes)
        data.extend_from_slice(&[0x12, 0x48, 0x12, 0x48, 0x12]);
        // Quality (all same)
        data.extend_from_slice(&[30; 10]);

        // CG:B,I tag with real CIGAR: 5M3D5M (3 operations)
        // Tag name: CG
        data.extend_from_slice(b"CG");
        // Tag type: B (array)
        data.push(b'B');
        // Array subtype: i (Int32, lowercase!)
        data.push(b'i');
        // Array count: 3
        data.extend_from_slice(&3u32.to_le_bytes());
        // Array elements (packed CIGAR operations)
        let real_op1: i32 = (5 << 4) | 0; // 5M (match)
        let real_op2: i32 = (3 << 4) | 2; // 3D (deletion)
        let real_op3: i32 = (5 << 4) | 0; // 5M (match)
        data.extend_from_slice(&real_op1.to_le_bytes());
        data.extend_from_slice(&real_op2.to_le_bytes());
        data.extend_from_slice(&real_op3.to_le_bytes());

        // Update block size
        let block_size = (data.len() - 4) as i32;
        data[0..4].copy_from_slice(&block_size.to_le_bytes());

        let result = parse_record(&data).expect("Should parse record with oversized CIGAR");

        // Debug: verify CG tag is present and correct
        let cg_tag = result.tags.get(b"CG").expect("Should parse CG tag").expect("CG tag should exist");
        match cg_tag.value {
            TagValue::Array(ArrayValue::Int32(ref values)) => {
                assert_eq!(values.len(), 3, "CG tag should have 3 values");
                // Verify the packed values
                assert_eq!(values[0], (5 << 4) | 0, "First CG value should be 5M");
                assert_eq!(values[1], (3 << 4) | 2, "Second CG value should be 3D");
                assert_eq!(values[2], (5 << 4) | 0, "Third CG value should be 5M");
            }
            _ => panic!("CG tag should be Int32 array"),
        }

        // Should have real CIGAR from CG tag, not placeholder
        assert_eq!(result.cigar.len(), 3, "Should have 3 CIGAR ops from CG tag");
        assert_eq!(result.cigar[0], CigarOp::Match(5), "First op should be 5M");
        assert_eq!(result.cigar[1], CigarOp::Deletion(3), "Second op should be 3D");
        assert_eq!(result.cigar[2], CigarOp::Match(5), "Third op should be 5M");
    }

    #[test]
    fn test_normal_cigar_not_affected() {
        // Test that normal CIGAR (not oversized pattern) is not affected
        let mut data = Vec::new();

        // Block size (will be updated)
        data.extend_from_slice(&0i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&100i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(60);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op = 2 (but NOT the oversized pattern)
        data.extend_from_slice(&2u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&10i32.to_le_bytes());
        // next_refID
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // next_pos
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");

        // Normal CIGAR: 5M5D (NOT the kS*N pattern)
        let op1: i32 = (5 << 4) | 0; // 5M
        let op2: i32 = (5 << 4) | 2; // 5D
        data.extend_from_slice(&op1.to_le_bytes());
        data.extend_from_slice(&op2.to_le_bytes());

        // Sequence (10 bases = 5 bytes)
        data.extend_from_slice(&[0x12, 0x48, 0x12, 0x48, 0x12]);
        // Quality
        data.extend_from_slice(&[30; 10]);

        // Update block size
        let block_size = (data.len() - 4) as i32;
        data[0..4].copy_from_slice(&block_size.to_le_bytes());

        let result = parse_record(&data).expect("Should parse record with normal CIGAR");

        // Should use original CIGAR, not look for CG tag
        assert_eq!(result.cigar.len(), 2, "Should have 2 CIGAR ops from record");
        assert_eq!(result.cigar[0], CigarOp::Match(5), "First op should be 5M");
        assert_eq!(result.cigar[1], CigarOp::Deletion(5), "Second op should be 5D");
    }

    #[test]
    fn test_oversized_pattern_without_cg_tag() {
        // Test that oversized pattern without CG tag uses placeholder CIGAR
        let seq_len = 10;
        let mut data = Vec::new();

        // Block size (will be updated)
        data.extend_from_slice(&0i32.to_le_bytes());
        // refID
        data.extend_from_slice(&0i32.to_le_bytes());
        // pos
        data.extend_from_slice(&100i32.to_le_bytes());
        // l_read_name
        data.push(5);
        // mapq
        data.push(60);
        // bin
        data.extend_from_slice(&0u16.to_le_bytes());
        // n_cigar_op = 2 (placeholder for oversized)
        data.extend_from_slice(&2u16.to_le_bytes());
        // flag
        data.extend_from_slice(&0u16.to_le_bytes());
        // l_seq
        data.extend_from_slice(&(seq_len as i32).to_le_bytes());
        // next_refID
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // next_pos
        data.extend_from_slice(&(-1i32).to_le_bytes());
        // tlen
        data.extend_from_slice(&0i32.to_le_bytes());
        // read_name
        data.extend_from_slice(b"read\0");

        // Placeholder CIGAR: kS*N (oversized pattern)
        let op1: i32 = (seq_len << 4) | 4; // 10S
        let op2: i32 = (100 << 4) | 3;      // 100N
        data.extend_from_slice(&op1.to_le_bytes());
        data.extend_from_slice(&op2.to_le_bytes());

        // Sequence (10 bases = 5 bytes)
        data.extend_from_slice(&[0x12, 0x48, 0x12, 0x48, 0x12]);
        // Quality
        data.extend_from_slice(&[30; 10]);

        // No CG tag - just use placeholder

        // Update block size
        let block_size = (data.len() - 4) as i32;
        data[0..4].copy_from_slice(&block_size.to_le_bytes());

        let result = parse_record(&data).expect("Should parse record even without CG tag");

        // Should use placeholder CIGAR when pattern matched but no CG tag
        assert_eq!(result.cigar.len(), 2, "Should have 2 CIGAR ops (placeholder)");
        assert_eq!(result.cigar[0], CigarOp::SoftClip(10), "First op should be 10S");
        assert_eq!(result.cigar[1], CigarOp::RefSkip(100), "Second op should be 100N");
    }
}
