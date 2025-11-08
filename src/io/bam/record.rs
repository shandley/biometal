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
use super::tags::{parse_tags, Tags};
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

    let quality = data[cursor..cursor + l_seq].to_vec();
    cursor += l_seq;

    // Optional tags (remaining bytes)
    let tags = if cursor < data.len() {
        parse_tags(&data[cursor..])?
    } else {
        Tags::new()
    };

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
}
