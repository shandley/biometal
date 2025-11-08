//! CIGAR (Compact Idiosyncratic Gapped Alignment Report) parsing.
//!
//! CIGAR strings describe how a read aligns to the reference, including
//! matches, insertions, deletions, and other operations.
//!
//! # BAM Format
//!
//! In BAM, CIGAR is stored as 32-bit integers:
//! - Low 4 bits: operation type (0-8)
//! - High 28 bits: operation length (0 to 268,435,455)
//!
//! # Operations
//!
//! - M: Match/mismatch (alignment match, can include mismatches)
//! - I: Insertion to reference
//! - D: Deletion from reference
//! - N: Skipped region from reference (intron for RNA-seq)
//! - S: Soft clipping (bases present in read, not in alignment)
//! - H: Hard clipping (bases not present in read)
//! - P: Padding (silent deletion from padded reference)
//! - =: Sequence match (bases match reference)
//! - X: Sequence mismatch (bases don't match reference)

use std::io;

/// CIGAR operation types.
///
/// Each operation describes a type of alignment event and its length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// Match or mismatch (M)
    Match(u32),
    /// Insertion to reference (I)
    Insertion(u32),
    /// Deletion from reference (D)
    Deletion(u32),
    /// Skipped region from reference (N)
    RefSkip(u32),
    /// Soft clipping (S)
    SoftClip(u32),
    /// Hard clipping (H)
    HardClip(u32),
    /// Padding (P)
    Padding(u32),
    /// Sequence match (=)
    SeqMatch(u32),
    /// Sequence mismatch (X)
    SeqMismatch(u32),
}

impl CigarOp {
    /// Get the operation count/length.
    ///
    /// Returns the number of bases consumed by this CIGAR operation.
    pub fn length(&self) -> u32 {
        match self {
            CigarOp::Match(len) => *len,
            CigarOp::Insertion(len) => *len,
            CigarOp::Deletion(len) => *len,
            CigarOp::RefSkip(len) => *len,
            CigarOp::SoftClip(len) => *len,
            CigarOp::HardClip(len) => *len,
            CigarOp::Padding(len) => *len,
            CigarOp::SeqMatch(len) => *len,
            CigarOp::SeqMismatch(len) => *len,
        }
    }

    /// Check if this operation has zero length.
    pub fn is_empty(&self) -> bool {
        self.length() == 0
    }

    /// Get the operation type as a character (for SAM format).
    pub fn as_char(&self) -> char {
        match self {
            CigarOp::Match(_) => 'M',
            CigarOp::Insertion(_) => 'I',
            CigarOp::Deletion(_) => 'D',
            CigarOp::RefSkip(_) => 'N',
            CigarOp::SoftClip(_) => 'S',
            CigarOp::HardClip(_) => 'H',
            CigarOp::Padding(_) => 'P',
            CigarOp::SeqMatch(_) => '=',
            CigarOp::SeqMismatch(_) => 'X',
        }
    }
}

impl std::fmt::Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.length(), self.as_char())
    }
}

/// Parse BAM CIGAR operations from binary format.
///
/// Each CIGAR operation is encoded as a 32-bit little-endian integer:
/// - Bits 0-3: operation type (0-8)
/// - Bits 4-31: operation length
///
/// # Arguments
///
/// * `data` - Raw CIGAR data (4 bytes per operation)
/// * `n_ops` - Number of CIGAR operations
///
/// # Returns
///
/// Vector of CIGAR operations.
///
/// # Errors
///
/// Returns error if:
/// - Data is too short for specified number of operations
/// - Invalid operation code encountered (not 0-8)
///
/// # Example
///
/// ```
/// use biometal::io::bam::{parse_cigar, CigarOp};
///
/// // 100M = 100 << 4 | 0 = 1600 = 0x00000640 (little-endian)
/// let data = vec![0x40, 0x06, 0x00, 0x00];
/// let cigar = parse_cigar(&data, 1).unwrap();
/// assert_eq!(cigar.len(), 1);
/// assert_eq!(cigar[0], CigarOp::Match(100));
/// ```
pub fn parse_cigar(data: &[u8], n_ops: usize) -> io::Result<Vec<CigarOp>> {
    // Check if we have enough data
    let required_bytes = n_ops * 4;
    if data.len() < required_bytes {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Insufficient CIGAR data: need {} bytes for {} operations, got {}",
                required_bytes,
                n_ops,
                data.len()
            ),
        ));
    }

    let mut ops = Vec::with_capacity(n_ops);

    for i in 0..n_ops {
        let offset = i * 4;

        // Read 32-bit little-endian integer
        let cigar_int = u32::from_le_bytes([
            data[offset],
            data[offset + 1],
            data[offset + 2],
            data[offset + 3],
        ]);

        // Extract operation length (high 28 bits) and type (low 4 bits)
        let length = cigar_int >> 4;
        let op_code = cigar_int & 0x0F;

        // Map operation code to CigarOp
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
                    format!("Invalid CIGAR operation code: {}", op_code),
                ))
            }
        };

        ops.push(cigar_op);
    }

    Ok(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn test_cigar_op_length() {
        assert_eq!(CigarOp::Match(100).length(), 100);
        assert_eq!(CigarOp::Insertion(5).length(), 5);
    }

    #[test]
    fn test_cigar_op_char() {
        assert_eq!(CigarOp::Match(100).as_char(), 'M');
        assert_eq!(CigarOp::Insertion(5).as_char(), 'I');
        assert_eq!(CigarOp::Deletion(3).as_char(), 'D');
    }

    #[test]
    fn test_cigar_op_display() {
        assert_eq!(format!("{}", CigarOp::Match(100)), "100M");
        assert_eq!(format!("{}", CigarOp::Insertion(5)), "5I");
    }

    #[test]
    fn test_parse_single_match() {
        // 100M = 100 << 4 | 0 = 1600 = 0x00000640
        let data = vec![0x40, 0x06, 0x00, 0x00];
        let cigar = parse_cigar(&data, 1).unwrap();
        assert_eq!(cigar.len(), 1);
        assert_eq!(cigar[0], CigarOp::Match(100));
    }

    #[test]
    fn test_parse_multiple_operations() {
        // 50M 5I 45M
        // 50M = 50 << 4 | 0 = 800 = 0x00000320
        // 5I = 5 << 4 | 1 = 81 = 0x00000051
        // 45M = 45 << 4 | 0 = 720 = 0x000002D0
        let data = vec![
            0x20, 0x03, 0x00, 0x00, // 50M
            0x51, 0x00, 0x00, 0x00, // 5I
            0xD0, 0x02, 0x00, 0x00, // 45M
        ];
        let cigar = parse_cigar(&data, 3).unwrap();
        assert_eq!(cigar.len(), 3);
        assert_eq!(cigar[0], CigarOp::Match(50));
        assert_eq!(cigar[1], CigarOp::Insertion(5));
        assert_eq!(cigar[2], CigarOp::Match(45));
    }

    #[test]
    fn test_parse_all_operations() {
        // Test all 9 operation types
        let data = vec![
            0x10, 0x00, 0x00, 0x00, // 1M
            0x11, 0x00, 0x00, 0x00, // 1I
            0x12, 0x00, 0x00, 0x00, // 1D
            0x13, 0x00, 0x00, 0x00, // 1N
            0x14, 0x00, 0x00, 0x00, // 1S
            0x15, 0x00, 0x00, 0x00, // 1H
            0x16, 0x00, 0x00, 0x00, // 1P
            0x17, 0x00, 0x00, 0x00, // 1=
            0x18, 0x00, 0x00, 0x00, // 1X
        ];
        let cigar = parse_cigar(&data, 9).unwrap();
        assert_eq!(cigar.len(), 9);
        assert_eq!(cigar[0], CigarOp::Match(1));
        assert_eq!(cigar[1], CigarOp::Insertion(1));
        assert_eq!(cigar[2], CigarOp::Deletion(1));
        assert_eq!(cigar[3], CigarOp::RefSkip(1));
        assert_eq!(cigar[4], CigarOp::SoftClip(1));
        assert_eq!(cigar[5], CigarOp::HardClip(1));
        assert_eq!(cigar[6], CigarOp::Padding(1));
        assert_eq!(cigar[7], CigarOp::SeqMatch(1));
        assert_eq!(cigar[8], CigarOp::SeqMismatch(1));
    }

    #[test]
    fn test_parse_empty_cigar() {
        let data = vec![];
        let cigar = parse_cigar(&data, 0).unwrap();
        assert_eq!(cigar.len(), 0);
    }

    #[test]
    fn test_insufficient_data_error() {
        let data = vec![0x10, 0x00]; // Only 2 bytes
        let result = parse_cigar(&data, 1); // Need 4 bytes
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_operation_code() {
        // Operation code 9 is invalid
        let data = vec![0x19, 0x00, 0x00, 0x00]; // 1 (invalid op 9)
        let result = parse_cigar(&data, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_large_length() {
        // Test maximum length (28 bits = 268,435,455)
        // 0xFFFFFFF0 = length 268435455, op 0 (Match)
        let data = vec![0xF0, 0xFF, 0xFF, 0xFF];
        let cigar = parse_cigar(&data, 1).unwrap();
        assert_eq!(cigar[0], CigarOp::Match(268_435_455));
    }

    // Property-based tests

    /// Helper to encode a CIGAR operation to binary format
    fn encode_cigar_op(op: &CigarOp) -> [u8; 4] {
        let length = op.length();
        let op_code = match op {
            CigarOp::Match(_) => 0,
            CigarOp::Insertion(_) => 1,
            CigarOp::Deletion(_) => 2,
            CigarOp::RefSkip(_) => 3,
            CigarOp::SoftClip(_) => 4,
            CigarOp::HardClip(_) => 5,
            CigarOp::Padding(_) => 6,
            CigarOp::SeqMatch(_) => 7,
            CigarOp::SeqMismatch(_) => 8,
        };
        let cigar_int = (length << 4) | op_code;
        cigar_int.to_le_bytes()
    }

    proptest! {
        #[test]
        fn prop_cigar_roundtrip(length in 1u32..10000u32, op_code in 0u32..=8u32) {
            // Create CIGAR op from length and op code
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
                _ => unreachable!(),
            };

            // Encode to binary
            let encoded = encode_cigar_op(&cigar_op);

            // Parse back
            let parsed = parse_cigar(&encoded, 1).unwrap();

            // Should roundtrip correctly
            prop_assert_eq!(parsed.len(), 1);
            prop_assert_eq!(parsed[0], cigar_op);
            prop_assert_eq!(parsed[0].length(), length);
        }

        #[test]
        fn prop_cigar_multiple_ops(
            ops in prop::collection::vec(0u32..=8u32, 1..10),
            lengths in prop::collection::vec(1u32..=1000u32, 1..10)
        ) {
            // Ensure ops and lengths have same length
            let ops = &ops[..ops.len().min(lengths.len())];
            let lengths = &lengths[..ops.len()];

            // Create CIGAR ops
            let cigar_ops: Vec<CigarOp> = ops.iter().zip(lengths.iter())
                .map(|(op_code, length)| match op_code {
                    0 => CigarOp::Match(*length),
                    1 => CigarOp::Insertion(*length),
                    2 => CigarOp::Deletion(*length),
                    3 => CigarOp::RefSkip(*length),
                    4 => CigarOp::SoftClip(*length),
                    5 => CigarOp::HardClip(*length),
                    6 => CigarOp::Padding(*length),
                    7 => CigarOp::SeqMatch(*length),
                    8 => CigarOp::SeqMismatch(*length),
                    _ => unreachable!(),
                })
                .collect();

            // Encode all ops
            let mut encoded = Vec::new();
            for op in &cigar_ops {
                encoded.extend_from_slice(&encode_cigar_op(op));
            }

            // Parse back
            let parsed = parse_cigar(&encoded, cigar_ops.len()).unwrap();

            // Should match
            prop_assert_eq!(parsed.len(), cigar_ops.len());
            for (parsed_op, original_op) in parsed.iter().zip(cigar_ops.iter()) {
                prop_assert_eq!(parsed_op, original_op);
            }
        }

        #[test]
        fn prop_cigar_display_format(length in 1u32..10000u32, op_code in 0u32..=8u32) {
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
                _ => unreachable!(),
            };

            let display_str = format!("{}", cigar_op);

            // Should start with length
            prop_assert!(display_str.starts_with(&length.to_string()));

            // Should end with correct character
            let expected_char = cigar_op.as_char();
            prop_assert!(display_str.ends_with(expected_char));
        }
    }
}
