//! 4-bit BAM sequence decoding.
//!
//! BAM stores sequences in 4-bit encoding (2 bases per byte) to save space.
//! Each nibble (4 bits) encodes one base using a lookup table.
//!
//! # Format
//!
//! - High nibble first, low nibble second within each byte
//! - 16 possible values (4 bases + ambiguity codes)
//! - Example: byte 0x12 -> bases at indices 1 ('A') and 2 ('C')
//!
//! # Evidence Base
//!
//! Phase 0 profiling showed sequence decoding is <6% CPU time,
//! so NEON optimization is deferred until profiling shows >=15% (Rule 1).
//! This implementation uses a simple scalar lookup table.
//!
//! See: `experiments/native-bam-implementation/DECISION_REVISED.md`

use std::io;

/// BAM 4-bit to ASCII base lookup table.
///
/// Encoding from SAM/BAM specification v1.6:
/// - 0 = '=' (match to reference, not used in practice)
/// - 1-4 = A, C, G, T
/// - 5-15 = IUPAC ambiguity codes (M, R, S, V, W, Y, H, K, D, B, N)
const SEQ_LOOKUP: [u8; 16] = [
    b'=', b'A', b'C', b'M', // 0-3
    b'G', b'R', b'S', b'V', // 4-7
    b'T', b'W', b'Y', b'H', // 8-11
    b'K', b'D', b'B', b'N', // 12-15
];

/// Decode a 4-bit encoded BAM sequence to ASCII.
///
/// BAM sequences are stored as 2 bases per byte (high nibble first).
/// This function unpacks the 4-bit encoding to standard ASCII bases.
///
/// # Arguments
///
/// * `data` - Packed 4-bit sequence data (ceil(length/2) bytes)
/// * `length` - Number of bases to decode
///
/// # Returns
///
/// Vector of ASCII bases (A, C, G, T, N, etc.)
///
/// # Errors
///
/// Returns error if `data` is too short for the specified `length`.
///
/// # Performance
///
/// This is a scalar implementation. NEON optimization deferred until
/// profiling shows sequence decoding is e15% CPU time (currently <6%).
///
/// # Example
///
/// ```
/// use biometal::io::bam::decode_sequence;
///
/// // Byte 0x12 encodes bases at indices 1 ('A') and 2 ('C')
/// let data = vec![0x12];
/// let sequence = decode_sequence(&data, 2).unwrap();
/// assert_eq!(sequence, b"AC");
/// ```
pub fn decode_sequence(data: &[u8], length: usize) -> io::Result<Vec<u8>> {
    // Check if we have enough data
    let required_bytes = length.div_ceil(2);
    if data.len() < required_bytes {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Insufficient sequence data: need {} bytes for {} bases, got {} (at offset 0)",
                required_bytes,
                length,
                data.len()
            ),
        ));
    }

    let mut sequence = Vec::with_capacity(length);

    for i in 0..length {
        let byte_idx = i / 2;
        let nibble = if i % 2 == 0 {
            // High nibble (first base in byte)
            data[byte_idx] >> 4
        } else {
            // Low nibble (second base in byte)
            data[byte_idx] & 0x0F
        };

        sequence.push(SEQ_LOOKUP[nibble as usize]);
    }

    Ok(sequence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn test_decode_single_base() {
        // 0x10 = high nibble 1 ('A'), low nibble 0 ('=')
        let data = vec![0x10];
        let seq = decode_sequence(&data, 1).unwrap();
        assert_eq!(seq, b"A");
    }

    #[test]
    fn test_decode_two_bases() {
        // 0x12 = 'A' and 'C'
        let data = vec![0x12];
        let seq = decode_sequence(&data, 2).unwrap();
        assert_eq!(seq, b"AC");
    }

    #[test]
    fn test_decode_acgt() {
        // Test all 4 standard bases
        // 0x12 = 'A' (1) + 'C' (2)
        // 0x48 = 'G' (4) + 'T' (8)
        let data = vec![0x12, 0x48];
        let seq = decode_sequence(&data, 4).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn test_decode_odd_length() {
        // Odd number of bases (5 bases = 3 bytes)
        // 0x12, 0x48, 0x10 = A, C, G, T, A
        let data = vec![0x12, 0x48, 0x10];
        let seq = decode_sequence(&data, 5).unwrap();
        assert_eq!(seq, b"ACGTA");
    }

    #[test]
    fn test_decode_ambiguity_codes() {
        // Test some IUPAC ambiguity codes
        // N (any base) = 15 = 0xF
        // 0xFF = N + N
        let data = vec![0xFF];
        let seq = decode_sequence(&data, 2).unwrap();
        assert_eq!(seq, b"NN");
    }

    #[test]
    fn test_decode_empty_sequence() {
        let data = vec![];
        let seq = decode_sequence(&data, 0).unwrap();
        assert_eq!(seq, b"");
    }

    #[test]
    fn test_insufficient_data_error() {
        let data = vec![0x12]; // Only 1 byte
        let result = decode_sequence(&data, 5); // Need 3 bytes for 5 bases
        assert!(result.is_err());
    }

    #[test]
    fn test_all_lookup_values() {
        // Verify all 16 lookup table entries
        let expected = b"=ACMGRSVTWYHKDBN";
        for i in 0..16 {
            assert_eq!(SEQ_LOOKUP[i], expected[i]);
        }
    }

    // Property-based tests

    /// Helper to encode a sequence to 4-bit format
    fn encode_sequence(bases: &[u8]) -> Vec<u8> {
        let mut encoded = Vec::with_capacity(bases.len().div_ceil(2));

        for chunk in bases.chunks(2) {
            let high = base_to_nibble(chunk[0]);
            let low = if chunk.len() > 1 {
                base_to_nibble(chunk[1])
            } else {
                0 // Padding for odd-length sequences
            };
            encoded.push((high << 4) | low);
        }

        encoded
    }

    /// Helper to convert ASCII base to 4-bit nibble
    fn base_to_nibble(base: u8) -> u8 {
        match base {
            b'=' => 0,
            b'A' => 1,
            b'C' => 2,
            b'M' => 3,
            b'G' => 4,
            b'R' => 5,
            b'S' => 6,
            b'V' => 7,
            b'T' => 8,
            b'W' => 9,
            b'Y' => 10,
            b'H' => 11,
            b'K' => 12,
            b'D' => 13,
            b'B' => 14,
            b'N' => 15,
            _ => 15, // Default to N for invalid bases
        }
    }

    proptest! {
        #[test]
        fn prop_sequence_roundtrip_acgt(sequence in "[ACGT]{1,1000}") {
            // Convert string to bytes
            let bases = sequence.as_bytes();
            let length = bases.len();

            // Encode to 4-bit
            let encoded = encode_sequence(bases);

            // Decode back
            let decoded = decode_sequence(&encoded, length).unwrap();

            // Should match original
            prop_assert_eq!(decoded, bases);
        }

        #[test]
        fn prop_sequence_roundtrip_all_bases(sequence in "[=ACMGRSVTWYHKDBN]{1,100}") {
            // Convert string to bytes
            let bases = sequence.as_bytes();
            let length = bases.len();

            // Encode to 4-bit
            let encoded = encode_sequence(bases);

            // Decode back
            let decoded = decode_sequence(&encoded, length).unwrap();

            // Should match original
            prop_assert_eq!(decoded, bases);
        }

        #[test]
        fn prop_sequence_length_invariant(length in 1usize..1000) {
            // Generate sequence of all A's
            let sequence = "A".repeat(length);
            let bases = sequence.as_bytes();

            // Encode
            let encoded = encode_sequence(bases);

            // Decoded length should match original length
            let decoded = decode_sequence(&encoded, length).unwrap();
            prop_assert_eq!(decoded.len(), length);
        }

        #[test]
        fn prop_sequence_odd_even_length(length in 1usize..500) {
            // Test both odd and even lengths
            for &len in &[length, length + 1] {
                let sequence = "G".repeat(len);
                let bases = sequence.as_bytes();
                let encoded = encode_sequence(bases);
                let decoded = decode_sequence(&encoded, len).unwrap();

                prop_assert_eq!(decoded.len(), len);
                prop_assert!(decoded.iter().all(|&b| b == b'G'));
            }
        }

        #[test]
        fn prop_sequence_encoded_size(length in 1usize..1000) {
            let sequence = "C".repeat(length);
            let bases = sequence.as_bytes();
            let encoded = encode_sequence(bases);

            // Encoded size should be ceil(length / 2)
            let expected_size = length.div_ceil(2);
            prop_assert_eq!(encoded.len(), expected_size);

            // Should be able to decode
            let decoded = decode_sequence(&encoded, length).unwrap();
            prop_assert_eq!(decoded.len(), length);
        }
    }
}
