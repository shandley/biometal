//! ARM NEON-optimized BAM sequence decoding.
//!
//! Implements 4-bit to ASCII decoding using ARM NEON SIMD instructions.
//! Expected 16-25× speedup vs scalar implementation (Rule 1).
//!
//! # Evidence Base
//!
//! - Sequence decoding is **30.2%** of total BAM parsing time (validated Nov 9, 2025)
//! - NEON optimization expected to provide **+38-40% overall BAM parsing speedup**
//! - Rule 1: Element-wise operations achieve 16-25× speedup with NEON
//!
//! See: `experiments/bam-simd-sequence-decoding/RESEARCH_LOG.md`

use std::arch::aarch64::*;
use std::io;

/// Decode a 4-bit encoded BAM sequence to ASCII using ARM NEON.
///
/// This function uses NEON SIMD instructions to decode BAM's 4-bit sequence
/// encoding (2 bases per byte) to standard ASCII bases in parallel.
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
/// Expected: 16-25× faster than scalar implementation (Rule 1)
/// - Processes 32 bases per NEON iteration (16-byte chunks)
/// - Uses `vqtbl1q_u8` for 16-entry table lookup
/// - Interleaves high/low nibbles with `vzip` instructions
///
/// # Example
///
/// ```ignore
/// use biometal::io::bam::sequence_neon::decode_sequence_neon;
///
/// // Byte 0x12 encodes bases at indices 1 ('A') and 2 ('C')
/// let data = vec![0x12];
/// let sequence = decode_sequence_neon(&data, 2).unwrap();
/// assert_eq!(sequence, b"AC");
/// ```
#[cfg(target_arch = "aarch64")]
#[inline]
pub fn decode_sequence_neon(data: &[u8], length: usize) -> io::Result<Vec<u8>> {
    // Check if we have enough data
    let required_bytes = length.div_ceil(2);
    if data.len() < required_bytes {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Insufficient sequence data: need {} bytes for {} bases, got {}",
                required_bytes,
                length,
                data.len()
            ),
        ));
    }

    let mut sequence: Vec<u8> = Vec::with_capacity(length);

    // BAM 4-bit to ASCII base lookup table (16 bytes)
    // Encoding from SAM/BAM specification v1.6
    // 16-byte alignment for optimal NEON load performance.
    // While vld1q_u8 handles unaligned data correctly, alignment
    // avoids potential performance penalty on older ARM cores.
    #[repr(align(16))]
    struct AlignedLookup([u8; 16]);

    let lookup_data = AlignedLookup([
        b'=', b'A', b'C', b'M', // 0-3
        b'G', b'R', b'S', b'V', // 4-7
        b'T', b'W', b'Y', b'H', // 8-11
        b'K', b'D', b'B', b'N', // 12-15
    ]);

    // Process 16 bytes (32 bases) at a time
    let full_chunks = length / 32;
    let mut offset = 0;

    // SAFETY: All NEON operations are safe because:
    // - Bounds checked before any pointer operations
    // - All pointer arithmetic stays within allocated capacity
    // - set_len() only called after data is written
    // - NEON instructions available on all ARM64 platforms
    unsafe {
        // Load lookup table into NEON register
        let lookup_table = vld1q_u8(lookup_data.0.as_ptr());

        // Constant for extracting low nibbles
        let nibble_mask = vdupq_n_u8(0x0F);

        for _ in 0..full_chunks {
            // Load 16 bytes (32 bases worth of packed data)
            let packed = vld1q_u8(data.as_ptr().add(offset));

            // Extract high nibbles (first 16 bases)
            let high_nibbles = vshrq_n_u8::<4>(packed);

            // Extract low nibbles (second 16 bases)
            let low_nibbles = vandq_u8(packed, nibble_mask);

            // Lookup bases for both high and low nibbles
            let bases_high = vqtbl1q_u8(lookup_table, high_nibbles);
            let bases_low = vqtbl1q_u8(lookup_table, low_nibbles);

            // Interleave high/low nibbles into final sequence order
            // vzip1q_u8: Interleaves lower 64 bits → [h0,l0,h1,l1,...,h7,l7]
            // vzip2q_u8: Interleaves upper 64 bits → [h8,l8,h9,l9,...,h15,l15]
            // Combined: 32 bases in correct nibble order
            let interleaved = vzip1q_u8(bases_high, bases_low);
            let interleaved2 = vzip2q_u8(bases_high, bases_low);

            // Store first 16 interleaved bases
            let dest_ptr = sequence.as_mut_ptr().add(sequence.len());
            vst1q_u8(dest_ptr, interleaved);
            sequence.set_len(sequence.len() + 16);

            // Store second 16 interleaved bases
            let dest_ptr2 = sequence.as_mut_ptr().add(sequence.len());
            vst1q_u8(dest_ptr2, interleaved2);
            sequence.set_len(sequence.len() + 16);

            offset += 16;
        }
    } // end unsafe block

    // Handle remaining bases (< 32)
    let remaining = length - (full_chunks * 32);

    if remaining > 0 {
        // Process remaining bases using scalar fallback
        // This is simpler than complex NEON masking for the tail
        for i in 0..remaining {
            let global_idx = full_chunks * 32 + i;
            let byte_idx = global_idx / 2;
            let nibble = if global_idx.is_multiple_of(2) {
                // High nibble
                data[byte_idx] >> 4
            } else {
                // Low nibble
                data[byte_idx] & 0x0F
            };

            // Use lookup table for consistency
            let base = lookup_data.0[nibble as usize];
            sequence.push(base);
        }
    }

    Ok(sequence)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_single_base() {
        // 0x10 = high nibble 1 ('A'), low nibble 0 ('=')
        let data = vec![0x10];
        let seq = decode_sequence_neon(&data, 1).unwrap();
        assert_eq!(seq, b"A");
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_two_bases() {
        // 0x12 = 'A' and 'C'
        let data = vec![0x12];
        let seq = decode_sequence_neon(&data, 2).unwrap();
        assert_eq!(seq, b"AC");
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_acgt() {
        // Test all 4 standard bases
        // 0x12 = 'A' (1) + 'C' (2)
        // 0x48 = 'G' (4) + 'T' (8)
        let data = vec![0x12, 0x48];
        let seq = decode_sequence_neon(&data, 4).unwrap();
        assert_eq!(seq, b"ACGT");
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_32_bases() {
        // Test exactly one NEON iteration (32 bases = 16 bytes)
        let data = vec![0x12; 16]; // All 'A' and 'C'
        let seq = decode_sequence_neon(&data, 32).unwrap();
        assert_eq!(seq, b"ACACACACACACACACACACACACACACACAC");
        assert_eq!(seq.len(), 32);
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_64_bases() {
        // Test exactly two NEON iterations (64 bases = 32 bytes)
        let data = vec![0x48; 32]; // All 'G' and 'T'
        let seq = decode_sequence_neon(&data, 64).unwrap();
        let expected = "GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT";
        assert_eq!(seq, expected.as_bytes());
        assert_eq!(seq.len(), 64);
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_100_bases() {
        // Test NEON path + scalar tail (100 = 96 NEON + 4 scalar)
        let data = vec![0x12; 50]; // All 'A' and 'C'
        let seq = decode_sequence_neon(&data, 100).unwrap();
        assert_eq!(seq.len(), 100);
        assert!(seq.iter().all(|&b| b == b'A' || b == b'C'));
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_odd_length() {
        // Odd number of bases (5 bases = 3 bytes)
        // 0x12, 0x48, 0x10 = A, C, G, T, A
        let data = vec![0x12, 0x48, 0x10];
        let seq = decode_sequence_neon(&data, 5).unwrap();
        assert_eq!(seq, b"ACGTA");
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_ambiguity_codes() {
        // Test some IUPAC ambiguity codes
        // N (any base) = 15 = 0xF
        // 0xFF = N + N
        let data = vec![0xFF];
        let seq = decode_sequence_neon(&data, 2).unwrap();
        assert_eq!(seq, b"NN");
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_decode_empty_sequence() {
        let data = vec![];
        let seq = decode_sequence_neon(&data, 0).unwrap();
        assert_eq!(seq, b"");
        assert_eq!(seq.len(), 0);
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_all_lookup_values() {
        // Verify all 16 lookup table entries
        let mut data = Vec::new();
        for i in 0..16 {
            data.push((i << 4) | i); // Same value in high and low nibble
        }

        let seq = decode_sequence_neon(&data, 32).unwrap();
        let expected = b"==AACCMMGGRRSSVVTTWWYYHHKKDDBBNN";

        // Check that we get pairs of each character
        for i in 0..16 {
            assert_eq!(seq[i * 2], expected[i * 2]);
            assert_eq!(seq[i * 2 + 1], expected[i * 2 + 1]);
        }
    }
}
