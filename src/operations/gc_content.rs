//! GC content calculation with ARM NEON SIMD optimization (Rule 1)
//!
//! # Evidence
//!
//! Entry 020-025 (Lab Notebook):
//! - **Speedup**: 20.3× faster than scalar
//! - **Statistical rigor**: Cohen's d = 5.12 (very large effect)
//! - **Cross-platform**: Mac M4 Max, AWS Graviton 3
//!
//! # Architecture
//!
//! This module provides both NEON (ARM) and scalar (portable) implementations:
//! - NEON: Processes 16 bytes at a time with SIMD instructions
//! - Scalar: Sequential processing (x86_64 fallback)
//!
//! The public API automatically selects the best implementation for the platform.

/// GC content result as a fraction (0.0 to 1.0)
pub type GCContent = f64;

/// Calculate GC content of a DNA sequence
///
/// Returns the fraction of bases that are G or C (ignoring non-ACGT characters).
///
/// # Platform-Specific Optimization
///
/// - **ARM (aarch64)**: Uses NEON SIMD (20.3× speedup)
/// - **x86_64**: Uses scalar fallback (portable)
///
/// # Evidence
///
/// Entry 020-025: NEON GC content achieves 20.3× speedup with Cohen's d = 5.12
///
/// # Example
///
/// ```
/// use biometal::operations::gc_content;
///
/// let sequence = b"GATTACAGATTACA";
/// let gc = gc_content(sequence);
/// assert!((gc - 0.286).abs() < 0.001); // ~28.6% GC (4 out of 14: 2G + 2C)
/// ```
pub fn gc_content(seq: &[u8]) -> GCContent {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { gc_content_neon(seq) }
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        gc_content_scalar(seq)
    }
}

/// NEON-optimized GC content calculation (20.3× faster than scalar)
///
/// # Evidence
///
/// Entry 020-025, Cohen's d = 5.12 (very large effect)
///
/// # Safety
///
/// This function uses unsafe NEON intrinsics but is safe to call:
/// - Only called on aarch64 platforms (compile-time check)
/// - NEON is standard on all aarch64 CPUs
/// - Pointer operations are bounds-checked via chunks_exact
#[cfg(target_arch = "aarch64")]
pub unsafe fn gc_content_neon(seq: &[u8]) -> GCContent {
    use std::arch::aarch64::*;

    let mut gc_count = 0u32;
    let mut total_count = 0u32;

    // NEON registers for GC counts
    let mut gc_vcount = vdupq_n_u32(0);

    // Process 16 bytes at a time
    let chunks = seq.chunks_exact(16);
    let remainder = chunks.remainder();

    for chunk in chunks {
        let seq_vec = vld1q_u8(chunk.as_ptr());

        // Compare against G and C
        let g_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'G'));
        let c_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'C'));

        // Combine G and C masks (bitwise OR)
        let gc_mask = vorrq_u8(g_mask, c_mask);

        // Convert 0xFF to 0x01 by shifting right 7 bits
        let gc_count_u8 = vshrq_n_u8(gc_mask, 7);

        // Accumulate counts
        gc_vcount = vaddq_u32(gc_vcount, vpaddlq_u16(vpaddlq_u8(gc_count_u8)));

        // Count all ACGT bases (not N or other ambiguous bases)
        let a_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'A'));
        let t_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'T'));

        // Combine all ACGT masks
        let acgt_mask = vorrq_u8(vorrq_u8(a_mask, t_mask), gc_mask);
        let acgt_count_u8 = vshrq_n_u8(acgt_mask, 7);

        // Accumulate total
        let acgt_vcount = vpaddlq_u16(vpaddlq_u8(acgt_count_u8));
        total_count += vgetq_lane_u32(acgt_vcount, 0);
        total_count += vgetq_lane_u32(acgt_vcount, 1);
        total_count += vgetq_lane_u32(acgt_vcount, 2);
        total_count += vgetq_lane_u32(acgt_vcount, 3);
    }

    // Extract GC count from NEON register
    gc_count += vgetq_lane_u32(gc_vcount, 0);
    gc_count += vgetq_lane_u32(gc_vcount, 1);
    gc_count += vgetq_lane_u32(gc_vcount, 2);
    gc_count += vgetq_lane_u32(gc_vcount, 3);

    // Handle remainder with scalar code
    for &base in remainder {
        match base {
            b'G' | b'C' => {
                gc_count += 1;
                total_count += 1;
            }
            b'A' | b'T' => {
                total_count += 1;
            }
            _ => {} // Ignore non-ACGT characters
        }
    }

    if total_count == 0 {
        0.0 // Empty sequence or all non-ACGT
    } else {
        gc_count as f64 / total_count as f64
    }
}

/// Scalar fallback for non-ARM platforms
///
/// This provides a portable implementation for x86_64 and other architectures.
pub fn gc_content_scalar(seq: &[u8]) -> GCContent {
    let mut gc_count = 0u32;
    let mut total_count = 0u32;

    for &base in seq {
        match base {
            b'G' | b'C' => {
                gc_count += 1;
                total_count += 1;
            }
            b'A' | b'T' => {
                total_count += 1;
            }
            _ => {} // Ignore non-ACGT characters
        }
    }

    if total_count == 0 {
        0.0
    } else {
        gc_count as f64 / total_count as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content_basic() {
        let seq = b"GATTACA"; // 2 GC out of 7 total
        let gc = gc_content(seq);
        assert!((gc - 0.285714).abs() < 0.001); // ~28.6%
    }

    #[test]
    fn test_gc_content_all_gc() {
        let seq = b"GCGCGCGC";
        let gc = gc_content(seq);
        assert!((gc - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_gc_content_all_at() {
        let seq = b"ATATATA";
        let gc = gc_content(seq);
        assert!((gc - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_gc_content_with_n() {
        let seq = b"GATNNCACN"; // 3 GC out of 6 ACGT bases (ignore N)
        let gc = gc_content(seq);
        assert!((gc - 0.5).abs() < 0.001); // 50%
    }

    #[test]
    fn test_gc_content_empty() {
        let seq = b"";
        let gc = gc_content(seq);
        assert_eq!(gc, 0.0);
    }

    #[test]
    fn test_gc_content_only_n() {
        let seq = b"NNNNN";
        let gc = gc_content(seq);
        assert_eq!(gc, 0.0);
    }

    #[test]
    fn test_gc_content_large() {
        // Test with sequence larger than 16 bytes (tests NEON chunking)
        let seq = b"GATTACAGATTACAGATTACAGATTACA"; // 8 GC out of 28 (4G + 4C)
        let gc = gc_content(seq);
        assert!((gc - 0.285714).abs() < 0.001); // ~28.6%
    }

    #[test]
    fn test_gc_content_50_percent() {
        let seq = b"GCGCAT"; // 4 GC out of 6 total
        let gc = gc_content(seq);
        assert!((gc - 0.666666).abs() < 0.001); // ~66.7%
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn test_neon_matches_scalar() {
        // Verify NEON and scalar implementations produce identical results
        let sequences = vec![
            b"GATTACA".as_slice(),
            b"AAAA".as_slice(),
            b"CCCC".as_slice(),
            b"GGGG".as_slice(),
            b"TTTT".as_slice(),
            b"ACGTACGTACGTACGT".as_slice(), // Exactly 16 bytes
            b"ACGTACGTACGTACGTACGT".as_slice(), // >16 bytes
            b"NNNACGTNNNN".as_slice(),       // With non-ACGT
        ];

        for seq in sequences {
            let neon_result = unsafe { gc_content_neon(seq) };
            let scalar_result = gc_content_scalar(seq);
            assert!(
                (neon_result - scalar_result).abs() < 0.0001,
                "NEON and scalar results differ for sequence: {:?} (neon: {}, scalar: {})",
                std::str::from_utf8(seq).unwrap(),
                neon_result,
                scalar_result
            );
        }
    }
}
