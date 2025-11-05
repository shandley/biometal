//! Masking operations for low-quality bases
//!
//! Provides operations to mask (replace with 'N') bases that fall below
//! quality thresholds. Unlike trimming, masking preserves read length,
//! which is important for some downstream analysis tools.
//!
//! # Design Principles
//!
//! 1. **Preserve read length**: Replace bases with 'N' instead of removing
//! 2. **Quality format**: Phred+33 (Illumina 1.8+, standard today)
//! 3. **In-place operation**: Mutates record for efficiency
//! 4. **Evidence-based**: Scalar-only (sequential operation)
//!
//! # Use Cases
//!
//! - **Alignment**: Some aligners handle 'N' specially
//! - **Variant calling**: Mark unreliable bases without changing coordinates
//! - **Length-dependent tools**: Tools that require fixed read length
//!
//! # Examples
//!
//! ```
//! use biometal::FastqRecord;
//! use biometal::operations::mask_low_quality;
//!
//! # fn main() -> biometal::Result<()> {
//! let mut record = FastqRecord::new(
//!     "read1".to_string(),
//!     b"ATGCATGC".to_vec(),
//!     b"IIII!!!!".to_vec(), // Q40 then Q0
//! );
//!
//! // Mask bases with Q<20 (replace with 'N')
//! mask_low_quality(&mut record, 20)?;
//! assert_eq!(record.sequence, b"ATGCNNNN"); // Low-quality bases masked
//! assert_eq!(record.quality, b"IIII!!!!"); // Quality unchanged
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::types::FastqRecord;

/// Mask (replace with 'N') bases below quality threshold
///
/// Replaces bases with quality scores below `min_quality` with 'N'. This
/// preserves read length while marking unreliable bases. Quality scores
/// are left unchanged.
///
/// # Arguments
///
/// * `record` - Mutable FASTQ record (modified in place)
/// * `min_quality` - Minimum Phred quality score (Phred+33 encoding)
///
/// # Returns
///
/// `Ok(())` on success, `Err` if sequence/quality lengths don't match
///
/// # Quality Score Format
///
/// - **Phred+33** (Illumina 1.8+, standard today)
/// - ASCII 33 = Q0, ASCII 73 = Q40
/// - Example: '!' = Q0, 'I' = Q40
///
/// # Performance
///
/// - Time: O(n) single pass through sequence
/// - Space: O(1) - in-place modification
/// - Evidence: Scalar-only (sequential operation, no NEON benefit)
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::mask_low_quality;
///
/// # fn main() -> biometal::Result<()> {
/// let mut record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"IIII!!!!".to_vec(),
/// );
///
/// mask_low_quality(&mut record, 20)?;
/// assert_eq!(record.sequence, b"ATGCNNNN");
///
/// // Quality scores are preserved
/// assert_eq!(record.quality, b"IIII!!!!");
///
/// // Read length preserved
/// assert_eq!(record.sequence.len(), 8);
/// # Ok(())
/// # }
/// ```
///
/// # Use Cases
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::mask_low_quality;
///
/// # fn main() -> biometal::Result<()> {
/// // Variant calling: Mark unreliable bases without changing positions
/// let mut read = FastqRecord::new(
///     "chr1:12345".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"II!III!I".to_vec(), // Mixed quality
/// );
///
/// mask_low_quality(&mut read, 20)?;
/// // Positions 2 and 6 masked (Q0) but coordinates preserved
/// assert_eq!(read.sequence, b"ATNCATNC");
///
/// // Alignment: Some aligners treat 'N' specially
/// // BWA-MEM, bowtie2 can soft-clip or report ambiguous bases
/// # Ok(())
/// # }
/// ```
pub fn mask_low_quality(record: &mut FastqRecord, min_quality: u8) -> Result<()> {
    if record.sequence.len() != record.quality.len() {
        return Err(BiometalError::InvalidRange(format!(
            "Sequence length ({}) doesn't match quality length ({})",
            record.sequence.len(),
            record.quality.len()
        )));
    }

    let phred_threshold = min_quality + 33; // Convert to Phred+33 encoding

    for (base, &qual) in record.sequence.iter_mut().zip(&record.quality) {
        if qual < phred_threshold {
            *base = b'N';
        }
    }

    Ok(())
}

/// Mask low-quality bases (non-mutating)
///
/// Returns a new record with low-quality bases masked. Use this when you
/// need to preserve the original record.
///
/// # Arguments
///
/// * `record` - Source FASTQ record (unchanged)
/// * `min_quality` - Minimum Phred quality score
///
/// # Returns
///
/// New `FastqRecord` with low-quality bases masked
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::mask_low_quality_copy;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGC".to_vec(),
///     b"I!!!".to_vec(),
/// );
///
/// let masked = mask_low_quality_copy(&record, 20)?;
/// assert_eq!(masked.sequence, b"ANNN");
///
/// // Original unchanged
/// assert_eq!(record.sequence, b"ATGC");
/// # Ok(())
/// # }
/// ```
pub fn mask_low_quality_copy(record: &FastqRecord, min_quality: u8) -> Result<FastqRecord> {
    let mut masked = record.clone();
    mask_low_quality(&mut masked, min_quality)?;
    Ok(masked)
}

/// Count masked bases in a sequence
///
/// Counts the number of 'N' bases in a sequence. Useful for QC metrics
/// after masking operations.
///
/// # Arguments
///
/// * `record` - FASTQ record to analyze
///
/// # Returns
///
/// Number of 'N' bases in the sequence
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::{mask_low_quality, count_masked_bases};
///
/// # fn main() -> biometal::Result<()> {
/// let mut record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"IIII!!!!".to_vec(),
/// );
///
/// assert_eq!(count_masked_bases(&record), 0);
///
/// mask_low_quality(&mut record, 20)?;
/// assert_eq!(count_masked_bases(&record), 4); // 4 bases masked
/// # Ok(())
/// # }
/// ```
pub fn count_masked_bases(record: &FastqRecord) -> usize {
    record.sequence.iter().filter(|&&b| b == b'N' || b == b'n').count()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Basic Masking Tests =====

    #[test]
    fn test_mask_low_quality_basic() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"IIII!!!!".to_vec(), // Q40 then Q0
        );

        mask_low_quality(&mut record, 20).unwrap();
        assert_eq!(record.sequence, b"ATGCNNNN");
        assert_eq!(record.quality, b"IIII!!!!"); // Quality unchanged
        assert_eq!(record.id, "read1");
    }

    #[test]
    fn test_mask_low_quality_no_masking_needed() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"IIII".to_vec(), // All Q40
        );

        mask_low_quality(&mut record, 20).unwrap();
        assert_eq!(record.sequence, b"ATGC"); // No masking
        assert_eq!(record.quality, b"IIII");
    }

    #[test]
    fn test_mask_low_quality_all_masked() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"!!!!".to_vec(), // All Q0
        );

        mask_low_quality(&mut record, 20).unwrap();
        assert_eq!(record.sequence, b"NNNN"); // All masked
        assert_eq!(record.quality, b"!!!!"); // Quality unchanged
    }

    #[test]
    fn test_mask_low_quality_mixed() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"I!I!I!I!".to_vec(), // Alternating Q40/Q0
        );

        mask_low_quality(&mut record, 20).unwrap();
        assert_eq!(record.sequence, b"ANGNANGN"); // Every other base masked (8 bases)
        assert_eq!(record.quality, b"I!I!I!I!");
    }

    #[test]
    fn test_mask_low_quality_threshold() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"5555".to_vec(), // All Q20 (ASCII 53 = 20 + 33)
        );

        // Threshold Q20: should NOT mask Q20 bases
        mask_low_quality(&mut record, 20).unwrap();
        assert_eq!(record.sequence, b"ATGC");

        // Threshold Q21: should mask Q20 bases
        let mut record2 = record.clone();
        mask_low_quality(&mut record2, 21).unwrap();
        assert_eq!(record2.sequence, b"NNNN");
    }

    #[test]
    fn test_mask_low_quality_preserves_length() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"I!I!I!I!".to_vec(),
        );

        let original_len = record.sequence.len();
        mask_low_quality(&mut record, 20).unwrap();
        assert_eq!(record.sequence.len(), original_len); // Length preserved
    }

    #[test]
    fn test_mask_low_quality_mismatched_lengths() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"AB".to_vec(), // Too short
        );

        assert!(mask_low_quality(&mut record, 20).is_err());
    }

    // ===== Copy Variant Tests =====

    #[test]
    fn test_mask_low_quality_copy_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"I!!!".to_vec(),
        );

        let masked = mask_low_quality_copy(&record, 20).unwrap();
        assert_eq!(masked.sequence, b"ANNN");

        // Original unchanged
        assert_eq!(record.sequence, b"ATGC");
    }

    #[test]
    fn test_mask_low_quality_copy_preserves_original() {
        let original = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"IIII!!!!".to_vec(),
        );

        let original_seq = original.sequence.clone();
        let _masked = mask_low_quality_copy(&original, 20).unwrap();

        // Verify original unchanged
        assert_eq!(original.sequence, original_seq);
    }

    // ===== Count Masked Bases Tests =====

    #[test]
    fn test_count_masked_bases_none() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"IIII".to_vec(),
        );

        assert_eq!(count_masked_bases(&record), 0);
    }

    #[test]
    fn test_count_masked_bases_some() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCNNNN".to_vec(),
            b"IIII!!!!".to_vec(),
        );

        assert_eq!(count_masked_bases(&record), 4);
    }

    #[test]
    fn test_count_masked_bases_all() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"NNNN".to_vec(),
            b"!!!!".to_vec(),
        );

        assert_eq!(count_masked_bases(&record), 4);
    }

    #[test]
    fn test_count_masked_bases_lowercase() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCnnnn".to_vec(),
            b"IIII!!!!".to_vec(),
        );

        assert_eq!(count_masked_bases(&record), 4); // Counts lowercase 'n' too
    }

    #[test]
    fn test_count_masked_bases_mixed() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATNGCNAN".to_vec(),
            b"IIIIIIII".to_vec(),
        );

        assert_eq!(count_masked_bases(&record), 3);
    }

    // ===== Edge Cases =====

    #[test]
    fn test_mask_empty_record() {
        let mut record = FastqRecord::new(
            "empty".to_string(),
            b"".to_vec(),
            b"".to_vec(),
        );

        assert!(mask_low_quality(&mut record, 20).is_ok());
        assert_eq!(count_masked_bases(&record), 0);
    }

    #[test]
    fn test_mask_preserves_existing_n() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATNCATGC".to_vec(),
            b"IIIIIIII".to_vec(),
        );

        mask_low_quality(&mut record, 20).unwrap();
        assert_eq!(record.sequence, b"ATNCATGC"); // N preserved
        assert_eq!(count_masked_bases(&record), 1);
    }

    // ===== Property-Based Tests =====

    #[cfg(test)]
    mod properties {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            /// Property: Masking preserves sequence length
            #[test]
            fn prop_mask_preserves_length(
                seq_len in 10usize..100,
                threshold in 10u8..30,
            ) {
                let seq = vec![b'A'; seq_len];
                let qual = vec![b'I'; seq_len];
                let mut record = FastqRecord::new("test".to_string(), seq, qual);

                let original_len = record.sequence.len();
                mask_low_quality(&mut record, threshold).unwrap();

                prop_assert_eq!(record.sequence.len(), original_len);
            }

            /// Property: Masking preserves quality scores
            #[test]
            fn prop_mask_preserves_quality(
                seq_len in 10usize..100,
                threshold in 10u8..30,
            ) {
                let seq = vec![b'A'; seq_len];
                let qual = vec![b'I'; seq_len];
                let original_qual = qual.clone();

                let mut record = FastqRecord::new("test".to_string(), seq, qual);
                mask_low_quality(&mut record, threshold).unwrap();

                prop_assert_eq!(record.quality, original_qual);
            }

            /// Property: All masked bases have low quality
            #[test]
            fn prop_masked_bases_are_low_quality(
                seq_len in 10usize..100,
                threshold in 10u8..30,
            ) {
                let seq = vec![b'A'; seq_len];
                let mut qual = vec![b'I'; seq_len]; // Start with high quality

                // Add some low-quality bases
                if seq_len > 10 {
                    for i in 0..5 {
                        qual[i] = b'!'; // Q0
                    }
                }

                let mut record = FastqRecord::new("test".to_string(), seq, qual.clone());
                mask_low_quality(&mut record, threshold).unwrap();

                let phred_threshold = threshold + 33;
                for (i, &base) in record.sequence.iter().enumerate() {
                    if base == b'N' {
                        prop_assert!(qual[i] < phred_threshold);
                    }
                }
            }

            /// Property: High-quality bases not masked
            #[test]
            fn prop_high_quality_not_masked(
                seq_len in 10usize..100,
                threshold in 10u8..30,
            ) {
                let mut seq = vec![b'A'; seq_len];
                let qual = vec![b'I'; seq_len]; // All high quality (Q40)

                // Add some variety to sequence
                for (i, base) in seq.iter_mut().enumerate() {
                    *base = [b'A', b'C', b'G', b'T'][i % 4];
                }

                let original_seq = seq.clone();
                let mut record = FastqRecord::new("test".to_string(), seq, qual);

                mask_low_quality(&mut record, threshold).unwrap();

                // No bases should be masked (all high quality)
                prop_assert_eq!(record.sequence, original_seq);
            }

            /// Property: Copy variant matches in-place
            #[test]
            fn prop_copy_matches_inplace(
                seq_len in 10usize..100,
                threshold in 10u8..30,
            ) {
                let seq = vec![b'A'; seq_len];
                let qual = vec![b'I'; seq_len];

                let record = FastqRecord::new("test".to_string(), seq.clone(), qual.clone());
                let mut record_mut = FastqRecord::new("test".to_string(), seq, qual);

                let copied = mask_low_quality_copy(&record, threshold).unwrap();
                mask_low_quality(&mut record_mut, threshold).unwrap();

                prop_assert_eq!(copied.sequence, record_mut.sequence);
            }

            /// Property: Count matches actual N count
            #[test]
            fn prop_count_matches_actual(
                seq_len in 10usize..100,
            ) {
                let mut seq = vec![b'A'; seq_len];

                // Add some N bases
                if seq_len > 10 {
                    for i in 0..5 {
                        seq[i] = b'N';
                    }
                }

                let qual = vec![b'I'; seq_len];
                let record = FastqRecord::new("test".to_string(), seq.clone(), qual);

                let count = count_masked_bases(&record);
                let expected = seq.iter().filter(|&&b| b == b'N' || b == b'n').count();

                prop_assert_eq!(count, expected);
            }
        }
    }
}
