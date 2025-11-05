//! Trimming operations for FASTQ records
//!
//! Provides fixed-position and quality-based trimming operations commonly
//! used in read preprocessing pipelines.
//!
//! # Design Principles
//!
//! 1. **Preserve sequence-quality coupling**: Trim both simultaneously
//! 2. **Quality format**: Phred+33 (Illumina 1.8+, standard today)
//! 3. **Zero-copy when possible**: Fixed trimming uses slicing
//! 4. **Evidence-based**: Scalar-only (no NEON benefit for sequential ops)
//!
//! # Trimming Types
//!
//! **Fixed-Position**: Trim N bases from start/end
//! - Fast, deterministic
//! - Common for removing adapters or barcodes
//!
//! **Quality-Based**: Trim until quality threshold met
//! - Data-driven, adaptive
//! - Common for removing low-quality tails
//!
//! **Sliding Window**: Trim using moving window average
//! - Robust to noise
//! - More sophisticated than single-base threshold
//!
//! # Examples
//!
//! ```
//! use biometal::FastqRecord;
//! use biometal::operations::{trim_quality_end, trim_start};
//!
//! # fn main() -> biometal::Result<()> {
//! let record = FastqRecord::new(
//!     "read1".to_string(),
//!     b"ATGCATGC".to_vec(),
//!     b"IIII!!!!".to_vec(), // Q40, then Q0
//! );
//!
//! // Fixed trimming: Remove first 2 bases
//! let trimmed = trim_start(&record, 2)?;
//! assert_eq!(trimmed.sequence, b"GCATGC");
//!
//! // Quality trimming: Remove low-quality tail (Q<20)
//! let qual_trimmed = trim_quality_end(&record, 20)?;
//! assert_eq!(qual_trimmed.sequence, b"ATGC"); // Removed low-quality tail
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::types::FastqRecord;

/// Validate that a record's sequence and quality lengths match
///
/// # Returns
///
/// `Ok(())` if lengths match, otherwise `Err` with detailed message
fn validate_record_alignment(record: &FastqRecord) -> Result<()> {
    if record.sequence.len() != record.quality.len() {
        return Err(BiometalError::InvalidRange(format!(
            "Sequence length ({}) doesn't match quality length ({})",
            record.sequence.len(),
            record.quality.len()
        )));
    }
    Ok(())
}

/// Trim N bases from the 5' end (start of read)
///
/// Removes the first `bases` nucleotides from the sequence and corresponding
/// quality scores. If `bases` exceeds the sequence length, returns an empty record.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `bases` - Number of bases to trim from start
///
/// # Returns
///
/// New `FastqRecord` with bases trimmed from the start
///
/// # Performance
///
/// - Time: O(n) where n = remaining length (copies subsequence)
/// - Space: O(n) for new record allocation
/// - Evidence: Scalar-only (simple slicing operation)
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::trim_start;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"ABCDEFGH".to_vec(),
/// );
///
/// let trimmed = trim_start(&record, 2)?;
/// assert_eq!(trimmed.sequence, b"GCATGC");
/// assert_eq!(trimmed.quality, b"CDEFGH");
///
/// // Trimming entire sequence returns empty record
/// let all_trimmed = trim_start(&record, 100)?;
/// assert_eq!(all_trimmed.sequence, b"");
/// # Ok(())
/// # }
/// ```
pub fn trim_start(record: &FastqRecord, bases: usize) -> Result<FastqRecord> {
    validate_record_alignment(record)?;

    let start = bases.min(record.sequence.len());

    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[start..].to_vec(),
        record.quality[start..].to_vec(),
    ))
}

/// Trim N bases from the 3' end (end of read)
///
/// Removes the last `bases` nucleotides from the sequence and corresponding
/// quality scores. If `bases` exceeds the sequence length, returns an empty record.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `bases` - Number of bases to trim from end
///
/// # Returns
///
/// New `FastqRecord` with bases trimmed from the end
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::trim_end;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"ABCDEFGH".to_vec(),
/// );
///
/// let trimmed = trim_end(&record, 2)?;
/// assert_eq!(trimmed.sequence, b"ATGCAT");
/// assert_eq!(trimmed.quality, b"ABCDEF");
/// # Ok(())
/// # }
/// ```
pub fn trim_end(record: &FastqRecord, bases: usize) -> Result<FastqRecord> {
    validate_record_alignment(record)?;

    let end = record.sequence.len().saturating_sub(bases);

    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[..end].to_vec(),
        record.quality[..end].to_vec(),
    ))
}

/// Trim bases from both ends
///
/// Removes `start_bases` from the 5' end and `end_bases` from the 3' end.
/// If total trimming exceeds sequence length, returns an empty record.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `start_bases` - Number of bases to trim from start
/// * `end_bases` - Number of bases to trim from end
///
/// # Returns
///
/// New `FastqRecord` with bases trimmed from both ends
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::trim_both;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"ABCDEFGH".to_vec(),
/// );
///
/// let trimmed = trim_both(&record, 2, 2)?;
/// assert_eq!(trimmed.sequence, b"GCAT");
/// assert_eq!(trimmed.quality, b"CDEF");
/// # Ok(())
/// # }
/// ```
pub fn trim_both(record: &FastqRecord, start_bases: usize, end_bases: usize) -> Result<FastqRecord> {
    validate_record_alignment(record)?;

    let start = start_bases.min(record.sequence.len());
    let end = record.sequence.len().saturating_sub(end_bases);
    let final_end = end.max(start); // Ensure start <= end

    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[start..final_end].to_vec(),
        record.quality[start..final_end].to_vec(),
    ))
}

/// Trim from 3' end until minimum quality threshold is met
///
/// Scans from the 3' end backward, trimming bases until a base with quality
/// ≥ `min_quality` is found. This removes low-quality tails common in
/// Illumina sequencing.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `min_quality` - Minimum Phred quality score (Phred+33 encoding)
///
/// # Returns
///
/// New `FastqRecord` with low-quality tail removed
///
/// # Quality Score Format
///
/// - **Phred+33** (Illumina 1.8+, standard today)
/// - ASCII 33 = Q0, ASCII 73 = Q40
/// - Example: '!' = Q0, 'I' = Q40
///
/// # Algorithm
///
/// 1. Start at end of quality array
/// 2. Scan backward
/// 3. Stop when quality[i] >= min_quality
/// 4. Return record[0..i+1]
///
/// # Empty Records
///
/// If all bases are below the quality threshold, returns an `Ok(FastqRecord)`
/// with empty sequence and quality vectors. Callers can check using
/// `record.is_empty()` if empty records need special handling.
///
/// # Evidence
///
/// - Scalar-only: Sequential scanning doesn't benefit from NEON
/// - Standard bioinformatics preprocessing step
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::trim_quality_end;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"IIII!!!!".to_vec(), // Q40 then Q0
/// );
///
/// // Trim bases with Q<20 (ASCII 53)
/// let trimmed = trim_quality_end(&record, 20)?;
/// assert_eq!(trimmed.sequence, b"ATGC");
/// assert_eq!(trimmed.quality, b"IIII");
///
/// // No trimming if all bases are high quality
/// let good = FastqRecord::new(
///     "read2".to_string(),
///     b"ATGC".to_vec(),
///     b"IIII".to_vec(),
/// );
/// let not_trimmed = trim_quality_end(&good, 20)?;
/// assert_eq!(not_trimmed.sequence, b"ATGC");
/// # Ok(())
/// # }
/// ```
pub fn trim_quality_end(record: &FastqRecord, min_quality: u8) -> Result<FastqRecord> {
    validate_record_alignment(record)?;

    // Find last position with quality >= threshold
    let phred_threshold = min_quality + 33; // Convert to Phred+33 encoding

    let mut end = 0; // Start with 0 (will return empty if no good bases found)
    for (i, &qual) in record.quality.iter().enumerate().rev() {
        if qual >= phred_threshold {
            end = i + 1;
            break;
        }
    }

    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[..end].to_vec(),
        record.quality[..end].to_vec(),
    ))
}

/// Trim from 5' end until minimum quality threshold is met
///
/// Scans from the 5' end forward, trimming bases until a base with quality
/// ≥ `min_quality` is found. Less common than 3' trimming but useful for
/// some protocols.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `min_quality` - Minimum Phred quality score (Phred+33 encoding)
///
/// # Returns
///
/// New `FastqRecord` with low-quality prefix removed
///
/// # Empty Records
///
/// If all bases are below the quality threshold, returns an `Ok(FastqRecord)`
/// with empty sequence and quality vectors. Callers can check using
/// `record.is_empty()` if empty records need special handling.
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::trim_quality_start;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"!!!!IIII".to_vec(), // Q0 then Q40
/// );
///
/// let trimmed = trim_quality_start(&record, 20)?;
/// assert_eq!(trimmed.sequence, b"ATGC");
/// assert_eq!(trimmed.quality, b"IIII");
/// # Ok(())
/// # }
/// ```
pub fn trim_quality_start(record: &FastqRecord, min_quality: u8) -> Result<FastqRecord> {
    validate_record_alignment(record)?;

    // Find first position with quality >= threshold
    let phred_threshold = min_quality + 33;

    let mut start = record.quality.len();
    for (i, &qual) in record.quality.iter().enumerate() {
        if qual >= phred_threshold {
            start = i;
            break;
        }
    }

    // If all bases are below threshold, return empty record
    if start == record.quality.len() {
        return Ok(FastqRecord::new(
            record.id.clone(),
            Vec::new(),
            Vec::new(),
        ));
    }

    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[start..].to_vec(),
        record.quality[start..].to_vec(),
    ))
}

/// Trim from both ends until minimum quality threshold is met
///
/// Applies quality trimming from both 5' and 3' ends. This is the most
/// common trimming operation in preprocessing pipelines.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `min_quality` - Minimum Phred quality score (Phred+33 encoding)
///
/// # Returns
///
/// New `FastqRecord` with low-quality ends removed
///
/// # Empty Records
///
/// If all bases are below the quality threshold, returns an `Ok(FastqRecord)`
/// with empty sequence and quality vectors. Callers can check using
/// `record.is_empty()` if empty records need special handling.
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::trim_quality_both;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"!!IIII!!".to_vec(), // Low quality at both ends
/// );
///
/// let trimmed = trim_quality_both(&record, 20)?;
/// // Trim from start removes first 2, trim from end removes last 2
/// assert_eq!(trimmed.sequence, b"GCAT");
/// assert_eq!(trimmed.quality, b"IIII");
/// # Ok(())
/// # }
/// ```
pub fn trim_quality_both(record: &FastqRecord, min_quality: u8) -> Result<FastqRecord> {
    // Validate alignment
    if record.sequence.len() != record.quality.len() {
        return Err(BiometalError::InvalidRange(format!(
            "Sequence length ({}) doesn't match quality length ({})",
            record.sequence.len(),
            record.quality.len()
        )));
    }

    let phred_threshold = min_quality + 33;

    // Find first good position from start
    let mut start = record.quality.len();
    for (i, &qual) in record.quality.iter().enumerate() {
        if qual >= phred_threshold {
            start = i;
            break;
        }
    }

    // If all low quality, return empty
    if start == record.quality.len() {
        return Ok(FastqRecord::new(
            record.id.clone(),
            Vec::new(),
            Vec::new(),
        ));
    }

    // Find last good position from end (searching from start position onward)
    let mut end = start + 1; // If we find no good bases from end, keep at least the first good base
    for (i, &qual) in record.quality[start..].iter().enumerate().rev() {
        if qual >= phred_threshold {
            end = start + i + 1;
            break;
        }
    }

    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[start..end].to_vec(),
        record.quality[start..end].to_vec(),
    ))
}

/// Trim using sliding window quality average
///
/// Uses a sliding window approach where the mean quality of a window must
/// exceed the threshold. This is more robust to noise than single-base
/// thresholding.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `min_quality` - Minimum mean Phred quality for window
/// * `window_size` - Size of sliding window
///
/// # Returns
///
/// New `FastqRecord` trimmed where window quality drops below threshold
///
/// # Algorithm
///
/// 1. Slide window from 3' end backward
/// 2. Calculate mean quality for each window
/// 3. Trim at first window with mean < threshold
/// 4. Return trimmed record
///
/// # Empty Records
///
/// If no window meets the quality threshold, returns an `Ok(FastqRecord)`
/// with empty sequence and quality vectors. Callers can check using
/// `record.is_empty()` if empty records need special handling.
///
/// # Evidence
///
/// - Scalar-only: Sequential scanning with accumulation
/// - Standard Trimmomatic/cutadapt approach
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::trim_quality_window;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"IIII!!!!".to_vec(), // Good quality then bad
/// );
///
/// // Window size 4, threshold Q20
/// let trimmed = trim_quality_window(&record, 20, 4)?;
/// assert_eq!(trimmed.sequence, b"ATGCAT"); // Trimmed when window average dropped
/// # Ok(())
/// # }
/// ```
pub fn trim_quality_window(record: &FastqRecord, min_quality: u8, window_size: usize) -> Result<FastqRecord> {
    validate_record_alignment(record)?;

    if window_size == 0 {
        return Err(BiometalError::InvalidRange(
            "Window size must be > 0".to_string()
        ));
    }

    if record.quality.len() < window_size {
        // If record is shorter than window, check entire record
        let mean_qual = record.quality.iter()
            .map(|&q| q as u32 - 33)
            .sum::<u32>() / record.quality.len() as u32;

        if mean_qual < min_quality as u32 {
            return Ok(FastqRecord::new(
                record.id.clone(),
                Vec::new(),
                Vec::new(),
            ));
        } else {
            return Ok(record.clone());
        }
    }

    let phred_threshold = min_quality as u32;

    // Scan from end backward with sliding window
    let mut end = record.quality.len();
    for i in (window_size..=record.quality.len()).rev() {
        let window_start = i - window_size;
        let window = &record.quality[window_start..i];

        // Calculate mean quality for window (Phred scores, not ASCII)
        let mean_qual = window.iter()
            .map(|&q| (q as u32).saturating_sub(33))
            .sum::<u32>() / window_size as u32;

        if mean_qual >= phred_threshold {
            end = i;
            break;
        }
    }

    // If no good window found, return empty
    if end == 0 {
        return Ok(FastqRecord::new(
            record.id.clone(),
            Vec::new(),
            Vec::new(),
        ));
    }

    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[..end].to_vec(),
        record.quality[..end].to_vec(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Fixed-Position Trimming Tests =====

    #[test]
    fn test_trim_start_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"ABCDEFGH".to_vec(),
        );

        let trimmed = trim_start(&record, 2).unwrap();
        assert_eq!(trimmed.sequence, b"GCATGC");
        assert_eq!(trimmed.quality, b"CDEFGH");
        assert_eq!(trimmed.id, "read1");
    }

    #[test]
    fn test_trim_start_excessive() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        let trimmed = trim_start(&record, 100).unwrap();
        assert_eq!(trimmed.sequence, b"");
        assert_eq!(trimmed.quality, b"");
    }

    #[test]
    fn test_trim_end_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"ABCDEFGH".to_vec(),
        );

        let trimmed = trim_end(&record, 2).unwrap();
        assert_eq!(trimmed.sequence, b"ATGCAT");
        assert_eq!(trimmed.quality, b"ABCDEF");
    }

    #[test]
    fn test_trim_both_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"ABCDEFGH".to_vec(),
        );

        let trimmed = trim_both(&record, 2, 2).unwrap();
        assert_eq!(trimmed.sequence, b"GCAT");
        assert_eq!(trimmed.quality, b"CDEF");
    }

    #[test]
    fn test_trim_both_excessive() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        let trimmed = trim_both(&record, 2, 3).unwrap();
        assert_eq!(trimmed.sequence, b"");
        assert_eq!(trimmed.quality, b"");
    }

    // ===== Quality-Based Trimming Tests =====

    #[test]
    fn test_trim_quality_end_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"IIII!!!!".to_vec(), // Q40 then Q0
        );

        // Trim Q<20 (ASCII 53 = '5')
        let trimmed = trim_quality_end(&record, 20).unwrap();
        assert_eq!(trimmed.sequence, b"ATGC");
        assert_eq!(trimmed.quality, b"IIII");
    }

    #[test]
    fn test_trim_quality_end_no_trim_needed() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"IIII".to_vec(), // All Q40
        );

        let trimmed = trim_quality_end(&record, 20).unwrap();
        assert_eq!(trimmed.sequence, b"ATGC");
        assert_eq!(trimmed.quality, b"IIII");
    }

    #[test]
    fn test_trim_quality_end_all_low() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"!!!!".to_vec(), // All Q0
        );

        let trimmed = trim_quality_end(&record, 20).unwrap();
        assert_eq!(trimmed.sequence, b"");
        assert_eq!(trimmed.quality, b"");
    }

    #[test]
    fn test_trim_quality_start_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"!!!!IIII".to_vec(), // Q0 then Q40
        );

        let trimmed = trim_quality_start(&record, 20).unwrap();
        assert_eq!(trimmed.sequence, b"ATGC");
        assert_eq!(trimmed.quality, b"IIII");
    }

    #[test]
    fn test_trim_quality_both_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"!!IIII!!".to_vec(), // Low at both ends
        );

        let trimmed = trim_quality_both(&record, 20).unwrap();
        // Trim from start: remove first 2 (!!) → "GCATGC" with "IIII!!"
        // Trim from end: remove last 2 (!!) → "GCAT" with "IIII"
        assert_eq!(trimmed.sequence, b"GCAT");
        assert_eq!(trimmed.quality, b"IIII");
    }

    // ===== Window-Based Trimming Tests =====

    #[test]
    fn test_trim_quality_window_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"IIII!!!!".to_vec(), // Good then bad
        );

        // Window of 4 will catch the quality drop
        // Window: IIII (mean Q40) -> III! (mean Q30) -> II!! (mean Q20) -> I!!! (mean Q10)
        let trimmed = trim_quality_window(&record, 20, 4).unwrap();
        // Should trim when window mean drops below 20
        assert_eq!(trimmed.sequence.len(), 6); // Keeps first 6, trims last 2
        assert_eq!(trimmed.sequence, b"ATGCAT");
    }

    #[test]
    fn test_trim_quality_window_short_sequence() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"AT".to_vec(),
            b"II".to_vec(),
        );

        // Window larger than sequence
        let trimmed = trim_quality_window(&record, 20, 10).unwrap();
        assert_eq!(trimmed.sequence, b"AT"); // High quality, keep all
    }

    #[test]
    fn test_trim_quality_window_zero_window() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"IIII".to_vec(),
        );

        assert!(trim_quality_window(&record, 20, 0).is_err());
    }

    // ===== Edge Cases =====

    #[test]
    fn test_trim_empty_record() {
        let record = FastqRecord::new(
            "empty".to_string(),
            b"".to_vec(),
            b"".to_vec(),
        );

        assert!(trim_start(&record, 0).is_ok());
        assert!(trim_end(&record, 0).is_ok());
        assert!(trim_quality_end(&record, 20).is_ok());
    }

    #[test]
    fn test_trim_mismatched_lengths() {
        let mut record = FastqRecord::new(
            "bad".to_string(),
            b"ATGC".to_vec(),
            b"AB".to_vec(), // Too short
        );

        assert!(trim_start(&record, 1).is_err());
        assert!(trim_quality_end(&record, 20).is_err());

        // Fix and verify
        record.quality = b"ABCD".to_vec();
        assert!(trim_start(&record, 1).is_ok());
    }

    // ===== Property-Based Tests =====

    #[cfg(test)]
    mod properties {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            /// Property: Trimming preserves sequence-quality alignment
            #[test]
            fn prop_trim_preserves_alignment(
                seq in "[ACGTN]{10,100}",
                trim_amount in 0usize..20,
            ) {
                let seq_bytes = seq.as_bytes();
                let qual = vec![b'I'; seq_bytes.len()];
                let record = FastqRecord::new("test".to_string(), seq_bytes.to_vec(), qual);

                let trimmed = trim_start(&record, trim_amount).unwrap();
                prop_assert_eq!(trimmed.sequence.len(), trimmed.quality.len());

                let trimmed_end = trim_end(&record, trim_amount).unwrap();
                prop_assert_eq!(trimmed_end.sequence.len(), trimmed_end.quality.len());
            }

            /// Property: Trimming never increases length
            #[test]
            fn prop_trim_never_increases_length(
                seq in "[ACGTN]{10,100}",
                trim_amount in 0usize..20,
            ) {
                let seq_bytes = seq.as_bytes();
                let qual = vec![b'I'; seq_bytes.len()];
                let record = FastqRecord::new("test".to_string(), seq_bytes.to_vec(), qual);

                let trimmed = trim_start(&record, trim_amount).unwrap();
                prop_assert!(trimmed.sequence.len() <= record.sequence.len());

                let trimmed_end = trim_end(&record, trim_amount).unwrap();
                prop_assert!(trimmed_end.sequence.len() <= record.sequence.len());
            }

            /// Property: Quality trimming result has all bases >= threshold
            #[test]
            fn prop_quality_trim_meets_threshold(
                seq_len in 10usize..100,
                threshold in 10u8..30,
            ) {
                // Create sequence with varying quality
                let seq = vec![b'A'; seq_len];
                let mut qual = vec![b'I'; seq_len]; // Start with high quality

                // Add low-quality tail
                if seq_len > 10 {
                    for i in (seq_len - 5)..seq_len {
                        qual[i] = b'!'; // Q0
                    }
                }

                let record = FastqRecord::new("test".to_string(), seq, qual);
                let trimmed = trim_quality_end(&record, threshold).unwrap();

                // All remaining bases should have quality >= threshold
                let phred_threshold = threshold + 33;
                for &q in &trimmed.quality {
                    prop_assert!(q >= phred_threshold);
                }
            }

            /// Property: Window trimming is deterministic
            #[test]
            fn prop_window_trim_deterministic(
                seq_len in 20usize..100,
                window_size in 2usize..10,
            ) {
                let seq = vec![b'A'; seq_len];
                let qual = vec![b'I'; seq_len];
                let record = FastqRecord::new("test".to_string(), seq, qual);

                let result1 = trim_quality_window(&record, 20, window_size).unwrap();
                let result2 = trim_quality_window(&record, 20, window_size).unwrap();

                prop_assert_eq!(result1.sequence, result2.sequence);
                prop_assert_eq!(result1.quality, result2.quality);
            }

            /// Property: Window size 1 should equal single-base threshold
            #[test]
            fn prop_window_size_one_equals_single_threshold(
                seq_len in 20usize..100,
                threshold in 10u8..30,
            ) {
                let seq = vec![b'A'; seq_len];
                let mut qual = vec![b'I'; seq_len];

                // Add low-quality tail
                if seq_len > 10 {
                    for i in (seq_len - 5)..seq_len {
                        qual[i] = b'!'; // Q0
                    }
                }

                let record = FastqRecord::new("test".to_string(), seq, qual);

                // Window size 1 should behave identically to single-base threshold
                let window_result = trim_quality_window(&record, threshold, 1).unwrap();
                let single_result = trim_quality_end(&record, threshold).unwrap();

                prop_assert_eq!(window_result.sequence, single_result.sequence,
                    "Window size 1 should match single-base trimming");
            }
        }
    }
}
