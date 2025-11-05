//! Record-level sequence operations
//!
//! Operations that work on complete `FastqRecord` structures, preserving
//! the relationship between sequence, quality scores, and metadata.
//!
//! # Design Principles
//!
//! 1. **Preserve sequence-quality coupling**: Operations maintain alignment
//! 2. **Zero-copy when possible**: Use slices for subsequence extraction
//! 3. **Validation**: Check bounds and ensure sequence/quality length match
//! 4. **Error handling**: Return `Result` for operations that can fail
//!
//! # Examples
//!
//! ```
//! use biometal::FastqRecord;
//! use biometal::operations::{extract_region, reverse_complement_record};
//!
//! # fn main() -> biometal::Result<()> {
//! let record = FastqRecord::new(
//!     "read1".to_string(),
//!     b"ATGCATGC".to_vec(),
//!     b"IIIIIIII".to_vec(),
//! );
//!
//! // Extract region (maintains sequence-quality alignment)
//! let extracted = extract_region(&record, 2, 6)?;
//! assert_eq!(extracted.sequence, b"GCAT");
//! assert_eq!(extracted.quality, b"IIII");
//!
//! // Reverse complement (sequence + quality reversed)
//! let rc_record = reverse_complement_record(&record);
//! assert_eq!(rc_record.sequence, b"GCATGCAT");
//! assert_eq!(rc_record.quality, b"IIIIIIII"); // Also reversed
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::operations::sequence::reverse_complement;
use crate::types::FastqRecord;

/// Extract a region from a FASTQ record (sequence + quality)
///
/// Creates a new record containing the subsequence from `start` to `end`
/// (exclusive), along with the corresponding quality scores.
///
/// # Arguments
///
/// * `record` - Source FASTQ record
/// * `start` - Start position (0-indexed, inclusive)
/// * `end` - End position (0-indexed, exclusive)
///
/// # Returns
///
/// New `FastqRecord` containing the extracted region
///
/// # Errors
///
/// Returns `BiometalError::InvalidRange` if:
/// - `start >= end`
/// - `end > sequence length`
/// - `sequence and quality lengths don't match`
///
/// # Performance
///
/// - Time: O(n) where n = extracted length (copies subsequence)
/// - Space: O(n) for new record allocation
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::extract_region;
///
/// # fn main() -> biometal::Result<()> {
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGCATGC".to_vec(),
///     b"IIIIIIII".to_vec(),
/// );
///
/// // Extract middle 4 bases
/// let extracted = extract_region(&record, 2, 6)?;
/// assert_eq!(extracted.id, "read1");
/// assert_eq!(extracted.sequence, b"GCAT");
/// assert_eq!(extracted.quality, b"IIII");
///
/// // Extract from start
/// let prefix = extract_region(&record, 0, 4)?;
/// assert_eq!(prefix.sequence, b"ATGC");
/// assert_eq!(prefix.quality, b"IIII");
///
/// // Extract to end
/// let suffix = extract_region(&record, 4, 8)?;
/// assert_eq!(suffix.sequence, b"ATGC");
/// assert_eq!(suffix.quality, b"IIII");
/// # Ok(())
/// # }
/// ```
pub fn extract_region(record: &FastqRecord, start: usize, end: usize) -> Result<FastqRecord> {
    // Validate inputs
    if record.sequence.len() != record.quality.len() {
        return Err(BiometalError::InvalidRange(format!(
            "Sequence length ({}) doesn't match quality length ({})",
            record.sequence.len(),
            record.quality.len()
        )));
    }

    if start >= end {
        return Err(BiometalError::InvalidRange(format!(
            "start ({}) >= end ({})",
            start, end
        )));
    }

    if end > record.sequence.len() {
        return Err(BiometalError::InvalidRange(format!(
            "end ({}) exceeds sequence length ({})",
            end,
            record.sequence.len()
        )));
    }

    // Extract subsequences
    Ok(FastqRecord::new(
        record.id.clone(),
        record.sequence[start..end].to_vec(),
        record.quality[start..end].to_vec(),
    ))
}

/// Reverse complement a FASTQ record (sequence + quality)
///
/// Creates a new record with:
/// - Sequence: Reverse complemented
/// - Quality: Reversed to maintain alignment with sequence
/// - ID: Preserved
///
/// # Arguments
///
/// * `record` - Source FASTQ record
///
/// # Returns
///
/// New `FastqRecord` with reverse complemented sequence and reversed quality
///
/// # Performance
///
/// - Time: O(n) where n = sequence length
/// - Space: O(n) for new record allocation
/// - Future: NEON optimization will provide 10-16× speedup
///
/// # Evidence
///
/// - Category 2 (Similar): Inherits from `reverse_complement()` optimization
/// - Expected NEON speedup: 10-16× (from Entry 020-025 patterns)
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::reverse_complement_record;
///
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGC".to_vec(),
///     b"ABCD".to_vec(),
/// );
///
/// let rc = reverse_complement_record(&record);
/// assert_eq!(rc.id, "read1");
/// assert_eq!(rc.sequence, b"GCAT");
/// assert_eq!(rc.quality, b"DCBA"); // Quality reversed to maintain alignment
///
/// // Involutive: RC(RC(x)) = x
/// let rc_rc = reverse_complement_record(&rc);
/// assert_eq!(rc_rc.sequence, record.sequence);
/// assert_eq!(rc_rc.quality, record.quality);
/// ```
pub fn reverse_complement_record(record: &FastqRecord) -> FastqRecord {
    let rc_seq = reverse_complement(&record.sequence);
    let rev_qual = record.quality.iter().rev().copied().collect();

    FastqRecord::new(record.id.clone(), rc_seq, rev_qual)
}

/// In-place reverse complement a FASTQ record
///
/// Mutates the record in place, avoiding allocation. Useful for large
/// records where allocation cost matters.
///
/// # Arguments
///
/// * `record` - Mutable FASTQ record (modified in place)
///
/// # Performance
///
/// - Time: O(n) where n = sequence length
/// - Space: O(1) - zero allocation (in-place transformation)
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::reverse_complement_record_inplace;
///
/// let mut record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGC".to_vec(),
///     b"ABCD".to_vec(),
/// );
///
/// reverse_complement_record_inplace(&mut record);
/// assert_eq!(record.sequence, b"GCAT");
/// assert_eq!(record.quality, b"DCBA");
///
/// // Involutive: applying twice returns to original
/// reverse_complement_record_inplace(&mut record);
/// assert_eq!(record.sequence, b"ATGC");
/// assert_eq!(record.quality, b"ABCD");
/// ```
pub fn reverse_complement_record_inplace(record: &mut FastqRecord) {
    // Reverse complement sequence in place
    crate::operations::sequence::reverse_complement_inplace(&mut record.sequence);

    // Reverse quality scores
    record.quality.reverse();
}

/// Get sequence length from a record
///
/// Convenience wrapper for checking record length.
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::sequence_length;
///
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGC".to_vec(),
///     b"IIII".to_vec(),
/// );
///
/// assert_eq!(sequence_length(&record), 4);
/// ```
pub fn sequence_length(record: &FastqRecord) -> usize {
    record.sequence.len()
}

/// Check if record meets length requirements
///
/// Validates that a record's sequence length falls within specified bounds.
///
/// # Arguments
///
/// * `record` - FASTQ record to check
/// * `min` - Minimum acceptable length (inclusive)
/// * `max` - Maximum acceptable length (inclusive)
///
/// # Returns
///
/// `true` if `min <= length <= max`, `false` otherwise
///
/// # Examples
///
/// ```
/// use biometal::FastqRecord;
/// use biometal::operations::meets_length_requirement;
///
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGC".to_vec(),
///     b"IIII".to_vec(),
/// );
///
/// assert!(meets_length_requirement(&record, 3, 10));  // 4 is in [3, 10]
/// assert!(!meets_length_requirement(&record, 5, 10)); // 4 < 5
/// assert!(!meets_length_requirement(&record, 1, 3));  // 4 > 3
/// ```
pub fn meets_length_requirement(record: &FastqRecord, min: usize, max: usize) -> bool {
    let len = record.sequence.len();
    len >= min && len <= max
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Extract Region Tests =====

    #[test]
    fn test_extract_region_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"ABCDEFGH".to_vec(),
        );

        let extracted = extract_region(&record, 2, 6).unwrap();
        assert_eq!(extracted.id, "read1");
        assert_eq!(extracted.sequence, b"GCAT");
        assert_eq!(extracted.quality, b"CDEF");
    }

    #[test]
    fn test_extract_region_full() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        let extracted = extract_region(&record, 0, 4).unwrap();
        assert_eq!(extracted.sequence, record.sequence);
        assert_eq!(extracted.quality, record.quality);
    }

    #[test]
    fn test_extract_region_start() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"ABCDEFGH".to_vec(),
        );

        let extracted = extract_region(&record, 0, 4).unwrap();
        assert_eq!(extracted.sequence, b"ATGC");
        assert_eq!(extracted.quality, b"ABCD");
    }

    #[test]
    fn test_extract_region_end() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"ABCDEFGH".to_vec(),
        );

        let extracted = extract_region(&record, 4, 8).unwrap();
        assert_eq!(extracted.sequence, b"ATGC");
        assert_eq!(extracted.quality, b"EFGH");
    }

    #[test]
    fn test_extract_region_single_base() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        let extracted = extract_region(&record, 1, 2).unwrap();
        assert_eq!(extracted.sequence, b"T");
        assert_eq!(extracted.quality, b"B");
    }

    #[test]
    fn test_extract_region_invalid_range() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        // start >= end
        assert!(extract_region(&record, 2, 2).is_err());
        assert!(extract_region(&record, 3, 2).is_err());

        // end > length
        assert!(extract_region(&record, 0, 5).is_err());
        assert!(extract_region(&record, 2, 10).is_err());
    }

    #[test]
    fn test_extract_region_mismatched_lengths() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"AB".to_vec(), // Too short
        );

        assert!(extract_region(&record, 0, 2).is_err());

        // Fix and test
        record.quality = b"ABCD".to_vec();
        assert!(extract_region(&record, 0, 2).is_ok());
    }

    // ===== Reverse Complement Record Tests =====

    #[test]
    fn test_reverse_complement_record_basic() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        let rc = reverse_complement_record(&record);
        assert_eq!(rc.id, "read1");
        assert_eq!(rc.sequence, b"GCAT");
        assert_eq!(rc.quality, b"DCBA"); // Quality reversed
    }

    #[test]
    fn test_reverse_complement_record_involutive() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCATGC".to_vec(),
            b"ABCDEFGH".to_vec(),
        );

        let rc = reverse_complement_record(&record);
        let rc_rc = reverse_complement_record(&rc);

        assert_eq!(rc_rc.sequence, record.sequence);
        assert_eq!(rc_rc.quality, record.quality);
        assert_eq!(rc_rc.id, record.id);
    }

    #[test]
    fn test_reverse_complement_record_with_n() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGCN".to_vec(),
            b"ABCDE".to_vec(),
        );

        let rc = reverse_complement_record(&record);
        assert_eq!(rc.sequence, b"NGCAT");
        assert_eq!(rc.quality, b"EDCBA");
    }

    #[test]
    fn test_reverse_complement_record_inplace_basic() {
        let mut record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        reverse_complement_record_inplace(&mut record);
        assert_eq!(record.sequence, b"GCAT");
        assert_eq!(record.quality, b"DCBA");
    }

    #[test]
    fn test_reverse_complement_record_inplace_involutive() {
        let original = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        let mut record = original.clone();
        reverse_complement_record_inplace(&mut record);
        reverse_complement_record_inplace(&mut record);

        assert_eq!(record.sequence, original.sequence);
        assert_eq!(record.quality, original.quality);
    }

    // ===== Length Operations Tests =====

    #[test]
    fn test_sequence_length() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        assert_eq!(sequence_length(&record), 4);
    }

    #[test]
    fn test_meets_length_requirement() {
        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"ABCD".to_vec(),
        );

        // Within range
        assert!(meets_length_requirement(&record, 3, 10));
        assert!(meets_length_requirement(&record, 4, 4)); // Exact match
        assert!(meets_length_requirement(&record, 1, 100));

        // Outside range
        assert!(!meets_length_requirement(&record, 5, 10)); // Too short
        assert!(!meets_length_requirement(&record, 1, 3));  // Too long
    }

    #[test]
    fn test_meets_length_requirement_edge_cases() {
        let empty = FastqRecord::new(
            "empty".to_string(),
            b"".to_vec(),
            b"".to_vec(),
        );

        assert!(meets_length_requirement(&empty, 0, 0));
        assert!(meets_length_requirement(&empty, 0, 10));
        assert!(!meets_length_requirement(&empty, 1, 10));
    }

    // ===== Property-Based Tests =====

    #[cfg(test)]
    mod properties {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            /// Property: extract_region preserves sequence-quality alignment
            #[test]
            fn prop_extract_preserves_alignment(
                seq in "[ACGTN]{10,100}",
                start in 0usize..50,
            ) {
                let seq_bytes = seq.as_bytes();
                let qual = vec![b'I'; seq_bytes.len()];
                let record = FastqRecord::new("test".to_string(), seq_bytes.to_vec(), qual);

                if start < record.sequence.len() {
                    let end = (start + 1).min(record.sequence.len());
                    let extracted = extract_region(&record, start, end).unwrap();
                    prop_assert_eq!(extracted.sequence.len(), extracted.quality.len());
                }
            }

            /// Property: RC(RC(record)) = record
            #[test]
            fn prop_reverse_complement_record_involutive(seq in "[ACGTN]{1,1000}") {
                let seq_bytes = seq.as_bytes();
                let qual = vec![b'I'; seq_bytes.len()];
                let record = FastqRecord::new("test".to_string(), seq_bytes.to_vec(), qual.clone());

                let rc = reverse_complement_record(&record);
                let rc_rc = reverse_complement_record(&rc);

                prop_assert_eq!(rc_rc.sequence, record.sequence);
                prop_assert_eq!(rc_rc.quality, record.quality);
            }

            /// Property: In-place matches allocating
            #[test]
            fn prop_inplace_matches_allocating(seq in "[ACGTN]{1,1000}") {
                let seq_bytes = seq.as_bytes();
                let qual = vec![b'I'; seq_bytes.len()];

                let record = FastqRecord::new("test".to_string(), seq_bytes.to_vec(), qual.clone());
                let rc_alloc = reverse_complement_record(&record);

                let mut record_inplace = record.clone();
                reverse_complement_record_inplace(&mut record_inplace);

                prop_assert_eq!(rc_alloc.sequence, record_inplace.sequence);
                prop_assert_eq!(rc_alloc.quality, record_inplace.quality);
            }

            /// Property: Extracted region is subset of original
            #[test]
            fn prop_extract_is_subset(
                seq in "[ACGTN]{10,100}",
                start in 0usize..50,
                len in 1usize..20,
            ) {
                let seq_bytes = seq.as_bytes();
                let qual = vec![b'I'; seq_bytes.len()];
                let record = FastqRecord::new("test".to_string(), seq_bytes.to_vec(), qual);

                if start < record.sequence.len() {
                    let end = (start + len).min(record.sequence.len());
                    if start < end {
                        let extracted = extract_region(&record, start, end).unwrap();
                        prop_assert_eq!(extracted.sequence, &record.sequence[start..end]);
                        prop_assert_eq!(extracted.quality, &record.quality[start..end]);
                    }
                }
            }
        }
    }
}
