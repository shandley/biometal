//! Genomic coordinate types and operations.
//!
//! This module provides types for working with genomic coordinates:
//! - [`GenomicInterval`]: A genomic interval (chromosome, start, end)
//! - [`Strand`]: DNA strand orientation (+, -, or unknown)
//!
//! # Coordinate System
//!
//! All genomic intervals use **0-based, half-open** coordinates `[start, end)`:
//! - Start position is inclusive (0-based)
//! - End position is exclusive
//! - Length = end - start
//!
//! This matches the BED format and is standard in bioinformatics.
//!
//! # Examples
//!
//! ```
//! use biometal::formats::primitives::GenomicInterval;
//!
//! // Create interval: chr1:100-200 (0-based, half-open)
//! let interval = GenomicInterval::new("chr1".to_string(), 100, 200)?;
//!
//! assert_eq!(interval.chrom, "chr1");
//! assert_eq!(interval.start, 100);
//! assert_eq!(interval.end, 200);
//! assert_eq!(interval.length(), 100);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ```
//! use biometal::formats::primitives::GenomicInterval;
//!
//! // Check overlap
//! let int1 = GenomicInterval::new("chr1".to_string(), 100, 200)?;
//! let int2 = GenomicInterval::new("chr1".to_string(), 150, 250)?;
//!
//! assert!(int1.overlaps(&int2));
//! assert!(int2.overlaps(&int1));  // Symmetric
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! ```
//! use biometal::formats::primitives::Strand;
//! use std::str::FromStr;
//!
//! // Parse strand
//! assert_eq!(Strand::from_str("+")?, Strand::Forward);
//! assert_eq!(Strand::from_str("-")?, Strand::Reverse);
//! assert_eq!(Strand::from_str(".")?, Strand::Unknown);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::formats::primitives::{FormatError, Result};
use std::fmt;
use std::str::FromStr;

/// A genomic interval with chromosome and coordinates.
///
/// Coordinates are **0-based, half-open** `[start, end)`:
/// - `start`: Inclusive start position (0-based)
/// - `end`: Exclusive end position
/// - Length: `end - start`
///
/// # Invariants
///
/// - `start < end` (enforced by constructor)
/// - Chromosome name is non-empty
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::GenomicInterval;
///
/// // Valid interval
/// let interval = GenomicInterval::new("chr1".to_string(), 100, 200)?;
/// assert_eq!(interval.length(), 100);
///
/// // Invalid interval (start >= end)
/// let result = GenomicInterval::new("chr1".to_string(), 200, 100);
/// assert!(result.is_err());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GenomicInterval {
    /// Chromosome or contig name.
    pub chrom: String,

    /// Start position (0-based, inclusive).
    pub start: u64,

    /// End position (0-based, exclusive).
    pub end: u64,
}

impl GenomicInterval {
    /// Creates a new genomic interval.
    ///
    /// # Arguments
    ///
    /// * `chrom` - Chromosome or contig name (must be non-empty)
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    ///
    /// # Errors
    ///
    /// Returns [`FormatError::InvalidInterval`] if `start >= end`.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::GenomicInterval;
    ///
    /// let interval = GenomicInterval::new("chr1".to_string(), 100, 200)?;
    /// assert_eq!(interval.length(), 100);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn new(chrom: String, start: u64, end: u64) -> Result<Self> {
        if start >= end {
            return Err(FormatError::InvalidInterval { start, end });
        }

        Ok(GenomicInterval { chrom, start, end })
    }

    /// Returns the length of this interval.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::GenomicInterval;
    ///
    /// let interval = GenomicInterval::new("chr1".to_string(), 100, 250)?;
    /// assert_eq!(interval.length(), 150);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[inline]
    pub fn length(&self) -> u64 {
        self.end - self.start
    }

    /// Checks if this interval overlaps with another interval.
    ///
    /// Two intervals overlap if they share at least one base position.
    /// Different chromosomes never overlap.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::GenomicInterval;
    ///
    /// let int1 = GenomicInterval::new("chr1".to_string(), 100, 200)?;
    /// let int2 = GenomicInterval::new("chr1".to_string(), 150, 250)?;
    /// let int3 = GenomicInterval::new("chr1".to_string(), 300, 400)?;
    ///
    /// assert!(int1.overlaps(&int2));  // Overlap
    /// assert!(!int1.overlaps(&int3)); // No overlap
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn overlaps(&self, other: &Self) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    /// Checks if this interval completely contains another interval.
    ///
    /// Interval A contains interval B if:
    /// - Same chromosome
    /// - A.start <= B.start
    /// - A.end >= B.end
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::GenomicInterval;
    ///
    /// let outer = GenomicInterval::new("chr1".to_string(), 100, 300)?;
    /// let inner = GenomicInterval::new("chr1".to_string(), 150, 200)?;
    ///
    /// assert!(outer.contains(&inner));
    /// assert!(!inner.contains(&outer));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn contains(&self, other: &Self) -> bool {
        self.chrom == other.chrom && self.start <= other.start && self.end >= other.end
    }
}

impl fmt::Display for GenomicInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/// DNA strand orientation.
///
/// Represents the orientation of a genomic feature:
/// - `Forward`: Plus strand (+)
/// - `Reverse`: Minus strand (-)
/// - `Unknown`: Strand not specified or unknown (.)
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::Strand;
/// use std::str::FromStr;
///
/// assert_eq!(Strand::from_str("+")?, Strand::Forward);
/// assert_eq!(Strand::from_str("-")?, Strand::Reverse);
/// assert_eq!(Strand::from_str(".")?, Strand::Unknown);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    /// Plus strand (+)
    Forward,

    /// Minus strand (-)
    Reverse,

    /// Unknown or unspecified strand (.)
    Unknown,
}

impl FromStr for Strand {
    type Err = FormatError;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "." => Ok(Strand::Unknown),
            _ => Err(FormatError::InvalidStrand(s.to_string())),
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genomic_interval_new() {
        let interval = GenomicInterval::new("chr1".to_string(), 100, 200).unwrap();
        assert_eq!(interval.chrom, "chr1");
        assert_eq!(interval.start, 100);
        assert_eq!(interval.end, 200);
    }

    #[test]
    fn test_genomic_interval_invalid() {
        // start == end
        let result = GenomicInterval::new("chr1".to_string(), 100, 100);
        assert!(result.is_err());

        // start > end
        let result = GenomicInterval::new("chr1".to_string(), 200, 100);
        assert!(result.is_err());
    }

    #[test]
    fn test_genomic_interval_length() {
        let interval = GenomicInterval::new("chr1".to_string(), 100, 200).unwrap();
        assert_eq!(interval.length(), 100);

        let interval = GenomicInterval::new("chr1".to_string(), 0, 1).unwrap();
        assert_eq!(interval.length(), 1);
    }

    #[test]
    fn test_genomic_interval_overlaps() {
        let int1 = GenomicInterval::new("chr1".to_string(), 100, 200).unwrap();
        let int2 = GenomicInterval::new("chr1".to_string(), 150, 250).unwrap();
        let int3 = GenomicInterval::new("chr1".to_string(), 300, 400).unwrap();
        let int4 = GenomicInterval::new("chr2".to_string(), 100, 200).unwrap();

        // Overlapping intervals
        assert!(int1.overlaps(&int2));
        assert!(int2.overlaps(&int1)); // Symmetric

        // Non-overlapping intervals
        assert!(!int1.overlaps(&int3));
        assert!(!int3.overlaps(&int1)); // Symmetric

        // Different chromosomes
        assert!(!int1.overlaps(&int4));
    }

    #[test]
    fn test_genomic_interval_contains() {
        let outer = GenomicInterval::new("chr1".to_string(), 100, 300).unwrap();
        let inner = GenomicInterval::new("chr1".to_string(), 150, 200).unwrap();
        let partial = GenomicInterval::new("chr1".to_string(), 150, 350).unwrap();
        let disjoint = GenomicInterval::new("chr1".to_string(), 400, 500).unwrap();

        // Contains
        assert!(outer.contains(&inner));
        assert!(!inner.contains(&outer));

        // Partial overlap
        assert!(!outer.contains(&partial));

        // Disjoint
        assert!(!outer.contains(&disjoint));
    }

    #[test]
    fn test_genomic_interval_display() {
        let interval = GenomicInterval::new("chr1".to_string(), 100, 200).unwrap();
        assert_eq!(interval.to_string(), "chr1:100-200");
    }

    #[test]
    fn test_strand_from_str() {
        assert_eq!(Strand::from_str("+").unwrap(), Strand::Forward);
        assert_eq!(Strand::from_str("-").unwrap(), Strand::Reverse);
        assert_eq!(Strand::from_str(".").unwrap(), Strand::Unknown);

        // Invalid strand
        let result = Strand::from_str("x");
        assert!(result.is_err());
    }

    #[test]
    fn test_strand_display() {
        assert_eq!(Strand::Forward.to_string(), "+");
        assert_eq!(Strand::Reverse.to_string(), "-");
        assert_eq!(Strand::Unknown.to_string(), ".");
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_interval_length_invariant(start in 0u64..100000, end in 0u64..100000) {
            if start < end {
                let interval = GenomicInterval::new("chr1".to_string(), start, end).unwrap();
                assert_eq!(interval.length(), end - start);
            }
        }

        #[test]
        fn test_interval_overlaps_symmetric(
            start1 in 0u64..10000,
            end1 in 0u64..10000,
            start2 in 0u64..10000,
            end2 in 0u64..10000,
        ) {
            // Ensure valid intervals
            let start1 = start1.min(end1);
            let end1 = end1.max(start1 + 1);
            let start2 = start2.min(end2);
            let end2 = end2.max(start2 + 1);

            let int1 = GenomicInterval::new("chr1".to_string(), start1, end1).unwrap();
            let int2 = GenomicInterval::new("chr1".to_string(), start2, end2).unwrap();

            // Overlaps should be symmetric
            assert_eq!(int1.overlaps(&int2), int2.overlaps(&int1));
        }

        #[test]
        fn test_interval_contains_implies_overlaps(
            start1 in 0u64..10000,
            end1 in 0u64..10000,
            start2 in 0u64..10000,
            end2 in 0u64..10000,
        ) {
            // Ensure valid intervals
            let start1 = start1.min(end1);
            let end1 = end1.max(start1 + 1);
            let start2 = start2.min(end2);
            let end2 = end2.max(start2 + 1);

            let int1 = GenomicInterval::new("chr1".to_string(), start1, end1).unwrap();
            let int2 = GenomicInterval::new("chr1".to_string(), start2, end2).unwrap();

            // If int1 contains int2, they must overlap
            if int1.contains(&int2) {
                assert!(int1.overlaps(&int2));
            }
        }
    }
}
