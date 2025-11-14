//! Shared primitives for bioinformatics file formats.
//!
//! This module provides reusable infrastructure for tab-delimited formats:
//! - Generic parsers and traits
//! - Field parsing utilities
//! - Genomic types (intervals, strands)
//! - Header/comment handling
//!
//! # Example: Using Genomic Types
//!
//! ```
//! use biometal::formats::primitives::{GenomicInterval, Strand};
//! use std::str::FromStr;
//!
//! // Create an interval
//! let interval = GenomicInterval::new("chr1".to_string(), 100, 200)?;
//! assert_eq!(interval.length(), 100);
//!
//! // Parse strand
//! let strand = Strand::from_str("+")?;
//! assert_eq!(strand, Strand::Forward);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use thiserror::Error;

pub mod genomic;
pub mod tab_delimited;

// Re-exports
pub use genomic::{GenomicInterval, Strand};
pub use tab_delimited::{TabDelimitedParser, TabDelimitedRecord};

/// Errors that can occur when parsing bioinformatics formats.
#[derive(Debug, Error)]
pub enum FormatError {
    /// Invalid number of tab-delimited fields.
    #[error("Invalid number of fields: expected {expected}, got {actual} at line {line}")]
    FieldCount {
        /// Expected number of fields
        expected: usize,
        /// Actual number of fields found
        actual: usize,
        /// Line number where error occurred
        line: usize,
    },

    /// Invalid field value.
    #[error("Invalid field '{field}' at line {line}: {reason}")]
    InvalidField {
        /// Field name
        field: String,
        /// Line number where error occurred
        line: usize,
        /// Reason for invalidity
        reason: String,
    },

    /// Parse error with context.
    #[error("Parse error at line {line}: {source}")]
    ParseError {
        /// Line number where error occurred
        line: usize,
        /// Underlying error
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    /// Invalid genomic interval (start >= end).
    #[error("Invalid genomic interval: start {start} >= end {end}")]
    InvalidInterval {
        /// Start position
        start: u64,
        /// End position
        end: u64,
    },

    /// Invalid strand specification.
    #[error("Invalid strand: {0} (expected '+', '-', or '.')")]
    InvalidStrand(String),

    /// I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// HTTP error.
    #[error("HTTP error: {0}")]
    #[cfg(feature = "network")]
    Http(#[from] reqwest::Error),
}

/// Result type for format operations.
pub type Result<T> = std::result::Result<T, FormatError>;
