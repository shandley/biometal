//! Error types for biometal

use thiserror::Error;

/// Result type alias for biometal operations
pub type Result<T> = std::result::Result<T, BiometalError>;

/// Error types that can occur in biometal
#[derive(Debug, Error)]
pub enum BiometalError {
    /// I/O error
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Invalid FASTQ format
    #[error("Invalid FASTQ format at line {line}: {msg}")]
    InvalidFastqFormat {
        /// Line number where error occurred
        line: usize,
        /// Error message
        msg: String,
    },

    /// Invalid FASTA format
    #[error("Invalid FASTA format at line {line}: {msg}")]
    InvalidFastaFormat {
        /// Line number where error occurred
        line: usize,
        /// Error message
        msg: String,
    },

    /// Compression/decompression error
    #[error("Compression error: {0}")]
    Compression(String),

    /// Paired-end read ID mismatch
    #[error("Paired-end read ID mismatch: R1={r1_id}, R2={r2_id}")]
    PairedEndMismatch {
        /// R1 read ID
        r1_id: String,
        /// R2 read ID
        r2_id: String,
    },

    /// Paired-end file length mismatch
    #[error("Paired-end files have different lengths")]
    PairedEndLengthMismatch,

    /// Invalid range or region
    #[error("Invalid range: {0}")]
    InvalidRange(String),

    /// Network error (Rule 6)
    #[cfg(feature = "network")]
    #[error("Network error: {0}")]
    Network(String),

    /// HTTP error (Rule 6)
    #[cfg(feature = "network")]
    #[error("HTTP error {status}: {url}")]
    Http {
        /// HTTP status code
        status: u16,
        /// URL that failed
        url: String,
    },

    /// Network timeout (Rule 6)
    #[cfg(feature = "network")]
    #[error("Network timeout after {seconds}s: {url}")]
    Timeout {
        /// Timeout duration in seconds
        seconds: u64,
        /// URL that timed out
        url: String,
    },

    /// Cache error
    #[cfg(feature = "network")]
    #[error("Cache error: {0}")]
    Cache(String),

    /// Metal GPU not available
    #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
    #[error("Metal GPU not available on this system")]
    MetalNotAvailable,
}
