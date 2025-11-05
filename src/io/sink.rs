//! Output destinations for streaming writes
//!
//! This module provides the `DataSink` abstraction, which is the write
//! counterpart to `DataSource`. It enables biometal to write to various
//! destinations while maintaining the same API.
//!
//! # Architecture
//!
//! `DataSink` mirrors `DataSource`:
//! - `DataSource::Local` (read) ↔ `DataSink::Local` (write)
//! - `DataSource::Http` (read) ↔ (future: `DataSink::S3` for write)
//! - (future: `DataSource::Sra`) ↔ N/A (SRA is read-only)
//!
//! # Example
//!
//! ```no_run
//! use biometal::io::DataSink;
//! use std::path::Path;
//!
//! // Write to local file
//! let sink = DataSink::from_path("output.fq.gz");
//!
//! // Write to stdout
//! let sink = DataSink::stdout();
//! ```

use std::path::{Path, PathBuf};

/// Output destination for streaming writes
///
/// This enum abstracts over different output destinations, allowing
/// writers to be agnostic to where data is being written.
///
/// # Memory Guarantees
///
/// All sinks maintain constant memory usage regardless of data volume.
#[derive(Debug, Clone)]
pub enum DataSink {
    /// Write to a local file path
    ///
    /// Compression format is auto-detected from file extension:
    /// - `.gz` → gzip compression
    /// - `.bgz` → bgzip compression (recommended for bioinformatics)
    /// - other → uncompressed
    Local(PathBuf),

    /// Write to standard output
    ///
    /// Useful for streaming pipelines:
    /// ```bash
    /// biometal filter input.fq.gz | biometal qc
    /// ```
    Stdout,

    // Future: Network destinations
    // S3(String),
    // Azure(String),
    // GCS(String),
}

impl DataSink {
    /// Create a sink from a file path
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::io::DataSink;
    ///
    /// let sink = DataSink::from_path("output.fq.gz");
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> Self {
        Self::Local(path.as_ref().to_path_buf())
    }

    /// Create a sink for standard output
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::io::DataSink;
    ///
    /// let sink = DataSink::stdout();
    /// ```
    pub fn stdout() -> Self {
        Self::Stdout
    }

    /// Get the file extension if this is a local file sink
    ///
    /// Used internally for compression auto-detection.
    pub(crate) fn extension(&self) -> Option<&str> {
        match self {
            Self::Local(path) => path.extension().and_then(|s| s.to_str()),
            Self::Stdout => None,
        }
    }

    /// Check if this sink represents a compressed output
    pub fn is_compressed(&self) -> bool {
        matches!(
            self.extension(),
            Some("gz") | Some("bgz") | Some("gzip")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_path() {
        let sink = DataSink::from_path("test.fq");
        match sink {
            DataSink::Local(path) => {
                assert_eq!(path, PathBuf::from("test.fq"));
            }
            _ => panic!("Expected Local variant"),
        }
    }

    #[test]
    fn test_stdout() {
        let sink = DataSink::stdout();
        matches!(sink, DataSink::Stdout);
    }

    #[test]
    fn test_extension_detection() {
        let sink = DataSink::from_path("test.fq.gz");
        assert_eq!(sink.extension(), Some("gz"));

        let sink = DataSink::from_path("test.fq.bgz");
        assert_eq!(sink.extension(), Some("bgz"));

        let sink = DataSink::from_path("test.fq");
        assert_eq!(sink.extension(), Some("fq"));

        let sink = DataSink::stdout();
        assert_eq!(sink.extension(), None);
    }

    #[test]
    fn test_is_compressed() {
        assert!(DataSink::from_path("test.fq.gz").is_compressed());
        assert!(DataSink::from_path("test.fq.bgz").is_compressed());
        assert!(!DataSink::from_path("test.fq").is_compressed());
        assert!(!DataSink::stdout().is_compressed());
    }
}
