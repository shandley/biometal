//! Paired-end FASTQ streaming (constant memory)
//!
//! # Overview
//!
//! Paired-end sequencing produces two files (R1 and R2) where records at the same
//! position represent forward and reverse reads from the same DNA fragment.
//!
//! This module provides synchronized streaming of paired FASTQ files with:
//! - Constant memory (~5 MB total for both files)
//! - Automatic mismatch detection (different record counts)
//! - Streaming architecture (no accumulation)
//!
//! # Example
//!
//! ```no_run
//! use biometal::PairedFastqStream;
//!
//! # fn main() -> biometal::Result<()> {
//! let paired = PairedFastqStream::from_paths("sample_R1.fq.gz", "sample_R2.fq.gz")?;
//!
//! for pair in paired {
//!     let (r1, r2) = pair?;
//!     // Process paired reads together
//!     println!("R1: {}, R2: {}", r1.id, r2.id);
//! }
//! # Ok(())
//! # }
//! ```

use crate::error::{BiometalError, Result};
use crate::io::{CompressedReader, DataSource, FastqStream};
use crate::types::FastqRecord;
use std::io::BufRead;
use std::path::Path;

/// Paired-end FASTQ stream iterator
///
/// Yields tuples of `(R1, R2)` records, ensuring both files are read synchronously.
/// Returns an error if the files have different numbers of records.
pub struct PairedFastqStream<R1: BufRead, R2: BufRead> {
    stream1: FastqStream<R1>,
    stream2: FastqStream<R2>,
    record_count: usize,
}

impl PairedFastqStream<CompressedReader, CompressedReader> {
    /// Create a paired stream from two file paths
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::PairedFastqStream;
    ///
    /// # fn main() -> biometal::Result<()> {
    /// let paired = PairedFastqStream::from_paths(
    ///     "sample_R1.fq.gz",
    ///     "sample_R2.fq.gz"
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_paths<P1: AsRef<Path>, P2: AsRef<Path>>(path1: P1, path2: P2) -> Result<Self> {
        let source1 = DataSource::from_path(path1);
        let source2 = DataSource::from_path(path2);
        Self::new(source1, source2)
    }

    /// Create a paired stream from two data sources
    ///
    /// This allows mixing local files and network sources.
    pub fn new(source1: DataSource, source2: DataSource) -> Result<Self> {
        let stream1 = FastqStream::new(source1)?;
        let stream2 = FastqStream::new(source2)?;
        Ok(Self::from_streams(stream1, stream2))
    }
}

impl<R1: BufRead, R2: BufRead> PairedFastqStream<R1, R2> {
    /// Create a paired stream from two existing FASTQ streams
    pub fn from_streams(stream1: FastqStream<R1>, stream2: FastqStream<R2>) -> Self {
        Self {
            stream1,
            stream2,
            record_count: 0,
        }
    }

    /// Get the number of paired records read so far
    pub fn records_read(&self) -> usize {
        self.record_count
    }
}

impl<R1: BufRead, R2: BufRead> Iterator for PairedFastqStream<R1, R2> {
    type Item = Result<(FastqRecord, FastqRecord)>;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.stream1.next(), self.stream2.next()) {
            // Both files have records - return pair
            (Some(Ok(r1)), Some(Ok(r2))) => {
                self.record_count += 1;
                Some(Ok((r1, r2)))
            }

            // Both files exhausted - normal completion
            (None, None) => None,

            // Mismatched record counts - R1 has more records
            (Some(_), None) => Some(Err(BiometalError::InvalidFastqFormat {
                line: self.record_count * 4 + 1,
                msg: format!(
                    "R1 file has more records than R2 (R1 continues after {} pairs)",
                    self.record_count
                ),
            })),

            // Mismatched record counts - R2 has more records
            (None, Some(_)) => Some(Err(BiometalError::InvalidFastqFormat {
                line: self.record_count * 4 + 1,
                msg: format!(
                    "R2 file has more records than R1 (R2 continues after {} pairs)",
                    self.record_count
                ),
            })),

            // R1 error
            (Some(Err(e)), _) => Some(Err(e)),

            // R2 error
            (_, Some(Err(e))) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn create_test_fastq(records: usize) -> Vec<u8> {
        let mut data = Vec::new();
        for i in 0..records {
            data.extend_from_slice(format!("@read_{}\n", i).as_bytes());
            data.extend_from_slice(b"ACGT\n");
            data.extend_from_slice(b"+\n");
            data.extend_from_slice(b"IIII\n");
        }
        data
    }

    #[test]
    fn test_paired_basic() {
        let r1_data = create_test_fastq(3);
        let r2_data = create_test_fastq(3);

        let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
        let stream2 = FastqStream::from_reader(Cursor::new(r2_data));

        let mut paired = PairedFastqStream::from_streams(stream1, stream2);

        // Should read 3 pairs
        let pairs: Result<Vec<_>> = paired.by_ref().collect();
        let pairs = pairs.unwrap();
        assert_eq!(pairs.len(), 3);

        // Verify record count tracking
        assert_eq!(paired.records_read(), 3);
    }

    #[test]
    fn test_paired_empty() {
        let r1_data = create_test_fastq(0);
        let r2_data = create_test_fastq(0);

        let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
        let stream2 = FastqStream::from_reader(Cursor::new(r2_data));

        let paired = PairedFastqStream::from_streams(stream1, stream2);
        let pairs: Vec<_> = paired.collect();
        assert_eq!(pairs.len(), 0);
    }

    #[test]
    fn test_paired_mismatch_r1_longer() {
        let r1_data = create_test_fastq(3); // 3 records
        let r2_data = create_test_fastq(2); // 2 records

        let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
        let stream2 = FastqStream::from_reader(Cursor::new(r2_data));

        let paired = PairedFastqStream::from_streams(stream1, stream2);
        let pairs: Vec<_> = paired.collect();

        // Should get 2 successful pairs, then 1 error
        assert_eq!(pairs.len(), 3);
        assert!(pairs[0].is_ok());
        assert!(pairs[1].is_ok());
        assert!(pairs[2].is_err());

        // Verify error message
        match &pairs[2] {
            Err(BiometalError::InvalidFastqFormat { msg, .. }) => {
                assert!(msg.contains("R1 file has more records"));
            }
            _ => panic!("Expected InvalidFastqFormat error"),
        }
    }

    #[test]
    fn test_paired_mismatch_r2_longer() {
        let r1_data = create_test_fastq(2); // 2 records
        let r2_data = create_test_fastq(3); // 3 records

        let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
        let stream2 = FastqStream::from_reader(Cursor::new(r2_data));

        let paired = PairedFastqStream::from_streams(stream1, stream2);
        let pairs: Vec<_> = paired.collect();

        // Should get 2 successful pairs, then 1 error
        assert_eq!(pairs.len(), 3);
        assert!(pairs[0].is_ok());
        assert!(pairs[1].is_ok());
        assert!(pairs[2].is_err());

        // Verify error message
        match &pairs[2] {
            Err(BiometalError::InvalidFastqFormat { msg, .. }) => {
                assert!(msg.contains("R2 file has more records"));
            }
            _ => panic!("Expected InvalidFastqFormat error"),
        }
    }

    #[test]
    fn test_paired_ids_match() {
        let r1_data = b"@read_1/1\nACGT\n+\nIIII\n@read_2/1\nACGT\n+\nIIII\n";
        let r2_data = b"@read_1/2\nACGT\n+\nIIII\n@read_2/2\nACGT\n+\nIIII\n";

        let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
        let stream2 = FastqStream::from_reader(Cursor::new(r2_data));

        let paired = PairedFastqStream::from_streams(stream1, stream2);

        for (i, pair) in paired.enumerate() {
            let (r1, r2) = pair.unwrap();
            assert_eq!(r1.id, format!("read_{}/1", i + 1));
            assert_eq!(r2.id, format!("read_{}/2", i + 1));
        }
    }

    // Property-based tests
    use proptest::prelude::*;

    proptest! {
        /// Test that paired reads with matching counts can be processed
        #[test]
        fn test_paired_matching_counts(
            records_count in 1..20usize,
        ) {
            let mut r1_data = Vec::new();
            let mut r2_data = Vec::new();

            for i in 0..records_count {
                let seq = "ACGT".repeat(25);
                let qual = "I".repeat(100);
                r1_data.extend_from_slice(format!("@read_{}/1\n{}\n+\n{}\n", i, seq, qual).as_bytes());
                r2_data.extend_from_slice(format!("@read_{}/2\n{}\n+\n{}\n", i, seq, qual).as_bytes());
            }

            let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
            let stream2 = FastqStream::from_reader(Cursor::new(r2_data));
            let paired = PairedFastqStream::from_streams(stream1, stream2);

            let pairs: Result<Vec<_>> = paired.collect();
            let pairs = pairs.unwrap();

            prop_assert_eq!(pairs.len(), records_count);
        }

        /// Test that mismatched counts are detected (R1 longer)
        #[test]
        fn test_paired_r1_longer(
            r1_count in 2..10usize,
        ) {
            let r2_count = r1_count - 1;

            let mut r1_data = Vec::new();
            let mut r2_data = Vec::new();

            for i in 0..r1_count {
                let seq = "ACGT".repeat(10);
                let qual = "I".repeat(40);
                r1_data.extend_from_slice(format!("@read_{}\n{}\n+\n{}\n", i, seq, qual).as_bytes());
            }
            for i in 0..r2_count {
                let seq = "ACGT".repeat(10);
                let qual = "I".repeat(40);
                r2_data.extend_from_slice(format!("@read_{}\n{}\n+\n{}\n", i, seq, qual).as_bytes());
            }

            let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
            let stream2 = FastqStream::from_reader(Cursor::new(r2_data));
            let paired = PairedFastqStream::from_streams(stream1, stream2);

            let results: Vec<_> = paired.collect();

            // Should have r2_count successful pairs, then 1 error
            prop_assert_eq!(results.len(), r1_count);
            prop_assert!(results.iter().take(r2_count).all(|r| r.is_ok()));
            prop_assert!(results[r2_count].is_err());
        }

        /// Test that mismatched counts are detected (R2 longer)
        #[test]
        fn test_paired_r2_longer(
            r2_count in 2..10usize,
        ) {
            let r1_count = r2_count - 1;

            let mut r1_data = Vec::new();
            let mut r2_data = Vec::new();

            for i in 0..r1_count {
                let seq = "ACGT".repeat(10);
                let qual = "I".repeat(40);
                r1_data.extend_from_slice(format!("@read_{}\n{}\n+\n{}\n", i, seq, qual).as_bytes());
            }
            for i in 0..r2_count {
                let seq = "ACGT".repeat(10);
                let qual = "I".repeat(40);
                r2_data.extend_from_slice(format!("@read_{}\n{}\n+\n{}\n", i, seq, qual).as_bytes());
            }

            let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
            let stream2 = FastqStream::from_reader(Cursor::new(r2_data));
            let paired = PairedFastqStream::from_streams(stream1, stream2);

            let results: Vec<_> = paired.collect();

            // Should have r1_count successful pairs, then 1 error
            prop_assert_eq!(results.len(), r2_count);
            prop_assert!(results.iter().take(r1_count).all(|r| r.is_ok()));
            prop_assert!(results[r1_count].is_err());
        }
    }
}
