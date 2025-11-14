//! Generic tab-delimited file parsing.
//!
//! This module provides generic infrastructure for parsing tab-delimited formats
//! like BED, GFF, VCF, and GFA. All these formats share common patterns:
//! - Tab-delimited fields
//! - Comment lines (starting with `#`)
//! - Line-based records
//!
//! # Design
//!
//! The [`TabDelimitedRecord`] trait defines the interface for parsing records.
//! The [`TabDelimitedParser`] provides a generic streaming parser that works
//! with any type implementing this trait.
//!
//! # Examples
//!
//! ```
//! use biometal::formats::primitives::{TabDelimitedRecord, TabDelimitedParser, Result};
//!
//! // Define a custom record type
//! #[derive(Debug, PartialEq)]
//! struct SimpleRecord {
//!     chrom: String,
//!     start: u64,
//!     end: u64,
//! }
//!
//! impl TabDelimitedRecord for SimpleRecord {
//!     fn from_line(line: &str) -> Result<Self> {
//!         let fields: Vec<_> = line.split('\t').collect();
//!         if fields.len() < 3 {
//!             return Err(biometal::formats::primitives::FormatError::FieldCount {
//!                 expected: 3,
//!                 actual: fields.len(),
//!                 line: 0,
//!             });
//!         }
//!
//!         Ok(SimpleRecord {
//!             chrom: fields[0].to_string(),
//!             start: fields[1].parse().unwrap(),
//!             end: fields[2].parse().unwrap(),
//!         })
//!     }
//!
//!     fn to_line(&self) -> String {
//!         format!("{}\t{}\t{}", self.chrom, self.start, self.end)
//!     }
//! }
//!
//! // Parse from string
//! let data = "chr1\t100\t200\nchr2\t300\t400\n";
//! let parser = TabDelimitedParser::<_, SimpleRecord>::new(data.as_bytes());
//!
//! let records: Vec<_> = parser.collect::<Result<_>>().unwrap();
//! assert_eq!(records.len(), 2);
//! assert_eq!(records[0].chrom, "chr1");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::formats::primitives::Result;
#[cfg(test)]
use crate::formats::primitives::FormatError;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::marker::PhantomData;
use std::path::Path;

/// Trait for types that can be parsed from tab-delimited lines.
///
/// Implement this trait to create custom parsers for tab-delimited formats.
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::{TabDelimitedRecord, Result};
///
/// #[derive(Debug, PartialEq)]
/// struct BedRecord {
///     chrom: String,
///     start: u64,
///     end: u64,
/// }
///
/// impl TabDelimitedRecord for BedRecord {
///     fn from_line(line: &str) -> Result<Self> {
///         let fields: Vec<_> = line.split('\t').collect();
///         Ok(BedRecord {
///             chrom: fields[0].to_string(),
///             start: fields[1].parse().unwrap(),
///             end: fields[2].parse().unwrap(),
///         })
///     }
///
///     fn to_line(&self) -> String {
///         format!("{}\t{}\t{}", self.chrom, self.start, self.end)
///     }
/// }
/// ```
pub trait TabDelimitedRecord: Sized {
    /// Parse a record from a tab-delimited line.
    ///
    /// The line should not include the trailing newline.
    ///
    /// # Errors
    ///
    /// Returns an error if the line is malformed or contains invalid data.
    fn from_line(line: &str) -> Result<Self>;

    /// Serialize this record to a tab-delimited line.
    ///
    /// The returned string should not include a trailing newline.
    fn to_line(&self) -> String;

    /// Expected number of tab-delimited fields.
    ///
    /// Returns `None` if the number of fields is variable.
    /// Default implementation returns `None`.
    fn expected_fields() -> Option<usize> {
        None
    }
}

/// Generic streaming parser for tab-delimited formats.
///
/// Parses records one at a time with constant memory usage.
/// Automatically skips:
/// - Empty lines
/// - Comment lines (starting with `#`)
///
/// # Type Parameters
///
/// - `R`: The underlying reader (anything implementing `Read`)
/// - `T`: The record type (must implement `TabDelimitedRecord`)
///
/// # Examples
///
/// ## Parse from file
///
/// ```no_run
/// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
///
/// # #[derive(Debug)]
/// # struct BedRecord { chrom: String, start: u64, end: u64 }
/// # impl TabDelimitedRecord for BedRecord {
/// #     fn from_line(line: &str) -> Result<Self> { todo!() }
/// #     fn to_line(&self) -> String { todo!() }
/// # }
/// # fn main() -> Result<()> {
/// let parser = TabDelimitedParser::<_, BedRecord>::from_path("data.bed")?;
///
/// for record in parser {
///     let record = record?;
///     // Process record
/// }
/// # Ok(())
/// # }
/// ```
///
/// ## Parse from compressed file
///
/// ```no_run
/// # use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
/// # #[derive(Debug)]
/// # struct BedRecord { chrom: String, start: u64, end: u64 }
/// # impl TabDelimitedRecord for BedRecord {
/// #     fn from_line(line: &str) -> Result<Self> { todo!() }
/// #     fn to_line(&self) -> String { todo!() }
/// # }
/// # fn main() -> Result<()> {
/// // Automatically decompresses with cloudflare_zlib
/// let parser = TabDelimitedParser::<_, BedRecord>::from_bgzip_path("data.bed.gz")?;
///
/// for record in parser {
///     let record = record?;
///     // Process record
/// }
/// # Ok(())
/// # }
/// ```
///
/// ## Parse from HTTP
///
/// ```no_run
/// # use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
/// # #[derive(Debug)]
/// # struct BedRecord { chrom: String, start: u64, end: u64 }
/// # impl TabDelimitedRecord for BedRecord {
/// #     fn from_line(line: &str) -> Result<Self> { todo!() }
/// #     fn to_line(&self) -> String { todo!() }
/// # }
/// # fn main() -> Result<()> {
/// let parser = TabDelimitedParser::<_, BedRecord>::from_url(
///     "https://example.com/data.bed"
/// )?;
///
/// for record in parser {
///     let record = record?;
///     // Process record
/// }
/// # Ok(())
/// # }
/// ```
pub struct TabDelimitedParser<R: Read, T: TabDelimitedRecord> {
    reader: BufReader<R>,
    line_buf: String,
    line_number: usize,
    _phantom: PhantomData<T>,
}

impl<R: Read, T: TabDelimitedRecord> TabDelimitedParser<R, T> {
    /// Creates a new parser from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct SimpleRecord { data: String }
    /// # impl TabDelimitedRecord for SimpleRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { Ok(SimpleRecord { data: line.to_string() }) }
    /// #     fn to_line(&self) -> String { self.data.clone() }
    /// # }
    /// let data = "field1\tfield2\tfield3\n";
    /// let parser = TabDelimitedParser::<_, SimpleRecord>::new(data.as_bytes());
    /// ```
    pub fn new(reader: R) -> Self {
        TabDelimitedParser {
            reader: BufReader::new(reader),
            line_buf: String::with_capacity(1024),
            line_number: 0,
            _phantom: PhantomData,
        }
    }

    /// Returns the current line number (1-based).
    ///
    /// Useful for error reporting.
    pub fn line_number(&self) -> usize {
        self.line_number
    }
}

impl<T: TabDelimitedRecord> TabDelimitedParser<File, T> {
    /// Creates a parser from a file path.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct BedRecord { chrom: String, start: u64, end: u64 }
    /// # impl TabDelimitedRecord for BedRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { todo!() }
    /// #     fn to_line(&self) -> String { todo!() }
    /// # }
    /// # fn main() -> Result<()> {
    /// let parser = TabDelimitedParser::<_, BedRecord>::from_path("data.bed")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self::new(file))
    }
}

impl<T: TabDelimitedRecord> TabDelimitedParser<MultiGzDecoder<File>, T> {
    /// Creates a parser from a gzip/bgzip-compressed file.
    ///
    /// Uses cloudflare_zlib backend for optimal decompression performance.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or is not valid gzip.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct BedRecord { chrom: String, start: u64, end: u64 }
    /// # impl TabDelimitedRecord for BedRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { todo!() }
    /// #     fn to_line(&self) -> String { todo!() }
    /// # }
    /// # fn main() -> Result<()> {
    /// let parser = TabDelimitedParser::<_, BedRecord>::from_bgzip_path("data.bed.gz")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_bgzip_path(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        let decoder = MultiGzDecoder::new(file);
        Ok(Self::new(decoder))
    }
}

#[cfg(feature = "network")]
impl<T: TabDelimitedRecord> TabDelimitedParser<Box<dyn Read>, T> {
    /// Creates a parser from an HTTP/HTTPS URL.
    ///
    /// Supports streaming from remote files without downloading.
    ///
    /// # Errors
    ///
    /// Returns an error if the URL is invalid or cannot be fetched.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use biometal::formats::primitives::{TabDelimitedParser, TabDelimitedRecord, Result};
    ///
    /// # #[derive(Debug)]
    /// # struct BedRecord { chrom: String, start: u64, end: u64 }
    /// # impl TabDelimitedRecord for BedRecord {
    /// #     fn from_line(line: &str) -> Result<Self> { todo!() }
    /// #     fn to_line(&self) -> String { todo!() }
    /// # }
    /// # fn main() -> Result<()> {
    /// let parser = TabDelimitedParser::<_, BedRecord>::from_url(
    ///     "https://example.com/data.bed"
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_url(url: &str) -> Result<Self> {
        let response = reqwest::blocking::get(url)?;
        let reader: Box<dyn Read> = Box::new(response);
        Ok(Self::new(reader))
    }
}

impl<R: Read, T: TabDelimitedRecord> Iterator for TabDelimitedParser<R, T> {
    type Item = Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buf.clear();

            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    self.line_number += 1;

                    // Trim trailing newline
                    let line = self.line_buf.trim_end();

                    // Skip empty lines
                    if line.is_empty() {
                        continue;
                    }

                    // Skip comments (lines starting with #)
                    if line.starts_with('#') {
                        continue;
                    }

                    // Parse record
                    return Some(T::from_line(line));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test record type
    #[derive(Debug, PartialEq)]
    struct TestRecord {
        chrom: String,
        start: u64,
        end: u64,
    }

    impl TabDelimitedRecord for TestRecord {
        fn from_line(line: &str) -> Result<Self> {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 3 {
                return Err(FormatError::FieldCount {
                    expected: 3,
                    actual: fields.len(),
                    line: 0,
                });
            }

            Ok(TestRecord {
                chrom: fields[0].to_string(),
                start: fields[1].parse().map_err(|e| FormatError::InvalidField {
                    field: "start".to_string(),
                    line: 0,
                    reason: format!("{}", e),
                })?,
                end: fields[2].parse().map_err(|e| FormatError::InvalidField {
                    field: "end".to_string(),
                    line: 0,
                    reason: format!("{}", e),
                })?,
            })
        }

        fn to_line(&self) -> String {
            format!("{}\t{}\t{}", self.chrom, self.start, self.end)
        }

        fn expected_fields() -> Option<usize> {
            Some(3)
        }
    }

    #[test]
    fn test_parse_basic() {
        let data = "chr1\t100\t200\nchr2\t300\t400\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].start, 100);
        assert_eq!(records[0].end, 200);
        assert_eq!(records[1].chrom, "chr2");
        assert_eq!(records[1].start, 300);
        assert_eq!(records[1].end, 400);
    }

    #[test]
    fn test_parse_skip_comments() {
        let data = "# This is a comment\nchr1\t100\t200\n# Another comment\nchr2\t300\t400\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].chrom, "chr2");
    }

    #[test]
    fn test_parse_skip_empty_lines() {
        let data = "chr1\t100\t200\n\n\nchr2\t300\t400\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
    }

    #[test]
    fn test_parse_mixed() {
        let data = "# Header\n\nchr1\t100\t200\n# Comment\n\nchr2\t300\t400\n\n# Trailer\n";
        let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].chrom, "chr2");
    }

    #[test]
    fn test_to_line_round_trip() {
        let original = TestRecord {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
        };

        let line = original.to_line();
        let parsed = TestRecord::from_line(&line).unwrap();

        assert_eq!(parsed, original);
    }

    #[test]
    fn test_line_number_tracking() {
        let data = "# Comment\nchr1\t100\t200\nchr2\t300\t400\n";
        let mut parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

        // Before first record
        assert_eq!(parser.line_number(), 0);

        // After first record (line 2, after comment)
        let _ = parser.next();
        assert_eq!(parser.line_number(), 2);

        // After second record
        let _ = parser.next();
        assert_eq!(parser.line_number(), 3);
    }
}
