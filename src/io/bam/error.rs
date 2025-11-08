//! BAM parsing error types.
//!
//! Provides structured error types for BAM parsing that enable:
//! - Precise error identification
//! - Pattern matching on specific error cases
//! - Better debugging with contextual information
//!
//! Inspired by noodles-bam's error hierarchy for production-grade error handling.

use std::{error, fmt, io};

/// Errors that can occur during BAM decoding.
///
/// This enum provides specific error variants for different parsing failures,
/// enabling precise error handling and better debugging messages.
///
/// # Example
///
/// ```
/// use biometal::io::bam::error::BamDecodeError;
///
/// fn handle_error(err: BamDecodeError) {
///     match err {
///         BamDecodeError::DuplicateTag { tag } => {
///             eprintln!("Found duplicate tag: {}{}", tag[0] as char, tag[1] as char);
///         }
///         BamDecodeError::InvalidReferenceId { value, field } => {
///             eprintln!("Invalid {} reference ID: {}", field, value);
///         }
///         _ => eprintln!("Other error: {}", err),
///     }
/// }
/// ```
#[derive(Debug)]
pub enum BamDecodeError {
    /// I/O error occurred during reading
    Io(io::Error),

    /// Invalid reference sequence ID (must be -1 or >= 0)
    InvalidReferenceId {
        /// The invalid reference ID value
        value: i32,
        /// Which field had the invalid ID ("read" or "mate")
        field: String,
    },

    /// Invalid alignment position
    InvalidPosition {
        /// The invalid position value
        value: i32,
        /// Offset in the data where error occurred
        offset: usize,
    },

    /// Invalid read name length (must be >= 1)
    InvalidReadNameLength {
        /// The invalid length value
        length: u8,
        /// Offset in the data where error occurred
        offset: usize,
    },

    /// Missing NUL terminator in string field
    MissingNulTerminator {
        /// Which field was missing the terminator
        field: String,
        /// Offset in the data where error occurred
        offset: usize,
    },

    /// Invalid UTF-8 in string field
    InvalidUtf8 {
        /// Which field had invalid UTF-8
        field: String,
        /// The underlying UTF-8 error
        source: std::string::FromUtf8Error,
    },

    /// Invalid tag type code
    InvalidTagType {
        /// The tag name
        tag: [u8; 2],
        /// The invalid type code
        type_code: u8,
    },

    /// Duplicate tag found (spec violation)
    DuplicateTag {
        /// The duplicate tag name
        tag: [u8; 2],
    },

    /// Array count too large for platform (overflow on 32-bit)
    ArrayCountOverflow {
        /// The count value that overflowed
        count: u32,
    },

    /// Invalid array subtype
    InvalidArraySubtype {
        /// The invalid subtype code
        subtype: u8,
    },

    /// CIGAR operation count too large
    CigarCountOverflow {
        /// The count value
        count: usize,
    },

    /// Invalid CIGAR operation code
    InvalidCigarOp {
        /// The invalid operation value
        value: u32,
    },

    /// Negative sequence length (invalid)
    NegativeSequenceLength {
        /// The negative length value
        length: i32,
    },

    /// Insufficient data for parsing
    UnexpectedEof {
        /// What was being parsed
        context: String,
        /// Expected number of bytes
        expected: usize,
        /// Actual number of bytes available
        actual: usize,
    },

    /// Invalid BAM magic bytes
    InvalidMagic {
        /// The actual bytes found
        actual: [u8; 4],
    },

    /// Generic invalid data error with context
    InvalidData {
        /// Description of what was invalid
        message: String,
    },
}

impl error::Error for BamDecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidUtf8 { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl fmt::Display for BamDecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {}", e),

            Self::InvalidReferenceId { value, field } => {
                write!(
                    f,
                    "Invalid {} reference ID: {} (must be -1 or >= 0)",
                    field, value
                )
            }

            Self::InvalidPosition { value, offset } => {
                write!(f, "Invalid position at offset {}: {}", offset, value)
            }

            Self::InvalidReadNameLength { length, offset } => {
                write!(
                    f,
                    "Invalid read name length at offset {}: {} (must be >= 1)",
                    offset, length
                )
            }

            Self::MissingNulTerminator { field, offset } => {
                write!(f, "Missing NUL terminator in {} at offset {}", field, offset)
            }

            Self::InvalidUtf8 { field, source } => {
                write!(f, "Invalid UTF-8 in {}: {}", field, source)
            }

            Self::InvalidTagType { tag, type_code } => {
                write!(
                    f,
                    "Invalid tag type for {}{}:: {}",
                    tag[0] as char, tag[1] as char, type_code
                )
            }

            Self::DuplicateTag { tag } => {
                write!(f, "Duplicate tag: {}{}", tag[0] as char, tag[1] as char)
            }

            Self::ArrayCountOverflow { count } => {
                write!(
                    f,
                    "Array count too large for platform: {} (exceeds usize::MAX)",
                    count
                )
            }

            Self::InvalidArraySubtype { subtype } => {
                write!(f, "Invalid array subtype: {}", subtype)
            }

            Self::CigarCountOverflow { count } => {
                write!(f, "CIGAR operation count too large: {}", count)
            }

            Self::InvalidCigarOp { value } => {
                write!(f, "Invalid CIGAR operation: {}", value)
            }

            Self::NegativeSequenceLength { length } => {
                write!(f, "Invalid negative sequence length: {}", length)
            }

            Self::UnexpectedEof {
                context,
                expected,
                actual,
            } => {
                write!(
                    f,
                    "Unexpected end of file while parsing {}: expected {} bytes, got {}",
                    context, expected, actual
                )
            }

            Self::InvalidMagic { actual } => {
                write!(
                    f,
                    "Invalid BAM magic bytes: expected [BAM\\x01], got {:?}",
                    actual
                )
            }

            Self::InvalidData { message } => {
                write!(f, "Invalid data: {}", message)
            }
        }
    }
}

impl From<io::Error> for BamDecodeError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<BamDecodeError> for io::Error {
    fn from(e: BamDecodeError) -> Self {
        match e {
            BamDecodeError::Io(io_err) => io_err,
            other => io::Error::new(io::ErrorKind::InvalidData, other),
        }
    }
}
