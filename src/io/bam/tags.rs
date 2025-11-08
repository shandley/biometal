//! BAM optional tags (auxiliary data).
//!
//! Optional tags store additional information about alignments such as:
//! - Edit distance (NM:i)
//! - Alignment score (AS:i)
//! - Secondary alignment status (SA:Z)
//! - Read group (RG:Z)
//! - MD string (MD:Z)
//! - Many others
//!
//! # Format
//!
//! Each tag is encoded as:
//! - 2 bytes: tag name (e.g., "NM")
//! - 1 byte: value type (A, i, f, Z, H, B)
//! - N bytes: value (format depends on type)
//!
//! # Tag Types
//!
//! - **A**: Printable character (1 byte)
//! - **i**: Signed integer (c, C, s, S, i, I - variable width)
//! - **f**: Float (4 bytes, IEEE 754)
//! - **Z**: Null-terminated string
//! - **H**: Hex string (stored as string)
//! - **B**: Array of numeric values
//!
//! # Phase 4 Implementation
//!
//! Phase 4 completes tag parsing with type-safe accessors, enabling:
//! - Individual tag extraction by name
//! - Type-safe value access
//! - Tag iteration
//! - Common tag convenience methods

use std::io;
use std::fmt;
use super::error::BamDecodeError;

/// Tag value types in BAM format.
///
/// Each tag has a type code and corresponding value encoding.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue {
    /// Character (A): Single printable character
    Char(u8),
    /// Integer (i): Signed integer (variable width)
    Int(i64),
    /// Float (f): IEEE 754 single-precision float
    Float(f32),
    /// String (Z): Null-terminated string
    String(String),
    /// Hex string (H): Hex-encoded data
    Hex(String),
    /// Array (B): Typed array of numbers
    Array(ArrayValue),
}

/// Array value types for tag arrays (B type).
#[derive(Debug, Clone, PartialEq)]
pub enum ArrayValue {
    /// Array of signed 8-bit integers
    Int8(Vec<i8>),
    /// Array of unsigned 8-bit integers
    UInt8(Vec<u8>),
    /// Array of signed 16-bit integers
    Int16(Vec<i16>),
    /// Array of unsigned 16-bit integers
    UInt16(Vec<u16>),
    /// Array of signed 32-bit integers
    Int32(Vec<i32>),
    /// Array of unsigned 32-bit integers
    UInt32(Vec<u32>),
    /// Array of 32-bit floats
    Float(Vec<f32>),
}

/// A single BAM tag with name and value.
#[derive(Debug, Clone, PartialEq)]
pub struct Tag {
    /// Two-character tag name (e.g., "NM", "AS", "RG")
    pub name: [u8; 2],
    /// Tag value
    pub value: TagValue,
}

impl Tag {
    /// Get tag name as a string slice.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tag;
    /// # use biometal::io::bam::TagValue;
    /// let tag = Tag {
    ///     name: *b"NM",
    ///     value: TagValue::Int(5),
    /// };
    /// assert_eq!(tag.name_str(), "NM");
    /// ```
    pub fn name_str(&self) -> &str {
        std::str::from_utf8(&self.name).unwrap_or("??")
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:", self.name_str())?;
        match &self.value {
            TagValue::Char(c) => write!(f, "A:{}", *c as char),
            TagValue::Int(i) => write!(f, "i:{}", i),
            TagValue::Float(fl) => write!(f, "f:{}", fl),
            TagValue::String(s) => write!(f, "Z:{}", s),
            TagValue::Hex(h) => write!(f, "H:{}", h),
            TagValue::Array(arr) => {
                write!(f, "B:")?;
                match arr {
                    ArrayValue::Int8(v) => write!(f, "c,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::UInt8(v) => write!(f, "C,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::Int16(v) => write!(f, "s,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::UInt16(v) => write!(f, "S,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::Int32(v) => write!(f, "i,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::UInt32(v) => write!(f, "I,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::Float(v) => write!(f, "f,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                }
            }
        }
    }
}

/// Container for BAM optional tags.
///
/// # Phase 4 Implementation
///
/// Now provides full tag parsing with type-safe accessors:
/// - Individual tag extraction by name
/// - Tag iteration
/// - Type-safe value access
/// - Common tag convenience methods
///
/// # Example
///
/// ```
/// # use biometal::io::bam::Tags;
/// # use biometal::io::bam::TagValue;
/// // Create tags from raw BAM data
/// let tags = Tags::from_raw(vec![
///     // NM:i:5 (edit distance)
///     b'N', b'M', b'i', 5, 0, 0, 0,
/// ]);
///
/// // Access tag by name
/// if let Some(tag) = tags.get(b"NM").unwrap() {
///     if let TagValue::Int(nm) = tag.value {
///         assert_eq!(nm, 5);
///     }
/// }
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Tags {
    /// Raw tag data from BAM record.
    ///
    /// Preserves all tag information for round-trip validation.
    /// Tags are parsed on-demand when accessed.
    data: Vec<u8>,
}

impl Tags {
    /// Create empty tags.
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// Create tags from raw BAM data.
    ///
    /// # Arguments
    ///
    /// * `data` - Raw tag bytes from BAM record
    ///
    /// # Phase 1 Note
    ///
    /// Currently just stores the raw bytes. Full parsing will be
    /// added in Phase 6 when we need tag-level access.
    pub fn from_raw(data: Vec<u8>) -> Self {
        Self { data }
    }

    /// Get the raw tag data.
    ///
    /// Useful for:
    /// - Differential testing (compare raw bytes)
    /// - Round-trip testing (preserve exact encoding)
    /// - Custom tag parsers
    pub fn as_raw(&self) -> &[u8] {
        &self.data
    }

    /// Check if tags are empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get the total size of tag data in bytes.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Get a specific tag by name.
    ///
    /// Returns `Ok(Some(tag))` if found, `Ok(None)` if not found, or `Err` if parsing fails.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::{Tags, TagValue};
    /// let tags = Tags::from_raw(vec![b'N', b'M', b'i', 5, 0, 0, 0]);
    ///
    /// let nm_tag = tags.get(b"NM").unwrap().unwrap();
    /// assert_eq!(nm_tag.name_str(), "NM");
    /// if let TagValue::Int(val) = nm_tag.value {
    ///     assert_eq!(val, 5);
    /// }
    /// ```
    pub fn get(&self, name: &[u8; 2]) -> io::Result<Option<Tag>> {
        let mut cursor = 0;

        while cursor < self.data.len() {
            // Read tag name (2 bytes)
            if cursor + 2 > self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Incomplete tag name at offset {}", cursor),
                ));
            }
            let tag_name = [self.data[cursor], self.data[cursor + 1]];
            cursor += 2;

            // Read tag type (1 byte)
            if cursor >= self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Missing tag type at offset {}", cursor),
                ));
            }
            let tag_type = self.data[cursor];
            cursor += 1;

            // Parse value based on type
            let (value, value_size) = parse_tag_value(&self.data[cursor..], tag_type)?;

            // Check if this is the tag we're looking for
            if &tag_name == name {
                return Ok(Some(Tag {
                    name: tag_name,
                    value,
                }));
            }

            cursor += value_size;
        }

        Ok(None)
    }

    /// Iterate over all tags.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tags;
    /// let tags = Tags::from_raw(vec![
    ///     b'N', b'M', b'i', 5, 0, 0, 0,
    ///     b'A', b'S', b'i', 100, 0, 0, 0,
    /// ]);
    ///
    /// for tag in tags.iter().unwrap() {
    ///     println!("{}: {:?}", tag.name_str(), tag.value);
    /// }
    /// ```
    pub fn iter(&self) -> io::Result<Vec<Tag>> {
        let mut tags = Vec::new();
        let mut seen_tags = std::collections::HashSet::new();
        let mut cursor = 0;

        while cursor < self.data.len() {
            // Read tag name (2 bytes)
            if cursor + 2 > self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Incomplete tag name at offset {}", cursor),
                ));
            }
            let tag_name = [self.data[cursor], self.data[cursor + 1]];
            cursor += 2;

            // Check for duplicate tags (BAM spec violation)
            if !seen_tags.insert(tag_name) {
                return Err(BamDecodeError::DuplicateTag { tag: tag_name }.into());
            }

            // Read tag type (1 byte)
            if cursor >= self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Missing tag type at offset {}", cursor),
                ));
            }
            let tag_type = self.data[cursor];
            cursor += 1;

            // Parse value
            let (value, value_size) = parse_tag_value(&self.data[cursor..], tag_type)?;
            cursor += value_size;

            tags.push(Tag {
                name: tag_name,
                value,
            });
        }

        Ok(tags)
    }

    /// Get an integer tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is an integer.
    pub fn get_i32(&self, name: &[u8; 2]) -> io::Result<Option<i32>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::Int(i) => Ok(Some(i as i32)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }

    /// Get a string tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is a string.
    pub fn get_string(&self, name: &[u8; 2]) -> io::Result<Option<String>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::String(s) => Ok(Some(s)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }
}

/// Parse tags from BAM record data.
///
/// # Phase 1 Implementation
///
/// Currently just copies the raw bytes. Full tag parsing will be
/// added in Phase 6.
///
/// # Arguments
///
/// * `data` - Raw tag data from end of BAM record
///
/// # Returns
///
/// Tags structure containing the raw data.
///
/// # Future Work
///
/// - Validate tag format (2-char name, type code, value)
/// - Parse individual tags into structured types
/// - Provide type-safe accessors
///
/// # Example
///
/// ```
/// use biometal::io::bam::Tags;
///
/// // Empty tags
/// let data = vec![];
/// let tags = Tags::from_raw(data);
/// assert!(tags.is_empty());
/// ```
pub fn parse_tags(data: &[u8]) -> io::Result<Tags> {
    // Phase 4: Store raw bytes (parsing happens on-demand via get/iter)
    Ok(Tags::from_raw(data.to_vec()))
}

/// Parse a single tag value from raw bytes.
///
/// Returns `(TagValue, bytes_consumed)`.
///
/// # Tag Type Codes
///
/// - `A` (65): Character (1 byte)
/// - `c` (99): Int8 (1 byte signed)
/// - `C` (67): UInt8 (1 byte unsigned)
/// - `s` (115): Int16 (2 bytes signed)
/// - `S` (83): UInt16 (2 bytes unsigned)
/// - `i` (105): Int32 (4 bytes signed)
/// - `I` (73): UInt32 (4 bytes unsigned)
/// - `f` (102): Float (4 bytes)
/// - `Z` (90): Null-terminated string
/// - `H` (72): Hex string (null-terminated)
/// - `B` (66): Typed array
fn parse_tag_value(data: &[u8], type_code: u8) -> io::Result<(TagValue, usize)> {
    match type_code {
        b'A' => {
            // Character (1 byte)
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for character tag",
                ));
            }
            Ok((TagValue::Char(data[0]), 1))
        }
        b'c' => {
            // Int8 (1 byte signed)
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for int8 tag",
                ));
            }
            Ok((TagValue::Int(data[0] as i8 as i64), 1))
        }
        b'C' => {
            // UInt8 (1 byte unsigned)
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for uint8 tag",
                ));
            }
            Ok((TagValue::Int(data[0] as i64), 1))
        }
        b's' => {
            // Int16 (2 bytes signed, little-endian)
            if data.len() < 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for int16 tag",
                ));
            }
            let value = i16::from_le_bytes([data[0], data[1]]);
            Ok((TagValue::Int(value as i64), 2))
        }
        b'S' => {
            // UInt16 (2 bytes unsigned, little-endian)
            if data.len() < 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for uint16 tag",
                ));
            }
            let value = u16::from_le_bytes([data[0], data[1]]);
            Ok((TagValue::Int(value as i64), 2))
        }
        b'i' => {
            // Int32 (4 bytes signed, little-endian)
            if data.len() < 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for int32 tag",
                ));
            }
            let value = i32::from_le_bytes([data[0], data[1], data[2], data[3]]);
            Ok((TagValue::Int(value as i64), 4))
        }
        b'I' => {
            // UInt32 (4 bytes unsigned, little-endian)
            if data.len() < 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for uint32 tag",
                ));
            }
            let value = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
            Ok((TagValue::Int(value as i64), 4))
        }
        b'f' => {
            // Float (4 bytes, little-endian)
            if data.len() < 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for float tag",
                ));
            }
            let value = f32::from_le_bytes([data[0], data[1], data[2], data[3]]);
            Ok((TagValue::Float(value), 4))
        }
        b'Z' | b'H' => {
            // Null-terminated string or hex string
            let null_pos = data.iter().position(|&b| b == 0).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Missing null terminator in string tag",
                )
            })?;

            let string = String::from_utf8(data[..null_pos].to_vec()).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid UTF-8 in string tag")
            })?;

            let value = if type_code == b'Z' {
                TagValue::String(string)
            } else {
                TagValue::Hex(string)
            };

            Ok((value, null_pos + 1))
        }
        b'B' => {
            // Typed array
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for array tag type",
                ));
            }

            let array_type = data[0];
            if data.len() < 5 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for array tag count",
                ));
            }

            let count_u32 = u32::from_le_bytes([data[1], data[2], data[3], data[4]]);
            let count = usize::try_from(count_u32)
                .map_err(|_| BamDecodeError::ArrayCountOverflow { count: count_u32 })?;
            let mut offset = 5;

            let array_value = match array_type {
                b'c' => {
                    // Int8 array
                    if data.len() < offset + count {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for int8 array",
                        ));
                    }
                    let values: Vec<i8> = data[offset..offset + count]
                        .iter()
                        .map(|&b| b as i8)
                        .collect();
                    offset += count;
                    ArrayValue::Int8(values)
                }
                b'C' => {
                    // UInt8 array
                    if data.len() < offset + count {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for uint8 array",
                        ));
                    }
                    let values = data[offset..offset + count].to_vec();
                    offset += count;
                    ArrayValue::UInt8(values)
                }
                b's' => {
                    // Int16 array
                    let bytes_needed = count * 2;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for int16 array",
                        ));
                    }
                    let values: Vec<i16> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 2;
                            i16::from_le_bytes([data[idx], data[idx + 1]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::Int16(values)
                }
                b'S' => {
                    // UInt16 array
                    let bytes_needed = count * 2;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for uint16 array",
                        ));
                    }
                    let values: Vec<u16> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 2;
                            u16::from_le_bytes([data[idx], data[idx + 1]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::UInt16(values)
                }
                b'i' => {
                    // Int32 array
                    let bytes_needed = count * 4;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for int32 array",
                        ));
                    }
                    let values: Vec<i32> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 4;
                            i32::from_le_bytes([data[idx], data[idx + 1], data[idx + 2], data[idx + 3]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::Int32(values)
                }
                b'I' => {
                    // UInt32 array
                    let bytes_needed = count * 4;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for uint32 array",
                        ));
                    }
                    let values: Vec<u32> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 4;
                            u32::from_le_bytes([data[idx], data[idx + 1], data[idx + 2], data[idx + 3]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::UInt32(values)
                }
                b'f' => {
                    // Float array
                    let bytes_needed = count * 4;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for float array",
                        ));
                    }
                    let values: Vec<f32> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 4;
                            f32::from_le_bytes([data[idx], data[idx + 1], data[idx + 2], data[idx + 3]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::Float(values)
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Invalid array type code: {}", array_type),
                    ));
                }
            };

            Ok((TagValue::Array(array_value), offset))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid tag type code: {}", type_code),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_tags() {
        let tags = Tags::new();
        assert!(tags.is_empty());
        assert_eq!(tags.len(), 0);
        assert_eq!(tags.as_raw(), &[] as &[u8]);
    }

    #[test]
    fn test_from_raw() {
        let data = vec![0x4E, 0x4D, 0x69, 0x01, 0x00, 0x00, 0x00]; // NM:i:1
        let tags = Tags::from_raw(data.clone());
        assert!(!tags.is_empty());
        assert_eq!(tags.len(), 7);
        assert_eq!(tags.as_raw(), &data);
    }

    #[test]
    fn test_parse_tags_empty() {
        let data = vec![];
        let tags = parse_tags(&data).unwrap();
        assert!(tags.is_empty());
    }

    #[test]
    fn test_parse_tags_preserves_data() {
        let data = vec![0x4E, 0x4D, 0x69, 0x01, 0x00, 0x00, 0x00]; // NM:i:1
        let tags = parse_tags(&data).unwrap();
        assert_eq!(tags.as_raw(), &data);
    }

    #[test]
    fn test_tags_clone() {
        let tags1 = Tags::from_raw(vec![1, 2, 3]);
        let tags2 = tags1.clone();
        assert_eq!(tags1, tags2);
    }

    #[test]
    fn test_tags_equality() {
        let tags1 = Tags::from_raw(vec![1, 2, 3]);
        let tags2 = Tags::from_raw(vec![1, 2, 3]);
        let tags3 = Tags::from_raw(vec![1, 2, 4]);

        assert_eq!(tags1, tags2);
        assert_ne!(tags1, tags3);
    }

    // Noodles-inspired robustness tests

    #[test]
    fn test_duplicate_tags_rejected() {
        let data = vec![
            b'N', b'M', b'i', 5, 0, 0, 0,  // NM:i:5
            b'N', b'M', b'i', 3, 0, 0, 0,  // NM:i:3 (duplicate!)
        ];
        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Duplicate tag"));
    }

    #[test]
    fn test_array_count_overflow() {
        // This test is primarily for 32-bit platforms, but good for completeness
        #[cfg(target_pointer_width = "32")]
        {
            let data = vec![
                b'T', b'S', b'B', b'I',  // TS:B:I (uint32 array)
                0xFF, 0xFF, 0xFF, 0xFF,  // count = u32::MAX (overflows usize on 32-bit)
            ];
            let tags = Tags::from_raw(data);
            let result = tags.iter();
            assert!(result.is_err());
            assert!(result.unwrap_err().to_string().contains("too large"));
        }
    }

    #[test]
    fn test_string_tag_missing_nul() {
        let data = vec![
            b'R', b'G', b'Z',  // RG:Z:...
            b'r', b'g', b'0',  // "rg0" (no NUL terminator!)
        ];
        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("null terminator"));
    }

    #[test]
    fn test_truncated_array() {
        let data = vec![
            b'T', b'S', b'B', b'I',  // TS:B:I (int32 array)
            0x03, 0x00, 0x00, 0x00,  // count = 3
            0x01, 0x00, 0x00, 0x00,  // [0] = 1
            0x02, 0x00, // [1] = incomplete! (only 2 bytes of 4)
        ];
        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Insufficient"));
    }
}
