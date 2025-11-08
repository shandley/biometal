//! BAM optional tags (auxiliary data).
//!
//! Optional tags store additional information about alignments such as:
//! - Edit distance (NM:i)
//! - Alignment score (AS:i)
//! - Secondary alignment status (SA:Z)
//! - Read group (RG:Z)
//! - Many others
//!
//! # Format
//!
//! Each tag is encoded as:
//! - 2 bytes: tag name (e.g., "NM")
//! - 1 byte: value type (A, i, f, Z, H, B)
//! - N bytes: value (format depends on type)
//!
//! # Phase 1 Simplification
//!
//! For Phase 1, we store tags as raw bytes and defer full parsing to Phase 6.
//! This allows us to focus on correctness of the core BAM parsing while
//! maintaining the ability to preserve tags for round-trip testing.
//!
//! # Evidence-Based Design
//!
//! This simplification follows `OPTIMIZATION_RULES.md` principles:
//! - Focus optimization effort on measured bottlenecks (BGZF: 66-80% CPU time)
//! - Defer optimizations until profiling proves value (Rule 1 threshold: >=15%)
//! - Maintain streaming architecture (Rule 5) without premature complexity

use std::io;

/// Container for BAM optional tags.
///
/// # Phase 1 Implementation
///
/// Currently stores tags as raw bytes. Full tag parsing (individual
/// tag extraction, type conversion) will be added in Phase 6.
///
/// # Future Work
///
/// - Tag iteration
/// - Type-safe value extraction
/// - Tag modification
/// - Common tag accessors (NM, AS, MD, etc.)
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Tags {
    /// Raw tag data from BAM record.
    ///
    /// This preserves all tag information for differential testing
    /// and round-trip validation without needing full tag parsing.
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
    // Phase 1: Just store raw bytes
    // Phase 6: Parse individual tags with type validation
    Ok(Tags::from_raw(data.to_vec()))
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
}
