//! BAM header parsing.
//!
//! The BAM header consists of:
//! 1. Magic bytes ("BAM\1")
//! 2. SAM header text (null-terminated string)
//! 3. Reference sequence dictionary
//!
//! # Format
//!
//! ```text
//! BAM Header:
//! - 4 bytes: Magic ("BAM\1")
//! - 4 bytes: SAM header text length (l_text, int32)
//! - l_text bytes: SAM header text
//! - 4 bytes: Number of reference sequences (n_ref, int32)
//! - For each reference:
//!   - 4 bytes: Reference name length (l_name, int32, includes null terminator)
//!   - l_name bytes: Reference name (null-terminated)
//!   - 4 bytes: Reference length (int32)
//! ```

use std::io::{self, Read};

/// BAM magic bytes.
const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

/// Reference sequence information.
///
/// Each reference sequence (chromosome/contig) has a name and length.
/// These are used to validate alignment positions and convert reference
/// IDs to names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    /// Reference sequence name (e.g., "chr1", "chrM")
    pub name: String,
    /// Reference sequence length in bases
    pub length: u32,
}

impl Reference {
    /// Create a new reference.
    pub fn new(name: String, length: u32) -> Self {
        Self { name, length }
    }
}

/// BAM file header.
///
/// Contains SAM header text and reference sequence information.
/// This metadata is required to interpret BAM records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    /// SAM header text (includes @HD, @SQ, @RG, @PG lines)
    pub text: String,
    /// Reference sequences (chromosomes/contigs)
    pub references: Vec<Reference>,
}

impl Header {
    /// Create a new header.
    pub fn new(text: String, references: Vec<Reference>) -> Self {
        Self { text, references }
    }

    /// Get reference by ID.
    ///
    /// Returns `None` if the reference ID is out of bounds.
    pub fn reference(&self, id: usize) -> Option<&Reference> {
        self.references.get(id)
    }

    /// Get reference name by ID.
    ///
    /// Returns `None` if the reference ID is out of bounds.
    pub fn reference_name(&self, id: usize) -> Option<&str> {
        self.reference(id).map(|r| r.name.as_str())
    }

    /// Get number of reference sequences.
    pub fn reference_count(&self) -> usize {
        self.references.len()
    }
}

/// Read and validate BAM magic bytes.
///
/// The first 4 bytes of a BAM file must be "BAM\1".
///
/// # Errors
///
/// Returns error if:
/// - Cannot read 4 bytes
/// - Magic bytes don't match "BAM\1"
pub fn read_magic<R: Read>(reader: &mut R) -> io::Result<()> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;

    if &magic != BAM_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Invalid BAM magic: expected {:?}, got {:?}",
                BAM_MAGIC, magic
            ),
        ));
    }

    Ok(())
}

/// Read SAM header text.
///
/// SAM header text is stored as:
/// - 4 bytes: length (int32)
/// - N bytes: UTF-8 text
///
/// # Errors
///
/// Returns error if:
/// - Cannot read length
/// - Cannot read text
/// - Text is not valid UTF-8
pub fn read_header_text<R: Read>(reader: &mut R) -> io::Result<String> {
    // Read length
    let mut len_bytes = [0u8; 4];
    reader.read_exact(&mut len_bytes)?;
    let len = i32::from_le_bytes(len_bytes);

    if len < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid SAM header length: {}", len),
        ));
    }

    let len = len as usize;

    // Read text
    let mut text_bytes = vec![0u8; len];
    reader.read_exact(&mut text_bytes)?;

    // Convert to UTF-8 string
    String::from_utf8(text_bytes).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid UTF-8 in SAM header: {}", e),
        )
    })
}

/// Read a single reference sequence.
///
/// Reference format:
/// - 4 bytes: name length (including null terminator)
/// - N bytes: name (null-terminated)
/// - 4 bytes: sequence length
///
/// # Errors
///
/// Returns error if:
/// - Cannot read name length
/// - Cannot read name
/// - Name is not valid UTF-8
/// - Cannot read sequence length
pub fn read_reference<R: Read>(reader: &mut R) -> io::Result<Reference> {
    // Read name length (includes null terminator)
    let mut name_len_bytes = [0u8; 4];
    reader.read_exact(&mut name_len_bytes)?;
    let name_len = i32::from_le_bytes(name_len_bytes);

    if name_len <= 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid reference name length: {}", name_len),
        ));
    }

    let name_len = name_len as usize;

    // Read name (null-terminated)
    let mut name_bytes = vec![0u8; name_len];
    reader.read_exact(&mut name_bytes)?;

    // Remove null terminator
    if name_bytes.last() != Some(&0) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Reference name not null-terminated",
        ));
    }
    name_bytes.pop();

    // Convert to UTF-8 string
    let name = String::from_utf8(name_bytes).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid UTF-8 in reference name: {}", e),
        )
    })?;

    // Read sequence length
    let mut len_bytes = [0u8; 4];
    reader.read_exact(&mut len_bytes)?;
    let length = i32::from_le_bytes(len_bytes);

    if length < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid reference length: {}", length),
        ));
    }

    Ok(Reference::new(name, length as u32))
}

/// Read all reference sequences.
///
/// Format:
/// - 4 bytes: number of references
/// - For each reference: read_reference()
///
/// # Errors
///
/// Returns error if:
/// - Cannot read number of references
/// - Cannot read any reference
pub fn read_references<R: Read>(reader: &mut R) -> io::Result<Vec<Reference>> {
    // Read number of references
    let mut count_bytes = [0u8; 4];
    reader.read_exact(&mut count_bytes)?;
    let count = i32::from_le_bytes(count_bytes);

    if count < 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid reference count: {}", count),
        ));
    }

    let count = count as usize;
    let mut references = Vec::with_capacity(count);

    for i in 0..count {
        let reference = read_reference(reader).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Error reading reference {}: {}", i, e),
            )
        })?;
        references.push(reference);
    }

    Ok(references)
}

/// Read complete BAM header.
///
/// Reads and validates:
/// - Magic bytes
/// - SAM header text
/// - Reference sequences
///
/// # Errors
///
/// Returns error if any component fails to parse.
///
/// # Example
///
/// ```no_run
/// use biometal::io::bam::header::read_header;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// # fn main() -> std::io::Result<()> {
/// let file = File::open("alignments.bam")?;
/// let mut reader = BufReader::new(file);
/// let header = read_header(&mut reader)?;
///
/// println!("References: {}", header.reference_count());
/// for (i, ref_seq) in header.references.iter().enumerate() {
///     println!("  {}: {} ({} bp)", i, ref_seq.name, ref_seq.length);
/// }
/// # Ok(())
/// # }
/// ```
pub fn read_header<R: Read>(reader: &mut R) -> io::Result<Header> {
    read_magic(reader)?;
    let text = read_header_text(reader)?;
    let references = read_references(reader)?;
    Ok(Header::new(text, references))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_magic_valid() {
        let data = b"BAM\x01extra data";
        let mut cursor = Cursor::new(data);
        assert!(read_magic(&mut cursor).is_ok());
    }

    #[test]
    fn test_read_magic_invalid() {
        let data = b"BAMX";
        let mut cursor = Cursor::new(data);
        assert!(read_magic(&mut cursor).is_err());
    }

    #[test]
    fn test_read_header_text_empty() {
        let data = vec![0, 0, 0, 0]; // length 0
        let mut cursor = Cursor::new(data);
        let text = read_header_text(&mut cursor).unwrap();
        assert_eq!(text, "");
    }

    #[test]
    fn test_read_header_text_simple() {
        let mut data = vec![5, 0, 0, 0]; // length 5
        data.extend_from_slice(b"@HD\tVN:1.6\n");
        // Wait, that's 11 bytes, let me fix
        let data = {
            let text = b"hello";
            let mut data = vec![5, 0, 0, 0]; // length 5
            data.extend_from_slice(text);
            data
        };
        let mut cursor = Cursor::new(data);
        let text = read_header_text(&mut cursor).unwrap();
        assert_eq!(text, "hello");
    }

    #[test]
    fn test_read_reference() {
        // Reference: "chr1" (length 5 including null) = 248956422 bp
        let mut data = vec![5, 0, 0, 0]; // name length = 5
        data.extend_from_slice(b"chr1\0"); // name with null terminator
        data.extend_from_slice(&248956422u32.to_le_bytes()); // length

        let mut cursor = Cursor::new(data);
        let reference = read_reference(&mut cursor).unwrap();
        assert_eq!(reference.name, "chr1");
        assert_eq!(reference.length, 248956422);
    }

    #[test]
    fn test_read_references_empty() {
        let data = vec![0, 0, 0, 0]; // count = 0
        let mut cursor = Cursor::new(data);
        let references = read_references(&mut cursor).unwrap();
        assert_eq!(references.len(), 0);
    }

    #[test]
    fn test_read_references_multiple() {
        let mut data = vec![2, 0, 0, 0]; // count = 2

        // Reference 1: "chr1"
        data.extend_from_slice(&5i32.to_le_bytes()); // name length
        data.extend_from_slice(b"chr1\0");
        data.extend_from_slice(&1000u32.to_le_bytes());

        // Reference 2: "chr2"
        data.extend_from_slice(&5i32.to_le_bytes()); // name length
        data.extend_from_slice(b"chr2\0");
        data.extend_from_slice(&2000u32.to_le_bytes());

        let mut cursor = Cursor::new(data);
        let references = read_references(&mut cursor).unwrap();
        assert_eq!(references.len(), 2);
        assert_eq!(references[0].name, "chr1");
        assert_eq!(references[0].length, 1000);
        assert_eq!(references[1].name, "chr2");
        assert_eq!(references[1].length, 2000);
    }

    #[test]
    fn test_read_full_header() {
        let mut data = Vec::new();

        // Magic
        data.extend_from_slice(b"BAM\x01");

        // SAM header text
        let header_text = "@HD\tVN:1.6\n";
        data.extend_from_slice(&(header_text.len() as i32).to_le_bytes());
        data.extend_from_slice(header_text.as_bytes());

        // References (1 reference)
        data.extend_from_slice(&1i32.to_le_bytes()); // count

        // Reference: "chr1"
        data.extend_from_slice(&5i32.to_le_bytes()); // name length
        data.extend_from_slice(b"chr1\0");
        data.extend_from_slice(&1000u32.to_le_bytes());

        let mut cursor = Cursor::new(data);
        let header = read_header(&mut cursor).unwrap();

        assert_eq!(header.text, "@HD\tVN:1.6\n");
        assert_eq!(header.references.len(), 1);
        assert_eq!(header.references[0].name, "chr1");
        assert_eq!(header.references[0].length, 1000);
    }

    #[test]
    fn test_header_reference_lookup() {
        let header = Header::new(
            String::from("@HD\tVN:1.6\n"),
            vec![
                Reference::new(String::from("chr1"), 1000),
                Reference::new(String::from("chr2"), 2000),
            ],
        );

        assert_eq!(header.reference_count(), 2);
        assert_eq!(header.reference_name(0), Some("chr1"));
        assert_eq!(header.reference_name(1), Some("chr2"));
        assert_eq!(header.reference_name(2), None);
    }
}
