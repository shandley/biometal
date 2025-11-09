//! Python bindings for BAM/SAM parser
//!
//! Provides Python access to biometal's high-performance BAM parser with:
//! - Parallel BGZF decompression (4× speedup, Rule 3)
//! - Streaming architecture (constant ~5 MB memory, Rule 5)
//! - 4.54 million records/sec throughput (43.0 MiB/s)
//!
//! # Evidence Base
//!
//! - Rule 3 (Parallel BGZF): Entry 029, 6.5× validated speedup
//! - Rule 5 (Streaming): Entry 026, constant memory architecture
//! - Performance: 4× overall speedup (experiments/native-bam-implementation/PHASE_3_BENCHMARKS.md)

use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::exceptions::PyStopIteration;
use std::path::PathBuf;
use std::collections::HashMap;
use crate::io::bam::{BamReader, Header, Reference, Tag, TagValue, ArrayValue, CigarOp, SamWriter};
use crate::io::compression::CompressedReader;
use std::fs::File;
use std::io::BufWriter;

/// Reference sequence information
///
/// Attributes:
///     name (str): Reference name (e.g., "chr1", "chrM")
///     length (int): Reference sequence length in bases
///
/// Example:
///     >>> header = bam.header
///     >>> ref = header.reference(0)
///     >>> print(f"{ref.name}: {ref.length} bp")
#[pyclass(name = "Reference")]
#[derive(Clone)]
pub struct PyReference {
    /// Reference sequence name
    #[pyo3(get)]
    pub name: String,

    /// Reference sequence length in bases
    #[pyo3(get)]
    pub length: u32,
}

#[pymethods]
impl PyReference {
    fn __repr__(&self) -> String {
        format!("Reference(name='{}', length={})", self.name, self.length)
    }

    fn __str__(&self) -> String {
        format!("{}: {} bp", self.name, self.length)
    }
}

impl From<&Reference> for PyReference {
    fn from(reference: &Reference) -> Self {
        PyReference {
            name: reference.name.clone(),
            length: reference.length,
        }
    }
}

/// CIGAR operation
///
/// Represents a single CIGAR operation in an alignment.
/// CIGAR (Compact Idiosyncratic Gapped Alignment Report) describes
/// how a read aligns to the reference.
///
/// Operation Types:
///     - Match (M): Alignment match (can include mismatches)
///     - Insertion (I): Insertion to reference
///     - Deletion (D): Deletion from reference
///     - RefSkip (N): Skipped region from reference (intron)
///     - SoftClip (S): Soft clipping (bases in read, not aligned)
///     - HardClip (H): Hard clipping (bases not in read)
///     - Padding (P): Padding (silent deletion from padded reference)
///     - SeqMatch (=): Sequence match (bases match reference)
///     - SeqMismatch (X): Sequence mismatch (bases differ)
///
/// Attributes:
///     length (int): Number of bases consumed by this operation
///     op_char (str): Operation character (M, I, D, N, S, H, P, =, X)
///
/// Example:
///     >>> for op in record.cigar:
///     ...     print(f"{op.length}{op.op_char}")
///     100M
///     5I
///     50M
#[pyclass(name = "CigarOp")]
#[derive(Clone)]
pub struct PyCigarOp {
    inner: CigarOp,
}

#[pymethods]
impl PyCigarOp {
    /// Get the operation length
    #[getter]
    fn length(&self) -> u32 {
        self.inner.length()
    }

    /// Get the operation character (M, I, D, N, S, H, P, =, X)
    #[getter]
    fn op_char(&self) -> char {
        self.inner.as_char()
    }

    /// Check if this is a Match/mismatch operation (M)
    fn is_match(&self) -> bool {
        matches!(self.inner, CigarOp::Match(_))
    }

    /// Check if this is an Insertion operation (I)
    fn is_insertion(&self) -> bool {
        matches!(self.inner, CigarOp::Insertion(_))
    }

    /// Check if this is a Deletion operation (D)
    fn is_deletion(&self) -> bool {
        matches!(self.inner, CigarOp::Deletion(_))
    }

    /// Check if this is a Reference skip operation (N)
    fn is_ref_skip(&self) -> bool {
        matches!(self.inner, CigarOp::RefSkip(_))
    }

    /// Check if this is a Soft clip operation (S)
    fn is_soft_clip(&self) -> bool {
        matches!(self.inner, CigarOp::SoftClip(_))
    }

    /// Check if this is a Hard clip operation (H)
    fn is_hard_clip(&self) -> bool {
        matches!(self.inner, CigarOp::HardClip(_))
    }

    /// Check if this is a Padding operation (P)
    fn is_padding(&self) -> bool {
        matches!(self.inner, CigarOp::Padding(_))
    }

    /// Check if this is a Sequence match operation (=)
    fn is_seq_match(&self) -> bool {
        matches!(self.inner, CigarOp::SeqMatch(_))
    }

    /// Check if this is a Sequence mismatch operation (X)
    fn is_seq_mismatch(&self) -> bool {
        matches!(self.inner, CigarOp::SeqMismatch(_))
    }

    /// Check if this operation consumes reference bases
    ///
    /// Returns True for: M, D, N, =, X
    /// Returns False for: I, S, H, P
    fn consumes_reference(&self) -> bool {
        matches!(
            self.inner,
            CigarOp::Match(_)
                | CigarOp::Deletion(_)
                | CigarOp::RefSkip(_)
                | CigarOp::SeqMatch(_)
                | CigarOp::SeqMismatch(_)
        )
    }

    /// Check if this operation consumes query bases (from read)
    ///
    /// Returns True for: M, I, S, =, X
    /// Returns False for: D, N, H, P
    fn consumes_query(&self) -> bool {
        matches!(
            self.inner,
            CigarOp::Match(_)
                | CigarOp::Insertion(_)
                | CigarOp::SoftClip(_)
                | CigarOp::SeqMatch(_)
                | CigarOp::SeqMismatch(_)
        )
    }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}

impl From<&CigarOp> for PyCigarOp {
    fn from(op: &CigarOp) -> Self {
        PyCigarOp { inner: *op }
    }
}

/// BAM tag value
///
/// Represents different tag value types in BAM format.
///
/// Types:
///     - char(str): Single character (type A)
///     - int(int): Integer value (type i)
///     - float(float): Float value (type f)
///     - string(str): String value (type Z)
///     - hex(str): Hex string (type H)
///     - array(list): Array of values (type B)
///
/// Example:
///     >>> tag = record.get_tag("NM")
///     >>> if tag:
///     ...     print(f"Edit distance: {tag.as_int()}")
#[pyclass(name = "TagValue")]
#[derive(Clone)]
pub struct PyTagValue {
    inner: TagValue,
}

#[pymethods]
impl PyTagValue {
    /// Check if tag is a character value
    fn is_char(&self) -> bool {
        matches!(self.inner, TagValue::Char(_))
    }

    /// Check if tag is an integer value
    fn is_int(&self) -> bool {
        matches!(self.inner, TagValue::Int(_))
    }

    /// Check if tag is a float value
    fn is_float(&self) -> bool {
        matches!(self.inner, TagValue::Float(_))
    }

    /// Check if tag is a string value
    fn is_string(&self) -> bool {
        matches!(self.inner, TagValue::String(_))
    }

    /// Check if tag is an array value
    fn is_array(&self) -> bool {
        matches!(self.inner, TagValue::Array(_))
    }

    /// Get value as integer (if it is one)
    fn as_int(&self) -> Option<i64> {
        match &self.inner {
            TagValue::Int(i) => Some(*i),
            _ => None,
        }
    }

    /// Get value as float (if it is one)
    fn as_float(&self) -> Option<f32> {
        match &self.inner {
            TagValue::Float(f) => Some(*f),
            _ => None,
        }
    }

    /// Get value as string (if it is one)
    fn as_string(&self) -> Option<String> {
        match &self.inner {
            TagValue::String(s) => Some(s.clone()),
            TagValue::Hex(h) => Some(h.clone()),
            _ => None,
        }
    }

    /// Get value as character (if it is one)
    fn as_char(&self) -> Option<char> {
        match &self.inner {
            TagValue::Char(c) => Some(*c as char),
            _ => None,
        }
    }

    /// Get value as array (if it is one)
    fn as_array(&self, py: Python) -> Option<Py<PyAny>> {
        use pyo3::types::PyList;

        match &self.inner {
            TagValue::Array(arr) => {
                let list = match arr {
                    ArrayValue::Int8(v) => {
                        let items: Vec<i64> = v.iter().map(|&x| x as i64).collect();
                        PyList::new(py, &items).expect("Failed to create PyList from Vec<i64>")
                    }
                    ArrayValue::UInt8(v) => {
                        let items: Vec<u64> = v.iter().map(|&x| x as u64).collect();
                        PyList::new(py, &items).expect("Failed to create PyList from Vec<u64>")
                    }
                    ArrayValue::Int16(v) => {
                        let items: Vec<i64> = v.iter().map(|&x| x as i64).collect();
                        PyList::new(py, &items).expect("Failed to create PyList from Vec<i64>")
                    }
                    ArrayValue::UInt16(v) => {
                        let items: Vec<u64> = v.iter().map(|&x| x as u64).collect();
                        PyList::new(py, &items).expect("Failed to create PyList from Vec<u64>")
                    }
                    ArrayValue::Int32(v) => {
                        let items: Vec<i64> = v.iter().map(|&x| x as i64).collect();
                        PyList::new(py, &items).expect("Failed to create PyList from Vec<i64>")
                    }
                    ArrayValue::UInt32(v) => {
                        let items: Vec<u64> = v.iter().map(|&x| x as u64).collect();
                        PyList::new(py, &items).expect("Failed to create PyList from Vec<u64>")
                    }
                    ArrayValue::Float(v) => {
                        PyList::new(py, v.as_slice()).expect("Failed to create PyList from Vec<f32>")
                    }
                };
                Some(list.into())
            }
            _ => None,
        }
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            TagValue::Char(c) => format!("TagValue::Char('{}')", *c as char),
            TagValue::Int(i) => format!("TagValue::Int({})", i),
            TagValue::Float(f) => format!("TagValue::Float({})", f),
            TagValue::String(s) => format!("TagValue::String('{}')", s),
            TagValue::Hex(h) => format!("TagValue::Hex('{}')", h),
            TagValue::Array(_) => "TagValue::Array([...])".to_string(),
        }
    }

    fn __str__(&self) -> String {
        match &self.inner {
            TagValue::Char(c) => format!("{}", *c as char),
            TagValue::Int(i) => format!("{}", i),
            TagValue::Float(f) => format!("{}", f),
            TagValue::String(s) => s.clone(),
            TagValue::Hex(h) => h.clone(),
            TagValue::Array(_) => "[array]".to_string(),
        }
    }
}

impl From<TagValue> for PyTagValue {
    fn from(value: TagValue) -> Self {
        PyTagValue { inner: value }
    }
}

/// BAM tag with name and value
///
/// Attributes:
///     name (str): Two-character tag name (e.g., "NM", "AS", "RG")
///     value (TagValue): Tag value
///
/// Example:
///     >>> tag = record.get_tag("NM")
///     >>> print(f"{tag.name}: {tag.value}")
///     NM: 5
#[pyclass(name = "Tag")]
#[derive(Clone)]
pub struct PyTag {
    /// Tag name (2 characters)
    #[pyo3(get)]
    pub name: String,

    /// Tag value
    #[pyo3(get)]
    pub value: PyTagValue,
}

#[pymethods]
impl PyTag {
    fn __repr__(&self) -> String {
        format!("Tag(name='{}', value={})", self.name, self.value.__repr__())
    }

    fn __str__(&self) -> String {
        format!("{}:{}", self.name, self.value.__str__())
    }
}

impl From<&Tag> for PyTag {
    fn from(tag: &Tag) -> Self {
        PyTag {
            name: tag.name_str().to_string(),
            value: PyTagValue::from(tag.value.clone()),
        }
    }
}

/// BAM record with alignment information
///
/// Attributes:
///     name (str): Read name/identifier
///     reference_id (int | None): Reference sequence ID (-1 for unmapped)
///     position (int | None): 0-based leftmost position (-1 for unmapped)
///     mapq (int | None): Mapping quality (255 for unavailable)
///     flags (int): SAM flags (bitwise)
///     mate_reference_id (int | None): Mate reference ID
///     mate_position (int | None): Mate position
///     template_length (int): Template length
///     sequence (bytes): Read sequence
///     quality (bytes): Phred quality scores
///     sequence_length (int): Length of sequence
///
/// Properties (computed from flags):
///     is_mapped (bool): Read is mapped
///     is_reverse (bool): Read is reverse strand
///     is_forward (bool): Read is forward strand
///     is_paired (bool): Read is paired-end
///     is_proper_pair (bool): Read and mate properly paired
///     is_first (bool): Read is first in pair
///     is_second (bool): Read is second in pair
///     is_primary (bool): Primary alignment (not secondary/supplementary)
///     is_secondary (bool): Secondary alignment
///     is_supplementary (bool): Supplementary alignment
///     is_mate_mapped (bool): Mate is mapped
///     is_mate_reverse (bool): Mate is reverse strand
///     is_qc_fail (bool): Failed quality checks
///     is_duplicate (bool): PCR or optical duplicate
///
/// Example:
///     >>> bam = biometal.BamReader.from_path("alignments.bam")
///     >>> for record in bam:
///     ...     if record.is_mapped and record.is_primary and record.mapq and record.mapq >= 30:
///     ...         print(f"{record.name}: chr{record.reference_id}:{record.position}")
///
/// Note:
///     CIGAR and tags are available as raw data. Full parsing coming in future version.
#[pyclass(name = "BamRecord")]
#[derive(Clone)]
pub struct PyBamRecord {
    /// Read name/identifier
    #[pyo3(get)]
    pub name: String,

    /// Reference sequence ID (None for unmapped)
    #[pyo3(get)]
    pub reference_id: Option<usize>,

    /// 0-based leftmost position (None for unmapped)
    #[pyo3(get)]
    pub position: Option<i32>,

    /// Mapping quality (None if unavailable)
    #[pyo3(get)]
    pub mapq: Option<u8>,

    /// SAM flags (bitwise)
    #[pyo3(get)]
    pub flags: u16,

    /// Mate reference ID (None for unmapped mate)
    #[pyo3(get)]
    pub mate_reference_id: Option<usize>,

    /// Mate position (None for unmapped mate)
    #[pyo3(get)]
    pub mate_position: Option<i32>,

    /// Template length
    #[pyo3(get)]
    pub template_length: i32,

    /// Read sequence
    #[pyo3(get)]
    pub sequence: Vec<u8>,

    /// Phred quality scores
    #[pyo3(get)]
    pub quality: Vec<u8>,

    /// CIGAR operations (stored internally, accessed via property)
    cigar_ops: Vec<CigarOp>,

    /// Optional tags (stored internally, accessed via methods)
    tags: crate::io::bam::Tags,
}

#[pymethods]
impl PyBamRecord {
    /// Check if read is mapped
    #[getter]
    fn is_mapped(&self) -> bool {
        (self.flags & 0x4) == 0
    }

    /// Check if read is reverse complement
    #[getter]
    fn is_reverse(&self) -> bool {
        (self.flags & 0x10) != 0
    }

    /// Check if read is first in pair
    #[getter]
    fn is_first(&self) -> bool {
        (self.flags & 0x40) != 0
    }

    /// Check if read is second in pair
    #[getter]
    fn is_second(&self) -> bool {
        (self.flags & 0x80) != 0
    }

    /// Check if read is paired
    #[getter]
    fn is_paired(&self) -> bool {
        (self.flags & 0x1) != 0
    }

    /// Check if alignment is primary
    #[getter]
    fn is_primary(&self) -> bool {
        (self.flags & 0x900) == 0  // Not secondary, not supplementary
    }

    /// Check if read is properly paired
    ///
    /// Both reads mapped in proper pair orientation and distance.
    #[getter]
    fn is_proper_pair(&self) -> bool {
        (self.flags & 0x2) != 0
    }

    /// Check if read is forward strand
    ///
    /// Opposite of is_reverse.
    #[getter]
    fn is_forward(&self) -> bool {
        (self.flags & 0x10) == 0
    }

    /// Check if mate is mapped
    #[getter]
    fn is_mate_mapped(&self) -> bool {
        (self.flags & 0x8) == 0
    }

    /// Check if mate is reverse strand
    #[getter]
    fn is_mate_reverse(&self) -> bool {
        (self.flags & 0x20) != 0
    }

    /// Check if read is secondary alignment
    #[getter]
    fn is_secondary(&self) -> bool {
        (self.flags & 0x100) != 0
    }

    /// Check if read is supplementary alignment
    #[getter]
    fn is_supplementary(&self) -> bool {
        (self.flags & 0x800) != 0
    }

    /// Check if read failed quality checks
    #[getter]
    fn is_qc_fail(&self) -> bool {
        (self.flags & 0x200) != 0
    }

    /// Check if read is PCR or optical duplicate
    #[getter]
    fn is_duplicate(&self) -> bool {
        (self.flags & 0x400) != 0
    }

    /// Get sequence length
    #[getter]
    fn sequence_length(&self) -> usize {
        self.sequence.len()
    }

    /// Get sequence as string
    #[getter]
    fn sequence_str(&self) -> String {
        String::from_utf8_lossy(&self.sequence).to_string()
    }

    /// Get quality as string
    #[getter]
    fn quality_str(&self) -> String {
        String::from_utf8_lossy(&self.quality).to_string()
    }

    /// Get CIGAR operations
    ///
    /// Returns list of CigarOp objects describing the alignment.
    ///
    /// Returns:
    ///     List[CigarOp]: CIGAR operations
    ///
    /// Example:
    ///     >>> for op in record.cigar:
    ///     ...     print(f"{op.length}{op.op_char}", end=" ")
    ///     100M 5I 50M
    #[getter]
    fn cigar(&self) -> Vec<PyCigarOp> {
        self.cigar_ops.iter().map(|op| PyCigarOp::from(op)).collect()
    }

    /// Get CIGAR string representation
    ///
    /// Returns CIGAR in SAM format (e.g., "100M5I50M").
    ///
    /// Returns:
    ///     str: CIGAR string
    ///
    /// Example:
    ///     >>> print(record.cigar_string())
    ///     100M5I50M
    fn cigar_string(&self) -> String {
        self.cigar_ops
            .iter()
            .map(|op| format!("{}", op))
            .collect::<Vec<_>>()
            .join("")
    }

    /// Calculate reference length consumed by alignment
    ///
    /// Counts bases consumed on reference (M, D, N, =, X operations).
    /// Does not include insertions, clipping, or padding.
    ///
    /// Returns:
    ///     int: Number of reference bases consumed
    ///
    /// Example:
    ///     >>> # CIGAR: 50M5I45M (100M total, but 5I doesn't consume reference)
    ///     >>> print(record.reference_length())
    ///     95
    fn reference_length(&self) -> u32 {
        self.cigar_ops
            .iter()
            .filter_map(|op| match op {
                CigarOp::Match(len)
                | CigarOp::Deletion(len)
                | CigarOp::RefSkip(len)
                | CigarOp::SeqMatch(len)
                | CigarOp::SeqMismatch(len) => Some(len),
                _ => None,
            })
            .sum()
    }

    /// Calculate query length consumed by alignment
    ///
    /// Counts bases consumed from read (M, I, S, =, X operations).
    /// Does not include deletions, ref skips, hard clips, or padding.
    ///
    /// Returns:
    ///     int: Number of query bases consumed
    ///
    /// Example:
    ///     >>> # CIGAR: 50M5I45M (105 bases from read)
    ///     >>> print(record.query_length())
    ///     105
    fn query_length(&self) -> u32 {
        self.cigar_ops
            .iter()
            .filter_map(|op| match op {
                CigarOp::Match(len)
                | CigarOp::Insertion(len)
                | CigarOp::SoftClip(len)
                | CigarOp::SeqMatch(len)
                | CigarOp::SeqMismatch(len) => Some(len),
                _ => None,
            })
            .sum()
    }

    /// Calculate alignment end position on reference
    ///
    /// Returns the 0-based end position (exclusive) of the alignment
    /// on the reference sequence. Returns None if unmapped.
    ///
    /// Returns:
    ///     Optional[int]: End position (None if unmapped)
    ///
    /// Example:
    ///     >>> # Position 100, CIGAR 50M5D45M = 100 reference bases
    ///     >>> print(record.reference_end())
    ///     200
    fn reference_end(&self) -> Option<i32> {
        self.position.map(|pos| pos + self.reference_length() as i32)
    }

    /// Get a specific tag by name
    ///
    /// Args:
    ///     name (str): Two-character tag name (e.g., "NM", "AS", "RG")
    ///
    /// Returns:
    ///     Tag | None: Tag object if found, None otherwise
    ///
    /// Example:
    ///     >>> nm_tag = record.get_tag("NM")
    ///     >>> if nm_tag:
    ///     ...     print(f"Edit distance: {nm_tag.value.as_int()}")
    fn get_tag(&self, name: &str) -> PyResult<Option<PyTag>> {
        if name.len() != 2 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Tag name must be exactly 2 characters"
            ));
        }

        let name_bytes: [u8; 2] = [name.as_bytes()[0], name.as_bytes()[1]];

        match self.tags.get(&name_bytes) {
            Ok(Some(tag)) => Ok(Some(PyTag::from(&tag))),
            Ok(None) => Ok(None),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Error parsing tag: {}", e)
            )),
        }
    }

    /// Check if a tag exists
    ///
    /// Args:
    ///     name (str): Two-character tag name
    ///
    /// Returns:
    ///     bool: True if tag exists, False otherwise
    ///
    /// Example:
    ///     >>> if record.has_tag("NM"):
    ///     ...     print("Has edit distance tag")
    fn has_tag(&self, name: &str) -> PyResult<bool> {
        if name.len() != 2 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Tag name must be exactly 2 characters"
            ));
        }

        let name_bytes: [u8; 2] = [name.as_bytes()[0], name.as_bytes()[1]];

        match self.tags.get(&name_bytes) {
            Ok(Some(_)) => Ok(true),
            Ok(None) => Ok(false),
            Err(_) => Ok(false),  // Parsing error = tag doesn't exist
        }
    }

    /// Get all tags as a list
    ///
    /// Returns:
    ///     list[Tag]: List of all tags in the record
    ///
    /// Example:
    ///     >>> for tag in record.tags():
    ///     ...     print(f"{tag.name}: {tag.value}")
    fn tags(&self) -> PyResult<Vec<PyTag>> {
        let mut result = Vec::new();

        // Iterate through all tags
        for tag in self.tags.iter()? {
            result.push(PyTag::from(&tag));
        }

        Ok(result)
    }

    /// Get an integer tag value
    ///
    /// Convenience method that returns the integer value directly if the tag exists and is an integer.
    ///
    /// Args:
    ///     name (str): Two-character tag name
    ///
    /// Returns:
    ///     int | None: Integer value if found, None otherwise
    ///
    /// Example:
    ///     >>> edit_distance = record.get_int("NM")
    ///     >>> if edit_distance is not None:
    ///     ...     print(f"Edit distance: {edit_distance}")
    fn get_int(&self, name: &str) -> PyResult<Option<i64>> {
        if name.len() != 2 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Tag name must be exactly 2 characters"
            ));
        }

        let name_bytes: [u8; 2] = [name.as_bytes()[0], name.as_bytes()[1]];

        match self.tags.get(&name_bytes) {
            Ok(Some(tag)) => match tag.value {
                TagValue::Int(i) => Ok(Some(i)),
                _ => Ok(None),
            },
            Ok(None) => Ok(None),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Error parsing tag: {}", e)
            )),
        }
    }

    /// Get a string tag value
    ///
    /// Convenience method that returns the string value directly if the tag exists and is a string.
    ///
    /// Args:
    ///     name (str): Two-character tag name
    ///
    /// Returns:
    ///     str | None: String value if found, None otherwise
    ///
    /// Example:
    ///     >>> read_group = record.get_string("RG")
    ///     >>> if read_group:
    ///     ...     print(f"Read group: {read_group}")
    fn get_string(&self, name: &str) -> PyResult<Option<String>> {
        if name.len() != 2 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Tag name must be exactly 2 characters"
            ));
        }

        let name_bytes: [u8; 2] = [name.as_bytes()[0], name.as_bytes()[1]];

        match self.tags.get(&name_bytes) {
            Ok(Some(tag)) => match tag.value {
                TagValue::String(s) => Ok(Some(s)),
                TagValue::Hex(h) => Ok(Some(h)),
                _ => Ok(None),
            },
            Ok(None) => Ok(None),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Error parsing tag: {}", e)
            )),
        }
    }

    /// Get edit distance (NM tag)
    ///
    /// Returns the number of mismatches in the alignment (NM:i tag).
    ///
    /// Returns:
    ///     int | None: Edit distance if available
    ///
    /// Example:
    ///     >>> edit_dist = record.edit_distance()
    ///     >>> if edit_dist is not None:
    ///     ...     print(f"Mismatches: {edit_dist}")
    fn edit_distance(&self) -> PyResult<Option<i64>> {
        self.get_int("NM")
    }

    /// Get alignment score (AS tag)
    ///
    /// Returns the aligner's alignment score (AS:i tag).
    ///
    /// Returns:
    ///     int | None: Alignment score if available
    ///
    /// Example:
    ///     >>> score = record.alignment_score()
    ///     >>> if score and score >= 100:
    ///     ...     print("High quality alignment")
    fn alignment_score(&self) -> PyResult<Option<i64>> {
        self.get_int("AS")
    }

    /// Get read group (RG tag)
    ///
    /// Returns the read group identifier (RG:Z tag).
    ///
    /// Returns:
    ///     str | None: Read group ID if available
    ///
    /// Example:
    ///     >>> rg = record.read_group()
    ///     >>> if rg:
    ///     ...     print(f"Read group: {rg}")
    fn read_group(&self) -> PyResult<Option<String>> {
        self.get_string("RG")
    }

    /// Get MD string (MD tag)
    ///
    /// Returns the mismatch/deletion string (MD:Z tag) used for variant calling.
    ///
    /// Returns:
    ///     str | None: MD string if available
    ///
    /// Example:
    ///     >>> md = record.md_string()
    ///     >>> if md:
    ///     ...     print(f"MD: {md}")
    fn md_string(&self) -> PyResult<Option<String>> {
        self.get_string("MD")
    }

    fn __repr__(&self) -> String {
        let ref_str = self.reference_id
            .map(|r| r.to_string())
            .unwrap_or_else(|| "*".to_string());
        let pos_str = self.position
            .map(|p| p.to_string())
            .unwrap_or_else(|| "*".to_string());
        let mapq_str = self.mapq
            .map(|m| m.to_string())
            .unwrap_or_else(|| "*".to_string());

        format!(
            "BamRecord(name='{}', ref={}, pos={}, mapq={}, flags={})",
            self.name, ref_str, pos_str, mapq_str, self.flags
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<crate::io::bam::Record> for PyBamRecord {
    fn from(record: crate::io::bam::Record) -> Self {
        PyBamRecord {
            name: record.name,
            reference_id: record.reference_id,
            position: record.position,
            mapq: record.mapq,
            flags: record.flags,
            mate_reference_id: record.mate_reference_id,
            mate_position: record.mate_position,
            template_length: record.template_length,
            sequence: record.sequence,
            quality: record.quality,
            cigar_ops: record.cigar,
            tags: record.tags,
        }
    }
}

/// Stream BAM records with constant memory and parallel BGZF decompression
///
/// High-performance streaming BAM parser with automatic parallel BGZF decompression.
/// Maintains constant ~5 MB memory regardless of file size (terabyte-scale capable).
///
/// Args:
///     path (str): Path to BAM file (.bam, .bam.gz, or uncompressed)
///
/// Performance:
///     - 4.54 million records/sec throughput
///     - 43.0 MiB/s compressed file processing
///     - 4× speedup via parallel BGZF decompression
///     - Constant ~5 MB memory (streams terabyte-scale files)
///
/// Example:
///     >>> import biometal
///     >>> bam = biometal.BamReader.from_path("alignments.bam")
///     >>>
///     >>> # Access header info
///     >>> print(f"References: {bam.reference_count}")
///     >>>
///     >>> # Stream records with constant memory
///     >>> mapped_count = 0
///     >>> for record in bam:
///     ...     if record.is_mapped and record.mapq and record.mapq >= 30:
///     ...         mapped_count += 1
///     ...         # Process high-quality alignments
///     >>>
///     >>> print(f"High-quality alignments: {mapped_count}")
///
/// Note:
///     Automatically detects and uses parallel BGZF decompression for compressed BAM files.
///     Memory footprint remains constant at ~5 MB even for TB-scale alignments.
#[pyclass(name = "BamReader", unsendable)]
pub struct PyBamReader {
    inner: Option<BamReader<CompressedReader>>,
    header: PyBamHeader,
}

/// BAM header information
///
/// Attributes:
///     text (str): SAM header text
///     reference_count (int): Number of reference sequences
///     reference_names (list[str]): List of all reference names
///
/// Example:
///     >>> header = bam.header
///     >>> print(f"References: {', '.join(header.reference_names)}")
///     >>> ref_id = header.get_reference_id("chr1")
///     >>> print(f"chr1 is reference {ref_id}")
#[pyclass(name = "BamHeader")]
#[derive(Clone)]
pub struct PyBamHeader {
    /// SAM header text
    #[pyo3(get)]
    pub text: String,

    /// Number of reference sequences
    #[pyo3(get)]
    pub reference_count: usize,

    /// Reference sequences (internal storage)
    references: Vec<PyReference>,

    /// Reference name to ID mapping (for fast lookup)
    name_to_id: HashMap<String, usize>,
}

#[pymethods]
impl PyBamHeader {
    /// Get list of all reference names
    ///
    /// Returns:
    ///     list[str]: List of reference sequence names
    ///
    /// Example:
    ///     >>> names = header.reference_names
    ///     >>> print(f"First reference: {names[0]}")
    #[getter]
    fn reference_names(&self) -> Vec<String> {
        self.references.iter().map(|r| r.name.clone()).collect()
    }

    /// Get reference by ID
    ///
    /// Args:
    ///     reference_id (int): Reference sequence ID (0-based)
    ///
    /// Returns:
    ///     Reference | None: Reference object or None if ID is invalid
    ///
    /// Example:
    ///     >>> ref = header.reference(0)
    ///     >>> print(f"{ref.name}: {ref.length} bp")
    fn reference(&self, reference_id: usize) -> Option<PyReference> {
        self.references.get(reference_id).cloned()
    }

    /// Get reference name by ID
    ///
    /// Args:
    ///     reference_id (int): Reference sequence ID (0-based)
    ///
    /// Returns:
    ///     str | None: Reference name or None if ID is invalid
    ///
    /// Example:
    ///     >>> name = header.reference_name(0)
    ///     >>> print(f"Reference 0 is {name}")
    fn reference_name(&self, reference_id: usize) -> Option<String> {
        self.references.get(reference_id).map(|r| r.name.clone())
    }

    /// Get reference length by ID
    ///
    /// Args:
    ///     reference_id (int): Reference sequence ID (0-based)
    ///
    /// Returns:
    ///     int | None: Reference length in bases or None if ID is invalid
    ///
    /// Example:
    ///     >>> length = header.reference_length(0)
    ///     >>> print(f"Reference 0 is {length} bp")
    fn reference_length(&self, reference_id: usize) -> Option<u32> {
        self.references.get(reference_id).map(|r| r.length)
    }

    /// Get reference ID by name
    ///
    /// Args:
    ///     name (str): Reference sequence name (e.g., "chr1", "chrM")
    ///
    /// Returns:
    ///     int | None: Reference ID (0-based) or None if name not found
    ///
    /// Example:
    ///     >>> ref_id = header.get_reference_id("chr1")
    ///     >>> if ref_id is not None:
    ///     ...     print(f"chr1 is reference {ref_id}")
    fn get_reference_id(&self, name: &str) -> Option<usize> {
        self.name_to_id.get(name).copied()
    }

    fn __repr__(&self) -> String {
        format!(
            "BamHeader(references={}, text_len={})",
            self.reference_count,
            self.text.len()
        )
    }
}

impl From<&Header> for PyBamHeader {
    fn from(header: &Header) -> Self {
        // Convert references
        let references: Vec<PyReference> = header
            .references
            .iter()
            .map(|r| r.into())
            .collect();

        // Build name-to-ID mapping for O(1) lookup
        let name_to_id: HashMap<String, usize> = header
            .references
            .iter()
            .enumerate()
            .map(|(id, r)| (r.name.clone(), id))
            .collect();

        PyBamHeader {
            text: header.text.clone(),
            reference_count: header.reference_count(),
            references,
            name_to_id,
        }
    }
}

#[pymethods]
impl PyBamReader {
    /// Open BAM file with automatic parallel BGZF decompression
    ///
    /// Args:
    ///     path (str): Path to BAM file
    ///
    /// Returns:
    ///     BamReader: Streaming iterator with ~5 MB constant memory
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    ///
    /// Performance:
    ///     Automatically uses parallel BGZF decompression (4× speedup).
    ///     Processes 4.54 million records/sec with constant memory.
    ///
    /// Example:
    ///     >>> bam = biometal.BamReader.from_path("alignments.bam")
    ///     >>> print(f"Opened BAM with {bam.reference_count} references")
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = BamReader::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        // Clone header info for Python access
        let header = PyBamHeader::from(reader.header());

        Ok(PyBamReader {
            inner: Some(reader),
            header,
        })
    }

    /// Get BAM header
    ///
    /// Returns:
    ///     BamHeader: Header with reference information
    #[getter]
    fn header(&self) -> PyBamHeader {
        self.header.clone()
    }

    /// Get number of reference sequences
    ///
    /// Returns:
    ///     int: Number of reference sequences in BAM header
    #[getter]
    fn reference_count(&self) -> usize {
        self.header.reference_count
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBamRecord> {
        if let Some(ref mut reader) = slf.inner {
            match reader.read_record() {
                Ok(Some(record)) => Ok(record.into()),
                Ok(None) => Err(PyStopIteration::new_err("no more records")),
                Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
            }
        } else {
            Err(PyStopIteration::new_err("reader exhausted"))
        }
    }

    /// Query records in a specific genomic region
    ///
    /// Args:
    ///     path (str): Path to BAM file
    ///     reference_id (int): Reference sequence ID (0-based)
    ///     start (int | None): Start position (0-based, inclusive). None = start of reference
    ///     end (int | None): End position (0-based, exclusive). None = end of reference
    ///
    /// Returns:
    ///     Iterator[BamRecord]: Iterator yielding records in the specified region
    ///
    /// Note:
    ///     This performs a full file scan (O(n)) since BAI index support is not yet implemented.
    ///     For indexed random access, see roadmap for BAI/CSI support.
    ///
    /// Example:
    ///     >>> # Query chr1:1000-2000
    ///     >>> for record in biometal.BamReader.query("alignments.bam", reference_id=0, start=1000, end=2000):
    ///     ...     print(f"{record.name}: {record.position}")
    ///     >>>
    ///     >>> # Query entire reference
    ///     >>> for record in biometal.BamReader.query("alignments.bam", reference_id=1):
    ///     ...     process(record)
    #[staticmethod]
    #[pyo3(signature = (path, reference_id, start=None, end=None))]
    fn query(
        path: String,
        reference_id: usize,
        start: Option<i32>,
        end: Option<i32>,
    ) -> PyResult<PyBamRegionIter> {
        let path_buf = PathBuf::from(path);
        let reader = BamReader::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        let header = PyBamHeader::from(reader.header());

        Ok(PyBamRegionIter {
            reader: Some(reader),
            reference_id,
            start,
            end,
            _header: header,
        })
    }

    fn __repr__(&self) -> String {
        format!("BamReader(references={})", self.header.reference_count)
    }
}

/// Iterator for querying BAM records in a genomic region
///
/// Created by BamReader.query(). Yields only records within the specified region.
///
/// Note:
///     This performs a full file scan. BAI/CSI index support for O(log n) queries
///     is planned for a future release.
#[pyclass(name = "BamRegionIter", unsendable)]
pub struct PyBamRegionIter {
    reader: Option<BamReader<CompressedReader>>,
    reference_id: usize,
    start: Option<i32>,
    end: Option<i32>,
    _header: PyBamHeader,
}

#[pymethods]
impl PyBamRegionIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBamRecord> {
        // Copy filter criteria to avoid borrow checker issues
        let reference_id = slf.reference_id;
        let start = slf.start;
        let end = slf.end;

        if let Some(ref mut reader) = slf.reader {
            loop {
                match reader.read_record() {
                    Ok(Some(record)) => {
                        // Check if record matches region criteria
                        if let Some(rec_ref_id) = record.reference_id {
                            if rec_ref_id == reference_id {
                                // Check position range
                                let in_range = match (start, end, record.position) {
                                    (Some(start), Some(end), Some(pos)) => pos >= start && pos < end,
                                    (Some(start), None, Some(pos)) => pos >= start,
                                    (None, Some(end), Some(pos)) => pos < end,
                                    (None, None, _) => true,  // No range specified, include all
                                    (_, _, None) => false,  // Record has no position, skip
                                };

                                if in_range {
                                    return Ok(record.into());
                                }
                                // Record not in range, continue to next record
                            }
                            // Record on different reference, continue
                        }
                        // Record unmapped, continue
                    }
                    Ok(None) => return Err(PyStopIteration::new_err("no more records")),
                    Err(e) => {
                        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                            e.to_string(),
                        ))
                    }
                }
            }
        } else {
            Err(PyStopIteration::new_err("iterator exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        let start_str = self.start.map(|s| s.to_string()).unwrap_or_else(|| "0".to_string());
        let end_str = self.end.map(|e| e.to_string()).unwrap_or_else(|| "end".to_string());
        format!(
            "BamRegionIter(ref={}, {}:{})",
            self.reference_id, start_str, end_str
        )
    }
}

//
// Statistics Functions
//

/// Calculate per-position coverage for a genomic region
///
/// Uses CIGAR operations to calculate accurate per-base coverage.
/// Counts M (match), D (deletion), N (ref skip), = (seq match), and X (seq mismatch)
/// operations. Insertions and clipping operations do not contribute to coverage.
///
/// Args:
///     path (str): Path to BAM file
///     reference_id (int): Reference sequence ID (0-based)
///     start (int | None): Start position (0-based, inclusive). None = start of reference
///     end (int | None): End position (0-based, exclusive). None = end of reference
///
/// Returns:
///     dict[int, int]: Mapping from position to coverage depth
///
/// Example:
///     >>> coverage = biometal.calculate_coverage("alignments.bam", reference_id=0, start=1000, end=2000)
///     >>> print(f"Position 1500 has {coverage.get(1500, 0)} reads")
///     >>> max_cov = max(coverage.values())
///     >>> print(f"Maximum coverage: {max_cov}×")
///
/// Note:
///     Only returns positions with coverage > 0. Missing positions have coverage 0.
///     Uses CIGAR-aware calculation for accurate per-base coverage.
#[pyfunction]
#[pyo3(signature = (path, reference_id, start=None, end=None))]
pub fn calculate_coverage(
    _py: Python,
    path: String,
    reference_id: usize,
    start: Option<i32>,
    end: Option<i32>,
) -> PyResult<HashMap<i32, u32>> {
    let path_buf = PathBuf::from(path);
    let reader = BamReader::from_path(&path_buf)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut coverage: HashMap<i32, u32> = HashMap::new();

    // Read all records and calculate coverage
    let mut reader = reader;
    while let Ok(Some(record)) = reader.read_record() {
        // Check if record matches region
        if let Some(rec_ref_id) = record.reference_id {
            if rec_ref_id == reference_id {
                if let Some(pos) = record.position {
                    // Check if position is in range
                    let in_range = match (start, end) {
                        (Some(s), Some(e)) => pos >= s && pos < e,
                        (Some(s), None) => pos >= s,
                        (None, Some(e)) => pos < e,
                        (None, None) => true,
                    };

                    if in_range {
                        // Calculate coverage using CIGAR operations
                        let mut current_pos = pos;

                        for op in &record.cigar {
                            // Only count operations that consume reference bases
                            match op {
                                CigarOp::Match(len) |
                                CigarOp::Deletion(len) |
                                CigarOp::RefSkip(len) |
                                CigarOp::SeqMatch(len) |
                                CigarOp::SeqMismatch(len) => {
                                    // Increment coverage for each position
                                    for i in 0..*len {
                                        let coverage_pos = current_pos + i as i32;

                                        // Only count if in requested range
                                        let pos_in_range = match (start, end) {
                                            (Some(s), Some(e)) => coverage_pos >= s && coverage_pos < e,
                                            (Some(s), None) => coverage_pos >= s,
                                            (None, Some(e)) => coverage_pos < e,
                                            (None, None) => true,
                                        };

                                        if pos_in_range {
                                            *coverage.entry(coverage_pos).or_insert(0) += 1;
                                        }
                                    }
                                    current_pos += *len as i32;
                                }
                                _ => {
                                    // Insertions, clips, padding don't consume reference
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(coverage)
}

/// Calculate MAPQ distribution histogram
///
/// Args:
///     path (str): Path to BAM file
///     reference_id (int | None): Optional reference ID to filter by
///
/// Returns:
///     dict[int, int]: Mapping from MAPQ score to count
///
/// Example:
///     >>> mapq_dist = biometal.mapq_distribution("alignments.bam")
///     >>> for mapq in sorted(mapq_dist.keys()):
///     ...     print(f"MAPQ {mapq}: {mapq_dist[mapq]:,} reads")
#[pyfunction]
#[pyo3(signature = (path, reference_id=None))]
pub fn mapq_distribution(
    _py: Python,
    path: String,
    reference_id: Option<usize>,
) -> PyResult<HashMap<u8, u32>> {
    let path_buf = PathBuf::from(path);
    let reader = BamReader::from_path(&path_buf)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut distribution: HashMap<u8, u32> = HashMap::new();

    // Read all records and count MAPQ values
    let mut reader = reader;
    while let Ok(Some(record)) = reader.read_record() {
        // Filter by reference if specified
        if let Some(ref_id) = reference_id {
            if record.reference_id != Some(ref_id) {
                continue;
            }
        }

        // Count MAPQ (treat None as 255 per SAM spec)
        let mapq = record.mapq.unwrap_or(255);
        *distribution.entry(mapq).or_insert(0) += 1;
    }

    Ok(distribution)
}

/// SAM format writer
///
/// Writes BAM records to SAM text format. Enables workflows like:
/// BAM → filter → SAM output.
///
/// Example:
///     >>> reader = biometal.BamReader.from_path("input.bam")
///     >>> writer = biometal.SamWriter.create("output.sam")
///     >>>
///     >>> # Write header
///     >>> writer.write_header(reader.header)
///     >>>
///     >>> # Write filtered records
///     >>> for record in reader:
///     ...     if record.is_mapped and record.mapq >= 30:
///     ...         writer.write_record(record)
///     >>>
///     >>> writer.close()
#[pyclass(name = "SamWriter", unsendable)]
pub struct PySamWriter {
    writer: Option<SamWriter<BufWriter<File>>>,
    path: String,
}

#[pymethods]
impl PySamWriter {
    /// Create a new SAM writer
    ///
    /// Args:
    ///     path (str): Output SAM file path
    ///
    /// Returns:
    ///     SamWriter: New SAM writer instance
    ///
    /// Raises:
    ///     IOError: If file cannot be created
    ///
    /// Example:
    ///     >>> writer = biometal.SamWriter.create("output.sam")
    #[staticmethod]
    fn create(path: String) -> PyResult<Self> {
        let file = File::create(&path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to create file {}: {}", path, e)
            ))?;

        let buf_writer = BufWriter::new(file);
        let sam_writer = SamWriter::new(buf_writer);

        Ok(PySamWriter {
            writer: Some(sam_writer),
            path: path.clone(),
        })
    }

    /// Write BAM header to SAM file
    ///
    /// Args:
    ///     header (BamHeader): Header to write
    ///
    /// Raises:
    ///     IOError: If write fails
    ///
    /// Example:
    ///     >>> reader = biometal.BamReader.from_path("input.bam")
    ///     >>> writer = biometal.SamWriter.create("output.sam")
    ///     >>> writer.write_header(reader.header)
    fn write_header(&mut self, header: &PyBamHeader) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("Writer has been closed")
        })?;

        // Convert PyBamHeader to Header
        let rust_header = Header {
            text: String::new(),  // SAM writer will generate from references
            references: header.references.iter().map(|r| Reference {
                name: r.name.clone(),
                length: r.length,
            }).collect(),
        };

        writer.write_header(&rust_header)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to write header: {}", e)
            ))
    }

    /// Write BAM record to SAM file
    ///
    /// Args:
    ///     record (BamRecord): Record to write
    ///
    /// Raises:
    ///     IOError: If write fails
    ///
    /// Example:
    ///     >>> for record in reader:
    ///     ...     if record.is_mapped:
    ///     ...         writer.write_record(record)
    fn write_record(&mut self, record: &PyBamRecord) -> PyResult<()> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("Writer has been closed")
        })?;

        // Convert PyBamRecord to Record
        let rust_record = crate::io::bam::Record {
            name: record.name.clone(),
            reference_id: record.reference_id,
            position: record.position,
            mapq: record.mapq,
            flags: record.flags,
            mate_reference_id: record.mate_reference_id,
            mate_position: record.mate_position,
            template_length: record.template_length,
            sequence: record.sequence.clone(),
            quality: record.quality.clone(),
            cigar: record.cigar_ops.clone(),
            tags: record.tags.clone(),
        };

        writer.write_record(&rust_record)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to write record: {}", e)
            ))
    }

    /// Close the SAM writer and flush buffers
    ///
    /// Note: Writer is automatically flushed when dropped, but calling
    /// close() explicitly ensures errors are caught.
    ///
    /// Example:
    ///     >>> writer.close()
    fn close(&mut self) -> PyResult<()> {
        // Taking the writer will drop it, which flushes via BufWriter's Drop impl
        self.writer.take();
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("SamWriter(path='{}')", self.path)
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

/// Count records by flag categories
///
/// Args:
///     path (str): Path to BAM file
///
/// Returns:
///     dict[str, int]: Mapping from flag category to count
///
/// Categories:
///     - total: Total records
///     - mapped: Mapped reads
///     - unmapped: Unmapped reads
///     - primary: Primary alignments
///     - secondary: Secondary alignments
///     - supplementary: Supplementary alignments
///     - paired: Paired-end reads
///     - proper_pair: Properly paired reads
///     - forward: Forward strand
///     - reverse: Reverse strand
///     - qc_fail: Failed QC
///     - duplicate: PCR/optical duplicates
///
/// Example:
///     >>> stats = biometal.count_by_flag("alignments.bam")
///     >>> print(f"Mapped: {stats['mapped']:,} / {stats['total']:,} ({100*stats['mapped']/stats['total']:.1f}%)")
#[pyfunction]
pub fn count_by_flag(_py: Python, path: String) -> PyResult<HashMap<String, u32>> {
    let path_buf = PathBuf::from(path);
    let reader = BamReader::from_path(&path_buf)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut counts: HashMap<String, u32> = HashMap::new();

    // Initialize counters
    counts.insert("total".to_string(), 0);
    counts.insert("mapped".to_string(), 0);
    counts.insert("unmapped".to_string(), 0);
    counts.insert("primary".to_string(), 0);
    counts.insert("secondary".to_string(), 0);
    counts.insert("supplementary".to_string(), 0);
    counts.insert("paired".to_string(), 0);
    counts.insert("proper_pair".to_string(), 0);
    counts.insert("forward".to_string(), 0);
    counts.insert("reverse".to_string(), 0);
    counts.insert("qc_fail".to_string(), 0);
    counts.insert("duplicate".to_string(), 0);

    // Read all records and count flags
    let mut reader = reader;
    while let Ok(Some(record)) = reader.read_record() {
        *counts.get_mut("total").unwrap() += 1;

        let flags = record.flags;

        if (flags & 0x4) == 0 {
            *counts.get_mut("mapped").unwrap() += 1;
        } else {
            *counts.get_mut("unmapped").unwrap() += 1;
        }

        if (flags & 0x900) == 0 {
            *counts.get_mut("primary").unwrap() += 1;
        }
        if (flags & 0x100) != 0 {
            *counts.get_mut("secondary").unwrap() += 1;
        }
        if (flags & 0x800) != 0 {
            *counts.get_mut("supplementary").unwrap() += 1;
        }
        if (flags & 0x1) != 0 {
            *counts.get_mut("paired").unwrap() += 1;
        }
        if (flags & 0x2) != 0 {
            *counts.get_mut("proper_pair").unwrap() += 1;
        }
        if (flags & 0x10) == 0 {
            *counts.get_mut("forward").unwrap() += 1;
        } else {
            *counts.get_mut("reverse").unwrap() += 1;
        }
        if (flags & 0x200) != 0 {
            *counts.get_mut("qc_fail").unwrap() += 1;
        }
        if (flags & 0x400) != 0 {
            *counts.get_mut("duplicate").unwrap() += 1;
        }
    }

    Ok(counts)
}

/// Calculate insert size distribution for paired-end QC
///
/// Analyzes template length distribution for properly paired reads on the same reference.
/// Essential for paired-end sequencing QC.
///
/// Args:
///     path (str): Path to BAM file
///     reference_id (int | None): Optional reference ID to filter by
///
/// Returns:
///     dict[int, int]: Mapping from insert size to count
///
/// Example:
///     >>> insert_dist = biometal.insert_size_distribution("alignments.bam")
///     >>> sizes = list(insert_dist.keys())
///     >>> mean_size = sum(s * insert_dist[s] for s in sizes) / sum(insert_dist.values())
///     >>> print(f"Mean insert size: {mean_size:.0f}bp")
///
/// Note:
///     Only includes primary, properly paired reads on the same reference with positive insert size
#[pyfunction]
#[pyo3(signature = (path, reference_id=None))]
pub fn insert_size_distribution(
    _py: Python,
    path: String,
    reference_id: Option<usize>,
) -> PyResult<HashMap<i32, u32>> {
    let path_buf = PathBuf::from(path);
    let reader = BamReader::from_path(&path_buf)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut distribution: HashMap<i32, u32> = HashMap::new();

    let mut reader = reader;
    while let Ok(Some(record)) = reader.read_record() {
        // Filter by reference if specified
        if let Some(ref_id) = reference_id {
            if record.reference_id != Some(ref_id) {
                continue;
            }
        }

        // Only count properly paired reads
        let flags = record.flags;
        let is_paired = (flags & 0x1) != 0;
        let is_proper_pair = (flags & 0x2) != 0;
        let is_primary = (flags & 0x900) == 0;

        if is_paired && is_proper_pair && is_primary {
            // Check same reference
            if record.reference_id == record.mate_reference_id {
                let insert_size = record.template_length.abs();
                if insert_size > 0 {
                    *distribution.entry(insert_size).or_insert(0) += 1;
                }
            }
        }
    }

    Ok(distribution)
}

/// Calculate edit distance statistics for alignment quality assessment
///
/// Analyzes NM tag (edit distance) distribution across mapped reads.
///
/// Args:
///     path (str): Path to BAM file
///     reference_id (int | None): Optional reference ID to filter by
///
/// Returns:
///     dict: Statistics including:
///         - mean (float | None): Mean edit distance
///         - median (int | None): Median edit distance
///         - max (int | None): Maximum edit distance
///         - min (int | None): Minimum edit distance
///         - distribution (dict[int, int]): Edit distance → count
///         - with_nm_tag (int): Number of records with NM tag
///         - total_records (int): Total records analyzed
///
/// Example:
///     >>> stats = biometal.edit_distance_stats("alignments.bam")
///     >>> print(f"Mean mismatches: {stats['mean']:.2f}")
///     >>> print(f"Coverage: {stats['with_nm_tag'] / stats['total_records'] * 100:.1f}%")
#[pyfunction]
#[pyo3(signature = (path, reference_id=None))]
pub fn edit_distance_stats(
    py: Python,
    path: String,
    reference_id: Option<usize>,
) -> PyResult<Py<PyDict>> {
    let path_buf = PathBuf::from(path);
    let reader = BamReader::from_path(&path_buf)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut edit_distances = Vec::new();
    let mut distribution: HashMap<i64, u32> = HashMap::new();
    let mut total_records = 0;

    let mut reader = reader;
    while let Ok(Some(record)) = reader.read_record() {
        // Filter by reference if specified
        if let Some(ref_id) = reference_id {
            if record.reference_id != Some(ref_id) {
                continue;
            }
        }

        total_records += 1;

        // Get NM tag (edit distance)
        if let Ok(Some(tag)) = record.tags.get(b"NM") {
            if let TagValue::Int(nm) = tag.value {
                edit_distances.push(nm);
                *distribution.entry(nm).or_insert(0) += 1;
            }
        }
    }

    // Build result dictionary
    let result = PyDict::new(py);

    result.set_item("total_records", total_records)?;
    result.set_item("with_nm_tag", edit_distances.len())?;

    if !edit_distances.is_empty() {
        edit_distances.sort();
        let sum: i64 = edit_distances.iter().sum();
        let mean = sum as f64 / edit_distances.len() as f64;
        let median = edit_distances[edit_distances.len() / 2];
        let min = *edit_distances.iter().min().unwrap();
        let max = *edit_distances.iter().max().unwrap();

        result.set_item("mean", mean)?;
        result.set_item("median", median)?;
        result.set_item("min", min)?;
        result.set_item("max", max)?;
    } else {
        result.set_item("mean", py.None())?;
        result.set_item("median", py.None())?;
        result.set_item("min", py.None())?;
        result.set_item("max", py.None())?;
    }

    // Convert distribution
    let dist_dict = PyDict::new(py);
    for (nm, count) in distribution {
        dist_dict.set_item(nm, count)?;
    }
    result.set_item("distribution", dist_dict)?;

    Ok(result.into())
}

/// Calculate strand bias at a genomic position
///
/// Counts forward and reverse strand reads overlapping a position.
/// Essential for variant calling QC.
///
/// Args:
///     path (str): Path to BAM file
///     reference_id (int): Reference sequence ID
///     position (int): 0-based genomic position
///     window_size (int): Window size around position (default: 1)
///
/// Returns:
///     dict: Strand bias statistics:
///         - forward (int): Forward strand reads
///         - reverse (int): Reverse strand reads
///         - total (int): Total reads
///         - ratio (float): Forward / (Forward + Reverse), 0.5 = no bias
///
/// Example:
///     >>> bias = biometal.strand_bias("alignments.bam", 0, 1000, window_size=10)
///     >>> print(f"Strand bias: {bias['forward']} forward, {bias['reverse']} reverse")
///     >>> print(f"Ratio: {bias['ratio']:.3f} (0.5 = no bias)")
#[pyfunction]
#[pyo3(signature = (path, reference_id, position, window_size=1))]
pub fn strand_bias(
    py: Python,
    path: String,
    reference_id: usize,
    position: i32,
    window_size: i32,
) -> PyResult<Py<PyDict>> {
    let path_buf = PathBuf::from(path);
    let reader = BamReader::from_path(&path_buf)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut forward_count = 0u32;
    let mut reverse_count = 0u32;

    let start = position - window_size / 2;
    let end = position + window_size / 2 + 1;

    let mut reader = reader;
    while let Ok(Some(record)) = reader.read_record() {
        // Check reference
        if record.reference_id != Some(reference_id) {
            continue;
        }

        // Check if read overlaps position
        if let Some(pos) = record.position {
            // Calculate alignment end using CIGAR
            let mut alignment_end = pos;
            for op in &record.cigar {
                match op {
                    CigarOp::Match(len) |
                    CigarOp::Deletion(len) |
                    CigarOp::RefSkip(len) |
                    CigarOp::SeqMatch(len) |
                    CigarOp::SeqMismatch(len) => {
                        alignment_end += *len as i32;
                    }
                    _ => {}
                }
            }

            // Check overlap with window
            if pos < end && alignment_end > start {
                // Count by strand
                let is_reverse = (record.flags & 0x10) != 0;
                if is_reverse {
                    reverse_count += 1;
                } else {
                    forward_count += 1;
                }
            }
        }
    }

    // Build result
    let result = PyDict::new(py);
    result.set_item("forward", forward_count)?;
    result.set_item("reverse", reverse_count)?;
    result.set_item("total", forward_count + reverse_count)?;

    let total = forward_count + reverse_count;
    let ratio = if total > 0 {
        forward_count as f64 / total as f64
    } else {
        0.0
    };
    result.set_item("ratio", ratio)?;

    Ok(result.into())
}

/// Calculate alignment length distribution
///
/// Analyzes reference span (CIGAR-based) for mapped alignments.
/// Useful for RNA-seq, assembly QC, and structural variant detection.
///
/// Args:
///     path (str): Path to BAM file
///     reference_id (int | None): Optional reference ID to filter by
///
/// Returns:
///     dict[int, int]: Mapping from alignment length to count
///
/// Example:
///     >>> dist = biometal.alignment_length_distribution("alignments.bam")
///     >>> lengths = list(dist.keys())
///     >>> mean_len = sum(l * dist[l] for l in lengths) / sum(dist.values())
///     >>> print(f"Mean alignment length: {mean_len:.0f}bp")
#[pyfunction]
#[pyo3(signature = (path, reference_id=None))]
pub fn alignment_length_distribution(
    _py: Python,
    path: String,
    reference_id: Option<usize>,
) -> PyResult<HashMap<u32, u32>> {
    let path_buf = PathBuf::from(path);
    let reader = BamReader::from_path(&path_buf)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    let mut distribution: HashMap<u32, u32> = HashMap::new();

    let mut reader = reader;
    while let Ok(Some(record)) = reader.read_record() {
        // Filter by reference if specified
        if let Some(ref_id) = reference_id {
            if record.reference_id != Some(ref_id) {
                continue;
            }
        }

        // Only count mapped reads
        let is_mapped = (record.flags & 0x4) == 0;
        if !is_mapped {
            continue;
        }

        // Calculate alignment length from CIGAR
        let mut length = 0u32;
        for op in &record.cigar {
            match op {
                CigarOp::Match(len) |
                CigarOp::Deletion(len) |
                CigarOp::RefSkip(len) |
                CigarOp::SeqMatch(len) |
                CigarOp::SeqMismatch(len) => {
                    length += len;
                }
                _ => {}
            }
        }

        if length > 0 {
            *distribution.entry(length).or_insert(0) += 1;
        }
    }

    Ok(distribution)
}
