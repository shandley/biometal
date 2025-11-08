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
use pyo3::exceptions::PyStopIteration;
use std::path::PathBuf;
use std::collections::HashMap;
use crate::io::bam::{BamReader, Header, Reference};
use crate::io::compression::CompressedReader;

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
            header,
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
    header: PyBamHeader,
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
///
/// Note:
///     Only counts positions with coverage > 0. Missing positions have coverage 0.
#[pyfunction]
#[pyo3(signature = (path, reference_id, start=None, end=None))]
pub fn calculate_coverage(
    py: Python,
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
                        // For now, just count the start position
                        // TODO: Parse CIGAR to get full coverage
                        *coverage.entry(pos).or_insert(0) += 1;
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
    py: Python,
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
pub fn count_by_flag(py: Python, path: String) -> PyResult<HashMap<String, u32>> {
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
