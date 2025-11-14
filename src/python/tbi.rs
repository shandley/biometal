//! Python bindings for TBI (Tabix Index)

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::formats::index::TbiIndex;
use crate::io::bam::index::VirtualOffset;

/// Python wrapper for TBI index
///
/// This class provides access to Tabix index functionality, enabling
/// fast random access to sorted, tab-delimited, BGZF-compressed files.
///
/// Example:
///     >>> import biometal
///     >>> # Load tabix index
///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
///     >>> print(f"Index has {len(index.references())} references")
///     >>>
///     >>> # Query region
///     >>> chunks = index.query("chr1", 1000000, 2000000)
///     >>> print(f"Found {len(chunks)} chunks to read")
///     >>>
///     >>> # Get metadata
///     >>> print(f"Format: {index.format()}")
///     >>> print(f"Sequence column: {index.col_seq()}")
#[pyclass(name = "TbiIndex")]
pub struct PyTbiIndex {
    inner: TbiIndex,
}

#[pymethods]
impl PyTbiIndex {
    /// Load TBI index from a file
    ///
    /// Args:
    ///     tbi_path: Path to .tbi index file
    ///
    /// Returns:
    ///     TbiIndex: Loaded index
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"Loaded {len(index.references())} references")
    #[staticmethod]
    fn from_path(tbi_path: String) -> PyResult<Self> {
        let index = TbiIndex::from_path(&tbi_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyTbiIndex { inner: index })
    }

    /// Query region and get chunks to read
    ///
    /// Returns list of (start_offset, end_offset) tuples representing
    /// file chunks that overlap the query region.
    ///
    /// Args:
    ///     ref_name: Reference sequence name (e.g., "chr1")
    ///     start: Start position (0-based, inclusive)
    ///     end: End position (0-based, exclusive)
    ///
    /// Returns:
    ///     list[tuple[int, int]]: List of (start, end) virtual file offset pairs
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> chunks = index.query("chr1", 1000000, 2000000)
    ///     >>> for start, end in chunks:
    ///     ...     print(f"Read chunk: {start:016x} - {end:016x}")
    fn query(&self, ref_name: String, start: u32, end: u32) -> PyResult<Vec<(u64, u64)>> {
        let chunks = self
            .inner
            .query(&ref_name, start, end)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;

        Ok(chunks
            .iter()
            .map(|chunk| (chunk.start.as_raw(), chunk.end.as_raw()))
            .collect())
    }

    /// Get file format type
    ///
    /// Returns:
    ///     str: Format type ("Generic", "SAM", or "VCF")
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"Format: {index.format()}")
    fn format(&self) -> String {
        format!("{:?}", self.inner.format())
    }

    /// Get column for sequence name (0-based)
    ///
    /// Returns:
    ///     int: Column index for sequence/chromosome name
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"Sequence column: {index.col_seq()}")
    fn col_seq(&self) -> i32 {
        self.inner.col_seq()
    }

    /// Get column for start position (0-based)
    ///
    /// Returns:
    ///     int: Column index for start position
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"Start column: {index.col_beg()}")
    fn col_beg(&self) -> i32 {
        self.inner.col_beg()
    }

    /// Get column for end position (0-based)
    ///
    /// Returns:
    ///     int: Column index for end position (0 if same as start)
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"End column: {index.col_end()}")
    fn col_end(&self) -> i32 {
        self.inner.col_end()
    }

    /// Get comment character for header lines
    ///
    /// Returns:
    ///     str: Comment character (usually '#')
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"Comment char: '{index.meta_char()}'")
    fn meta_char(&self) -> String {
        self.inner.meta_char().to_string()
    }

    /// Get number of header lines to skip
    ///
    /// Returns:
    ///     int: Number of lines to skip at start of file
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"Skip {index.skip_lines()} header lines")
    fn skip_lines(&self) -> i32 {
        self.inner.skip_lines()
    }

    /// Get list of reference sequence names
    ///
    /// Returns:
    ///     list[str]: Reference sequence names in order
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> for ref_name in index.references():
    ///     ...     print(ref_name)
    fn references(&self) -> Vec<String> {
        self.inner
            .references()
            .iter()
            .map(|r| r.name.clone())
            .collect()
    }

    /// Check if a reference exists in the index
    ///
    /// Args:
    ///     ref_name: Reference sequence name
    ///
    /// Returns:
    ///     bool: True if reference exists
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> if index.contains("chr1"):
    ///     ...     print("chr1 found")
    fn contains(&self, ref_name: String) -> bool {
        self.inner.get_reference(&ref_name).is_some()
    }

    /// Get information about a specific reference
    ///
    /// Args:
    ///     ref_name: Reference sequence name
    ///
    /// Returns:
    ///     dict: Dictionary with 'name', 'n_bins', 'n_intervals'
    ///           Returns None if reference not found
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> info = index.get_info("chr1")
    ///     >>> print(f"chr1: {info['n_bins']} bins, {info['n_intervals']} intervals")
    fn get_info(&self, py: Python, ref_name: String) -> Option<Py<PyDict>> {
        self.inner.get_reference(&ref_name).map(|reference| {
            let dict = PyDict::new(py);
            dict.set_item("name", &reference.name).unwrap();
            dict.set_item("n_bins", reference.bins.len()).unwrap();
            dict.set_item("n_intervals", reference.intervals.len())
                .unwrap();
            dict.into()
        })
    }

    /// Get number of references in index
    ///
    /// Returns:
    ///     int: Number of reference sequences
    ///
    /// Example:
    ///     >>> index = biometal.TbiIndex.from_path("data.vcf.gz.tbi")
    ///     >>> print(f"Index has {len(index)} references")
    fn __len__(&self) -> usize {
        self.inner.references().len()
    }

    /// String representation
    fn __repr__(&self) -> String {
        format!(
            "TbiIndex(format={:?}, references={})",
            self.inner.format(),
            self.inner.references().len()
        )
    }
}
