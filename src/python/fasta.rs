//! Python bindings for FASTA format and FAI index

use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;

use crate::io::fasta::{FaiIndex, FastaStream};

/// Python wrapper for FAI index
///
/// This class provides access to FASTA index functionality, enabling
/// fast random access to sequences in FASTA files.
///
/// Example:
///     >>> import biometal
///     >>> # Build index from FASTA file
///     >>> index = biometal.FaiIndex.build("genome.fa")
///     >>> index.write("genome.fa.fai")
///     >>>
///     >>> # Load existing index
///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
///     >>> print(f"Index contains {len(index)} sequences")
///     >>>
///     >>> # Fetch entire sequence
///     >>> chr1 = index.fetch("chr1", "genome.fa")
///     >>> print(f"chr1: {len(chr1)} bp")
///     >>>
///     >>> # Fetch region (0-based, half-open)
///     >>> region = index.fetch_region("chr1", 1000, 2000, "genome.fa")
///     >>> print(f"chr1:1000-2000: {region}")
#[pyclass(name = "FaiIndex")]
pub struct PyFaiIndex {
    inner: FaiIndex,
}

#[pymethods]
impl PyFaiIndex {
    /// Build FAI index from a FASTA file
    ///
    /// Args:
    ///     fasta_path: Path to FASTA file
    ///
    /// Returns:
    ///     FaiIndex: Built index
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.build("genome.fa")
    ///     >>> index.write("genome.fa.fai")
    #[staticmethod]
    fn build(fasta_path: String) -> PyResult<Self> {
        let index = FaiIndex::build(&fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyFaiIndex { inner: index })
    }

    /// Load FAI index from a file
    ///
    /// Args:
    ///     fai_path: Path to .fai index file
    ///
    /// Returns:
    ///     FaiIndex: Loaded index
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
    ///     >>> print(f"Loaded {len(index)} sequences")
    #[staticmethod]
    fn from_path(fai_path: String) -> PyResult<Self> {
        let index = FaiIndex::from_path(&fai_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyFaiIndex { inner: index })
    }

    /// Write index to a file
    ///
    /// Args:
    ///     fai_path: Path to write .fai file
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.build("genome.fa")
    ///     >>> index.write("genome.fa.fai")
    fn write(&self, fai_path: String) -> PyResult<()> {
        self.inner
            .write(&fai_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
    }

    /// Fetch entire sequence by name
    ///
    /// Args:
    ///     name: Sequence name
    ///     fasta_path: Path to FASTA file
    ///
    /// Returns:
    ///     str: Full sequence
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
    ///     >>> chr1 = index.fetch("chr1", "genome.fa")
    ///     >>> print(f"chr1: {len(chr1)} bp")
    fn fetch(&self, name: String, fasta_path: String) -> PyResult<String> {
        self.inner
            .fetch(&name, &fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    }

    /// Fetch a region of a sequence (0-based, half-open interval [start, end))
    ///
    /// Args:
    ///     name: Sequence name
    ///     start: Start position (0-based, inclusive)
    ///     end: End position (0-based, exclusive)
    ///     fasta_path: Path to FASTA file
    ///
    /// Returns:
    ///     str: Sequence region
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
    ///     >>> region = index.fetch_region("chr1", 1000, 2000, "genome.fa")
    ///     >>> print(f"Region: {region}")  # 1000 bases
    fn fetch_region(
        &self,
        name: String,
        start: u64,
        end: u64,
        fasta_path: String,
    ) -> PyResult<String> {
        self.inner
            .fetch_region(&name, start, end, &fasta_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    }

    /// Get number of sequences in index
    ///
    /// Returns:
    ///     int: Number of sequences
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
    ///     >>> print(f"Index contains {len(index)} sequences")
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Get sequence names in index
    ///
    /// Returns:
    ///     list[str]: List of sequence names in order
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
    ///     >>> for name in index.sequence_names():
    ///     ...     print(name)
    fn sequence_names(&self) -> Vec<String> {
        self.inner.sequence_names.clone()
    }

    /// Get information about a specific sequence
    ///
    /// Args:
    ///     name: Sequence name
    ///
    /// Returns:
    ///     dict: Dictionary with 'length', 'offset', 'line_bases', 'line_width'
    ///           Returns None if sequence not found
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
    ///     >>> info = index.get_info("chr1")
    ///     >>> print(f"chr1 length: {info['length']} bp")
    fn get_info(&self, py: Python, name: String) -> Option<Py<PyDict>> {
        self.inner.get(&name).map(|entry| {
            let dict = PyDict::new(py);
            dict.set_item("name", &entry.name).unwrap();
            dict.set_item("length", entry.length).unwrap();
            dict.set_item("offset", entry.offset).unwrap();
            dict.set_item("line_bases", entry.line_bases).unwrap();
            dict.set_item("line_width", entry.line_width).unwrap();
            if let Some(qual_offset) = entry.qual_offset {
                dict.set_item("qual_offset", qual_offset).unwrap();
            }
            dict.into()
        })
    }

    /// Check if a sequence exists in the index
    ///
    /// Args:
    ///     name: Sequence name
    ///
    /// Returns:
    ///     bool: True if sequence exists
    ///
    /// Example:
    ///     >>> index = biometal.FaiIndex.from_path("genome.fa.fai")
    ///     >>> if index.contains("chr1"):
    ///     ...     print("chr1 found")
    fn contains(&self, name: String) -> bool {
        self.inner.get(&name).is_some()
    }

    /// String representation
    fn __repr__(&self) -> String {
        format!("FaiIndex(sequences={})", self.inner.len())
    }
}
