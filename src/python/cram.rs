//! Python bindings for CRAM format (reference-based compressed alignment)

use pyo3::prelude::*;
use pyo3::exceptions::PyStopIteration;
use crate::io::cram::{CramReader, CramRecords};
use crate::python::bam::{PyBamRecord, PyBamHeader};
use std::path::PathBuf;
use std::fs::File;
use std::io::BufReader;

/// Read CRAM files with ARM NEON optimizations
///
/// Native ARM-optimized CRAM 3.0/3.1 parser with streaming architecture.
/// CRAM files are 30-60% smaller than BAM through reference-based compression.
///
/// Args:
    ///     path (str): Path to CRAM file (.cram)
///
/// Features:
///     - ARM NEON optimized (9× base counting speedup)
///     - Streaming architecture (constant ~5 MB memory)
///     - Multi-codec support (gzip, bzip2, lzma, rANS 4x16)
///     - BAM-compatible output (same Record format)
///     - Optional reference FASTA for sequence reconstruction
///
/// Example:
///     >>> import biometal
///     >>> cram = biometal.CramReader.from_path("alignments.cram")
///     >>>
///     >>> # Access header info (BAM-compatible)
///     >>> print(f"CRAM version: {cram.major_version}.{cram.minor_version}")
///     >>> print(f"References: {cram.header.reference_count}")
///     >>>
///     >>> # Stream records with constant memory
///     >>> for record in cram:
///     ...     if record.is_mapped and record.mapq and record.mapq >= 30:
///     ...         print(f"{record.name}: {record.position}")
///     ...         # Process high-quality alignments
///     >>>
///     >>> # With reference FASTA for full sequence reconstruction
///     >>> cram = biometal.CramReader.from_path_with_reference(
///     ...     "alignments.cram",
///     ...     "reference.fa"
///     ... )
///     >>> for record in cram:
///     ...     print(f"{record.name}: {len(record.sequence)} bp")
///
/// Note:
///     - Returns same PyBamRecord format as BamReader for compatibility
///     - Memory footprint remains constant at ~5 MB
///     - Unique ARM-native CRAM implementation
#[pyclass(name = "CramReader", unsendable)]
pub struct PyCramReader {
    /// Iterator over CRAM records
    inner: Option<CramRecords<BufReader<File>>>,
    /// BAM-compatible header (for Python access)
    header: PyBamHeader,
    /// CRAM major version
    major_version: u8,
    /// CRAM minor version
    minor_version: u8,
}

#[pymethods]
impl PyCramReader {
    /// Open CRAM file for streaming
    ///
    /// Args:
    ///     path (str): Path to CRAM file
    ///
    /// Returns:
    ///     CramReader: Streaming iterator with ~5 MB constant memory
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    ///
    /// Example:
    ///     >>> cram = biometal.CramReader.from_path("sample.cram")
    ///     >>> print(f"Opened CRAM {cram.major_version}.{cram.minor_version}")
    ///     >>> print(f"References: {cram.header.reference_count}")
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = CramReader::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        // Extract metadata before consuming reader into iterator
        let header = PyBamHeader::from(reader.header());
        let major_version = reader.major_version();
        let minor_version = reader.minor_version();

        // Convert to iterator (consumes the reader)
        let records = reader.records();

        Ok(PyCramReader {
            inner: Some(records),
            header,
            major_version,
            minor_version,
        })
    }

    /// Open CRAM file with reference FASTA for sequence reconstruction
    ///
    /// Args:
    ///     cram_path (str): Path to CRAM file
    ///     reference_path (str): Path to reference FASTA file (e.g., "hg38.fa")
    ///
    /// Returns:
    ///     CramReader: Streaming iterator with full sequence reconstruction
    ///
    /// Raises:
    ///     IOError: If files cannot be opened
    ///     ValueError: If file formats are invalid
    ///
    /// Note:
    ///     Requires reference FASTA index file (.fai) in same directory.
    ///     Use `samtools faidx reference.fa` to create if missing.
    ///
    /// Example:
    ///     >>> cram = biometal.CramReader.from_path_with_reference(
    ///     ...     "alignments.cram",
    ///     ...     "hg38.fa"
    ///     ... )
    ///     >>> for record in cram:
    ///     ...     print(f"{record.name}: {len(record.sequence)} bp")
    #[staticmethod]
    fn from_path_with_reference(cram_path: String, reference_path: String) -> PyResult<Self> {
        let cram_path_buf = PathBuf::from(cram_path);
        let reference_path_buf = PathBuf::from(reference_path);

        let reader = CramReader::from_path_with_reference(&cram_path_buf, &reference_path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        // Extract metadata before consuming reader into iterator
        let header = PyBamHeader::from(reader.header());
        let major_version = reader.major_version();
        let minor_version = reader.minor_version();

        // Convert to iterator (consumes the reader)
        let records = reader.records();

        Ok(PyCramReader {
            inner: Some(records),
            header,
            major_version,
            minor_version,
        })
    }

    /// Get CRAM header (BAM-compatible format)
    ///
    /// Returns:
    ///     BamHeader: Header with reference information
    ///
    /// Example:
    ///     >>> header = cram.header
    ///     >>> print(f"References: {', '.join(header.reference_names)}")
    #[getter]
    fn header(&self) -> PyBamHeader {
        self.header.clone()
    }

    /// Get CRAM major version
    ///
    /// Returns:
    ///     int: CRAM major version (2 or 3)
    ///
    /// Example:
    ///     >>> print(f"CRAM version: {cram.major_version}.{cram.minor_version}")
    #[getter]
    fn major_version(&self) -> u8 {
        self.major_version
    }

    /// Get CRAM minor version
    ///
    /// Returns:
    ///     int: CRAM minor version
    ///
    /// Example:
    ///     >>> print(f"CRAM version: {cram.major_version}.{cram.minor_version}")
    #[getter]
    fn minor_version(&self) -> u8 {
        self.minor_version
    }

    /// Get number of reference sequences
    ///
    /// Returns:
    ///     int: Number of reference sequences in header
    ///
    /// Example:
    ///     >>> print(f"References: {cram.reference_count}")
    #[getter]
    fn reference_count(&self) -> usize {
        self.header.reference_count
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyBamRecord> {
        if let Some(ref mut records) = slf.inner {
            match records.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(PyStopIteration::new_err("iterator exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "CramReader(version={}.{}, references={})",
            self.major_version, self.minor_version, self.header.reference_count
        )
    }
}
