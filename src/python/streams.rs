//! Python wrappers for streaming classes

use pyo3::prelude::*;
use pyo3::exceptions::PyStopIteration;
use std::path::PathBuf;
use crate::io::{FastqStream, FastaStream, CompressedReader};
use crate::python::records::{PyFastqRecord, PyFastaRecord};

/// Stream FASTQ records with constant memory
///
/// Streaming iterator that processes FASTQ files one record at a time,
/// maintaining constant ~5 MB memory regardless of file size.
///
/// Args:
///     path (str): Path to FASTQ file (.fq, .fastq, .fq.gz, .fastq.gz)
///
/// Example:
///     >>> stream = biometal.FastqStream.from_path("data.fq.gz")
///     >>> for record in stream:
///     ...     gc = biometal.gc_content(record.sequence)
///     ...     print(f"{record.id}: {gc:.2%}")
///
/// Note:
///     Memory footprint remains constant at ~5 MB even for TB-scale files.
///     This enables analysis on consumer hardware without downloading.
#[pyclass(name = "FastqStream", unsendable)]
pub struct PyFastqStream {
    inner: Option<FastqStream<CompressedReader>>,
}

#[pymethods]
impl PyFastqStream {
    /// Create FASTQ stream from file path
    ///
    /// Args:
///         path (str): Path to FASTQ file
    ///
    /// Returns:
    ///     FastqStream: Streaming iterator
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let stream = FastqStream::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastqStream {
            inner: Some(stream),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyFastqRecord> {
        if let Some(ref mut stream) = slf.inner {
            match stream.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(PyStopIteration::new_err("stream exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "FastqStream(...)".to_string()
    }
}

/// Stream FASTA records with constant memory
///
/// Streaming iterator that processes FASTA files one record at a time,
/// maintaining constant ~5 MB memory regardless of file size.
///
/// Args:
///     path (str): Path to FASTA file (.fa, .fasta, .fa.gz, .fasta.gz)
///
/// Example:
///     >>> stream = biometal.FastaStream.from_path("genome.fa.gz")
///     >>> for record in stream:
///     ...     length = len(record.sequence)
///     ...     print(f"{record.id}: {length} bp")
///
/// Note:
///     Memory footprint remains constant at ~5 MB even for entire genomes.
#[pyclass(name = "FastaStream", unsendable)]
pub struct PyFastaStream {
    inner: Option<FastaStream<CompressedReader>>,
}

#[pymethods]
impl PyFastaStream {
    /// Create FASTA stream from file path
    ///
    /// Args:
    ///     path (str): Path to FASTA file
    ///
    /// Returns:
    ///     FastaStream: Streaming iterator
    ///
    /// Raises:
    ///     IOError: If file cannot be opened
    ///     ValueError: If file format is invalid
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let stream = FastaStream::from_path(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

        Ok(PyFastaStream {
            inner: Some(stream),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyFastaRecord> {
        if let Some(ref mut stream) = slf.inner {
            match stream.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(PyStopIteration::new_err("stream exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "FastaStream(...)".to_string()
    }
}
