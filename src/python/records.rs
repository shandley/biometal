//! Python wrappers for FASTQ and FASTA records

use pyo3::prelude::*;
use crate::types::{FastqRecord, FastaRecord};

/// FASTQ record with sequence, quality, and metadata
///
/// Attributes:
///     id (str): Read identifier
///     sequence (bytes): DNA/RNA sequence
///     quality (bytes): Phred quality scores
///
/// Example:
///     >>> record = stream.next()
///     >>> print(record.id)
///     >>> gc = biometal.gc_content(record.sequence)
#[pyclass(name = "FastqRecord")]
#[derive(Clone)]
pub struct PyFastqRecord {
    /// Read identifier (e.g., "@SRR390728.1")
    #[pyo3(get)]
    pub id: String,
    /// DNA/RNA sequence as bytes
    #[pyo3(get)]
    pub sequence: Vec<u8>,
    /// Phred quality scores as bytes (ASCII-encoded)
    #[pyo3(get)]
    pub quality: Vec<u8>,
}

#[pymethods]
impl PyFastqRecord {
    fn __repr__(&self) -> String {
        format!(
            "FastqRecord(id='{}', seq_len={}, qual_len={})",
            self.id,
            self.sequence.len(),
            self.quality.len()
        )
    }

    fn __str__(&self) -> String {
        format!("@{}\n{}\n+\n{}",
            self.id,
            String::from_utf8_lossy(&self.sequence),
            String::from_utf8_lossy(&self.quality)
        )
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
}

impl From<FastqRecord> for PyFastqRecord {
    fn from(record: FastqRecord) -> Self {
        PyFastqRecord {
            id: record.id,
            sequence: record.sequence,
            quality: record.quality,
        }
    }
}

/// FASTA record with sequence and metadata
///
/// Attributes:
///     id (str): Sequence identifier
///     sequence (bytes): DNA/RNA/protein sequence
///
/// Example:
///     >>> record = stream.next()
///     >>> print(record.id)
///     >>> gc = biometal.gc_content(record.sequence)
#[pyclass(name = "FastaRecord")]
#[derive(Clone)]
pub struct PyFastaRecord {
    /// Sequence identifier (e.g., ">chr1")
    #[pyo3(get)]
    pub id: String,
    /// DNA/RNA/protein sequence as bytes
    #[pyo3(get)]
    pub sequence: Vec<u8>,
}

#[pymethods]
impl PyFastaRecord {
    fn __repr__(&self) -> String {
        format!(
            "FastaRecord(id='{}', seq_len={})",
            self.id,
            self.sequence.len()
        )
    }

    fn __str__(&self) -> String {
        format!(">{}\n{}",
            self.id,
            String::from_utf8_lossy(&self.sequence)
        )
    }

    /// Get sequence as string
    #[getter]
    fn sequence_str(&self) -> String {
        String::from_utf8_lossy(&self.sequence).to_string()
    }
}

impl From<FastaRecord> for PyFastaRecord {
    fn from(record: FastaRecord) -> Self {
        PyFastaRecord {
            id: record.id,
            sequence: record.sequence,
        }
    }
}
