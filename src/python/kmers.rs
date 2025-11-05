//! K-mer extraction utilities for machine learning preprocessing

use pyo3::prelude::*;

/// Extract overlapping k-mers from sequence
///
/// Extracts all k-mers of specified length using a sliding window.
/// Useful for BERT/transformer model preprocessing.
///
/// Args:
///     sequence (bytes): DNA/RNA sequence
///     k (int): K-mer length (typically 3-8 for BERT models)
///
/// Returns:
///     list[str]: List of k-mers as strings
///
/// Example:
///     >>> kmers = biometal.extract_kmers(b"ATGCAT", k=3)
///     >>> print(kmers)
///     ['ATG', 'TGC', 'GCA', 'CAT']
///
/// Note:
///     For BERT preprocessing, typical k values are:
///     - DNABert: k=3, k=4, k=5, k=6
///     - Enformer: k=6
///     - Nucleotide Transformer: k=6
///
/// Performance:
///     This operation is memory-efficient and processes sequences
///     in a single pass. For large datasets, combine with streaming:
///
///     >>> stream = biometal.FastqStream.from_path("data.fq.gz")
///     >>> for record in stream:
///     ...     kmers = biometal.extract_kmers(record.sequence, k=6)
///     ...     # Feed to BERT model
#[pyfunction(name = "extract_kmers")]
pub fn py_extract_kmers(sequence: &[u8], k: usize) -> PyResult<Vec<String>> {
    if k == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "k must be greater than 0"
        ));
    }

    if k > sequence.len() {
        return Ok(Vec::new());
    }

    let mut kmers = Vec::with_capacity(sequence.len().saturating_sub(k - 1));

    for window in sequence.windows(k) {
        kmers.push(String::from_utf8_lossy(window).to_string());
    }

    Ok(kmers)
}

/// Extract non-overlapping k-mers from sequence
///
/// Extracts k-mers without overlap, stepping by k each time.
/// Useful for k-mer counting and composition analysis.
///
/// Args:
///     sequence (bytes): DNA/RNA sequence
///     k (int): K-mer length
///
/// Returns:
///     list[str]: List of non-overlapping k-mers
///
/// Example:
///     >>> kmers = biometal.extract_kmers_non_overlapping(b"ATGCATGC", k=4)
///     >>> print(kmers)
///     ['ATGC', 'ATGC']
///
/// Note:
///     Partial k-mers at the end are discarded.
#[pyfunction(name = "extract_kmers_non_overlapping")]
pub fn py_extract_kmers_non_overlapping(sequence: &[u8], k: usize) -> PyResult<Vec<String>> {
    if k == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "k must be greater than 0"
        ));
    }

    let mut kmers = Vec::with_capacity(sequence.len() / k);

    for chunk in sequence.chunks_exact(k) {
        kmers.push(String::from_utf8_lossy(chunk).to_string());
    }

    Ok(kmers)
}
