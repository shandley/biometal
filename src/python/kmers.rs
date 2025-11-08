//! Python bindings for k-mer operations
//!
//! Provides Python access to biometal's evidence-based k-mer operations.
//! All operations follow Entry 034 findings (scalar-only optimal).

use pyo3::prelude::*;
use pyo3::types::PyDict;
use crate::operations::kmer;

/// Extract overlapping k-mers from sequence
///
/// Uses evidence-based scalar implementation (Entry 034: parallel provides
/// only 2.2× for large datasets, scalar is optimal for most use cases).
///
/// Args:
///     sequence (bytes): DNA/RNA sequence
///     k (int): K-mer length (typically 3-8 for BERT models)
///
/// Returns:
///     list[bytes]: List of k-mers as byte strings
///
/// Example:
///     >>> kmers = biometal.extract_kmers(b"ATGCAT", k=3)
///     >>> print(kmers)
///     [b'ATG', b'TGC', b'GCA', b'CAT']
///
/// Note:
///     For DNABert preprocessing, typical k values are:
///     - DNABert: k=3, k=4, k=5, k=6
///     - Enformer: k=6
///     - Nucleotide Transformer: k=6
///
/// Performance:
///     This operation is scalar-only (Entry 034: NEON provides no benefit).
///     For large datasets (10K+ sequences), use KmerExtractor with parallel=True
///     for 2.2× speedup.
///
///     >>> # Streaming + constant memory
///     >>> stream = biometal.FastqStream.from_path("data.fq.gz")
///     >>> for record in stream:
///     ...     kmers = biometal.extract_kmers(record.sequence, k=6)
///     ...     # Feed to BERT model
#[pyfunction(name = "extract_kmers")]
pub fn py_extract_kmers(sequence: &[u8], k: usize) -> PyResult<Vec<Vec<u8>>> {
    if k == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "k must be greater than 0"
        ));
    }

    Ok(kmer::extract_kmers(sequence, k))
}

/// Extract minimizers from sequence (scalar-only)
///
/// Minimizers are the lexicographically smallest k-mer (by hash) in each
/// sliding window. Used for efficient indexing (minimap2-style).
///
/// Args:
///     sequence (bytes): DNA/RNA sequence
///     k (int): K-mer length
///     w (int): Window size (number of k-mers per window)
///
/// Returns:
///     list[dict]: List of minimizers with 'position', 'hash', 'kmer' keys
///
/// Example:
///     >>> minimizers = biometal.extract_minimizers(b"ATGCATGCATGC", k=3, w=5)
///     >>> for m in minimizers:
///     ...     print(f"Position {m['position']}: {m['kmer']}")
///
/// Note:
///     Entry 034 validated minimap2's scalar design (parallel provides only
///     1.26× - below ≥5× threshold).
#[pyfunction(name = "extract_minimizers")]
pub fn py_extract_minimizers(
    py: Python,
    sequence: &[u8],
    k: usize,
    w: usize
) -> PyResult<Vec<Py<PyDict>>> {
    if k == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "k must be greater than 0"
        ));
    }

    if w == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "w must be greater than 0"
        ));
    }

    let minimizers = kmer::extract_minimizers(sequence, k, w);

    let mut result = Vec::with_capacity(minimizers.len());
    for minimizer in minimizers {
        let dict = PyDict::new(py);
        dict.set_item("position", minimizer.position)?;
        dict.set_item("hash", minimizer.hash)?;
        dict.set_item("kmer", minimizer.kmer(sequence))?;
        result.push(dict.into());
    }

    Ok(result)
}

/// K-mer spectrum (frequency counting, scalar-only)
///
/// Counts k-mer frequencies across sequences. This operation is intentionally
/// scalar-only - Entry 034 found that parallelization causes HashMap contention,
/// making it SLOWER (0.95-1.88× inconsistent).
///
/// Args:
///     sequences (list[bytes]): DNA/RNA sequences to count k-mers from
///     k (int): K-mer length
///
/// Returns:
///     dict: K-mer frequencies as {kmer: count}
///
/// Example:
///     >>> sequences = [b"ATGCAT", b"GCATGC"]
///     >>> spectrum = biometal.kmer_spectrum(sequences, k=3)
///     >>> print(spectrum)
///     {b'ATG': 2, b'TGC': 2, b'GCA': 2, b'CAT': 2}
///
/// Note:
///     DO NOT attempt to parallelize - Entry 034 proved it makes performance worse!
#[pyfunction(name = "kmer_spectrum")]
pub fn py_kmer_spectrum(
    py: Python,
    sequences: Vec<Vec<u8>>,
    k: usize
) -> PyResult<Py<PyDict>> {
    if k == 0 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "k must be greater than 0"
        ));
    }

    let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();
    let spectrum = kmer::kmer_spectrum(&seq_refs, k);

    let dict = PyDict::new(py);
    for (kmer, count) in spectrum {
        dict.set_item(kmer, count)?;
    }

    Ok(dict.into())
}

/// K-mer extractor with optional parallelization (2.2× for large datasets)
///
/// Provides scalar (default) and parallel k-mer extraction. Parallel extraction
/// provides 2.19-2.38× speedup for large datasets (Entry 034).
///
/// Args:
///     parallel (bool): Enable parallel extraction (default: False)
///     threads (int): Number of threads (default: 4, capped at 4 per Entry 034)
///
/// Example:
///     >>> # Default (scalar) - simple, fast for most use cases
///     >>> extractor = biometal.KmerExtractor()
///     >>> kmers = extractor.extract(sequences, k=6)
///
///     >>> # Parallel (opt-in) - 2.2× faster for large datasets (10K+ sequences)
///     >>> extractor = biometal.KmerExtractor(parallel=True, threads=4)
///     >>> kmers = extractor.extract(large_sequences, k=6)
///
/// Note:
///     Parallel is only used for ≥1000 sequences (amortize overhead).
///     Threads are automatically capped at 4 (Entry 034: optimal).
#[pyclass(name = "KmerExtractor")]
pub struct PyKmerExtractor {
    extractor: kmer::KmerExtractor,
}

#[pymethods]
impl PyKmerExtractor {
    /// Create a new k-mer extractor
    ///
    /// Args:
    ///     parallel (bool): Enable parallel extraction (default: False)
    ///     threads (int): Number of threads (default: 4)
    ///
    /// Returns:
    ///     KmerExtractor: Configured extractor
    #[new]
    #[pyo3(signature = (parallel=false, threads=4))]
    fn new(parallel: bool, threads: usize) -> Self {
        let extractor = if parallel {
            kmer::KmerExtractor::with_parallel(threads)
        } else {
            kmer::KmerExtractor::new()
        };

        Self { extractor }
    }

    /// Extract k-mers from multiple sequences
    ///
    /// Args:
    ///     sequences (list[bytes]): DNA sequences
    ///     k (int): K-mer size
    ///
    /// Returns:
    ///     list[bytes]: Flattened list of all k-mers from all sequences
    ///
    /// Example:
    ///     >>> extractor = biometal.KmerExtractor()
    ///     >>> sequences = [b"ATGCAT", b"GCATGC"]
    ///     >>> kmers = extractor.extract(sequences, k=3)
    fn extract(&self, sequences: Vec<Vec<u8>>, k: usize) -> PyResult<Vec<Vec<u8>>> {
        if k == 0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "k must be greater than 0"
            ));
        }

        let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();
        Ok(self.extractor.extract(&seq_refs, k))
    }

    fn __repr__(&self) -> String {
        // Note: Fields are private, just return a generic repr
        "KmerExtractor()".to_string()
    }
}
