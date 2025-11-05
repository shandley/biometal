//! Python wrappers for ARM NEON-accelerated operations

use pyo3::prelude::*;
use crate::operations;

/// Calculate GC content of DNA/RNA sequence
///
/// Uses ARM NEON SIMD for 20.3x speedup on ARM platforms.
/// Automatically falls back to scalar code on x86_64.
///
/// Args:
///     sequence (bytes): DNA/RNA sequence
///
/// Returns:
///     float: GC content as fraction (0.0 to 1.0)
///
/// Example:
///     >>> gc = biometal.gc_content(b"ATGCATGC")
///     >>> print(f"GC content: {gc:.2%}")
///     GC content: 50.00%
///
/// Performance:
///     - ARM (NEON): ~5,954 Kseq/s (20.3x speedup)
///     - x86_64 (scalar): ~294 Kseq/s
#[pyfunction(name = "gc_content")]
pub fn py_gc_content(sequence: &[u8]) -> PyResult<f64> {
    Ok(operations::gc_content(sequence))
}

/// Count bases in DNA/RNA sequence
///
/// Uses ARM NEON SIMD for 16.7x speedup on ARM platforms.
/// Automatically falls back to scalar code on x86_64.
///
/// Args:
///     sequence (bytes): DNA/RNA sequence
///
/// Returns:
///     dict: Base counts as {"A": int, "C": int, "G": int, "T": int}
///
/// Example:
///     >>> counts = biometal.count_bases(b"ATGCATGC")
///     >>> print(counts)
///     {'A': 2, 'C': 2, 'G': 2, 'T': 2}
///
/// Performance:
///     - ARM (NEON): ~5,254 Kseq/s (16.7x speedup)
///     - x86_64 (scalar): ~315 Kseq/s
#[pyfunction(name = "count_bases")]
pub fn py_count_bases(py: Python, sequence: &[u8]) -> PyResult<Py<pyo3::PyAny>> {
    let counts = operations::count_bases(sequence);

    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("A", counts[0])?;
    dict.set_item("C", counts[1])?;
    dict.set_item("G", counts[2])?;
    dict.set_item("T", counts[3])?;

    Ok(dict.into())
}

/// Calculate mean Phred quality score
///
/// Uses ARM NEON SIMD for 25.1x speedup on ARM platforms.
/// Automatically falls back to scalar code on x86_64.
///
/// Args:
///     quality (bytes): Phred quality scores (ASCII-encoded)
///
/// Returns:
///     float: Mean quality score
///
/// Example:
///     >>> quality = record.quality
///     >>> mean_q = biometal.mean_quality(quality)
///     >>> if mean_q > 30.0:
///     ...     print("High quality read")
///
/// Performance:
///     - ARM (NEON): ~6,143 Kseq/s (25.1x speedup)
///     - x86_64 (scalar): ~245 Kseq/s
#[pyfunction(name = "mean_quality")]
pub fn py_mean_quality(quality: &[u8]) -> PyResult<f64> {
    Ok(operations::mean_quality(quality))
}
