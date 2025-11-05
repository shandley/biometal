//! Python bindings for biometal
//!
//! This module provides PyO3 wrappers for biometal's core functionality,
//! enabling Python users to leverage ARM-native streaming and NEON acceleration.

use pyo3::prelude::*;

mod records;
mod streams;
mod operations;
mod kmers;

pub use records::*;
pub use streams::*;
pub use operations::*;
pub use kmers::*;

/// biometal: ARM-native bioinformatics library
///
/// Features:
/// - Streaming architecture (constant ~5 MB memory)
/// - ARM NEON acceleration (16-25x speedup)
/// - Network streaming (analyze without downloading)
/// - Production-grade (121 tests, Grade A+)
///
/// Example:
///     >>> import biometal
///     >>> stream = biometal.FastqStream.from_path("data.fq.gz")
///     >>> for record in stream:
///     ...     gc = biometal.gc_content(record.sequence)
///     ...     print(f"{record.id}: GC={gc:.2%}")
#[pymodule]
fn biometal(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Register record types
    m.add_class::<PyFastqRecord>()?;
    m.add_class::<PyFastaRecord>()?;

    // Register streaming classes
    m.add_class::<PyFastqStream>()?;
    m.add_class::<PyFastaStream>()?;

    // Register operations
    m.add_function(wrap_pyfunction!(py_gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(py_count_bases, m)?)?;
    m.add_function(wrap_pyfunction!(py_mean_quality, m)?)?;

    // Register k-mer utilities
    m.add_function(wrap_pyfunction!(py_extract_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(py_extract_kmers_non_overlapping, m)?)?;

    // Module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", "ARM-native bioinformatics library with streaming architecture")?;

    Ok(())
}
