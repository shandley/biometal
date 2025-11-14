//! Python bindings for biometal
//!
//! This module provides PyO3 wrappers for biometal's core functionality,
//! enabling Python users to leverage ARM-native streaming and NEON acceleration.

use pyo3::prelude::*;

mod records;
mod streams;
mod operations;
mod kmers;
mod sequence;
mod record_ops;
mod trimming;
mod masking;
mod bam;
mod bed;
mod gfa;
mod vcf;
mod gff;
mod fasta;
mod tbi;

pub use records::*;
pub use streams::*;
pub use operations::*;
pub use kmers::*;
pub use sequence::*;
pub use record_ops::*;
pub use trimming::*;
pub use masking::*;
pub use bam::*;
pub use bed::*;
pub use gfa::*;
pub use vcf::*;
pub use gff::*;
pub use fasta::*;
pub use tbi::*;

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
    m.add_class::<PyBamReader>()?;

    // Register BAM types
    m.add_class::<PyBamRecord>()?;
    m.add_class::<PyBamHeader>()?;
    m.add_class::<PyReference>()?;
    m.add_class::<PyCigarOp>()?;
    m.add_class::<PyBamRegionIter>()?;
    m.add_class::<PyBaiIndex>()?;
    m.add_class::<PyBamIndexedRegionIter>()?;
    m.add_class::<PyTag>()?;
    m.add_class::<PyTagValue>()?;
    m.add_class::<PySamWriter>()?;

    // Register BED format types
    m.add_class::<PyBed3Record>()?;
    m.add_class::<PyBed6Record>()?;
    m.add_class::<PyBed12Record>()?;
    m.add_class::<PyBed3Stream>()?;
    m.add_class::<PyBed6Stream>()?;
    m.add_class::<PyBed12Stream>()?;

    // Register GFA format types
    m.add_class::<PyGfaSegment>()?;
    m.add_class::<PyGfaLink>()?;
    m.add_class::<PyGfaPath>()?;
    m.add_class::<PyGfaStream>()?;

    // Register VCF format types
    m.add_class::<PyVcfHeader>()?;
    m.add_class::<PyVcfRecord>()?;
    m.add_class::<PyVcfStream>()?;

    // Register GFF3 format types
    m.add_class::<PyGff3Record>()?;
    m.add_class::<PyGff3Stream>()?;

    // Register FASTA index
    m.add_class::<PyFaiIndex>()?;

    // Register TBI index
    m.add_class::<PyTbiIndex>()?;

    // Register operations
    m.add_function(wrap_pyfunction!(py_gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(py_count_bases, m)?)?;
    m.add_function(wrap_pyfunction!(py_mean_quality, m)?)?;

    // Register k-mer utilities
    m.add_function(wrap_pyfunction!(py_extract_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(py_extract_minimizers, m)?)?;
    m.add_function(wrap_pyfunction!(py_kmer_spectrum, m)?)?;

    // Register k-mer extractor class
    m.add_class::<PyKmerExtractor>()?;

    // Register sequence operations (Phase 4)
    m.add_function(wrap_pyfunction!(py_reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(py_complement, m)?)?;
    m.add_function(wrap_pyfunction!(py_reverse, m)?)?;
    m.add_function(wrap_pyfunction!(py_is_valid_dna, m)?)?;
    m.add_function(wrap_pyfunction!(py_is_valid_rna, m)?)?;
    m.add_function(wrap_pyfunction!(py_count_invalid_bases, m)?)?;

    // Register record operations (Phase 4)
    m.add_function(wrap_pyfunction!(py_extract_region, m)?)?;
    m.add_function(wrap_pyfunction!(py_reverse_complement_record, m)?)?;
    m.add_function(wrap_pyfunction!(py_sequence_length, m)?)?;
    m.add_function(wrap_pyfunction!(py_meets_length_requirement, m)?)?;
    m.add_function(wrap_pyfunction!(py_to_fasta_record, m)?)?;

    // Register trimming operations (Phase 4)
    m.add_function(wrap_pyfunction!(py_trim_start, m)?)?;
    m.add_function(wrap_pyfunction!(py_trim_end, m)?)?;
    m.add_function(wrap_pyfunction!(py_trim_both, m)?)?;
    m.add_function(wrap_pyfunction!(py_trim_quality_end, m)?)?;
    m.add_function(wrap_pyfunction!(py_trim_quality_start, m)?)?;
    m.add_function(wrap_pyfunction!(py_trim_quality_both, m)?)?;
    m.add_function(wrap_pyfunction!(py_trim_quality_window, m)?)?;

    // Register masking operations (Phase 4)
    m.add_function(wrap_pyfunction!(py_mask_low_quality, m)?)?;
    m.add_function(wrap_pyfunction!(py_count_masked_bases, m)?)?;

    // BAM statistics functions
    m.add_function(wrap_pyfunction!(calculate_coverage, m)?)?;
    m.add_function(wrap_pyfunction!(mapq_distribution, m)?)?;
    m.add_function(wrap_pyfunction!(count_by_flag, m)?)?;
    m.add_function(wrap_pyfunction!(insert_size_distribution, m)?)?;
    m.add_function(wrap_pyfunction!(edit_distance_stats, m)?)?;
    m.add_function(wrap_pyfunction!(strand_bias, m)?)?;
    m.add_function(wrap_pyfunction!(alignment_length_distribution, m)?)?;

    // Module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", "ARM-native bioinformatics library with streaming architecture")?;

    Ok(())
}
