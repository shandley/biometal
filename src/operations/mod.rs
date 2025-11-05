//! ARM NEON-optimized operations and sequence manipulation primitives
//!
//! This module provides:
//! - Rule 1 (ARM NEON SIMD): 16-25× speedup on ARM platforms with scalar fallback
//! - Sequence manipulation primitives: reverse complement, trimming, masking, validation
//! - Record-level operations: extract regions, reverse complement records
//! - Trimming operations: fixed-position and quality-based trimming
//! - Masking operations: replace low-quality bases with 'N'
//!
//! # Organization
//!
//! - `base_counting`, `gc_content`, `quality_filter`: NEON-optimized statistics
//! - `sequence`: Core sequence transformations (reverse complement, etc.)
//! - `record_ops`: Record-level operations (extract_region, etc.)
//! - `trimming`: Fixed and quality-based trimming operations
//! - `masking`: Quality-based masking operations

pub mod base_counting;
pub mod gc_content;
pub mod masking;
pub mod quality_filter;
pub mod record_ops;
pub mod sequence;
pub mod trimming;

pub use base_counting::count_bases;
pub use gc_content::{gc_content, gc_content_scalar};
pub use quality_filter::{mean_quality, mean_quality_scalar, passes_quality_filter};

// Sequence manipulation primitives
pub use sequence::{
    complement, complement_inplace, count_invalid_bases, is_valid_dna, is_valid_rna, reverse,
    reverse_complement, reverse_complement_inplace, reverse_inplace,
};

// Record-level operations
pub use record_ops::{
    extract_region, meets_length_requirement, reverse_complement_record,
    reverse_complement_record_inplace, sequence_length,
};

// Trimming operations
pub use trimming::{
    trim_both, trim_end, trim_quality_both, trim_quality_end, trim_quality_start,
    trim_quality_window, trim_start,
};

// Masking operations
pub use masking::{count_masked_bases, mask_low_quality, mask_low_quality_copy};
