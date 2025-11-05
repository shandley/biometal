//! ARM NEON-optimized operations and sequence manipulation primitives
//!
//! This module provides:
//! - Rule 1 (ARM NEON SIMD): 16-25× speedup on ARM platforms with scalar fallback
//! - Sequence manipulation primitives: reverse complement, trimming, validation
//!
//! # Organization
//!
//! - `base_counting`, `gc_content`, `quality_filter`: NEON-optimized statistics
//! - `sequence`: Core sequence transformations (reverse complement, etc.)

pub mod base_counting;
pub mod gc_content;
pub mod quality_filter;
pub mod sequence;

pub use base_counting::count_bases;
pub use gc_content::{gc_content, gc_content_scalar};
pub use quality_filter::{mean_quality, mean_quality_scalar, passes_quality_filter};

// Sequence manipulation primitives
pub use sequence::{
    complement, complement_inplace, count_invalid_bases, is_valid_dna, is_valid_rna, reverse,
    reverse_complement, reverse_complement_inplace, reverse_inplace,
};
