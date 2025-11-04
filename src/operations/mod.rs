//! ARM NEON-optimized operations
//!
//! This module implements Rule 1 (ARM NEON SIMD) from OPTIMIZATION_RULES.md,
//! providing 16-25× speedup on ARM platforms with automatic scalar fallback.

pub mod base_counting;
pub mod gc_content;
pub mod quality_filter;

pub use base_counting::count_bases;
pub use gc_content::{gc_content, gc_content_scalar};
pub use quality_filter::{mean_quality, mean_quality_scalar, passes_quality_filter};
