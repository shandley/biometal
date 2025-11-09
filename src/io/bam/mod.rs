//! Native BAM (Binary Alignment Map) parser for biometal.
//!
//! This module provides a streaming, ARM-optimized BAM parser with:
//! - Constant memory footprint (~5 MB, Rule 5)
//! - Parallel BGZF decompression (6.5x speedup, Rule 3)
//! - ARM NEON sequence decoding (4.62x speedup, Rule 1) - v1.5.0
//! - Full-stack ARM-native ownership
//! - Spec-compliant parsing
//!
//! # Performance (v1.5.0)
//!
//! - **Overall throughput**: 55.1 MiB/s compressed BAM (5× faster than scalar baseline)
//! - **Record processing**: 5.82 million records/sec
//! - **Memory footprint**: Constant ~5 MB (streams terabyte-scale alignments)
//!
//! # Architecture
//!
//! The BAM implementation follows biometal's evidence-based design:
//! - **Primary optimization**: Parallel BGZF (captures 66-80% CPU bottleneck)
//! - **Secondary optimization**: ARM NEON sequence decoding (30.2% CPU time, +27.5% improvement)
//! - **Streaming-first**: Iterator interface, constant memory
//! - **ARM-native**: NEON optimizations validated by profiling (≥15% CPU time threshold)
//!
//! # Evidence Base
//!
//! Initial profiling (Phase 0, 100K records, flamegraph analysis):
//! - BGZF decompression: 66-80% CPU time → Targeted with parallel decompression
//! - Sequence decoding: <6% CPU time initially
//!
//! Post-BGZF optimization (v1.5.0, microbenchmark validation):
//! - Sequence decoding: 30.2% CPU time (exposed after BGZF optimization)
//! - NEON implementation: 4.62× speedup (4-8× typical for memory-bound operations)
//! - Overall improvement: +27.5% faster BAM parsing
//!
//! # Optimization Rules Applied
//!
//! This implementation follows evidence-based rules from `OPTIMIZATION_RULES.md`:
//!
//! - **Rule 1** (ARM NEON SIMD): Sequence decoding validated at 30.2% CPU time (v1.5.0)
//!   - Memory-bound operation: 4.62× speedup (vs 16-25× for compute-bound)
//!   - Uses `vqtbl1q_u8` for 16-entry vector table lookup
//! - **Rule 3** (Parallel BGZF): 6.5× speedup via 8-block parallel decompression
//! - **Rule 5** (Streaming): Iterator-based design maintains ~5 MB constant memory footprint
//!
//! See: `experiments/bam-simd-sequence-decoding/FINDINGS.md` for v1.5.0 NEON validation
//!
//! # Production Status (v1.5.0)
//!
//! Complete implementation:
//! - ✅ Header parsing (magic bytes, SAM text, references)
//! - ✅ Record parsing (all fields: position, MAPQ, FLAGS, sequence, quality, CIGAR)
//! - ✅ ARM NEON sequence decoding (4.62× faster, v1.5.0)
//! - ✅ CIGAR parsing (all 9 operation types)
//! - ✅ Tag parsing (typed values, convenience accessors v1.4.0)
//! - ✅ Parallel BGZF decompression (6.5× speedup)
//! - ✅ Streaming iterator interface
//! - ✅ Python bindings (v1.3.0+)
//!
//! # Example
//!
//! ```no_run
//! use biometal::io::bam::BamReader;
//! use std::fs::File;
//! use std::io::BufReader;
//!
//! # fn main() -> std::io::Result<()> {
//! let file = File::open("alignments.bam")?;
//! let reader = BufReader::new(file);
//! let mut bam = BamReader::new(reader)?;
//!
//! println!("Header: {} references", bam.header().reference_count());
//!
//! // Stream records with constant memory
//! for record in bam.records() {
//!     let record = record?;
//!     if let Some(pos) = record.position {
//!         println!("Read {} at position {}", record.name, pos);
//!     }
//! }
//! # Ok(())
//! # }
//! ```
//!
//! # Opening Files
//!
//! ```no_run
//! use biometal::io::bam::BamReader;
//!
//! # fn main() -> biometal::Result<()> {
//! // Convenience method
//! let mut bam = BamReader::from_path("alignments.bam")?;
//!
//! // Iterate over records
//! for record in bam.records() {
//!     let record = record?;
//!     println!("{}", record.name);
//! }
//! # Ok(())
//! # }
//! ```

pub mod cigar;
pub mod error;
pub mod header;
pub mod record;
pub mod reader;
pub mod sam_writer;
pub mod sequence;
pub mod tags;

// Re-export main types for convenience
pub use cigar::{CigarOp, parse_cigar};
pub use error::BamDecodeError;
pub use header::{Header, Reference};
pub use record::{Record, parse_record};
pub use reader::{BamReader, Records};
pub use sam_writer::SamWriter;
pub use sequence::decode_sequence;
pub use tags::{Tags, Tag, TagValue, ArrayValue, parse_tags};
