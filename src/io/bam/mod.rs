//! Native BAM (Binary Alignment Map) parser for biometal.
//!
//! This module provides a streaming, ARM-optimized BAM parser with:
//! - Constant memory footprint (~5 MB, Rule 5)
//! - Parallel BGZF decompression (6.5x speedup, Rule 3) - Phase 2
//! - Full-stack ARM-native ownership
//! - Spec-compliant parsing
//!
//! # Architecture
//!
//! The BAM implementation follows biometal's evidence-based design:
//! - **Primary optimization**: Parallel BGZF (captures 66-80% CPU bottleneck)
//! - **Streaming-first**: Iterator interface, constant memory
//! - **ARM-native**: NEON optimizations where profiling shows value (>=15% CPU time)
//!
//! # Evidence Base
//!
//! Phase 0 profiling (100K records, flamegraph analysis):
//! - BGZF decompression: 66-80% CPU time (PRIMARY TARGET)
//! - Record parsing: 6-13% CPU time
//! - Sequence decoding: <6% CPU time (NEON deferred until proven needed)
//!
//! # Optimization Rules Applied
//!
//! This implementation follows evidence-based rules from `OPTIMIZATION_RULES.md`:
//!
//! - **Rule 1** (NEON threshold): Sequence decoding <6% CPU time, below >=15% threshold for NEON optimization
//! - **Rule 3** (Parallel BGZF): Phase 2 will add 6.5x speedup via parallel decompression
//! - **Rule 5** (Streaming): Iterator-based design maintains ~5 MB constant memory footprint
//!
//! See: `experiments/native-bam-implementation/DECISION_REVISED.md` for Phase 0 profiling details
//!
//! # Phase 1 Status
//!
//! Currently implements:
//! -  Header parsing (magic bytes, SAM text, references)
//! -  Record parsing (all fields: position, MAPQ, FLAGS, sequence, quality, CIGAR)
//! -  4-bit sequence decoding (scalar, NEON deferred)
//! -  CIGAR parsing (all 9 operation types)
//! -  Tags (raw bytes, full parsing in Phase 6)
//! -  Streaming iterator interface
//!
//! Phase 2 will add:
//! - Parallel BGZF decompression (THE BIG WIN: 6.5� � ~4-5� overall)
//! - Integration with biometal compression module
//! - Benchmarking vs noodles
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
//! # fn main() -> std::io::Result<()> {
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
pub mod header;
pub mod record;
pub mod reader;
pub mod sequence;
pub mod tags;

// Re-export main types for convenience
pub use cigar::{CigarOp, parse_cigar};
pub use header::{Header, Reference};
pub use record::{Record, parse_record};
pub use reader::{BamReader, Records};
pub use sequence::decode_sequence;
pub use tags::{Tags, parse_tags};
