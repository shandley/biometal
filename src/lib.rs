//! biometal: ARM-native bioinformatics library with streaming architecture
//!
//! # Overview
//!
//! biometal enables analysis of terabyte-scale genomics datasets on consumer hardware
//! through streaming architecture, ARM NEON optimization, and network streaming.
//!
//! ## Key Features
//!
//! - **Streaming**: Constant ~5 MB memory regardless of dataset size
//! - **ARM-Native**: 16-25× speedup using NEON SIMD
//! - **Network Streaming**: Analyze without downloading (HTTP/SRA)
//! - **Intelligent I/O**: 16.3× speedup (parallel bgzip + mmap)
//! - **Evidence-Based**: Every optimization validated experimentally
//!
//! ## Quick Start
//!
//! ```no_run
//! use biometal::FastqStream;
//!
//! # fn main() -> biometal::Result<()> {
//! // Stream FASTQ from file (constant memory)
//! let stream = FastqStream::from_path("large.fq.gz")?;
//!
//! for record in stream {
//!     let record = record?;
//!     // Process one record at a time
//! }
//! # Ok(())
//! # }
//! ```
//!
//! ## Evidence Base
//!
//! biometal's design is grounded in 1,357 experiments (40,710 measurements, N=30):
//! - Full methodology: <https://github.com/shandley/apple-silicon-bio-bench>
//! - Optimization rules: See `OPTIMIZATION_RULES.md`
//!
//! ## Module Organization
//!
//! - [`io`]: Streaming parsers (FASTQ, FASTA, compression, network)
//! - [`operations`]: ARM NEON-optimized operations
//! - [`optimization`]: Auto-detection and platform tuning

#![warn(missing_docs)]
#![warn(rustdoc::missing_crate_level_docs)]

pub mod error;
pub mod io;
pub mod operations;
pub mod optimization;
pub mod types;

// Python bindings (Week 5-6)
#[cfg(feature = "python")]
pub mod python;

// Re-export commonly used types
pub use error::{BiometalError, Result};
pub use io::{FastaStream, FastqStream, PairedFastqStream};
pub use types::{FastaRecord, FastqRecord};

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
