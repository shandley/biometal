//! I/O module: Streaming parsers and compression
//!
//! This module implements streaming architecture with constant memory (~5 MB)
//! regardless of dataset size, following Rules 2-6 from OPTIMIZATION_RULES.md.

pub mod compression;
mod fasta;
mod fastq;
mod paired;
pub mod sink;

pub use compression::{decompress_bgzip_parallel, CompressedReader, CompressedWriter, DataSource, MMAP_THRESHOLD};
pub use fasta::FastaStream;
pub use fastq::{FastqStream, FastqWriter};
pub use paired::PairedFastqStream;
pub use sink::DataSink;

// Week 3-4: Network streaming (Rule 6)
#[cfg(feature = "network")]
pub mod network;
#[cfg(feature = "network")]
pub use network::{HttpClient, HttpReader};

// Week 3-4: SRA integration (Rule 6)
#[cfg(feature = "network")]
pub mod sra;
#[cfg(feature = "network")]
pub use sra::{is_sra_accession, sra_to_url};

// Native BAM implementation (Phase 1-6)
// ARM-optimized BAM parser with parallel BGZF decompression
// See experiments/native-bam-implementation/ for design and profiling
pub mod bam;
pub use bam::BamReader;
