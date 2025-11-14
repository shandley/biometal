//! FASTA format support: streaming parser and FAI index
//!
//! This module provides:
//! - Streaming FASTA parser with constant memory (Rule 5)
//! - FAI (FASTA Index) for random access to sequences
//!
//! # Basic Usage
//!
//! ```no_run
//! use biometal::io::fasta::FastaStream;
//!
//! // Stream through FASTA file
//! let stream = FastaStream::from_path("genome.fa.gz")?;
//! for record in stream {
//!     let record = record?;
//!     println!("{}: {} bp", record.id, record.sequence.len());
//! }
//! # Ok::<(), biometal::error::BiometalError>(())
//! ```
//!
//! # Indexed Access (with FAI)
//!
//! ```no_run
//! use biometal::io::fasta::{FaiIndex, FastaStream};
//!
//! // Load index
//! let index = FaiIndex::from_path("genome.fa.fai")?;
//!
//! // Query specific sequence by name
//! let seq = index.fetch("chr1", "genome.fa")?;
//! println!("chr1: {} bp", seq.len());
//!
//! // Query subsequence (chr1:1000-2000)
//! let subseq = index.fetch_region("chr1", 1000, 2000, "genome.fa")?;
//! println!("Region: {} bp", subseq.len());
//! # Ok::<(), biometal::error::BiometalError>(())
//! ```

mod parser;
pub mod index;

pub use parser::FastaStream;
pub use index::FaiIndex;
