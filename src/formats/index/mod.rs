//! Index formats for genomic data files
//!
//! This module provides support for index formats that enable random access
//! to genomic data files:
//!
//! - **TBI (Tabix)**: Index for tab-delimited files (BED, VCF, GFF3)
//! - **CSI**: Coordinate-sorted index (successor to TBI, supports larger chromosomes)
//!
//! # Overview
//!
//! Tabix indexes enable O(log n) region queries on sorted, tab-delimited,
//! BGZF-compressed genomic files. They work with:
//! - BED (genomic intervals)
//! - VCF (variant calls)
//! - GFF3 (gene features)
//! - Generic tab-delimited files
//!
//! # Example
//!
//! ```no_run
//! use biometal::formats::index::TbiIndex;
//!
//! # fn main() -> biometal::Result<()> {
//! // Load tabix index
//! let index = TbiIndex::from_path("variants.vcf.gz.tbi")?;
//!
//! // Query region
//! let chunks = index.query("chr1", 1000000, 2000000)?;
//! println!("Found {} chunks for region", chunks.len());
//! # Ok(())
//! # }
//! ```

pub mod tbi;

pub use tbi::TbiIndex;
