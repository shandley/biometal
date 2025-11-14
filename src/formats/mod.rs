//! Bioinformatics file format parsers.
//!
//! This module provides streaming-first parsers for common bioinformatics formats:
//! - **BED**: Browser Extensible Data (genomic intervals)
//! - **GFA**: Graphical Fragment Assembly (assembly graphs, pangenomes)
//! - **VCF**: Variant Call Format (genetic variants)
//! - **GFF/GTF**: General Feature Format (gene annotations)
//!
//! All parsers implement:
//! - **Streaming architecture**: Constant memory usage regardless of file size
//! - **Iterator-based API**: Process records one at a time
//! - **Network streaming**: Read from HTTP/HTTPS URLs
//! - **Compression support**: Transparent gzip/bgzip decompression
//! - **Error handling**: All operations return `Result` for robust error handling
//!
//! # Design Principles
//!
//! ## 1. Streaming-First
//!
//! All parsers are designed for constant memory usage:
//!
//! ```rust,ignore
//! // Bad: Loads entire file into memory
//! let records: Vec<BedRecord> = parser.collect()?;
//!
//! // Good: Constant memory streaming
//! for record in parser {
//!     let record = record?;
//!     // Process one record at a time
//! }
//! ```
//!
//! ## 2. Evidence-Based Optimization
//!
//! Optimizations are applied based on profiling and evidence:
//! - **NEON SIMD**: Applied to CPU-bound operations after profiling
//! - **cloudflare_zlib**: Used for all gzip/bgzip decompression
//! - **BufReader**: Buffered I/O for all file operations
//!
//! ## 3. Production Quality
//!
//! - No `panic!` in library code (all errors use `Result`)
//! - Property-based testing for correctness
//! - Cross-platform support (Mac ARM, Linux ARM, x86_64)
//! - Comprehensive documentation with examples
//!
//! # Module Organization
//!
//! - [`primitives`]: Shared infrastructure for tab-delimited formats
//!   - Generic parsers, field utilities, genomic types
//! - Format-specific modules (coming soon):
//!   - `bed`: BED format parser
//!   - `gfa`: GFA format parser
//!   - `vcf`: VCF format parser
//!   - `gff`: GFF3/GTF format parser

pub mod primitives;

// Re-export commonly used types
pub use primitives::{
    FormatError, GenomicInterval, Strand, TabDelimitedParser, TabDelimitedRecord,
};
