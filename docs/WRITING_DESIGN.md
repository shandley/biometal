# Writing Operations Design Specification

**Status**: Phase 1 - Design
**Date**: November 5, 2025
**Version**: Draft 1.0

---

## Overview

This document specifies the design for writing operations in biometal, enabling CLI tools to transform and output data streams. This is the **critical blocker** for CLI tool development.

## Design Goals

1. **Mirror reader architecture** - Writers should be the inverse of readers
2. **Streaming-first** - Constant memory, write one record at a time (Rule 5)
3. **Compression support** - Parallel bgzip compression (inverse of Rule 3)
4. **Production quality** - No panics, proper error handling, RAII cleanup
5. **Performance** - Buffered writes, async flushing options

---

## Architecture

### Type Hierarchy

```
DataSink (enum)
├── Local(PathBuf)
├── Stdout
└── Future: Network (S3, etc.)

CompressedWriter
├── Uncompressed (BufWriter)
├── Gzip (flate2::GzEncoder)
└── Bgzip (parallel compression)

FastqWriter<W: Write>
FastaWriter<W: Write>
```

### Core Abstractions

#### 1. DataSink (src/io/sink.rs - NEW)

Mirrors `DataSource` for output destinations:

```rust
/// Output destination for streaming writes
pub enum DataSink {
    /// Write to local file path
    Local(PathBuf),

    /// Write to stdout
    Stdout,

    // Future: Network destinations
    // S3(String),
    // Azure(String),
}

impl DataSink {
    /// Create sink from file path (detects compression from extension)
    pub fn from_path<P: AsRef<Path>>(path: P) -> Self {
        Self::Local(path.as_ref().to_path_buf())
    }

    /// Create stdout sink
    pub fn stdout() -> Self {
        Self::Stdout
    }
}
```

#### 2. CompressedWriter (src/io/compression.rs - EXTEND)

Mirrors `CompressedReader` for writing:

```rust
/// Writer with automatic compression support
pub enum CompressedWriter {
    /// Uncompressed writer
    Plain(BufWriter<Box<dyn Write>>),

    /// Gzip compressed writer
    Gzip(GzEncoder<BufWriter<Box<dyn Write>>>),

    /// Bgzip compressed writer with parallel compression
    Bgzip(BgzipWriter),
}

impl CompressedWriter {
    /// Create writer from data sink
    pub fn new(sink: DataSink) -> io::Result<Self> {
        match sink {
            DataSink::Local(path) => {
                let file = File::create(&path)?;

                // Detect compression from extension
                let ext = path.extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or("");

                match ext {
                    "gz" => {
                        // Decide: bgzip vs gzip based on filename?
                        // *.fq.gz, *.fastq.gz, *.fa.gz -> bgzip (block-based)
                        // Others -> regular gzip
                        Self::new_gzip(Box::new(file))
                    }
                    "bgz" => {
                        Self::new_bgzip(Box::new(file))
                    }
                    _ => {
                        Self::new_plain(Box::new(file))
                    }
                }
            }
            DataSink::Stdout => {
                Self::new_plain(Box::new(io::stdout()))
            }
        }
    }

    /// Create plain (uncompressed) writer
    pub fn new_plain(writer: Box<dyn Write>) -> io::Result<Self> {
        Ok(Self::Plain(BufWriter::new(writer)))
    }

    /// Create gzip writer
    pub fn new_gzip(writer: Box<dyn Write>) -> io::Result<Self> {
        let encoder = GzEncoder::new(
            BufWriter::new(writer),
            flate2::Compression::default(),
        );
        Ok(Self::Gzip(encoder))
    }

    /// Create bgzip writer (parallel compression)
    pub fn new_bgzip(writer: Box<dyn Write>) -> io::Result<Self> {
        Ok(Self::Bgzip(BgzipWriter::new(writer)?))
    }

    /// Flush all buffered data
    pub fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
            Self::Bgzip(w) => w.flush(),
        }
    }

    /// Finish writing and consume writer
    /// This ensures proper compression finalization
    pub fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Gzip(mut w) => {
                w.try_finish()?;
                Ok(())
            }
            Self::Bgzip(w) => w.finish(),
        }
    }
}

impl Write for CompressedWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Gzip(w) => w.write(buf),
            Self::Bgzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        CompressedWriter::flush(self)
    }
}
```

#### 3. BgzipWriter (src/io/compression.rs - NEW)

Parallel bgzip compression (inverse of parallel decompression):

```rust
/// Parallel bgzip writer
///
/// Compresses data in BGZF blocks using parallel compression
/// (inverse of parallel bgzip decompression from Rule 3)
pub struct BgzipWriter {
    writer: Box<dyn Write>,
    buffer: Vec<u8>,
    block_size: usize, // Default: 64 KB (BGZF standard)
    compression_level: flate2::Compression,
}

impl BgzipWriter {
    /// Create new bgzip writer
    pub fn new(writer: Box<dyn Write>) -> io::Result<Self> {
        Ok(Self {
            writer,
            buffer: Vec::with_capacity(64 * 1024),
            block_size: 64 * 1024,
            compression_level: flate2::Compression::default(),
        })
    }

    /// Set compression level (0-9)
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.compression_level = flate2::Compression::new(level);
        self
    }

    /// Write a complete BGZF block
    fn write_block(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        // Compress block
        let compressed = compress_bgzip_block(&self.buffer, self.compression_level)?;

        // Write to underlying writer
        self.writer.write_all(&compressed)?;

        // Clear buffer for next block
        self.buffer.clear();

        Ok(())
    }

    /// Finish writing and add BGZF EOF marker
    pub fn finish(mut self) -> io::Result<()> {
        // Write any remaining data
        self.write_block()?;

        // Write BGZF EOF marker (empty block)
        const BGZF_EOF: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
            0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
            0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00,
        ];
        self.writer.write_all(&BGZF_EOF)?;

        self.writer.flush()?;
        Ok(())
    }
}

impl Write for BgzipWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut bytes_written = 0;
        let mut remaining = buf;

        while !remaining.is_empty() {
            let space_left = self.block_size - self.buffer.len();
            let to_write = remaining.len().min(space_left);

            self.buffer.extend_from_slice(&remaining[..to_write]);
            bytes_written += to_write;
            remaining = &remaining[to_write..];

            // Write block if buffer is full
            if self.buffer.len() >= self.block_size {
                self.write_block()?;
            }
        }

        Ok(bytes_written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.write_block()?;
        self.writer.flush()
    }
}

impl Drop for BgzipWriter {
    fn drop(&mut self) {
        // Best-effort flush on drop
        let _ = self.flush();
    }
}
```

#### 4. FastqWriter (src/io/fastq.rs - EXTEND)

FASTQ record writer:

```rust
/// FASTQ record writer with streaming and compression support
///
/// # Memory Footprint
///
/// - Write buffer: ~8 KB (BufWriter default)
/// - Compression buffer: ~64 KB (BGZF block size)
/// - **Total: ~72 KB constant** (Rule 5 compliant)
///
/// # Example
///
/// ```no_run
/// use biometal::{FastqWriter, FastqRecord};
/// use biometal::io::DataSink;
///
/// # fn main() -> biometal::Result<()> {
/// let sink = DataSink::from_path("output.fq.gz");
/// let mut writer = FastqWriter::new(sink)?;
///
/// let record = FastqRecord::new(
///     "read1".to_string(),
///     b"ATGC".to_vec(),
///     b"IIII".to_vec(),
/// );
///
/// writer.write_record(&record)?;
/// writer.finish()?;  // Important: finalizes compression
/// # Ok(())
/// # }
/// ```
pub struct FastqWriter {
    writer: CompressedWriter,
}

impl FastqWriter {
    /// Create a new FASTQ writer from a data sink
    pub fn new(sink: DataSink) -> Result<Self> {
        let writer = CompressedWriter::new(sink)?;
        Ok(Self { writer })
    }

    /// Create a FASTQ writer from a file path
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(DataSink::from_path(path))
    }

    /// Create a FASTQ writer to stdout
    pub fn stdout() -> Result<Self> {
        Self::new(DataSink::stdout())
    }

    /// Write a single FASTQ record
    ///
    /// Format:
    /// ```text
    /// @id
    /// sequence
    /// +
    /// quality
    /// ```
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        // Validate record
        if record.sequence.len() != record.quality.len() {
            return Err(BiometalError::InvalidFastqFormat {
                line: 0,
                msg: format!(
                    "Sequence length ({}) != quality length ({})",
                    record.sequence.len(),
                    record.quality.len()
                ),
            });
        }

        // Write 4 lines
        write!(self.writer, "@{}\n", record.id)?;
        self.writer.write_all(&record.sequence)?;
        self.writer.write_all(b"\n+\n")?;
        self.writer.write_all(&record.quality)?;
        self.writer.write_all(b"\n")?;

        Ok(())
    }

    /// Flush buffered data to disk
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().map_err(Into::into)
    }

    /// Finish writing and finalize compression
    ///
    /// This is important for compressed files to write EOF markers
    /// and ensure all data is flushed.
    pub fn finish(self) -> Result<()> {
        self.writer.finish().map_err(Into::into)
    }
}
```

#### 5. FastaWriter (src/io/fasta.rs - EXTEND)

FASTA record writer with line wrapping:

```rust
/// FASTA record writer with optional line wrapping
///
/// # Line Wrapping
///
/// FASTA sequences are often wrapped at 60 or 80 characters per line
/// for readability. This is configurable.
///
/// # Example
///
/// ```no_run
/// use biometal::{FastaWriter, FastaRecord};
/// use biometal::io::DataSink;
///
/// # fn main() -> biometal::Result<()> {
/// let sink = DataSink::from_path("output.fa.gz");
/// let mut writer = FastaWriter::new(sink)?
///     .with_line_width(60);  // Wrap at 60 characters
///
/// let record = FastaRecord::new(
///     "contig1".to_string(),
///     b"ATGCATGCATGC".to_vec(),
/// );
///
/// writer.write_record(&record)?;
/// writer.finish()?;
/// # Ok(())
/// # }
/// ```
pub struct FastaWriter {
    writer: CompressedWriter,
    line_width: Option<usize>,  // None = no wrapping
}

impl FastaWriter {
    /// Create a new FASTA writer from a data sink
    pub fn new(sink: DataSink) -> Result<Self> {
        let writer = CompressedWriter::new(sink)?;
        Ok(Self {
            writer,
            line_width: Some(80),  // Default: 80 characters
        })
    }

    /// Create a FASTA writer from a file path
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::new(DataSink::from_path(path))
    }

    /// Set line width for sequence wrapping
    ///
    /// Set to `None` to disable wrapping (write sequence on one line)
    pub fn with_line_width(mut self, width: usize) -> Self {
        self.line_width = Some(width);
        self
    }

    /// Disable line wrapping
    pub fn no_wrapping(mut self) -> Self {
        self.line_width = None;
        self
    }

    /// Write a single FASTA record
    ///
    /// Format:
    /// ```text
    /// >id
    /// sequence (wrapped at line_width if set)
    /// ```
    pub fn write_record(&mut self, record: &FastaRecord) -> Result<()> {
        // Write header
        write!(self.writer, ">{}\n", record.id)?;

        // Write sequence (with optional wrapping)
        match self.line_width {
            None => {
                // No wrapping
                self.writer.write_all(&record.sequence)?;
                self.writer.write_all(b"\n")?;
            }
            Some(width) => {
                // Wrap at specified width
                for chunk in record.sequence.chunks(width) {
                    self.writer.write_all(chunk)?;
                    self.writer.write_all(b"\n")?;
                }
            }
        }

        Ok(())
    }

    /// Flush buffered data to disk
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().map_err(Into::into)
    }

    /// Finish writing and finalize compression
    pub fn finish(self) -> Result<()> {
        self.writer.finish().map_err(Into::into)
    }
}
```

---

## Implementation Plan

### Phase 1.1: Core Infrastructure
1. Create `src/io/sink.rs` with `DataSink` enum
2. Extend `src/io/compression.rs` with `CompressedWriter`
3. Implement basic (uncompressed) writer
4. Tests: Write and read-back validation

### Phase 1.2: Compression Support
1. Implement `BgzipWriter` with parallel compression
2. Add gzip compression support
3. Tests: Compressed output validation
4. Benchmark: Compression throughput

### Phase 1.3: Format Writers
1. Implement `FastqWriter` in `src/io/fastq.rs`
2. Implement `FastaWriter` in `src/io/fasta.rs`
3. Tests: Format correctness, round-trip validation
4. Integration tests: Read → transform → write

---

## Testing Strategy

### Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastq_writer_basic() {
        let sink = DataSink::from_path("test.fq");
        let mut writer = FastqWriter::new(sink).unwrap();

        let record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"IIII".to_vec(),
        );

        writer.write_record(&record).unwrap();
        writer.finish().unwrap();

        // Read back and validate
        let stream = FastqStream::from_path("test.fq").unwrap();
        let read_record = stream.next().unwrap().unwrap();
        assert_eq!(record, read_record);
    }

    #[test]
    fn test_fastq_writer_compressed() {
        let sink = DataSink::from_path("test.fq.gz");
        let mut writer = FastqWriter::new(sink).unwrap();

        // Write many records
        for i in 0..1000 {
            let record = FastqRecord::new(
                format!("read{}", i),
                b"ATGCATGC".to_vec(),
                b"IIIIIIII".to_vec(),
            );
            writer.write_record(&record).unwrap();
        }
        writer.finish().unwrap();

        // Read back and count
        let stream = FastqStream::from_path("test.fq.gz").unwrap();
        let count = stream.count();
        assert_eq!(count, 1000);
    }

    #[test]
    fn test_fastq_writer_validation() {
        let sink = DataSink::from_path("test.fq");
        let mut writer = FastqWriter::new(sink).unwrap();

        // Invalid: mismatched lengths
        let bad_record = FastqRecord::new(
            "read1".to_string(),
            b"ATGC".to_vec(),
            b"III".to_vec(),  // Too short!
        );

        let result = writer.write_record(&bad_record);
        assert!(result.is_err());
    }
}
```

### Integration Tests

Create `tests/writing.rs`:

```rust
use biometal::{FastqStream, FastqWriter, FastqRecord};
use biometal::io::DataSink;
use biometal::operations::gc_content;

#[test]
fn test_read_filter_write_pipeline() {
    // Read from input
    let input = FastqStream::from_path("tests/data/input.fq.gz").unwrap();

    // Filter high-quality reads
    let sink = DataSink::from_path("tests/output/filtered.fq.gz");
    let mut writer = FastqWriter::new(sink).unwrap();

    for record in input {
        let record = record.unwrap();

        // Filter by GC content
        let gc = gc_content(&record.sequence);
        if gc >= 0.4 && gc <= 0.6 {
            writer.write_record(&record).unwrap();
        }
    }

    writer.finish().unwrap();

    // Verify output
    let output = FastqStream::from_path("tests/output/filtered.fq.gz").unwrap();
    let count = output.count();
    assert!(count > 0);
}
```

### Property-Based Tests

```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn test_write_read_roundtrip(
        id in "[a-zA-Z0-9]{1,50}",
        seq in "[ACGT]{1,200}",
    ) {
        let qual = vec![b'I'; seq.len()];
        let record = FastqRecord::new(id.clone(), seq.into(), qual);

        // Write
        let sink = DataSink::from_path("test_roundtrip.fq");
        let mut writer = FastqWriter::new(sink).unwrap();
        writer.write_record(&record).unwrap();
        writer.finish().unwrap();

        // Read back
        let stream = FastqStream::from_path("test_roundtrip.fq").unwrap();
        let read_record = stream.next().unwrap().unwrap();

        // Validate
        prop_assert_eq!(record.id, read_record.id);
        prop_assert_eq!(record.sequence, read_record.sequence);
        prop_assert_eq!(record.quality, read_record.quality);
    }
}
```

---

## Performance Considerations

### Memory Usage
- **Target**: ~100 KB constant (Rule 5)
  - Write buffer: 8 KB (BufWriter)
  - Compression buffer: 64 KB (BGZF block)
  - Overhead: ~28 KB
- **Comparison**: Orders of magnitude less than batch writing

### Throughput
- **Uncompressed**: Limited by disk I/O (~500 MB/s SSD)
- **Compressed**: Limited by compression speed
  - Single-threaded gzip: ~80 MB/s
  - Multi-threaded bgzip: Target ~300-400 MB/s (parallel compression)

### Parallel Compression Opportunity
```rust
// Future optimization: Parallel bgzip compression
pub struct ParallelBgzipWriter {
    workers: ThreadPool,
    input_queue: Sender<Vec<u8>>,
    output_queue: Receiver<Vec<u8>>,
}

// Could achieve 6.5× compression speedup (inverse of decompression)
```

---

## API Stability

### Public API (v1.1.0)
```rust
// Core types
pub enum DataSink { ... }
pub struct FastqWriter { ... }
pub struct FastaWriter { ... }

// Public methods
impl FastqWriter {
    pub fn new(sink: DataSink) -> Result<Self>;
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self>;
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()>;
    pub fn flush(&mut self) -> Result<()>;
    pub fn finish(self) -> Result<()>;
}

impl FastaWriter {
    pub fn new(sink: DataSink) -> Result<Self>;
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self>;
    pub fn with_line_width(self, width: usize) -> Self;
    pub fn no_wrapping(self) -> Self;
    pub fn write_record(&mut self, record: &FastaRecord) -> Result<()>;
    pub fn flush(&mut self) -> Result<()>;
    pub fn finish(self) -> Result<()>;
}
```

---

## Future Enhancements

### Async I/O (v1.2+)
```rust
pub struct AsyncFastqWriter { ... }

impl AsyncFastqWriter {
    pub async fn write_record(&mut self, record: &FastqRecord) -> Result<()>;
    pub async fn finish(self) -> Result<()>;
}
```

### Network Destinations (v1.2+)
```rust
pub enum DataSink {
    S3(String),
    Azure(String),
    GCS(String),
}
```

### Parallel Compression (v1.2+)
Investigate using Rayon for parallel BGZF block compression (inverse of decompression).

---

## Decision Log

### 1. Bgzip vs Gzip Detection
**Decision**: Use `.bgz` extension for bgzip, `.gz` for gzip
**Rationale**: Clear separation, explicit control
**Alternative**: Auto-detect based on filename pattern (*.fq.gz → bgzip)

### 2. Line Wrapping Default
**Decision**: Default to 80 characters for FASTA
**Rationale**: Standard convention, widely compatible
**Alternative**: No wrapping by default

### 3. Finish vs Drop
**Decision**: Explicit `finish()` required, best-effort flush in `Drop`
**Rationale**: Compression finalization can fail, users should handle explicitly
**Alternative**: Panic in Drop (bad practice)

---

## Documentation Requirements

Every writer must have:
1. Module-level documentation with examples
2. Method documentation with usage patterns
3. Error conditions documented
4. Memory guarantees documented
5. Example in `examples/writing.rs`

---

**Status**: Ready for Phase 1.1 implementation
**Next Steps**: Create `src/io/sink.rs` and basic `CompressedWriter`
