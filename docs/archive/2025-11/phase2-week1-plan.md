# Phase 2: Week 1 - Format Primitives Foundation

**Week**: 1 of 16 (Nov 13-20, 2025)
**Focus**: Design and implement shared format infrastructure
**Goal**: Reusable primitives for tab-delimited formats (BED, GFF, VCF, GFA)

---

## Context

Before implementing specific formats (BED, GFA, VCF, GFF3), we need shared infrastructure. All these formats share common patterns:
- Tab-delimited fields
- Header/comment lines (starting with `#`)
- Field validation and parsing
- Genomic coordinates (chromosome, position, strand)

Building primitives once → faster format implementation later.

---

## Week 1 Deliverables

### 1. Module Structure

**Create**:
```
src/formats/
├── mod.rs               # Public API, re-exports
├── primitives/
│   ├── mod.rs           # Primitives module
│   ├── tab_delimited.rs # Generic tab-delimited parser
│   ├── header.rs        # Header/comment handling
│   ├── fields.rs        # Field parsing utilities
│   └── genomic.rs       # Genomic types (Strand, Interval, etc.)
└── README.md            # Format library overview
```

### 2. Core Types and Traits

#### GenomicInterval (genomic.rs)
```rust
/// Genomic interval (0-based, half-open [start, end))
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicInterval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl GenomicInterval {
    pub fn new(chrom: String, start: u64, end: u64) -> Result<Self>;
    pub fn length(&self) -> u64;
    pub fn overlaps(&self, other: &Self) -> bool;
    pub fn contains(&self, other: &Self) -> bool;
}
```

#### Strand (genomic.rs)
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl FromStr for Strand {
    // Parse '+', '-', '.'
}
```

#### TabDelimitedRecord (tab_delimited.rs)
```rust
/// Generic tab-delimited record
pub trait TabDelimitedRecord: Sized {
    /// Parse from tab-delimited line
    fn from_line(line: &str) -> Result<Self>;

    /// Serialize to tab-delimited line
    fn to_line(&self) -> String;

    /// Expected number of fields (None = variable)
    fn expected_fields() -> Option<usize> { None }
}
```

#### TabDelimitedParser (tab_delimited.rs)
```rust
/// Generic streaming parser for tab-delimited formats
pub struct TabDelimitedParser<R: Read, T: TabDelimitedRecord> {
    reader: BufReader<R>,
    line_buf: String,
    line_number: usize,
    _phantom: PhantomData<T>,
}

impl<R: Read, T: TabDelimitedRecord> TabDelimitedParser<R, T> {
    pub fn new(reader: R) -> Self;
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self>;
    pub fn from_url(url: &str) -> Result<Self>;  // HTTP streaming
    pub fn from_bgzip_path(path: impl AsRef<Path>) -> Result<Self>;  // bgzip support
}

impl<R: Read, T: TabDelimitedRecord> Iterator for TabDelimitedParser<R, T> {
    type Item = Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buf.clear();
            match self.reader.read_line(&mut self.line_buf) {
                Ok(0) => return None,
                Ok(_) => {
                    self.line_number += 1;

                    // Skip empty lines
                    if self.line_buf.trim().is_empty() {
                        continue;
                    }

                    // Skip comments (handled by HeaderParser)
                    if self.line_buf.starts_with('#') {
                        continue;
                    }

                    return Some(T::from_line(&self.line_buf));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}
```

#### HeaderParser (header.rs)
```rust
/// Parse header/comment lines from tab-delimited formats
pub struct HeaderParser<R: Read> {
    reader: BufReader<R>,
    line_buf: String,
}

impl<R: Read> HeaderParser<R> {
    pub fn new(reader: R) -> Self;

    /// Parse all header lines (lines starting with '#')
    /// Returns (headers, remaining_reader) for data parsing
    pub fn parse_headers(mut self) -> Result<(Vec<String>, R)>;
}
```

#### Field Parsing Utilities (fields.rs)
```rust
/// Parse required field
pub fn parse_required<T: FromStr>(
    field: &str,
    field_name: &str,
) -> Result<T>
where
    T::Err: std::error::Error + Send + Sync + 'static;

/// Parse optional field ('.', '*', or empty → None)
pub fn parse_optional<T: FromStr>(
    field: &str,
    field_name: &str,
) -> Result<Option<T>>
where
    T::Err: std::error::Error + Send + Sync + 'static;

/// Parse key=value pairs (e.g., GFF attributes, VCF INFO)
pub fn parse_attributes(attr_str: &str) -> Result<HashMap<String, String>>;

/// Split tab-delimited line with validation
pub fn split_fields(line: &str, expected: Option<usize>) -> Result<Vec<&str>>;
```

### 3. Error Handling

```rust
// src/formats/primitives/mod.rs

#[derive(Debug, thiserror::Error)]
pub enum FormatError {
    #[error("Invalid number of fields: expected {expected}, got {actual}")]
    FieldCount { expected: usize, actual: usize },

    #[error("Invalid field '{field}': {reason}")]
    InvalidField { field: String, reason: String },

    #[error("Parse error at line {line}: {source}")]
    ParseError { line: usize, source: Box<dyn std::error::Error + Send + Sync> },

    #[error("Invalid genomic interval: start {start} >= end {end}")]
    InvalidInterval { start: u64, end: u64 },

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("HTTP error: {0}")]
    Http(#[from] reqwest::Error),
}

pub type Result<T> = std::result::Result<T, FormatError>;
```

---

## Implementation Steps

### Session 1: Module Structure + Types (2-3 hours)

1. **Create module structure**
   - `src/formats/mod.rs`
   - `src/formats/primitives/mod.rs`
   - Empty modules for each primitive

2. **Implement genomic types** (`genomic.rs`)
   - `GenomicInterval` struct
   - `Strand` enum
   - Tests for overlaps, contains, parsing

3. **Design errors** (`primitives/mod.rs`)
   - `FormatError` enum with thiserror
   - `Result<T>` type alias

### Session 2: Tab-Delimited Parser (2-3 hours)

1. **Implement `TabDelimitedRecord` trait** (`tab_delimited.rs`)
   - Trait definition
   - Documentation with examples

2. **Implement `TabDelimitedParser`** (`tab_delimited.rs`)
   - Generic streaming parser
   - Iterator implementation
   - Skip comments/empty lines
   - Line number tracking for errors

3. **Add constructors**
   - `from_path()`
   - `from_url()` (HTTP streaming)
   - `from_bgzip_path()` (compressed)

### Session 3: Field Parsing + Headers (2-3 hours)

1. **Implement field utilities** (`fields.rs`)
   - `parse_required()`
   - `parse_optional()`
   - `parse_attributes()`
   - `split_fields()`

2. **Implement `HeaderParser`** (`header.rs`)
   - Parse lines starting with `#`
   - Return headers + remaining reader

3. **Integration tests**
   - Test with sample tab-delimited data
   - Test header parsing
   - Test field parsing

### Session 4: Testing + Documentation (1-2 hours)

1. **Property tests** (proptest)
   - Round-trip: `GenomicInterval` serialization
   - Field parsing edge cases
   - Strand parsing

2. **Documentation**
   - API docs with examples
   - Module-level overview
   - Usage patterns

3. **Examples**
   - `examples/format_primitives_demo.rs`

---

## Testing Strategy

### Unit Tests

```rust
// src/formats/primitives/genomic.rs

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genomic_interval_length() {
        let interval = GenomicInterval::new("chr1".to_string(), 100, 200).unwrap();
        assert_eq!(interval.length(), 100);
    }

    #[test]
    fn test_genomic_interval_invalid() {
        let result = GenomicInterval::new("chr1".to_string(), 200, 100);
        assert!(result.is_err());
    }

    #[test]
    fn test_strand_parsing() {
        assert_eq!(Strand::from_str("+").unwrap(), Strand::Forward);
        assert_eq!(Strand::from_str("-").unwrap(), Strand::Reverse);
        assert_eq!(Strand::from_str(".").unwrap(), Strand::Unknown);
    }
}
```

### Property Tests

```rust
// src/formats/primitives/genomic.rs

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_interval_overlaps_symmetric(
            start1 in 0u64..10000,
            end1 in 0u64..10000,
            start2 in 0u64..10000,
            end2 in 0u64..10000,
        ) {
            let start1 = start1.min(end1);
            let end1 = end1.max(start1 + 1);
            let start2 = start2.min(end2);
            let end2 = end2.max(start2 + 1);

            let int1 = GenomicInterval::new("chr1".to_string(), start1, end1).unwrap();
            let int2 = GenomicInterval::new("chr1".to_string(), start2, end2).unwrap();

            // Overlaps should be symmetric
            assert_eq!(int1.overlaps(&int2), int2.overlaps(&int1));
        }
    }
}
```

### Integration Tests

```rust
// tests/format_primitives.rs

use biometal::formats::primitives::*;

#[derive(Debug, PartialEq)]
struct TestRecord {
    chrom: String,
    start: u64,
    end: u64,
}

impl TabDelimitedRecord for TestRecord {
    fn from_line(line: &str) -> Result<Self> {
        let fields = split_fields(line, Some(3))?;
        Ok(TestRecord {
            chrom: fields[0].to_string(),
            start: parse_required(fields[1], "start")?,
            end: parse_required(fields[2], "end")?,
        })
    }

    fn to_line(&self) -> String {
        format!("{}\t{}\t{}", self.chrom, self.start, self.end)
    }

    fn expected_fields() -> Option<usize> {
        Some(3)
    }
}

#[test]
fn test_tab_delimited_parser() {
    let data = "chr1\t100\t200\nchr2\t300\t400\n";
    let parser = TabDelimitedParser::<_, TestRecord>::new(data.as_bytes());

    let records: Vec<_> = parser.collect::<Result<_>>().unwrap();

    assert_eq!(records.len(), 2);
    assert_eq!(records[0].chrom, "chr1");
    assert_eq!(records[1].start, 300);
}
```

---

## Success Criteria

**Code**:
- [ ] `src/formats/primitives/` module created
- [ ] `GenomicInterval` implemented with tests
- [ ] `Strand` implemented with tests
- [ ] `TabDelimitedRecord` trait defined
- [ ] `TabDelimitedParser` implemented (streaming)
- [ ] `HeaderParser` implemented
- [ ] Field parsing utilities implemented
- [ ] HTTP streaming support (`from_url`)
- [ ] Compression support (`from_bgzip_path`)

**Testing**:
- [ ] Unit tests for genomic types
- [ ] Property tests for interval operations
- [ ] Integration tests for tab-delimited parsing
- [ ] Integration tests for header parsing
- [ ] Test with sample data files

**Documentation**:
- [ ] API docs for all public items
- [ ] Module-level overview
- [ ] Usage examples in docs
- [ ] Example program (`examples/format_primitives_demo.rs`)

**Quality**:
- [ ] No panics (all errors use `Result`)
- [ ] No clippy warnings
- [ ] Formatted with rustfmt
- [ ] Tests pass on Mac ARM, Linux ARM, x86_64

---

## Optimization Considerations

### Apply Proven Optimizations

**Streaming (Rule 5)**: ✅
- Iterator-based parsing
- Line buffer reuse (`line_buf.clear()`)
- Constant memory regardless of file size

**NEON (Rule 1)**: ⏸️ Defer
- Profile first to find hotspots
- Field splitting likely not CPU-bound
- Add NEON later if profiling shows benefit

**cloudflare_zlib**: ✅
- `from_bgzip_path()` uses cloudflare_zlib backend
- Already implemented for BAM/SAM

**Network Streaming (Rule 6)**: ✅
- `from_url()` provides HTTP/HTTPS support
- Reuses reqwest infrastructure

### What NOT to Do

❌ **Don't parallelize parsing**
- Conflicts with streaming architecture
- Tab-delimited parsing is fast enough

❌ **Don't use mmap**
- Minimal benefit for text parsing
- `BufReader` is sufficient

❌ **Don't optimize speculatively**
- Implement, profile, then optimize hotspots

---

## Next Week Preview

**Week 2: BED Format** (validates primitives)
- Use `TabDelimitedRecord` trait
- Implement BED3, BED6, BED12 variants
- Test with real BED files
- Benchmark against existing tools

---

## Questions Before Starting?

1. **Module naming**: `src/formats/` or `src/parsers/`?
   - **Recommendation**: `src/formats/` (more general, includes writing)

2. **Network streaming**: Blocking (reqwest::blocking) or async?
   - **Recommendation**: Blocking for now (simpler, matches current BAM code)

3. **Compression**: Auto-detect gzip vs bgzip?
   - **Recommendation**: Separate methods (`from_gzip_path`, `from_bgzip_path`)

4. **Error handling**: Custom error type per format or shared?
   - **Recommendation**: Shared `FormatError` in primitives, format-specific errors as needed

---

**Status**: Ready to begin
**Estimated Duration**: 8-12 hours across multiple sessions
**Deliverable**: Reusable format primitives for BED, GFA, VCF, GFF3
