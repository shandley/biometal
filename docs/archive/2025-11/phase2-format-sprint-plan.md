# Phase 2: Format Library Sprint

**Date**: November 13, 2025
**Decision**: Build comprehensive bioinformatics format library with streaming-first design
**Duration**: 12-16 weeks
**Status**: 🚀 APPROVED - Ready to begin

---

## Strategic Context

### What Changed

**Original Plan**: Phase 2 = Rules 3+4 implementation for 16× speedup
**Reality Check**: Rules 3+4 not viable (see RULES_3_4_REALITY_CHECK.md)
**New Understanding**: biometal is **already fast** - the value is in building comprehensive, robust primitives

### Core Insight

> "Part of what I wanted to build was a robust rust library and speed and other optimizations have always been secondary (although exciting and fun!). I really just want my own software library and primitives to build my own software. Full stack from bottom up."

**What We've Already Built** (v1.7.0):
- ✅ Streaming architecture (constant memory, terabyte-scale)
- ✅ Competitive performance (92 MiB/s BAM, 16-25× NEON)
- ✅ Production quality (403 tests, 2 audit reports)
- ✅ Evidence-based design (ASBB-validated optimizations)
- ✅ Formats: FASTQ, FASTA, BAM, SAM, BAI

**What's Missing**:
- Comprehensive format coverage (VCF, BED, GFF, GFA, etc.)
- User adoption and real-world validation
- Your own projects built on biometal primitives

### The New Direction

**Focus**: Build a comprehensive bioinformatics format library
- Streaming-first design for all formats
- Apply proven optimizations (NEON, cloudflare_zlib)
- Production quality (tests, docs, benchmarks)
- Don't chase speculative optimizations

**Goal**: Full-stack bioinformatics primitives library in Rust

---

## Format Strategy

### Tier 1: High-Value Formats (Phase 2 Focus)

**Week 1-2: Format Primitives** (Foundation)
- Generic tab-delimited streaming parser
- Header/comment handling
- Field validation framework
- Record trait abstraction

**Week 3-4: BED Format** (Simplest)
- BED3, BED6, BED12 variants
- Interval operations
- Validates primitives

**Week 5-7: GFA Format** (Graph Assembly)
- S/L/P/W/C/H line types
- Path extraction/manipulation
- Assembly graph analysis/statistics
- General-purpose (not limited to single use case)

**Week 8-10: VCF Format** (Variants)
- Header parsing (##INFO, ##FORMAT, etc.)
- Record parsing (REF, ALT, QUAL, INFO, genotypes)
- High-impact for variant calling workflows

**Week 11-12: GFF3 Format** (Annotation)
- 9-column format
- Attributes parsing
- Gene/transcript features

**Week 13-14: Utilities**
- FAI index (FASTA random access)
- TBI index (Tabix for VCF/BED/GFF)
- Format conversion helpers

**Week 15-16: Integration**
- Cross-format examples
- Python bindings for new formats
- Documentation
- v2.0.0 release

### Tier 2: Future Formats (Post-Phase 2)

- **PAF**: Minimap2 pairwise alignment format
- **GTF**: Gene transfer format (similar to GFF)
- **CRAM**: Reference-based compression (complex)
- **BCF**: Binary VCF
- **narrowPeak/broadPeak**: ChIP-seq peaks

### Tier 3: Specialized Formats (Community-Driven)

- **BigWig/BigBed**: Visualization
- **MAF**: Multiple alignment
- **rGFA**: Reference GFA
- **HAL**: Hierarchical alignment (cactus)

---

## Format Development Framework

### For Each Format, Apply Proven Optimizations

#### 1. Streaming Architecture (Rule 5) ✅

**Always**:
- Iterator-based API (`impl Iterator<Item = Result<Record>>`)
- Constant memory (no `Vec<Record>` accumulation)
- Lazy parsing (parse on demand)
- Buffered I/O (`BufReader`)

**Example**:
```rust
pub struct BedParser<R: Read> {
    reader: BufReader<R>,
    line_buf: String,  // Reused buffer
}

impl<R: Read> Iterator for BedParser<R> {
    type Item = Result<BedRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line_buf.clear();
        match self.reader.read_line(&mut self.line_buf) {
            Ok(0) => None,
            Ok(_) => Some(self.parse_line()),
            Err(e) => Some(Err(e.into())),
        }
    }
}
```

#### 2. ARM NEON Optimization (Rule 1) ✅

**When Applicable**:
- String operations (field extraction, splitting)
- Numeric parsing (quality scores, coverage values)
- Pattern matching (CIGAR parsing, attribute parsing)

**Decision Framework**:
- Profile first (is this a hotspot?)
- Check if NEON helps (sequence operations > text parsing)
- Implement scalar fallback always

**Example** (GFA sequence parsing):
```rust
// GFA S line: S<tab>segment_name<tab>ACGT...<tab>...
#[cfg(target_arch = "aarch64")]
fn parse_sequence_neon(seq: &[u8]) -> Result<SegmentSequence> {
    // NEON base counting, validation
}

#[cfg(not(target_arch = "aarch64"))]
fn parse_sequence_scalar(seq: &[u8]) -> Result<SegmentSequence> {
    // Portable fallback
}
```

#### 3. cloudflare_zlib (Rule 1 variant) ✅

**When Applicable**:
- BGZF-compressed formats (VCF.gz, BED.gz via bgzip)
- Compressed text formats

**Pattern**:
```rust
use flate2::read::MultiGzDecoder;

pub fn from_bgzip_path(path: impl AsRef<Path>) -> Result<Self> {
    let file = File::open(path)?;
    let decoder = MultiGzDecoder::new(file);
    Self::from_reader(BufReader::new(decoder))
}
```

#### 4. Network Streaming (Rule 6) ✅

**Always Provide**:
- HTTP/HTTPS support
- SRA support (if applicable)
- Range request support (with TBI/CSI indices)

**Pattern**:
```rust
pub fn from_url(url: &str) -> Result<Self> {
    let response = reqwest::blocking::get(url)?;
    Self::from_reader(BufReader::new(response))
}
```

#### 5. Indexing (When Needed)

**For Random Access Formats**:
- TBI (Tabix): VCF, BED, GFF
- CSI: Large genomes
- FAI: FASTA

**Pattern** (learned from BAI):
```rust
pub struct BedIndex {
    // Coordinate-based indexing
}

impl BedParser<R> {
    pub fn query(&mut self, region: &str) -> Result<impl Iterator<Item = Result<BedRecord>>> {
        // Use index for fast seeking
    }
}
```

### What NOT to Do (Avoid Optimization Trap)

❌ **Don't implement optimizations without profiling**
- Profile first, optimize hotspots
- Most format parsing is I/O-bound, not CPU-bound

❌ **Don't chase parallel parsing (Rule 3)**
- Conflicts with streaming architecture
- Text parsing is fast enough

❌ **Don't use mmap speculatively (Rule 4)**
- Only ~1% benefit for compressed files
- Standard BufReader is fine

❌ **Don't batch-accumulate records**
- Breaks streaming promise
- Use iterators, not `Vec<Record>`

✅ **Do apply proven optimizations**:
- Streaming (Rule 5)
- NEON for hotspots (Rule 1)
- cloudflare_zlib for compression
- Network streaming (Rule 6)

✅ **Do validate with property tests**:
- Round-trip tests (parse → write → parse)
- NEON == scalar correctness
- Edge case handling

✅ **Do benchmark (but don't obsess)**:
- N=30 statistical rigor
- Compare to existing tools (if available)
- Document performance characteristics

---

## GFA Format Specifics

### Why GFA Is Important

**Use Cases**:
1. **Pangenome graphs**: T2T, 1000 Genomes graph reference
2. **Assembly graphs**: Flye, Canu, hifiasm output
3. **Variation graphs**: vg toolkit
4. **Read overlap graphs**: Long-read assembly

**Current Rust Ecosystem**:
- `gfa` crate: Basic, unmaintained (last update 2020)
- `handlegraph` crate: Trait-based, complex
- **Gap**: No streaming-first, production-quality GFA parser

### GFA Format Structure

**Line Types** (tab-delimited):
```
H  Header    H  VN:Z:1.0  ...
S  Segment   S  seg1  ACGT  LN:i:4  ...
L  Link      L  seg1  +  seg2  -  4M  ...
P  Path      P  path1  seg1+,seg2-  4M,5M  ...
W  Walk      W  sample1  0  chrom1  100  200  >seg1<seg2  ...
C  Containment  (less common)
```

**Your Project Needs**:
- Path extraction/manipulation (P/W lines)
- Assembly graph analysis/statistics (S/L lines, graph properties)

**General GFA Library Should Support**:
1. **Streaming parsing**: Process line-by-line (constant memory)
2. **Line type filtering**: Iterator adapters for S/L/P/W
3. **Graph construction**: Optional in-memory graph from stream
4. **Path operations**: Extract, filter, traverse paths
5. **Statistics**: Segment count, N50, connectivity, etc.
6. **Writing**: Serialize back to GFA format

### GFA API Design (Proposal)

```rust
// Streaming API (constant memory)
pub struct GfaParser<R: Read> {
    reader: BufReader<R>,
}

pub enum GfaLine {
    Header(Header),
    Segment(Segment),
    Link(Link),
    Path(Path),
    Walk(Walk),
    Containment(Containment),
}

pub struct Segment {
    pub name: String,
    pub sequence: Option<Vec<u8>>,  // Optional (can be '*')
    pub length: Option<usize>,      // LN tag
    pub tags: HashMap<String, Tag>, // Optional tags
}

pub struct Link {
    pub from: String,
    pub from_orient: Orientation,
    pub to: String,
    pub to_orient: Orientation,
    pub overlap: String,            // CIGAR string
}

pub struct Path {
    pub name: String,
    pub segments: Vec<(String, Orientation)>,
    pub overlaps: Vec<String>,      // CIGARs
}

impl<R: Read> Iterator for GfaParser<R> {
    type Item = Result<GfaLine>;
    // Streaming line-by-line
}

// Filtered iterators
impl<R: Read> GfaParser<R> {
    pub fn segments(self) -> impl Iterator<Item = Result<Segment>> {
        self.filter_map(|line| match line {
            Ok(GfaLine::Segment(s)) => Some(Ok(s)),
            Err(e) => Some(Err(e)),
            _ => None,
        })
    }

    pub fn paths(self) -> impl Iterator<Item = Result<Path>> {
        // Filter for paths
    }
}

// Optional in-memory graph (for your project)
pub struct GfaGraph {
    segments: HashMap<String, Segment>,
    links: Vec<Link>,
    paths: HashMap<String, Path>,
}

impl GfaGraph {
    pub fn from_parser<R: Read>(parser: GfaParser<R>) -> Result<Self> {
        // Build graph in memory (for analysis)
    }

    // Graph analysis methods
    pub fn segment_count(&self) -> usize;
    pub fn n50(&self) -> usize;
    pub fn connected_components(&self) -> Vec<Vec<String>>;
    pub fn extract_path(&self, path_name: &str) -> Option<Vec<u8>>;
}

// Writing
pub struct GfaWriter<W: Write> {
    writer: BufWriter<W>,
}

impl<W: Write> GfaWriter<W> {
    pub fn write_header(&mut self, header: &Header) -> Result<()>;
    pub fn write_segment(&mut self, segment: &Segment) -> Result<()>;
    pub fn write_path(&mut self, path: &Path) -> Result<()>;
}
```

**Design Philosophy**:
- **Streaming-first**: Low-level API is iterator-based
- **Opt-in graph**: Build in-memory graph when needed
- **Flexible**: Filter to specific line types
- **Your project**: Can use both streaming (large graphs) and in-memory (analysis)
- **General-purpose**: Supports all GFA use cases (assembly, pangenome, variation)

---

## Testing Strategy

### For Each Format

**1. Property-Based Tests** (proptest)
```rust
proptest! {
    #[test]
    fn test_bed_round_trip(
        chrom in "[a-z]{3,10}",
        start in 0u64..1000000,
        end in 0u64..1000000,
    ) {
        let start = start.min(end);
        let end = end.max(start + 1);

        let record = BedRecord { chrom, start, end, ..Default::default() };
        let serialized = record.to_string();
        let parsed = BedRecord::from_str(&serialized)?;

        prop_assert_eq!(record, parsed);
    }
}
```

**2. Real-World Data Tests**
- Download example files from public databases
- Parse and validate (e.g., GENCODE GTF, 1000 Genomes VCF)
- Compare against reference implementations

**3. Edge Case Tests**
- Empty files
- Malformed records
- Missing fields
- Comment-only files
- Unicode in headers

**4. Performance Benchmarks** (N=30)
- Large file parsing (1-10 GB files)
- Memory usage validation (constant ~5 MB)
- Comparison to existing tools (if available)

**5. NEON Correctness** (if applicable)
```rust
#[test]
fn test_neon_equals_scalar() {
    let test_data = generate_test_sequences();

    for seq in test_data {
        let neon_result = parse_sequence_neon(seq);
        let scalar_result = parse_sequence_scalar(seq);
        assert_eq!(neon_result, scalar_result);
    }
}
```

---

## Documentation Requirements

### For Each Format

**1. API Documentation**
- Doctests for every public function
- Usage examples (basic, intermediate, advanced)
- Performance characteristics

**2. Format Guide** (docs/)
- Format specification reference
- Common operations (examples)
- Performance tips
- Migration from other libraries

**3. Tutorial Notebooks** (notebooks/)
- Getting started
- Real-world workflows
- Integration with other formats

**4. Benchmarks** (benchmarks/)
- Parsing performance
- Memory usage
- Comparison to existing tools

---

## Success Criteria

### Per-Format Checklist

**Implementation**:
- [ ] Streaming parser (`impl Iterator<Item = Result<Record>>`)
- [ ] Writer (if write support desired)
- [ ] Index support (if applicable: TBI, CSI, FAI)
- [ ] Network streaming (HTTP/SRA)
- [ ] Compression support (gzip/bgzip)
- [ ] NEON optimization for hotspots (if profiling shows benefit)

**Testing**:
- [ ] Property-based tests (round-trip, edge cases)
- [ ] Real-world data validation
- [ ] Edge case coverage (malformed, empty, unicode)
- [ ] NEON == scalar correctness (if applicable)
- [ ] Benchmarks (N=30, memory validation)

**Documentation**:
- [ ] API docs with examples
- [ ] Format guide (docs/FORMATNAME.md)
- [ ] Tutorial notebook (notebooks/NN_formatname.ipynb)
- [ ] Performance benchmarks documented

**Quality**:
- [ ] No panics in library code (Result everywhere)
- [ ] Cross-platform validated (Mac ARM, Linux ARM, x86_64)
- [ ] Python bindings (PyO3)
- [ ] Zero clippy warnings

### Overall Phase 2 Success

**Formats Implemented**: BED, GFA, VCF, GFF3, FAI, TBI
**Timeline**: 12-16 weeks
**Tests**: 500+ total (including new formats)
**Documentation**: 60,000+ words (add 20k to existing 40k)
**Release**: v2.0.0 - Comprehensive format library

---

## First Steps (This Session)

### 1. Format Primitives Design (Week 1-2 Foundation)

**Immediate Tasks**:
1. ✅ Document format sprint plan (this file)
2. Design generic tab-delimited parser
3. Design header/comment handling
4. Design field validation framework
5. Create shared types (GenomicInterval, Strand, etc.)

**Deliverables**:
- `src/formats/mod.rs` - Format primitives module
- `src/formats/primitives/` - Shared infrastructure
  - `tab_delimited.rs` - Generic parser
  - `header.rs` - Header/comment handling
  - `fields.rs` - Field validation
  - `genomic.rs` - Shared genomic types
- Tests for primitives

### 2. BED Format (Week 3-4)

**After primitives are solid**:
- `src/formats/bed/` - BED implementation
- Tests with real BED files
- Documentation
- Python bindings

### 3. GFA Format (Week 5-7)

**Your project focus**:
- `src/formats/gfa/` - GFA implementation
- Streaming + optional in-memory graph
- Path extraction/manipulation
- Assembly statistics
- Tests, docs, Python bindings

---

## Optimization Decision Framework

### When Considering an Optimization

**Ask these questions**:

1. **Is this a proven optimization from OPTIMIZATION_RULES.md?**
   - ✅ Streaming (Rule 5): Always apply
   - ✅ NEON (Rule 1): Apply to hotspots after profiling
   - ✅ cloudflare_zlib: Apply to compressed formats
   - ✅ Network streaming (Rule 6): Always provide
   - ❌ Parallel parsing (Rule 3): Conflicts with streaming
   - ❌ mmap (Rule 4): Minimal benefit, don't use

2. **Have I profiled and found this is a hotspot?**
   - ✅ Profile shows >10% time here: Consider optimization
   - ❌ No profile data: Don't optimize yet

3. **Does this break streaming architecture?**
   - ✅ Maintains constant memory: OK to implement
   - ❌ Requires accumulating records: Don't implement

4. **Is there a scalar fallback?**
   - ✅ ARM + x86_64 both work: OK to implement
   - ❌ ARM-only: Add scalar fallback first

**If uncertain**: Default to streaming + standard library. Optimize later if profiling shows need.

---

## Communication and Planning

### Weekly Updates

**Each week, document**:
1. Format(s) worked on
2. Progress (tests passing, docs written)
3. Performance results (if benchmarked)
4. Decisions made
5. Next week's plan

**File**: `PHASE2_WEEKLY_UPDATES.md`

### Format-Specific Decisions

**When facing design choices**, document:
- Question/choice
- Options considered
- Decision made
- Rationale
- Reference to evidence (if applicable)

**File per format**: `docs/decisions/FORMATNAME_DECISIONS.md`

---

## Project Integration

### Your GFA Project

**Phase 2 Timeline Alignment**:
- Week 1-2: Format primitives (foundation)
- Week 3-4: BED (validation)
- **Week 5-7: GFA (your project can start)**

**By Week 7**:
- biometal GFA parser complete
- Your project can use: streaming parsing, path extraction, graph analysis
- General-purpose enough for others' pangenome/assembly use cases

**Your Project Benefits**:
- Streaming-first (handle huge assembly graphs)
- Path manipulation (extract/filter paths)
- Graph statistics (N50, connectivity, etc.)
- Production quality (tested, documented)
- Python bindings (if needed)

---

## Summary

**Decision**: Build comprehensive format library (BED, GFA, VCF, GFF3, utilities)
**Duration**: 12-16 weeks
**Philosophy**:
- Apply proven optimizations (streaming, NEON, compression)
- Don't chase speculative optimizations
- Production quality (tests, docs, benchmarks)
- Support your projects + general community use

**Next Steps**:
1. Design format primitives (tab-delimited, headers, fields)
2. Implement BED (validates primitives)
3. Implement GFA (your project + general-purpose)
4. Continue with VCF, GFF3, utilities

**Goal**: v2.0.0 - Full-stack bioinformatics primitives library in Rust

---

**Status**: 🚀 Ready to begin
**First Task**: Design format primitives (this session)
**Expected Completion**: March 2026 (16 weeks from Nov 13, 2025)
