# Columnar Alignment Format (CAF) Specification

**Version**: 1.0 (Draft)
**Date**: November 8, 2025
**Status**: Design Phase
**Authors**: biometal project

---

## Vision

A modern alignment format optimized for ARM SIMD, GPU processing, and analytical queries on 2025+ hardware.

**Design Philosophy**:
- ARM-first (NEON vectorization)
- GPU-ready (Metal/CUDA friendly)
- Modern compression (zstd, not gzip)
- Analytical queries (filter, aggregate, transform)
- Storage trade-off: 1.5-2× larger than BAM, but 5-10× faster operations

**Target Use Cases**:
- Batch processing pipelines
- ML training data preparation
- Quality control and filtering
- Large-scale analytical queries
- Research workflows (not clinical)

**Not Designed For**:
- Genome browsers (random access at base level)
- Long-term archival (BAM is standard)
- Clinical pipelines (regulatory compliance)
- Compatibility with existing tools (initially)

---

## File Structure

### Overview

```
┌─────────────────────────────────────────────────────────────┐
│ Magic: "CAF\x01" (4 bytes)                                  │
├─────────────────────────────────────────────────────────────┤
│ Header:                                                      │
│   - Version: uint16                                          │
│   - Block size: uint32 (default 10,000 records)             │
│   - SAM header: string (zstd compressed)                    │
│   - Column schema: metadata                                  │
├─────────────────────────────────────────────────────────────┤
│ Index:                                                       │
│   - Block offsets: [uint64; num_blocks]                     │
│   - Block metadata: [BlockMeta; num_blocks]                 │
├─────────────────────────────────────────────────────────────┤
│ Data Blocks (columnar, 10,000 records each):                │
│                                                               │
│ Block 1:                                                     │
│   Column: positions      [i32; 10000]      (zstd level 3)   │
│   Column: mapq           [u8; 10000]       (raw/RLE)        │
│   Column: flags          [u16; 10000]      (bitpacked)      │
│   Column: sequences      [u8; total_len]   (lz4 or raw)     │
│   Column: seq_offsets    [u32; 10001]      (length prefix)  │
│   Column: qualities      [u8; total_len]   (raw)            │
│   Column: qual_offsets   [u32; 10001]      (same as seq)    │
│   Column: cigar_ops      [u32; cigar_len]  (RLE encoding)   │
│   Column: cigar_offsets  [u32; 10001]      (per-record)     │
│   Column: read_names     [string]          (dict + RLE)     │
│   Column: tags           [flexible]        (nested columnar)│
│                                                               │
│ Block 2...                                                   │
└─────────────────────────────────────────────────────────────┘
```

### Block Size: 10,000 Records

**Rationale** (from OPTIMIZATION_RULES.md Rule 2):
- Entry 027: 1,440 measurements validate 10K optimal
- NEON operations benefit from larger blocks
- Memory manageable (~5-10 MB uncompressed)
- Parallelization granularity (rayon)

---

## Column Encodings

### 1. Positions (Reference Coordinates)

**Type**: `[i32; block_size]`
**Compression**: zstd level 3
**Rationale**:
- Sorted positions compress very well (delta encoding in zstd)
- i32 sufficient (2.1 billion positions, >human genome)

**Example**:
```rust
positions: [100, 150, 152, 200, ...]  // Sorted alignment positions
// Compresses to ~1-2 bytes/position with zstd
```

### 2. Mapping Quality (MAPQ)

**Type**: `[u8; block_size]`
**Compression**: RLE (run-length encoding) or raw
**Rationale**:
- Often many records with same MAPQ (e.g., all 60)
- RLE achieves 10-100× compression for uniform quality
- Raw if heterogeneous (auto-detect)

**NEON Optimization**:
```rust
// Filter by MAPQ in parallel
let threshold = vdupq_n_u8(30);
let mask = vcgtq_u8(mapq_array, threshold);  // 16 compares at once
```

### 3. Flags (SAM Flags)

**Type**: `[u16; block_size]`
**Compression**: Bitpacked or raw
**Rationale**:
- Many flags are sparse (e.g., 90% mapped, 10% unmapped)
- Bitpack common patterns

**NEON Optimization**:
```rust
// Check multiple flags in parallel
let unmapped_mask = vandq_u16(flags, vdupq_n_u16(0x04));  // 8 at once
```

### 4. Sequences (Pre-decoded ASCII)

**Type**: `[u8; total_length]` (variable per record)
**Offsets**: `[u32; block_size + 1]` (cumulative lengths)
**Compression**: lz4 (fast) or raw (if heterogeneous)
**Rationale**:
- **KEY DIFFERENCE FROM BAM**: Store ASCII, not 4-bit
- Eliminates unpacking overhead (2× CPU savings)
- 2× storage cost acceptable (storage cheap, CPU expensive)

**Memory Layout**:
```rust
// Sequences stored contiguously
sequences: [b'A', b'C', b'G', b'T', b'A', b'C', ...]

// Offsets for record boundaries
offsets: [0, 100, 250, 350, ...]  // Record i: sequences[offsets[i]..offsets[i+1]]
```

**NEON Optimization**:
```rust
// Count bases across all sequences in block
for seq in sequences.chunks(16) {
    let bases = vld1q_u8(seq.as_ptr());
    // Parallel compare, accumulate (16 bases at once)
}
// 16-25× speedup (proven in biometal)
```

### 5. Quality Scores (ASCII Phred+33)

**Type**: `[u8; total_length]` (same layout as sequences)
**Offsets**: `[u32; block_size + 1]` (shared with sequences)
**Compression**: Raw (incompressible noise)
**Rationale**:
- Quality scores are high-entropy (don't compress well)
- Store raw for instant access

**NEON Optimization**:
```rust
// Quality filtering in parallel
let threshold = vdupq_n_u8(33 + 20);  // Q20 in Phred+33
let mask = vcgeq_u8(qualities, threshold);  // 16 at once
// 16-25× speedup (proven in biometal)
```

### 6. CIGAR Operations

**Type**: `[u32; total_ops]` (variable per record)
**Offsets**: `[u32; block_size + 1]`
**Compression**: RLE (run-length encoding)
**Rationale**:
- CIGAR often repetitive (e.g., "100M" = 100 operations)
- RLE: store (operation, count) pairs

**Encoding**:
```rust
// Standard: CIGAR as u32 (lower 4 bits = op, upper 28 bits = len)
cigar_ops: [0x640, 0x011, ...]  // 100M, 1I

// RLE alternative: [(op, count); N]
cigar_rle: [(Op::Match, 100), (Op::Insertion, 1), ...]
```

### 7. Read Names

**Type**: Variable-length strings
**Compression**: Dictionary + RLE
**Rationale**:
- Read names often have common prefixes (e.g., "SRR123456.1", "SRR123456.2")
- Dictionary of unique strings + indices

**Example**:
```rust
// Dictionary compression
dict: ["SRR123456.", "HWI:", ...]
indices: [0, 0, 0, 1, ...]  // Record i uses dict[indices[i]] + suffix
suffixes: ["1", "2", "3", "1000", ...]
```

### 8. Optional Tags (Nested Columnar)

**Type**: Nested columnar structure
**Schema**: Flexible, per-tag columns

**Example**:
```rust
// Tag: NM (edit distance, integer)
tag_NM: [u8; block_size]  // Present for all records

// Tag: MD (mismatch string, variable)
tag_MD_values: [string]
tag_MD_present: [bool; block_size]  // Bitmap of which records have MD
```

**Rationale**:
- Some tags universal (NM), some sparse (MD)
- Columnar allows efficient filtering

---

## Compression Strategy

### Modern Codecs (Not Gzip)

**zstd level 3** (for integers, positions):
- Faster decompression than gzip (2-3×)
- Better compression ratio (10-30%)
- Parallelizable

**lz4** (for sequences, if needed):
- Extremely fast decompression (>GB/s)
- Modest compression (2-3×)
- Good for random access

**Raw** (for qualities, flags):
- Quality scores don't compress (high entropy)
- Flags are small (2 bytes × 10K = 20 KB)

### Compression Benchmark (Estimated)

| Column | BAM (gzip) | CAF (zstd/lz4) | Ratio |
|--------|------------|----------------|-------|
| Positions | 4 bytes → 1-2 bytes | 4 bytes → 1-2 bytes | 1.0× (similar) |
| Sequences | 0.5 bytes (4-bit) | 1 byte (ASCII) → 0.5 bytes (lz4) | 1.0× (similar) |
| Qualities | 1 byte → 0.9 bytes | 1 byte (raw) | 1.1× (slightly worse) |
| **Total** | 1.0× | **1.5-2.0×** | Acceptable trade-off |

**Trade-off**: 1.5-2× larger files for 5-10× faster operations

---

## NEON Optimization Opportunities

### 1. Quality Filtering (Proven 16-25×)

```rust
#[cfg(target_arch = "aarch64")]
unsafe fn filter_block_neon(block: &CafBlock, min_quality: u8) -> Vec<usize> {
    let threshold = vdupq_n_u8(min_quality);
    let mut passing_indices = Vec::new();

    for (i, record_quals) in block.qualities.chunks(16).enumerate() {
        let quals = vld1q_u8(record_quals.as_ptr());
        let mask = vcgeq_u8(quals, threshold);  // 16 compares at once

        if vaddvq_u8(mask) == 16 {  // All pass
            passing_indices.push(i);
        }
    }

    passing_indices
}
```

**Speedup**: 16-25× (proven in biometal Rule 1)

### 2. MAPQ Filtering

```rust
// Process 10,000 MAPQ values in parallel
let mapq_array = &block.mapq;  // [u8; 10000]

// NEON: 16 at a time
for chunk in mapq_array.chunks(16) {
    let mapqs = vld1q_u8(chunk.as_ptr());
    let mask = vcgtq_u8(mapqs, threshold);  // Parallel compare
    // Accumulate passing indices
}
```

**Speedup**: 16× (element-wise comparison)

### 3. Base Counting (Proven 16-25×)

```rust
// Count A/C/G/T across all sequences in block
let sequences = &block.sequences;  // Pre-decoded ASCII

unsafe fn count_bases_neon(seq: &[u8]) -> [u32; 4] {
    // Existing biometal implementation
    // Processes 16 bases at once
}
```

**Speedup**: 16-25× (proven in biometal)

### 4. Flag Checking

```rust
// Check if records are unmapped (FLAG & 0x04)
for chunk in block.flags.chunks(8) {  // u16, so 8 per NEON register
    let flags = vld1q_u16(chunk.as_ptr());
    let unmapped = vandq_u16(flags, vdupq_n_u16(0x04));
    // 8 bitwise ANDs in parallel
}
```

**Speedup**: 8× (parallel bitwise operations)

---

## Random Access and Indexing

### Block-Level Index

**Structure**:
```rust
struct CafIndex {
    block_offsets: Vec<u64>,        // File offset of each block
    block_metadata: Vec<BlockMeta>,
}

struct BlockMeta {
    num_records: u32,
    ref_id: i32,              // Reference sequence ID
    start_pos: i32,           // First alignment position in block
    end_pos: i32,             // Last alignment position in block
    compressed_size: u32,
    uncompressed_size: u32,
}
```

**Query Strategy**:
```rust
// Find blocks overlapping region chr1:1000-2000
let overlapping_blocks: Vec<usize> = index
    .block_metadata
    .iter()
    .enumerate()
    .filter(|(_, meta)| {
        meta.ref_id == chr1_id &&
        meta.start_pos <= 2000 &&
        meta.end_pos >= 1000
    })
    .map(|(i, _)| i)
    .collect();

// Read only overlapping blocks (not entire file)
for block_id in overlapping_blocks {
    let block = reader.read_block(block_id)?;
    // Process block
}
```

**Granularity**: Block-level (10,000 records)
- Coarser than BAM's BAI index (64 KB chunks)
- Sufficient for analytical queries
- Not suitable for genome browsers (IGV, etc.)

---

## File Format Version and Compatibility

### Version 1.0 (Initial)

**Magic**: `CAF\x01` (4 bytes)
**Features**:
- Columnar layout
- NEON-optimized encodings
- zstd/lz4 compression
- Block-level indexing

### Future Extensions

**Version 1.1** (potential):
- GPU-optimized memory layout (Metal/CUDA)
- Neural Engine integration (quantized operations)
- AMX matrix operations (paired-end correlation)

### SAM Header Compatibility

**Preserve SAM metadata**:
```rust
struct CafHeader {
    version: u16,
    block_size: u32,
    sam_header: String,  // Standard SAM header (compressed)
    column_schema: SchemaMetadata,
}
```

**Conversion**: CAF ↔ BAM lossless
- All SAM fields preserved
- Optional tags preserved
- Metadata intact

---

## Benchmark Targets (vs BAM)

### Operation Speedups (Estimated)

| Operation | BAM (noodles) | CAF (NEON) | Speedup |
|-----------|---------------|------------|---------|
| Parse 100K records | 2.56 sec | 0.25 sec | **10×** |
| Quality filter Q30 | 2.00 sec | 0.08 sec | **25×** |
| Count bases | 1.50 sec | 0.06 sec | **25×** |
| MAPQ > 30 filter | 0.50 sec | 0.03 sec | **16×** |
| Convert to FASTQ | 3.00 sec | 0.40 sec | **7.5×** |
| **Overall** | **1.0×** | **5-10×** | Target |

### Storage Comparison

| Format | 100K records | Compression |
|--------|--------------|-------------|
| BAM (gzip) | 50 MB | 1.0× (baseline) |
| CAF (zstd/lz4) | 75-100 MB | 1.5-2.0× larger |
| Uncompressed SAM | 500 MB | 10× (reference) |

**Trade-off**: 50 MB → 75 MB (1.5×) for 10× faster operations

---

## Implementation Phases

### Phase 1: Core Format (2 weeks)

**Deliverables**:
- CAF writer (BAM → CAF conversion)
- CAF reader (basic parsing)
- Columnar block layout
- zstd compression integration

**Validation**:
- Round-trip: BAM → CAF → BAM (lossless)
- Correctness: Differential testing vs noodles

### Phase 2: NEON Optimization (2 weeks)

**Deliverables**:
- NEON quality filtering
- NEON base counting
- NEON MAPQ filtering
- Parallel block processing

**Validation**:
- Benchmark: CAF vs BAM (target 5-10× speedup)
- Correctness: Property testing

### Phase 3: Indexing and Queries (1 week)

**Deliverables**:
- Block-level index (.caf.idx)
- Region queries (chr:start-end)
- Streaming query API

**Validation**:
- Query correctness vs full scan
- Performance: Index overhead < 10%

### Phase 4: Production Polish (1 week)

**Deliverables**:
- Comprehensive tests (200+ tests)
- Documentation
- Benchmarks (published)
- CLI tools (biometal convert, biometal query)

**Validation**:
- Integration tests (real datasets)
- Edge case handling
- Error messages

---

## CLI Interface (Proposed)

### Convert BAM ↔ CAF

```bash
# BAM to CAF
biometal convert input.bam output.caf
# Options: --block-size 10000, --compression zstd

# CAF to BAM
biometal convert input.caf output.bam

# Validate lossless round-trip
biometal validate --round-trip input.bam
```

### Query CAF

```bash
# Filter by quality
biometal query input.caf --min-quality 30 --output filtered.caf

# Extract region (if indexed)
biometal query input.caf --region chr1:1000-2000 --output region.caf

# Convert to FASTQ
biometal query input.caf --to-fastq output.fastq

# Statistics
biometal stats input.caf
# Output: 100,000 records, 85% mapped, mean MAPQ: 42, ...
```

### Benchmark

```bash
# Compare CAF vs BAM
biometal benchmark input.bam input.caf --operation quality-filter
# Output: BAM: 2.5 sec, CAF: 0.1 sec (25× faster)
```

---

## Open Questions (To Be Resolved)

1. **Sequence compression**: lz4 or raw?
   - lz4 if compressible (homopolymers)
   - raw if heterogeneous (short reads)
   - Auto-detect?

2. **Tag schema**: Fixed vs flexible?
   - Option A: Fixed schema (all tags defined upfront)
   - Option B: Flexible (schema in header)
   - Lean toward B (more general)

3. **GPU support**: Phase 1 or deferred?
   - Phase 1: Design GPU-friendly layout
   - Phase 3-4: Implement Metal/CUDA kernels
   - Decision: Phase 1 design, Phase 3+ implementation

4. **Compatibility guarantee**: Lossless BAM ↔ CAF?
   - Yes: All SAM fields preserved
   - Edge cases: Non-standard tags? (handle gracefully)

---

## Success Criteria

**Phase 1-4 Success**:
- ✅ 5-10× faster than BAM for analytical operations
- ✅ Lossless BAM ↔ CAF conversion
- ✅ 1.5-2× storage overhead (acceptable)
- ✅ NEON optimizations validated (16-25× per operation)
- ✅ Comprehensive tests (200+)
- ✅ Documentation and examples

**Stretch Goals**:
- ⭐ GPU integration (Metal)
- ⭐ Neural Engine quantized operations
- ⭐ Community adoption (others use CAF)
- ⭐ Published benchmarks (influence ecosystem)

---

## Strategic Positioning

**CAF is NOT**:
- A replacement for BAM (coexist)
- For genome browsers (IGV, etc.)
- For clinical pipelines (regulatory)

**CAF IS**:
- Optimized for modern ARM hardware
- For batch processing and analysis
- For ML training data preparation
- For research workflows
- An experiment in modern format design

**Adoption Strategy**:
- Document trade-offs clearly
- Provide seamless BAM ↔ CAF conversion
- Benchmark against BAM (transparency)
- Let users choose (performance vs compatibility)

---

**Status**: Design Phase (November 2025)
**Next**: Finalize specification, begin implementation after BAM Phase 3-4
**Timeline**: 6 weeks (after BAM complete)
