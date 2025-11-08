# Columnar Alignment Format (CAF) Experiment

**Status**: Design Phase (Pre-implementation)
**Start Date**: After BAM Phase 3-4 (estimated ~Dec 2025)
**Timeline**: 6 weeks
**Type**: Format Innovation Experiment

---

## Overview

CAF is a modern, ARM-native columnar alignment format designed to leverage 2025+ hardware capabilities (NEON SIMD, GPU, Neural Engine) that weren't available when BAM was designed in 2009.

**Key Innovation**: Replace BAM's row-oriented, 4-bit compressed layout with columnar arrays optimized for ARM NEON vectorization.

---

## The Problem with BAM (2009 Design)

BAM was optimized for 2009 constraints:
- **Disk expensive, CPU cheap**: 4-bit sequence encoding (saves storage, burns CPU)
- **gzip compression**: Standard in 2009, now outdated vs zstd/lz4
- **Row-oriented**: Variable-length records limit SIMD vectorization
- **Random access**: Optimized for genome browsers, not analytical queries

**Result**: Modern ARM hardware underutilized, significant CPU overhead for simple operations.

---

## CAF Design Philosophy

Optimize for **2025+ hardware**:
- ✅ ARM NEON: 128-bit SIMD (16 operations/cycle)
- ✅ Modern codecs: zstd (2-3× faster than gzip)
- ✅ Columnar layout: Arrays enable NEON vectorization
- ✅ Pre-decoded sequences: ASCII storage (no 4-bit unpacking overhead)

**Trade-off**: 1.5-2× larger files for **5-10× faster operations**

---

## Expected Performance

### Operation Speedups (Target)

| Operation | BAM (noodles) | CAF (NEON) | Speedup |
|-----------|---------------|------------|---------|
| Parse 100K records | 2.56 sec | 0.25 sec | **10×** |
| Quality filter Q30 | 2.00 sec | 0.08 sec | **25×** |
| Count bases | 1.50 sec | 0.06 sec | **25×** |
| MAPQ > 30 filter | 0.50 sec | 0.03 sec | **16×** |
| Overall | 1.0× | **5-10×** | Target |

### Storage Comparison

| Format | 100K records | Trade-off |
|--------|--------------|-----------|
| BAM (gzip) | 50 MB | Baseline |
| CAF (zstd/lz4) | 75-100 MB | **1.5-2× larger** |

**Acceptable**: Storage is cheap (~$10/TB), CPU time is valuable

---

## Key Design Features

### 1. Columnar Block Layout

**Block size**: 10,000 records (Rule 2 from OPTIMIZATION_RULES.md)

```
Block structure:
  positions:    [i32; 10000]      (zstd compressed)
  mapq:         [u8; 10000]       (raw or RLE)
  flags:        [u16; 10000]      (bitpacked)
  sequences:    [u8; total_len]   (ASCII, lz4 compressed)
  qualities:    [u8; total_len]   (raw)
  cigar_ops:    [u32; cigar_len]  (RLE)
  read_names:   [string]          (dictionary)
  tags:         [nested columnar] (flexible schema)
```

### 2. Pre-decoded Sequences (No 4-bit)

**BAM**: 4-bit encoding (unpack every time)
```
[0x12, 0x48] → decode → "ACGT"  # Overhead on every access
```

**CAF**: ASCII storage (ready to use)
```
[b'A', b'C', b'G', b'T']  # Zero unpacking overhead
```

**Result**: 2× storage, but **eliminates CPU bottleneck**

### 3. NEON Optimization Opportunities

**Quality Filtering** (proven 16-25× in biometal):
```rust
// Process 16 quality scores in parallel
let quals = vld1q_u8(block.qualities.as_ptr());
let mask = vcgeq_u8(quals, threshold);  // 16 compares at once
```

**MAPQ Filtering**:
```rust
// Process 16 MAPQ values in parallel
let mapqs = vld1q_u8(block.mapq.as_ptr());
let mask = vcgtq_u8(mapqs, threshold);  // Parallel comparison
```

**Base Counting** (proven 16-25× in biometal):
```rust
// Count A/C/G/T across all sequences
// NEON processes 16 bases at once
```

### 4. Modern Compression

- **zstd level 3**: Positions, metadata (faster than gzip, better ratio)
- **lz4**: Sequences (extremely fast decompression, >GB/s)
- **Raw**: Quality scores (incompressible, high entropy)

---

## Implementation Phases

### Phase 1: Core Format (2 weeks)
- Columnar block layout
- BAM → CAF conversion (lossless)
- CAF → BAM conversion (lossless)
- zstd/lz4 compression

### Phase 2: NEON Optimization (2 weeks)
- NEON quality filtering
- NEON base counting
- NEON MAPQ filtering
- Benchmark: CAF vs BAM

### Phase 3: Indexing (1 week)
- Block-level index (.caf.idx)
- Region queries (chr:start-end)
- Streaming query API

### Phase 4: Production Polish (1 week)
- 200+ tests
- Documentation
- CLI tools (convert, query, stats)
- Published benchmarks

---

## Use Cases

### ✅ CAF is Great For:
- Batch processing pipelines
- ML training data preparation
- Quality control and filtering
- Large-scale analytical queries
- Research workflows

### ❌ CAF is NOT For:
- Genome browsers (IGV, etc.) - use BAM
- Long-term archival - use BAM
- Clinical pipelines - use BAM
- Sharing with others (initially) - convert to BAM

---

## Success Criteria

**Phase 1-4 Complete**:
- ✅ Lossless BAM ↔ CAF conversion
- ✅ ≥5× speedup for analytical operations
- ✅ 1.5-2× storage overhead (acceptable)
- ✅ Comprehensive tests (200+)
- ✅ Documentation and examples

**Stretch Goals**:
- ⭐ ≥10× speedup (exceptional)
- ⭐ GPU integration (Metal)
- ⭐ Community adoption
- ⭐ Published benchmarks

---

## Files

```
experiments/columnar-alignment-format/
├── README.md                  # This file
├── PROPOSAL.md                # 6-week implementation plan
├── CAF_SPECIFICATION.md       # Format design (v1.0 draft)
├── IMPLEMENTATION_LOG.md      # Development diary (future)
├── BENCHMARKS.md              # Performance results (future)
└── LESSONS_LEARNED.md         # Retrospective (future)
```

---

## Strategic Positioning

**CAF is NOT**:
- A replacement for BAM
- For universal compatibility
- For clinical pipelines

**CAF IS**:
- Optimized for modern ARM hardware
- For analytical performance (5-10×)
- An experiment in format design
- A research tool

**Adoption Strategy**:
- Document trade-offs clearly
- Provide seamless BAM ↔ CAF conversion
- Let users choose (performance vs compatibility)
- "Okay if nobody uses it" - innovation for its own sake :)

---

## Relationship to BAM Implementation

**Sequence**:
1. Complete BAM implementation (Phase 0-6, 9-12 weeks)
2. Start CAF development (Phase 1-4, 6 weeks)

**BAM provides**:
- Format understanding (specifications, edge cases)
- Reference for correctness (differential testing)
- Conversion foundation (BAM ↔ CAF)
- Baseline for benchmarks

**CAF extends**:
- Columnar layout (vs row-oriented)
- Modern compression (zstd vs gzip)
- NEON optimization (extensive vs limited)
- Analytical performance (5-10× faster)

---

## Current Status

**Phase**: Design (Pre-implementation)
**Documents**:
- ✅ PROPOSAL.md (complete)
- ✅ CAF_SPECIFICATION.md (draft v1.0)
- ⏳ Implementation (starts after BAM Phase 3-4)

**Next**: Continue BAM implementation, revisit CAF design in ~6-8 weeks

---

## Questions/Discussion

For questions about CAF design, see:
- `CAF_SPECIFICATION.md` - Technical format details
- `PROPOSAL.md` - Implementation plan and timeline
- `../native-bam-implementation/` - BAM implementation (foundation)

---

**Date**: November 8, 2025
**Status**: Design complete, implementation pending
**Expected Start**: ~December 2025 (after BAM Phase 3-4)
