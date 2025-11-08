# Phase 0 Research: BAM/SAM/CRAM Format Analysis

**Date**: November 8, 2025
**Status**: In Progress
**Duration**: Week 1 (7 days)

---

## Objectives

1. **Format Deep Dive**: Understand BAM/SAM/CRAM specifications thoroughly
2. **noodles Analysis**: Study noodles-bam architecture and implementation
3. **Complexity Mapping**: Identify edge cases and implementation challenges
4. **ARM Opportunities**: Find NEON optimization opportunities

---

## SAM/BAM/CRAM Specification Analysis

### Specification Source

**SAM/BAM v1.6**: https://samtools.github.io/hts-specs/SAMv1.pdf
**Reading Status**: Starting

### Format Overview

#### SAM (Sequence Alignment/Map) Format

**Structure**:
- Header section (optional, starts with '@')
- Alignment section (tab-delimited text)

**Header Lines**:
- `@HD` - Header line (version, sort order)
- `@SQ` - Reference sequence dictionary
- `@RG` - Read group
- `@PG` - Program
- `@CO` - Comment

**Alignment Fields** (11 mandatory):
1. QNAME - Query template name
2. FLAG - Bitwise flags (16 bits)
3. RNAME - Reference sequence name
4. POS - 1-based leftmost mapping position
5. MAPQ - Mapping quality (Phred-scaled)
6. CIGAR - CIGAR string (alignment operations)
7. RNEXT - Reference name of mate/next read
8. PNEXT - Position of mate/next read
9. TLEN - Observed template length
10. SEQ - Segment sequence
11. QUAL - ASCII Phred+33 quality scores

**Optional Fields**: TAG:TYPE:VALUE format
- Types: A (char), i (integer), f (float), Z (string), H (hex), B (array)

#### BAM (Binary Alignment/Map) Format

**Structure**:
```
BAM = BGZF-compressed blocks
Each block contains:
  - Header (text SAM header)
  - Reference sequences (n_ref * ref_data)
  - Alignment records (binary encoded)
```

**Binary Record Structure**:
```
block_size: int32
refID: int32
pos: int32
bin_mq_nl: uint32 (bin:16, mapq:8, l_read_name:8)
flag_nc: uint32 (flag:16, n_cigar_op:16)
l_seq: int32
next_refID: int32
next_pos: int32
tlen: int32
read_name: char[l_read_name]
cigar: uint32[n_cigar_op]
seq: uint8[(l_seq+1)/2] (4-bit encoding: =ACMGRSVTWYHKDBN)
qual: char[l_seq]
tags: (tag_data)*
```

**CIGAR Encoding**:
- Each operation: 4-byte uint32
- Lower 4 bits: operation (M, I, D, N, S, H, P, =, X)
- Upper 28 bits: length

**Sequence Encoding** (4-bit):
```
=:0, A:1, C:2, M:3, G:4, R:5, S:6, V:7,
T:8, W:9, Y:10, H:11, K:12, D:13, B:14, N:15
```

#### CRAM Format

**Overview**:
- Reference-based compression (requires reference genome)
- More complex than BAM (additional compression layers)
- Much higher compression ratio

**Complexity**: HIGH (Phase 5 - optional)

---

## noodles Architecture Analysis

### Repository

**URL**: https://github.com/zaeleus/noodles
**Version**: Latest cloned (Nov 8, 2025)
**Status**: In Progress

### Module Structure Analysis

**noodles-bam** architecture:
```
noodles-bam/src/
├── record/
│   ├── codec/
│   │   ├── decoder/
│   │   │   ├── sequence.rs       ← 4-bit decoding (NEON target!)
│   │   │   ├── cigar/op.rs       ← Bitwise ops (fast already)
│   │   │   ├── quality_scores.rs ← Quality extraction
│   │   │   └── ...
│   │   └── encoder/
│   │       └── ... (mirror structure)
│   ├── sequence.rs
│   └── ...
├── reader.rs                      ← Main BAM reader
├── writer.rs
└── ...
```

**Key observations**:
- Clean codec architecture (separate encode/decode)
- Zero-copy where possible (mutable record reuse)
- Generic over `Read` trait (streaming-friendly)
- BGZF handled separately (noodles-bgzf)

### Implementation Patterns Analysis

#### 1. Sequence Decoding (sequence.rs:36-59)

**Current implementation**:
```rust
pub(super) fn read_sequence(
    src: &mut &[u8],
    sequence: &mut Sequence,
    base_count: usize,
) -> Result<(), DecodeError> {
    let len = base_count.div_ceil(2);
    let (buf, rest) = src.split_at_checked(len)?;
    *src = rest;

    // CRITICAL: Iterator + match-based lookup (line 49-55)
    let bases = buf
        .iter()
        .flat_map(|&b| [decode_base(b >> 4), decode_base(b)]);

    let dst = sequence.as_mut();
    dst.clear();
    dst.extend(bases);
    dst.truncate(base_count);

    Ok(())
}

fn decode_base(n: u8) -> u8 {
    match n & 0x0f {
        0 => b'=', 1 => b'A', 2 => b'C', 3 => b'M',
        4 => b'G', 5 => b'R', 6 => b'S', 7 => b'V',
        8 => b'T', 9 => b'W', 10 => b'Y', 11 => b'H',
        12 => b'K', 13 => b'D', 14 => b'B', 15 => b'N',
        _ => unreachable!(),
    }
}
```

**🎯 PRIMARY NEON OPPORTUNITY IDENTIFIED**:
- Current: Sequential match-based lookup (1 nibble/iteration)
- NEON: Parallel table lookup (`vqtbl1q_u8`)
- Pattern: Process 16 packed bytes → 32 bases in parallel
- Estimated speedup: **16-32× on ARM** (proven pattern from biometal Rule 1)
- **Confidence**: HIGH (this is exactly the pattern that NEON excels at)

#### 2. CIGAR Parsing (cigar/op.rs:36-40)

**Current implementation**:
```rust
pub(crate) fn decode_op(n: u32) -> Result<Op, DecodeError> {
    let kind = decode_kind(n).map_err(DecodeError::InvalidKind)?;
    let len = usize::try_from(n >> 4).map_err(DecodeError::InvalidLength)?;
    Ok(Op::new(kind, len))
}
```

**Analysis**:
- Simple bitwise operations (already 1-2 CPU cycles)
- NEON unlikely to help (no data parallelism)
- Estimated speedup: **1.0×** (no benefit)
- **Confidence**: HIGH (too simple for SIMD)

#### 3. Memory Allocation Patterns

**Record reuse pattern** (from format-integration experiment):
```rust
let mut record = Record::default();  // Allocate once
while reader.read_record(&mut record)? != 0 {
    // Reuse same buffer - constant memory
}
```

**Observation**:
- Excellent zero-copy design
- Mutable record reuse prevents allocation per-record
- biometal will adopt same pattern

#### 4. Error Handling

**Pattern**: Result-based, clear error types
```rust
pub enum DecodeError {
    UnexpectedEof,
    InvalidLength(num::TryFromIntError),
}
```

**Observation**:
- Clean, composable errors
- biometal will use similar approach (BiometalError)

### Performance Characteristics (Preliminary)

**From format-integration experiment**:
- noodles direct: 39.1 Melem/s (baseline)
- Architecture-agnostic (no ARM-specific code detected)

**Expected hotspots** (to be validated by profiling):
1. **Sequence decoding** (decode_base match) - 🎯 NEON target
2. **BGZF decompression** (noodles-bgzf) - biometal already has parallel version (6.5×)
3. **Memory allocation** - minimized via record reuse
4. **CIGAR parsing** - likely negligible (simple bitops)

---

## Complexity Mapping

### Edge Cases to Document

**SAM Parsing**:
- [ ] Malformed CIGAR strings
- [ ] Missing optional fields
- [ ] Invalid FLAG combinations
- [ ] Reference sequence mismatches
- [ ] Quality string length != sequence length
- [ ] Unicode in read names/tags

**BAM Binary Parsing**:
- [ ] Truncated records
- [ ] Invalid BGZF blocks
- [ ] Corrupted compression
- [ ] Integer overflow in positions
- [ ] Invalid sequence encoding (4-bit)
- [ ] Tag type validation

**CRAM** (Phase 5):
- [ ] Reference genome availability
- [ ] Compression dictionary handling
- [ ] Container/slice structure
- [ ] Encoding schemes

---

## ARM NEON Optimization Opportunities

### Hypothesis: ≥3 NEON opportunities with ≥3× potential each

**Potential Areas** (to be validated):

#### 1. Sequence Decoding (4-bit → ASCII)

**Current approach** (likely):
```rust
for byte in packed_seq {
    let base1 = LOOKUP_TABLE[(byte >> 4) as usize];
    let base2 = LOOKUP_TABLE[(byte & 0x0F) as usize];
}
```

**NEON potential**:
- Load 16 packed bytes (32 bases) with `vld1q_u8`
- Parallel table lookup with `vqtbl1q_u8`
- Estimated speedup: **8-16×** (Rule 1 pattern)

**Validation needed**: Profile noodles to confirm this is a hotspot

#### 2. Quality Score Filtering

**Pattern**: Filter records by mean quality ≥ threshold

**NEON potential**:
- We already have quality_filter_neon in biometal
- Integration: Parse BAM → extract quality → NEON filter
- Estimated speedup: **16-25×** (proven in biometal)

**Note**: This is more integration than new NEON work

#### 3. CIGAR String Parsing

**Current approach** (likely):
```rust
for cigar_op in cigar_data {
    let op = cigar_op & 0x0F;
    let len = cigar_op >> 4;
    // Process operation
}
```

**NEON potential**: UNCERTAIN
- CIGAR operations are variable-length, irregular
- May not map well to SIMD (data-dependent branching)
- Estimated speedup: **1.2-2×?** (low confidence)

**Validation needed**: Profile to see if CIGAR parsing is even a bottleneck

#### 4. Base Counting (for statistics)

**Pattern**: Count A/C/G/T/N across alignment

**NEON potential**:
- We already have count_bases_neon in biometal
- Integration: Parse BAM seq → NEON count
- Estimated speedup: **16-25×** (proven)

**Note**: Again, integration not new optimization

#### 5. FLAG Bit Manipulation

**Pattern**: Check multiple FLAG bits (unmapped, duplicate, etc.)

**NEON potential**: LOW
- Bitwise operations on single uint16 are already fast
- SIMD unlikely to help (scalar is 1-2 cycles)
- Estimated speedup: **1.0×** (no benefit)

#### 6. Record Validation

**Pattern**: Validate record consistency (seq length == qual length, etc.)

**NEON potential**: MEDIUM
- Parallel length comparisons
- Parallel character validation
- Estimated speedup: **3-5×?** (uncertain)

**Validation needed**: Profile validation overhead

---

## Profiling Plan

### Setup

**Tools**:
```bash
# macOS Instruments (CPU profiling)
instruments -t "Time Profiler" target/release/noodles_benchmark

# cargo-flamegraph
cargo install cargo-flamegraph
cargo flamegraph --bench noodles_bam_read
```

**Test Data** (to acquire):
- 1000 Genomes BAM (Illumina short reads)
- PacBio HiFi BAM (long reads)
- ONT nanopore BAM (very long reads)
- Different aligners: BWA, minimap2, STAR

### Profiling Experiments

**Experiment 1: Baseline noodles Performance**
```rust
// Measure throughput: records/sec, MB/sec
let mut reader = bam::io::Reader::new(File::open(path)?);
let mut record = Record::default();
let start = Instant::now();
let mut count = 0;

while reader.read_record(&mut record)? != 0 {
    count += 1;
}

let elapsed = start.elapsed();
println!("Throughput: {} records/sec", count as f64 / elapsed.as_secs_f64());
```

**Experiment 2: Hotspot Identification**
- Run with Instruments Time Profiler
- Identify top 10 functions by CPU time
- **Goal**: Find functions consuming ≥10% CPU time

**Experiment 3: ARM vs x86 Comparison**
- Run same benchmark on Mac ARM (M1/M2/M3)
- Run on x86_64 (via Docker or cloud)
- **Hypothesis**: If noodles is architecture-agnostic, ARM should not show significant advantage
- **Expected**: Similar throughput on both (within 10-20%)
- **If ARM is already 2× faster**: NEON optimization less valuable
- **If x86 is faster or equal**: Strong NEON opportunity

**Experiment 4: Memory Allocation Profiling**
```bash
# Use Instruments Allocations template
instruments -t "Allocations" target/release/noodles_benchmark
```
- **Goal**: Identify allocation hotspots
- **Target**: Functions allocating ≥1 MB/sec

---

## Architecture Design (Pending)

### To be designed after research phase:

1. **Module structure** (src/io/bam/, src/io/sam/)
2. **Streaming API** (consistent with FastqStream)
3. **NEON integration points**
4. **Error handling strategy**
5. **Testing approach** (differential vs noodles)

---

## Research Schedule

### Day 1-2: Specification Study (Nov 8-9)
- [ ] Read SAM/BAM v1.6 spec completely (50 pages)
- [ ] Document format structure
- [ ] Identify complexity hotspots
- [ ] Map edge cases

### Day 3-4: noodles Analysis (Nov 10-11)
- [ ] Clone noodles repository
- [ ] Study noodles-bam architecture
- [ ] Trace record parsing flow
- [ ] Document implementation patterns
- [ ] Extract edge case handling

### Day 5-6: ARM Profiling (Nov 12-13)
- [ ] Download test BAM files (diverse)
- [ ] Set up profiling infrastructure
- [ ] Run baseline noodles benchmarks
- [ ] Identify hotspots (Instruments)
- [ ] ARM vs x86 comparison
- [ ] Memory allocation analysis

### Day 7: Synthesis (Nov 14)
- [ ] Consolidate findings
- [ ] Design ARM-first architecture
- [ ] Write NEON_OPPORTUNITIES.md
- [ ] Create DECISION.md (GO/NO-GO)

---

## Success Criteria (Phase 0)

**GO to Phase 1 if**:
- ✅ Identified ≥2 NEON opportunities with ≥3× potential each
- ✅ Profiling shows clear ARM optimization path
- ✅ Architecture design is clean and feasible
- ✅ No major technical blockers discovered

**NO-GO (pivot to noodles) if**:
- ❌ No clear NEON opportunities found
- ❌ noodles already ARM-optimized (unlikely but check)
- ❌ Complexity far exceeds estimates
- ❌ Fundamental architectural issues

---

## Notes and Observations

### [Date: Nov 8] - Day 1 Complete

**Completed**:
- ✅ SAM/BAM v1.6 specification analysis (WebFetch)
- ✅ noodles source code cloned and analyzed
- ✅ Primary NEON opportunity identified (4-bit sequence decoding)
- ✅ NEON_OPPORTUNITIES.md created (comprehensive analysis)
- ✅ Paper planning complete (dual-format strategy documented)

**Initial observations**:
- BAM binary format is complex but well-documented
- Sequence encoding (4-bit) is a **clear NEON target** (high confidence)
- CIGAR parsing may not be SIMD-friendly (data-dependent)
- Quality filtering can leverage existing biometal ops
- Main question: Are these hotspots or negligible overhead?

**Key insight**: Profiling is critical - don't assume bottlenecks, measure them

**Primary NEON Opportunity** (from source analysis):
```rust
// noodles sequence.rs:49-55 - Sequential match-based lookup
fn decode_base(n: u8) -> u8 {
    match n & 0x0f {
        0 => b'=', 1 => b'A', 2 => b'C', 3 => b'M',
        4 => b'G', 5 => b'R', 6 => b'S', 7 => b'V',
        8 => b'T', 9 => b'W', 10 => b'Y', 11 => b'H',
        12 => b'K', 13 => b'D', 14 => b'B', 15 => b'N',
        _ => unreachable!(),
    }
}
```

**Estimated speedup** (to be validated by profiling):
- Operation-level: 3-6× (NEON table lookup vs sequential match)
- Overall BAM parsing: 2-3× (combined with parallel BGZF at 6.5×)

**Next steps** (Days 2-7):
- Days 2-4: Complete specification deep dive, document edge cases
- Days 5-6: Create profiling benchmark, measure actual bottlenecks
- Day 7: GO/NO-GO decision based on profiling data

---

## Links and Resources

### Specifications
- SAM/BAM v1.6: https://samtools.github.io/hts-specs/SAMv1.pdf
- CRAM v3.0: https://samtools.github.io/hts-specs/CRAMv3.pdf

### Reference Implementations
- noodles: https://github.com/zaeleus/noodles
- samtools (C): https://github.com/samtools/samtools
- HTSlib (C): https://github.com/samtools/htslib

### Test Data Sources
- 1000 Genomes: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/
- GIAB (benchmarking): https://www.nist.gov/programs-projects/genome-bottle

---

**Status**: Day 1 - Specification study in progress
**Next**: Complete SAM/BAM v1.6 spec reading, begin documenting format structure
