# BAM Phase 0 - Day 1 Summary

**Date**: November 8, 2025
**Status**: Day 1 Complete ✅
**Next**: Continue with profiling (Days 2-7)

---

## What We Accomplished Today

### 1. Dual-Format Strategy Established ✅

**Decision**: Pursue both BAM (compatibility) + CAF (performance)
- BAM: ARM-optimized, spec-compliant (2-3× speedup)
- CAF: Columnar, NEON-first (5-10× speedup)
- Publication goal: "Dual-Format Strategy for ARM-Native Bioinformatics"

### 2. Complete Paper Planning ✅

**Created**:
- `experiments/alignment-formats-paper/PAPER_OUTLINE.md` (15 KB)
- `experiments/alignment-formats-paper/BENCHMARKING_PLAN.md` (18 KB)
- `experiments/alignment-formats-paper/README.md` (9 KB)

**Timeline**: 24 weeks to submission (Q2 2026)
- Weeks 1-12: BAM implementation
- Weeks 13-18: CAF implementation
- Weeks 19-20: Comprehensive benchmarking
- Weeks 21-24: Paper writing & submission

### 3. BAM Phase 0 Research (Day 1) ✅

**Documents Created**:
- `PROPOSAL.md` (14 KB) - 9-12 week implementation plan
- `PHASE_0_RESEARCH.md` (14 KB) - Format & noodles analysis
- `NEON_OPPORTUNITIES.md` (10 KB) - Optimization targets

**Key Findings**:

**A. Format Analysis** (SAM/BAM v1.6):
```
BAM Structure:
- Magic: "BAM\1"
- Header: SAM header (text, zstd compressed)
- References: Dictionary of sequences
- Alignments: BGZF-compressed records

Record encoding:
- Positions: int32
- MAPQ: uint8
- FLAGS: uint16
- Sequences: 4-bit encoding (2 bases/byte) ← NEON TARGET
- Quality: uint8 array (Phred+33)
- CIGAR: uint32 array (op + length)
```

**B. noodles Source Analysis**:

Primary NEON opportunity identified:
```rust
// File: noodles/noodles-bam/src/record/codec/decoder/sequence.rs:49-81

// Current: Sequential match-based 4-bit decoding
fn decode_base(n: u8) -> u8 {
    match n & 0x0f {
        0 => b'=', 1 => b'A', 2 => b'C', 3 => b'M',
        4 => b'G', 5 => b'R', 6 => b'S', 7 => b'V',
        8 => b'T', 9 => b'W', 10 => b'Y', 11 => b'H',
        12 => b'K', 13 => b'D', 14 => b'B', 15 => b'N',
        _ => unreachable!(),
    }
}

let bases = buf.iter()
    .flat_map(|&b| [decode_base(b >> 4), decode_base(b)]);
```

**NEON Optimization Strategy**:
```rust
#[cfg(target_arch = "aarch64")]
unsafe fn decode_sequence_neon(packed: &[u8], output: &mut [u8]) {
    let lookup = vld1q_u8(LOOKUP_TABLE.as_ptr());

    for chunk in packed.chunks_exact(16) {
        let packed_bytes = vld1q_u8(chunk.as_ptr());

        // Process 16 packed bytes → 32 bases in parallel
        let high_nibbles = vshrq_n_u8(packed_bytes, 4);
        let low_nibbles = vandq_u8(packed_bytes, vdupq_n_u8(0x0F));

        let bases_high = vqtbl1q_u8(lookup, high_nibbles);  // 16 lookups at once
        let bases_low = vqtbl1q_u8(lookup, low_nibbles);    // 16 lookups at once

        // Interleave and store 32 decoded bases
    }
}
```

**Estimated speedup**:
- **Operation-level**: 3-6× (NEON parallel lookup vs sequential match)
- **Overall BAM parsing**: 2-3× (combined with parallel BGZF at 6.5×)

**C. Additional NEON Opportunities**:
- Parallel BGZF decompression: 6.5× (already have in biometal)
- Quality filtering: 16-25× (existing biometal ops, integration)
- Base counting: 16-25× (existing biometal ops, integration)
- MAPQ filtering: 16× (NEON parallel compare)

---

## Documents Created (Today)

### BAM Experiment (native-bam-implementation/)
1. **PROPOSAL.md** (14 KB) - Complete 9-12 week implementation plan
2. **PHASE_0_RESEARCH.md** (14 KB) - Format analysis, noodles study
3. **NEON_OPPORTUNITIES.md** (10 KB) - Detailed optimization analysis

### CAF Experiment (columnar-alignment-format/)
4. **PROPOSAL.md** (14 KB) - 6-week CAF implementation plan
5. **CAF_SPECIFICATION.md** (17 KB) - Complete format design (v1.0 draft)
6. **README.md** (7 KB) - Overview and use cases

### Paper Planning (alignment-formats-paper/)
7. **PAPER_OUTLINE.md** (15 KB) - Complete paper structure (4,000 words)
8. **BENCHMARKING_PLAN.md** (18 KB) - Comprehensive benchmark strategy
9. **README.md** (9 KB) - Paper overview and timeline

### Strategic Documents
10. **ALIGNMENT_FORMATS_ROADMAP.md** (15 KB) - 24-week roadmap

**Total**: ~120 KB of comprehensive documentation, 10 major documents

---

## Experiment Registry Updated

`.experiments.toml` now tracks:
1. **native-bam-implementation** (research phase, Week 1)
2. **columnar-alignment-format** (design phase, starts Week 13)
3. **alignment-formats-paper** (planning phase, 24-week timeline)

---

## Next Steps (Continuing Phase 0)

### Days 2-4: Specification Deep Dive
**Goal**: Document BAM format edge cases, complexity

**Tasks**:
- [ ] Read complete SAM/BAM v1.6 specification (50 pages)
- [ ] Document edge cases (unmapped reads, supplementary alignments, etc.)
- [ ] Map CIGAR complexity (which operations matter)
- [ ] Understand optional tags (which are common)

**Deliverable**: Updated PHASE_0_RESEARCH.md with edge cases section

### Days 5-6: ARM Profiling (CRITICAL)
**Goal**: Validate NEON opportunity hypothesis

**Approach**:
```rust
// Create profiling benchmark using noodles
// Measure % CPU time in decode_base / read_sequence

use noodles_bam as bam;
use std::time::Instant;

fn profile_parse(path: &str) {
    let mut reader = bam::io::Reader::new(File::open(path)?);
    let header = reader.read_header()?;

    let mut record = bam::Record::default();
    let mut count = 0;

    let start = Instant::now();
    while reader.read_record(&mut record)? != 0 {
        count += 1;
    }
    let elapsed = start.elapsed();

    println!("Parsed {} records in {:?}", count, elapsed);
    println!("Throughput: {} Mrec/s", count as f64 / elapsed.as_secs_f64() / 1e6);
}
```

**Profiling tools**:
- macOS Instruments (Time Profiler)
- cargo-flamegraph
- Manual timing measurements

**Critical validation**:
- Is `decode_base` ≥15% of CPU time? (GO if yes)
- Is `read_sequence` a bottleneck? (GO if yes)
- What other hotspots exist? (document)

**Test data** (to download/generate):
- Small: 10K Illumina records (~1 MB BAM)
- Medium: 100K Illumina records (~50 MB BAM)

**GO Criteria**:
- Sequence decoding ≥15% CPU time → High confidence 2× overall
- Sequence decoding ≥30% CPU time → Likely 3-4× overall

**NO-GO Risk**:
- Sequence decoding <10% CPU time → May not reach 2× threshold
- I/O dominates ≥80% → Limited optimization opportunity

### Day 7: GO/NO-GO Decision
**Goal**: Make evidence-based decision

**Decision Document** (`DECISION.md`):
```markdown
# Phase 0 GO/NO-GO Decision

**Date**: November 15, 2025
**Decision**: GO / NO-GO

## Profiling Results

- Sequence decoding: X% of CPU time
- BGZF decompression: Y% of CPU time
- Other: Z%

## NEON Opportunities Validated

1. Sequence decoding: X% → Estimated A× speedup
2. BGZF parallel: Y% → 6.5× speedup (proven)
3. Combined estimate: B× overall

## Decision Rationale

**GO if**: Overall estimate ≥2×, clear NEON path
**NO-GO if**: Estimate <2×, no clear optimization

## Outcome: [GO/NO-GO]

[Rationale based on evidence]
```

**If GO**: Proceed to Phase 1 (Minimal BAM Reader)
**If NO-GO**: Pivot to noodles integration (document learnings)

---

## Current Status

**Week 1, Day 1**: ✅ Complete
- Specification analysis ✅
- noodles source analysis ✅
- NEON opportunity identified ✅
- Paper planning complete ✅

**Week 1, Days 2-7**: Pending
- Specification deep dive (Days 2-4)
- ARM profiling validation (Days 5-6)
- GO/NO-GO decision (Day 7)

---

## Strategic Context

This BAM implementation is **Part 1 of dual-format paper**:

**Paper**: "Dual-Format Strategy for ARM-Native Bioinformatics: BAM and CAF"
**Components**:
1. ARM-native BAM (2-3× speedup, spec-compliant)
2. CAF (5-10× speedup, columnar format)
3. Lossless BAM ↔ CAF conversion
4. Decision framework (when to use each)
5. Comprehensive benchmarks (N=30)

**Publication target**: Bioinformatics (Oxford) or BMC Bioinformatics
**Expected submission**: Q2 2026 (~6 months)

---

## Key Insights (Day 1)

1. **NEON opportunity is real**: 4-bit sequence decoding is perfect NEON pattern (table lookup, data parallel)

2. **Profiling is critical**: Must validate actual CPU time % (Days 5-6)

3. **Dual-format strategy is strong**: BAM (compatibility) + CAF (performance) tells complete story

4. **Paper-driven development**: Having publication goal focuses implementation priorities

5. **Evidence-based approach**: Every decision based on measurements (profiling, benchmarks)

---

## Profiling Infrastructure (Day 1 Continuation)

**Created**: Profiling tool and test data infrastructure

### Profiling Tool ✅

**Location**: `profiling/`

```rust
// profiling/src/main.rs - Measures noodles performance
fn profile_parse(path: &str) -> io::Result<()> {
    let mut bam_reader = bam::io::Reader::new(reader);
    let mut record = bam::Record::default();

    // Parse all records, measure throughput
    while bam_reader.read_record(&mut record)? != 0 {
        // Access sequence and quality to ensure decoding
        let seq = record.sequence();
        let qual = record.quality_scores();
    }

    // Report: records/sec, Mbp/sec
}
```

**Files Created**:
- `profiling/Cargo.toml` - noodles-bam 0.68 dependencies
- `profiling/src/main.rs` - Profiling binary (compiles ✅)
- `profiling/README.md` - Complete usage guide
- `profiling/get_test_data.sh` - Download 1000 Genomes BAM
- `profiling/generate_test_bam.sh` - Generate synthetic BAM
- `profiling/PROFILING_STATUS.md` - Status and next steps

**Build**: ✅ Compiles successfully
```bash
cd profiling
cargo build --release
./target/release/bam-profiling <file.bam>
```

### Next Steps (Days 5-6)

**Step 1**: Install samtools (if generating synthetic data)
```bash
brew install samtools
```

**Step 2**: Get test data
```bash
./get_test_data.sh           # Download 1000 Genomes
# OR
./generate_test_bam.sh 100000  # Generate 100K records
```

**Step 3**: Run profiling
```bash
instruments -t "Time Profiler" ./target/release/bam-profiling ../test-data/test.bam
```

**Step 4**: Analyze results
- Look for CPU time in `decode_base` and `read_sequence`
- **GO Criteria**: ≥15% CPU time in sequence decoding
- **NO-GO Risk**: <10% CPU time

**Step 5**: Document findings in `PHASE_0_RESEARCH.md`

---

**Status**: Excellent progress! Day 1 complete, profiling infrastructure ready
**Next session**: Continue with specification deep dive (Days 2-4) and profiling (Days 5-6)
**Decision point**: November 15, 2025 (Day 7)
