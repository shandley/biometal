# ARM NEON Optimization Opportunities

**Date**: November 8, 2025 (Day 1 Analysis)
**Status**: Preliminary (profiling validation pending)
**Based on**: noodles source code analysis

---

## Executive Summary

**Hypothesis**: Native BAM implementation can achieve ≥2× ARM speedup via NEON

**Initial Findings**:
- ✅ **1 HIGH-CONFIDENCE NEON opportunity identified** (sequence decoding)
- ✅ **2 INTEGRATION opportunities** (leverage existing biometal ops)
- ⚠️ **Profiling needed** to validate actual bottlenecks

---

## Opportunity #1: 4-bit Sequence Decoding (HIGH CONFIDENCE)

### Current Implementation (noodles)

**File**: `noodles-bam/src/record/codec/decoder/sequence.rs:36-81`

```rust
pub(super) fn read_sequence(
    src: &mut &[u8],
    sequence: &mut Sequence,
    base_count: usize,
) -> Result<(), DecodeError> {
    let len = base_count.div_ceil(2);
    let (buf, rest) = src.split_at_checked(len)?;
    *src = rest;

    // Sequential processing: 1 nibble per iteration
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

### NEON Optimization Strategy

**Pattern**: Table lookup with SIMD

```rust
#[cfg(target_arch = "aarch64")]
unsafe fn decode_sequence_neon(packed: &[u8], output: &mut Vec<u8>) {
    use std::arch::aarch64::*;

    // Lookup table for 4-bit -> ASCII base
    const LOOKUP: [u8; 16] = [
        b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V',
        b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
    ];

    let lookup_table = vld1q_u8(LOOKUP.as_ptr());

    // Process 16 packed bytes (32 bases) at a time
    for chunk in packed.chunks_exact(16) {
        let packed_bytes = vld1q_u8(chunk.as_ptr());

        // Extract high nibbles (first 16 bases)
        let high_nibbles = vshrq_n_u8(packed_bytes, 4);
        let bases_high = vqtbl1q_u8(lookup_table, high_nibbles);

        // Extract low nibbles (next 16 bases)
        let low_nibbles = vandq_u8(packed_bytes, vdupq_n_u8(0x0F));
        let bases_low = vqtbl1q_u8(lookup_table, low_nibbles);

        // Interleave and store
        let interleaved_low = vzip1q_u8(bases_high, bases_low);
        let interleaved_high = vzip2q_u8(bases_high, bases_low);

        vst1q_u8(output.as_mut_ptr().add(output.len()), interleaved_low);
        vst1q_u8(output.as_mut_ptr().add(output.len() + 16), interleaved_high);

        output.set_len(output.len() + 32);
    }

    // Handle remainder (< 16 bytes) with scalar fallback
}
```

### Performance Estimate

**Current (noodles scalar)**:
- Processing: 1 nibble/iteration (sequential match)
- Estimated: ~2-4 cycles/base (branch misprediction + match overhead)
- For 100 bp read: ~200-400 cycles

**NEON optimized**:
- Processing: 32 bases in parallel (16-byte chunk)
- Operations: 1 load + 2 shifts + 2 table lookups + 2 interleaves + 2 stores
- Estimated: ~20 cycles per 32 bases = **0.625 cycles/base**
- For 100 bp read: ~63 cycles

**Speedup**: **3-6× for sequence decoding operation**

**Caveat**: This is per-operation speedup. Overall BAM parsing speedup depends on what % of time is spent in sequence decoding (profiling needed).

### Implementation Complexity

**Effort**: LOW-MEDIUM
- Pattern already proven in biometal (base counting uses similar NEON intrinsics)
- Reference implementation available (noodles source)
- Clear specification (SAM/BAM v1.6)

**Risks**:
- Edge case: Odd-length sequences (last byte has unused nibble)
- Correctness: Must match noodles exactly (differential testing)

---

## Opportunity #2: Quality Score Operations (INTEGRATION)

### Pattern

BAM records contain Phred+33 quality scores (same as FASTQ)

### Optimization

**Leverage existing biometal operations**:
- `quality_filter_neon`: 16-25× speedup (proven)
- `calculate_mean_quality_neon`: Similar speedup expected

### Implementation

```rust
pub struct BamRecord {
    // ... other fields
    quality_scores: Vec<u8>,
}

impl BamRecord {
    pub fn mean_quality(&self) -> f64 {
        // Use biometal's ARM-optimized operation
        biometal::operations::calculate_mean_quality(&self.quality_scores)
    }

    pub fn passes_quality_filter(&self, min_quality: u8) -> bool {
        // Use biometal's NEON filter
        biometal::operations::quality_filter(&self.quality_scores, min_quality)
    }
}
```

### Performance Estimate

**Integration-only** (not new NEON work):
- Speedup: 16-25× for quality operations
- Overall impact: Depends on user workflow (filtering-heavy workloads benefit most)

---

## Opportunity #3: Base Counting (INTEGRATION)

### Pattern

Common analysis: Count A/C/G/T/N distribution across alignments

### Optimization

**Leverage biometal's count_bases_neon**:
- Proven 16-25× speedup
- Works on decoded sequences

### Implementation

```rust
impl BamRecord {
    pub fn base_counts(&self) -> [u32; 4] {
        // After decoding sequence to ASCII
        biometal::operations::count_bases(&self.sequence)
    }
}
```

### Performance Estimate

**Integration-only**:
- Speedup: 16-25× for base counting
- Overall impact: Depends on analysis type

---

## Opportunity #4: BGZF Decompression (ALREADY HAVE)

### Pattern

BAM files use BGZF (Blocked GZip Format) compression

### Optimization

**Use biometal's parallel BGZF**:
- Already implemented: `decompress_bgzip_parallel`
- Proven 6.5× speedup

### Implementation

```rust
pub struct BamReader<R: Read> {
    bgzf_reader: biometal::io::CompressedReader<R>,
    // ... rest
}
```

### Performance Estimate

**Already implemented**:
- Speedup: 6.5× for decompression
- Impact: Significant for large files (I/O-bound scenarios)

---

## Opportunities REJECTED

### CIGAR String Parsing

**Analysis**: Already optimal
- Simple bitwise operations (1-2 cycles)
- No data parallelism opportunity
- NEON would add overhead

**Decision**: No NEON optimization

### FLAG Bit Checking

**Analysis**: Trivial cost
- Single uint16 bitwise AND (1 cycle)
- No benefit from SIMD

**Decision**: No optimization needed

---

## Overall Speedup Estimate

### Pessimistic Scenario

**Assumptions**:
- Sequence decoding: 20% of parse time → 3× faster → **1.4× overall**
- BGZF decompression: 30% of time → 6.5× faster → **1.6× overall**
- Combined: **~2.2× total speedup**

### Optimistic Scenario

**Assumptions**:
- Sequence decoding: 40% of parse time → 5× faster → **2.5× overall**
- BGZF decompression: 30% of time → 6.5× faster → **1.6× overall**
- Quality operations (if used): 10% of time → 20× faster → **1.15× overall**
- Combined: **~4-5× total speedup**

### Realistic Estimate (Pending Profiling)

**Expected**: **2-3× overall ARM speedup**
- Primary: Sequence decoding + BGZF (proven patterns)
- Secondary: Quality operations (user workflow-dependent)

**Threshold**: ≥2× speedup required (PROPOSAL.md:362)
**Status**: ✅ **LIKELY ACHIEVABLE** (pending validation)

---

## Validation Plan

### Phase 0 Profiling (This Week)

**Objective**: Validate sequence decoding is actual bottleneck

**Method**:
1. Profile noodles with real BAM files (Instruments Time Profiler)
2. Measure % CPU time in `decode_base` / `read_sequence`
3. Validate assumptions about hotspots

**GO Criteria**:
- Sequence decoding ≥15% of parse time → High-confidence 2× overall
- Sequence decoding ≥30% of parse time → Likely 3-4× overall

**NO-GO Risk**:
- Sequence decoding <10% of parse time → May not reach 2× threshold
- I/O dominates ≥80% → Limited optimization opportunity

### Phase 1 Prototype (Week 2-3)

**Objective**: Implement NEON sequence decoder

**Validation**:
- Correctness: Differential testing vs noodles (100% match required)
- Performance: Isolated benchmark (expect 3-6× for operation)

### Phase 3 Integration (Week 5-6)

**Objective**: Measure end-to-end speedup

**Validation**:
- Full BAM parsing benchmark
- Target: ≥2× overall ARM vs noodles scalar
- Stretch: ≥3× overall

---

## Risk Assessment

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Sequence decoding not hotspot | Medium | High | **Profile first** (this week) |
| NEON speedup <3× per operation | Low | Medium | Pattern proven in biometal |
| Overall speedup <2× | Low-Medium | Medium | BGZF gives 1.6× floor already |
| Edge cases break NEON | Low | Low | Comprehensive differential testing |

### Strategic Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Profiling reveals no clear target | Low | High | Pivot to noodles integration |
| Implementation complexity exceeds estimate | Medium | Low | Phased approach with decision points |

---

## Phase 0 Decision Criteria

**GO to Phase 1 if**:
- ✅ Sequence decoding ≥15% of CPU time (profiling)
- ✅ NEON prototype shows ≥3× operation speedup
- ✅ No architectural blockers

**NO-GO (pivot to noodles) if**:
- ❌ Sequence decoding <10% of CPU time
- ❌ Fundamental NEON implementation issues

**Current Status**: Day 1
- ✅ Source analysis complete
- ✅ NEON opportunity identified (high confidence)
- ⏳ Profiling pending (Day 5-6)
- ⏳ Decision point: Day 7

---

## Next Steps

### This Week (Phase 0)

**Day 2-4**: Complete specification study
**Day 5-6**: Profile noodles on ARM
- Download diverse BAM files (Illumina, PacBio, ONT)
- Instruments Time Profiler
- Identify actual hotspots
- Measure sequence decoding % of time

**Day 7**: Decision
- If profiling validates: **GO to Phase 1**
- If profiling invalidates: **Pivot to noodles integration**

---

**Status**: Day 1 - Source analysis complete, high-confidence NEON target identified
**Next**: Profiling validation (Day 5-6)
**Decision Point**: End of Week 1
