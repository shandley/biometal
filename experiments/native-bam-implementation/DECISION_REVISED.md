# Phase 0 GO Decision (Revised)

**Date**: November 8, 2025
**Status**: **GO** (with strategic clarification)
**Decision Point**: End of Phase 0 (Day 7)

---

## Executive Summary

**Decision**: **GO** for native BAM implementation

**Primary Optimization**: Parallel BGZF decompression (captures 66-80% bottleneck)

**Strategic Rationale**: Full-stack ownership + ARM-native architecture aligns with biometal's mission. Speedup comes primarily from parallel BGZF (~6.5×), which we get regardless of parsing approach. Native parsing provides strategic value (no external dependencies, complete control, ARM-first design).

**Expected Speedup**: ~4-5× overall (primarily from parallel BGZF)

---

## Key Finding from Profiling

**CPU Time Distribution** (100K records, flamegraph):

| Category | CPU % | Optimization Strategy |
|----------|-------|----------------------|
| **BGZF decompression** | **66-80%** | ✅ Parallel BGZF (proven 6.5× speedup) |
| Record parsing | 6-13% | Native implementation (full control) |
| Sequence decoding | <6% | Defer NEON optimization (not bottleneck) |
| Memory operations | ~6% | Standard optimizations |

**Critical Insight**: BGZF decompression is the bottleneck, and this bottleneck exists **regardless** of whether we use native or noodles parsing.

---

## Decision Analysis

### Initial Assessment (Speed-focused)

**Original thinking**: NEON sequence decoding is insufficient (<6% CPU time), therefore NO-GO for native implementation.

**Conclusion**: Use noodles + parallel BGZF for faster time-to-market.

### Strategic Reassessment (Full-stack ownership)

**User insight**: If speedup is similar either way (both get ~4-5× from parallel BGZF), strategic considerations become decisive.

**Key question**: Do we want full-stack ownership or external dependency?

### Strategic Arguments FOR Native Implementation

**1. Full-Stack Ownership**
- biometal controls entire pipeline (parsing + decompression)
- No dependency on noodles maintenance/versioning
- Consistent architecture (all biometal, all ARM-native)

**2. ARM-First Design Philosophy**
- Even if parsing is only 6-13% CPU time, we can optimize it
- Aligns with biometal's mission (ARM-native bioinformatics)
- Future-proof for ARM-specific optimizations

**3. Complete Control**
- Understand every line of code
- Can add biometal-specific features (streaming, memory layout)
- Can optimize any part of the stack as needed

**4. No External Dependencies**
- Self-contained (no noodles versioning/breaking changes)
- Simpler deployment (one library, not two)
- Better for LMIC users (fewer dependencies)

**5. Foundation for CAF**
- Deep BAM format understanding informs CAF design
- Complete format expertise in-house
- Stronger paper narrative (native BAM → revolutionary CAF)

### Speedup Comparison

**Option A: noodles + parallel BGZF**
- Speedup: ~4-5× (from parallel BGZF)
- Effort: 500-1,000 LOC, 2-3 weeks
- Strategic: External dependency (noodles)

**Option B: native + parallel BGZF**
- Speedup: ~4-5× (from parallel BGZF) + marginal parsing gains
- Effort: 7,000-9,000 LOC, 9-12 weeks
- Strategic: Full-stack ownership (all biometal)

**Key realization**: Speedup is similar (~4-5× either way), so strategic value is the differentiator.

---

## Revised Decision

**GO** for native BAM implementation with these priorities:

### Phase 1-2: Minimal Native Parser + Parallel BGZF (Weeks 2-5)

**Focus**: Correctness + capturing the big bottleneck

```rust
// Native BAM implementation
pub struct BamReader {
    // Our parallel BGZF (6.5× proven speedup)
    bgzf: biometal::io::compression::ParallelBgzfReader,

    // Native parsing (full-stack ownership)
    parser: NativeBamRecordParser,

    // Streaming architecture (constant ~5 MB memory)
    buffer: RecordBuffer,
}
```

**Goals**:
- ✅ Native parsing (full biometal ownership)
- ✅ Parallel BGZF integration (captures 66-80% bottleneck)
- ✅ Spec-compliant (works with existing tools)
- ✅ Streaming architecture (constant memory)
- ⏸️ NEON optimizations (defer until profiling shows they matter)

### Phase 3: Profile Again (Week 6)

**Why re-profile**: With BGZF 6.5× faster, bottlenecks may shift

**Questions to answer**:
- Is sequence decoding now a bigger % of CPU time?
- Are there new bottlenecks revealed?
- Is NEON sequence optimization worth effort now?

**Expected**: Even with fast BGZF, sequence likely still <15% (not worth NEON effort)

### Phase 4-6: Production Polish (Weeks 7-12)

**Focus**: Completeness, edge cases, optimization (if needed)

- Complete BAM feature set
- Edge case handling (unmapped, supplementary, etc.)
- Comprehensive testing (differential vs noodles)
- NEON optimizations IF profiling shows they matter
- Documentation, examples

---

## Why This is Better Than Original NO-GO

### Original NO-GO Logic

**Focused on**: Optimization efficiency (1.05× for 7,000 LOC is poor ROI)

**Missed**: Strategic value (full-stack ownership, ARM-native philosophy)

### Revised GO Logic

**Recognizes**:
1. Speedup comes primarily from parallel BGZF (~6.5×) regardless of parsing approach
2. Native vs noodles parsing has similar overall speedup (~4-5× both)
3. Strategic considerations (full-stack, no dependencies) are decisive
4. biometal's mission is ARM-native bioinformatics (native aligns with mission)

**Result**: Native implementation provides strategic value for similar speedup effort

---

## BGZF: No Alternatives for BAM

**Question raised**: Are there alternatives to BGZF we should consider?

**Answer**: NO for BAM (BGZF is spec-required)

**Why BGZF is mandatory**:
1. **Spec compliance**: BAM specification requires BGZF compression
2. **Random access**: BAM index (BAI/CSI) assumes BGZF block structure
3. **Tooling compatibility**: All BAM tools expect BGZF format
4. **Block independence**: Each 64KB block decompresses independently

**Our optimization**: Parallelize BGZF decompression (proven 6.5× speedup)

**For CAF**: We CAN use modern compression (zstd, lz4) - this is a key differentiator!

---

## Updated Implementation Plan

### Phase 1: Minimal BAM Reader (Weeks 2-3)

**Goal**: Basic BAM parsing, spec-compliant

**Deliverables**:
- BAM header parsing (SAM header, reference sequences)
- Record parsing (positions, MAPQ, FLAGS, sequences, qualities, CIGAR)
- Streaming iterator interface
- Differential testing harness (vs noodles)

**Focus**: Correctness over speed

### Phase 2: Parallel BGZF Integration (Weeks 3-4)

**Goal**: Integrate biometal's parallel BGZF (THE BIG WIN)

**Deliverables**:
- Replace sequential BGZF with parallel version
- Benchmark: baseline vs parallel
- Validate ~6.5× BGZF speedup translates to ~4-5× overall
- Profile to confirm BGZF bottleneck is resolved

**Expected**: ~4-5× overall speedup

### Phase 3: Streaming Architecture (Week 4)

**Goal**: Constant memory (~5 MB) regardless of file size

**Deliverables**:
- Block-based processing (Rule 2: 10,000 records)
- Record reuse pattern (zero-copy where possible)
- Memory profiling (validate constant footprint)

**Target**: ~5 MB memory for any BAM file size

### Phase 4: ARM NEON Optimization (IF NEEDED) (Weeks 5-6)

**Goal**: Profile-driven optimization

**Approach**:
1. Profile native + parallel BGZF implementation
2. Identify remaining bottlenecks (if any)
3. Apply NEON optimizations ONLY if >15% CPU time
4. Re-measure, validate improvement

**Expected**: Probably no NEON needed (BGZF was the bottleneck)

### Phase 5: SAM Support (Optional) (Week 7)

**Goal**: Support uncompressed SAM files

**Deliverables**:
- SAM text parsing
- Unified interface (same API for BAM/SAM)

**Benefit**: Complete format coverage

### Phase 6: Production Polish (Weeks 8-12)

**Goal**: Production-ready quality

**Deliverables**:
- Edge case handling (200+ test cases)
- Complete documentation
- Examples, benchmarks
- Error handling, validation
- CLI tools (if needed)

**Quality**: Ready for biometal 2.0 release

---

## Expected Outcomes

### Performance

**Baseline** (noodles):
- Throughput: ~8.7 Mrec/sec (measured)
- Bandwidth: ~873 Mbp/sec

**Native + Parallel BGZF** (estimated):
- Throughput: ~35-40 Mrec/sec (4-5× improvement)
- Bandwidth: ~3.5-4.0 Gbp/sec
- Memory: ~5 MB constant (Rule 5)

**Source of speedup**:
- Parallel BGZF: 6.5× on 70% CPU time = ~4.5× overall
- Native parsing optimizations: marginal (5-10%)
- **Combined**: ~4-5× total

### Strategic Value

**Full-Stack Ownership**:
- ✅ No external dependencies (noodles)
- ✅ Complete ARM-native stack
- ✅ biometal controls entire pipeline
- ✅ Foundation for CAF design

**ARM Optimization**:
- ✅ Parallel BGZF (ARM-optimized threading)
- ✅ Native parsing (can add ARM-specific optimizations later)
- ✅ Streaming architecture (memory-efficient for all ARM devices)

**Mission Alignment**:
- ✅ Aligns with biometal's ARM-native philosophy
- ✅ Demonstrates full-stack capability
- ✅ Positions biometal as complete solution

---

## Paper Impact

### Revised Paper Structure

**Title**: "Dual-Format Strategy for ARM-Native Bioinformatics: Native BAM Optimization and CAF"

**BAM Section**:
1. **Profiling Analysis** (Evidence-based approach)
   - Identified BGZF as bottleneck (66-80% CPU time)
   - Measured sequence decoding (<6%, not primary target)

2. **Native Implementation** (Full-stack ownership)
   - Built ARM-native BAM parser
   - Integrated parallel BGZF (6.5× proven speedup)
   - Achieved ~4-5× overall speedup

3. **Strategic Decision** (Full-stack vs integration)
   - Chose native for complete control
   - Achieved similar speedup (~4-5×) either way
   - Demonstrates ARM-first design philosophy

**CAF Section** (unchanged):
- Columnar format, modern compression
- 5-10× speedup for analytical operations
- Freedom from legacy format constraints

**Paper Strength**: Shows both pragmatic optimization (target real bottlenecks) AND strategic design (full-stack ownership)

---

## Timeline

**Total**: 9-12 weeks (as originally planned)

| Phase | Weeks | Focus | Deliverable |
|-------|-------|-------|-------------|
| Phase 1 | 2-3 | Minimal parser | Correct BAM parsing |
| Phase 2 | 3-4 | Parallel BGZF | ~4-5× speedup |
| Phase 3 | 4 | Streaming | ~5 MB constant memory |
| Phase 4 | 5-6 | Profile/optimize | Evidence-based NEON (if needed) |
| Phase 5 | 7 | SAM support | Complete format coverage |
| Phase 6 | 8-12 | Production | 200+ tests, docs |

**Followed by**:
- CAF implementation (Weeks 13-18)
- Benchmarking (Weeks 19-20)
- Paper writing (Weeks 21-24)

---

## Success Criteria

### Phase 0 (Complete ✅)

- ✅ Profiling data collected (100K records)
- ✅ Bottleneck identified (BGZF 66-80%)
- ✅ GO/NO-GO decision made (GO with strategic rationale)

### Phase 1-2 (Weeks 2-5)

- ✅ Native BAM parsing working
- ✅ Parallel BGZF integrated
- ✅ ~4-5× speedup achieved
- ✅ Differential testing passes

### Phase 3 (Week 6)

- ✅ Constant memory (~5 MB)
- ✅ Profile confirms no new bottlenecks
- ✅ Streaming architecture validated

### Phase 4-6 (Weeks 7-12)

- ✅ 200+ tests passing
- ✅ Complete documentation
- ✅ Production-ready quality
- ✅ Ready for biometal 2.0

---

## Conclusion

**Phase 0 SUCCESS**: Profiling revealed actual bottleneck (BGZF, not sequence decoding)

**Decision**: **GO** for native BAM implementation

**Rationale**:
1. ✅ Full-stack ownership (strategic value)
2. ✅ Parallel BGZF captures bottleneck (~4-5× speedup)
3. ✅ Aligns with biometal's ARM-native mission
4. ✅ Similar speedup to alternative (~4-5× either way)
5. ✅ Foundation for CAF design

**Timeline**: 9-12 weeks (as planned)

**Next**: Begin Phase 1 (Minimal BAM Reader, Week 2)

---

**Signed**: biometal development team
**Date**: November 8, 2025
**Phase 0 Status**: Complete ✅
**Decision**: GO ✅
**Primary Optimization**: Parallel BGZF (66-80% bottleneck) ✅
**Strategic Value**: Full-stack ARM-native ownership ✅
