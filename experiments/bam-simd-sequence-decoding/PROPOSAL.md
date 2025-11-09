# Experiment: BAM SIMD Sequence Decoding

**Status**: 🔬 Proposed
**Date**: November 9, 2025
**Researcher**: Claude (AI assistant)
**Duration**: 2-3 weeks (estimated)

---

## Executive Summary

Implement ARM NEON SIMD acceleration for BAM 4-bit to ASCII sequence decoding. Expected 16-25× speedup could translate to 10-15% overall BAM parsing improvement now that BGZF bottleneck has been optimized.

---

## Hypothesis

ARM NEON-accelerated 4-bit to ASCII sequence decoding can achieve **16-25× speedup** vs scalar implementation, translating to **10-15% overall BAM parsing speedup** (now that BGZF has been optimized from 70% → ~30% of CPU time).

---

## Motivation

### Current State (Post-BGZF Optimization)

After v1.4.0's parallel BGZF optimization:
- **BGZF decompression**: 70% → ~30% (4× faster, but still 30% of time)
- **Sequence decoding**: <6% → **likely ~10-15%** (proportionally larger)
- **Record parsing**: 6-13% → **likely ~15-20%**
- **Tag parsing**: NEW in v1.4.0, adds ~2-3%

**Conclusion**: Sequence decoding is now likely ≥15% of CPU time, meeting Rule 1 threshold for SIMD optimization.

### Evidence Base

**OPTIMIZATION_RULES.md Rule 1**:
- **When**: Element-wise operations (complexity 0.30-0.40), ≥15% CPU time
- **Expected**: 16-25× speedup using ARM NEON SIMD
- **Evidence**: Entry 020-025 (307 experiments, 9,210 measurements)
- **Platform**: Mac ARM (M1/M2/M3/M4) optimized

**Current Implementation** (scalar):
- Simple lookup table per nibble
- Loop through sequence byte-by-byte
- ~6% CPU time in Phase 0 profiling (pre-BGZF optimization)

**NEON Potential**:
- Process 16 bytes at once (32 bases)
- Use VTBL (vector table lookup) for 4-bit → ASCII
- Cache-friendly, memory-bound operation

---

## Research Questions

### Primary Questions

1. **Q1**: What is sequence decoding's CPU time percentage after BGZF optimization?
   - **Method**: Profile 100K record parsing with Instruments (macOS)
   - **Target**: Validate ≥15% threshold for Rule 1

2. **Q2**: What speedup does NEON achieve for isolated sequence decoding?
   - **Method**: Microbenchmark NEON vs scalar (N=30, criterion)
   - **Target**: 16-25× speedup (Rule 1 expectation)

3. **Q3**: What overall BAM parsing speedup does NEON sequence decoding provide?
   - **Method**: Full BAM parsing benchmark (100K records, N=30)
   - **Target**: 10-15% overall improvement

### Secondary Questions

4. **Q4**: Does NEON maintain correctness for all IUPAC ambiguity codes?
   - **Method**: Property-based testing (proptest, 10K cases)
   - **Target**: 100% correctness

5. **Q5**: What is memory bandwidth utilization for sequence decoding?
   - **Method**: Performance counters (Apple AMX monitoring)
   - **Target**: Understand bottleneck (compute vs memory)

---

## Methodology

### Phase 1: Profiling and Validation (Days 1-3)

**Goal**: Validate that sequence decoding is now ≥15% of CPU time.

**Tasks**:
1. Profile current BAM parser with Instruments
2. Measure CPU time breakdown:
   - BGZF decompression (expected ~30%)
   - Sequence decoding (expected ~10-15%)
   - Record parsing (expected ~15-20%)
   - Tag parsing (expected ~2-3%)
   - Other overhead (expected ~30-35%)
3. Document findings in `RESEARCH_LOG.md`

**Deliverable**: Profiling report with CPU time breakdown

### Phase 2: NEON Prototype Implementation (Days 4-10)

**Goal**: Implement ARM NEON 4-bit to ASCII decoder.

**Approach**: VTBL (Vector Table Lookup) Strategy

```rust
#[cfg(target_arch = "aarch64")]
unsafe fn decode_sequence_neon(data: &[u8], length: usize) -> Vec<u8> {
    use std::arch::aarch64::*;

    // Create lookup table in NEON register
    // BAM encoding: =ACMGRSVTWYHKDBN (indices 0-15)
    let lookup_table = vld1q_u8([
        b'=', b'A', b'C', b'M',
        b'G', b'R', b'S', b'V',
        b'T', b'W', b'Y', b'H',
        b'K', b'D', b'B', b'N',
    ].as_ptr());

    // Process 16 bytes (32 bases) at a time
    for chunk in data.chunks(16) {
        // Load 16 bytes
        let packed = vld1q_u8(chunk.as_ptr());

        // Extract high nibbles (first 16 bases)
        let high_nibbles = vshrq_n_u8(packed, 4);
        let high_bases = vqtbl1q_u8(lookup_table, high_nibbles);

        // Extract low nibbles (second 16 bases)
        let low_nibbles = vandq_u8(packed, vdupq_n_u8(0x0F));
        let low_bases = vqtbl1q_u8(lookup_table, low_nibbles);

        // Interleave high and low bases
        // ... (implementation detail)
    }
}
```

**Fallback**: Maintain scalar implementation for x86_64

**Tasks**:
1. Implement NEON decoder in `src/io/bam/sequence_neon.rs`
2. Add platform-specific dispatch:
   ```rust
   pub fn decode_sequence(data: &[u8], length: usize) -> io::Result<Vec<u8>> {
       #[cfg(target_arch = "aarch64")]
       { unsafe { decode_sequence_neon(data, length) } }

       #[cfg(not(target_arch = "aarch64"))]
       { decode_sequence_scalar(data, length) }
   }
   ```
3. Implement scalar fallback (existing implementation)
4. Add microbenchmark in `benches/sequence_decode.rs`

**Deliverable**: NEON implementation with benchmarks

### Phase 3: Validation and Benchmarking (Days 11-14)

**Goal**: Validate correctness and measure performance.

**Correctness Testing**:
1. Property-based tests (10K random sequences)
2. Edge cases: empty, odd length, all ambiguity codes
3. Comparison: NEON vs scalar output must match exactly

**Performance Benchmarking**:
1. **Microbenchmark** (isolated sequence decoding):
   - Input: 100K bases (50KB packed data)
   - N=30 samples, criterion framework
   - Compare: NEON vs scalar
   - Expected: 16-25× speedup

2. **Integration Benchmark** (full BAM parsing):
   - Input: 100K record BAM file (synthetic_100000.bam)
   - N=30 samples, criterion framework
   - Compare: With NEON vs without NEON
   - Expected: 10-15% overall speedup

**Deliverable**: Benchmark results and validation report

### Phase 4: Documentation (Days 15-17)

**Goal**: Document findings and integrate into codebase.

**Tasks**:
1. Update `FINDINGS.md` with results
2. Update `OPTIMIZATION_RULES.md` if new insights
3. Update BAM performance documentation (`docs/BAM_PERFORMANCE.md`)
4. Add to CHANGELOG.md as experimental feature
5. Update `.experiments.toml` registry

**Deliverable**: Complete documentation

---

## Go/No-Go Criteria

### Go Criteria (Proceed to Integration)

**Must Have** (all required):
1. ✅ **Q1 Validated**: Sequence decoding is ≥15% of CPU time
2. ✅ **Q2 Achieved**: NEON achieves ≥10× speedup vs scalar (microbenchmark)
3. ✅ **Q3 Achieved**: Overall BAM parsing improves by ≥5%
4. ✅ **Q4 Validated**: 100% correctness (all property tests pass)

**Nice to Have** (optional):
- Q5: Memory bandwidth analysis completed
- NEON speedup ≥16× (Rule 1 expected range)
- Overall speedup ≥10%

### No-Go Criteria (Abandon or Defer)

**Abandon if**:
1. ❌ Sequence decoding is <10% CPU time (too small to optimize)
2. ❌ NEON speedup is <5× (not worth complexity)
3. ❌ Overall BAM parsing improvement is <2% (negligible)
4. ❌ Correctness issues cannot be resolved

**Defer if**:
- Other higher-priority bottlenecks identified (e.g., tag parsing >20%)
- Platform compatibility issues arise
- Time-box exceeded (3 weeks)

---

## Success Metrics

### Primary Metrics

| Metric | Baseline | Target | Stretch Goal |
|--------|----------|--------|--------------|
| Sequence decode speedup | 1× (scalar) | 10-16× | 16-25× |
| BAM parsing speedup | 43.0 MiB/s | 47-49 MiB/s (+10%) | 49-50 MiB/s (+15%) |
| Record throughput | 4.4M rec/s | 4.8M rec/s (+10%) | 5.0M rec/s (+15%) |
| Correctness | 100% | 100% | 100% |

### Secondary Metrics

| Metric | Target |
|--------|--------|
| Memory usage | No regression (constant ~5 MB) |
| Code complexity | <200 LOC for NEON implementation |
| Test coverage | ≥95% for new code |
| Platform support | ARM + x86_64 fallback |

---

## Risks and Mitigation

### Risk 1: Sequence Decoding Still <15% CPU Time

**Probability**: Medium (30%)
**Impact**: High (blocks SIMD optimization)

**Mitigation**:
- Profile first (Phase 1) before implementing NEON
- If <15%, document findings and defer optimization
- Consider alternative targets (e.g., tag parsing)

### Risk 2: NEON Speedup Lower Than Expected

**Probability**: Low (15%)
**Impact**: Medium (reduces overall benefit)

**Mitigation**:
- Memory-bound operations may not achieve 16-25× (compute-bound expectation)
- Even 8-10× speedup would provide 5-8% overall improvement
- Document actual speedup and adjust expectations

### Risk 3: Platform Compatibility Issues

**Probability**: Low (10%)
**Impact**: Low (fallback to scalar)

**Mitigation**:
- Always maintain scalar fallback
- Test on multiple ARM platforms (M1, M2, M3, M4)
- Document platform-specific behavior

### Risk 4: Correctness Issues with IUPAC Codes

**Probability**: Low (10%)
**Impact**: High (blocks integration)

**Mitigation**:
- Extensive property-based testing
- Manual validation of lookup table
- Comparison testing (NEON vs scalar must match)

---

## Timeline

| Phase | Duration | Tasks |
|-------|----------|-------|
| Phase 1: Profiling | 3 days | Profile, validate ≥15% threshold |
| Phase 2: Implementation | 7 days | NEON decoder, microbenchmarks |
| Phase 3: Validation | 4 days | Property tests, benchmarks |
| Phase 4: Documentation | 3 days | Findings, integration docs |
| **Total** | **17 days** | **(~2.5 weeks)** |

**Buffer**: +3 days for unexpected issues (total: 3 weeks)

---

## Expected Outcomes

### Best Case (Go + Stretch Goals)

- Sequence decoding: 20-25× NEON speedup
- Overall BAM parsing: 12-15% faster (49-50 MiB/s, 5.0M rec/s)
- Rule 1 validated for memory-bound operations
- New entry added to OPTIMIZATION_RULES.md
- Production-ready feature in v1.5.0

### Likely Case (Go)

- Sequence decoding: 12-16× NEON speedup
- Overall BAM parsing: 8-10% faster (46-47 MiB/s, 4.8M rec/s)
- Validated approach for SIMD in BAM operations
- Experimental feature, documented for future work

### Worst Case (No-Go)

- Sequence decoding still <10% CPU time (deferred)
- Document findings: BGZF still dominates, focus elsewhere
- Recommend: BAM indexing (region queries) or VCF parser
- No code integrated, clean no-go decision

---

## Resources Required

### Hardware

- M3 MacBook Pro (ARM64) - Primary development
- Access to M1/M2/M4 for cross-validation (optional)

### Software

- Rust 1.82+ (NEON intrinsics)
- Instruments (macOS profiling)
- Criterion (benchmarking, N=30)
- Proptest (property-based testing)

### Data

- `synthetic_100000.bam` (100K records, ~946 KB)
- Located at: `experiments/native-bam-implementation/test-data/`

### Time

- 2-3 weeks (17-20 days)
- ~3-4 hours/day sustained focus

---

## Related Work

### biometal

- **OPTIMIZATION_RULES.md Rule 1**: ARM NEON SIMD (16-25× speedup)
- **Entry 020-025**: 307 experiments, 9,210 measurements (base counting, GC content, quality filtering)
- **Phase 0 Profiling**: BGZF 66-80%, sequence <6%, record 6-13%
- **v1.4.0**: Parallel BGZF (4× speedup), shifts bottleneck proportions

### External

- **noodles**: Pure Rust genomics library (reference for format correctness)
- **samtools**: C implementation (reference for performance comparison)
- **htslib**: C library with SIMD in some paths (not documented)

---

## Next Steps

1. **Approve Proposal**: Review and approve this experiment
2. **Phase 1**: Profile current BAM parser (validate ≥15% threshold)
3. **Go/No-Go Decision**: Based on Phase 1 profiling results
4. **Phase 2-4**: Implement, validate, document (if go)
5. **Integration**: Merge into main codebase (if successful)

---

**Proposal Status**: ✅ Ready for Review
**Recommended Action**: Approve and proceed to Phase 1 (Profiling)
