# simd-minimizers-analysis: Status & Next Steps

**Date**: November 7, 2025 (FINAL UPDATE)
**Status**: ✅ **COMPLETE** - SIMD Integration Validated and Documented

---

## ✅ FINAL RESULTS (Nov 7, 2025)

### Actual Path Taken: Library Integration (Not Porting)

**Decision (Nov 7)**: Instead of porting ~500 LOC from simd-minimizers, we integrated
the library as a dependency. This approach was faster and leverages ongoing maintenance
from the original authors.

**Implementation**:
- Added simd-minimizers v2.2 as optional dependency
- Created `extract_minimizers_simd()` wrapper function
- Recomputed ntHash values for API compatibility
- Opt-in via `features = ["simd"]` in Cargo.toml

### Performance Results (Entry 036-E)

**Validated benchmarks** (Mac M1, RUSTFLAGS="-C target-cpu=native", N=20):

| Sequence Size | Fast (MiB/s) | SIMD (MiB/s) | Speedup | Status |
|---------------|--------------|--------------|---------|--------|
| 100bp         | 225          | 92           | 0.41×   | Regression |
| 1Kbp          | 222          | 183          | 0.82×   | Regression |
| 10Kbp         | 225          | 195          | 0.87×   | Regression |
| **100Kbp**    | **122**      | **184**      | **1.51×** | **Crossover** |
| **1Mbp**      | **99.4**     | **193**      | **1.94×** | **Production** ✓ |
| **10Mbp**     | **98.2**     | **147**      | **1.50×** | **Production** ✓ |

### Gap from Original Expectations

**Expected** (from original analysis):
- 100-200× speedup based on Entry 036 baseline (3.7 Mbp/s)
- Direct port of ntHash + two stacks + SIMD

**Achieved** (actual implementation):
- **1.5-1.9× speedup** for production workloads (≥100Kbp)
- **0.41-0.87× regression** for small sequences (<100Kbp)
- Library integration (not port)

**Root cause of gap**:
1. **Double hashing overhead**: SIMD positions + ntHash recomputation (~2× cost)
2. **Chunking overhead**: Splits sequence into 8 chunks (dominates for small inputs)
3. **Canonical overhead**: Strand-aware selection adds complexity

### Decision: Accept and Document

✅ **GO** - 1.5-1.9× is valuable for production genomics (≥100Kbp sequences)

**Rationale**:
- Meaningful improvement for real-world workloads (bacterial genomes, contigs)
- Clear crossover point (~100Kbp) guides users effectively
- Opt-in feature flag minimizes risk
- Evidence-based documentation prevents misuse
- Gap from literature is explained and acceptable

### Deliverables

- ✅ `extract_minimizers_simd()` in `src/operations/kmer.rs:530`
- ✅ Benchmark: `benches/minimizer_simd.rs`
- ✅ Test: `test_extract_minimizers_simd_correctness`
- ✅ Documentation: Entry 036-E results and usage guidelines
- ✅ Entry 036-E summary: ENTRY_036_E_SIMD_VALIDATION.md

### Lessons Learned

1. **Library integration > porting**: Faster implementation, ongoing maintenance
2. **Double hashing matters**: API compatibility has measurable performance cost
3. **Scaling is non-linear**: SIMD overhead dominates for small inputs
4. **Evidence-based shipping**: 1.5-1.9× is valuable, don't wait for perfect

---

## Historical Record: Original Plan (Nov 6, 2025)

Below is the original analysis and plan from November 6. We ultimately took a
different approach (library integration instead of porting), but the research
and baseline validation remain valuable.

### Original Executive Summary

✅ **GO DECISION VALIDATED** - Entry 036 baseline reveals 221× speedup opportunity (not 4-8×!)

**Original baseline**: 3.7 Mbp/s (k=21, w=11, 100K sequences)
**Original target (realistic)**: 370 Mbp/s (100× speedup)
**Original target (exceptional)**: 555 Mbp/s (150× speedup)

**Note**: These targets were based on the assumption of directly porting the SIMD
implementation. The library integration approach achieved 1.5-1.9× instead due to
double hashing overhead

---

## Progress Review: What We've Accomplished

### ✅ Day 1: Source Code Analysis (Nov 4)
**Deliverables**:
- PROPOSAL.md: Hypothesis, go/no-go criteria, 1-week timeline
- Cloned rust-seq/simd-minimizers repository
- Analyzed algorithm structure (ntHash + two stacks + SIMD)
- Identified NEON compatibility (portable SIMD via packed-seq + wide)

**Key Finding**: Algorithm is clear, well-documented, MIT licensed (low integration risk)

### ✅ Day 2: Benchmarking (Nov 5)
**Deliverables**:
- Built simd-minimizers on Mac M-series (ARM NEON)
- Created custom benchmark (stable Rust compatible)
- Measured performance: 820.62 Mbp/s forward, 640.53 Mbp/s canonical

**Key Finding**: 8-16× faster than Entry 034 estimates (validates GO decision)

### ✅ Day 3: Algorithm Deep Dive & Streaming Adaptation (Nov 6)
**Deliverables**:
- Analyzed ntHash mechanism (vectorizable: table lookup + rotate + XOR)
- Analyzed two stacks algorithm (O(1) amortized sliding minimum)
- Created streaming_adaptation_prototype.rs (block-based, 10K blocks)
- Assessed memory trade-offs (97-99.99% reduction)

**Key Finding**: Block-based streaming preserves SIMD speedup while maintaining O(1) memory

### ✅ GO DECISION (Nov 6)
**Deliverables**:
- GO_DECISION.md: All 5 criteria met with high confidence
- Integration roadmap: 3 phases, 3 weeks to v1.3.0

**Key Finding**: Trade-off is acceptable (25% slower for 99.99% memory reduction)

### ✅ Entry 036 Baseline (Nov 6) - CRITICAL VALIDATION
**Deliverables**:
- Rigorous baseline measurement (N=100, 1,600 measurements)
- Statistical validation (95% CI, CV < 2%)
- ASBB Entry 036 complete and committed
- ENTRY_035_BASELINE_RESULTS.md: Comprehensive analysis

**Key Finding**: Baseline is 221× slower than SimdMinimizers (not 8-16×!)
- Entry 034 overestimated by 10-20×
- Revised speedup potential: 100-200× (not 4-8×!)
- Opportunity is 12-25× larger than originally estimated

### ✅ Applicability Analysis (Nov 6)
**Deliverables**:
- APPLICABILITY_ANALYSIS.md: Evidence-based evaluation
- Confirmed ntHash + two stacks applies ONLY to minimizers

**Key Finding**: Other operations already optimally optimized (16-25× NEON for compute-bound, scalar for memory/data-structure-bound)

---

## Current State

### Completed
- ✅ Source code analysis (Day 1)
- ✅ Performance benchmarking (Day 2)
- ✅ Algorithm deep dive (Day 3)
- ✅ GO decision made (all 5 criteria met)
- ✅ Rigorous baseline established (Entry 036)
- ✅ Applicability validated (minimizers only)
- ✅ All documentation complete and committed

### Ready to Begin
- 🎯 Phase 1: Core Implementation (ntHash + two stacks ports)

### Pending (Sequential)
- ⏳ Phase 2: Validation (Entry 036-B benchmark)
- ⏳ Phase 3: Release (v1.3.0 preparation)

---

## Decision Point: Start Phase 1 Implementation?

### Option 1: Begin Phase 1 Now ⭐ RECOMMENDED

**Tasks** (Week 1, ~20-30 hours):
1. Port ntHash from seq-hash crate (~300 LOC)
   - Adapt table lookup for NEON intrinsics style
   - Implement rolling hash state management
   - Unit tests for correctness (compare to seq-hash output)

2. Port two stacks from simd-minimizers (~200 LOC)
   - Ring buffer implementation (prefix/suffix minimums)
   - O(1) amortized sliding minimum logic
   - Property-based tests (compare to naive O(w) scan)

3. Implement block-based streaming wrapper (~100 LOC)
   - 10K block size (Rule 2 from OPTIMIZATION_RULES.md)
   - Overlap handling (k+w-1 bytes for boundary minimizers)
   - Integration with FastqStream

4. Integration testing
   - Correctness validation (matches naive implementation)
   - Performance validation (≥50× speedup minimum)
   - Memory validation (O(1) streaming maintained)

**Deliverables**:
- `src/operations/nthash.rs` (~300 LOC)
- `src/operations/sliding_min.rs` (~200 LOC)
- Updated `src/operations/kmer.rs` (~100 LOC)
- Tests passing (unit + integration + property-based)

**Timeline**: 1 week (Nov 7-14, 2025)

**Success Criteria**:
- Correctness: Output matches naive implementation exactly
- Performance: ≥50× speedup (conservative threshold)
- Memory: Constant ~5 MB regardless of sequence size
- Tests: All passing (unit, integration, property-based)

**Pros**:
- ✅ Clear path forward (algorithm understood)
- ✅ High confidence (Entry 036 baseline, GO decision)
- ✅ Low risk (~500 LOC, well-documented algorithms)
- ✅ Dramatic improvement (100-200× expected)

**Cons**:
- ⚠️ 1 week engineering effort
- ⚠️ Requires NEON adaptation (tested patterns from Entry 020-025)

**Recommendation**: **YES - Begin Phase 1 now**

### Option 2: Defer Phase 1

**Rationale for deferring**:
- ❓ Other priorities (none identified)
- ❓ Wait for additional validation (already rigorous with Entry 036)
- ❓ Uncertainty about approach (algorithm is clear)

**Cons**:
- ❌ Delays 100-200× improvement
- ❌ No active blockers or risks
- ❌ All validation complete

**Recommendation**: **NO - No reason to defer**

---

## Recommended Next Steps

### Immediate (This Week)

**1. Phase 1 Day 1: ntHash Port (Nov 7, ~8 hours)**

Start with ntHash because:
- Independent of two stacks (can develop/test separately)
- Well-documented in seq-hash crate
- Can validate correctness against seq-hash output

**Tasks**:
- [ ] Create `src/operations/nthash.rs`
- [ ] Port forward hash logic from seq-hash
- [ ] Port reverse-complement hash logic
- [ ] Implement rolling update (in_out_mapper)
- [ ] Adapt for NEON intrinsics (if beneficial)
- [ ] Unit tests (compare output to seq-hash)

**Validation criteria**:
- Output matches seq-hash for all test cases
- Handles canonical hashing correctly
- Rolling updates preserve hash values

**Phase 1 Day 2: Two Stacks Port (Nov 8, ~8 hours)**

**Tasks**:
- [ ] Create `src/operations/sliding_min.rs`
- [ ] Port ring buffer implementation
- [ ] Port prefix/suffix minimum logic
- [ ] Implement sliding minimum tracker
- [ ] Property-based tests (compare to naive O(w) scan)

**Validation criteria**:
- Output matches naive sliding minimum
- O(1) amortized performance (profile vs O(w))
- Handles edge cases (w=1, w=len, empty sequences)

**Phase 1 Day 3: Block-Based Streaming (Nov 9, ~4 hours)**

**Tasks**:
- [ ] Update `src/operations/kmer.rs`
- [ ] Implement `StreamingMinimizerExtractor`
- [ ] Block buffer with overlap handling
- [ ] Boundary minimizer merging
- [ ] Integration with FastqStream

**Validation criteria**:
- Streaming interface works with FastqStream
- Memory stays constant (profile with valgrind/heaptrack)
- Boundary minimizers handled correctly

**Phase 1 Day 4-5: Integration & Testing (Nov 10-11, ~8 hours)**

**Tasks**:
- [ ] Integration tests (end-to-end FastqStream → minimizers)
- [ ] Performance validation (≥50× speedup minimum)
- [ ] Cross-platform testing (Mac ARM, x86_64 fallback)
- [ ] Property-based tests (fuzzing with proptest)
- [ ] Documentation updates

**Validation criteria**:
- All tests passing
- Performance meets ≥50× threshold (conservative)
- Memory usage < 1 GB for human genome
- Documentation complete

### Week 2: Validation (Entry 036-B)

**Phase 2 Day 1-2: Entry 036-B Benchmark (Nov 12-13)**

**Tasks**:
- [ ] Run Entry 036-B benchmark (same 16 configs as Entry 036)
- [ ] Statistical comparison (speedup, Cohen's d, 95% CI)
- [ ] Generate comparison plots
- [ ] Document findings in Entry 036-B

**Success Validation**:
- ≥50× speedup: SUCCESS (conservative threshold) ✅
- ≥100× speedup: EXCEPTIONAL (realistic target) ✅
- ≥150× speedup: OUTSTANDING (stretch goal) ✅

**Entry 036-B deliverables**:
- ASBB entry with full statistical comparison
- Speedup validation (Entry 036-B / Entry 036)
- Cohen's d effect size (expected >> 2.0)
- 95% CI non-overlapping validation

**Phase 2 Day 3: Cross-Platform Validation (Nov 14)**

**Tasks**:
- [ ] Test on AWS Graviton (Linux ARM)
- [ ] Test x86_64 fallback (scalar implementation)
- [ ] Validate performance within 10% across ARM platforms

**Validation criteria**:
- Mac ARM: 100× speedup (optimized, Entry 036 baseline)
- Graviton: 80-100× speedup (portable NEON, Entry 021 precedent)
- x86_64: 1× (scalar fallback, no SIMD)

### Week 3: Release Preparation (v1.3.0)

**Phase 3 Day 1-2: Documentation & Python (Nov 15-16)**

**Tasks**:
- [ ] Update CHANGELOG.md (v1.3.0)
- [ ] Algorithm explanation in docs/
- [ ] Performance characteristics documentation
- [ ] Python bindings for new minimizer interface
- [ ] Usage examples

**Phase 3 Day 3: Release (Nov 17)**

**Tasks**:
- [ ] Final testing (all platforms)
- [ ] Version bump (Cargo.toml, pyproject.toml)
- [ ] Tag release (v1.3.0)
- [ ] Publish to crates.io
- [ ] Publish to PyPI
- [ ] Community announcement

---

## Risk Assessment

### Technical Risks (All LOW)

**Risk 1**: Block-based speedup lower than 100×
- **Likelihood**: Low (SIMD preserved within blocks, Entry 036 baseline rigorous)
- **Impact**: Medium (may achieve 50-80× instead of 100×)
- **Mitigation**: Still exceeds ≥50× conservative threshold
- **Fallback**: 50× speedup is still publication-quality improvement

**Risk 2**: Integration complexity underestimated
- **Likelihood**: Low (algorithms clear, ~500 LOC, MIT licensed)
- **Impact**: Low (extend Phase 1 by 2-3 days)
- **Mitigation**: Incremental development (ntHash → two stacks → streaming)
- **Fallback**: Complete Phase 1 in 10 days instead of 7

**Risk 3**: NEON adaptation challenges
- **Likelihood**: Low (ntHash operations are NEON-friendly: table lookup, rotate, XOR)
- **Impact**: Low (can use portable SIMD initially, optimize later)
- **Mitigation**: Use packed-seq patterns (proven in simd-minimizers)
- **Fallback**: Portable SIMD provides 50-75% of hand-coded NEON performance

### Strategic Risks (All LOW)

**Risk 4**: Entry 036-B validation shows <50× speedup
- **Likelihood**: Very Low (SimdMinimizers achieved 221×, we're targeting 100×)
- **Impact**: Medium (would need to investigate bottleneck)
- **Mitigation**: Incremental validation during Phase 1 (catch issues early)
- **Fallback**: Even 25× speedup would be significant improvement

**Overall Risk Level**: **LOW** (high confidence, clear path, rigorous validation)

---

## Timeline Summary

**Total Duration**: 3 weeks (Nov 7 - Nov 17, 2025)

**Week 1** (Nov 7-11): Phase 1 Implementation
- Days 1-2: ntHash + two stacks ports
- Days 3-5: Streaming wrapper + integration testing
- Deliverable: Working implementation with tests passing

**Week 2** (Nov 12-14): Phase 2 Validation
- Days 1-2: Entry 036-B benchmark (compare to Entry 036)
- Day 3: Cross-platform validation (Graviton, x86_64)
- Deliverable: Statistical validation of ≥100× speedup

**Week 3** (Nov 15-17): Phase 3 Release
- Days 1-2: Documentation + Python bindings
- Day 3: Release v1.3.0 (crates.io, PyPI)
- Deliverable: biometal v1.3.0 with 100× faster minimizers

---

## Success Metrics

### Quantitative (Entry 036-B validation)

- [ ] ≥50× speedup (conservative threshold)
- [ ] ≥100× speedup (realistic target) ← Primary
- [ ] ≥150× speedup (exceptional target)
- [ ] Memory usage < 1 GB for human genome
- [ ] O(1) memory scaling (constant for all input sizes)
- [ ] Performance within 10% across Mac ARM/Graviton

### Qualitative

- [ ] Code passes rust-code-quality-reviewer
- [ ] Documentation explains algorithm clearly
- [ ] Integration with FastqStream is seamless
- [ ] Community feedback is positive
- [ ] Publication-quality evidence (Cohen's d >> 2.0)

---

## Recommendation

### 🎯 START PHASE 1 IMMEDIATELY

**Why now**:
1. ✅ All validation complete (Entry 036 baseline, GO decision)
2. ✅ High confidence (221× measured vs 100× target)
3. ✅ Clear path (algorithm understood, ~500 LOC)
4. ✅ Low risk (MIT licensed, incremental development)
5. ✅ Dramatic impact (100-200× improvement for minimizers)

**Next action**:
- Create `src/operations/nthash.rs`
- Begin porting ntHash from seq-hash crate
- Target: Day 1 complete by Nov 7 EOD

**Expected outcome**:
- Phase 1: Working implementation (Nov 7-11)
- Phase 2: ≥100× speedup validated (Nov 12-14)
- Phase 3: biometal v1.3.0 released (Nov 17)

---

## Questions for Decision

1. **Timeline**: Is 3-week timeline acceptable? (Nov 7-17)
   - Alternative: Can extend to 4 weeks if needed
   - Recommendation: 3 weeks is realistic

2. **Scope**: Implement minimizers only? (not k-mer extraction/spectrum)
   - Evidence: APPLICABILITY_ANALYSIS.md shows minimizers only
   - Recommendation: Minimizers only (evidence-based)

3. **Testing rigor**: N=100 for Entry 036-B? (same as Entry 036)
   - Alternative: Could use N=30 (ASBB standard)
   - Recommendation: N=100 for consistency with Entry 036

4. **Release target**: v1.3.0 or v2.0.0?
   - v1.3.0: Feature addition (minimizer optimization)
   - v2.0.0: Major version (if API changes)
   - Recommendation: v1.3.0 (no breaking API changes expected)

**No blockers identified** - Ready to proceed with Phase 1 implementation.

---

**Status**: ✅ READY FOR PHASE 1
**Confidence**: HIGH (rigorous validation, clear path, low risk)
**Next Action**: Begin ntHash port (Day 1, Nov 7)
**Expected Timeline**: 3 weeks to v1.3.0 release
