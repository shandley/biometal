# Experiment: SIMD Minimizers Performance Analysis

**Status**: Research Phase
**Start Date**: November 6, 2025
**Owner**: Scott Handley
**Timeline**: 1 week (time-boxed)

---

## Hypothesis

SimdMinimizers' rolling hash + SIMD argmin approach can achieve 4-9.5× speedup for minimizer extraction, contradicting our Entry 034 findings (1.02-1.26× with NEON).

---

## Background

### Our Current State (Entry 034)
- **Implementation**: FNV-1a hash + HashMap-based minimizer extraction
- **NEON speedup**: 1.02-1.26× (below ≥5× threshold)
- **Conclusion**: Data-structure-bound (HashMap dominates 50-60% of runtime)
- **Decision**: Scalar-only implementation

### SimdMinimizers Claims (bioRxiv Jan 2025)
- **Paper**: "SimdMinimizers: Computing random minimizers, fast"
- **Authors**: Ragnar Groot Koerkamp, Igor Martayan (ETH Zurich, ENS Rennes)
- **Implementation**: Rust library `rust-seq/simd-minimizers`
- **Performance**:
  - 9.5× speedup for w=5 (small window)
  - 4.5× speedup for w=19 (large window)
  - 16× faster than existing `minimizer-iter` crate
  - Human genome (3.2 Gbp): 4.1 seconds (forward), 6.0 seconds (canonical)

### The Discrepancy

**Critical Question:** How did they achieve 9.5× when we got 1.26×?

**Hypotheses:**
1. **Different algorithm**: Rolling hash (vectorizable) vs FNV-1a (sequential)
2. **Decoupled operations**: SIMD-friendly parts (hash, argmin) separated from data structures
3. **Batch processing**: Compute all hashes first, then find minimizers
4. **Memory layout**: Cache-friendly access patterns
5. **Our testing was flawed**: We may have coupled too many operations

---

## Research Questions

### Primary Questions
1. **What algorithm do they use?** Rolling hash? Polynomial hash? Other?
2. **How is SIMD applied?** Which operations are vectorized?
3. **What's the data flow?** When does HashMap/deduplication happen?
4. **Can we replicate their speedup?** Is it architecture-specific or generalizable?

### Secondary Questions
5. **What's the memory footprint?** Constant or scales with sequence length?
6. **Does it work on ARM NEON?** Or AVX2-only?
7. **What are the limitations?** k/w parameter ranges, sequence properties?
8. **How does it compare to minimap2?** minimap2 uses scalar minimizers

---

## Success Criteria (Go/No-Go)

### GO Criteria (Pursue Integration)
- [ ] **Speedup**: Can demonstrate ≥4× speedup with rolling hash + SIMD approach
- [ ] **Constant memory**: Maintains streaming architecture (no O(n) buffers)
- [ ] **ARM compatible**: Works on Apple Silicon with NEON
- [ ] **Understandable**: Clear algorithmic advantage we can implement
- [ ] **Evidence-based**: Can experimentally validate speedup (Entry 035)

### NO-GO Criteria (Stay with Current)
- [ ] **Architecture-specific**: Only works on AVX2/AVX-512, not NEON
- [ ] **Memory trade-off**: Requires O(n) buffering (breaks streaming)
- [ ] **Complex**: Implementation complexity >> benefit
- [ ] **Marginal improvement**: <3× speedup (not worth rewrite)
- [ ] **Incorrect comparison**: Their benchmark doesn't match our use case

---

## Methodology

### Phase 1: Code Analysis (Days 1-2)
**Objective**: Understand their implementation

1. **Clone repository**
   ```bash
   git clone https://github.com/rust-seq/simd-minimizers
   cd simd-minimizers
   ```

2. **Study source code**
   - Core algorithm: `src/lib.rs`
   - SIMD implementation: Look for `#[cfg(target_arch)]` blocks
   - Benchmark code: `benches/` directory
   - Test suite: Understand correctness validation

3. **Read paper thoroughly**
   - bioRxiv: https://www.biorxiv.org/content/10.1101/2025.01.27.634998v1
   - Blog post: https://curiouscoding.nl/posts/simd-minimizers/
   - SEA 2025 proceedings: Full technical details

4. **Document architecture**
   - Data structures used
   - Control flow diagram
   - Memory access patterns
   - SIMD instruction usage

### Phase 2: Benchmarking (Day 3)
**Objective**: Validate their performance claims

1. **Install and build**
   ```bash
   cargo build --release
   cargo test --release
   ```

2. **Run their benchmarks**
   ```bash
   cargo bench
   ```

3. **Custom benchmarks** (comparable to Entry 034)
   - Same input sequences (E. coli, random DNA)
   - Same k/w parameters (k=15, w=10 typical)
   - Same platform (Mac M-series)
   - Measure: throughput, memory usage

4. **Compare to our implementation**
   - Load Entry 034 data
   - Run side-by-side comparison
   - Analyze discrepancy sources

### Phase 3: Algorithm Analysis (Days 4-5)
**Objective**: Understand the technique

1. **Identify SIMD-friendly operations**
   - Which operations use SIMD intrinsics?
   - Where is the 9.5× speedup coming from?
   - What's scalar vs vectorized?

2. **Trace data flow**
   ```
   Input sequence → [Step 1] → [Step 2] → ... → Minimizers

   Identify:
   - Where hashing happens (SIMD?)
   - Where minimum finding happens (SIMD?)
   - Where deduplication happens (HashMap?)
   - Memory allocation points
   ```

3. **Compare to our Entry 034 approach**
   ```
   Our approach:
   sequence → sliding windows → FNV-1a hash → HashMap insert

   Their approach:
   sequence → ??? → ??? → minimizers

   Find the difference!
   ```

4. **Understand rolling hash**
   - What's the formula?
   - Why is it vectorizable?
   - How does it compare to FNV-1a?

### Phase 4: Feasibility Assessment (Days 6-7)
**Objective**: Determine if we should integrate

1. **NEON compatibility check**
   - Do they have ARM NEON implementation?
   - If not, can we port AVX2 → NEON?
   - Complexity estimate

2. **Integration complexity**
   - Lines of code to add
   - Breaking changes required?
   - Test coverage needed

3. **Performance projection**
   - If we port to biometal, expected speedup?
   - Memory footprint impact?
   - Streaming compatibility?

4. **Decision matrix**
   | Factor | Weight | Score | Weighted |
   |--------|--------|-------|----------|
   | Speedup | 40% | ? | ? |
   | NEON support | 30% | ? | ? |
   | Complexity | 20% | ? | ? |
   | Memory | 10% | ? | ? |
   | **Total** | 100% | | **?** |

   **GO if total score ≥ 70/100**

---

## Deliverables

### Research Phase Documents
1. **RESEARCH_LOG.md** - Daily findings, experiments, insights
2. **ARCHITECTURE.md** - Their implementation architecture
3. **COMPARISON.md** - Side-by-side with our Entry 034
4. **BENCHMARK_RESULTS.md** - Performance data

### Decision Document
5. **FINDINGS.md** - Final analysis and go/no-go recommendation

### If GO Decision
6. **IMPLEMENTATION_PLAN.md** - Detailed integration roadmap
7. **Entry 035** (ASBB) - Experimental validation of new approach

### If NO-GO Decision
6. **RATIONALE.md** - Why we're staying with current approach
7. Update Entry 034 with additional context

---

## Timeline

**Week 1: Research & Decision**
- Day 1-2: Code analysis, understand algorithm
- Day 3: Benchmark validation, performance comparison
- Day 4-5: Deep algorithmic analysis, identify technique
- Day 6-7: Feasibility assessment, go/no-go decision

**Total time commitment**: 1 week (5-7 days)

**Decision deadline**: November 13, 2025

---

## Risks & Mitigations

### Risk 1: Their speedup is AVX2-specific
**Mitigation**: Focus on understanding the *algorithmic* advantage, not just SIMD. Even if we can't match 9.5×, understanding rolling hash may improve our scalar implementation.

### Risk 2: Apples-to-oranges comparison
**Mitigation**: Run identical benchmarks with same input data. Document any differences in methodology.

### Risk 3: Time sink (over 1 week)
**Mitigation**: Strict time-box. If not clear by Day 7, make decision with available data.

### Risk 4: Implementation too complex
**Mitigation**: Complexity is a valid NO-GO criterion. Document finding and move on.

---

## Expected Outcomes

### Best Case (GO)
- Understand technique completely
- 4-9× speedup achievable on ARM NEON
- Clear integration path
- Entry 035 planned for full validation
- v1.3.0 feature

### Good Case (GO with caveats)
- Understand technique
- 2-4× speedup achievable (worth it)
- Moderate complexity
- Phased integration approach

### Acceptable Case (NO-GO, learned something)
- Understand why they got better speedup
- Identify what we did right in Entry 034
- Document edge case or use-case difference
- No action needed, current approach validated

### Worst Case (NO-GO, wasted time)
- Can't replicate their speedup
- Implementation is impenetrable
- No clear learnings
- 1 week lost

**Probability assessment**:
- Best/Good case: 60-70%
- Acceptable case: 20-30%
- Worst case: 5-10%

---

## Success Metrics

### Quantitative
- [ ] Replicate their claimed speedup (±20%)
- [ ] Identify specific SIMD operations providing speedup
- [ ] Benchmark on our hardware (Mac M-series)
- [ ] Measure memory footprint

### Qualitative
- [ ] Understand algorithmic difference from Entry 034
- [ ] Assess NEON portability
- [ ] Evaluate integration complexity
- [ ] Make confident go/no-go decision

---

## References

**SimdMinimizers**:
- Paper: https://www.biorxiv.org/content/10.1101/2025.01.27.634998v1
- GitHub: https://github.com/rust-seq/simd-minimizers
- Blog: https://curiouscoding.nl/posts/simd-minimizers/
- SEA 2025: https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2025.20

**Our Work**:
- Entry 034: `apple-silicon-bio-bench/lab-notebook/2025-11/20251106-034-EXPERIMENT-kmer-operations.md`
- K-mer module: `src/operations/kmer.rs`
- Minimizer benchmarks: `benches/kmer.rs`

**Background**:
- Farrar (2007): Striped SIMD (alignment context)
- minimap2: Uses scalar minimizers (validation point)
- Parasail: SIMD pairwise alignment (related technique)

---

## Notes

- This experiment follows evidence-based methodology
- Clear go/no-go criteria prevent scope creep
- Time-boxed to 1 week (learned from sra-decoder)
- Hypothesis is testable and specific
- Outcome documented regardless of result

**Key principle**: It's okay to say NO-GO if evidence doesn't support integration. The point is to make an *informed* decision.

---

**Status**: Ready to begin
**Next action**: Clone repository and start Day 1 code analysis
