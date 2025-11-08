# Format Integration Experiment: Research Log

**Experiment**: BAM Wrapper Prototype
**Start Date**: November 7, 2025
**Status**: Day 1 Complete

---

## Day 1: Prototype Implementation (Nov 7, 2025)

### Objective

Create lightweight wrapper around noodles-bam to validate hybrid integration strategy.

### Approach

1. **Dependency Selection**: Use noodles-bam 0.67 (pure Rust, well-maintained)
2. **API Design**: Iterator pattern matching FastqStream/FastaStream
3. **Feature Flag**: Optional "formats" feature to avoid bloating default build
4. **Error Handling**: Integrate with BiometalError for consistency

### Implementation Details

#### Initial Challenges

**Challenge 1: noodles API Discovery**
- **Issue**: Initial uncertainty about correct noodles API patterns
- **Investigation**: Used WebFetch to examine noodles-bam and noodles-sam docs
- **Finding**: noodles uses `.records()` iterator pattern + `read_record()` method
- **Decision**: Implement Iterator trait manually using `read_record()` for full control

**Challenge 2: Type Resolution**
- **Issue**: `Reader::new(R)` wraps input in BGZF reader → type mismatch
- **Error**: `expected Reader<R>, found Reader<Reader<R>>`
- **Solution**: Created type alias `BamReaderInner<R> = BamIoReader<BgzfReader<R>>`
- **Lesson**: BAM files are inherently BGZF-compressed; wrapper must account for this

**Challenge 3: Version Conflicts**
- **Issue**: Direct dependency on noodles-sam 0.66 conflicted with noodles-bam's 0.64
- **Error**: `two different versions of crate noodles_sam`
- **Investigation**: `cargo tree` revealed noodles-bam uses 0.64.0
- **Solution**: Match exact version (noodles-sam 0.64) to avoid conflicts
- **Lesson**: Let transitive dependencies guide version selection for wrapper crates

#### Final Implementation

**File Structure**:
```
src/io/bam.rs           # BamReader<R> implementation
examples/bam_reader_example.rs  # Usage demonstration
Cargo.toml              # Dependencies + formats feature
```

**Dependencies Added**:
```toml
noodles-bam = { version = "0.67", optional = true }
noodles-sam = { version = "0.64", optional = true }  # Match noodles-bam
noodles-bgzf = { version = "0.33", optional = true }
```

**API Pattern**:
```rust
pub struct BamReader<R: Read> {
    reader: BamReaderInner<R>,  // = BamIoReader<BgzfReader<R>>
    header: noodles_sam::Header,
}

impl<R: Read> Iterator for BamReader<R> {
    type Item = Result<Record>;
    // Streaming iteration, constant memory
}
```

### Results

#### ✅ Success Criteria Met

**Quantitative** (Partial):
- ✅ Wrapper compiles and runs
- ⏳ Performance overhead < 5% (needs benchmarking)
- ⏳ Memory usage equivalent to noodles (needs profiling)
- ✅ Error handling integrates cleanly

**Qualitative**:
- ✅ API feels natural and consistent (Iterator pattern)
- ✅ Documentation is clear (examples provided)
- ✅ Integration path is obvious (feature flag + re-exports)

#### Code Quality

**Compilation**: Clean build with `--features formats`
```
Finished `dev` profile in 0.82s
```

**Tests**: 1 API test passes
```
test io::bam::tests::test_bam_reader_api_compiles ... ok
```

**Example**: Compiles with 1 harmless warning (unused import in demo code)

### Key Decisions

**Decision 1: Generic over Read trait**
- **Rationale**: Enables file, network stream, memory buffer usage
- **Tradeoff**: Slightly more complex types vs maximum flexibility
- **Outcome**: Consistent with FastqStream design

**Decision 2: Escape hatches for advanced usage**
- **Methods**: `inner_mut()`, `into_inner()`
- **Rationale**: Users can drop to noodles for operations not exposed by wrapper
- **Benefit**: Zero-cost abstraction - no performance penalty for advanced users

**Decision 3: Feature-gated dependencies**
- **Rationale**: BAM support is niche; don't bloat default builds
- **Implementation**: `formats` feature flag
- **Benefit**: Users opt-in to format support

### Observations

#### What Worked Well

1. **noodles Integration**: Smooth, no major API friction
2. **Type System**: Rust's type system caught mismatches early (version conflicts, Reader wrapping)
3. **Cargo Ecosystem**: `cargo tree` invaluable for dependency resolution
4. **Documentation**: noodles docs sufficient for integration

#### What Could Be Improved

1. **Version Discovery**: Trial-and-error for matching versions (could be faster with better tooling)
2. **Type Complexity**: BGZF reader wrapping adds cognitive load (but unavoidable for BAM)
3. **Testing**: Need real BAM files for end-to-end validation

### Lessons Learned

1. **Wrapper Pattern Works**: Minimal code (~150 LOC) for significant value (unified API)
2. **Version Matching Critical**: Must match transitive dependency versions exactly
3. **Type Aliases Help**: Simplify complex generic types (BamReaderInner)
4. **Escape Hatches Essential**: Don't force users through wrapper for all operations

### Next Actions (Day 2)

**Priority 1: Performance Validation** (Quantitative)
- [ ] Create benchmark comparing wrapper vs direct noodles
- [ ] Measure with 100K BAM records
- [ ] Target: < 5% overhead
- [ ] **GO/NO-GO Decision Point**

**Priority 2: Real-World Testing**
- [ ] Obtain small test BAM file (~1MB)
- [ ] Verify end-to-end iteration
- [ ] Test error handling with malformed input
- [ ] Memory profiling

**Priority 3: Documentation**
- [ ] Document findings in FINDINGS.md
- [ ] Update FILE_FORMAT_INTEGRATION_ANALYSIS.md
- [ ] Recommendation for next formats (VCF, GFF, BED)

### Time Tracking

- **Planning**: 30 min (PROPOSAL.md creation)
- **Dependency Setup**: 15 min (Cargo.toml, version resolution)
- **Implementation**: 90 min (API design, type resolution, error handling)
- **Testing/Validation**: 30 min (tests, example, compilation)
- **Documentation**: 45 min (PROPOSAL updates, example documentation)
- **Total**: ~3.5 hours (well within Day 1 budget)

### Preliminary Assessment

**Hypothesis Validation**: ✅ Wrapper approach appears viable

**Evidence**:
- Clean API integration (Iterator pattern works)
- Minimal code (~150 LOC for full wrapper)
- No performance measurement yet, but no obvious bottlenecks
- User experience matches biometal patterns (FastqStream consistency)

**Confidence**: 75% GO (pending benchmark validation)

**Risk**: If overhead > 5%, may need optimization pass or direct noodles recommendation

---

## Next Session Prep

**Environment**:
- ✅ Code compiles cleanly
- ✅ Tests pass
- ⏳ Need BAM test file for benchmarking
- ⏳ Need to create benchmark harness

**Questions to Answer**:
1. What is the actual performance overhead?
2. Does memory usage stay constant during iteration?
3. How does error handling perform with malformed BAM?
4. Should we proceed to VCF/GFF/BED wrappers?

**Files to Review**:
- `src/io/bam.rs` (implementation)
- `examples/bam_reader_example.rs` (usage patterns)
- `experiments/format-integration/PROPOSAL.md` (objectives)

---

## Day 2: Performance Validation (Nov 7, 2025)

### Objective

Measure wrapper performance overhead and make GO/NO-GO decision.

### Approach

1. **Benchmark Infrastructure**: Created criterion-based benchmarks
2. **Synthetic Data**: Generated minimal BAM files (100, 1K, 10K records)
3. **Three Benchmark Groups**:
   - `bam_wrapper`: biometal wrapper throughput
   - `bam_direct`: noodles baseline throughput
   - `bam_overhead_comparison`: Head-to-head comparison

### Implementation

**Challenge**: Initial noodles Writer API issues
- Error: `write_alignment_record` doesn't exist → use `write_record`
- Error: `finish` method requires trait import → use `try_finish()`
- Solution: Simplified to minimal record generation (empty records valid for overhead testing)

**Benchmark Design**:
```rust
// Baseline (direct noodles)
let mut reader = noodles_bam::io::Reader::new(cursor);
let mut record = Record::default();
while reader.read_record(&mut record)? != 0 { count += 1; }

// Wrapper (biometal)
let reader = BamReader::from_reader(cursor)?;
for result in reader {
    let _record = result?;
    count += 1;
}
```

### Results

#### Quantitative Performance

| Records | Baseline (µs) | Wrapper (µs) | Overhead | Status |
|---------|--------------|-------------|----------|--------|
| 100     | 7.72         | 8.71        | **+12.8%** | ❌ NO-GO |
| 1,000   | 31.77        | 49.83       | **+56.9%** | ❌ NO-GO |
| 10,000  | 255.67       | 360.11      | **+40.8%** | ❌ NO-GO |

**Throughput Impact**:
- 100: 12.96 → 11.48 Melem/s (-11.4%)
- 1K: 31.48 → 20.07 Melem/s (-36.3%)
- 10K: 39.11 → 27.77 Melem/s (-29.0%)

**Statistical Confidence**: N=100 samples per benchmark, minimal outliers, consistent across runs

#### Success Criteria Evaluation

- ❌ **Target**: < 5% overhead (GO) → FAILED (12.8-56.9%)
- ❌ **Acceptable**: 5-10% overhead (CONDITIONAL) → FAILED
- ❌ **Blocker**: > 10% overhead (NO-GO) → **TRIGGERED**

**Decision**: **CLEAR NO-GO** for naive wrapper approach

### Root Cause Analysis

**Primary Bottleneck: Record Allocation**

```rust
fn next(&mut self) -> Option<Self::Item> {
    let mut record = Record::default();  // ← Allocates EVERY iteration
    match self.reader.read_record(&mut record) { ... }
}
```

**Architecture Mismatch**:
- **noodles pattern**: Pre-allocate record, reuse buffer (zero-allocation)
- **Iterator pattern**: Yields owned records (allocation per iteration)
- **Conflict**: Can't have Iterator AND zero-copy simultaneously

**Secondary Factors**:
- Error wrapping: `BiometalError` creation + `format!()` macro
- Double wrapping: `Option<Result<Record>>`
- Minor: Virtual dispatch through Iterator trait

**Conclusion**: Iterator API ownership semantics fundamentally conflict with record reuse pattern.

### Key Findings

1. **Performance**: 40% overhead unacceptable for production use
2. **API**: Wrapper API is elegant and consistent (qualitative success)
3. **Integration**: noodles integrates cleanly (no technical blockers)
4. **Architecture**: Iterator pattern inappropriate for record reuse libraries

### Lessons Learned

#### 1. API Elegance ≠ Performance

**Trade-off**:
- ✅ Iterator API: Ergonomic, consistent, familiar
- ❌ 40% overhead: Unacceptable for bioinformatics workloads

**Lesson**: Sometimes direct library usage is better than "unified" wrapper.

#### 2. Benchmark Early = Huge ROI

**Time Investment**:
- Day 1 (prototype): 3.5 hours
- Day 2 (benchmark): 2 hours
- **Total**: 5.5 hours to NO-GO decision

**Alternative Timeline**:
- Without benchmarking: 2-3 weeks of integration before discovering overhead
- **ROI**: ~25:1 (time saved by early validation)

**Conclusion**: Experiment process worked perfectly - caught dead-end early.

#### 3. Mature Libraries Are Fast

**Observation**: noodles achieves 38-39 Melem/s throughput at 10K records.

**Implication**: Hard to improve on specialized, mature libraries.

**Strategy**: Integrate, don't replicate.

### Decision Matrix

**Three Options Identified**:

**Option A: Direct noodles Usage** (RECOMMENDED)
- ✅ Zero overhead (by definition)
- ✅ Full API access
- ✅ Less maintenance
- ❌ API inconsistency with biometal

**Option B: Optimized Wrapper**
- ✅ Could achieve 2-5% overhead with record reuse API
- ❌ More complex (dual API pattern)
- ❌ 1-2 weeks additional development
- ⚠️ Needs validation

**Option C: "Convenience API" Wrapper**
- ✅ Useful for prototyping
- ❌ 40% overhead is a liability
- ❌ Reputation risk for "performance" library

### Revised Strategy

**Original Hypothesis**: Wrap all formats for unified API

**Revised Strategy** (evidence-based):

**Tier 1** (Performance-Critical): Direct usage
- BAM/SAM/CRAM: Use noodles directly
- Provide: Documentation, examples, re-exports

**Tier 2** (Moderate Overhead OK): Lightweight wrappers
- VCF/GFF/BED: Wrapper acceptable (less performance-sensitive)
- Provide: Wrapper + escape hatches

**Tier 3** (Native Implementation): biometal-specific
- Formats with NEON potential: Native implementation
- Example: Quality score parsing with SIMD

### Next Actions

**Immediate (Day 3)**:
1. Choose Option A, B, or C (recommendation: **A**)
2. Update FILE_FORMAT_INTEGRATION_ANALYSIS.md
3. Document "Using noodles with biometal" patterns
4. Remove or mark experimental: src/io/bam.rs

**Short-term (Week 1-2)**:
1. Evaluate VCF wrapper potential
2. Research GFF/BED integration
3. Create noodles integration examples

### Time Tracking

**Day 2 Breakdown**:
- Benchmark creation: 60 min
- API fixes (Writer issues): 30 min
- Running benchmarks: 15 min
- Analysis + documentation: 45 min
- **Total Day 2**: 2.5 hours

**Cumulative**: 6 hours (Day 1 + Day 2)

### Experimental Process Evaluation

**What Worked**:
- ✅ Time-boxed approach caught NO-GO early
- ✅ Clear quantitative thresholds (5%/10%)
- ✅ Comprehensive benchmarking (multiple sizes, N=100)
- ✅ Complete documentation (PROPOSAL → RESEARCH_LOG → FINDINGS)

**What Could Improve**:
- Earlier API research (record reuse pattern)
- Prototype both patterns (Iterator vs buffer reuse) on Day 1
- Memory profiling (deferred due to early NO-GO)

**Overall Assessment**: **Excellent** - Experiment worked as designed, delivered clear decision with evidence.

---

**Status**: Day 2 Complete - NO-GO Decision Made
**Next**: Finalize format integration strategy (Day 3)
