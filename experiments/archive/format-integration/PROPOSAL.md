# Format Integration Experiment: BAM Wrapper Prototype

**Date**: November 7, 2025
**Status**: Proposal
**Type**: Technical Validation Experiment

---

## Hypothesis

A lightweight wrapper around noodles BAM reader can provide:
1. Unified API consistent with biometal patterns
2. Performance overhead < 5%
3. Maintained streaming characteristics
4. Path to cross-format operations

## Research Questions

1. **Overhead**: What performance penalty does wrapping introduce?
2. **API Design**: Can we create intuitive, consistent API?
3. **Integration**: How well does noodles integrate with our error handling?
4. **Streaming**: Is constant memory maintained?

## Success Criteria

### Quantitative
- ✅ Wrapper compiles and runs
- ✅ Performance overhead < 5% vs direct noodles
- ✅ Memory usage equivalent to noodles
- ✅ Error handling integrates cleanly

### Qualitative
- ✅ API feels natural and consistent
- ✅ Documentation is clear
- ✅ Integration path is obvious

## Go/No-Go Decision

### GO if:
- Overhead < 5%
- API design is clean
- No major integration issues

### NO-GO if:
- Overhead > 10%
- API requires complex workarounds
- Integration breaks streaming

## Timeline

- **Day 1** (Nov 7): Prototype implementation
  - Add noodles dependency
  - Create RecordReader trait
  - Implement BAM wrapper
  - Basic validation

- **Day 2** (Nov 8): Benchmarking and validation
  - Performance comparison (wrapper vs direct)
  - Memory profiling
  - Integration testing
  - Decision point

## Expected Outcome

**Hypothesis**: Wrapper adds < 5% overhead

**Evidence needed**:
- Benchmark: Read 100K BAM records (wrapper vs direct noodles)
- Profile: Memory usage during streaming
- Code review: API ergonomics

## Deliverables

- [x] `src/io/bam.rs` - BAM wrapper implementation
- [x] Cargo.toml - noodles dependencies (0.67/0.64/0.33)
- [x] Example code - `examples/bam_reader_example.rs`
- [x] Test validation - compilation and API tests pass
- [ ] Real BAM file testing (needs test data)
- [ ] Benchmark results (overhead measurement)
- [ ] GO/NO-GO decision document

## Progress Update (Day 1 - Nov 7, 2025)

### ✅ Completed

**Implementation** (src/io/bam.rs):
- Created BamReader<R> wrapper with streaming Iterator API
- Integrated biometal error handling (BiometalError)
- Generic over Read trait (works with files, network streams, etc.)
- Escape hatch methods: `inner_mut()`, `into_inner()` for advanced usage
- Type-safe wrapper over noodles BGZF-compressed BAM reader

**Dependencies** (Cargo.toml):
- noodles-bam 0.67 (BAM I/O)
- noodles-sam 0.64 (header types, matched to noodles-bam)
- noodles-bgzf 0.33 (BGZF compression)
- Optional "formats" feature flag (no bloat for default builds)

**Validation**:
- ✅ Compiles successfully
- ✅ Tests pass (1 API test)
- ✅ Example compiles (demonstrates usage patterns)
- ✅ Clean API consistent with FastqStream

### 🔍 Observations

**API Design** (Qualitative Success Criterion):
- Iterator pattern matches FastqStream/FastaStream ✅
- Consistent error handling via Result<Record> ✅
- Clean from_path() / from_reader() constructors ✅
- No panic-inducing unwraps ✅

**Integration** (Success Criterion):
- noodles integrates cleanly with biometal types ✅
- Version resolution straightforward (match noodles-bam's dependencies) ✅
- BGZF reader composition handled via type alias ✅

### ⏭️ Next Steps (Day 2)

**Performance Validation** (Quantitative):
- [ ] Create benchmark: wrapper vs direct noodles
- [ ] Measure overhead with 100K records
- [ ] Target: < 5% overhead (GO threshold)
- [ ] Blocker: < 10% overhead (NO-GO threshold)

**Real-World Testing**:
- [ ] Obtain test BAM file (small, ~1MB)
- [ ] Verify record iteration works end-to-end
- [ ] Test error handling with malformed input
- [ ] Memory profiling (constant usage validation)

**Decision Point**:
- If overhead < 5%: GO → Document and recommend expansion
- If 5-10%: CONDITIONAL → Identify optimization opportunities
- If > 10%: NO-GO → Consider alternative approach or direct noodles

---

## Final Status (Day 3 - Nov 7, 2025)

### ✅ Experiment Complete

**Decision**: **NO-GO** for wrapper, **GO** for direct noodles integration

**Phase 1 Implementation** (Completed):
- ✅ Removed experimental wrapper code
- ✅ Cleaned up dependencies
- ✅ Created comprehensive integration guide (docs/NOODLES_INTEGRATION.md)
- ✅ Build verified

**Deliverables**:
- [x] Comprehensive analysis (FINDINGS.md, OPTION_D_ANALYSIS.md, ECOSYSTEM_ANALYSIS.md)
- [x] Integration guide (docs/NOODLES_INTEGRATION.md)
- [x] Evidence-based recommendation
- [x] Clean codebase (wrapper removed)

**Time Investment**: 8 hours total (Day 1-3)
**ROI**: Excellent - validated approach early, prevented weeks of wasted effort

---

**Final Status**: Complete - Integration strategy implemented
