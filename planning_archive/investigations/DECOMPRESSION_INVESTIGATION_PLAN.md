# Decompression Library Investigation Plan

**Date**: November 11, 2025
**Goal**: Validate 2-3× decompression speedup to address actual bottleneck (98.7% of time)
**Effort**: 20-30 hours investigation
**Decision**: Go/No-Go based on validated results

---

## Discovery: flate2 Supports zlib-ng Backend!

**Current State**:
- Using `flate2 = "1.1.5"` with default `rust_backend` feature (pure Rust `miniz_oxide`)
- Decompression: 4.37s for 544 MB file (98.7% of total time)

**Available backends in flate2**:
- ✅ `default = [rust_backend]` - Pure Rust miniz_oxide (CURRENT)
- ✅ `zlib-ng` - High-performance zlib-ng backend (TARGET #1)
- ✅ `zlib` - Standard system zlib
- ✅ `cloudflare_zlib` - Cloudflare's optimized zlib
- ✅ `zlib-rs` - Pure Rust zlib reimplementation

**Key Insight**: Can test zlib-ng with ZERO code changes - just modify Cargo.toml!

---

## Investigation Strategy

### Phase 1: flate2 with zlib-ng Backend (PRIORITY - Easiest)

**Effort**: 2-4 hours
**Code changes**: NONE (just Cargo.toml feature flag)
**Risk**: Very low

**Steps**:
1. Modify `Cargo.toml` to use `flate2` with `zlib-ng` feature
2. Re-run existing benchmarks (rule4_mmap_validation)
3. Compare decompression time: Current 4.37s vs zlib-ng
4. **Decision**: If >1.5× speedup → DONE (integrate), otherwise proceed to Phase 2

**Expected Outcome**:
- If zlib-ng provides 2-3× speedup: Investigation complete!
- If zlib-ng provides <1.5× speedup: Try libdeflate (Phase 2)

---

### Phase 2: libdeflate Integration (IF Phase 1 insufficient)

**Effort**: 15-20 hours
**Code changes**: Moderate (different API, streaming adaptation)
**Risk**: Medium (integration complexity)

**Steps**:
1. Create `libdeflate` wrapper for BGZF decompression
2. Benchmark with N=30 samples
3. Compare vs zlib-ng results
4. Assess integration complexity
5. **Decision**: If combined benefit >1.5× → Integrate best performer

**Challenges**:
- libdeflate API is different from flate2 (buffer-based, not streaming)
- Need to adapt to biometal's streaming architecture
- May conflict with Rule 5 (constant memory streaming)

---

### Phase 3: Alternative Backends (IF Phase 1+2 insufficient)

**Effort**: 5-10 hours
**Options**:
- Cloudflare's zlib (via flate2 feature)
- zlib-rs (pure Rust, via flate2 feature)
- Intel ISA-L igzip (x86_64 only, not relevant for M4 Max)

**Decision**: Test if no clear winner from Phases 1+2

---

## Benchmark Plan

### Test Matrix

**Files** (same as Rule 4 investigation):
- Small: 5.4 MB (`large_100k_150bp.fq.gz`)
- Medium: 54 MB (`vlarge_1m_150bp.fq.gz`)
- Large: 544 MB (`huge_10m_150bp.fq.gz`)

**Backends to test**:
1. **Baseline**: flate2 with `rust_backend` (current)
2. **Target #1**: flate2 with `zlib-ng` (easiest)
3. **Target #2**: libdeflate (if needed)
4. **Optional**: flate2 with `cloudflare_zlib` (if needed)

**Metrics**:
- Decompression time (primary)
- Total time (I/O + decompression)
- Memory usage (ensure no regression)
- Cross-platform implications

**Sample Size**: N=30 per configuration (statistical rigor)

---

## Success Criteria

### Go Decision (Integrate)
- ✅ Decompression speedup ≥1.5× (validated with N=30)
- ✅ No memory regression (stays ~5 MB constant)
- ✅ Works on ARM (M4 Max validation)
- ✅ Integration effort reasonable (<40 hours)
- ✅ Cross-platform compatible (Mac ARM, Linux ARM, x86_64)

### No-Go Decision (Accept Current Performance)
- ❌ Decompression speedup <1.5× (DAG pruning threshold)
- ❌ Memory regression (violates Rule 5)
- ❌ Integration too complex (>60 hours effort)
- ❌ Platform-specific (doesn't work on all targets)

---

## Implementation Plan (IF Go)

### Week 1-2: Investigation (Current)
1. **Day 1**: Phase 1 testing (flate2 + zlib-ng)
2. **Day 2**: Phase 2 testing (libdeflate) IF needed
3. **Day 3**: Analysis and decision
4. **Day 4**: Documentation

### Week 3-4: Integration (IF Go)
1. Update Cargo.toml dependencies
2. Validate on M4 Max (Mac ARM)
3. Test on Linux ARM (GitHub Actions Graviton)
4. Test on x86_64 (GitHub Actions ubuntu-latest)
5. Memory profiling (ensure no regression)
6. Update documentation

---

## Expected Outcomes

### Optimistic (2-3× speedup achieved)
- Decompression: 4.37s → 1.5-2.2s
- Total time: 4.42s → 1.5-2.2s
- Overall speedup: **2.0-2.9×**
- BAM parsing: 55 MiB/s → **110-160 MiB/s**
- **Status**: EXCEEDS original Phase 2A target (16.3× was unrealistic, but 2-3× is real!)

### Realistic (1.5-2× speedup achieved)
- Decompression: 4.37s → 2.2-2.9s
- Total time: 4.42s → 2.2-2.9s
- Overall speedup: **1.5-2.0×**
- BAM parsing: 55 MiB/s → **82-110 MiB/s**
- **Status**: Solid improvement, worth integrating

### Pessimistic (< 1.5× speedup)
- Decompression: 4.37s → >2.9s
- Overall speedup: **<1.5×**
- **Status**: DAG prune, accept current performance, move to horizontal expansion

---

## Risk Analysis

### Phase 1 Risks (flate2 + zlib-ng)

**Low Risk**:
- ✅ Zero code changes (just Cargo.toml)
- ✅ Same API (drop-in replacement)
- ✅ Quick to test (2-4 hours)
- ✅ Easy to revert if doesn't work

**Unknowns**:
- ⚠️ ARM NEON optimizations in zlib-ng (might be x86_64 focused)
- ⚠️ macOS compatibility (should work but untested)
- ⚠️ Cross-platform build complexity

### Phase 2 Risks (libdeflate)

**Medium Risk**:
- ⚠️ API differences (buffer-based vs streaming)
- ⚠️ Integration effort (15-20 hours)
- ⚠️ May conflict with streaming architecture (Rule 5)
- ⚠️ Maintenance burden (two decompression paths)

**Mitigation**:
- Only pursue if Phase 1 insufficient
- Thorough testing before integration
- Clear documentation of trade-offs

---

## Timeline

### Immediate (Today)
- [x] Research available libraries ✓
- [x] Discover flate2 supports zlib-ng ✓
- [x] Create investigation plan ✓
- [ ] Phase 1: Test flate2 + zlib-ng (IN PROGRESS)

### Day 1-2 (November 11-12)
- [ ] Phase 1: Benchmark zlib-ng with N=30
- [ ] Phase 1: Analyze results
- [ ] Phase 1: Decision (Go/No-Go/Try Phase 2)

### Day 3-4 (November 13-14) - IF Phase 2 needed
- [ ] Phase 2: Benchmark libdeflate with N=30
- [ ] Phase 2: Assess integration complexity
- [ ] Phase 2: Final decision

### Day 5-7 (November 15-17) - IF Go
- [ ] Integration and cross-platform testing
- [ ] Documentation
- [ ] Commit and push

---

## Deliverables

### Investigation Phase (This Week)
1. **Benchmark results**: N=30 samples per backend
2. **Analysis document**: `DECOMPRESSION_INVESTIGATION_FINDINGS.md`
3. **Decision**: Go/No-Go with evidence
4. **Updated roadmap**: Reflect decision

### Integration Phase (Next Week) - IF Go
1. **Code changes**: Cargo.toml (minimal) or integration (moderate)
2. **Cross-platform validation**: Mac ARM, Linux ARM, x86_64
3. **Performance validation**: 2-3× speedup confirmed
4. **Documentation**: OPTIMIZATION_RULES.md update

---

## Next Actions

**Immediate**:
1. Modify Cargo.toml to test flate2 + zlib-ng backend
2. Re-run rule4_mmap_validation benchmark
3. Compare decompression times
4. Quick decision: Does zlib-ng provide ≥1.5×?

**If Yes**:
- Document findings
- Plan integration
- Test cross-platform

**If No**:
- Proceed to Phase 2 (libdeflate)
- OR accept current performance
- OR try other backends

---

## Document Status

**Created**: November 11, 2025
**Purpose**: Investigation plan for decompression bottleneck
**Status**: Phase 1 starting
**Next Update**: After Phase 1 benchmark results

**Key Decision**: Start with easiest option (flate2 + zlib-ng) before complex integration

---

**Remember**: We're addressing the actual bottleneck (98.7% decompression time). Even 1.5× improvement here is worth more than 10× improvement on I/O (1.3% of time). Amdahl's Law favors optimizing the dominant portion!
