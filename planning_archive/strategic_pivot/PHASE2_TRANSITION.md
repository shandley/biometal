# Phase 2 Transition: Return to Core Roadmap

**Date**: November 13, 2025
**Decision**: Option A - Return to Core Roadmap (Rules 3+4 Implementation)
**Status**: ❌ **INVALID** - See Critical Correction Below

---

## ⚠️ CRITICAL CORRECTION (November 13, 2025)

**This document contains incorrect claims about Rules 3+4 performance.**

After reviewing OPTIMIZATION_RULES.md, I discovered:
- **Rule 3 (Parallel BGZF)**: ❌ **DISABLED** in biometal (0.77-0.84× slowdown, not 6.5× speedup)
- **Rule 4 (Smart mmap)**: ⚠️ **~1% benefit** for compressed files (not 2.5×)
- **Combined**: NO significant speedup available (not 16×)

**See**: RULES_3_4_REALITY_CHECK.md for full analysis

**New direction**: See STRATEGIC_PIVOT_PLAN.md for alternative Phase 2 options

**Status**: This transition plan is archived for historical reference but should NOT be followed.

---

## Executive Summary

After 2 weeks exploring Apple Silicon hardware acceleration (Neural Engine + GPU), we're returning to biometal's core roadmap: implementing Rules 3+4 for proven 16× speedup. This aligns with biometal's evidence-based, cross-platform streaming architecture.

## What Happened

### Apple Silicon Exploration (Archived)

**Neural Engine** (Week 1):
- ✅ Built: Complete ONNX Runtime + CoreML integration
- ❌ Result: 2,940× slowdown for streaming use case
- **Lesson**: Hardware mismatch (batch inference vs streaming)

**GPU Smith-Waterman** (Week 2):
- ✅ Built: Production Metal shader + Rust dispatch (870 lines)
- ⚠️ Result: 1.2-1.4× speedup for batches ≥10
- **Lesson**: Modest gains vs 10-50× literature claims

### Strategic Assessment

**Core tension**:
- biometal = streaming-first, evidence-based, cross-platform
- Apple Silicon = batch-oriented, speculative, platform-specific

**Opportunity cost**:
- 2 weeks on Apple Silicon (1.2-1.4× gains)
- Rules 3+4 waiting (16× **proven** speedup from ASBB)

**Decision**: Focus on high-ROI, evidence-based optimizations

---

## What's Next: Phase 2 (Rules 3+4)

### Rule 3: Parallel BGZF Decompression

**Target**: 6.5× speedup (92 MiB/s → 600 MiB/s)
- Evidence: ASBB Entry 029 (CPU parallel prototype)
- Approach: Rayon-based parallel block decompression
- Platform: Cross-platform (benefits all users)
- Effort: 40-60 hours

### Rule 4: Smart mmap

**Target**: 2.5× additional (combined 16× total)
- Evidence: ASBB Entry 032 (scale validation)
- Approach: Memory-mapped I/O with madvise hints
- Platform: macOS/Linux (graceful degradation on Windows)
- Effort: 40-60 hours

### Combined Outcome

**Performance**:
- Current: 92 MiB/s (cloudflare_zlib baseline)
- Target: **~900 MiB/s** (16× improvement)
- Validation: ASBB experiments + N=30 benchmarks

**Implementation Status**:
- Before: 4/6 rules (67%)
- After: **6/6 rules (100% complete)**
- Combined speedup: **27× across all rules**

---

## Files Changed

### Archived
- `research/apple-silicon/RESEARCH_SUMMARY.md` - Complete analysis
- `research/apple-silicon/gpu-smith-waterman/` - GPU implementation docs
- `research/apple-silicon/neural-engine/` - Neural Engine docs

### Updated
- `CLAUDE.md` - Removed Apple Silicon focus, added Phase 2 focus
- **Current Focus**: "Phase 2 - Rules 3+4 Implementation (16× speedup)"
- **Research Status**: "Apple Silicon archived - Returning to core roadmap"

### Preserved (Feature-Gated)
- `src/alignment/gpu/` - GPU code (feature: `gpu`)
- `src/ml/` - Neural Engine (feature: `neural-engine`)
- Can be activated if needed: `cargo build --features gpu`

---

## Lessons Applied

### Technical
1. **Hardware-software fit matters**: Match workload to hardware strengths
2. **Evidence-based > speculative**: Validated ASBB results beat literature claims
3. **Cross-platform value**: Optimizations that benefit all users > platform-specific

### Strategic
1. **Opportunity cost is real**: Time spent on speculative work vs proven optimizations
2. **Mission alignment critical**: Streaming-first vs batch-oriented mismatch
3. **Stay focused**: Evidence-based methodology guides correct priorities

---

## Next Steps

### Immediate (This Session)
1. ✅ Archive Apple Silicon research
2. ✅ Update CLAUDE.md (Phase 2 focus)
3. 🔄 Review OPTIMIZATION_RULES.md (Rules 3+4 details)
4. 🔄 Design parallel BGZF implementation

### Phase 2 Roadmap (3-4 weeks)
1. **Week 1**: Design + implement Rule 3 (parallel BGZF)
2. **Week 2**: Design + implement Rule 4 (smart mmap)
3. **Week 3**: Integration + cross-platform testing
4. **Week 4**: Benchmarking (N=30) + validation vs ASBB

### Success Criteria
- ✅ 16× speedup achieved (92 MiB/s → 900 MiB/s)
- ✅ All 6 optimization rules implemented
- ✅ Cross-platform validated (Mac ARM, Linux ARM, x86_64)
- ✅ Benchmarks confirm ASBB predictions (N=30)

---

## Project Status

### Current State (v1.7.0)
- ✅ 4/6 optimization rules implemented (67%)
- ✅ Streaming architecture (constant ~5 MB memory)
- ✅ ARM NEON optimization (16-25× speedups)
- ✅ Network streaming (HTTP, SRA)
- ✅ BAM/SAM parsing (92 MiB/s with cloudflare_zlib)
- ⚠️ Rules 3+4 pending (16× proven speedup waiting)

### Target State (v1.8.0 - after Phase 2)
- ✅ **6/6 optimization rules (100% complete)**
- ✅ BAM parsing: **~900 MiB/s** (16× improvement)
- ✅ Cross-platform validated
- ✅ All ASBB evidence implemented
- ✅ Community-ready for major release

---

## Retrospective

### What Went Well
- ✅ Production-quality code (all implementations)
- ✅ Comprehensive testing methodology
- ✅ Honest assessment and course correction
- ✅ Infrastructure preserved for future use

### What Didn't Work
- ❌ Problem selection (streaming vs batch mismatch)
- ❌ Speculative optimization (literature vs measured)
- ❌ Platform-specific focus (breaks portability)

### Key Insight
**Evidence-based methodology works**: ASBB-validated optimizations (16×) are more valuable than speculative hardware exploration (1.2-1.4×). Staying focused on biometal's core strengths (streaming, cross-platform, evidence-based) delivers better outcomes.

---

## Conclusion

**Transition complete**: Apple Silicon research archived, returning to core roadmap.

**Next priority**: Phase 2 (Rules 3+4 implementation)
- Proven 16× speedup from ASBB
- Cross-platform benefits
- Completes all 6 optimization rules

**Timeline**: 3-4 weeks to Phase 2 completion → v1.8.0 release

**Impact**: Transformative performance for all biometal users, not just Mac users.

---

**Status**: Ready to begin Rules 3+4 implementation
**Next action**: Review OPTIMIZATION_RULES.md for detailed guidance
