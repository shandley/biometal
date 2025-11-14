# Apple Silicon Research (November 4-13, 2025)

**Period**: November 4-13, 2025 (2 weeks)
**Location**: research/apple-silicon/
**Status**: ✅ Archived - Returning to core roadmap

## What Was Built

**Neural Engine (Week 1)**:
- ✅ Complete ONNX Runtime + CoreML integration
- ✅ Quality prediction model (PyTorch → ONNX → CoreML)
- ❌ Result: 2,940× **slowdown** for streaming use case
- **Finding**: Neural Engine optimized for batch inference, not per-read streaming

**GPU Smith-Waterman (Week 2)**:
- ✅ Metal compute shader (340 lines)
- ✅ Rust GPU dispatch (430 lines)
- ✅ All tests passing (430 tests)
- ⚠️ Result: 1.2-1.4× speedup for batches ≥10
- **Finding**: Modest gains vs 10-50× from literature (needs anti-diagonal parallelization)

## Strategic Decision (November 13, 2025)

**Decision**: Archive Apple Silicon research

**Rationale**:
- Apple Silicon (batch-oriented) vs biometal's streaming-first architecture = poor fit
- Modest results: Neural Engine (2,940× slowdown), GPU (1.2-1.4× speedup)
- Platform-specific (Mac-only) breaks cross-platform promise
- **Better path**: Focus on cross-platform, evidence-based features

**Critical Discovery (Same Day)**:
After archiving Apple Silicon, reviewed OPTIMIZATION_RULES.md and discovered:
- ❌ Rule 3 (Parallel BGZF): **DISABLED** (0.77-0.84× slowdown, not 6.5× speedup)
- ⚠️ Rule 4 (Smart mmap): **~1% benefit** (not 2.5×, CPU-bound workload)
- **Result**: Original Phase 2 plan (Rules 3+4 for 16×) is invalid

**See**: RULES_3_4_REALITY_CHECK.md for full analysis, STRATEGIC_PIVOT_PLAN.md for options

## Lessons Learned

1. **Hardware-software fit matters**: Neural Engine/GPU excel at batch processing, not streaming
2. **Read documentation carefully**: Summary tables can be outdated, check detailed sections
3. **Context-dependent optimizations**: Entry 029's 6.5× works for unbounded memory, not streaming
4. **Platform-specific code has costs**: Breaks portability, smaller user base
5. **Amdahl's Law applies**: 2.5× I/O speedup on 1.3% of time = ~1% overall

**Value Preserved**:
- ✅ Infrastructure reusable for future use cases (adapter detection, variant calling)
- ✅ Code quality demonstrates capabilities (all production-ready)
- ✅ GPU available via feature flag (`--features gpu`)

**Status**: Archived November 13, 2025. See research/apple-silicon/RESEARCH_SUMMARY.md for full analysis.
