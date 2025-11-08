# SRA Decoder Experiment

**Status**: Research (Day 1 of 14)
**Hypothesis**: ARM NEON-optimized SRA decoder can achieve ≥10× speedup vs scalar
**Timeline**: Nov 5-15, 2025 (2-week research phase)

---

## Quick Overview

This experiment investigates whether we can build a native ARM-optimized decoder for NCBI's Sequence Read Archive (SRA) format that outperforms the existing SRA Toolkit.

**Why?** SRA stores petabases of sequencing data, but requires a complex C++ toolkit to decode. An ARM-native streaming decoder would:
- Enable 5 MB memory footprint (vs downloading 50-200 MB files)
- Provide 10-40× speedup with NEON SIMD
- Integrate cleanly with biometal's streaming architecture
- Democratize access for LMIC researchers

**Risk**: SRA format may be too complex or change too frequently. **Time-boxed to 2 weeks** to validate feasibility.

---

## Current Status

### Week 1 Progress

- [x] **Day 1**: Experiment setup, proposal complete
- [ ] **Day 2-3**: Format analysis, scalar baseline
- [ ] **Day 4-6**: NEON optimization
- [ ] **Day 7-8**: Benchmarking (N=30)
- [ ] **Day 9-10**: Validation and analysis
- [ ] **Day 11-12**: Go/No-Go decision

### Key Metrics (To Be Measured)

| Metric | Target | Current | Status |
|--------|--------|---------|--------|
| Speedup (NEON vs scalar) | ≥10× | TBD | 🟡 Pending |
| Memory footprint | ~5 MB | TBD | 🟡 Pending |
| Correctness | 100% | TBD | 🟡 Pending |
| Code maintainability | Clean | TBD | 🟡 Pending |

---

## Technical Approach

### Phase 1: Scalar Baseline (Days 1-3)
1. Study SRA format specification
2. Implement scalar base unpacking (2-bit → 8-bit)
3. Implement quality score decompression
4. Validate against SRA Toolkit output

### Phase 2: NEON Optimization (Days 4-8)
1. NEON-optimized base unpacking (hot path)
2. Profile for additional optimization opportunities
3. Block-based processing (10K records)
4. Comprehensive benchmarks (N=30)

### Phase 3: Decision (Days 9-12)
1. Statistical analysis (mean, 95% CI, Cohen's d)
2. Correctness validation (1000+ reads)
3. Make Go/No-Go decision
4. Document findings

---

## Hot Paths Identified

### 1. Base Unpacking (2-bit → 8-bit) - PRIMARY TARGET
- **Operation**: Convert packed bases (00=A, 01=C, 10=G, 11=T) to ASCII
- **NEON potential**: Process 32 bases per instruction
- **Expected speedup**: 16-25× (Rule 1)
- **CPU time**: ~40% of total

### 2. Quality Score Decompression
- **Operation**: Inflate compressed quality strings
- **Approach**: Use existing flate2 (already optimized)
- **CPU time**: ~30% of total

### 3. Read Name Reconstruction
- **Priority**: LOWER (I/O dominated)
- **CPU time**: ~30% of total

---

## Success Criteria

### ✅ GO Decision (proceed to full implementation)
- Speedup: **≥10×** (NEON vs scalar)
- Memory: **Constant ~5 MB**
- Correctness: **100% match** with SRA Toolkit
- Integration: **Clean API** matching FastqStream

### 🤔 MAYBE Decision (deeper investigation)
- Speedup: 5-10× (marginal benefit)
- Some edge cases fail (scope assessment needed)

### ❌ NO-GO Decision (wrap existing toolkit)
- Speedup: **<5×** (not worth custom implementation)
- Correctness: Fails on common inputs
- Complexity: Format too complex to reverse-engineer

---

## Development

### Build & Test

```bash
# From experiments/sra-decoder/
cargo build
cargo test
cargo bench
```

### Run Benchmarks (When Ready)

```bash
# Run full benchmark suite (N=30)
cargo bench -- --sample-size 30

# View results
open target/criterion/report/index.html
```

### Validate Correctness (When Ready)

```bash
# Compare against SRA Toolkit
fastq-dump SRR390728 > reference.fastq
cargo run --example validate -- SRR390728 reference.fastq
```

---

## Files

- `PROPOSAL.md` - Full experiment proposal ⭐ START HERE
- `RESEARCH_LOG.md` - Daily progress log (updated daily)
- `FINDINGS.md` - Final results (filled out at end)
- `src/lib.rs` - Scalar baseline + dispatch
- `src/neon.rs` - ARM NEON optimizations
- `benches/benchmark.rs` - Performance benchmarks

---

## Resources

**SRA Format:**
- [NCBI SRA Format Docs](https://github.com/ncbi/sra-tools/wiki)
- [sra-tools source code](https://github.com/ncbi/sra-tools)
- [SRA Handbook](https://www.ncbi.nlm.nih.gov/books/NBK47540/)

**Evidence Base:**
- [OPTIMIZATION_RULES.md](../../OPTIMIZATION_RULES.md) - 6 rules from 1,357 experiments
- [ASBB Lab Notebook](https://github.com/shandley/apple-silicon-bio-bench)

**Test Data:**
- SRR390728 (E. coli, 195 MB, well-characterized)
- SRR1553425 (Human, longer reads)

---

## Publication Plan

### If Successful (≥10× speedup)
- **Target**: Bioinformatics journal or ISMB/RECOMB conference
- **Title**: "ARM-Native SRA Decoder: Democratizing Access to Petabases of Sequencing Data"
- **Timeline**: Submit Q1 2026

### If Moderately Successful (5-10× speedup)
- **Target**: bioRxiv preprint + blog post
- **Value**: Documents approach for others

### If Unsuccessful (<5× speedup)
- **Target**: Technical report
- **Value**: Prevents others from wasting time on dead-end
- **Still valuable**: Shows evidence-based process works

---

## Contact

**Owner**: Scott Handley
**Status**: Research Phase
**Decision Date**: Nov 15, 2025

---

**Remember**: This is research. Failure is an option. We learn from both success and failure.
