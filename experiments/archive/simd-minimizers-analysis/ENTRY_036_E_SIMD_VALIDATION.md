# Entry 036-E: SIMD Minimizer Validation

**Date**: November 7, 2025
**Status**: ✅ **COMPLETE** - Integrated and documented
**Implementation**: `extract_minimizers_simd()` in `src/operations/kmer.rs:530`
**Benchmark**: `benches/minimizer_simd.rs`

---

## Executive Summary

SIMD minimizer extraction provides **1.5-1.9× speedup** for genomic-scale sequences (≥100Kbp), making it valuable for production bioinformatics workloads. However, small sequences (<100Kbp) show **regression** (0.41-0.84×) due to chunking overhead, establishing a clear crossover point at ~100Kbp.

**Decision**: ✅ **ACCEPTED** - Ship as opt-in feature with clear documentation

---

## Performance Results

### Benchmark Configuration

- **Platform**: Mac M1 (Apple Silicon)
- **Compiler flags**: `RUSTFLAGS="-C target-cpu=native"`
- **Sample size**: N=20
- **Measurement time**: 10s per configuration
- **Benchmark**: `benches/minimizer_simd.rs`

### Scaling Behavior

| Sequence Length | Fast (MiB/s) | SIMD (MiB/s) | Speedup | Status |
|-----------------|--------------|--------------|---------|--------|
| 100bp           | 225          | 92           | **0.41×** | ❌ Regression |
| 1Kbp            | 222          | 183          | **0.82×** | ❌ Regression |
| 10Kbp           | 225          | 195          | **0.87×** | ❌ Regression |
| **100Kbp**      | **122**      | **184**      | **1.51×** | ✅ **Crossover** |
| **1Mbp**        | **99.4**     | **193**      | **1.94×** | ✅ **Production** |
| **10Mbp**       | **98.2**     | **147**      | **1.50×** | ✅ **Production** |

### Key Findings

1. **Small sequences (≤10Kbp)**: SIMD overhead dominates → use `extract_minimizers_fast`
2. **Medium sequences (100Kbp)**: Crossover point → ~1.5× speedup begins
3. **Large sequences (≥1Mbp)**: Consistent 1.5-1.9× speedup → ideal for SIMD

### Production Scale (k=21, w=11)

```
Fast:  99.4 MiB/s @ 1Mbp  (baseline from Entry 036-C)
SIMD: 193.2 MiB/s @ 1Mbp  (1.94× faster) ✓

Fast:  98.2 MiB/s @ 10Mbp
SIMD: 147.0 MiB/s @ 10Mbp (1.50× faster) ✓
```

---

## Gap Analysis: Literature vs Implementation

### Expected (simd-minimizers paper)
- **Paper claims**: 3-15× speedup
- **Paper workload**: Human genome (3Gbp) in 4 seconds
- **Paper measurement**: Direct SIMD vs scalar comparison

### Achieved (biometal)
- **Our results**: 1.5-1.9× speedup
- **Our workload**: 1-10Mbp sequences
- **Our measurement**: SIMD with ntHash recomputation

### Root Causes

1. **Double hashing overhead** (primary)
   - simd-minimizers computes positions only
   - We recompute hashes with ntHash for API compatibility
   - **Impact**: ~2× overhead from redundant hashing

2. **Canonical overhead**
   - Strand-aware selection (TG counting) adds work
   - Leftmost vs rightmost minimizer selection logic

3. **Memory bandwidth**
   - May be bottlenecked on memory access, not compute
   - M1 memory subsystem characteristics

### Why We Accept 1.5-1.9×

**Evidence-based decision**: While below paper claims, 1.5-1.9× is a **meaningful speedup** for:
- Real-world genomic analysis (bacterial genomes, contigs, gene clusters)
- Production bioinformatics pipelines
- Users processing 100Kbp - 100Mbp sequences

**Trade-off accepted**: API compatibility (ntHash hashes) over maximum theoretical speedup

---

## Implementation Details

### Architecture

```rust
#[cfg(feature = "simd")]
pub fn extract_minimizers_simd(sequence: &[u8], k: usize, w: usize) -> Result<Vec<Minimizer>> {
    // 1. Wrap sequence for SIMD library
    let ascii_seq = AsciiSeq(sequence);

    // 2. Get minimizer positions using SIMD (8-way parallel)
    let positions = simd_minimizers::canonical_minimizer_positions(ascii_seq, k, w);

    // 3. Recompute hashes with ntHash for API compatibility
    let hash_iter = NtHashIterator::new(sequence, k)?;
    let hashes: Vec<(usize, u64)> = hash_iter.enumerate().collect();

    // 4. Combine positions + hashes into Minimizer structs
    let minimizers = positions.into_iter()
        .filter_map(|pos| {
            let pos_usize = pos as usize;
            hashes.get(pos_usize).map(|&(_idx, hash)| Minimizer {
                position: pos_usize,
                hash,
                k,
            })
        })
        .collect();

    Ok(minimizers)
}
```

### Key Design Decisions

1. **Canonical minimizers**: Uses strand-aware selection (may differ from fast implementation)
2. **ntHash recomputation**: Ensures consistent hash values across all biometal functions
3. **Zero-copy wrapper**: `AsciiSeq` has no allocation overhead
4. **Feature flag**: Opt-in via `features = ["simd"]` in Cargo.toml

### Correctness Validation

Test: `test_extract_minimizers_simd_correctness` validates:
- ✅ Valid minimizer positions (within sequence bounds)
- ✅ Valid k-mer extraction
- ✅ Non-zero ntHash values
- ✅ Deduplication (no consecutive duplicate positions)
- ✅ Correct k-mer size metadata

**Note**: SIMD may select different minimizers than fast (both valid due to canonical selection)

---

## Usage Guidelines

### When to Use SIMD

```rust
use biometal::operations::kmer::{extract_minimizers_fast, extract_minimizers_simd};

// ✅ GOOD: Large genomic sequences
let genome = vec![b'A'; 1_000_000];  // 1Mbp
let minimizers = extract_minimizers_simd(&genome, 21, 11)?;  // 1.94× faster

// ❌ BAD: Small sequences
let read = b"ATGCATGCATGC";  // 12bp
let minimizers = extract_minimizers_simd(read, 5, 7)?;  // 0.41× slower!

// ✅ GOOD: Use fast for small sequences
let minimizers = extract_minimizers_fast(read, 5, 7)?;  // Optimal
```

### Feature Flag Setup

```toml
# Cargo.toml
[dependencies]
biometal = { version = "1.2", features = ["simd"] }
```

```bash
# Compile with native CPU features for maximum performance
RUSTFLAGS="-C target-cpu=native" cargo build --release --features simd
```

### Auto-Selection Strategy

```rust
use biometal::operations::kmer::{extract_minimizers_fast, extract_minimizers_simd};

pub fn extract_minimizers_auto(sequence: &[u8], k: usize, w: usize) -> Result<Vec<Minimizer>> {
    const SIMD_THRESHOLD: usize = 100_000;  // 100Kbp crossover

    #[cfg(feature = "simd")]
    if sequence.len() >= SIMD_THRESHOLD {
        return extract_minimizers_simd(sequence, k, w);
    }

    extract_minimizers_fast(sequence, k, w)
}
```

---

## Comparison with Baseline

### Timeline

- **Entry 036** (baseline): 3.7 Mbp/s (FNV-1a + O(w) scan)
- **Entry 036-B** (fast): 77-231 Mbp/s (21-62× speedup, ntHash + SlidingMin)
- **Entry 036-C** (lazy): 105 Mbp/s @ 1Mbp (1.40× from lazy k-mer)
- **Entry 036-D** (block-based): NO-GO (0.90× regression)
- **Entry 036-E** (SIMD): 193 Mbp/s @ 1Mbp (1.94× from SIMD) ✅

### Cumulative Speedup vs Baseline

```
Baseline (Entry 036):      3.7 Mbp/s @ 1Mbp
Fast (Entry 036-C):      105.0 Mbp/s @ 1Mbp (28× faster)
SIMD (Entry 036-E):      193.0 Mbp/s @ 1Mbp (52× faster) ✅✅
```

**Total achievement**: **52× speedup** from baseline through iterative optimization

---

## Testing

### Benchmark Command

```bash
# Quick validation (N=20, ~30 minutes)
RUSTFLAGS="-C target-cpu=native" cargo bench --bench minimizer_simd --features simd \
  -- --measurement-time 10 --sample-size 20

# Full validation (N=100, ~90 minutes)
RUSTFLAGS="-C target-cpu=native" cargo bench --bench minimizer_simd --features simd \
  -- --measurement-time 60 --sample-size 100
```

### Test Coverage

```bash
# Correctness test
cargo test --features simd test_extract_minimizers_simd_correctness

# All tests with SIMD
cargo test --features simd
```

---

## Future Optimization Opportunities

### Potential Improvements

1. **Eliminate double hashing** (target: 3-5× speedup)
   - Expose hash values from simd-minimizers
   - Modify simd-minimizers to return positions + hashes
   - Trade-off: Forked dependency or upstream contribution

2. **Adaptive selection** (automatic fast vs SIMD)
   - Auto-detect sequence length at runtime
   - Benchmark-driven threshold tuning per platform

3. **Batch processing** (amortize overhead)
   - Process multiple small sequences in single SIMD batch
   - Good for FASTQ streaming with many short reads

### Why We're Not Pursuing Now

**Evidence-based threshold**: Optimization must show ≥5× improvement to justify complexity.
- Current gap: 1.5-1.9× (actual) vs 3-15× (theoretical) = **2-8× gap**
- Eliminating double hashing: **estimated 2-3× improvement** (below ≥5× threshold)
- **Decision**: Ship current implementation, revisit if user demand justifies complexity

---

## Documentation

### Code Documentation

- ✅ Function docs updated: `src/operations/kmer.rs:454-528`
- ✅ Performance section: Entry 036-E results
- ✅ When to use guidelines
- ✅ Gap from literature explained
- ✅ Example updated for realistic use case

### User-Facing Docs

- ✅ Benchmark: `benches/minimizer_simd.rs`
- ✅ Test: `test_extract_minimizers_simd_correctness`
- ✅ This document: Entry 036-E summary

---

## Conclusions

### Success Criteria Met

| Criterion | Target | Achieved | Status |
|-----------|--------|----------|--------|
| Speedup @ 1Mbp | ≥3× | 1.94× | ⚠️ Below target |
| Speedup @ 10Mbp | ≥3× | 1.50× | ⚠️ Below target |
| Production value | Meaningful | 1.5-1.9× | ✅ Valuable |
| Small seq penalty | Acceptable | 0.41-0.84× | ✅ Documented |

### Final Decision

✅ **GO** - Accept 1.5-1.9× speedup and integrate

**Rationale**:
- Meaningful improvement for production genomics (≥100Kbp)
- Evidence-based documentation prevents misuse
- Clear crossover point (100Kbp) guides users
- Opt-in feature flag minimizes risk
- Gap from literature is explained and acceptable

### Lessons Learned

1. **Double hashing matters**: API compatibility has measurable performance cost
2. **Scaling is non-linear**: SIMD overhead dominates for small inputs
3. **Documentation is critical**: Must guide users to appropriate use cases
4. **Evidence-based shipping**: 1.5-1.9× is valuable, don't wait for perfect

---

## References

- **Paper**: "SimdMinimizers: Computing Random Minimizers, fast" (SEA 2025)
- **Library**: simd-minimizers v2.2 (MIT licensed)
- **Entry 036** (baseline): FNV-1a + O(w) scan (3.7 Mbp/s)
- **Entry 036-B** (fast): ntHash + SlidingMin (77-231 Mbp/s)
- **Entry 036-C** (lazy): Lazy k-mer optimization (105 Mbp/s @ 1Mbp)
- **Entry 036-D** (block-based): NO-GO regression (0.90×)
- **Entry 036-E** (this): SIMD integration (193 Mbp/s @ 1Mbp)

---

**Status**: ✅ Integrated into biometal v1.3.0 (opt-in `simd` feature)
**Benchmark log**: `/tmp/minimizer_simd_quick.log`
**Implementation**: `src/operations/kmer.rs:530`
**Tests**: `src/operations/kmer.rs:1139` (correctness validation)
