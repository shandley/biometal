# ntHash + Two Stacks Applicability Analysis

**Question**: Does the 221× speedup potential from Entry 036 apply to our other operations?

**Answer**: **No - This is minimizer-specific.**

---

## Why Minimizers Achieve 221× Speedup

The speedup comes from three specific algorithmic improvements:

1. **ntHash** (rolling hash):
   - Vectorizable operations (table lookup, rotate, XOR)
   - No data dependencies like FNV-1a
   - ~8× faster hash computation

2. **Two stacks** (sliding minimum):
   - O(1) amortized (vs O(w) naive scan)
   - Each element touched ≤2 times
   - Eliminates redundant comparisons

3. **8-way SIMD parallelism**:
   - Process 8 windows simultaneously
   - Packed hash (16-bit) + position (16-bit)
   - Full utilization of NEON lanes

**Key insight**: These improvements target the **specific bottlenecks in minimizer extraction**:
- Heavy hashing workload (60% of runtime)
- Sliding window minimum finding (25% of runtime)
- Sequential processing (no parallelism)

---

## Applicability to Other Operations

### 1. Base Counting, GC Content, Quality Filter ❌ No

**Current optimization**: NEON SIMD (Rule 1, Entry 020-025)
- Base counting: **16.7× speedup** (Cohen's d = 4.82)
- Quality filter: **25.1× speedup** (Cohen's d = 5.14)

**Why ntHash + two stacks doesn't apply**:
- ❌ No hashing involved (element-wise operations)
- ❌ No sliding minimum computation
- ❌ Already faster than minimizer target (16-25× vs our 100-200× minimizer goal)

**Bottleneck**: Compute-bound → NEON vectorization is optimal

**Conclusion**: **Keep existing NEON implementation** (already superior)

---

### 2. Sequence Transforms (reverse_complement, complement, reverse) ❌ No

**Current optimization**: Scalar-only (Entry 033)
- reverse_complement: **1.03× NEON** (negligible benefit)

**Why ntHash + two stacks doesn't apply**:
- ❌ No hashing involved (simple byte transformations)
- ❌ No sliding window operations
- ❌ Memory-bound, not compute-bound

**Bottleneck**: Memory bandwidth → No optimization helps

**Conclusion**: **Keep scalar implementation** (evidence-based, Entry 033)

---

### 3. K-mer Extraction ❌ No

**Current optimization**: Scalar + optional parallel-4t (Entry 034)
- Parallel-4t: **2.19-2.38× speedup** (borderline threshold)

**Why ntHash + two stacks doesn't apply**:
- ❌ No sliding minimum computation (just extraction)
- ⚠️ ntHash could help with hashing, but Entry 035 showed **0.67-0.94× (SLOWER!)**
- ❌ Bottleneck is Vec allocations and memory operations (30-40% of runtime)

**Entry 035 evidence**:
- Tested ntHash for k-mer operations
- Result: **SLOWER** than simple FNV-1a on ASCII
- Reason: Conversion overhead + algorithm complexity

**Bottleneck**: Data-structure-bound (allocations, memory) → Parallel provides modest benefit

**Conclusion**: **Keep scalar + optional parallel-4t** (Entry 034 evidence)

---

### 4. K-mer Spectrum (Frequency Counting) ❌ No

**Current optimization**: Scalar-only (Entry 034)
- Parallel: **0.95-1.88× (INCONSISTENT, sometimes SLOWER!)**

**Why ntHash + two stacks doesn't apply**:
- ❌ No sliding minimum computation
- ⚠️ ntHash could help with hashing (50-60% of runtime)
- ❌ Bottleneck is HashMap updates (30-40% of runtime, sequential)

**Critical finding** (Entry 034):
- Parallelization causes **HashMap contention**
- Thread merging causes cache thrashing
- Result: Sometimes **slower** than scalar

**Bottleneck**: HashMap-bound (sequential data structure) → No parallelization helps

**Conclusion**: **Keep scalar-only** (Entry 034 evidence - parallel makes it worse!)

---

### 5. Minimizers ✅ YES (Only applicable operation!)

**Current optimization**: Scalar-only (Entry 034, Entry 036 baseline)
- Entry 034 (pilot): 1.02-1.26× max with NEON/Parallel
- Entry 036 (rigorous): 1.7-5.5 Mbp/s baseline

**Why ntHash + two stacks applies**:
- ✅ Heavy hashing workload (60% of runtime) → ntHash helps
- ✅ Sliding window minimum finding (25% of runtime) → Two stacks helps
- ✅ Sequential processing → 8-way SIMD parallelism helps

**Evidence**:
- SimdMinimizers: 820 Mbp/s (full SIMD)
- Our baseline: 3.7 Mbp/s (Entry 036)
- Potential: **221× speedup**

**Bottleneck**: Hash + sliding-minimum-bound → ntHash + two stacks + SIMD is optimal

**Conclusion**: **Integrate ntHash + two stacks** (GO decision, Entry 036 validates)

---

## Evidence Summary

| Operation | Current | Speedup | ntHash + Two Stacks? | Rationale |
|-----------|---------|---------|----------------------|-----------|
| Base counting | NEON | 16.7× | ❌ No | Already faster, no hashing/sliding-min |
| Quality filter | NEON | 25.1× | ❌ No | Already faster, no hashing/sliding-min |
| reverse_complement | Scalar | 1.03× | ❌ No | Memory-bound, no improvement possible |
| K-mer extraction | Scalar + opt. parallel | 2.2× | ❌ No | Entry 035: ntHash is SLOWER (0.67-0.94×) |
| K-mer spectrum | Scalar | 1× | ❌ No | HashMap-bound, parallel makes it worse |
| **Minimizers** | **Scalar** | **1×** | ✅ **YES** | **221× potential, GO decision made** |

---

## Strategic Recommendation

### Selective Optimization (Evidence-Based)

**DO implement ntHash + two stacks for**:
- ✅ Minimizer extraction ONLY

**DO NOT retrofit to**:
- ❌ Base counting, GC content, quality filter (already 16-25× with NEON)
- ❌ Sequence transforms (memory-bound, Entry 033 evidence)
- ❌ K-mer extraction (Entry 035: ntHash is slower)
- ❌ K-mer spectrum (HashMap-bound, parallel is worse)

### Why Selective Optimization Is Correct

**1. Different bottlenecks require different solutions**:
- Compute-bound (base counting) → NEON vectorization ✅
- Memory-bound (reverse_complement) → No optimization helps ✅
- Data-structure-bound (k-mer extraction) → Modest parallel benefit ✅
- HashMap-bound (k-mer spectrum) → Scalar is optimal ✅
- Hash+sliding-min-bound (minimizers) → ntHash + two stacks ✅

**2. Evidence validates selective approach**:
- Entry 020-025: NEON provides 16-25× for element-wise ops
- Entry 033: NEON provides 1.03× for sequence transforms (negligible)
- Entry 034: Parallel provides 2.2× for k-mer extraction, 0.95× for spectrum
- Entry 035: ntHash is 0.67-0.94× for simple k-mers (SLOWER!)
- Entry 036: Baseline is 221× slower than SimdMinimizers for minimizers

**3. Complexity vs benefit**:
- ntHash port: ~300 LOC
- Two stacks port: ~200 LOC
- Total: ~500 LOC

Only justified where evidence shows ≥5× improvement (Rule from Phase 4).

**4. Aligns with biometal's core principle**:
> "Evidence-based optimization: Every rule validated experimentally"

We don't optimize for the sake of optimization - we optimize where evidence proves benefit.

---

## Answer to User's Question

> "Will we need to adopt this to all of our previously coded operations?"

**No.** The 221× speedup is **minimizer-specific** because:

1. **Other operations are already optimized appropriately**:
   - Base counting: 16.7× NEON (Entry 020-025)
   - Quality filter: 25.1× NEON (Entry 020-025)
   - These are already faster than our minimizer target!

2. **ntHash + two stacks doesn't apply to different bottlenecks**:
   - Element-wise ops: Already solved by NEON (Rule 1)
   - Memory-bound ops: No optimization helps (Entry 033)
   - Data-structure-bound ops: ntHash is slower (Entry 035)
   - HashMap-bound ops: Scalar is optimal (Entry 034)

3. **Evidence says selective optimization is correct**:
   - Entry 035 explicitly tested ntHash for k-mers: SLOWER
   - Entry 033 tested NEON for transforms: negligible (1.03×)
   - Entry 034 tested parallel for spectrum: worse (contention)

4. **We should implement ntHash + two stacks ONLY for minimizers**:
   - This is where the 221× opportunity exists
   - Clear evidence from Entry 036 baseline
   - GO decision already made with high confidence

---

## Conclusion

**The right strategy is selective optimization**:
- ✅ Minimizers: Integrate ntHash + two stacks (100-200× projected)
- ✅ Base operations: Keep NEON (16-25× proven)
- ✅ K-mer extraction: Keep scalar + optional parallel (2.2× proven)
- ✅ K-mer spectrum: Keep scalar (parallel is worse)
- ✅ Sequence transforms: Keep scalar (NEON doesn't help)

**This is evidence-based design at its best** - optimize where evidence proves benefit, don't over-engineer where it doesn't help.

**biometal remains**:
- **Fast where it matters** (minimizers: 100-200×, base counting: 16×, quality filter: 25×)
- **Simple where appropriate** (transforms: scalar, spectrum: scalar)
- **Evidence-driven** (every optimization validated experimentally)

---

**Recommendation**: Proceed with Phase 1 implementation for **minimizers only**. Do not retrofit ntHash + two stacks to other operations.
