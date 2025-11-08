# SIMD Acceleration Roadmap for biometal

**Date**: November 6, 2025
**Status**: Strategic Planning Document
**Purpose**: Chart path for expanding biometal's SIMD capabilities

---

## Current State Assessment

### ✅ What We Have (Excellent Foundation)

**Element-Wise Operations** - World-class NEON performance:
- Base counting: **16.7× speedup** (Entry 020, Cohen's d = 4.82)
- GC content: **20.3× speedup** (Entry 021)
- Quality filtering: **25.1× speedup** (Entry 022)
- Complexity scoring: **18.2× speedup** (Entry 023)

**Technical Excellence:**
- Clean NEON intrinsics implementation
- Portable scalar fallbacks
- Evidence-based design (1,357 experiments)
- Production-ready (347 tests)

**Coverage Map:**
```
SIMD Operations Spectrum
│
├─ Element-wise (16-25×) ✅ DONE
│  ├─ Base counting
│  ├─ GC content
│  ├─ Quality scores
│  └─ Complexity
│
├─ Data-structure-bound (1-2×) ✅ DONE (Scalar-only)
│  ├─ K-mer extraction
│  ├─ Minimizers
│  └─ K-mer spectrum
│
└─ Sequence alignment (??×) ❌ NOT STARTED
   ├─ Local (Smith-Waterman)
   ├─ Global (Needleman-Wunsch)
   └─ Semi-global
```

### ⚠️ What We're Missing (Significant Gaps)

**1. Alignment Operations** - The Big Gap
- No Smith-Waterman (local alignment)
- No Needleman-Wunsch (global alignment)
- No semi-global alignment
- **This is where Parasail dominates**

**2. Advanced K-mer Operations** - Recent Innovations
- No SIMD minimizers (SimdMinimizers just published Jan 2025)
- No advanced k-mer schemes (mod-minimizers, syncmers)
- Our minimizers are scalar-only

**3. Pattern Matching** - Unexplored
- No SIMD string search
- No motif finding acceleration
- No regex-like operations

---

## Competitive Landscape: SIMD in Genomics (2024-2025)

### SimdMinimizers (Jan 2025) - **Direct Competition**

**What They Built:**
- Rust library (`rust-seq/simd-minimizers`)
- AVX2 + NEON support (cross-platform)
- Random minimizer algorithm with SIMD

**Performance:**
- **9.5× faster** than scalar for w=5 (small window)
- **4.5× faster** for w=19 (large window)
- **16× faster** than existing minimizer-iter crate
- Processes 3.2 Gbp human genome in **4.1 seconds** (forward)
- Processes 3.2 Gbp human genome in **6.0 seconds** (canonical)

**Technical Approach:**
```rust
// Pseudo-code of their approach
pub struct SimdMinimizer<S: Simd> {
    // Generic over SIMD instruction set (AVX2/NEON)
    _phantom: PhantomData<S>,
}

impl SimdMinimizer {
    // Core innovation: SIMD-accelerated minimizer scanning
    fn compute_minimizers_simd(seq: &[u8], k: usize, w: usize) -> Vec<Minimizer> {
        // Use SIMD to:
        // 1. Compute rolling hashes in parallel
        // 2. Find minimum in sliding window with SIMD min operations
        // 3. Track positions with vectorized comparisons
    }
}
```

**Critical Insight:**
- They achieved SIMD speedup on **data-structure-bound** operation (minimizers)
- We said minimizers are 1.02-1.26× in Entry 034
- **They got 4-10× - we need to understand why**

**Hypothesis on Discrepancy:**
1. **Different algorithm:** They may use rolling hash vs our FNV-1a hash
2. **Vectorized minimum finding:** SIMD can efficiently find min in window
3. **Batching strategy:** Process multiple windows simultaneously
4. **Cache optimization:** Better memory access patterns

**Action Item:** Deep dive into their GitHub repo to understand technique.

---

### Parasail (2016, Still State-of-Art) - **Gold Standard**

**What They Built:**
- Comprehensive alignment library (C)
- SSE2, SSE4.1, AVX2, AltiVec, NEON support
- All three alignment types (local, global, semi-global)

**Performance:**
- **136 GCUPS** (Giga Cell Updates Per Second)
- Dual Xeon E5-2670 (24 cores)
- **Highest reported** for striped approach

**Key Innovation: Farrar's Striped Smith-Waterman (2007)**

Traditional approach:
```
Query:    A T G C A T
          ↓ ↓ ↓ ↓ ↓ ↓
Database: A T C C G T

DP matrix filled column by column
Limited parallelism
```

Striped approach:
```
Stripe across diagonal:
  A T G C A T
A [0 1 2 3 4 5]  ← Process these 6 cells in parallel
T [1 . . . . .]
C [2 . . . . .]
G [3 . . . . .]

SIMD register holds 8-16 cells simultaneously
Massive parallelism (6-8× speedup reported)
```

**Technical Details:**

```c
// Pseudo-code of striped SW
void striped_smith_waterman(char* query, char* db) {
    // Query sequence is loaded into SIMD vectors
    __m128i query_vec[query_len / 16];

    // Process database in striped fashion
    for (int i = 0; i < db_len; i++) {
        __m128i db_char = _mm_set1_epi8(db[i]);

        // Compare query vector against single db char
        // All 16 comparisons happen in parallel
        __m128i match = _mm_cmpeq_epi8(query_vec[0], db_char);

        // Update DP scores with SIMD max/add operations
        h_vec = _mm_max_epi8(
            _mm_add_epi8(h_vec, match_score),
            zero_vec
        );
    }
}
```

**Why This Works:**
- Alignment is compute-intensive (not data-structure-bound)
- Embarrassingly parallel across query positions
- Perfect fit for SIMD (element-wise max/add operations)

**Performance Characteristics:**

| Algorithm | Scalar | SIMD | Speedup | Source |
|-----------|--------|------|---------|--------|
| Striped SW | 22 GCUPS | 136 GCUPS | **6.2×** | Parasail paper |
| Our base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7×** | Entry 020 |
| Our GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3×** | Entry 021 |

**Parasail vs biometal NEON:**
- Parasail: 6× speedup (alignment-specific)
- biometal: 16-25× speedup (element-wise ops)
- Different operation types, both excellent

---

### SSW Library (2013) - **Extended Striped SW**

**Innovation Over Farrar:**
- Returns **alignment information** (not just score)
- Traceback implementation
- Production-ready library

**Adoption:**
- MOSAIK (read mapper)
- SCISSORS (split-read mapper)
- TANGRAM (MEI detector)
- RZMBLR (read-overlap graphs)

**Critical for biometal:** If we implement alignment, we need:
1. Score (like Parasail)
2. Alignment info (like SSW)
3. Both are essential for production use

---

### Block Aligner (2022) - **Adaptive SIMD**

**Innovation:**
- Adaptive block sizes based on sequence similarity
- x86 SSE2/AVX2, ARM NEON, WebAssembly SIMD
- Flexible scoring (PAM, BLOSUM, custom)

**Approach:**
```
Instead of fixed stripe width:
┌─────────────────┐
│ Adaptive blocks │  ← Size varies based on alignment difficulty
│  [████]          │  ← Small block for easy region
│  [██████████]    │  ← Large block for hard region
└─────────────────┘
```

**Relevance:** If we implement alignment, adaptive approach may be superior to fixed striping.

---

## Gap Analysis: Why We're Behind

### 1. **Minimizers: SimdMinimizers Got 9.5×, We Got 1.26×**

**Our Entry 034 Finding:**
```
Minimizers NEON: 1.02-1.26× speedup (below ≥5× threshold)
Conclusion: Scalar-only optimal
```

**Their Finding (Jan 2025):**
```
Minimizers SIMD: 9.5× speedup for w=5, 4.5× for w=19
Human genome (3.2 Gbp): 4.1 seconds
```

**Critical Questions:**
1. **Different algorithm?** Rolling hash vs FNV-1a?
2. **Better vectorization?** Finding minimum in window with SIMD?
3. **Data layout?** Memory access patterns optimized?
4. **Implementation quality?** Did we miss optimization opportunity?

**Action Required:**
- Clone `rust-seq/simd-minimizers`
- Study their implementation
- Re-evaluate our Entry 034 conclusion
- Consider: Were we right about data-structure-bound, but wrong about *which* data structure matters?

**Hypothesis:**
- Our Entry 034 tested **FNV-1a hash + HashMap**
- SimdMinimizers may use **rolling hash + SIMD minimum** (no HashMap during scan)
- **The minimizer *selection* can be SIMD-friendly** even if final storage is HashMap

**Revised Understanding:**
```rust
// Our approach (Entry 034):
for window in sequence.windows(w) {
    let mut min_hash = u64::MAX;
    for kmer in window.kmers(k) {
        let hash = fnv1a_hash(kmer);  // ← Scalar hash
        if hash < min_hash {           // ← Scalar comparison
            min_hash = hash;
        }
    }
    minimizers.insert(min_hash);      // ← HashMap (data-structure-bound)
}

// SimdMinimizers approach (hypothesis):
let hashes = compute_rolling_hashes_simd(sequence);  // ← SIMD hashing
for window in hashes.windows(w) {
    let min_idx = simd_argmin(window);  // ← SIMD find minimum index
    minimizers.push(min_idx);           // ← Vec::push (not HashMap)
}
// ↑ This separates SIMD-friendly part (finding min) from data structure
```

**Key Insight:** They may have **decoupled** the SIMD-friendly operations (rolling hash, argmin) from the data-structure-bound operations (deduplication).

---

### 2. **Alignment: We Have Nothing**

**Gap:**
- No Smith-Waterman
- No Needleman-Wunsch
- No semi-global

**Why This Matters:**
- Alignment is **foundational** in bioinformatics
- Most workflows require it (read mapping, variant calling, assembly)
- SIMD provides **6-10× speedup** (proven)
- **High ROI** for implementation

**Complexity:**
- Smith-Waterman: Medium (well-understood, striped approach proven)
- Implementation: ~2,000-3,000 LOC for production-ready version
- Testing: Requires extensive validation (many edge cases)
- Evidence: Would need Entry 035+ experiments

**Technical Debt:**
- No alignment means we can't be comprehensive library
- Limits adoption severely
- rust-bio has it, we don't

---

### 3. **Cross-Platform SIMD: We're ARM-Only**

**Current State:**
```rust
#[cfg(target_arch = "aarch64")]
pub unsafe fn operation_neon(input: &[u8]) -> Result { /* NEON */ }

#[cfg(not(target_arch = "aarch64"))]
pub fn operation_scalar(input: &[u8]) -> Result { /* Scalar */ }
```

**Gap:** No AVX2 implementation for x86_64

**SimdMinimizers approach:**
```rust
// Generic over SIMD instruction set
pub trait Simd {
    fn compare(&self, a: Self::Vector, b: Self::Vector) -> Self::Mask;
    fn min(&self, a: Self::Vector, b: Self::Vector) -> Self::Vector;
    // ...
}

impl Simd for Avx2 { /* AVX2 intrinsics */ }
impl Simd for Neon { /* NEON intrinsics */ }
```

**Benefit:**
- x86_64 users get speedup too (not just ARM)
- Broader adoption
- Fair comparison (both platforms optimized)

**Challenge:**
- Increased complexity
- More testing burden
- Maintenance overhead

**Strategic Decision:**
- **Stay ARM-focused?** (democratization mission)
- **Add AVX2?** (broader adoption, fairer comparisons)

---

## Strategic Options for biometal

### Option 1: **Double Down on Current Strengths** (Low Risk)

**Focus:** Polish what we have
- Perfect the element-wise operations (already 16-25×)
- Add comprehensive benchmarks vs BioPython/rust-bio
- Publish evidence base
- Educational content expansion

**Pros:**
- ✅ Low risk (we're already excellent here)
- ✅ Clear differentiation (best element-wise performance)
- ✅ Publishable now
- ✅ Aligns with evidence-based mission

**Cons:**
- ❌ Narrow scope (preprocessing only)
- ❌ Doesn't address alignment gap
- ❌ Falls behind SimdMinimizers on minimizers
- ❌ Limited adoption potential

**Timeline:** 1-2 months
**Outcome:** Solid niche tool, academic publication

---

### Option 2: **Fix Minimizers with SimdMinimizers Approach** (Medium Risk)

**Focus:** Match SimdMinimizers performance
- Study their implementation
- Implement rolling hash + SIMD argmin
- Re-benchmark (target: 4-5× speedup)
- Update Entry 034 with revised findings

**Pros:**
- ✅ Competitive with cutting-edge (Jan 2025 paper)
- ✅ Improves k-mer operations significantly
- ✅ Still within our data-structure expertise
- ✅ Evidence-based approach (re-evaluate Entry 034)

**Cons:**
- ⚠️ Admits we missed optimization in Entry 034
- ⚠️ Medium complexity (new algorithm)
- ❌ Still doesn't solve alignment gap

**Timeline:** 1-2 months (research + implementation)
**Outcome:** Best-in-class minimizers, matches SimdMinimizers

**Evidence Requirements:**
- Entry 035: Re-evaluate minimizers with rolling hash
- Compare: FNV-1a vs rolling hash vs SIMD minimizer selection
- Target: 4-5× speedup (match SimdMinimizers)

---

### Option 3: **Implement SIMD Alignment** (High Risk, High Reward)

**Focus:** Add striped Smith-Waterman
- Implement Farrar's striped approach
- NEON-optimized (follow ARM focus)
- Local, global, semi-global variants
- Comprehensive testing

**Pros:**
- ✅ **Closes major gap** (alignment is foundational)
- ✅ High-impact feature (enables many workflows)
- ✅ Proven technique (Farrar 2007, Parasail validation)
- ✅ Aligns with SIMD expertise
- ✅ 6-10× speedup proven achievable

**Cons:**
- ❌ High complexity (~2,000-3,000 LOC)
- ❌ Extensive testing required
- ❌ 2-3 months implementation time
- ❌ Requires Entry 035+ experiments for evidence base

**Timeline:** 2-3 months (research + implementation + testing)
**Outcome:** Comprehensive library, major capability boost

**Technical Approach:**

```rust
// High-level design
pub mod alignment {
    pub struct StripedSmithWaterman {
        query: Vec<u8>,
        scoring: ScoringMatrix,
    }

    impl StripedSmithWaterman {
        #[cfg(target_arch = "aarch64")]
        pub unsafe fn align_neon(&self, target: &[u8]) -> Alignment {
            // Farrar's striped approach with ARM NEON
            // Process 16 query positions in parallel
            // Expected: 6-10× speedup
        }

        #[cfg(not(target_arch = "aarch64"))]
        pub fn align_scalar(&self, target: &[u8]) -> Alignment {
            // Fallback implementation
        }
    }

    pub struct Alignment {
        pub score: i32,
        pub query_start: usize,
        pub query_end: usize,
        pub target_start: usize,
        pub target_end: usize,
        pub cigar: Option<String>,  // Optional: SSW-style traceback
    }
}
```

**Evidence Requirements:**
- Entry 035: Striped Smith-Waterman NEON implementation
- Benchmark: Query lengths 50-500bp, various target sizes
- Compare: Scalar vs NEON, target 6-10× speedup
- Validate: Correctness with known alignments

---

### Option 4: **Hybrid: Fix Minimizers + Add Alignment** (Highest Risk/Reward)

**Focus:** Address both gaps simultaneously
- Phase 1: Fix minimizers (1 month)
- Phase 2: Add alignment (2 months)
- Comprehensive SIMD library

**Pros:**
- ✅ **Comprehensive solution** (matches rust-bio scope)
- ✅ Competitive with all SIMD tools
- ✅ Strong differentiation (evidence-based + SIMD)
- ✅ Maximum impact potential

**Cons:**
- ❌ **Very high risk** (3+ months, complex)
- ❌ Extensive testing burden
- ❌ May delay publication
- ❌ Resource-intensive (just you?)

**Timeline:** 3-4 months
**Outcome:** Comprehensive bioinformatics library with SIMD

---

## Recommendation: Phased Approach

### Phase 1 (Immediate: 0-1 month)
**Publish What We Have**

1. **Academic paper** - Evidence-based optimization methodology
   - Target: Genome Biology or Bioinformatics
   - Focus: 1,357 experiments, statistical rigor
   - Case study: Element-wise operations (16-25× proven)

2. **Release v1.2.0** as production-ready
   - Current capabilities well-documented
   - Limitations acknowledged
   - Clear roadmap

**Rationale:** Secure academic credit NOW before expanding scope

---

### Phase 2 (Next: 1-2 months)
**Fix Minimizers (Match SimdMinimizers)**

1. **Deep dive** into SimdMinimizers implementation
   - Clone repo: `rust-seq/simd-minimizers`
   - Understand rolling hash + SIMD argmin approach
   - Identify why they got 9.5× vs our 1.26×

2. **Implement** improved minimizers
   - Rolling hash for SIMD-friendly computation
   - SIMD minimum finding in windows
   - Maintain scalar fallback

3. **Evidence** - Entry 035
   - Benchmark: Rolling hash vs FNV-1a
   - Target: 4-5× speedup (match their performance)
   - Validate: Human genome in ~5 seconds

4. **Release v1.3.0** - Improved k-mer operations

**Rationale:**
- Medium complexity, high value
- Keeps us competitive with cutting-edge
- Stays within our data-structure expertise

---

### Phase 3 (Future: 2-3 months)
**Add SIMD Alignment (If Community Demands)**

1. **Research** Farrar's striped Smith-Waterman
   - Study Parasail implementation
   - Understand striping technique thoroughly
   - Design ARM NEON-specific optimizations

2. **Implement** core alignment
   - Smith-Waterman (local) first
   - NEON-optimized striped approach
   - Target: 6-10× speedup

3. **Expand** alignment types
   - Needleman-Wunsch (global)
   - Semi-global variants
   - Traceback for CIGAR strings

4. **Evidence** - Entry 036-038
   - Benchmark: Various query/target sizes
   - Validate: Correctness tests
   - Compare: Parasail performance

5. **Release v2.0.0** - Comprehensive SIMD library

**Rationale:**
- Only if Phase 1-2 generate community interest
- High complexity, but proven technique
- Transforms biometal from niche to comprehensive

---

## Technical Deep Dives

### A. How SimdMinimizers Achieves 9.5× Speedup

Based on their paper and our Entry 034:

**Key Innovation: Separate SIMD-Friendly Operations**

```rust
// Traditional approach (our Entry 034):
// Everything coupled together
for window in sequence.windows(w) {
    let mut min_hash = u64::MAX;
    let mut min_pos = 0;

    for (i, kmer) in window.kmers(k).enumerate() {
        let hash = fnv1a(kmer);  // ← Scalar (not vectorizable)
        if hash < min_hash {      // ← Scalar comparison
            min_hash = hash;
            min_pos = i;
        }
    }

    minimizers.insert(min_pos, min_hash);  // ← HashMap (data-structure-bound)
}

// SimdMinimizers approach:
// Step 1: Compute ALL hashes with SIMD
let all_hashes: Vec<u64> = compute_rolling_hashes_simd(sequence);
// ↑ This is SIMD-friendly (parallel hash computation)

// Step 2: Find minimums in windows with SIMD
for window in all_hashes.windows(w) {
    let (min_idx, min_hash) = simd_find_minimum(window);
    // ↑ SIMD min reduction (hardware-optimized)

    minimizers.push((min_idx, min_hash));
    // ↑ Simple Vec::push (no HashMap during scan)
}

// Step 3: Deduplicate (data-structure-bound part is isolated)
minimizers.dedup();  // ← This part is still slow, but happens once
```

**Why This Works:**

1. **Rolling hash is SIMD-friendly:**
   ```
   hash[i+1] = (hash[i] - seq[i] * POW) * BASE + seq[i+k]
   ↑ All arithmetic, can vectorize 16 positions at once
   ```

2. **Finding minimum is SIMD-friendly:**
   ```neon
   vminq_u64(vec1, vec2)  // ← Single NEON instruction
   ```

3. **Data structure operations isolated:**
   - Deduplication happens once at end
   - Not in critical path
   - 10% of total time vs 50% in our approach

**Our Mistake in Entry 034:**
- We tested FNV-1a hash (not vectorizable)
- We coupled hash + HashMap operations (inseparable)
- We didn't try rolling hash + batch processing

**Lesson:** Sometimes the *algorithm choice* determines SIMD-friendliness, not just the operation type.

---

### B. Striped Smith-Waterman Implementation Guide

**Core Concept:**

Traditional DP:
```
Query:     A T G C
Database:  A T C G
           ↓
Fill column-by-column (sequential)
```

Striped DP:
```
Process diagonal stripe:
Query positions:  0 1 2 3
                 ┌─┬─┬─┬─┐
                 │A│T│G│C│  ← Load into SIMD register
                 └─┴─┴─┴─┘
                  ▼ ▼ ▼ ▼
Database char: A → Compare all 4 positions in parallel
              T → Compare all 4 positions in parallel
              C → Compare all 4 positions in parallel
```

**NEON Implementation Pattern:**

```rust
pub unsafe fn striped_smith_waterman_neon(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringMatrix,
) -> i32 {
    let stripe_width = 16; // NEON processes 16 bytes at once
    let n_stripes = (query.len() + stripe_width - 1) / stripe_width;

    // Initialize DP vectors (horizontal and vertical)
    let mut h_vecs: Vec<uint8x16_t> = vec![vdupq_n_u8(0); n_stripes];
    let mut e_vecs: Vec<uint8x16_t> = vec![vdupq_n_u8(0); n_stripes];
    let mut f_vecs: Vec<uint8x16_t> = vec![vdupq_n_u8(0); n_stripes];

    let mut max_score = vdupq_n_u8(0);

    // Process target sequence
    for &target_char in target {
        let target_vec = vdupq_n_u8(target_char);

        for stripe_idx in 0..n_stripes {
            // Load query stripe
            let query_stripe = load_query_stripe(query, stripe_idx, stripe_width);

            // Match/mismatch scores (SIMD comparison)
            let match_mask = vceqq_u8(query_stripe, target_vec);
            let match_score = vbslq_u8(
                match_mask,
                vdupq_n_u8(scoring.match_score as u8),
                vdupq_n_u8(scoring.mismatch_penalty as u8)
            );

            // Compute H (diagonal + match)
            let h_diag = if stripe_idx > 0 {
                shift_right(h_vecs[stripe_idx - 1], 1)
            } else {
                vdupq_n_u8(0)
            };

            h_vecs[stripe_idx] = vqaddq_u8(h_diag, match_score);

            // Compute E (insertion)
            e_vecs[stripe_idx] = vmaxq_u8(
                vqsubq_u8(h_vecs[stripe_idx], vdupq_n_u8(scoring.gap_open as u8)),
                vqsubq_u8(e_vecs[stripe_idx], vdupq_n_u8(scoring.gap_extend as u8))
            );

            // Compute F (deletion)
            f_vecs[stripe_idx] = vmaxq_u8(
                vqsubq_u8(h_vecs[stripe_idx], vdupq_n_u8(scoring.gap_open as u8)),
                vqsubq_u8(f_vecs[stripe_idx], vdupq_n_u8(scoring.gap_extend as u8))
            );

            // Take maximum (Smith-Waterman max)
            h_vecs[stripe_idx] = vmaxq_u8(
                vmaxq_u8(h_vecs[stripe_idx], e_vecs[stripe_idx]),
                vmaxq_u8(f_vecs[stripe_idx], vdupq_n_u8(0))
            );

            // Track maximum score across all cells
            max_score = vmaxq_u8(max_score, h_vecs[stripe_idx]);
        }
    }

    // Extract maximum score from SIMD register
    extract_max_u8(max_score) as i32
}

unsafe fn extract_max_u8(vec: uint8x16_t) -> u8 {
    let mut max = 0u8;
    for i in 0..16 {
        let val = vgetq_lane_u8(vec, i);
        if val > max {
            max = val;
        }
    }
    max
}
```

**Expected Performance:**
- Scalar: ~100-200 MCUPS (Million Cell Updates Per Second)
- NEON: ~600-1,200 MCUPS (6-10× speedup)
- Parasail baseline: 136 GCUPS on dual Xeon (our target is lower but proportional)

**Testing Strategy:**
1. **Correctness:** Known alignment benchmarks (standard test suite)
2. **Performance:** Query lengths 50-500bp, various target sizes
3. **Edge cases:** Empty sequences, identical sequences, no match
4. **Scoring matrices:** PAM, BLOSUM, custom

---

### C. Cross-Platform SIMD: Generic SIMD Trait

If we want to support both NEON and AVX2:

```rust
pub trait SimdOps {
    type Vector;
    type Mask;

    // Core operations
    fn load(data: &[u8]) -> Self::Vector;
    fn store(vec: Self::Vector, dest: &mut [u8]);
    fn compare_eq(a: Self::Vector, b: Self::Vector) -> Self::Mask;
    fn add(a: Self::Vector, b: Self::Vector) -> Self::Vector;
    fn max(a: Self::Vector, b: Self::Vector) -> Self::Vector;
    fn blend(mask: Self::Mask, a: Self::Vector, b: Self::Vector) -> Self::Vector;
}

#[cfg(target_arch = "aarch64")]
pub struct NeonOps;

#[cfg(target_arch = "aarch64")]
impl SimdOps for NeonOps {
    type Vector = uint8x16_t;
    type Mask = uint8x16_t;

    fn load(data: &[u8]) -> Self::Vector {
        unsafe { vld1q_u8(data.as_ptr()) }
    }

    fn compare_eq(a: Self::Vector, b: Self::Vector) -> Self::Mask {
        unsafe { vceqq_u8(a, b) }
    }

    // ... implement all operations
}

#[cfg(target_arch = "x86_64")]
pub struct Avx2Ops;

#[cfg(target_arch = "x86_64")]
impl SimdOps for Avx2Ops {
    type Vector = __m256i;
    type Mask = __m256i;

    fn load(data: &[u8]) -> Self::Vector {
        unsafe { _mm256_loadu_si256(data.as_ptr() as *const __m256i) }
    }

    fn compare_eq(a: Self::Vector, b: Self::Vector) -> Self::Mask {
        unsafe { _mm256_cmpeq_epi8(a, b) }
    }

    // ... implement all operations
}

// Then operations are generic:
pub fn count_bases_simd<S: SimdOps>(seq: &[u8]) -> [u32; 4] {
    // Works with both NEON and AVX2
    let a_vec = S::load(...);
    let comparison = S::compare_eq(a_vec, ...);
    // ...
}
```

**Pros:**
- Single codebase for both platforms
- Fair performance comparisons
- Broader adoption

**Cons:**
- Additional complexity
- More testing surface
- May compromise per-platform optimization

**Decision:** Consider for v2.0, focus on ARM for v1.x (mission alignment).

---

## Concrete Next Steps (Action Plan)

### Week 1-2: Research Phase

**Tasks:**
1. Clone and study `rust-seq/simd-minimizers`
   - Read their paper thoroughly
   - Understand rolling hash implementation
   - Analyze SIMD argmin technique

2. Benchmark their library
   - Install and run on your M-series Mac
   - Compare performance to our Entry 034 results
   - Understand real-world performance

3. Deep dive Parasail
   - Clone and build Parasail
   - Study striped SW implementation
   - Understand NEON specifics (if they have ARM version)

**Deliverable:** `docs/SIMD_RESEARCH_FINDINGS.md` with:
- SimdMinimizers technique analysis
- Parasail implementation notes
- Revised understanding of Entry 034

---

### Week 3-4: Decision Point

**Decision Matrix:**

| Option | Complexity | Timeline | Impact | Risk |
|--------|------------|----------|--------|------|
| Publish as-is | Low | 2 weeks | Medium | Low |
| Fix minimizers | Medium | 1-2 months | Medium-High | Medium |
| Add alignment | High | 2-3 months | High | High |
| Both | Very High | 3-4 months | Very High | Very High |

**Recommendation:**
1. **Publish evidence base** (academic paper) - 2 weeks
2. **Fix minimizers** (v1.3.0) - 1-2 months
3. **Decide on alignment** based on community feedback

---

### Month 2-3: Implementation (If Minimizers)

**Phase 1: Rolling Hash Implementation**
```rust
// New module: src/operations/rolling_hash.rs
pub struct RollingHash {
    // Implement rolling hash with SIMD support
}

// Update: src/operations/kmer.rs
pub fn extract_minimizers_fast(seq: &[u8], k: usize, w: usize) -> Vec<Minimizer> {
    // Use rolling hash + SIMD argmin
    // Target: 4-5× speedup over current implementation
}
```

**Phase 2: Benchmarking**
- Entry 035: Rolling hash vs FNV-1a
- Target: Match SimdMinimizers (4-5×)
- Validate: Human genome in ~5 seconds

**Phase 3: Integration**
- Update Python bindings
- Add benchmarks
- Documentation
- Release v1.3.0

---

## Success Metrics

### Short-term (3 months)
- ✅ Academic paper submitted (evidence base)
- ✅ Minimizers performance: 4-5× speedup (match SimdMinimizers)
- ✅ v1.3.0 released with improved k-mer operations

### Medium-term (6 months)
- ✅ Paper accepted and published
- ✅ Community adoption (100+ GitHub stars)
- ✅ Decision made on alignment implementation
- ✅ v2.0.0 roadmap finalized

### Long-term (12 months)
- ✅ Alignment implementation (if pursued): 6-10× speedup
- ✅ Comprehensive SIMD library (element-wise + k-mer + alignment)
- ✅ Competitive with rust-bio in scope
- ✅ Best-in-class ARM performance

---

## Conclusion

**Where We Stand:**
- ✅ **Excellent** at element-wise operations (16-25×)
- ⚠️ **Behind** on minimizers (1.26× vs SimdMinimizers' 9.5×)
- ❌ **Missing** alignment entirely

**Strategic Path Forward:**
1. **Publish evidence base NOW** (secure academic credit)
2. **Fix minimizers** (match cutting-edge performance)
3. **Add alignment if community demands** (high-risk, high-reward)

**Core Insight:**
SimdMinimizers showed us that **algorithm design** determines SIMD-friendliness more than operation type. We can achieve better k-mer performance by:
- Using rolling hash (vectorizable)
- Separating SIMD-friendly ops (argmin) from data-structure-bound ops (deduplication)
- Batch processing

**Final Recommendation:**
Take phased approach - secure what we have (publish), then expand capabilities based on community feedback. Don't over-commit before validating demand.

---

**Document Status:** Strategic planning complete
**Next Action:** Weeks 1-2 research phase (SimdMinimizers deep dive)
**Decision Point:** Week 4 (publish vs expand vs both)
