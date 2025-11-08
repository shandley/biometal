# Phase 4: Sequence Manipulation Primitives

**Status**: Planning (Not Started)
**Version**: v1.1.0 (post v1.0.0)
**Date**: November 5, 2025

---

## Mission

Provide foundational sequence manipulation primitives that CLI tools can compose to implement features like:
- Adapter trimming
- Quality-based read filtering
- Read orientation correction
- Subsequence extraction
- Target region isolation

**Philosophy**: Library provides primitives, CLI tools compose them.

---

## Architecture Principle

**Bad** (monolithic):
```rust
// ❌ Library shouldn't provide full workflow
pub fn trim_adapters_and_filter(input: &Path, output: &Path, params: TrimParams)
```

**Good** (composable primitives):
```rust
// ✅ Library provides building blocks
let rc = reverse_complement(&seq);
let trimmed = trim_quality(&record, 20);
let extracted = extract_region(&record, 10, 50);

// CLI tool composes these primitives into workflows
```

---

## Primitive Categories

### 1. Core Sequence Operations
Basic DNA/RNA sequence manipulations

### 2. Record Transformations
Operations that modify FastqRecord structures

### 3. Quality-Based Operations
Operations driven by quality scores

### 4. Validation & Statistics
Helper functions for correctness checking

---

## Complete Primitive Inventory

### Category 1: Core Sequence Operations

#### 1.1 Reverse Complement
```rust
/// Reverse complement a DNA sequence
///
/// # Examples
/// ```
/// use biometal::operations::reverse_complement;
///
/// let seq = b"ATGC";
/// let rc = reverse_complement(seq);
/// assert_eq!(rc, b"GCAT");
/// ```
///
/// # Evidence
/// - Scalar: No evidence needed (standard algorithm)
/// - NEON: Category 2 (Similar to base_counting, potential 10-16× speedup)
///
/// # Performance
/// - Scalar: O(n) single pass
/// - NEON: O(n/16) with 16-byte SIMD operations
pub fn reverse_complement(seq: &[u8]) -> Vec<u8>

/// In-place reverse complement (avoids allocation)
pub fn reverse_complement_inplace(seq: &mut [u8])
```

**Priority**: HIGH
**Use Cases**:
- Read orientation correction
- Adapter matching (must check both orientations)
- Paired-end read processing
- Reference genome queries

**Evidence Category**:
- **Scalar**: No validation needed (standard algorithm)
- **NEON**: Category 2 (Similar) - Quick criterion benchmark (N=10)
- **Expected NEON Speedup**: 10-16× (similar to base_counting)

**Implementation Notes**:
- Complement table: A↔T, G↔C, handle ambiguous bases (N, R, Y, etc.)
- NEON: Can process 16 bases at a time with lookup tables
- Edge cases: Non-DNA characters, empty sequences

---

#### 1.2 Reverse Only
```rust
/// Reverse a sequence without complementing
///
/// # Examples
/// ```
/// let seq = b"ATGC";
/// let rev = reverse(seq);
/// assert_eq!(rev, b"CGTA");
/// ```
pub fn reverse(seq: &[u8]) -> Vec<u8>

/// In-place reverse
pub fn reverse_inplace(seq: &mut [u8])
```

**Priority**: MEDIUM
**Use Cases**: Rarely needed alone, but useful for testing

---

#### 1.3 Complement Only
```rust
/// Complement a sequence without reversing
///
/// # Examples
/// ```
/// let seq = b"ATGC";
/// let comp = complement(seq);
/// assert_eq!(comp, b"TACG");
/// ```
pub fn complement(seq: &[u8]) -> Vec<u8>
pub fn complement_inplace(seq: &mut [u8])
```

**Priority**: LOW
**Use Cases**: Rare, but completes the API

---

#### 1.4 Validate Sequence
```rust
/// Check if sequence contains only valid DNA bases
///
/// # Examples
/// ```
/// assert!(is_valid_dna(b"ATGCN"));
/// assert!(!is_valid_dna(b"ATGCX"));
/// ```
pub fn is_valid_dna(seq: &[u8]) -> bool
pub fn is_valid_rna(seq: &[u8]) -> bool

/// Count invalid characters
pub fn count_invalid_bases(seq: &[u8]) -> usize
```

**Priority**: MEDIUM
**Use Cases**: Input validation, quality control

---

### Category 2: Subsequence Operations

#### 2.1 Extract Subsequence
```rust
/// Extract subsequence by position (zero-copy when possible)
///
/// # Examples
/// ```
/// let seq = b"ATGCATGC";
/// let sub = subsequence(seq, 2, 6);
/// assert_eq!(sub, b"GCAT");
/// ```
///
/// # Evidence
/// Category 1 (Built-in): Uses Rust slice operations
pub fn subsequence(seq: &[u8], start: usize, end: usize) -> &[u8]
```

**Priority**: HIGH
**Evidence**: Category 1 (built-in Rust slicing, no validation needed)

---

#### 2.2 Extract with Quality
```rust
/// Extract region from FastqRecord (sequence + quality)
///
/// # Examples
/// ```
/// let record = FastqRecord::new("read1", b"ATGCATGC", b"IIIIIIII");
/// let extracted = extract_region(&record, 2, 6);
/// assert_eq!(extracted.sequence, b"GCAT");
/// assert_eq!(extracted.quality, b"IIII");
/// ```
pub fn extract_region(record: &FastqRecord, start: usize, end: usize) -> Result<FastqRecord>
```

**Priority**: HIGH
**Use Cases**:
- Target region extraction
- Barcode/UMI extraction
- Read trimming

**Validation**:
- Check bounds
- Ensure sequence/quality length match
- Handle edge cases (start=end, etc.)

---

### Category 3: Trimming Operations

#### 3.1 Fixed-Position Trimming
```rust
/// Trim N bases from 5' end (start of read)
///
/// # Examples
/// ```
/// let seq = b"ATGCATGC";
/// let trimmed = trim_start(seq, 2);
/// assert_eq!(trimmed, b"GCATGC");
/// ```
pub fn trim_start(seq: &[u8], bases: usize) -> &[u8]

/// Trim N bases from 3' end (end of read)
pub fn trim_end(seq: &[u8], bases: usize) -> &[u8]

/// Trim from both ends
pub fn trim_both(seq: &[u8], start_bases: usize, end_bases: usize) -> &[u8]
```

**Priority**: HIGH
**Evidence**: Category 1 (built-in slicing)

---

#### 3.2 Quality-Based Trimming
```rust
/// Trim from 3' end until minimum quality threshold met
///
/// # Algorithm
/// Scans from 3' end backward, stops when quality ≥ threshold
///
/// # Examples
/// ```
/// let record = FastqRecord::new("read1", b"ATGCATGC", b"IIII!!!!"); // Q40, Q0
/// let trimmed = trim_quality_end(&record, 20); // Min Q=20
/// assert_eq!(trimmed.sequence, b"ATGC"); // Trimmed low-quality tail
/// ```
///
/// # Evidence
/// Scalar only (no NEON benefit for sequential quality scanning)
pub fn trim_quality_end(record: &FastqRecord, min_quality: u8) -> Result<FastqRecord>

/// Trim from 5' end until quality threshold met
pub fn trim_quality_start(record: &FastqRecord, min_quality: u8) -> Result<FastqRecord>

/// Trim from both ends
pub fn trim_quality_both(record: &FastqRecord, min_quality: u8) -> Result<FastqRecord>
```

**Priority**: HIGH
**Evidence**: Scalar only (no NEON optimization needed)

**Algorithm**:
```
For 3' trimming:
  1. Start at end of quality array
  2. Scan backward
  3. Stop when quality[i] >= threshold
  4. Return record[0..i+1]
```

**Use Cases**:
- Standard preprocessing step
- Illumina reads often have quality drop at 3' end
- Improves downstream analysis (assembly, alignment)

---

#### 3.3 Sliding Window Quality Trimming
```rust
/// Trim using sliding window approach
///
/// # Algorithm
/// - Use window of size N (e.g., 4 bases)
/// - Calculate mean quality of window
/// - Trim when window quality < threshold
///
/// # Examples
/// ```
/// let record = FastqRecord::new("read1", b"ATGCATGC", b"IIII!!!!");
/// let trimmed = trim_quality_window(&record, 20, 4); // threshold=20, window=4
/// ```
///
/// # Evidence
/// Scalar (no NEON benefit for mean calculation in small windows)
pub fn trim_quality_window(
    record: &FastqRecord,
    min_quality: u8,
    window_size: usize
) -> Result<FastqRecord>
```

**Priority**: MEDIUM
**Use Cases**: More robust than single-base trimming, handles noise better

---

### Category 4: Masking Operations

#### 4.1 Quality-Based Masking
```rust
/// Mask (replace with 'N') bases below quality threshold
///
/// # Examples
/// ```
/// let mut record = FastqRecord::new("read1", b"ATGC", b"I!I!");
/// mask_low_quality(&mut record, 20);
/// assert_eq!(record.sequence, b"ANAC"); // Masked Q<20
/// ```
pub fn mask_low_quality(record: &mut FastqRecord, min_quality: u8)

/// Return masked copy (non-mutating)
pub fn mask_low_quality_copy(record: &FastqRecord, min_quality: u8) -> FastqRecord
```

**Priority**: MEDIUM
**Use Cases**:
- Keep read length but mark unreliable bases
- Some aligners handle 'N' specially
- Quality-aware variant calling

---

### Category 5: Orientation Operations

#### 5.1 Record-Level Reverse Complement
```rust
/// Reverse complement an entire FastqRecord (sequence + quality)
///
/// # Examples
/// ```
/// let record = FastqRecord::new("read1", b"ATGC", b"IIII");
/// let rc_record = reverse_complement_record(&record);
/// assert_eq!(rc_record.sequence, b"GCAT");
/// assert_eq!(rc_record.quality, b"IIII"); // Quality reversed too
/// ```
///
/// # Evidence
/// Scalar: Composition of reverse_complement + quality reverse
/// NEON: Category 2 (inherits from reverse_complement speedup)
pub fn reverse_complement_record(record: &FastqRecord) -> FastqRecord
pub fn reverse_complement_record_inplace(record: &mut FastqRecord)
```

**Priority**: HIGH
**Use Cases**:
- Paired-end read processing (R2 is reverse complement of insert)
- Strand-specific protocols
- Read orientation correction

---

### Category 6: Length Operations

#### 6.1 Length Filtering Helpers
```rust
/// Check if record meets length criteria
///
/// # Examples
/// ```
/// let record = FastqRecord::new("read1", b"ATGC", b"IIII");
/// assert!(meets_length_requirement(&record, 3, 10));
/// assert!(!meets_length_requirement(&record, 5, 10));
/// ```
pub fn meets_length_requirement(record: &FastqRecord, min: usize, max: usize) -> bool

/// Get sequence length (convenience wrapper)
pub fn sequence_length(record: &FastqRecord) -> usize
```

**Priority**: LOW (trivial, but completes API)

---

## API Design Principles

### 1. Zero-Copy When Possible
```rust
// ✅ Good: Returns slice (zero-copy)
pub fn trim_start(seq: &[u8], bases: usize) -> &[u8]

// ⚠️  Acceptable: Must allocate (reverse complement)
pub fn reverse_complement(seq: &[u8]) -> Vec<u8>

// ✅ Better: Offer in-place variant
pub fn reverse_complement_inplace(seq: &mut [u8])
```

### 2. Error Handling
```rust
// Operations that can fail return Result
pub fn extract_region(record: &FastqRecord, start: usize, end: usize) -> Result<FastqRecord>

// Operations that can't fail return directly
pub fn trim_start(seq: &[u8], bases: usize) -> &[u8] // Clamps to valid range
```

### 3. Consistent Naming
- `trim_*`: Remove bases from ends
- `extract_*`: Get subsequence
- `mask_*`: Replace bases with 'N'
- `*_inplace`: Mutating variant
- `*_record`: Operates on FastqRecord (not just sequence)

### 4. Performance Variants
```rust
// Standard (allocating)
pub fn reverse_complement(seq: &[u8]) -> Vec<u8>

// In-place (no allocation)
pub fn reverse_complement_inplace(seq: &mut [u8])

// NEON-optimized (internal)
#[cfg(target_arch = "aarch64")]
unsafe fn reverse_complement_neon(seq: &[u8]) -> Vec<u8>

// Public API (dispatches to best implementation)
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    #[cfg(target_arch = "aarch64")]
    { unsafe { reverse_complement_neon(seq) } }

    #[cfg(not(target_arch = "aarch64"))]
    { reverse_complement_scalar(seq) }
}
```

---

## Testing Strategy

### 1. Property-Based Testing (proptest)
```rust
proptest! {
    #[test]
    fn test_reverse_complement_involutive(seq in "[ACGT]{1,1000}") {
        // RC(RC(seq)) = seq
        let rc = reverse_complement(seq.as_bytes());
        let rc_rc = reverse_complement(&rc);
        prop_assert_eq!(rc_rc, seq.as_bytes());
    }

    #[test]
    fn test_trim_preserves_validity(
        seq in "[ACGT]{10,100}",
        start in 0usize..10,
        end in 0usize..10
    ) {
        let trimmed = trim_both(seq.as_bytes(), start, end);
        prop_assert!(is_valid_dna(trimmed));
    }
}
```

### 2. Unit Tests
- Edge cases (empty sequences, length=1, all same base)
- Known transformations (specific inputs → expected outputs)
- Complement table correctness (A→T, G→C, etc.)
- Boundary conditions (start=0, end=length, etc.)

### 3. Integration Tests
- CLI-ready patterns (read → trim → filter → write)
- Round-trip correctness
- Large-scale tests (10K+ records)

### 4. Benchmark Tests (criterion, N=10)
- Only for NEON optimization candidates
- Compare scalar vs NEON
- Measure speedup against threshold (e.g., ≥5× for worthwhile)

---

## Evidence Classification

### Category 1 (Mirror/Built-in) - No Validation Needed
- `subsequence()` - Rust slicing
- `trim_*()` family - Rust slicing
- `extract_region()` - Composition of validated operations

### Category 2 (Similar) - Quick Validation (N=10)
- `reverse_complement()` NEON - Similar to base_counting
- `reverse_complement_record()` NEON - Inherits from above

### Category 3 (Novel) - Full Validation (N=30)
- None expected for Phase 4 (all operations are standard algorithms)

### Scalar Only (No NEON Benefit)
- Quality-based trimming (sequential scanning)
- Quality-based masking (sequential processing)
- Validation functions (early-exit patterns)

---

## Implementation Order

### Week 1: Core Primitives (Days 1-7)
1. `reverse_complement()` - scalar + tests
2. `reverse_complement_inplace()` - in-place variant
3. `complement()` + `reverse()` - building blocks
4. `is_valid_dna()` - validation helper
5. Property tests for reversibility

**Success Criteria**: All core sequence ops work, 100% test coverage

---

### Week 2: Record Operations (Days 8-14)
1. `extract_region()` - sequence + quality extraction
2. `trim_start()`, `trim_end()`, `trim_both()` - fixed trimming
3. `reverse_complement_record()` - full record transformation
4. Integration tests showing CLI patterns

**Success Criteria**: Can compose operations in CLI-ready workflows

---

### Week 3: Quality Operations (Days 15-21)
1. `trim_quality_end()` - 3' quality trimming
2. `trim_quality_start()` - 5' quality trimming
3. `trim_quality_both()` - bidirectional trimming
4. `trim_quality_window()` - sliding window approach
5. `mask_low_quality()` - quality-based masking

**Success Criteria**: Quality-aware operations match Trimmomatic behavior

---

### Week 4: Optimization & Polish (Days 22-28)
1. Benchmark with criterion (N=10)
2. NEON optimization for `reverse_complement()` if justified
3. Documentation review
4. Performance guide for users
5. Update CHANGELOG.md

**Success Criteria**:
- All operations meet performance targets
- Documentation complete
- Ready for v1.1.0 release

---

## File Structure

```
src/
├── lib.rs                    # Re-export operations
└── operations/
    ├── mod.rs                # Module organization
    ├── base_counting.rs      # Existing
    ├── gc_content.rs         # Existing
    ├── quality_filter.rs     # Existing
    ├── sequence.rs           # NEW: Core sequence ops
    │   ├── reverse_complement()
    │   ├── complement()
    │   ├── reverse()
    │   └── is_valid_dna()
    ├── trimming.rs           # NEW: Trimming operations
    │   ├── trim_start()
    │   ├── trim_end()
    │   ├── trim_quality_*()
    │   └── trim_quality_window()
    ├── masking.rs            # NEW: Masking operations
    │   └── mask_low_quality()
    └── record_ops.rs         # NEW: Record-level operations
        ├── extract_region()
        ├── reverse_complement_record()
        └── meets_length_requirement()

tests/
└── sequence_integration.rs   # NEW: Integration tests
```

---

## Success Criteria

### Functional
- ✅ All 20+ primitives implemented
- ✅ Property-based tests pass (1000+ generated cases each)
- ✅ Integration tests show CLI-ready patterns
- ✅ Round-trip correctness verified
- ✅ Edge cases handled (empty, length=1, etc.)

### Performance
- ✅ Scalar implementations correct
- ✅ NEON optimization only if benchmark justifies (≥5× speedup)
- ✅ Memory usage constant (Rule 5)
- ✅ No performance regressions in existing code

### Quality
- ✅ All operations documented with examples
- ✅ Evidence citations for optimizations
- ✅ Error handling complete (no panics)
- ✅ API consistent with biometal patterns

### Deliverable
- ✅ Ready for v1.1.0 release
- ✅ Updated README with new capabilities
- ✅ CHANGELOG.md reflects new primitives
- ✅ Migration guide for users (if needed)

---

## Open Questions for User

### Question 1: RNA Support?
Should we include RNA operations (U instead of T)?

**Options**:
- A: DNA only (simpler, most common)
- B: DNA + RNA (more complete)
- C: Generic nucleotide operations (handles both)

**Recommendation**: Start with DNA (A), add RNA if requested

---

### Question 2: Ambiguous Base Handling?
How to handle IUPAC ambiguity codes (N, R, Y, W, S, K, etc.)?

**Options**:
- A: Strict (only ACGT, error on ambiguous)
- B: Pass-through (preserve ambiguous codes)
- C: Configurable (parameter for behavior)

**Recommendation**: Pass-through (B), document behavior

---

### Question 3: Quality Score Format?
Which quality score format(s) to support?

**Current**: Phred+33 (Illumina 1.8+, most common)

**Options**:
- A: Phred+33 only (simplest)
- B: Auto-detect format (Phred+33 vs Phred+64)
- C: Configurable via parameter

**Recommendation**: Phred+33 only (A), document assumption

---

### Question 4: NEON Optimization Priority?
When should we add NEON optimizations?

**Options**:
- A: Immediately for `reverse_complement()` (known hot path)
- B: After benchmarks show need (evidence-based)
- C: Never (let users profile first)

**Recommendation**: B (wait for evidence)

---

### Question 5: CLI Examples?
Should we create example CLI tools showing primitive composition?

**Options**:
- A: Yes, in `examples/` directory
- B: No, let users create their own
- C: External repository (separate project)

**Recommendation**: A (demonstrates intended use)

---

## Next Steps

**Before Implementation**:
1. ✅ Review this plan with user
2. ⏳ Answer open questions
3. ⏳ Finalize API signatures
4. ⏳ Create initial test cases
5. ⏳ Set up benchmark framework (criterion)

**Start Implementation**:
1. Create `src/operations/sequence.rs`
2. Implement `reverse_complement()` scalar
3. Add proptest for reversibility
4. Verify tests pass

---

**Plan Status**: ⏳ AWAITING USER APPROVAL

