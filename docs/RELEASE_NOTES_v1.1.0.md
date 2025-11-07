# biometal v1.1.0 - K-mer Operations & Complexity Scoring

**Release Date**: November 6, 2025
**Grade**: A+ (rust-code-quality-reviewer)

---

## 🧬 What's New

### K-mer Operations (ASBB Entry 034)

Production-ready k-mer operations with evidence-based design showing that k-mer operations are **data-structure-bound** (hash+HashMap), not compute-bound.

**New Functions**:
- `extract_kmers(sequence, k)` - Simple overlapping k-mer extraction
- `kmer_iter(sequence, k)` - Zero-copy streaming iterator
- `extract_minimizers(sequence, k, w)` - minimap2-style sketching
- `kmer_spectrum(sequences, k)` - K-mer frequency counting
- `KmerExtractor::with_parallel(threads)` - Opt-in 2.2× speedup for large datasets

**Python Bindings**:
```python
import biometal

# Simple extraction
kmers = biometal.extract_kmers(b"ATGCAT", k=3)
# [b'ATG', b'TGC', b'GCA', b'CAT']

# Minimizers (minimap2-style)
minimizers = biometal.extract_minimizers(b"ATGCATGC", k=3, w=5)

# Frequency spectrum
spectrum = biometal.kmer_spectrum([b"ATGCAT", b"GCATGC"], k=3)

# Parallel extraction (2.2× speedup for ≥1000 sequences)
extractor = biometal.KmerExtractor(parallel=True, threads=4)
kmers = extractor.extract(large_sequences, k=6)
```

**Performance** (Entry 034 findings):
- **Minimizers**: Scalar-only optimal (1.02-1.26× parallel, below ≥5× threshold)
- **Spectrum**: Scalar-only optimal (parallel sometimes SLOWER: 0.95-1.88×)
- **Extraction**: 2.19-2.38× with Parallel-4t for ≥1000 sequences

**Use Cases**:
- BERT/DNABert preprocessing (k=3-6)
- minimap2-style indexing
- De Bruijn graph assembly
- K-mer frequency analysis

---

### Complexity Scoring

Shannon entropy calculation for sequence complexity analysis.

**New Function**:
- `complexity_score(sequence)` - Returns 0.0-1.0 complexity score
- NEON-optimized (reuses `count_bases()` for 16.7× speedup)

**Python Binding**:
```python
import biometal

# Filter low-complexity regions (metagenomics QC)
sequence = b"ATGCATGCATGCATGC"
if biometal.complexity_score(sequence) > 0.5:
    # High enough complexity - process this sequence
    pass
```

**Use Cases**:
- Metagenomics QC (filter low-complexity regions)
- Homopolymer detection
- Read quality assessment

---

### Record Operations Extension

**New Function**:
- `to_fasta_record(fastq_record)` - Convert FASTQ to FASTA (drop quality scores)

---

## 🔧 Code Quality Improvements

Addressed **all 8 issues** from rust-code-quality-reviewer:

**HIGH Priority** (1):
- ✅ Removed unused HashMap import in Python bindings

**MEDIUM Priority** (4):
- ✅ Fixed unwrap() in documentation example
- ✅ Exposed `PARALLEL_THRESHOLD` constant (1000 sequences)
- ✅ Added `will_use_parallel()` method for API transparency
- ✅ Added 6 property-based tests validating k-mer invariants

**LOW Priority** (3):
- ✅ Improved thread pool error handling (fallback to scalar on failure)
- ✅ Removed unused `complexity_score_neon()` function
- ✅ Added comprehensive k-mer benchmarks

---

## 🧪 Testing

**Test Count**: 260 passing (254 unit/integration + 6 property-based)

**New Property-Based Tests** (using proptest):
- `prop_kmer_count` - Validates k-mer count formula
- `prop_iter_matches_extract` - Iterator/extract equivalence
- `prop_spectrum_sum` - Spectrum frequencies sum correctly
- `prop_minimizer_uniqueness` - Consecutive minimizers differ
- `prop_parallel_matches_scalar` - Parallel extraction correctness
- `prop_edge_cases_no_panic` - Robustness testing

**New Benchmarks**:
- `benches/kmer_operations.rs` (182 LOC, 6 benchmark groups)
  - K-mer extraction across k values
  - Sequence sizes (100 bp - 100K bp)
  - Minimizer extraction
  - Spectrum counting
  - Parallel vs scalar comparison
  - Operations comparison at 150bp

---

## 📊 Evidence Validation

**Entry 034 Key Findings**:
- K-mer operations are data-structure-bound (not compute-bound)
- Hash computation: 50-60% of runtime (sequential, can't vectorize)
- HashMap operations: 30-40% of runtime (thread contention with parallel)
- Base validation: 5-10% (only NEON-friendly part, minimal impact)

**Validates Existing Tools**:
- minimap2's scalar minimizer design (empirically confirmed)
- Identifies 2.2× speedup opportunity for DNABert k-mer extraction

---

## 📦 Installation

### Rust

```toml
[dependencies]
biometal = "1.1"
```

### Python

```bash
pip install biometal-rs==1.1.0
```

> **Note**: Package name is `biometal-rs` on PyPI, but import as `biometal` in Python.

---

## 🔗 Links

- **Crates.io**: https://crates.io/crates/biometal
- **PyPI**: https://pypi.org/project/biometal-rs/
- **Documentation**: https://docs.rs/biometal
- **Repository**: https://github.com/shandley/biometal
- **ASBB Entry 034**: [K-mer Operations Analysis](https://github.com/shandley/apple-silicon-bio-bench/blob/main/lab-notebook/2025-11/20251106-034-EXPERIMENT-kmer-operations.md)

---

## 🙏 Acknowledgments

This release implements evidence-based k-mer operations validated through systematic benchmarking (ASBB Entry 034: 1,357+ experiments, 40,710 measurements, N=30, 95% CI).

---

**Full Changelog**: https://github.com/shandley/biometal/blob/main/CHANGELOG.md#110---2025-11-06
