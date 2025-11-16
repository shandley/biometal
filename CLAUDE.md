# biometal: Claude Development Guide

**Project**: biometal - ARM-native bioinformatics library
**Latest Release**: v1.11.0 (November 16, 2025)
**Current Focus**: Phase 2 COMPLETE - Planning next phase
**Phase**: Post-Phase 2 - Determining next direction

---

## ⚠️ CRITICAL CONSTRAINTS (READ THIS FIRST!)

### DO NOT Recommend Performance Optimization (Rules 3+4)

❌ **Rule 3 (Parallel BGZF)**: 0.77-0.84× **SLOWDOWN** - PERMANENTLY DISABLED
- Conflicts with streaming architecture (Rule 5)
- Experimentally validated as harmful (Nov 2025)
- See STRATEGIC_PIVOT_PLAN.md for details

❌ **Rule 4 (Smart mmap)**: ~1% benefit - NOT WORTH IT
- Minimal impact for CPU-bound decompression
- Complexity not justified by tiny gains

### Phase 2 Complete: FORMAT LIBRARY Achieved

✅ **What biometal HAS**: 14 formats with READ+WRITE, 4 READ-only formats, 3 index formats
✅ **What biometal is NOT doing**: Performance optimization (Rules 3+4 disabled)

**Performance is already strong**: 92 MiB/s BAM parsing, 16-25× NEON speedup, constant 5 MB memory

**See**: docs/archive/2025-11/phase2-format-sprint-plan.md for completed plan

---

## Mission

Democratize bioinformatics by enabling 5TB dataset analysis on consumer hardware through:
- **Streaming architecture**: Constant ~5 MB memory (not load-all)
- **ARM-native performance**: 16-25× NEON speedup
- **Network streaming**: Analyze without downloading
- **Evidence-based optimization**: Every rule validated experimentally

**Target**: LMIC researchers, small labs, students, field researchers, ML practitioners

---

## Core Principles

### 1. Evidence-Based Design

Every optimization comes from validated experimental results (apple-silicon-bio-bench):
- Follow OPTIMIZATION_RULES.md: 6 rules from 1,357 experiments (N=30)
- Don't guess: Reference ASBB evidence for optimization decisions
- Document rationale: Link implementations to specific rules/entries

### 2. Streaming-First Architecture

Always design for constant memory:
- Bad: `Vec<FastqRecord>` (accumulates in memory)
- Good: `Iterator<Item = FastqRecord>` (constant memory)
- Target: ~5 MB regardless of dataset size (Rule 5)

### 3. ARM-Native with Portable Fallback

Always provide both ARM and fallback implementations:
```rust
#[cfg(target_arch = "aarch64")]
pub fn operation_neon(input: &[u8]) -> Result { /* 16-25× faster */ }

#[cfg(not(target_arch = "aarch64"))]
pub fn operation_scalar(input: &[u8]) -> Result { /* x86_64 fallback */ }

pub fn operation(input: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { operation_neon(input) }
    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(input) }
}
```

Platform priority: Mac ARM → Linux ARM (Graviton) → x86_64 fallback

### 4. Production Quality

- Use `Result<T, BiometalError>` (no panics in library)
- Document every public API with examples
- Property-based testing (proptest)
- Benchmarks (criterion, N=30)
- No `unwrap()` or `expect()` in library code

---


## Recent Research

### Apple Silicon Exploration (November 4-13, 2025) - [Full Report](docs/archive/2025-11/apple-silicon-research.md)
Explored Neural Engine and GPU acceleration for bioinformatics workloads. Results: 1.2-1.4× speedup for GPU batch processing, 2,940× slowdown for Neural Engine streaming. Finding: Architectural mismatch between batch-oriented hardware and streaming-first design. Status: Archived to research/apple-silicon/ (infrastructure preserved, feature-gated).

---

### Compression Optimization (v1.7.0) - [Details](docs/archive/2025-11/compression-optimization.md)
Switched to cloudflare_zlib backend: 1.67× decompression, 2.29× compression speedup. BAM parsing improved: 55 MiB/s → 92 MiB/s (+67%). See [planning_archive/investigations/BACKEND_COMPARISON_FINDINGS.md](planning_archive/investigations/BACKEND_COMPARISON_FINDINGS.md) for full analysis.

---

## Current Status (v1.11.0)

### Released Features
- **FASTQ/FASTA** streaming parsers (constant memory)
- **ARM NEON** operations (base counting, GC content, quality filtering)
- **Sequence manipulation** (reverse_complement, trimming, masking)
- **K-mer operations** (extraction, minimizers, spectrum)
- **Network streaming** (HTTP, SRA)
- **BAM/SAM parser** (v1.4.0, production-ready)
  - Full BAM parsing (header, records, CIGAR, tags, sequences)
  - ARM NEON sequence decoding (+27.5% BAM parsing speedup)
  - SAM writing for downstream tools
  - Robustness features (oversized CIGAR, malformed record handling)
- **CRAM reader** (v1.12.0, production-ready)
  - Full CRAM 3.0/3.1 decoding (reference-based compression)
  - Native ARM-optimized implementation (zero external dependencies except codecs)
  - Multi-codec support (gzip, bzip2, lzma, rANS 4x16)
  - Streaming architecture (constant ~5 MB memory)
  - 46 CRAM tests passing, validated with 30K+ real-world records
  - ARM NEON optimizations (9× base counting, ~10% overall parsing improvement)
- **Python bindings** (PyO3 0.27, 100+ functions, **complete format coverage**)
  - **READ + WRITE**: FASTQ, FASTA, BAM, SAM, BED (3/6/12 + narrowPeak), GFA, VCF, GFF3, GTF, PAF
  - **READ only**: CRAM, GenBank, BLAST tabular
  - CIGAR operations, BAM/SAM headers, tags
  - Index support (BAI, CSI, FAI, TBI) with region queries
  - Statistics functions (coverage, insert size, MAPQ, edit distance, strand bias)
  - Sequence operations (reverse complement, trimming, masking, validation)
  - K-mer operations (extraction, minimizers, spectrum)
  - PyO3 registration bug resolved (Nov 16, 2025) - BamWriter, SamReader, GfaWriter now working
- **Format Library (READ + WRITE)**: Production-ready streaming parsers and writers
  - **BED (BED3/6/12 + narrowPeak)**: Genomic intervals and ChIP-seq peaks (v1.8.0, v1.10.0)
  - **GFA (Graphical Fragment Assembly)**: Assembly graphs (v1.8.0)
  - **VCF (Variant Call Format)**: Genetic variants, VCF 4.2 spec (v1.8.0)
  - **GFF3 (General Feature Format)**: Hierarchical gene annotations (v1.8.0)
  - **GTF (Gene Transfer Format)**: RNA-seq gene annotations (v1.10.0)
  - **PAF (Pairwise mApping Format)**: minimap2 long-read alignments (v1.10.0)
- **Format Library (READ-only)**: Specialized formats
  - **GenBank**: NCBI annotated sequences with features (LOCUS, FEATURES, ORIGIN) (v1.11.0)
  - **BLAST tabular**: Sequence alignment results (outfmt 6/7, 12 standard columns) (v1.11.0)
  - **Python optimizations**: 50-60% memory reduction per record (v1.10.0)
  - **Property-based testing**: 23 tests validating format invariants
  - **Real-world validation**: 6 integration tests with ENCODE, UCSC, Ensembl, 1000 Genomes
  - **Writer support**: All formats support streaming write with compression (Nov 16, 2025)
- **Index Support (v1.9.0)**: Efficient random access
  - **FAI (FASTA Index)**: O(1) sequence/region access, samtools-compatible
  - **TBI (Tabix Index)**: O(log n) region queries on BGZF files (VCF, BED, GFF3)
  - 100-1000× speedup vs sequential scanning
  - Integration examples: 3 Rust + 1 Python notebook
- **Tests**: 670 total (669 passing + 1 ignored, 100% pass rate)
  - 669 library tests (includes 46 CRAM, 9 GenBank/BLAST tests with NEON optimizations)
  - See CHANGELOG.md for full breakdown

### Optimization Rules Implemented

| Rule | Feature | Status | Impact |
|------|---------|--------|--------|
| **Rule 1** | ARM NEON SIMD | ✅ v1.0.0 | 16-25× speedup |
| **Rule 2** | Block-based processing | ✅ v1.0.0 | Preserves NEON gains |
| **Rule 3** | Parallel BGZF | ❌ Disabled | Conflicts with streaming |
| **Rule 4** | Smart mmap | ⏳ Optional | ~1% for compressed files |
| **Rule 5** | Constant-memory streaming | ✅ v1.0.0 | 99.5% memory reduction |
| **Rule 6** | Network streaming | ✅ v1.0.0 | Enables remote analysis |

**Current**: 4/6 rules implemented (Rules 1, 2, 5, 6)
**Rule 3**: Not viable for streaming architecture (see OPTIMIZATION_RULES.md)
**Rule 4**: Minimal benefit (~1%) for CPU-bound decompression workloads

### Distribution
- **PyPI**: biometal-rs v1.11.0 (pip install biometal-rs)
- **crates.io**: biometal v1.11.0 (cargo add biometal)

---


## Phase 2 Complete: Format Library (v1.8.0 - v1.11.0)

**Status**: ✅ **COMPLETE** (November 16, 2025)

### Achievements
- **14 formats with READ + WRITE**: FASTQ, FASTA, BAM, SAM, BED (3/6/12 + narrowPeak), GFA, VCF, GFF3, GTF, PAF
- **4 formats READ-only**: CRAM, GenBank, BLAST tabular, CSI index
- **3 index formats**: BAI, FAI, TBI (random access, 100-1000× speedup)
- **670 tests**: 669 passing + 1 ignored (100% pass rate)
  - Property-based testing (23 tests validating format invariants)
  - Real-world validation (6 integration tests with ENCODE, UCSC, Ensembl, 1000 Genomes)
- **Python bindings**: Complete feature parity (100+ functions)
  - All read/write formats supported
  - Index support with region queries
  - Statistics and sequence operations
- **Performance**: 92 MiB/s BAM parsing, 16-25× NEON speedup, constant 5 MB memory

### Timeline
**Duration**: 3.5 months (August - November 2025)
**Key Milestones**:
- v1.8.0 (Nov 13): BED, GFA, VCF, GFF3 parsers
- v1.9.0 (Nov 14): FAI/TBI indices + file writing support
- v1.10.0 (Nov 14): GTF, PAF, narrowPeak + Python optimizations
- v1.11.0 (Nov 16): GenBank, BLAST parsers + code quality (Grade A)

**See**: [docs/archive/2025-11/phase2-format-sprint-plan.md](docs/archive/2025-11/phase2-format-sprint-plan.md) for full plan

---

## Recent Milestones (November 2025)

- **v1.11.0** (Nov 16): GenBank + BLAST parsers, code quality Grade A
- **v1.10.0** (Nov 14): GTF, PAF, narrowPeak formats + Python optimizations (50-60% memory reduction)
- **v1.9.0** (Nov 14): FAI/TBI indices + file writing support (FASTQ, FASTA, BED, GTF, VCF)
- **Tab-delimited writers** (Nov 16): BED, GFF3, GTF writers with compression + validation
- **PyO3 bug fix** (Nov 16): Resolved registration issue affecting BamWriter, SamReader, GfaWriter
- **CRAM Python bindings** (Nov 15): PyCramReader with full streaming support

**See**: [docs/archive/2025-11/sessions/](docs/archive/2025-11/sessions/) for detailed session notes

---


## Next Steps: Post-Phase 2

**Phase 2 Status**: ✅ COMPLETE (v1.11.0, November 16, 2025)

### Option 1: Code Quality & Polish (1-2 weeks, **RECOMMENDED**)
**Why**: Makes biometal production-ready for community adoption
- Comprehensive benchmarking (N=30 validation for all formats)
- Cross-platform testing (AWS Graviton ARM, Linux x86_64)
- Documentation improvements (user guides, tutorials, API docs)
- Performance profiling and optimization opportunities
- Blog post/paper draft showcasing Phase 2 achievements

### Option 2: Advanced Format Support (2-4 weeks)
**Why**: Comprehensive format coverage
- **BCF** (Binary VCF) - 30-40 hours, high community value
- **CRAM writing** (complex, 40-60 hours, lower priority)
- Additional indices (CSI completion, other formats)
**Risk**: Diminishing returns, most common formats done

### Option 3: GPU/ML Research (3-6 months, **HIGH NOVELTY**)
**Why**: High research impact, unique to Apple Silicon
- **Smith-Waterman GPU acceleration** (10-50× potential speedup)
- **Streaming BERT data loaders** (constant memory for genomic ML)
- **Quality-aware tokenization** (NOVEL - publishable research)
- **Neural Engine integration** (Apple Silicon exclusive)
**See**: PROJECT_TODOS.md for detailed 24-week roadmap

### Recommendation
**Phase 3**: Option 1 (Polish, 1-2 weeks) → Option 3 (GPU/ML, 3-6 months) → v2.0 release
**Rationale**: Polish ensures quality foundation, GPU/ML has highest research impact

**Decision Point**: Which option to pursue next session?

---

## Key Documentation

### For Users
- **📘 User Guide**: docs/USER_GUIDE.md - Comprehensive onboarding (installation → optimization)
- **📓 Tutorials**: notebooks/ - 7 Jupyter notebooks (including BAI indexed queries)
- **⚡ Performance Guide**: docs/PERFORMANCE_OPTIMIZATION_GUIDE.md - Maximize performance
- **📊 Benchmarks**: benchmarks/comparison/BENCHMARK_COMPARISON.md - vs samtools/pysam

### For Developers
- **📐 Architecture**: docs/ARCHITECTURE.md - Technical design
- **🔬 Optimization Rules**: OPTIMIZATION_RULES.md - Evidence-based optimization (6 rules)
- **🐍 Python Bindings**: docs/PYTHON.md - Python-specific details
- **🧬 BAM API**: docs/BAM_API.md - Complete BAM/SAM parser reference

### For Planning
- **📈 Phase 1 Progress**: PHASE1_PROGRESS_REPORT.md - Current consolidation status
- **🗺️ Strategic Analysis**: NEXT_STEPS_ANALYSIS.md - Long-term roadmap (3-phase, 14 weeks)
- **📝 Changelog**: CHANGELOG.md - Version history

---

## Development Workflow

### Session Guidelines

**What to Emphasize**:
- Evidence-based design (follow OPTIMIZATION_RULES.md)
- Streaming-first architecture (constant memory)
- ARM-native with portable fallback
- Production quality (error handling, docs, tests)

**What NOT to Do**:
- ❌ **NEVER recommend Rules 3+4 performance optimization** (disabled, see top of document)
- ❌ **NEVER suggest parallel BGZF or mmap for performance gains**
- ❌ Don't make up optimization parameters (refer to evidence)
- ❌ Don't accumulate records in memory (streaming only)
- ❌ Don't panic in library code (use Result)
- ❌ Don't implement ARM-only code without scalar fallback

**Decision Framework**:
When implementing features:
1. Check OPTIMIZATION_RULES.md for relevant rule
2. Follow the implementation pattern for that rule
3. Link to evidence (lab notebook entry)

When evaluating optimizations:
1. Is this validated in ASBB experiments?
2. If yes: Which rule/entry documents it?
3. If no: Suggest validating first or using proven approach

### Testing Strategy

**Property-Based Testing**:
```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn test_base_counting_matches_naive(seq in "[ACGT]{1,1000}") {
        let neon_result = count_bases_neon(seq.as_bytes());
        let naive_result = count_bases_naive(seq.as_bytes());
        prop_assert_eq!(neon_result, naive_result);
    }
}
```

**Benchmarking**:
```rust
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_operation(c: &mut Criterion) {
    let data = generate_test_data(100_000);
    c.bench_function("operation_neon", |b| {
        b.iter(|| operation_neon(&data))
    });
}

criterion_group!(benches, bench_operation);
criterion_main!(benches);
```

### Python Bindings Development

**⚠️ CRITICAL: maturin develop Issue**

When adding or updating Python bindings, `maturin develop` does **NOT** properly update the `.so` file in site-packages, even though it reports success. This was discovered on Nov 14, 2025 during FastqWriter/FastaWriter implementation.

**The Problem**:
- `maturin develop --release --features python` builds successfully
- Reports "✏️ Setting installed package as editable"
- BUT: The `.so` file in site-packages is NOT updated
- Result: New Python classes don't appear in the module

**The Solution** (Required after every Python binding change):

**Option 1: Use the helper script** (recommended):
```bash
# Builds, copies .so, and verifies the class exists
./.claude/scripts/update_python_bindings.sh YourNewClass
```

**Option 2: Manual steps**:
```bash
# 1. Build with maturin
maturin develop --release --features python

# 2. MANUALLY copy the .so file (REQUIRED!)
cp target/release/libbiometal.dylib \
   $(python3 -c "import site; print(site.getusersitepackages())")/biometal/biometal.cpython-314-darwin.so

# 3. Verify the update worked
python3 -c "import biometal; print(hasattr(biometal, 'YourNewClass'))"
```

**Quick verification script**:
```bash
# Check .so file timestamp
stat -f "Modified: %Sm" $(python3 -c "import site; print(site.getsitepackages()[0])")/biometal/biometal.cpython-314-darwin.so

# Should match or be newer than target build
stat -f "Modified: %Sm" target/release/libbiometal.dylib
```

**Checklist for Python Binding Updates**:
1. ✅ Implement PyO3 wrapper class with `#[pyclass(name = "ClassName", unsendable)]`
2. ✅ Add methods with `#[pymethods]`
3. ✅ Register in `src/python/mod.rs` with `m.add_class::<PyClassName>()?`
4. ✅ Build: `maturin develop --release --features python`
5. ✅ **MANUALLY copy .so file** (critical step!)
6. ✅ Verify: `python3 -c "import biometal; print(hasattr(biometal, 'ClassName'))"`
7. ✅ Test round-trip functionality

**Example Python Binding Pattern** (FastaWriter):
```rust
// src/python/streams.rs
#[pyclass(name = "FastaWriter", unsendable)]
pub struct PyFastaWriter {
    inner: Option<FastaWriter>,
}

#[pymethods]
impl PyFastaWriter {
    #[staticmethod]
    fn create(path: String) -> PyResult<Self> {
        let writer = FastaWriter::create(&PathBuf::from(path))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(PyFastaWriter { inner: Some(writer) })
    }

    fn write_record(&mut self, record: &PyFastaRecord) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            writer.write_record(&record.to_fasta_record())
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyIOError, _>("writer already finished"))
        }
    }
}

// src/python/mod.rs
m.add_class::<PyFastaWriter>()?;
```

---

## Performance Expectations

### Current Performance (v1.10.0)

| Operation | Scalar | Optimized | Speedup |
|-----------|--------|-----------|---------|
| Base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7× (NEON)** |
| GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3× (NEON)** |
| Quality filter | 245 Kseq/s | 6,143 Kseq/s | **25.1× (NEON)** |
| BAM parsing | ~11 MiB/s | 92.0 MiB/s | **8.4× (BGZF + NEON + cloudflare_zlib)** |
| BAM indexed query | O(n) full scan | O(log n) indexed | **1.68-500× (scales with file size)** |

**✅ Performance is COMPLETE**: No major optimizations remain viable
- ❌ Rule 3 (Parallel BGZF): 0.77-0.84× slowdown (conflicts with streaming)
- ❌ Rule 4 (Smart mmap): ~1% benefit (not worth complexity)

**Current Focus**: File format writing support (NOT performance optimization)

---

## Competitive Position (Validated Week 2)

### vs samtools

| Metric | biometal | samtools | Advantage |
|--------|----------|----------|-----------|
| Sequential parsing | 55.1 MiB/s | ~45-50 MiB/s | ✅ Competitive |
| Indexed queries | 1.68-500× | ~1.2-1.5× | ✅ **Superior** |
| Memory | **5 MB** | 20-50 MB | ✅ **10× lower** |
| ARM NEON | **4-25× speedup** | None | ✅ **Exclusive** |

### vs pysam

| Metric | biometal | pysam | Advantage |
|--------|----------|-------|-----------|
| Python performance | ~45 MiB/s | ~30-40 MiB/s | ✅ 1.5-2× faster |
| Memory | **5 MB** | 50 MB-1 GB | ✅ **10-200× lower** |
| API | Streaming | Context managers | ✅ Simpler |
| ARM NEON | **4-25× speedup** | None | ✅ **Exclusive** |

**Production Use Cases**:
- ✅ Large-file targeted analysis (indexed queries)
- ✅ Memory-constrained environments
- ✅ ARM infrastructure (Apple Silicon, Graviton)
- ✅ Terabyte-scale streaming

---

## Quick Reference

### Evidence Base
- **1,357 experiments**, 40,710 measurements (N=30)
- Source: apple-silicon-bio-bench
- Rules: 6 optimization rules (OPTIMIZATION_RULES.md)

### Platform Support
1. **Mac ARM** (M1/M2/M3/M4): 16-25× NEON speedup (optimized)
2. **Linux ARM** (Graviton): 6-10× NEON speedup (portable)
3. **x86_64**: 1× scalar fallback (portable)

### File Formats
**READ + WRITE** (Rust + Python):
- ✅ FASTQ, FASTA (v1.0.0 read, v1.9.0 write)
- ✅ BAM, SAM (v1.4.0 read, v1.7.0-v1.8.0 write)
- ✅ BED (BED3/6/12 + narrowPeak) (v1.8.0 read, v1.10.0+ write)
- ✅ GFA (Segment/Link/Path) (v1.8.0 read/write)
- ✅ VCF (VCF 4.2) (v1.8.0 read/write)
- ✅ GFF3 (hierarchical annotations) (v1.8.0 read, v1.10.0+ write)
- ✅ GTF (RNA-seq annotations) (v1.10.0 read/write)
- ✅ PAF (minimap2 alignments) (v1.10.0 read/write)

**READ only** (Rust + Python):
- ✅ CRAM (v1.12.0) - Production-ready decoder, ARM NEON optimized
- ✅ GenBank (v1.11.0) - NCBI annotated sequences (LOCUS, FEATURES, ORIGIN)
- ✅ BLAST tabular (v1.11.0) - Sequence alignment results (outfmt 6/7, 12 columns)

**Index support** (READ only):
- ✅ BAI, CSI (BAM/CRAM indices) (v1.6.0)
- ✅ FAI (FASTA index) (v1.9.0)
- ✅ TBI (Tabix index for VCF/BED/GFF) (v1.9.0)

**Not planned**:
- ❌ BCF (deferred - lower priority)
- ❌ CRAM writing (complex, lower priority)

### Tests
- **670 total tests** (669 passing + 1 ignored, 100% pass rate)
  - 46 CRAM tests with NEON optimizations
  - 9 GenBank/BLAST tests (v1.11.0)
  - 24 tab-delimited writer tests (BED: 3, GFF3: 9, GTF: 12)
  - Property-based tests (format invariants)
  - Real-world integration tests (ENCODE, UCSC, Ensembl, 1000 Genomes, pUC19, E. coli)
  - See CHANGELOG.md for full breakdown

---

## Session Checklist

When starting a new session:
- [ ] Review current phase status (PHASE1_PROGRESS_REPORT.md)
- [ ] Check CHANGELOG.md for recent changes
- [ ] Review relevant optimization rules (OPTIMIZATION_RULES.md)
- [ ] Check open issues/PRs if community-facing work

When implementing features:
- [ ] Follow evidence-based design (link to ASBB entry)
- [ ] Use streaming architecture (constant memory)
- [ ] Provide ARM + fallback implementations
- [ ] Add property-based tests
- [ ] Benchmark with criterion (N=30)
- [ ] Document with examples

When wrapping up:
- [ ] Update CHANGELOG.md
- [ ] Run full test suite
- [ ] Update relevant planning documents
- [ ] Document any decisions made

---

**Last Updated**: November 16, 2025 (v1.11.0 - Phase 2 COMPLETE)
**Current Phase**: Post-Phase 2 - Planning next direction
**Latest Release**: v1.11.0 (GenBank + BLAST parsers, Grade A code quality)
**Phase 2 Achievement**: 14 formats (READ+WRITE), 4 formats (READ-only), 3 indices, 670 tests, 100+ Python functions
**Next Decision**: Option 1 (Polish) vs Option 2 (Advanced Formats) vs Option 3 (GPU/ML Research)
