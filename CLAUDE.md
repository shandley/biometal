# biometal: Claude Development Guide

**Project**: biometal - ARM-native bioinformatics library
**Latest Release**: v1.9.0 (November 14, 2025)
**Current Focus**: Indexed random access (FAI and TBI for efficient region queries)
**Phase**: Phase 2 Format Library Sprint - Week 2 complete (indices + integration examples)

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

## Current Status (v1.9.0)

### Released Features
- **FASTQ/FASTA** streaming parsers (constant memory)
- **ARM NEON** operations (base counting, GC content, quality filtering)
- **Sequence manipulation** (reverse_complement, trimming, masking)
- **K-mer operations** (extraction, minimizers, spectrum)
- **Network streaming** (HTTP, SRA)
- **BAM/SAM parser** (v1.4.0, November 8, 2025 - production-ready)
  - Full BAM parsing (header, records, CIGAR, tags, sequences)
  - ARM NEON sequence decoding (+27.5% BAM parsing speedup)
  - SAM writing for downstream tools
  - Robustness features (oversized CIGAR, malformed record handling)
  - 70 tests passing (integration complete)
- **Python bindings** (PyO3 0.27, 70+ functions)
  - Full BAM/FASTQ/FASTA support
  - CIGAR operations, SAM writing
  - BED, GFA, VCF, GFF3 format support
  - FAI and TBI index support
- **Format Library (v1.8.0)**: Production-ready parsers for genomic annotation and assembly formats
  - **BED (Browser Extensible Data)**: Genomic intervals (BED3/6/12 support)
  - **GFA (Graphical Fragment Assembly)**: Assembly graphs (Segment/Link/Path records)
  - **VCF (Variant Call Format)**: Genetic variants (VCF 4.2 spec compliance)
  - **GFF3 (General Feature Format)**: Hierarchical gene annotations
  - **Python gzip support**: All 4 formats support `.gz` files in Python
  - **Property-based testing**: 23 tests validating format invariants
  - **Real-world validation**: 6 integration tests with ENCODE, UCSC, Ensembl, 1000 Genomes data
- **Index Support (v1.9.0)**: Efficient random access to genomic data
  - **FAI (FASTA Index)**: O(1) random access to sequences and regions
    - Build, load, and query FASTA indices
    - Compatible with samtools faidx format
    - ~200 bytes per sequence in memory
  - **TBI (Tabix Index)**: O(log n) region queries on BGZF-compressed files
    - Support for VCF, BED, GFF3 indices
    - Hierarchical binning for fast queries
    - 100-1000× speedup vs sequential scanning
    - Compatible with samtools tabix format
  - **Integration examples**: 3 Rust examples + 1 Python notebook
  - **29 tests**: FAI and TBI unit + integration tests
- **Tests**: 890 passing (670 unit + 211 doc + 23 property-based + 6 real-world integration)

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
- **PyPI**: biometal-rs v1.9.0 (pip install biometal-rs)
- **crates.io**: biometal v1.9.0 (cargo add biometal)

---


## Current Work: Phase 2 Format Library Sprint

**Status**: ✅ **Week 2 Complete** - Index support shipped in v1.9.0

### Completed (v1.9.0)
- ✅ BED format (BED3/6/12 variants, genomic intervals) - v1.8.0
- ✅ GFA format (assembly graphs, Segment/Link/Path) - v1.8.0
- ✅ VCF format (VCF 4.2 spec, variant calling) - v1.8.0
- ✅ GFF3 format (hierarchical gene annotations) - v1.8.0
- ✅ Python gzip support for all 4 formats - v1.8.0
- ✅ Property-based testing framework (23 tests) - v1.8.0
- ✅ Real-world data validation (6 integration tests) - v1.8.0
- ✅ **FAI (FASTA Index)**: Build, load, query sequences/regions - v1.9.0
- ✅ **TBI (Tabix Index)**: Region queries on BGZF files - v1.9.0
- ✅ **Integration examples**: 3 Rust + 1 Python notebook - v1.9.0
- ✅ **Python bindings**: FAI and TBI support - v1.9.0

### In Progress
- Additional formats (GTF, PAF, narrowPeak) - Week 3-6
- Cross-format integration and workflows - Week 7-10
- Python enhancements and documentation - Week 11-14
- Performance benchmarking and v2.0.0 release - Week 15-16

### Strategic Context

**Decision Made**: Format Library Expansion (STRATEGIC_PIVOT_PLAN.md Option 1)
- Build useful bioinformatics primitives with proven Apple Silicon optimizations
- Focus on streaming architecture (constant memory) + ARM NEON acceleration
- Evidence-based design (all optimizations validated via apple-silicon-bio-bench)

**Historical Context**:
- Phase 1 consolidation complete (documentation, benchmarking, quality)
- Apple Silicon research complete (modest results, archived)
- Rules 3+4 found non-viable for streaming architecture

**See**: [PHASE2_FORMAT_LIBRARY_SPRINT.md](PHASE2_FORMAT_LIBRARY_SPRINT.md) for detailed 16-week plan

---

## Future Work

**Status**: Pending strategic direction decision

**Options** (see STRATEGIC_PIVOT_PLAN.md):
- Format expansion (CRAM, VCF, BCF, CSI)
- Community building (blog, social media, adoption)
- Quality assurance (testing, cross-platform validation)
- Maintenance mode (wait for community feedback)

**Decision Required**: Which direction should biometal pursue?

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
- Don't make up optimization parameters (refer to evidence)
- Don't accumulate records in memory (streaming only)
- Don't panic in library code (use Result)
- Don't implement ARM-only code without scalar fallback

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

---

## Performance Expectations

### Current Performance (v1.7.0)

| Operation | Scalar | Optimized | Speedup |
|-----------|--------|-----------|---------|
| Base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7× (NEON)** |
| GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3× (NEON)** |
| Quality filter | 245 Kseq/s | 6,143 Kseq/s | **25.1× (NEON)** |
| BAM parsing | ~11 MiB/s | 92.0 MiB/s | **8.4× (BGZF + NEON + cloudflare_zlib)** |
| BAM indexed query | O(n) full scan | O(log n) indexed | **1.68-500× (scales with file size)** |

**Note**: No major performance optimizations remain viable:
- Rule 3 (Parallel BGZF): Conflicts with streaming architecture (0.77-0.84× slowdown)
- Rule 4 (Smart mmap): ~1% benefit for CPU-bound decompression

**Focus**: Feature expansion (CRAM, VCF) and community adoption

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
- ✅ FASTQ, FASTA (v1.0.0)
- ✅ BAM, SAM (v1.4.0)
- ✅ BAI index (v1.6.0)
- ✅ BED (BED3/6/12) (v1.8.0)
- ✅ GFA (Segment/Link/Path) (v1.8.0)
- ✅ VCF (VCF 4.2) (v1.8.0)
- ✅ GFF3 (hierarchical annotations) (v1.8.0)
- ⏳ FAI, TBI indices (Phase 2 Week 2-3)
- ⏳ GTF, PAF, narrowPeak (Phase 2 Week 4-8)
- ❌ CRAM, BCF (future)

### Tests
- **860 passing** (100% pass rate)
  - 649 unit tests (core + formats)
  - 211 documentation tests
  - 23 property-based tests (format invariants)
  - 6 real-world integration tests (ENCODE, UCSC, Ensembl, 1000 Genomes)

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

**Last Updated**: November 14, 2025 (Phase 2 Week 1 complete, v1.8.0 released)
**Current Phase**: Phase 2 Format Library Sprint (4 formats shipped, 12-15 weeks remaining)
**Next Milestone**: Format utilities (FAI, TBI indices) - Week 2-3
