# Changelog

All notable changes to biometal will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-11-05

### 🎉 First Stable Release

Production-ready ARM-native bioinformatics library with streaming architecture and evidence-based optimization.

**Grade**: A+ (rust-code-quality-reviewer)
**Tests**: 121 passing (87 unit + 7 integration + 27 doc)
**Evidence Base**: 1,357 experiments, 40,710 measurements (N=30, 95% CI)

### Added

#### Core Features
- **Streaming FASTQ/FASTA parsers** with constant ~5 MB memory footprint
- **ARM NEON SIMD acceleration** (16-25× speedup on Mac ARM, Graviton)
  - Base counting: 16.7× speedup (5,254 Kseq/s vs 315 Kseq/s)
  - GC content: 20.3× speedup (5,954 Kseq/s vs 294 Kseq/s)
  - Quality filtering: 25.1× speedup (6,143 Kseq/s vs 245 Kseq/s)
- **Network streaming** with HTTP/HTTPS support (analyze without downloading)
- **SRA integration** (stream directly from NCBI SRA)
- **Python bindings** (PyO3 0.27, Python 3.9-3.14 support)

#### Optimizations (Evidence-Based)
- **Rule 1**: ARM NEON SIMD (Entry 020-025, 16-25× speedup)
- **Rule 2**: Block-based processing (Entry 027, 10K records preserves NEON)
- **Rule 3**: Parallel bgzip decompression (Entry 029, 6.5× speedup)
- **Rule 4**: Smart mmap for files ≥50 MB (Entry 032, 2.5× additional)
- **Rule 5**: Streaming architecture (Entry 026, 99.5% memory reduction)
- **Rule 6**: Network streaming (Entry 028, 264-352× I/O dominance)

#### Network Streaming Features
- Smart LRU cache (configurable 1 MB - 10 GB, default 50 MB)
- Background prefetching (hides network latency)
- Bounded worker pool (4 workers, prevents resource exhaustion)
- Request deduplication (eliminates duplicate fetches)
- Graceful cache poisoning recovery
- Automatic retry with exponential backoff

#### Python API
- `FastqStream.from_path()` - Stream FASTQ files
- `FastaStream.from_path()` - Stream FASTA files
- `gc_content()` - Calculate GC content (ARM NEON accelerated)
- `count_bases()` - Count A/C/G/T bases (ARM NEON accelerated)
- `mean_quality()` - Calculate mean Phred quality (ARM NEON accelerated)
- `extract_kmers()` - K-mer extraction for ML preprocessing

#### Documentation
- Comprehensive README with installation, usage, and examples
- Python installation guide and usage examples
- Jupyter notebook with interactive examples
- API documentation with performance characteristics
- Evidence links to experimental validation (OPTIMIZATION_RULES.md)
- Architecture documentation (docs/ARCHITECTURE.md)
- Performance tuning guide (docs/PERFORMANCE_TUNING.md)

### Cross-Platform Support

**Tested and Validated** (November 2025):
- ✅ **Mac ARM** (M1/M2/M3/M4): 121/121 tests pass, 16-25× NEON speedup (OPTIMIZED)
- ✅ **AWS Graviton 3**: 121/121 tests pass, 6-10× NEON speedup (PORTABLE)
- ✅ **x86_64 Intel**: 118/118 tests pass, 1× scalar fallback (PORTABLE)

**Strategy**: Optimized for Mac ARM (consumer hardware democratization), other platforms supported with correct, production-ready code.

### Performance

#### Memory Efficiency
- 99.5% memory reduction vs batch processing
- Constant ~5 MB footprint (even for TB-scale files)
- Streaming enables analysis on consumer hardware

#### Throughput (Mac M3 Max)
- Base counting: 5,254 Kseq/s (16.7× vs scalar)
- GC content: 5,954 Kseq/s (20.3× vs scalar)
- Quality filtering: 6,143 Kseq/s (25.1× vs scalar)

#### I/O Optimization
- Parallel bgzip: 6.5× speedup
- Smart mmap (≥50 MB): 2.5× additional speedup
- Combined: 16.3× I/O improvement

### Code Quality

**All 8 quality improvements completed**:
1. ✅ SRA URL pattern fix (CRITICAL)
2. ✅ SRA format limitation docs (HIGH)
3. ✅ Bounded thread pool for prefetch (HIGH)
4. ✅ Graceful cache poisoning recovery (MEDIUM)
5. ✅ Cache size validation (LOW)
6. ✅ Request deduplication (MEDIUM)
7. ✅ Enhanced documentation (LOW)
8. ✅ Additional test coverage (LOW)

**Quality Metrics**:
- Zero clippy warnings
- Zero unsafe blocks (outside NEON intrinsics)
- Complete error handling (no panics in library code)
- Property-based testing with proptest
- Comprehensive benchmarks with criterion (N=30)

### Fixed

- PyO3 compatibility with Python 3.14 (upgraded to PyO3 0.27)
- Unnecessary reference deref in compression.rs:98
- Deprecated io::Error pattern → io::Error::other (Rust 1.77+)
- Missing documentation for Python record fields

### Changed

- Updated `PyObject` → `Py<PyAny>` (PyO3 0.27 best practice)
- Platform strategy: Mac ARM optimized, others portable (documented)

### Dependencies

- PyO3: 0.22 → 0.27 (Python 3.14 support)
- All other dependencies at stable versions

## Evidence Base

All optimizations validated through [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench):
- 1,357 experiments
- 40,710 measurements (N=30)
- 95% confidence intervals
- Cohen's d effect sizes
- Lab notebook with 33 entries

See [OPTIMIZATION_RULES.md](OPTIMIZATION_RULES.md) for detailed evidence links.

## Publications (In Preparation)

1. **DAG Framework**: BMC Bioinformatics
2. **biometal Library**: Bioinformatics (Application Note) or JOSS
3. **Four-Pillar Democratization**: GigaScience

## Acknowledgments

Built with evidence-based optimization principles and rigorous experimental validation.

**Reviewed by**: rust-code-quality-reviewer (Grade A+)
**Mission**: Democratizing bioinformatics compute for LMIC researchers, small labs, students, and field researchers.

---

[1.0.0]: https://github.com/shandley/biometal/releases/tag/v1.0.0
