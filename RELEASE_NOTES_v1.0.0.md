# biometal v1.0.0 - Production Release 🎉

**Released**: November 5, 2025

Production-ready ARM-native bioinformatics library with streaming architecture and evidence-based optimization.

**Grade**: A+ (rust-code-quality-reviewer)
**Tests**: 121 passing (87 unit + 7 integration + 27 doc)
**Evidence Base**: 1,357 experiments, 40,710 measurements (N=30, 95% CI)

---

## 🚀 Installation

```bash
pip install biometal-rs
```

**Requirements**: Python 3.9+ (tested through 3.14)

---

## ✨ Core Features

### Streaming FASTQ/FASTA Parsers
- Constant ~5 MB memory footprint regardless of dataset size
- Analyze 5TB datasets on consumer hardware without downloading
- 99.5% memory reduction compared to batch processing

### ARM NEON SIMD Acceleration
16-25× speedup on Mac ARM (M1/M2/M3/M4), automatic fallback on x86_64:

| Operation | Scalar | NEON | Speedup |
|-----------|--------|------|---------|
| Base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7×** |
| GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3×** |
| Quality filtering | 245 Kseq/s | 6,143 Kseq/s | **25.1×** |

### Network Streaming
- Stream directly from HTTP/HTTPS sources (analyze without downloading)
- SRA integration (stream from NCBI SRA)
- Smart LRU caching with background prefetching
- Analyze terabyte-scale datasets on laptops

### Python Bindings
- PyO3 0.27 (Python 3.9-3.14 support)
- Streaming API: `FastqStream.from_path()`, `FastaStream.from_path()`
- Operations: `gc_content()`, `count_bases()`, `mean_quality()`, `extract_kmers()`
- Full ARM NEON acceleration from Python

---

## 🌍 Cross-Platform Support

**Tested and Validated** (November 2025):

| Platform | Performance | Status |
|----------|-------------|---------|
| **Mac ARM** (M1/M2/M3/M4) | 16-25× NEON speedup | ✅ Optimized |
| **AWS Graviton 3** | 6-10× NEON speedup | ✅ Portable |
| **Linux x86_64** | 1× scalar fallback | ✅ Portable |

**Strategy**: Optimized for Mac ARM (consumer hardware democratization), other platforms supported with correct, production-ready code.

---

## 📊 Evidence-Based Optimization

All optimizations validated through [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench):

- **Rule 1**: ARM NEON SIMD (Entry 020-025, 16-25× speedup)
- **Rule 2**: Block-based processing (Entry 027, 10K records preserves NEON)
- **Rule 3**: Parallel bgzip decompression (Entry 029, 6.5× speedup)
- **Rule 4**: Smart mmap for files ≥50 MB (Entry 032, 2.5× additional)
- **Rule 5**: Streaming architecture (Entry 026, 99.5% memory reduction)
- **Rule 6**: Network streaming (Entry 028, 264-352× I/O dominance)

See [OPTIMIZATION_RULES.md](https://github.com/shandley/biometal/blob/main/OPTIMIZATION_RULES.md) for detailed evidence links.

---

## 🎯 Democratizing Bioinformatics

biometal addresses four barriers that lock researchers out of genomics:

### 1. Economic Barrier
- **Problem**: Most tools require $50K+ servers
- **Solution**: Consumer ARM laptops ($1,400) deliver production performance
- **Impact**: Small labs and LMIC researchers can compete

### 2. Environmental Barrier
- **Problem**: HPC clusters consume massive energy
- **Solution**: ARM efficiency inherent in architecture
- **Impact**: Reduced carbon footprint for genomics research

### 3. Portability Barrier
- **Problem**: Vendor lock-in (x86-only, cloud-only tools)
- **Solution**: Works across ARM ecosystem (Mac, Graviton, Ampere, RPi)
- **Impact**: No platform dependencies, true portability

### 4. Data Access Barrier ⭐
- **Problem**: 5TB datasets require 5TB storage + days to download
- **Solution**: Network streaming with smart caching
- **Impact**: Analyze 5TB datasets on 24GB laptops without downloading

---

## 📚 Quick Start

### Python Example

```python
import biometal

# Stream FASTQ from local file (constant memory)
stream = biometal.FastqStream.from_path("large_dataset.fq.gz")

for record in stream:
    # ARM NEON-optimized operations (16-25× faster)
    gc = biometal.gc_content(record.sequence)
    counts = biometal.count_bases(record.sequence)
    mean_q = biometal.mean_quality(record.quality)

    print(f"{record.id}: GC={gc:.2%}, Mean Q={mean_q:.1f}")
    # Memory stays constant at ~5 MB
```

### Network Streaming Example

```python
import biometal

# Stream from NCBI SRA (no download required!)
stream = biometal.FastqStream.from_sra("SRR390728")  # E. coli dataset

for record in stream:
    # Process 40 MB dataset with only ~5 MB memory
    gc = biometal.gc_content(record.sequence)
    print(f"{record.id}: GC={gc:.2%}")
```

---

## 🏆 Code Quality

### All 8 Quality Improvements Completed

1. ✅ SRA URL pattern fix (CRITICAL)
2. ✅ SRA format limitation docs (HIGH)
3. ✅ Bounded thread pool for prefetch (HIGH)
4. ✅ Graceful cache poisoning recovery (MEDIUM)
5. ✅ Cache size validation (LOW)
6. ✅ Request deduplication (MEDIUM)
7. ✅ Enhanced documentation (LOW)
8. ✅ Additional test coverage (LOW)

### Quality Metrics
- Zero clippy warnings
- Zero unsafe blocks (outside NEON intrinsics)
- Complete error handling (no panics in library code)
- Property-based testing with proptest
- Comprehensive benchmarks with criterion (N=30)

---

## 📦 Available Wheels

This release includes pre-built wheels for:
- macOS ARM (M1/M2/M3/M4)
- macOS x86_64 (Intel Macs)
- Linux x86_64

Plus source distribution for other platforms.

**Note**: Linux ARM (Graviton, Raspberry Pi) wheels will be added in v1.0.1 after resolving cross-compilation setup. Linux ARM users can build from source in the meantime.

---

## 🔗 Resources

- **Documentation**: https://docs.rs/biometal
- **Repository**: https://github.com/shandley/biometal
- **Evidence Base**: https://github.com/shandley/apple-silicon-bio-bench
- **PyPI**: https://pypi.org/project/biometal-rs/

---

## 📝 Citation

If you use biometal in your research, please cite:

```bibtex
@software{biometal2025,
  author = {Handley, Scott},
  title = {biometal: ARM-native bioinformatics with streaming architecture},
  year = {2025},
  version = {1.0.0},
  url = {https://github.com/shandley/biometal}
}
```

---

## 🙏 Acknowledgments

Built with evidence-based optimization principles and rigorous experimental validation.

**Reviewed by**: rust-code-quality-reviewer (Grade A+)
**Mission**: Democratizing bioinformatics compute for LMIC researchers, small labs, students, and field researchers.

---

**Full Changelog**: https://github.com/shandley/biometal/blob/main/CHANGELOG.md
