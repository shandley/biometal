<div align="center">
  <img src="biometal_logo.png" alt="biometal logo" width="200"/>
</div>

<h1 align="center">biometal</h1>

<p align="center">
<strong>ARM-native bioinformatics library with streaming architecture and evidence-based optimization</strong>
</p>

<p align="center">
<a href="https://crates.io/crates/biometal"><img src="https://img.shields.io/crates/v/biometal.svg" alt="Crates.io"></a>
<a href="https://docs.rs/biometal"><img src="https://docs.rs/biometal/badge.svg" alt="Documentation"></a>
<a href="https://pypi.org/project/biometal-rs/"><img src="https://img.shields.io/pypi/v/biometal-rs.svg" alt="PyPI"></a>
<a href="https://pypi.org/project/biometal-rs/"><img src="https://img.shields.io/pypi/pyversions/biometal-rs.svg" alt="Python"></a>
<a href="https://deepwiki.com/shandley/biometal"><img src="https://deepwiki.com/badge.svg" alt="Ask DeepWiki"></a>
<a href="https://github.com/shandley/biometal#license"><img src="https://img.shields.io/crates/l/biometal.svg" alt="License"></a>
</p>

---

## What Makes biometal Different?

Stream data directly from networks and analyze terabyte-scale datasets on consumer hardware without downloading.

- **Constant ~5 MB memory** regardless of dataset size (99.5% reduction)
- **16-25× speedup** using ARM NEON SIMD on Apple Silicon
- **Network streaming** from HTTP/HTTPS sources (no download needed)
- **Evidence-based** design (1,357 experiments, 40,710 measurements)

---

## Quick Start

### Installation

**Rust:**
```toml
[dependencies]
biometal = "1.2"
```

**Python:**
```bash
pip install biometal-rs  # Install
python -c "import biometal; print(biometal.__version__)"  # Test
```

> **Note**: Package is `biometal-rs` on PyPI, but imports as `biometal` in Python.

### Basic Usage

**Rust:**
```rust
use biometal::FastqStream;

// Stream FASTQ with constant memory (~5 MB)
let stream = FastqStream::from_path("dataset.fq.gz")?;

for record in stream {
    let record = record?;
    // Process one record at a time
}
```

**Python:**
```python
import biometal

# Stream FASTQ with constant memory (~5 MB)
stream = biometal.FastqStream.from_path("dataset.fq.gz")

for record in stream:
    # ARM NEON accelerated (16-25× speedup)
    gc = biometal.gc_content(record.sequence)
    counts = biometal.count_bases(record.sequence)
    mean_q = biometal.mean_quality(record.quality)
```

---

## 📚 Documentation

- **📖 [Comprehensive Docs](https://deepwiki.com/shandley/biometal)** - DeepWiki AI-assisted documentation
- **📓 [Interactive Tutorials](notebooks/)** - Jupyter notebooks with real workflows
- **🦀 [API Reference](https://docs.rs/biometal)** - Full Rust documentation
- **🐍 [Python Guide](docs/PYTHON.md)** - Python-specific documentation
- **📐 [Architecture](docs/ARCHITECTURE.md)** - Technical design details
- **📊 [Benchmarks](docs/BENCHMARKS.md)** - Performance analysis
- **❓ [FAQ](FAQ.md)** - Frequently asked questions

---

## 📓 Interactive Tutorials

Learn biometal through hands-on Jupyter notebooks (4 complete, ~2 hours):

| Notebook | Duration | Topics |
|----------|----------|--------|
| [01. Getting Started](notebooks/01_getting_started.ipynb) | 15-20 min | Streaming, GC content, quality analysis |
| [02. Quality Control](notebooks/02_quality_control_pipeline.ipynb) | 30-40 min | Trimming, filtering, masking (v1.2.0) |
| [03. K-mer Analysis](notebooks/03_kmer_analysis.ipynb) | 30-40 min | ML preprocessing, DNABert (v1.1.0) |
| [04. Network Streaming](notebooks/04_network_streaming.ipynb) | 30-40 min | HTTP streaming, public data (v1.0.0) |

👉 **[Browse all tutorials →](notebooks/README.md)**

---

## 🚀 Key Features

### Streaming Architecture
- **Constant ~5 MB memory** regardless of dataset size
- Analyze 5TB datasets on laptops without downloading
- 99.5% memory reduction vs. traditional approaches

### ARM-Native Performance
- **16-25× speedup** using ARM NEON SIMD
- Optimized for Apple Silicon (M1/M2/M3/M4)
- Automatic scalar fallback on x86_64

### Network Streaming
- Stream directly from HTTP/HTTPS (no download)
- Smart LRU caching + background prefetching
- Access public data (ENA, S3, GCS, Azure)

### Operations Library
- **Core operations**: GC content, base counting, quality scores
- **K-mer operations**: Extraction, minimizers, spectrum (v1.1.0)
- **QC operations**: Trimming, filtering, masking (v1.2.0)
- **40+ Python functions** for bioinformatics workflows

---

## Performance Highlights

| Operation | Scalar | ARM NEON | Speedup |
|-----------|--------|----------|---------|
| Base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7×** |
| GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3×** |
| Quality filter | 245 Kseq/s | 6,143 Kseq/s | **25.1×** |

| Dataset Size | Traditional | biometal | Reduction |
|--------------|-------------|----------|-----------|
| 100K sequences | 134 MB | 5 MB | 96.3% |
| 1M sequences | 1,344 MB | 5 MB | 99.5% |
| **5TB dataset** | **5,000 GB** | **5 MB** | **99.9999%** |

---

## Platform Support

| Platform | Performance | Tests | Status |
|----------|-------------|-------|--------|
| **Mac ARM** (M1-M4) | **16-25× speedup** | ✅ 347/347 | Optimized |
| **AWS Graviton** | 6-10× speedup | ✅ 347/347 | Portable |
| **Linux x86_64** | 1× (scalar) | ✅ 344/344 | Portable |

---

## Evidence-Based Design

biometal's design is grounded in comprehensive experimental validation:

- **1,357 experiments** (40,710 measurements, N=30)
- **Statistical rigor** (95% CI, Cohen's d effect sizes)
- **Full methodology**: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)
- **6 optimization rules** documented in [OPTIMIZATION_RULES.md](OPTIMIZATION_RULES.md)

---

## Roadmap

**v1.0.0** (Released Nov 5, 2025) ✅ - Core library + network streaming
**v1.1.0** (Released Nov 6, 2025) ✅ - K-mer operations
**v1.2.0** (Released Nov 6, 2025) ✅ - Python bindings for Phase 4 QC

**Future** (Community Driven):
- Extended operations (alignment, assembly)
- Additional formats (BAM/SAM, VCF)
- Metal GPU acceleration (Mac-specific)

See [CHANGELOG.md](CHANGELOG.md) for detailed release notes.

---

## Mission: Democratizing Bioinformatics

biometal addresses barriers that lock researchers out of genomics:

1. **Economic**: Consumer ARM laptops ($1,400) deliver production performance
2. **Environmental**: ARM efficiency reduces carbon footprint
3. **Portability**: Works across ARM ecosystem (Mac, Graviton, Ampere, RPi)
4. **Data Access**: Analyze 5TB datasets on 24GB laptops without downloading

---

## Example Use Cases

### Quality Control Pipeline

```python
import biometal

stream = biometal.FastqStream.from_path("raw_reads.fq.gz")

for record in stream:
    # Trim low-quality ends
    trimmed = biometal.trim_quality_window(record, min_quality=20, window_size=4)

    # Length filter
    if biometal.meets_length_requirement(trimmed, min_len=50, max_len=150):
        # Mask remaining low-quality bases
        masked = biometal.mask_low_quality(trimmed, min_quality=20)

        # Check masking rate
        mask_rate = biometal.count_masked_bases(masked) / len(masked.sequence)
        if mask_rate < 0.1:
            # Pass QC - process further
            pass
```

### K-mer Extraction for ML

```python
import biometal

# Extract k-mers for DNABert preprocessing
stream = biometal.FastqStream.from_path("dataset.fq.gz")

for record in stream:
    # Extract overlapping k-mers (k=6 typical for DNABert)
    kmers = biometal.extract_kmers(record.sequence, k=6)

    # Format for transformer models
    kmer_string = " ".join(kmer.decode() for kmer in kmers)

    # Feed to DNABert - constant memory!
    model.process(kmer_string)
```

### Network Streaming

```python
import biometal

# Stream from HTTP without downloading
# Works with ENA, S3, GCS, Azure public data
url = "https://example.com/dataset.fq.gz"
stream = biometal.FastqStream.from_path(url)

for record in stream:
    # Analyze directly - no download needed!
    # Memory: constant ~5 MB
    gc = biometal.gc_content(record.sequence)
```

---

## FAQ

**Q: Why `biometal-rs` on PyPI but `biometal` everywhere else?**
A: The `biometal` name was taken on PyPI, so we use `biometal-rs` for installation. You still import as `import biometal`.

**Q: What platforms are supported?**
A: Mac ARM (optimized), Linux ARM/x86_64 (portable). Pre-built wheels for common platforms. See [docs/CROSS_PLATFORM_TESTING.md](docs/CROSS_PLATFORM_TESTING.md).

**Q: Why ARM-native?**
A: To democratize bioinformatics by enabling world-class performance on consumer hardware ($1,400 MacBooks vs. $50,000 servers).

More questions? See [FAQ.md](FAQ.md)

---

## Contributing

We welcome contributions! See [CLAUDE.md](CLAUDE.md) for development guidelines.

biometal is built on evidence-based optimization - new features should:
1. Have clear use cases
2. Be validated experimentally (when adding optimizations)
3. Maintain platform portability
4. Follow [OPTIMIZATION_RULES.md](OPTIMIZATION_RULES.md)

---

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

---

## Citation

If you use biometal in your research:

```bibtex
@software{biometal2025,
  author = {Handley, Scott},
  title = {biometal: ARM-native bioinformatics with streaming architecture},
  year = {2025},
  url = {https://github.com/shandley/biometal}
}
```

For the experimental methodology:
```bibtex
@misc{asbb2025,
  author = {Handley, Scott},
  title = {Apple Silicon Bio Bench: Systematic Hardware Characterization},
  year = {2025},
  url = {https://github.com/shandley/apple-silicon-bio-bench}
}
```

---

<p align="center">
<strong>Status:</strong> v1.2.0 - Python Phase 4 Bindings 🎉<br>
<strong>Released:</strong> November 6, 2025<br>
<strong>Tests:</strong> 347 passing (260 library + 87 doc)<br>
<strong>Python Functions:</strong> 40+ (core ops + k-mers + Phase 4)<br>
<strong>Evidence Base:</strong> 1,357 experiments, 40,710 measurements
</p>
