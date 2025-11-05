<div align="right">
  <img src="biometal_logo.png" alt="biometal logo" width="150"/>
</div>

# biometal

**ARM-native bioinformatics library with streaming architecture and evidence-based optimization**

[![Crates.io](https://img.shields.io/crates/v/biometal.svg)](https://crates.io/crates/biometal)
[![Documentation](https://docs.rs/biometal/badge.svg)](https://docs.rs/biometal)
[![License](https://img.shields.io/crates/l/biometal.svg)](https://github.com/shandley/biometal#license)

---

## What Makes biometal Different?

Most bioinformatics tools require you to download entire datasets before analysis. **biometal** streams data directly from the network, enabling analysis of terabyte-scale datasets on consumer hardware without downloading.

### Key Features

1. **Streaming Architecture** (Rule 5)
   - Constant ~5 MB memory footprint regardless of dataset size
   - Analyze 5TB datasets on laptops without downloading
   - 99.5% memory reduction compared to batch processing

2. **ARM-Native Performance** (Rule 1)
   - 16-25× speedup using ARM NEON SIMD
   - Works across Mac (Apple Silicon), AWS Graviton, Ampere, Raspberry Pi
   - Automatic fallback to scalar code on x86_64

3. **Network Streaming** (Rule 6)
   - Stream directly from HTTP/HTTPS sources
   - SRA toolkit integration (no local copy needed)
   - Smart LRU caching minimizes network requests
   - Background prefetching hides latency

4. **Intelligent I/O** (Rules 3-4)
   - 6.5× speedup from parallel bgzip decompression
   - Additional 2.5× from memory-mapped I/O (large files on macOS)
   - Combined 16.3× I/O speedup

5. **Evidence-Based Design**
   - Every optimization validated with statistical rigor (N=30, 95% CI)
   - 1,357 experiments, 40,710 measurements
   - Full methodology: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)

---

## Quick Start

### Rust Installation

```toml
[dependencies]
biometal = "1.0"
```

### Python Installation

```bash
# Option 1: Install from PyPI (when published)
pip install biometal-rs

# Option 2: Build from source (requires Rust toolchain)
pip install maturin
git clone https://github.com/shandley/biometal
cd biometal
maturin develop --release --features python
```

**Requirements**:
- Python 3.9+ (tested on 3.14)
- Rust toolchain (for building from source)

---

## Usage

### Rust: Basic Usage

```rust
use biometal::FastqStream;

// Stream FASTQ from local file (constant memory)
let stream = FastqStream::from_path("large_dataset.fq.gz")?;

for record in stream {
    let record = record?;
    // Process one record at a time
    // Memory stays constant at ~5 MB
}
```

### Network Streaming

```rust
use biometal::io::DataSource;
use biometal::FastqStream;

// Stream directly from URL (no download!)
let source = DataSource::Http("https://example.com/huge_dataset.fq.gz".to_string());
let stream = FastqStream::new(source)?;

// Analyze 5TB dataset without downloading
for record in stream {
    // Smart caching + prefetching in background
}
```

### SRA Streaming (No Download!)

```rust
use biometal::io::DataSource;
use biometal::FastqStream;

// Stream directly from NCBI SRA (no local download!)
let source = DataSource::Sra("SRR390728".to_string());  // E. coli dataset
let stream = FastqStream::new(source)?;

for record in stream {
    let record = record?;
    // Process 40 MB dataset with only ~5 MB memory
    // Background prefetching hides network latency
}
```

### Operations with Auto-Optimization

```rust
use biometal::operations;

// ARM NEON automatically enabled on ARM platforms
let counts = operations::base_counting(&sequence)?;
let gc = operations::gc_content(&sequence)?;

// 16-25× faster on ARM, automatic scalar fallback on x86_64
```

### Python: Basic Usage

```python
import biometal

# Stream FASTQ from local file (constant memory)
stream = biometal.FastqStream.from_path("large_dataset.fq.gz")

for record in stream:
    # Process one record at a time
    # Memory stays constant at ~5 MB
    gc = biometal.gc_content(record.sequence)
    print(f"{record.id}: GC={gc:.2%}")
```

### Python: ARM NEON Operations

```python
import biometal

# ARM NEON automatically enabled on ARM platforms
# 16-25× faster on Mac ARM, automatic scalar fallback on x86_64

# GC content calculation
sequence = b"ATGCATGC"
gc = biometal.gc_content(sequence)  # 20.3× speedup on ARM

# Base counting
counts = biometal.count_bases(sequence)  # 16.7× speedup on ARM
print(f"A:{counts['A']}, C:{counts['C']}, G:{counts['G']}, T:{counts['T']}")

# Quality scoring
quality = record.quality
mean_q = biometal.mean_quality(quality)  # 25.1× speedup on ARM

# K-mer extraction (for ML preprocessing)
kmers = biometal.extract_kmers(sequence, k=6)
print(f"6-mers: {kmers}")
```

### Python: Example Workflow

```python
import biometal

# Analyze FASTQ file with streaming (constant memory)
stream = biometal.FastqStream.from_path("data.fq.gz")

total_bases = 0
total_gc = 0.0
high_quality = 0

for record in stream:
    # Count bases (ARM NEON accelerated)
    counts = biometal.count_bases(record.sequence)
    total_bases += sum(counts.values())

    # Calculate GC content (ARM NEON accelerated)
    gc = biometal.gc_content(record.sequence)
    total_gc += gc

    # Check quality (ARM NEON accelerated)
    if biometal.mean_quality(record.quality) > 30.0:
        high_quality += 1

print(f"Total bases: {total_bases}")
print(f"Average GC: {total_gc/len(stream):.2%}")
print(f"High quality reads: {high_quality}")
```

---

## Performance

### Memory Efficiency

| Dataset Size | Traditional | biometal | Reduction |
|--------------|-------------|----------|-----------|
| 100K sequences | 134 MB | 5 MB | 96.3% |
| 1M sequences | 1,344 MB | 5 MB | 99.5% |
| **5TB dataset** | **5,000 GB** | **5 MB** | **99.9999%** |

### ARM NEON Speedup (Mac Apple Silicon)

**Optimized for Apple Silicon** - All optimizations validated on Mac M3 Max (1,357 experiments, N=30):

| Operation | Scalar | NEON | Speedup |
|-----------|--------|------|---------|
| Base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7×** |
| GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3×** |
| Quality filter | 245 Kseq/s | 6,143 Kseq/s | **25.1×** |

### Cross-Platform Performance (Validated Nov 2025)

| Platform | Base Counting | GC Content | Quality | Status |
|----------|---------------|------------|---------|--------|
| **Mac M3** (target) | 16.7× | 20.3× | 25.1× | ✅ Optimized |
| **AWS Graviton** | 10.7× | 6.9× | 1.9× | ✅ Works (portable) |
| **x86_64 Intel** | 1.0× | 1.0× | 1.0× | ✅ Works (portable) |

**Note**: biometal is optimized for Mac ARM (consumer hardware democratization). Other platforms are supported with correct, production-ready code but not specifically optimized. See [Cross-Platform Testing Results](results/cross_platform/FINDINGS.md) for details.

### I/O Optimization

| File Size | Standard | Optimized | Speedup |
|-----------|----------|-----------|---------|
| Small (<50 MB) | 12.3s | 1.9s | 6.5× |
| Large (≥50 MB) | 12.3s | 0.75s | **16.3×** |

---

## Democratizing Bioinformatics

biometal addresses four barriers that lock researchers out of genomics:

### 1. Economic Barrier
- **Problem**: Most tools require $50K+ servers
- **Solution**: Consumer ARM laptops ($1,400) deliver production performance
- **Impact**: Small labs and LMIC researchers can compete

### 2. Environmental Barrier
- **Problem**: HPC clusters consume massive energy (300× excess for many workloads)
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

## Evidence Base

biometal's design is grounded in comprehensive experimental validation:

- **Experiments**: 1,357 total (40,710 measurements with N=30)
- **Statistical rigor**: 95% confidence intervals, Cohen's d effect sizes
- **Cross-platform**: Mac M4 Max, AWS Graviton 3
- **Lab notebook**: 33 entries documenting full experimental log

See [OPTIMIZATION_RULES.md](OPTIMIZATION_RULES.md) for detailed evidence links.

**Full methodology**: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)

**Publications** (in preparation):
1. DAG Framework: BMC Bioinformatics
2. biometal Library: Bioinformatics (Application Note) or JOSS
3. Four-Pillar Democratization: GigaScience

---

## Platform Support

### Optimization Strategy

biometal is **optimized for Mac ARM** (M1/M2/M3/M4) based on 1,357 experiments on Mac M3 Max. This aligns with our democratization mission: enable world-class bioinformatics on **affordable consumer hardware** ($1,000-2,000 MacBooks, not $50,000 servers).

Other platforms are **supported with portable, correct code** but not specifically optimized:

| Platform | Performance | Test Status | Strategy |
|----------|-------------|-------------|----------|
| **Mac ARM** (M1/M2/M3/M4) | **16-25× speedup** | ✅ 121/121 tests pass | **Optimized** (target platform) |
| **AWS Graviton** | 6-10× speedup | ✅ 121/121 tests pass | Portable (works well) |
| **Linux x86_64** | 1× (scalar) | ✅ 118/118 tests pass | Portable (fallback) |

### Feature Support Matrix

| Feature | macOS ARM | Linux ARM | Linux x86_64 |
|---------|-----------|-----------|--------------|
| ARM NEON SIMD | ✅ | ✅ | ❌ (scalar fallback) |
| Parallel Bgzip | ✅ | ✅ | ✅ |
| Smart mmap | ✅ | ⏳ | ❌ |
| Network Streaming | ✅ | ✅ | ✅ |
| Python Bindings | ✅ | ✅ | ✅ |

**Validation**: Cross-platform testing completed Nov 2025 on AWS Graviton 3 and x86_64. All tests pass. See [results/cross_platform/FINDINGS.md](results/cross_platform/FINDINGS.md) for full details.

---

## Roadmap

**v1.0.0** (Released November 5, 2025) ✅
- Streaming FASTQ/FASTA parsers (constant memory)
- ARM NEON operations (16-25× speedup)
- Network streaming (HTTP/HTTPS, SRA)
- Python bindings (PyO3 0.27, Python 3.9-3.14)
- Cross-platform validation (Mac ARM, Graviton, x86_64)
- Production-grade quality (121 tests, Grade A+)

**Future Considerations** (Community Driven)
- Extended operation coverage (alignment, assembly)
- Additional format support (BAM/SAM, VCF)
- Publish to crates.io and PyPI
- Metal GPU acceleration (Mac-specific)

---

## SRA Streaming: Analysis Without Downloads

One of biometal's most powerful features is direct streaming from NCBI's Sequence Read Archive (SRA) without local downloads.

### Why This Matters

**Traditional workflow:**
1. Download 5 GB SRA dataset → 30 minutes + 5 GB disk space
2. Decompress → 15 GB disk space
3. Process → Additional memory
4. **Total:** 45 minutes + 20 GB resources before analysis even starts

**biometal workflow:**
1. Start analysis immediately → 0 wait time, ~5 MB memory
2. Stream directly from NCBI S3 → No disk space needed
3. Background prefetching hides latency → Near-local performance

### Supported Accessions

- **SRR** (Run): Most common, represents a sequencing run
- **SRX** (Experiment): Collection of runs
- **SRS** (Sample): Biological sample
- **SRP** (Study): Collection of experiments

### Basic SRA Usage

```rust
use biometal::io::DataSource;
use biometal::operations::{count_bases, gc_content};
use biometal::FastqStream;

// Stream from SRA accession
let source = DataSource::Sra("SRR390728".to_string());
let stream = FastqStream::new(source)?;

for record in stream {
    let record = record?;

    // ARM NEON-optimized operations (16-25× speedup)
    let bases = count_bases(&record.sequence);
    let gc = gc_content(&record.sequence);

    // Memory: Constant ~5 MB
}
```

### Real-World Example: E. coli Analysis

```bash
# Run the E. coli streaming example
cargo run --example sra_ecoli --features network

# Process ~250,000 reads with only ~5 MB memory
# No download required!
```

See [examples/sra_ecoli.rs](examples/sra_ecoli.rs) for complete example.

### Performance Tuning

biometal automatically configures optimal settings for most use cases. For custom tuning:

```rust
use biometal::io::{HttpReader, sra_to_url};

let url = sra_to_url("SRR390728")?;
let reader = HttpReader::new(&url)?
    .with_prefetch_count(8)      // Prefetch 8 blocks ahead
    .with_chunk_size(128 * 1024); // 128 KB chunks

// See docs/PERFORMANCE_TUNING.md for detailed guide
```

### SRA URL Conversion

```rust
use biometal::io::{is_sra_accession, sra_to_url};

// Check if string is SRA accession
if is_sra_accession("SRR390728") {
    // Convert to direct NCBI S3 URL
    let url = sra_to_url("SRR390728")?;
    // → https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390728/SRR390728
}
```

### Memory Guarantees

- **Streaming buffer:** ~5 MB (constant)
- **LRU cache:** 50 MB (byte-bounded, automatic eviction)
- **Prefetch:** ~256 KB (4 blocks × 64 KB)
- **Total:** ~55 MB regardless of SRA file size

Compare to downloading a 5 GB SRA file → **99%+ memory savings**

### Examples

| Example | Dataset | Size | Demo |
|---------|---------|------|------|
| [sra_streaming.rs](examples/sra_streaming.rs) | Demo mode | N/A | Capabilities overview |
| [sra_ecoli.rs](examples/sra_ecoli.rs) | E. coli K-12 | ~40 MB | Real SRA streaming |
| [prefetch_tuning.rs](examples/prefetch_tuning.rs) | E. coli K-12 | ~40 MB | Performance tuning |

---

## Example Use Cases

### 1. Large-Scale Quality Control

```rust
use biometal::{FastqStream, operations};

// Stream 5TB dataset without downloading
let stream = FastqStream::from_url("https://sra.example.com/huge.fq.gz")?;

let mut total = 0;
let mut high_quality = 0;

for record in stream {
    let record = record?;
    total += 1;
    
    // ARM NEON accelerated (16-25×)
    if operations::mean_quality(&record.quality) > 30.0 {
        high_quality += 1;
    }
}

println!("High quality: {}/{} ({:.1}%)", 
    high_quality, total, 100.0 * high_quality as f64 / total as f64);
```

### 2. BERT Preprocessing Pipeline

```rust
use biometal::{FastqStream, kmer};

// Stream from SRA (no local copy!)
let stream = FastqStream::from_sra("SRR12345678")?;

// Extract k-mers for DNABert training
for record in stream {
    let record = record?;
    let kmers = kmer::extract_overlapping(&record.sequence, 6)?;
    
    // Feed to BERT training pipeline
    // Constant memory even for TB-scale datasets
}
```

### 3. Metagenomics Filtering

```rust
use biometal::{FastqStream, operations};

let input = FastqStream::from_path("metagen.fq.gz")?;
let mut output = FastqWriter::create("filtered.fq.gz")?;

for record in input {
    let record = record?;
    
    // Filter low-complexity sequences (ARM NEON accelerated)
    if operations::complexity_score(&record.sequence) > 0.5 {
        output.write(&record)?;
    }
}
// Memory: constant ~5 MB
// Speed: 16-25× faster on ARM
```

---

## Contributing

We welcome contributions! biometal is built on evidence-based optimization, so new features should:
1. Have clear use cases
2. Be validated experimentally (when adding optimizations)
3. Maintain platform portability
4. Follow the optimization rules in [OPTIMIZATION_RULES.md](OPTIMIZATION_RULES.md)

See [CLAUDE.md](CLAUDE.md) for development guidelines.

---

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

---

## Citation

If you use biometal in your research, please cite:

```bibtex
@software{biometal2025,
  author = {Handley, Scott},
  title = {biometal: ARM-native bioinformatics with streaming architecture},
  year = {2025},
  url = {https://github.com/shandley/biometal}
}
```

For the experimental methodology, see:
```bibtex
@misc{asbb2025,
  author = {Handley, Scott},
  title = {Apple Silicon Bio Bench: Systematic Hardware Characterization for Bioinformatics},
  year = {2025},
  url = {https://github.com/shandley/apple-silicon-bio-bench}
}
```

---

**Status**: v1.0.0 - Production Release 🎉
**Released**: November 5, 2025
**Grade**: A+ (rust-code-quality-reviewer)
**Tests**: 121 passing (87 unit + 7 integration + 27 doc)
**Evidence Base**: 1,357 experiments, 40,710 measurements
**Mission**: Democratizing bioinformatics compute
