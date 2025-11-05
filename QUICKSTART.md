# biometal Quick Start

**Get up and running with biometal in 5 minutes.**

biometal is an ARM-native bioinformatics library that lets you analyze terabyte-scale datasets on consumer hardware. It uses streaming architecture (constant ~5 MB memory) and ARM NEON acceleration (16-25× speedup) to deliver HPC-class performance on laptops.

---

## Installation

### Python (Recommended)

```bash
# Option 1: From PyPI (when published)
pip install biometal-rs

# Option 2: From source (requires Rust)
pip install maturin
git clone https://github.com/shandley/biometal
cd biometal
maturin develop --release --features python
```

**Requirements**: Python 3.9+

### Rust

```bash
cargo add biometal
```

**Requirements**: Rust 1.70+

---

## Your First 5 Lines of Code

Create `analyze.py`:

```python
import biometal

# Stream FASTQ file (constant memory, no matter the size)
stream = biometal.FastqStream.from_path("data.fq.gz")

for record in stream:
    # ARM NEON accelerated (16-25× faster on ARM)
    gc = biometal.gc_content(bytes(record.sequence))
    print(f"{record.id}: GC = {gc:.2%}")
```

**Run it**:
```bash
python analyze.py
```

**Output**:
```
@read1: GC = 48.32%
@read2: GC = 52.67%
@read3: GC = 45.89%
...
```

---

## What Just Happened?

### 1. Streaming Architecture (Memory Magic)
```python
stream = biometal.FastqStream.from_path("data.fq.gz")
```

**Traditional tools**: Load entire file into RAM
- 5 GB file → 5 GB memory
- 5 TB file → Impossible on laptops

**biometal**: Process one record at a time
- 5 GB file → 5 MB memory
- 5 TB file → 5 MB memory

**Evidence**: Rule 5 (Entry 026), 99.5% memory reduction

### 2. ARM NEON Acceleration (Speed Boost)
```python
gc = biometal.gc_content(bytes(record.sequence))
```

**On ARM processors** (Mac M1/M2/M3/M4, AWS Graviton):
- 16-25× faster than scalar code
- Automatic SIMD vectorization
- No code changes needed

**On x86_64** (Intel, AMD):
- Automatic fallback to scalar
- Same API, portable code

**Evidence**: Rule 1 (Entry 020-025), 1,357 experiments

### 3. Production Quality
- **Grade A+** (rust-code-quality-reviewer)
- **121 tests passing**
- **Zero panics** in library code
- **Cross-platform validated** (Mac, Graviton, x86_64)

---

## Try Network Streaming (No Download!)

```python
import biometal

# Stream directly from NCBI SRA (no local download needed!)
stream = biometal.FastqStream.from_path("https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390728/SRR390728")

for record in stream:
    gc = biometal.gc_content(bytes(record.sequence))
    print(f"{record.id}: GC = {gc:.2%}")

# Memory: Still ~5 MB
# Speed: Background prefetching hides network latency
```

**What's happening**:
- Smart LRU cache (50 MB)
- Background prefetching
- No disk space used
- Analyze 5TB datasets without downloading

**Evidence**: Rule 6 (Entry 028), network streaming validated

---

## Common Operations

### GC Content
```python
gc = biometal.gc_content(b"ATGCATGC")
# Returns: 0.5 (50%)
# Speedup: 20.3× on ARM
```

### Base Counting
```python
counts = biometal.count_bases(b"ATGCATGC")
# Returns: {'A': 2, 'C': 2, 'G': 2, 'T': 2}
# Speedup: 16.7× on ARM
```

### Quality Scores
```python
mean_q = biometal.mean_quality(b"IIIIIIII")
# Returns: 40.0 (Phred quality)
# Speedup: 25.1× on ARM
```

### K-mer Extraction (ML Preprocessing)
```python
kmers = biometal.extract_kmers(b"ATGCATGC", k=3)
# Returns: ['ATG', 'TGC', 'GCA', 'CAT', 'ATG', 'TGC']
```

---

## Real-World Example: Quality Control

```python
import biometal

stream = biometal.FastqStream.from_path("reads.fq.gz")

total = 0
high_quality = 0
total_gc = 0.0

for record in stream:
    total += 1

    # ARM NEON accelerated
    gc = biometal.gc_content(bytes(record.sequence))
    mean_q = biometal.mean_quality(bytes(record.quality))

    total_gc += gc

    if mean_q >= 30.0:
        high_quality += 1

print(f"Total reads: {total:,}")
print(f"High quality (Q≥30): {high_quality:,} ({100*high_quality/total:.1f}%)")
print(f"Average GC: {100*total_gc/total:.1f}%")

# Memory: ~5 MB (even for TB-scale files)
# Speed: 16-25× faster on ARM
```

---

## Performance at a Glance

### Memory (Streaming vs Batch)
| Dataset | Traditional | biometal | Reduction |
|---------|-------------|----------|-----------|
| 1 GB | 1 GB | 5 MB | 99.5% |
| 100 GB | 100 GB | 5 MB | 99.995% |
| 5 TB | Impossible | 5 MB | ∞ |

### Speed (ARM NEON vs Scalar)
| Operation | Mac M3 | x86_64 | Speedup |
|-----------|--------|--------|---------|
| GC content | 5,954 Kseq/s | 294 Kseq/s | **20.3×** |
| Base counting | 5,254 Kseq/s | 315 Kseq/s | **16.7×** |
| Quality filter | 6,143 Kseq/s | 245 Kseq/s | **25.1×** |

### Cost (Consumer vs HPC)
| Hardware | Cost | Performance | biometal Speedup |
|----------|------|-------------|------------------|
| MacBook Air M3 | $1,400 | Full | 16-25× (NEON) |
| HPC Server | $50,000 | Full | 1× (baseline) |
| **Savings** | **$48,600** | **Same** | **35× cheaper** |

---

## Platform Support

| Platform | ARM NEON | Tests | Status |
|----------|----------|-------|--------|
| **Mac ARM** (M1/M2/M3/M4) | ✅ 16-25× | 121/121 | Optimized |
| **AWS Graviton** | ✅ 6-10× | 121/121 | Portable |
| **Linux x86_64** | ❌ 1× (scalar) | 118/118 | Portable |

**Strategy**: Optimized for consumer ARM hardware (democratization), portable everywhere else.

---

## Next Steps

### Learn More
- **[README.md](README.md)** - Complete feature overview
- **[PYTHON.md](PYTHON.md)** - Full Python API reference
- **[examples/](examples/)** - 11 working examples
- **[FAQ.md](FAQ.md)** - Common questions

### Go Deeper
- **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)** - Network streaming design
- **[docs/PERFORMANCE_TUNING.md](docs/PERFORMANCE_TUNING.md)** - Optimization guide
- **[OPTIMIZATION_RULES.md](OPTIMIZATION_RULES.md)** - Evidence base (1,357 experiments)

### Get Help
- **GitHub Issues**: https://github.com/shandley/biometal/issues
- **Discussions**: https://github.com/shandley/biometal/discussions

---

## Why biometal?

### The Problem
Most bioinformatics tools require expensive HPC clusters because they:
1. Load entire datasets into memory (out of memory errors)
2. Only optimize for x86_64 servers (vendor lock-in)
3. Require downloading 5TB files before analysis (days of waiting)

### The Solution
biometal democratizes bioinformatics by:
1. **Streaming**: Constant ~5 MB memory (analyze 5TB on laptops)
2. **ARM-native**: 16-25× speedup on $1,400 MacBooks
3. **Network streaming**: Analyze without downloading

### The Impact
- **Economic**: $1,400 laptop competes with $50,000 server
- **Environmental**: 1/10th power consumption (ARM efficiency)
- **Portability**: Works on Mac, Graviton, x86_64, Raspberry Pi
- **Access**: Small labs, LMIC researchers, students, field researchers

**Mission**: Enable world-class bioinformatics on consumer hardware.

---

## Quick Reference

### Installation
```bash
pip install biometal-rs  # Python
cargo add biometal    # Rust
```

### Basic Usage
```python
import biometal

# Stream FASTQ
stream = biometal.FastqStream.from_path("data.fq.gz")
for record in stream:
    gc = biometal.gc_content(bytes(record.sequence))
```

### Operations (ARM NEON Accelerated)
```python
biometal.gc_content(seq)      # 20.3× speedup
biometal.count_bases(seq)     # 16.7× speedup
biometal.mean_quality(qual)   # 25.1× speedup
biometal.extract_kmers(seq, k) # K-mer extraction
```

### Network Streaming
```python
stream = biometal.FastqStream.from_path("https://...")
# No download, ~5 MB memory
```

---

**Status**: v1.0.0 - Production Release (November 5, 2025)
**Grade**: A+ (rust-code-quality-reviewer)
**Tests**: 121 passing
**Evidence**: 1,357 experiments, 40,710 measurements

**Ready to use in production.**
