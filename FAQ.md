# Frequently Asked Questions (FAQ)

Common questions and troubleshooting for biometal.

---

## Installation

### Q: How do I install biometal for Python?

**A:** Two options:

```bash
# Option 1: From PyPI (easiest, when published)
pip install biometal-rs

# Option 2: From source (requires Rust toolchain)
pip install maturin
git clone https://github.com/shandley/biometal
cd biometal
maturin develop --release --features python
```

**Requirements**: Python 3.9+

---

### Q: Do I need Rust to use biometal with Python?

**A:** Not once we publish to PyPI (coming soon). For now, yes - you need the Rust toolchain to build from source.

**Install Rust**:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

After PyPI publication, you'll be able to `pip install biometal-rs` without Rust.

---

### Q: Installation fails with "cargo not found"

**A:** You need the Rust toolchain. Install it:

```bash
# macOS/Linux
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Then reload your shell
source $HOME/.cargo/env

# Verify
cargo --version
```

Then retry: `maturin develop --release --features python`

---

### Q: Build fails with "linker 'cc' not found"

**A:** You need a C compiler (for system dependencies).

**macOS**:
```bash
xcode-select --install
```

**Ubuntu/Debian**:
```bash
sudo apt update
sudo apt install build-essential
```

**CentOS/RHEL**:
```bash
sudo yum groupinstall "Development Tools"
```

---

### Q: How do I know if biometal is using ARM NEON?

**A:** Check at runtime:

```python
import biometal

# If on ARM, NEON is automatic
# No configuration needed

# To verify speedup, run benchmarks:
# cargo bench --features python
```

**ARM platforms** (NEON enabled):
- Mac M1/M2/M3/M4
- AWS Graviton
- Ampere Altra
- Raspberry Pi 4+ (64-bit mode)

**x86_64 platforms** (scalar fallback):
- Intel processors
- AMD processors

---

## Common Errors

### Q: "list object cannot be converted to PyBytes"

**A:** Record sequences are returned as `list[int]` (Python representation of Rust `Vec<u8>`). Convert to `bytes` for operations:

```python
# ✗ Wrong
gc = biometal.gc_content(record.sequence)

# ✓ Correct
seq_bytes = bytes(record.sequence)
gc = biometal.gc_content(seq_bytes)
```

**Why?** PyO3 bindings expose Rust types directly. Use `bytes()` to convert.

---

### Q: "No such file or directory" when opening file

**A:** Check the path and file exists:

```python
import os

path = "data.fq.gz"

# Check if file exists
if not os.path.exists(path):
    print(f"File not found: {path}")
    print(f"Current directory: {os.getcwd()}")
else:
    stream = biometal.FastqStream.from_path(path)
```

**Common issue**: Relative paths when script is in different directory.

**Solution**: Use absolute paths:
```python
import os
path = os.path.abspath("data.fq.gz")
stream = biometal.FastqStream.from_path(path)
```

---

### Q: "Invalid gzip magic" error

**A:** biometal expects gzip-compressed files (`.gz`). If you have uncompressed FASTQ:

```bash
# Compress the file
gzip data.fq

# Now use data.fq.gz
```

Or biometal can read uncompressed files too:
```python
stream = biometal.FastqStream.from_path("data.fq")  # No .gz
```

---

### Q: "Stream ended unexpectedly" or "Truncated gzip"

**A:** The file may be corrupted or incomplete. Verify:

```bash
# Check file integrity
gunzip -t data.fq.gz

# If corrupted, re-download or re-generate
```

---

## Performance Questions

### Q: How much faster is biometal on ARM?

**A:** Measured on Mac M3 Max (N=30):

| Operation | Speedup vs Scalar |
|-----------|-------------------|
| GC content | 20.3× |
| Base counting | 16.7× |
| Quality filtering | 25.1× |

**Evidence**: OPTIMIZATION_RULES.md (1,357 experiments)

---

### Q: What's the performance on x86_64?

**A:** Scalar fallback (1× baseline). biometal is **optimized for ARM**, but works correctly on x86_64.

**Cross-platform validated** (Nov 2025):
- **Mac ARM**: 16-25× speedup (optimized)
- **AWS Graviton**: 6-10× speedup (portable)
- **x86_64**: 1× baseline (portable)

**Strategy**: Consumer ARM democratization (not x86 optimization).

---

### Q: How much memory does biometal use?

**A:** Constant ~5 MB regardless of file size.

**Comparison**:
| Dataset Size | Traditional | biometal | Reduction |
|--------------|-------------|----------|-----------|
| 1 GB | 1 GB | 5 MB | 99.5% |
| 100 GB | 100 GB | 5 MB | 99.995% |
| 5 TB | Impossible | 5 MB | ∞ |

**Evidence**: Rule 5 (Entry 026), streaming architecture

---

### Q: Can I process files larger than RAM?

**A:** Yes! That's the whole point of streaming architecture.

```python
# This works on a 16 GB laptop
stream = biometal.FastqStream.from_path("5TB_dataset.fq.gz")
for record in stream:
    # Process one record at a time
    # Memory stays at ~5 MB
    pass
```

**Traditional tools**: Load entire file → Out of memory
**biometal**: Stream records → Constant memory

---

### Q: Is network streaming slower than local files?

**A:** Not significantly, thanks to smart caching and prefetching.

**Features**:
- LRU cache (50 MB, configurable)
- Background prefetching (hides latency)
- Request deduplication

**Evidence**: Rule 6 (Entry 028), I/O dominates 264-352×

---

## Usage Questions

### Q: Can I use biometal with pandas?

**A:** Yes! Example:

```python
import biometal
import pandas as pd

stream = biometal.FastqStream.from_path("data.fq.gz")

data = []
for record in stream:
    seq_bytes = bytes(record.sequence)
    qual_bytes = bytes(record.quality)

    data.append({
        'read_id': record.id,
        'length': len(record.sequence),
        'gc_content': biometal.gc_content(seq_bytes),
        'mean_quality': biometal.mean_quality(qual_bytes),
    })

df = pd.DataFrame(data)
print(df.describe())
```

See [PYTHON.md](PYTHON.md) for more integration examples.

---

### Q: How do I filter reads by quality?

**A:**

```python
import biometal

stream = biometal.FastqStream.from_path("reads.fq.gz")

high_quality = []
for record in stream:
    qual_bytes = bytes(record.quality)
    mean_q = biometal.mean_quality(qual_bytes)

    if mean_q >= 30.0:
        high_quality.append(record)

print(f"High quality reads: {len(high_quality)}")
```

---

### Q: Can I write filtered reads to a new file?

**A:** Yes, but you need to handle writing manually (biometal is read-only for now):

```python
import biometal
import gzip

stream = biometal.FastqStream.from_path("input.fq.gz")

with gzip.open("output.fq.gz", "wt") as out:
    for record in stream:
        qual_bytes = bytes(record.quality)
        mean_q = biometal.mean_quality(qual_bytes)

        if mean_q >= 30.0:
            # Write FASTQ format
            out.write(f"@{record.id}\n")
            out.write(f"{record.sequence_str}\n")
            out.write("+\n")
            out.write(f"{record.quality_str}\n")
```

---

### Q: How do I stream from SRA without downloading?

**A:**

```python
import biometal

# Stream directly from NCBI SRA (S3 bucket)
url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390728/SRR390728"
stream = biometal.FastqStream.from_path(url)

for record in stream:
    # Process without downloading
    # Memory: ~5 MB
    pass
```

See [examples/sra_streaming.rs](examples/sra_streaming.rs) for more details.

---

### Q: Can I extract k-mers for machine learning?

**A:** Yes, biometal includes k-mer utilities:

```python
import biometal

seq = b"ATGCATGC"

# Overlapping k-mers
kmers = biometal.extract_kmers(seq, k=3)
# Returns: ['ATG', 'TGC', 'GCA', 'CAT', 'ATG', 'TGC']

# Non-overlapping k-mers
kmers = biometal.extract_kmers_non_overlapping(seq, k=4)
# Returns: ['ATGC', 'ATGC']
```

Perfect for BERT preprocessing pipelines.

---

## Platform-Specific

### Q: Does biometal work on Raspberry Pi?

**A:** Yes, on 64-bit ARM mode with NEON support.

**Requirements**:
- Raspberry Pi 4 or newer
- 64-bit OS (e.g., Raspberry Pi OS 64-bit)
- Rust 1.70+

**Performance**: Slower than Mac M-series, but still benefits from NEON.

---

### Q: Does biometal work on AWS Graviton?

**A:** Yes! Cross-platform validated on Graviton 3 (Nov 2025).

**Performance**: 6-10× NEON speedup (portable, not optimized)
**Tests**: 121/121 passing

See [docs/CROSS_PLATFORM_TESTING.md](docs/CROSS_PLATFORM_TESTING.md) for details.

---

### Q: Does biometal work on Intel Macs (x86_64)?

**A:** Yes, with scalar fallback (no NEON). Performance will be baseline (1×).

**Tests**: 118/118 passing (3 ARM-specific tests skipped)

**Recommendation**: Upgrade to Apple Silicon for 16-25× speedup.

---

### Q: Can I use biometal on Windows?

**A:** Not officially tested yet. Cross-compilation should work, but ARM NEON benefits require ARM Windows (Surface X, etc.).

**Status**: Community contributions welcome for Windows testing.

---

## Advanced Usage

### Q: How do I tune performance for my workload?

**A:** See [docs/PERFORMANCE_TUNING.md](docs/PERFORMANCE_TUNING.md) for detailed guide.

**Quick tips**:
- Use `.gz` files (parallel bgzip decompression)
- Files ≥50 MB benefit from smart mmap (macOS)
- Network streaming: Tune prefetch count and cache size

---

### Q: What's the overhead of streaming vs batch?

**A:** Negligible. Block-based processing (10K records) preserves NEON speedup.

**Evidence**: Rule 2 (Entry 027), 82-86% overhead avoided with blocks

---

### Q: Can I contribute to biometal?

**A:** Yes! See [CLAUDE.md](CLAUDE.md) for development guidelines.

**Key principles**:
1. Evidence-based optimization (validate with experiments)
2. Streaming-first architecture (constant memory)
3. ARM-native with portable fallback
4. Production quality (tests, docs, error handling)

---

## Troubleshooting

### Q: My code is slow. What's wrong?

**Check**:
1. **ARM platform?** x86_64 uses scalar fallback (1× baseline)
2. **Using bytes()?** Operations require `bytes`, not `list[int]`
3. **File compressed?** `.gz` files enable parallel decompression
4. **Memory pressure?** Check system resources (Activity Monitor/htop)

**Profile**:
```python
import time

start = time.time()
# Your code here
elapsed = time.time() - start

print(f"Elapsed: {elapsed:.2f}s")
```

---

### Q: Tests are failing. What should I check?

**A:**

```bash
# Run full test suite
cargo test --all-features

# Run Python tests
maturin develop --release --features python
python test_python_bindings.py

# Check platform
uname -m  # Should show aarch64 (ARM) or x86_64
```

**If tests fail**:
1. Check Rust version: `rustc --version` (need 1.70+)
2. Check Python version: `python --version` (need 3.9+)
3. Report issue: https://github.com/shandley/biometal/issues

---

### Q: How do I report a bug?

**A:**

1. **Check FAQ** (this file) first
2. **Check existing issues**: https://github.com/shandley/biometal/issues
3. **Create new issue** with:
   - biometal version (`biometal.__version__`)
   - Platform (`uname -m`, `uname -s`)
   - Python version (`python --version`)
   - Minimal reproducible example
   - Error message (full traceback)

---

## Integration

### Q: Can I use biometal with BioPython?

**A:** Yes, but biometal is designed to replace BioPython's parsing (faster, less memory).

**Comparison**:
| Feature | BioPython | biometal |
|---------|-----------|----------|
| Parsing | `SeqIO.parse()` | `FastqStream.from_path()` |
| Memory | Loads all | Constant ~5 MB |
| Speed (ARM) | 1× | 16-25× |
| Streaming | No | Yes |

**Use case**: Use biometal for parsing, BioPython for downstream analysis.

---

### Q: Can I use biometal with samtools/bwa?

**A:** biometal reads FASTQ/FASTA. For BAM/SAM, use samtools directly.

**Pipeline example**:
```bash
# Option 1: biometal preprocessing → samtools alignment
python preprocess.py input.fq.gz > filtered.fq
samtools view -bS filtered.fq > aligned.bam

# Option 2: Use biometal for QC, samtools for alignment
python qc_report.py input.fq.gz  # biometal
bwa mem ref.fa input.fq.gz > aligned.sam  # existing tools
```

**Future**: BAM/SAM support planned (community-driven).

---

## Getting Help

### Q: Where can I get help?

**A:**

1. **Check documentation**:
   - [QUICKSTART.md](QUICKSTART.md) - 5-minute start
   - [README.md](README.md) - Complete overview
   - [PYTHON.md](PYTHON.md) - Python API reference
   - [FAQ.md](FAQ.md) - This file

2. **Search issues**: https://github.com/shandley/biometal/issues

3. **Ask questions**: https://github.com/shandley/biometal/discussions

4. **Report bugs**: https://github.com/shandley/biometal/issues/new

---

### Q: Is there a Slack/Discord for biometal?

**A:** Not yet. Use GitHub Discussions for now:
https://github.com/shandley/biometal/discussions

---

### Q: How do I cite biometal?

**A:**

```bibtex
@software{biometal2025,
  author = {Handley, Scott},
  title = {biometal: ARM-native bioinformatics with streaming architecture},
  year = {2025},
  version = {1.0.0},
  url = {https://github.com/shandley/biometal}
}
```

See [README.md](README.md) for full citation.

---

## Didn't Find Your Answer?

**Open an issue**: https://github.com/shandley/biometal/issues/new

**Start a discussion**: https://github.com/shandley/biometal/discussions/new

**We're here to help!**

---

**Status**: v1.0.0 - Production Release (November 5, 2025)
**Tests**: 121 passing
**Grade**: A+ (rust-code-quality-reviewer)
**Platforms**: Mac ARM (optimized), Graviton/x86_64 (portable)
