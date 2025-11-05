# biometal Python Bindings

Python interface for biometal - ARM-native bioinformatics library with streaming architecture and NEON acceleration.

## Features

- **Streaming architecture**: Constant ~5 MB memory (analyze TB-scale files)
- **ARM NEON acceleration**: 16-25× speedup on ARM processors
- **Zero-copy operations**: Efficient data processing
- **Pythonic API**: Familiar iterator patterns
- **Type hints**: Full IDE support (coming soon)

## Installation

### From Source

```bash
# Install maturin
pip install maturin

# Build and install in development mode
maturin develop --release

# Or build a wheel
maturin build --release
pip install target/wheels/biometal_rs-*.whl
```

### Requirements

- Python 3.8+
- Rust 1.70+ (for building from source)
- ARM processor for NEON acceleration (optional, falls back to scalar)

## Quick Start

```python
import biometal

# Stream FASTQ file
stream = biometal.FastqStream.from_path("data.fq.gz")

for record in stream:
    # Convert to bytes for operations
    seq_bytes = bytes(record.sequence)
    qual_bytes = bytes(record.quality)

    # NEON-accelerated operations
    gc = biometal.gc_content(seq_bytes)
    mean_q = biometal.mean_quality(qual_bytes)
    counts = biometal.count_bases(seq_bytes)

    print(f"{record.id}: GC={gc:.2%}, Q={mean_q:.1f}")
```

## API Reference

### Streaming Classes

#### `FastqStream`

Stream FASTQ records with constant memory.

```python
stream = biometal.FastqStream.from_path("data.fq.gz")

for record in stream:
    print(record.id)              # Read identifier
    print(record.sequence)        # Sequence as list[int]
    print(record.sequence_str)    # Sequence as string
    print(record.quality)         # Quality as list[int]
    print(record.quality_str)     # Quality as string
```

**Supported formats**: `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz`

#### `FastaStream`

Stream FASTA records with constant memory.

```python
stream = biometal.FastaStream.from_path("genome.fa.gz")

for record in stream:
    print(record.id)              # Sequence identifier
    print(record.sequence)        # Sequence as list[int]
    print(record.sequence_str)    # Sequence as string
```

**Supported formats**: `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`

### NEON-Accelerated Operations

All operations use ARM NEON SIMD for 16-25× speedup. Automatically falls back to scalar code on x86_64.

#### `gc_content(sequence: bytes) -> float`

Calculate GC content (20.3× speedup on ARM).

```python
gc = biometal.gc_content(b"ATGCATGC")
print(f"GC content: {gc:.2%}")  # GC content: 50.00%
```

#### `count_bases(sequence: bytes) -> dict`

Count bases (16.7× speedup on ARM).

```python
counts = biometal.count_bases(b"ATGCATGC")
print(counts)  # {'A': 2, 'C': 2, 'G': 2, 'T': 2}
```

#### `mean_quality(quality: bytes) -> float`

Calculate mean Phred quality score (25.1× speedup on ARM).

```python
mean_q = biometal.mean_quality(b"IIIIIIII")
print(f"Mean Q: {mean_q:.1f}")  # Mean Q: 40.0
```

### K-mer Extraction

#### `extract_kmers(sequence: bytes, k: int) -> list[str]`

Extract overlapping k-mers.

```python
kmers = biometal.extract_kmers(b"ATGCATGC", 3)
print(kmers)  # ['ATG', 'TGC', 'GCA', 'CAT', 'ATG', 'TGC']
```

#### `extract_kmers_non_overlapping(sequence: bytes, k: int) -> list[str]`

Extract non-overlapping k-mers.

```python
kmers = biometal.extract_kmers_non_overlapping(b"ATGCATGC", 4)
print(kmers)  # ['ATGC', 'ATGC']
```

## Examples

### Quality Filtering

```python
stream = biometal.FastqStream.from_path("reads.fq.gz")

high_quality = []
for record in stream:
    qual_bytes = bytes(record.quality)
    mean_q = biometal.mean_quality(qual_bytes)

    if mean_q >= 30.0:
        high_quality.append(record)

print(f"High quality reads: {len(high_quality)}")
```

### GC Content Distribution

```python
stream = biometal.FastqStream.from_path("reads.fq.gz")

gc_values = []
for record in stream:
    seq_bytes = bytes(record.sequence)
    gc = biometal.gc_content(seq_bytes)
    gc_values.append(gc)

import statistics
print(f"Mean GC: {statistics.mean(gc_values):.2%}")
print(f"Stdev: {statistics.stdev(gc_values):.2%}")
```

### K-mer Counting for ML

```python
from collections import Counter

stream = biometal.FastqStream.from_path("reads.fq.gz")

all_kmers = []
for record in stream:
    seq_bytes = bytes(record.sequence)
    kmers = biometal.extract_kmers(seq_bytes, 3)
    all_kmers.extend(kmers)

# Get k-mer frequencies
kmer_counts = Counter(all_kmers)
print(f"Most common k-mers:")
for kmer, count in kmer_counts.most_common(10):
    print(f"  {kmer}: {count}")
```

### Integration with Pandas

```python
import pandas as pd
import biometal

stream = biometal.FastqStream.from_path("reads.fq.gz")

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

# Plot
import matplotlib.pyplot as plt
df['gc_content'].hist(bins=20)
plt.xlabel('GC Content')
plt.ylabel('Frequency')
plt.title('GC Content Distribution')
plt.show()
```

## Performance

### ARM NEON Speedups

Measured on Apple M3 Max vs scalar baseline (N=30):

| Operation | Speedup | Throughput (ARM) | Throughput (x86_64) |
|-----------|---------|------------------|---------------------|
| GC content | 20.3× | ~5,954 Kseq/s | ~294 Kseq/s |
| Base counting | 16.7× | ~5,254 Kseq/s | ~315 Kseq/s |
| Quality filtering | 25.1× | ~6,143 Kseq/s | ~245 Kseq/s |

### Memory Usage

- **Streaming**: Constant ~5 MB (regardless of file size)
- **Parallel decompression**: ~1 MB bounded
- **Total footprint**: ~6 MB for TB-scale files

## Examples

See the `examples/` directory:

- `basic_usage.py`: Command-line examples
- `biometal_demo.ipynb`: Jupyter notebook with visualizations

## Testing

Run the test suite:

```bash
python test_python_bindings.py
```

Expected output:
```
============================================================
biometal Python Bindings Test Suite
============================================================
✓ Version test passed
✓ FASTQ streaming test passed
✓ FASTA streaming test passed
✓ GC content test passed
✓ Base counting test passed
✓ Mean quality test passed
✓ K-mer extraction test passed
✓ Integration test passed
============================================================
✓ ALL TESTS PASSED
============================================================
```

## Advanced Usage

### Batch Processing

Process multiple files efficiently:

```python
import glob
import biometal

files = glob.glob("data/*.fq.gz")

for file in files:
    stream = biometal.FastqStream.from_path(file)

    total_reads = 0
    total_bases = 0

    for record in stream:
        total_reads += 1
        total_bases += len(record.sequence)

    print(f"{file}: {total_reads:,} reads, {total_bases:,} bases")
```

### Filtering and Writing

Filter reads and write to new file:

```python
import gzip
import biometal

# Read from stream
stream = biometal.FastqStream.from_path("input.fq.gz")

# Write filtered reads
with gzip.open("filtered.fq.gz", "wt") as out:
    for record in stream:
        seq_bytes = bytes(record.sequence)
        qual_bytes = bytes(record.quality)

        # Filter by quality and GC content
        mean_q = biometal.mean_quality(qual_bytes)
        gc = biometal.gc_content(seq_bytes)

        if mean_q >= 30.0 and 0.4 <= gc <= 0.6:
            # Write FASTQ record
            out.write(f"@{record.id}\\n")
            out.write(f"{record.sequence_str}\\n")
            out.write("+\\n")
            out.write(f"{record.quality_str}\\n")
```

## Troubleshooting

### "list object cannot be converted to PyBytes"

Record sequences are returned as `list[int]` (Vec<u8> in Rust). Convert to bytes for operations:

```python
# ✗ Wrong
gc = biometal.gc_content(record.sequence)

# ✓ Correct
seq_bytes = bytes(record.sequence)
gc = biometal.gc_content(seq_bytes)
```

### "Invalid gzip magic" error

biometal expects gzip-compressed files. If you have uncompressed files, compress them first:

```bash
gzip data.fq
```

Or use the compressed format directly from upstream tools.

### Building on x86_64

biometal builds on x86_64 but uses scalar fallback (no NEON). Performance will be ~16-25× slower than ARM, but still functional.

## Contributing

See `CLAUDE.md` for development guidelines and evidence-based optimization rules.

## License

MIT OR Apache-2.0 (dual-licensed)

## Citation

If you use biometal in your research, please cite:

```bibtex
@software{biometal2025,
  author = {Handley, Scott},
  title = {biometal: ARM-native bioinformatics library},
  year = {2025},
  url = {https://github.com/shandley/biometal}
}
```

## Links

- Repository: https://github.com/shandley/biometal
- Documentation: https://docs.rs/biometal
- Issues: https://github.com/shandley/biometal/issues
