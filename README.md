<div align="right">
  <img src="biometal_logo.png" alt="biometal logo" width="150"/>
</div>

# biometal

**ARM-native bioinformatics library with streaming architecture and evidence-based optimization**

[![Crates.io](https://img.shields.io/crates/v/biometal.svg)](https://crates.io/crates/biometal)
[![Documentation](https://docs.rs/biometal/badge.svg)](https://docs.rs/biometal)
[![PyPI](https://img.shields.io/pypi/v/biometal-rs.svg)](https://pypi.org/project/biometal-rs/)
[![Python](https://img.shields.io/pypi/pyversions/biometal-rs.svg)](https://pypi.org/project/biometal-rs/)
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
biometal = "1.2"
```

### Python Installation

```bash
# Install from PyPI
pip install biometal-rs

# Then import as 'biometal'
python -c "import biometal; print(biometal.__version__)"
```

> **Note**: The package name is `biometal-rs` on PyPI (the `biometal` name was already taken), but you import it as `biometal` in your Python code. See [FAQ](#faq) for details.

**Alternative - Build from source**:
```bash
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

### K-mer Operations (Evidence-Based)

biometal provides k-mer operations optimized based on ASBB Entry 034 findings.

**Key finding**: K-mer operations are **data-structure-bound** (hash+HashMap), not compute-bound. Unlike element-wise operations (base counting, GC content), k-mers spend 50-60% of runtime on hash computation and 30-40% on data structure operations. Therefore, NEON/GPU provide no benefit.

#### Rust: K-mer Operations

```rust
use biometal::operations::kmer::{extract_kmers, extract_minimizers, kmer_spectrum, KmerExtractor};

// 1. Simple k-mer extraction (scalar-only, optimal)
let sequence = b"ATGCATGCATGC";
let kmers = extract_kmers(sequence, 6);  // Returns Vec<Vec<u8>>

// 2. Minimizers (minimap2-style, scalar-only)
let minimizers = extract_minimizers(sequence, 6, 10);  // k=6, w=10
for minimizer in minimizers {
    println!("Position {}: {:?}", minimizer.position, minimizer.kmer);
}

// 3. K-mer spectrum (frequency counting, scalar-only)
let sequences = vec![b"ATGCAT".as_ref(), b"GCATGC".as_ref()];
let spectrum = kmer_spectrum(&sequences, 3);  // HashMap<Vec<u8>, usize>

// 4. Parallel extraction (opt-in for large datasets, 2.2× speedup)
let extractor = KmerExtractor::with_parallel(4);  // 4 threads (optimal per Entry 034)
let large_dataset: Vec<&[u8]> = /* 10K+ sequences */;
let kmers = extractor.extract(&large_dataset, 6);  // 2.2× faster
```

#### Python: K-mer Operations

```python
import biometal

# 1. Simple k-mer extraction (scalar-only, optimal)
sequence = b"ATGCATGCATGC"
kmers = biometal.extract_kmers(sequence, k=6)  # Returns list[bytes]
print(f"Extracted {len(kmers)} k-mers")

# 2. Minimizers (minimap2-style, scalar-only)
minimizers = biometal.extract_minimizers(sequence, k=6, w=10)
for m in minimizers:
    print(f"Position {m['position']}: {m['kmer']}")

# 3. K-mer spectrum (frequency counting, scalar-only)
sequences = [b"ATGCAT", b"GCATGC"]
spectrum = biometal.kmer_spectrum(sequences, k=3)  # Returns dict
print(f"Unique k-mers: {len(spectrum)}")

# 4. Parallel extraction (opt-in for large datasets, 2.2× speedup)
extractor = biometal.KmerExtractor(parallel=True, threads=4)
large_dataset = [...]  # 10K+ sequences
kmers = extractor.extract(large_dataset, k=6)  # 2.2× faster
```

**Evidence (Entry 034)**:
- **Minimizers**: 1.02-1.26× (NEON/Parallel) → Scalar-only
- **K-mer Spectrum**: 0.95-1.88× (sometimes SLOWER with parallel!) → Scalar-only
- **K-mer Extraction**: 2.19-2.38× (Parallel-4t) → Opt-in parallel

This validates minimap2's scalar design and identifies a 2.2× optimization opportunity for DNABert preprocessing.

### Sequence Manipulation Operations (Phase 4)

biometal provides comprehensive sequence manipulation primitives for read processing pipelines. All operations maintain production quality with proper error handling.

#### Python: Sequence Operations

```python
import biometal

# 1. Reverse complement (standard molecular biology operation)
sequence = b"ATGCATGC"
rc = biometal.reverse_complement(sequence)  # b"GCATGCAT"

# 2. Complement only (preserves 5'→3' orientation)
comp = biometal.complement(sequence)  # b"TACGTACG"

# 3. Reverse only (no complementation)
rev = biometal.reverse(sequence)  # b"CGTAGCTA"

# 4. Sequence validation
if biometal.is_valid_dna(sequence):
    print("Valid DNA sequence")

if biometal.is_valid_rna(b"AUGCAUGC"):
    print("Valid RNA sequence")

# 5. Count invalid bases (for QC)
invalid_count = biometal.count_invalid_bases(sequence)
print(f"Invalid bases: {invalid_count}")
```

#### Python: Record Operations

```python
import biometal

stream = biometal.FastqStream.from_path("data.fq.gz")

for record in stream:
    # 1. Extract region [start, end)
    region = biometal.extract_region(record, start=10, end=50)
    print(f"Extracted 40bp region: {len(region.sequence)}bp")

    # 2. Reverse complement record (preserves quality alignment)
    rc_record = biometal.reverse_complement_record(record)
    # Both sequence AND quality are reversed

    # 3. Get sequence length
    length = biometal.sequence_length(record)

    # 4. Length filtering
    if biometal.meets_length_requirement(record, min_len=50, max_len=150):
        print("Read passes length filter")

    # 5. Convert FASTQ → FASTA (drops quality scores)
    fasta = biometal.to_fasta_record(record)
    print(f">{fasta.id}")
    print(fasta.sequence_str)

    break
```

#### Python: Quality-Based Trimming

```python
import biometal

stream = biometal.FastqStream.from_path("data.fq.gz")

for record in stream:
    # 1. Fixed position trimming
    trimmed = biometal.trim_start(record, bases=10)  # Remove first 10bp
    trimmed = biometal.trim_end(record, bases=5)     # Remove last 5bp
    trimmed = biometal.trim_both(record, start_bases=10, end_bases=5)

    # 2. Quality-based trimming (Phred+33, Q20 = 99% accuracy)
    trimmed = biometal.trim_quality_end(record, min_quality=20)
    trimmed = biometal.trim_quality_start(record, min_quality=20)
    trimmed = biometal.trim_quality_both(record, min_quality=20)

    # 3. Sliding window trimming (Trimmomatic-style)
    # Trim when 4bp window average drops below Q20
    trimmed = biometal.trim_quality_window(record, min_quality=20, window_size=4)

    # 4. QC pipeline: trim + length filter
    trimmed = biometal.trim_quality_both(record, min_quality=20)
    if len(trimmed.sequence) >= 50:  # Keep only ≥50bp after trimming
        print(f"Pass QC: {len(trimmed.sequence)}bp after trimming")

    break
```

#### Python: Quality-Based Masking

```python
import biometal

stream = biometal.FastqStream.from_path("data.fq.gz")

for record in stream:
    # Mask low-quality bases with 'N' (preserves length unlike trimming)
    masked = biometal.mask_low_quality(record, min_quality=20)

    # Count masked bases (for QC metrics)
    n_count = biometal.count_masked_bases(masked)
    mask_rate = n_count / len(masked.sequence)

    # Quality filter: reject if >10% masked
    if mask_rate < 0.1:
        print(f"Pass QC: {mask_rate*100:.1f}% masked")
    else:
        print(f"Fail QC: {mask_rate*100:.1f}% masked")

    break
```

#### Python: Complete QC Pipeline

```python
import biometal

# Quality control pipeline: trim → filter → mask
stream = biometal.FastqStream.from_path("raw_reads.fq.gz")

passed = 0
failed_quality = 0
failed_length = 0

for record in stream:
    # Step 1: Quality-based trimming (Q20, Trimmomatic-style)
    trimmed = biometal.trim_quality_window(record, min_quality=20, window_size=4)

    # Step 2: Length filter (keep 50-150bp)
    if not biometal.meets_length_requirement(trimmed, min_len=50, max_len=150):
        failed_length += 1
        continue

    # Step 3: Mask remaining low-quality bases
    masked = biometal.mask_low_quality(trimmed, min_quality=20)

    # Step 4: Final QC check (<10% masked)
    mask_rate = biometal.count_masked_bases(masked) / len(masked.sequence)
    if mask_rate > 0.1:
        failed_quality += 1
        continue

    passed += 1
    # Write to output or process further

print(f"QC Results:")
print(f"  Passed: {passed}")
print(f"  Failed (length): {failed_length}")
print(f"  Failed (quality): {failed_quality}")
```

**Use Cases**:
- **Trimming**: Remove low-quality ends before alignment (preserves high-quality core)
- **Masking**: Variant calling pipelines (preserves read structure for alignment)
- **Region extraction**: Extract specific genomic windows or features
- **Reverse complement**: Convert reads to correct strand orientation
- **FASTQ→FASTA**: Convert after quality filtering for downstream tools

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
- Published to crates.io and PyPI

**v1.1.0** (Released November 6, 2025) ✅
- K-mer operations (extraction, minimizers, spectrum)
- Shannon entropy complexity scoring
- Parallel k-mer extraction (opt-in, 2.2× speedup)
- Python bindings for k-mer operations
- Evidence-based design (Entry 034)
- 260 tests, Grade A+

**v1.2.0** (Released November 6, 2025) ✅
- Python bindings for Phase 4 sequence operations (20 functions)
- Complete QC pipelines (trim → filter → mask)
- Trimmomatic-compatible sliding window trimming
- Quality-based masking for variant calling
- 347 tests (260 library + 87 doc), Grade A

**Future Considerations** (Community Driven)
- Python examples & tutorials (Jupyter notebooks)
- Performance benchmarking vs existing tools
- Extended operation coverage (alignment, assembly)
- Additional format support (BAM/SAM, VCF)
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

### 2. BERT Preprocessing Pipeline (DNABert/ML)

```rust
use biometal::{FastqStream, operations::kmer};
use biometal::io::DataSource;

// Stream from SRA (no local copy!)
let source = DataSource::Sra("SRR12345678".to_string());
let stream = FastqStream::new(source)?;

// Extract k-mers for DNABert training
for record in stream {
    let record = record?;

    // Extract overlapping k-mers (Entry 034: scalar-only optimal)
    let kmers = kmer::extract_kmers(&record.sequence, 6);

    // Feed to BERT training pipeline immediately
    // Constant memory even for TB-scale datasets (~5 MB)
}
```

**Python equivalent**:
```python
import biometal

stream = biometal.FastqStream.from_path("dataset.fq.gz")

for record in stream:
    # Extract k-mers for DNABert (k=3, 4, 5, or 6 typical)
    kmers = biometal.extract_kmers(record.sequence, k=6)

    # Feed to model - constant memory!
    model.process(kmers)
```

**For large batches** (10K+ sequences), use parallel extraction:
```python
# Opt-in parallel for 2.2× speedup (Entry 034)
extractor = biometal.KmerExtractor(parallel=True, threads=4)

sequences = [record.sequence for record in batch]
kmers = extractor.extract(sequences, k=6)  # 2.2× faster
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

## FAQ

### Why is the package called `biometal-rs` on PyPI but `biometal` everywhere else?

The `biometal` name was already taken on PyPI when we published v1.0.0, so we used `biometal-rs` (following the Rust convention). However:

- **GitHub repository**: `shandley/biometal`
- **Python import**: `import biometal` (not `biometal_rs`)
- **Rust crate**: `biometal`
- **PyPI package**: `biometal-rs` (install name only)

This means you install with:
```bash
pip install biometal-rs
```

But use it as:
```python
import biometal  # Not biometal_rs!
```

This is a common pattern for Rust-based Python packages and provides the best user experience (clean import name).

### What platforms are supported?

**Pre-built wheels available for**:
- macOS ARM (M1/M2/M3/M4) - Optimized with NEON (16-25× speedup)
- macOS x86_64 (Intel Macs) - Scalar fallback
- Linux x86_64 - Scalar fallback

**Coming soon**:
- Linux ARM (Graviton, Raspberry Pi) - Will be added in v1.0.1

**Build from source**: All other platforms can build from the source distribution (requires Rust toolchain).

### Does it work on Windows?

Currently untested. Building from source may work with the Rust toolchain installed, but we haven't validated it. Community contributions for Windows support are welcome!

### Why ARM-native? What about x86_64?

biometal is designed to democratize bioinformatics by enabling world-class performance on **consumer hardware**. Modern ARM laptops (like MacBooks with M-series chips) cost $1,400 vs $50,000+ for traditional HPC servers.

**Performance philosophy**:
- **Mac ARM** (M1/M2/M3/M4): Optimized target - 16-25× NEON speedup
- **Other platforms**: Correct, production-ready code with scalar fallback

The library works great on x86_64 (all tests pass), it's just not specifically optimized for it. Our mission is enabling field researchers, students, and small labs in LMICs to do cutting-edge work on affordable hardware.

### How do I get support?

- **Bug reports**: [GitHub Issues](https://github.com/shandley/biometal/issues)
- **Questions**: [GitHub Discussions](https://github.com/shandley/biometal/discussions)
- **Documentation**: [https://docs.rs/biometal](https://docs.rs/biometal)

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

**Status**: v1.2.0 - Python Phase 4 Bindings 🎉
**Released**: November 6, 2025
**Grade**: A (rust-code-quality-reviewer compatible)
**Tests**: 347 passing (260 library + 87 doc)
**Python Functions**: 40+ (core ops + k-mers + Phase 4)
**Evidence Base**: 1,357 experiments, 40,710 measurements
**Mission**: Democratizing bioinformatics compute
