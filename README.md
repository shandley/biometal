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

## 🎉 NEW in v1.11.0: GenBank & BLAST Parsers

biometal now supports **14+ bioinformatics file formats** with production-ready streaming parsers:

**Sequences & Reads**:
- **FASTQ/FASTA**: Read sequences with quality scores
- **BAM/SAM**: Binary alignment format with indexing (BAI)
- **CRAM**: Reference-based compression (ARM NEON optimized)
- **GenBank**: NCBI annotated sequences with features **[NEW v1.11.0]**

**Annotations & Features**:
- **BED/narrowPeak**: Genomic intervals and ChIP-seq peaks (ENCODE)
- **GFF3**: Hierarchical gene features (genes, mRNAs, exons, CDS)
- **GTF**: Gene annotations for RNA-seq (GENCODE, Ensembl)

**Variants & Alignments**:
- **VCF**: Genetic variants (SNPs, indels, structural variants)
- **PAF**: minimap2 pairwise alignments (long-read analysis)
- **BLAST tabular**: Sequence alignment results (outfmt 6/7) **[NEW v1.11.0]**

**Graphs & Assembly**:
- **GFA**: Assembly graphs (pangenomes, read overlap graphs)

**Indices**:
- **FAI**: FASTA index for O(1) sequence lookup
- **TBI**: Tabix index for O(log n) region queries
- **BAI/CSI**: BAM/CRAM indices for random access

All formats support:
- ✅ Streaming architecture (constant ~5 MB memory)
- ✅ **Read AND Write** capabilities (FASTQ, FASTA, BED, GFF3, GTF) **[NEW in Nov 2025]**
- ✅ Automatic compression/decompression (`.gz`, `.bgz` files)
- ✅ Python bindings with optimized memory usage
- ✅ Real-world validation (ENCODE, UCSC, Ensembl, 1000 Genomes)

---

## Quick Start

### Installation

**Rust:**
```toml
[dependencies]
biometal = "1.10"
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

### 📝 NEW: Writing Genomic Data (Nov 2025)

biometal now supports **writing** for tab-delimited formats with automatic compression:

**Write BED intervals:**
```rust
use biometal::formats::bed_writer::BedWriter;
use biometal::formats::bed::Bed6Record;

let mut writer = BedWriter::create("peaks.bed.gz")?;  // Auto-compresses
let record = Bed6Record { /* ... */ };
writer.write_bed6(&record)?;
writer.finish()?;  // IMPORTANT: Flush data
```

**Write GFF3 annotations:**
```rust
use biometal::formats::gff_writer::Gff3Writer;

let mut writer = Gff3Writer::create("genes.gff3.gz")?;
writer.write_record(&record)?;  // Auto-writes header
writer.finish()?;
```

**Write GTF for RNA-seq:**
```rust
use biometal::formats::gtf_writer::GtfWriter;

let mut writer = GtfWriter::create("transcripts.gtf.gz")?;
writer.write_record(&record)?;  // Validates required attributes
writer.finish()?;
```

All writers support:
- ✅ Automatic `.gz`/`.bgz` compression (fast cloudflare_zlib backend)
- ✅ Streaming to `stdout` for pipelines
- ✅ Comprehensive validation before writing
- ✅ Constant memory (write terabyte-scale files)

---

## 📚 Documentation

### Start Here
- **📘 [User Guide](docs/USER_GUIDE.md)** - Comprehensive guide: installation, core concepts, common workflows, troubleshooting, and migration from pysam/samtools **(NEW - v1.6.0)**

### In-Depth Resources
- **📖 [DeepWiki AI Docs](https://deepwiki.com/shandley/biometal)** - AI-assisted documentation with Q&A
- **📓 [Interactive Tutorials](notebooks/)** - Jupyter notebooks with real workflows
- **🦀 [API Reference](https://docs.rs/biometal)** - Full Rust documentation
- **🐍 [Python Guide](docs/PYTHON.md)** - Python-specific documentation
- **🧬 [BAM API Reference](docs/BAM_API.md)** - Complete BAM/SAM parser API (v1.4.0)
- **⚡ [BAM Performance Guide](docs/BAM_PERFORMANCE.md)** - Benchmarks and optimization (v1.4.0)
- **📐 [Architecture](docs/ARCHITECTURE.md)** - Technical design details
- **❓ [FAQ](FAQ.md)** - Frequently asked questions

---

## 📓 Interactive Tutorials

Learn biometal through hands-on Jupyter notebooks (5 complete, ~2.5 hours):

| Notebook | Duration | Topics |
|----------|----------|--------|
| [01. Getting Started](notebooks/01_getting_started.ipynb) | 15-20 min | Streaming, GC content, quality analysis |
| [02. Quality Control](notebooks/02_quality_control_pipeline.ipynb) | 30-40 min | Trimming, filtering, masking (v1.2.0) |
| [03. K-mer Analysis](notebooks/03_kmer_analysis.ipynb) | 30-40 min | ML preprocessing, DNABert (v1.1.0) |
| [04. Network Streaming](notebooks/04_network_streaming.ipynb) | 30-40 min | HTTP streaming, public data (v1.0.0) |
| [05. BAM Alignment Analysis](notebooks/05_bam_alignment_analysis.ipynb) | 30-40 min | BAM parsing, 4× speedup, filtering (v1.2.0+) |
| [06. BAM Production Workflows](notebooks/06_bam_production_workflows.ipynb) | 45-60 min | Tag parsing, QC statistics, production pipelines (v1.4.0) |

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
- **BAM/SAM parser**: Production-ready with 8.4× speedup via parallel BGZF + NEON + cloudflare_zlib
  - 5.82 million records/sec throughput
  - 92.0 MiB/s compressed file processing (+67% from cloudflare_zlib in v1.7.0)
  - Constant ~5 MB memory (streams terabyte-scale alignments)
  - **Python bindings (v1.3.0)**: CIGAR operations, SAM writing, alignment metrics
  - **Production polish (v1.4.0)**: Tag convenience methods, statistics functions
    - 6 tag accessors: `edit_distance()`, `alignment_score()`, `read_group()`, etc.
    - 4 statistics functions: `insert_size_distribution()`, `edit_distance_stats()`, `strand_bias()`, `alignment_length_distribution()`
  - **NEON optimization (v1.5.0)**: ARM SIMD sequence decoding (4.62× faster)
  - **BAI index (v1.6.0)**: Indexed region queries with 1.68-500× speedup
    - O(log n) random access to BAM files
    - Near-zero overhead (<1ms index loading)
    - Speedup scales with file size (10-500× for 1-10 GB files)
- **Format Library (v1.8.0)**: Production-ready parsers for genomic annotation and assembly formats
  - **BED (Browser Extensible Data)**: Genomic intervals with streaming architecture
    - BED3/6/12 format support
    - 0-based half-open coordinate system
    - Constant memory (~5 MB) for terabyte-scale peak files
  - **GFA (Graphical Fragment Assembly)**: Assembly graph format
    - Segment, Link, Path record types
    - Graph connectivity validation
    - Streaming architecture for large assembly graphs
  - **VCF (Variant Call Format)**: Genetic variant data
    - VCF 4.2 specification compliance
    - Header parsing with sample/contig/INFO extraction
    - SNP/indel classification
    - Multi-allelic variant support
  - **GFF3 (General Feature Format)**: Hierarchical gene annotations
    - 1-based inclusive coordinate system
    - Parent-child relationship tracking (gene → mRNA → exon/CDS)
    - Attribute parsing with convenience methods
    - Coordinate conversion to BED (0-based)
  - **Testing**: 23 property-based tests + 6 real-world integration tests
  - **Python bindings**: Full streaming API for all formats
- **60+ Python functions** for bioinformatics workflows

---

## Performance Highlights

| Operation | Scalar | Optimized | Speedup |
|-----------|--------|-----------|---------|
| Base counting | 315 Kseq/s | 5,254 Kseq/s | **16.7× (NEON)** |
| GC content | 294 Kseq/s | 5,954 Kseq/s | **20.3× (NEON)** |
| Quality filter | 245 Kseq/s | 6,143 Kseq/s | **25.1× (NEON)** |
| **BAM parsing** | **~11 MiB/s** | **92.0 MiB/s** | **8.4× (BGZF + NEON + cloudflare_zlib v1.7.0)** |

| Dataset Size | Traditional | biometal | Reduction |
|--------------|-------------|----------|-----------|
| 100K sequences | 134 MB | 5 MB | 96.3% |
| 1M sequences | 1,344 MB | 5 MB | 99.5% |
| **5TB dataset** | **5,000 GB** | **5 MB** | **99.9999%** |

**📊 [Comprehensive Benchmark Comparison vs samtools/pysam →](benchmarks/comparison/BENCHMARK_COMPARISON.md)**

---

## Platform Support

| Platform | Performance | Tests | Status |
|----------|-------------|-------|--------|
| **Mac ARM** (M1-M4) | **16-25× speedup** | ✅ 670/670 | Optimized |
| **AWS Graviton** | 6-10× speedup | ✅ 670/670 | Portable |
| **Linux x86_64** | 1× (scalar) | ✅ 670/670 | Portable |

*Test count: 669 library tests passing + 1 ignored (670 total, 100% pass rate) + 23 property-based tests*

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
**BAM/SAM** (Integrated Nov 8, 2025) ✅ - Native streaming alignment parser with parallel BGZF (4× speedup)
**v1.3.0** (Released Nov 9, 2025) ✅ - Python BAM bindings with CIGAR operations and SAM writing
**v1.4.0** (Released Nov 9, 2025) ✅ - BAM tag convenience methods and statistics functions
**v1.5.0** (Released Nov 9, 2025) ✅ - ARM NEON sequence decoding (+27.5% BAM parsing speedup)
**v1.6.0** (Released Nov 10, 2025) ✅ - BAI index support (indexed region queries, 1.68-500× speedup)
**v1.7.0** (Released Nov 13, 2025) ✅ - cloudflare_zlib backend (1.67× decompression, 2.29× compression speedups)
**v1.8.0** (Released Nov 13, 2025) ✅ - Format library (BED, GFA, VCF, GFF3) with property-based testing
  - 4 production-ready format parsers with streaming architecture
  - 23 property-based tests + 6 real-world integration tests
  - Tested against ENCODE, UCSC, Ensembl, 1000 Genomes data
  - Full Python bindings for all formats

**Next** (Planned):
- CSI index support (for references >512 Mbp)
- Extended tag parsing (full type support)
- Additional alignment statistics
- Community feedback & benchmarking

**Future** (Community Driven):
- Extended operations (alignment, assembly)
- Additional formats (BCF, CRAM)
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

### BAM Alignment Analysis (v1.4.0)

```python
import biometal

# Stream BAM file with constant memory (~5 MB)
reader = biometal.BamReader.from_path("alignments.bam")

for record in reader:
    # Access alignment details
    print(f"{record.name}: MAPQ={record.mapq}, pos={record.position}")

    # NEW v1.4.0: Tag convenience methods
    edit_dist = record.edit_distance()  # NM tag
    align_score = record.alignment_score()  # AS tag
    read_group = record.read_group()  # RG tag
    print(f"  Edit distance: {edit_dist}, Score: {align_score}, RG: {read_group}")

    # CIGAR operations (v1.3.0)
    for op in record.cigar:
        if op.is_insertion() and op.length >= 5:
            print(f"  Found {op.length}bp insertion")

# NEW v1.4.0: Built-in statistics functions
# Insert size distribution (paired-end QC)
dist = biometal.insert_size_distribution("alignments.bam")
print(f"Mean insert size: {sum(s*c for s,c in dist.items())/sum(dist.values()):.1f}bp")

# Edit distance statistics (alignment quality)
stats = biometal.edit_distance_stats("alignments.bam")
print(f"Mean edit distance: {stats['mean']:.2f} mismatches/read")

# Strand bias (variant calling QC)
bias = biometal.strand_bias("alignments.bam", reference_id=0, position=1000)
print(f"Strand bias at chr1:1000: {bias['ratio']:.2f}:1")

# Alignment length distribution (RNA-seq QC)
lengths = biometal.alignment_length_distribution("alignments.bam")
print(f"Intron-spanning reads: {sum(c for l,c in lengths.items() if l > 1000)}")
```

### BAI Indexed Region Queries (v1.6.0)

```python
import biometal

# Load BAI index for fast random access
index = biometal.BaiIndex.from_path("alignments.bam.bai")

# Query specific genomic region (1.68× faster than full scan for small files)
# Speedup increases dramatically with file size (10-500× for 1-10 GB files)
for record in biometal.BamReader.query_region(
    "alignments.bam",
    index,
    "chr1",
    1000000,  # start position
    2000000   # end position
):
    # Only reads overlapping region are returned
    if record.is_mapped and record.mapq >= 30:
        print(f"{record.name}: {record.position}-{record.reference_end()}")

# Reuse index for multiple queries (index loading: <1ms overhead)
regions = [
    ("chr1", 1000000, 2000000),
    ("chr1", 5000000, 6000000),
    ("chr2", 100000, 200000),
]

for chrom, start, end in regions:
    count = sum(1 for _ in biometal.BamReader.query_region(
        "alignments.bam", index, chrom, start, end
    ))
    print(f"{chrom}:{start}-{end}: {count} reads")

# Full workflow: Coverage calculation for specific region
from collections import defaultdict

coverage = defaultdict(int)
for record in biometal.BamReader.query_region(
    "alignments.bam", index, "chr1", 1000, 2000
):
    if record.is_mapped and record.position is not None:
        # Calculate coverage from CIGAR
        pos = record.position
        for op in record.cigar:
            if op.consumes_reference():
                for i in range(op.length):
                    coverage[pos] += 1
                    pos += 1

print(f"Mean coverage: {sum(coverage.values())/len(coverage):.1f}×")
```

**Performance Characteristics:**
- Index loading: < 1ms (negligible overhead)
- Small region query (1 Kbp): ~11 ms vs 18 ms full scan (1.68× speedup)
- Speedup scales with file size:
  - 1 MB file: 1.7× speedup
  - 100 MB file: 10-20× speedup
  - 1 GB file: 50-100× speedup
  - 10 GB file: 200-500× speedup

### Format Library: BED/GFA/VCF/GFF3 (v1.8.0)

```python
import biometal

# BED: Parse genomic intervals (ChIP-seq peaks, gene annotations)
stream = biometal.Bed6Stream.from_path("peaks.bed.gz")
for record in stream:
    print(f"{record.chrom}:{record.start}-{record.end} score={record.score}")
    length = record.length()
    if length > 1000:
        print(f"  Long peak: {length}bp")

# GFA: Parse assembly graphs (genome assembly, pangenomes)
stream = biometal.GfaStream.from_path("assembly.gfa")
segments = []
for record in stream:
    if isinstance(record, biometal.GfaSegment):
        segments.append(record)
        print(f"Segment {record.name}: {len(record.sequence)}bp")

# VCF: Parse genetic variants (SNPs, indels)
stream = biometal.VcfStream.from_path("variants.vcf.gz")
header = stream.header()  # Note: header() not parse_header()
print(f"VCF version: {header.fileformat}, Samples: {len(header.samples)}")

for variant in stream:
    if variant.quality and variant.quality > 30:
        print(f"{variant.chrom}:{variant.pos} {variant.reference}→{variant.alternate[0]}")
        if variant.is_snp():
            print(f"  SNP with quality {variant.quality}")

# GFF3: Parse hierarchical gene annotations (genes, mRNAs, exons, CDS)
stream = biometal.Gff3Stream.from_path("annotations.gff3.gz")
for feature in stream:
    if feature.feature_type == "gene":
        gene_id = feature.get_id()
        length = feature.length()  # 1-based inclusive coordinates
        print(f"Gene {gene_id}: {length}bp on {feature.strand}")

    elif feature.feature_type == "exon":
        parent = feature.get_parent()
        # Note: interval() method not available in Python bindings
        # Use feature.start and feature.end directly (1-based inclusive)
        print(f"  Exon of {parent}: {feature.start}-{feature.end}")
```

**Format Library Features:**
- **Streaming architecture**: Constant ~5 MB memory for all formats
- **Production-ready**: Tested against real ENCODE, UCSC, Ensembl, 1000 Genomes data
- **Property-based testing**: 23 tests validating format invariants (round-trip parsing, coordinate systems, specification compliance)
- **Real-world validation**: 6 integration tests with production files (61,547 GFF3 features, 1,000 UCSC genes, 10 VCF variants)
- **Python bindings**: Full streaming API with Pythonic interfaces

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
<strong>Status:</strong> v1.11.0 released 🚀 (Phase 2 COMPLETE)<br>
<strong>Latest:</strong> GenBank + BLAST parsers, Grade A code quality (Nov 16, 2025)<br>
<strong>Tests:</strong> 670 total (669 passing + 1 ignored, 100% pass rate) + 23 property-based<br>
<strong>Performance:</strong> 5.82M records/sec, 92.0 MiB/s throughput, 50-60% Python memory reduction<br>
<strong>Python Functions:</strong> 100+ (14 formats READ+WRITE, 4 READ-only, 3 indices)<br>
<strong>Evidence Base:</strong> 1,357 experiments, 40,710 measurements
</p>
