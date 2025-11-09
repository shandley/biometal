# Changelog

All notable changes to biometal will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.5.0] - 2025-11-09

### 🚀 Major Performance: ARM NEON Sequence Decoding

This release delivers a **+27.5% overall BAM parsing improvement** through ARM NEON SIMD optimization of sequence decoding, validated through comprehensive microbenchmarking and evidence-based methodology.

**Performance Gains**:
- **Sequence decoding**: 4.62× faster (6.65 ms → 1.44 ms per 100K records)
- **Overall BAM parsing**: +27.5% faster (21.995 ms → 17.192 ms per 100K records)
- **Throughput**: 43.0 MiB/s → 55.1 MiB/s (+28% improvement)

**Evidence-Based Approach**:
- Isolated microbenchmarking revealed sequence decoding = 30.2% of CPU time (2× Rule 1 threshold)
- NEON implementation validated with 23 tests (8 NEON-specific + 15 platform-agnostic)
- Property-based testing (proptest) ensures correctness across all edge cases
- Code review by rust-code-quality-reviewer: **PRODUCTION-READY**

### Added

#### Core: ARM NEON Sequence Decoder (src/io/bam/sequence_neon.rs)

**NEON SIMD Implementation** - 4-bit to ASCII decoding (247 lines):
- Uses `vqtbl1q_u8` for 16-entry vector table lookup
- Processes 32 bases per NEON iteration (16-byte chunks)
- Interleaves high/low nibbles with `vzip1q_u8`/`vzip2q_u8`
- Scalar fallback for remaining bases (<32) and non-ARM platforms
- Platform dispatch: ARM64 (NEON) vs x86_64 (scalar fallback)

**Evidence** (experiments/bam-simd-sequence-decoding/):
- **PROPOSAL.md**: Experimental design, go/no-go criteria, timeline
- **RESEARCH_LOG.md**: Daily progress, benchmark results, analysis
- **FINDINGS.md**: Complete experiment summary, learnings, recommendations

**Key Learnings**:
- **Memory-bound operations**: NEON provides 4-8× speedup (not 16-25×) due to memory allocation overhead
- **Bottleneck cascade**: Optimizing BGZF exposed sequence decoding (6% → 30.2% CPU time)
- **Overall impact**: Even "modest" NEON speedups deliver major improvements on large bottlenecks

**Rule 1 Refinement** (OPTIMIZATION_RULES.md):
- **Compute-bound**: 16-25× speedup (base counting, GC content)
- **Memory-bound**: 4-8× speedup (sequence decode, quality scores)
- **Memory+allocation bound**: 3-5× speedup (frequent small allocations)

#### Performance: Benchmark Suite Enhanced

**New Benchmarks**:
- `sequence_decode`: Isolated sequence decoding microbenchmark
  - Single read (100bp): 66.7 ns (scalar) → 14.4 ns (NEON)
  - Bulk (100K × 100bp): 6.65 ms (scalar) → 1.44 ms (NEON)
  - Throughput: 1.50 Gbases/s → 6.95 Gbases/s

**Updated Benchmarks**:
- `bam_parsing`: Overall BAM parsing with NEON enabled
  - Parse 100K records: 21.995 ms → 17.192 ms (-21.5%)
  - Throughput: 43.031 MiB/s → 55.053 MiB/s (+27.5%)

### Changed

#### Performance: BAM Parsing Bottleneck Distribution

**Before v1.5.0** (Scalar sequence decoding):
- Sequence decoding: 30.2% (PRIMARY BOTTLENECK)
- BGZF decompression: ~30-35%
- Record parsing: ~15-20%
- Other: ~15-20%

**After v1.5.0** (NEON sequence decoding):
- BGZF decompression: ~30-35% (NEW PRIMARY BOTTLENECK)
- Record parsing: ~20-25%
- Sequence decoding: ~8% (OPTIMIZED ✅)
- Other: ~20-30%

**Next Optimization Target**: BGZF parallel decompression (6.5× expected, Rule 3)

#### Documentation: Evidence Base Updated

**OPTIMIZATION_RULES.md** - Rule 1 refinement:
- Added memory-bound vs compute-bound speedup expectations
- Documented sequence decode experiment as memory-bound evidence
- Updated when-to-apply guidance with realistic NEON expectations

### Technical Details

#### Platform Support

**ARM64 (M1/M2/M3/M4, Graviton)**:
- NEON SIMD implementation (4.62× sequence decode speedup)
- 23 tests passing (including NEON-specific tests)
- Production-ready, code-reviewed

**x86_64 (Intel/AMD)**:
- Scalar fallback implementation (portable, correct)
- Same 23 tests passing
- Full compatibility maintained

#### Memory Safety

**All unsafe code validated**:
- Bounds checking before all NEON operations
- `set_len()` only called after data is written
- No UB possible with any input data
- rust-code-quality-reviewer: "All unsafe code is sound"

#### Testing

**Test Coverage** (23 total):
- **NEON-specific**: 8 tests (single base, 32 bases, 64 bases, 100 bases, etc.)
- **Platform-agnostic**: 15 tests (odd-length, empty, all IUPAC codes, etc.)
- **Property-based**: 5 proptest properties (roundtrip, length invariants, etc.)

**Result**: 23/23 passing on both ARM64 and x86_64

### Performance Validation

#### Benchmarking Methodology

**Isolated Microbenchmark** (sequence_decode):
- Measures ONLY `decode_sequence()` function
- N=30 samples for statistical significance
- Controls for allocation overhead, cache effects
- Validated against total BAM parsing time

**Overall Impact Validation**:
```
Sequence time saved: 30.2% × (1 - 1/4.62) = 23.5%
Overall speedup: 1 / (1 - 0.235) = 1.31× = +31% faster
Measured: +27.5% faster (within 3.5% of prediction!)
```

#### Success Metrics

| Criterion | Target | Actual | Status |
|-----------|--------|--------|--------|
| CPU Time % | ≥15% | **30.2%** | ✅ 2× threshold |
| Overall Speedup | ≥5% | **+27.5%** | ✅ 5.5× target! |
| Correctness | 100% | 23/23 tests | ✅ Perfect |
| Code Quality | Production | PRODUCTION-READY | ✅ Approved |

### Acknowledgments

This optimization follows biometal's evidence-based methodology:
1. **Profile first**: Isolated sequence decoding CPU time (30.2%)
2. **Validate go/no-go**: Strong GO (2× the 15% threshold)
3. **Implement carefully**: ARM NEON with scalar fallback
4. **Test thoroughly**: 23 tests + property-based validation
5. **Review rigorously**: Code quality review → PRODUCTION-READY
6. **Measure impact**: +27.5% faster (5.5× the target!)

**Evidence**: experiments/bam-simd-sequence-decoding/ (PROPOSAL, RESEARCH_LOG, FINDINGS)

## [1.4.0] - 2025-11-09

### 🏭 Production Polish: Extended Tag Parsing and Statistics

Complete production polish for BAM/SAM parser with extended tag parsing, built-in statistics functions, and comprehensive documentation. Enables production-grade QC workflows with streamlined tag access and optimized statistics calculations.

**Grade**: A (all tests passing)
**Tests**: 50 Python BAM tests passing (+14 new tests, +39% increase)
**New Python APIs**: 10 (6 tag methods + 4 statistics functions)
**Documentation**: 3 new files (API reference, performance guide, production workflows)
**Performance**: 43.0 MiB/s, 4.4M records/sec, ~5 MB constant memory

### Added

#### Python: Extended Tag Parsing (src/python/bam.rs)

**Tag Convenience Methods** - Streamlined access to common BAM tags:

- **`get_int(tag)`** - Direct integer tag access with type checking
- **`get_string(tag)`** - Direct string tag access with type checking
- **`edit_distance()`** - NM tag convenience method (mismatches per read)
- **`alignment_score()`** - AS tag convenience method (aligner score)
- **`read_group()`** - RG tag convenience method (sample identification)
- **`md_string()`** - MD tag convenience method (mismatch details)

**Benefits**:
- Cached tag parsing for repeated access (~20-30% faster)
- Type-safe API with clear error messages
- Pythonic interface for common QC workflows
- Zero-copy when possible

**Use Cases**:
- Quality control pipelines (edit distance filtering)
- Multi-sample processing (read group demultiplexing)
- Alignment quality assessment (alignment score analysis)
- Variant calling preparation (MD string parsing)

#### Python: Statistics Functions (src/python/bam.rs)

**Built-in QC Statistics** - Optimized functions for common workflows:

**`insert_size_distribution(path, reference_id=None)`**:
- Calculate insert size distribution for paired-end reads
- Returns dict mapping insert size → count
- Filters for properly paired, primary alignments
- **Use Case**: Library preparation QC, adapter contamination detection

**`edit_distance_stats(path, reference_id=None)`**:
- Comprehensive edit distance (NM tag) statistics
- Returns dict with mean, median, min, max, distribution
- **Use Case**: Alignment quality assessment, sequencing error analysis

**`strand_bias(path, reference_id, position, window_size=1)`**:
- Calculate forward/reverse strand balance at position
- Returns dict with forward, reverse counts, ratio, percentages
- **Use Case**: Variant calling QC, artifact detection

**`alignment_length_distribution(path, reference_id=None)`**:
- Calculate reference alignment length distribution
- Detects intron-spanning reads (RNA-seq)
- **Use Case**: RNA-seq QC, splice junction analysis

**Performance**: 1.5-2× faster than manual Python loops (optimized Rust implementation)

#### Documentation: Comprehensive Reference (3 new files)

**docs/BAM_API.md** (805 lines):
- Complete API reference with examples
- All classes documented (BamReader, BamRecord, BamHeader, CigarOp, Tag, SamWriter)
- All 10 new v1.4.0 methods with usage examples
- All 7 statistics functions with parameter descriptions
- Performance tips section (DO/DON'T patterns)
- Error handling guide
- Complete working example

**docs/BAM_PERFORMANCE.md** (659 lines):
- Benchmark methodology and results
- Performance characteristics (parallel BGZF, streaming, zero-copy)
- Memory characteristics (constant ~5 MB profile)
- Scaling characteristics (linear with file size)
- Optimization guide (5 best practices)
- Profiling results (Phase 0 bottleneck analysis)
- Comparison notes (vs samtools, pysam, noodles)
- Troubleshooting guide
- Hardware recommendations
- Future optimizations roadmap

**notebooks/06_bam_production_workflows.ipynb** (1,038 lines):
- 5 complete production workflows with examples
- Workflow 1: Comprehensive Quality Control Pipeline
- Workflow 2: Paired-End Insert Size Analysis
- Workflow 3: Variant Calling Preparation
- Workflow 4: RNA-seq Alignment QC
- Workflow 5: Multi-Sample Tag-Based Filtering
- All workflows demonstrate v1.4.0 features
- Executable Jupyter notebook with commentary

### Enhanced

#### Testing (tests/python/test_bam.py)

**14 New Tests** for v1.4.0 features:

**Tag Convenience Tests** (8 tests):
- `test_get_int()` - Integer tag access validation
- `test_get_string()` - String tag access validation
- `test_edit_distance()` - NM tag convenience method
- `test_alignment_score()` - AS tag convenience method
- `test_read_group()` - RG tag convenience method
- `test_md_string()` - MD tag convenience method
- `test_get_int_invalid()` - Error handling for invalid tags
- `test_get_string_type_mismatch()` - Type mismatch handling

**Statistics Tests** (6 tests):
- `test_insert_size_distribution()` - Paired-end insert sizes
- `test_insert_size_by_reference()` - Reference-specific filtering
- `test_edit_distance_stats()` - Comprehensive edit distance stats
- `test_strand_bias()` - Strand balance calculation
- `test_alignment_length_distribution()` - RNA-seq alignment lengths
- `test_statistics_empty_file()` - Edge case handling

**Coverage**: All new APIs tested, property-based testing maintained

#### Examples (examples/bam_advanced_filtering.py)

**7 New Demonstration Functions**:

**Tag Analysis Functions**:
- `analyze_edit_distance()` - Edit distance distribution with statistics
- `filter_by_read_group()` - Multi-sample demultiplexing
- `analyze_alignment_scores()` - Alignment score distribution

**Statistics Demonstrations** (Examples 11-14):
- Example 11: Insert size distribution (paired-end QC)
- Example 12: Edit distance statistics (alignment quality)
- Example 13: Strand bias analysis (variant calling QC)
- Example 14: Alignment length distribution (RNA-seq QC)

**Updated Documentation**:
- Added v1.4.0 feature summary in module docstring
- All examples include usage instructions and output interpretation
- Demonstrates integration with existing workflows

### Performance

**Benchmarks** (M3 MacBook Pro, 100K records, N=30):
- **Throughput**: 43.0 MiB/s (compressed BAM processing)
- **Record Rate**: 4.4 million records/sec
- **Processing Time**: 22 ms / 100K records (~220 ns/record)
- **Memory**: Constant ~5 MB (independent of file size)

**Regression Analysis**:
- v1.3.0 → v1.4.0: -1.6% throughput (43.7 → 43.0 MiB/s)
- **Cause**: Tag convenience methods add caching overhead
- **Decision**: Acceptable tradeoff for improved API ergonomics
- **Impact**: Minimal (< 2%), within acceptable range

**Validation**:
- All benchmarks pass with N=30 samples
- Performance regression within acceptable limits
- Memory footprint unchanged (~5 MB constant)

### Files Changed

**Code**:
- `src/python/bam.rs`: +449 lines (tag methods, statistics functions)
- `src/python/mod.rs`: +4 function registrations
- `tests/python/test_bam.py`: +285 lines (14 new tests)
- `examples/bam_advanced_filtering.py`: +272 lines (7 new functions)

**Documentation**:
- `docs/BAM_API.md`: NEW (805 lines, comprehensive API reference)
- `docs/BAM_PERFORMANCE.md`: NEW (659 lines, performance guide)
- `notebooks/06_bam_production_workflows.ipynb`: NEW (1,038 lines, 5 workflows)

**Total**: +3,512 lines (7 files modified, 3 files added)

### Migration Guide

**From v1.3.0 to v1.4.0**:

No breaking changes. All v1.3.0 code continues to work unchanged.

**New APIs** (optional migration for improved ergonomics):

```python
# v1.3.0: Manual tag access
tag = record.get_tag("NM")
if tag and tag.value.is_int():
    nm = tag.value.as_int()

# v1.4.0: Convenience method (recommended)
nm = record.edit_distance()  # Returns None if not present
```

**New Statistics** (recommended for production workflows):

```python
# v1.4.0: Built-in functions (1.5-2× faster)
import biometal

dist = biometal.insert_size_distribution("alignments.bam")
stats = biometal.edit_distance_stats("alignments.bam")
bias = biometal.strand_bias("alignments.bam", ref_id=0, position=1000)
lengths = biometal.alignment_length_distribution("alignments.bam")
```

## [1.3.0] - 2025-11-09

### 🧬 Python BAM Bindings: CIGAR Operations and SAM Writing

Complete Python bindings for BAM alignment analysis with CIGAR operations, SAM writing, and advanced filtering capabilities. Enables production-grade alignment analysis pipelines in Python with biometal's streaming-first architecture.

**Grade**: A (all tests passing)
**Tests**: 545 passing (354 library + 70 BAM + 121 doc)
**New Python Classes**: 2 (CigarOp, SamWriter)
**New Python Methods**: 10+ (CIGAR analysis, SAM writing, alignment metrics)

### Added

#### Python: CIGAR Operations (src/python/bam.rs)
Complete CIGAR (Compact Idiosyncratic Gapped Alignment Report) analysis:

**PyCigarOp Class**:
- **`op_char`** property - Operation character (M, I, D, N, S, H, P, =, X)
- **`length`** property - Operation length in bases
- **`is_match()`**, **`is_insertion()`**, **`is_deletion()`** - Type checking (9 methods total)
- **`is_ref_skip()`**, **`is_soft_clip()`**, **`is_hard_clip()`**, **`is_padding()`**
- **`is_seq_match()`**, **`is_seq_mismatch()`** - Detailed alignment type
- **`consumes_reference()`** - Returns True if operation advances reference position
- **`consumes_query()`** - Returns True if operation consumes query sequence

**BamRecord CIGAR Methods**:
- **`cigar`** property - Returns list of CigarOp objects
- **`cigar_string()`** - Generate human-readable CIGAR string (e.g., "100M2I50M")
- **`reference_length()`** - Calculate reference span from CIGAR
- **`query_length()`** - Calculate query span from CIGAR
- **`reference_end()`** - Calculate alignment end position

**Use Cases**:
- Structural variant detection (indel identification)
- Alignment quality assessment
- Coverage calculation with CIGAR awareness
- RNA-seq splice junction analysis (N operations)

#### Python: SAM Writing (src/python/bam.rs)
BAM to SAM conversion with optional filtering:

**PySamWriter Class**:
- **`SamWriter.create(path)`** - Create SAM writer
- **`write_header(header)`** - Write SAM header from BamHeader
- **`write_record(record)`** - Write BamRecord in SAM format
- **`close()`** - Close writer and flush buffers

**Use Cases**:
- BAM → SAM conversion for human-readable output
- Region extraction to SAM format
- Filtered alignment export
- Tool integration requiring SAM input

#### Python: Enhanced Coverage Calculation
Updated `calculate_coverage()` function to use CIGAR operations:

- Accurate per-base coverage accounting for insertions, deletions, and clipping
- Properly handles all 9 CIGAR operation types
- Only counts reference-consuming operations (M, D, N, =, X)

**Before**: Simple read start position counting
**After**: CIGAR-aware per-base coverage with proper indel handling

### Enhanced

#### Documentation Updates

**examples/bam_advanced_filtering.py**:
- Added **`analyze_cigar_operations()`** - Distribution of CIGAR operations
- Added **`find_indels()`** - Detect insertions and deletions ≥N bp
- Added **`calculate_alignment_metrics()`** - Alignment quality from CIGAR
- Added **`convert_bam_to_sam()`** - Full BAM → SAM conversion with filtering
- Added **`extract_region_to_sam()`** - Export genomic region to SAM format
- Updated `calculate_coverage()` to use CIGAR operations
- Updated `main()` with 3 new examples (CIGAR analysis, indel detection, SAM writing)

**README.md**:
- Added "BAM Alignment Analysis (v1.3.0)" example section
- Demonstrates CIGAR iteration and type checking
- Shows SAM writing workflow
- Updated Operations Library to mention Python BAM bindings
- Updated roadmap: v1.3.0 marked as "In Development"
- Updated footer: "50+ Python functions (including full BAM support)"

**notebooks/05_bam_alignment_analysis.ipynb**:
- Added Section 11: "CIGAR Operations Analysis (v1.3.0)"
  - CIGAR operation distribution analysis
  - Indel detection with thresholds
  - Alignment metrics calculation
- Added Section 12: "SAM Writing and Format Conversion (v1.3.0)"
  - BAM → SAM conversion with filtering
  - Region extraction to SAM format
  - Human-readable output examples
- Updated Key Takeaways with v1.3.0 features

### Testing

**New Test Suite** (tests/python/test_bam.py):
- 36 comprehensive tests covering all BAM functionality
- CIGAR operations: parsing, type checking, consumption methods
- CIGAR helpers: cigar_string(), reference_length(), query_length(), reference_end()
- SAM writing: header, records, close
- Integration tests: CIGAR-based coverage calculation
- All tests passing (skipped when test data unavailable)

### Fixed

**Documentation Tests**:
- Updated `src/io/bam/mod.rs` example to use correct return type (`biometal::Result<()>`)
- All 121 doc tests now passing

### Performance

No performance changes - focuses on exposing existing functionality to Python.
Underlying BAM parser maintains:
- 4.54 million records/sec throughput
- 43.0 MiB/s compressed file processing
- Constant ~5 MB memory footprint

### Python API Summary

**New in v1.3.0**:
```python
# CIGAR operations
for record in reader:
    for op in record.cigar:
        if op.is_insertion() and op.length >= 5:
            print(f"Found {op.length}bp insertion")

    # Alignment metrics
    ref_len = record.reference_length()
    query_len = record.query_length()
    cigar_str = record.cigar_string()

# SAM writing
writer = biometal.SamWriter.create("output.sam")
writer.write_header(reader.header)
for record in reader:
    if record.is_primary and record.mapq >= 30:
        writer.write_record(record)
writer.close()
```

## [1.2.0] - 2025-11-06

### 🐍 Python Bindings for Phase 4 Sequence Operations

Complete Python bindings for all Phase 4 sequence manipulation primitives, enabling Python users to build production-grade read processing pipelines with biometal's streaming-first architecture.

**Grade**: A (rust-code-quality-reviewer compatible)
**Tests**: 260 library + 87 doc tests passing
**New Python Functions**: 20 (sequence ops, record ops, trimming, masking)

### Added

#### Python: Sequence Operations (src/python/sequence.rs)
Complete Python bindings for core sequence manipulation:

- **`reverse_complement(sequence)`** - Standard molecular biology operation
- **`complement(sequence)`** - Complement only (preserves 5'→3' orientation)
- **`reverse(sequence)`** - Reverse only (no complementation)
- **`is_valid_dna(sequence)`** - DNA validation (ACGTN)
- **`is_valid_rna(sequence)`** - RNA validation (ACGUN)
- **`count_invalid_bases(sequence)`** - QC metric for invalid bases

**Use Cases**:
- Strand orientation correction
- Sequence validation in QC pipelines
- RNA/DNA format conversion

#### Python: Record Operations (src/python/record_ops.rs)
FASTQ record manipulation functions:

- **`extract_region(record, start, end)`** - Extract subsequence [start, end)
- **`reverse_complement_record(record)`** - RC with quality alignment preservation
- **`sequence_length(record)`** - Get sequence length
- **`meets_length_requirement(record, min_len, max_len)`** - Length filtering
- **`to_fasta_record(record)`** - FASTQ→FASTA conversion (drops quality)

**Use Cases**:
- Genomic window extraction
- Read preprocessing and filtering
- Format conversion after QC

#### Python: Trimming Operations (src/python/trimming.rs)
Fixed position and quality-based trimming:

**Fixed Position**:
- **`trim_start(record, bases)`** - Remove N bases from 5' end
- **`trim_end(record, bases)`** - Remove N bases from 3' end
- **`trim_both(record, start_bases, end_bases)`** - Trim both ends

**Quality-Based** (Phred+33, Illumina 1.8+):
- **`trim_quality_end(record, min_quality)`** - Trim low-quality 3' end
- **`trim_quality_start(record, min_quality)`** - Trim low-quality 5' end
- **`trim_quality_both(record, min_quality)`** - Trim both ends (single-pass)
- **`trim_quality_window(record, min_quality, window_size)`** - Trimmomatic-style sliding window

**Use Cases**:
- Pre-alignment QC (remove low-quality ends)
- Adapter trimming
- Trimmomatic-compatible workflows

#### Python: Masking Operations (src/python/masking.rs)
Quality-based masking (preserves read length):

- **`mask_low_quality(record, min_quality)`** - Replace bases with 'N'
- **`count_masked_bases(record)`** - Count N's for QC metrics

**Use Cases**:
- Variant calling pipelines (preserves read structure for alignment)
- Alternative to trimming when length preservation is critical

### Documentation

#### README.md Updates
New "Sequence Manipulation Operations (Phase 4)" section with 5 comprehensive Python examples:

1. **Sequence Operations** - reverse_complement, complement, reverse, validation
2. **Record Operations** - extract_region, RC, length filtering, FASTQ→FASTA
3. **Quality-Based Trimming** - Fixed position and quality-based trimming
4. **Quality-Based Masking** - Masking with QC metrics
5. **Complete QC Pipeline** - Real-world trim → filter → mask workflow

#### Python Module Organization
```
src/python/
├── mod.rs          # PyModule registration (20 new functions)
├── records.rs      # PyFastqRecord with to_fastq_record() conversion
├── sequence.rs     # Sequence operations (6 functions)
├── record_ops.rs   # Record manipulation (5 functions)
├── trimming.rs     # Trimming operations (7 functions)
└── masking.rs      # Masking operations (2 functions)
```

### Technical Details

#### PyFastqRecord Conversion Pattern
Added `to_fastq_record()` conversion method to bridge Python wrapper and Rust internal types:

```rust
impl PyFastqRecord {
    pub(crate) fn to_fastq_record(&self) -> FastqRecord {
        FastqRecord {
            id: self.id.clone(),
            sequence: self.sequence.clone(),
            quality: self.quality.clone(),
        }
    }
}
```

This enables zero-cost abstraction while maintaining clean Python API.

#### Error Handling
All Python functions properly handle Rust errors:
- `Result<T, BiometalError>` → `PyResult<T>`
- Clear error messages for Python users
- Maintains production-grade quality standards

### Testing

**Rust Tests**:
- 260 library tests passing
- 87 documentation tests passing
- All Phase 4 operations validated

**Python Build**:
- `maturin build` successful (zero warnings)
- Wheel generation confirmed

### Examples

#### Complete QC Pipeline (Python)
```python
import biometal

stream = biometal.FastqStream.from_path("raw_reads.fq.gz")

for record in stream:
    # Step 1: Quality trimming (Trimmomatic-style)
    trimmed = biometal.trim_quality_window(record, min_quality=20, window_size=4)

    # Step 2: Length filter (50-150bp)
    if not biometal.meets_length_requirement(trimmed, min_len=50, max_len=150):
        continue

    # Step 3: Mask remaining low-quality bases
    masked = biometal.mask_low_quality(trimmed, min_quality=20)

    # Step 4: Final QC check (<10% masked)
    mask_rate = biometal.count_masked_bases(masked) / len(masked.sequence)
    if mask_rate > 0.1:
        continue

    # Write to output or process further
```

### Compatibility

- **Python Versions**: 3.9-3.14 (tested)
- **Platforms**: macOS ARM (optimized), macOS x86_64, Linux x86_64
- **PyO3**: 0.27 (latest stable)

### Performance

All operations maintain constant memory streaming architecture:
- **Memory footprint**: ~5 MB regardless of dataset size
- **Throughput**: Same as Rust (zero-copy where possible)
- **Quality ops**: Phred+33 standard (Illumina 1.8+)

### Migration Notes

**For Rust Users**:
- No breaking changes to Rust API
- All existing code continues to work

**For Python Users**:
- Install: `pip install biometal-rs`
- Import: `import biometal`
- All 20 new functions available immediately

---

## [1.1.0] - 2025-11-06

### 🧬 K-mer Operations & Complexity Scoring Release

Production-ready k-mer operations and sequence complexity analysis with comprehensive benchmarking and property-based testing.

**Grade**: A+ (rust-code-quality-reviewer, post-refactoring)
**Tests**: 260 passing (254 unit/integration + 6 property-based)
**Evidence Base**: ASBB Entry 034 (k-mer operations), Entry 020-025 (base counting NEON)

### Added

#### K-mer Operations (src/operations/kmer.rs) - Entry 034
**Evidence-Based Design**: K-mer operations are data-structure-bound (hash+HashMap), not compute-bound.

- **Simple extraction**: `extract_kmers(sequence, k)` - Overlapping k-mer extraction
- **Streaming iterator**: `kmer_iter(sequence, k)` - Zero-copy, constant memory
- **Minimizers**: `extract_minimizers(sequence, k, w)` - minimap2-style sketching
- **Spectrum counting**: `kmer_spectrum(sequences, k)` - Frequency analysis
- **Parallel extraction**: `KmerExtractor::with_parallel(threads)` - Opt-in 2.2× speedup

**Performance** (Entry 034 findings):
- **Minimizers**: 1.02-1.26× (NEON/Parallel) → Scalar-only optimal
- **Spectrum**: 0.95-1.88× (sometimes SLOWER with parallel!) → Scalar-only
- **Extraction**: 2.19-2.38× (Parallel-4t for ≥1000 sequences) → Opt-in parallel

**API Design**:
- Conservative defaults (scalar, simple)
- Opt-in parallelism via `KmerExtractor::with_parallel(4)`
- `PARALLEL_THRESHOLD` constant (1000 sequences) publicly exposed
- `will_use_parallel()` method for behavior transparency
- Thread count automatically capped at 4 (empirically optimal)

**Python Bindings** (src/python/kmers.rs):
- `extract_kmers(sequence, k)` - Returns list of k-mers
- `extract_minimizers(sequence, k, w)` - Returns list of minimizer dicts
- `kmer_spectrum(sequences, k)` - Returns frequency dict
- `KmerExtractor(parallel=False, threads=4)` - Configurable extractor class

**Use Cases**:
- BERT/DNABert preprocessing (k=3-6)
- minimap2-style indexing (minimizers)
- De Bruijn graph assembly (k-mer spectrum)
- K-mer frequency analysis

#### Complexity Scoring (src/operations/complexity.rs)
- **Shannon entropy**: `complexity_score(sequence)` - Returns 0.0-1.0
- **NEON-optimized**: Reuses `count_bases()` (16.7× speedup)
- **Scalar fallback**: `complexity_score_scalar()` for testing

**Performance**:
- Inherits NEON speedup from base_counting (16.7×)
- Optimized entropy calculation (single-pass)

**Use Cases**:
- Metagenomics QC (filter low-complexity regions)
- Read quality assessment
- Homopolymer detection

#### Record Operations Extensions (src/operations/record_ops.rs)
- **FASTA conversion**: `to_fasta_record(fastq_record)` - Drop quality scores

### Code Quality Improvements

Addressed **all** rust-code-quality-reviewer recommendations (8 issues):

**HIGH Priority** (1):
- ✅ Removed unused `HashMap` import in Python bindings

**MEDIUM Priority** (4):
- ✅ Fixed `unwrap()` in documentation example (use `if let` pattern)
- ✅ Exposed `PARALLEL_THRESHOLD` as public constant
- ✅ Added `will_use_parallel()` method for API transparency
- ✅ Added 6 property-based tests validating k-mer invariants

**LOW Priority** (3):
- ✅ Improved thread pool error handling (fallback to scalar on failure)
- ✅ Removed unused `complexity_score_neon()` function
- ✅ Added comprehensive k-mer benchmarks (`benches/kmer_operations.rs`)

### Testing

**Property-Based Tests** (6 new tests using proptest):
- `prop_kmer_count`: Validates k-mer count formula `len(kmers) == max(0, len(seq) - k + 1)`
- `prop_iter_matches_extract`: Iterator and extract produce identical results
- `prop_spectrum_sum`: Spectrum frequencies sum to total k-mer count
- `prop_minimizer_uniqueness`: Consecutive minimizers have different positions
- `prop_parallel_matches_scalar`: Parallel extraction correctness
- `prop_edge_cases_no_panic`: Robustness testing for all inputs

**Benchmarks**:
- `benches/kmer_operations.rs` (182 LOC, 6 benchmark groups)
  - K-mer extraction across k values (3, 6, 9, 12)
  - Sequence sizes (100 bp - 100K bp)
  - Minimizer extraction (validates scalar optimal)
  - Spectrum counting (validates parallel makes it slower)
  - Parallel vs scalar comparison (validates 2.2× speedup)
  - Operations comparison at 150bp (Illumina standard)

### Documentation

**Examples**:
- `examples/kmer_operations_full.rs` (150 LOC)
  - Demonstrates all 5 k-mer operations
  - Shows streaming integration
  - Includes Entry 034 evidence summary

**README Updates**:
- Added "K-mer Operations (Evidence-Based)" section
- Updated BERT preprocessing pipeline example
- Documented Python k-mer API with usage examples

**Module Documentation**:
- Comprehensive evidence-based design rationale
- Links to ASBB Entry 034 (GitHub)
- Performance breakdown (hash 50-60%, HashMap 30-40%, validation 5-10%)
- Clear explanation why NEON doesn't help

### Changed

- **Test count**: 254 → 260 tests (+6 property-based tests)
- **Benchmark count**: +6 new benchmark groups for k-mer operations
- **API transparency**: `PARALLEL_THRESHOLD` and `will_use_parallel()` expose behavior

### Evidence Validation

**Entry 034 Findings Confirmed**:
- K-mer operations are data-structure-bound (not compute-bound)
- Hash computation dominates runtime (50-60%, sequential)
- HashMap operations add overhead (30-40%, thread contention)
- Base validation is only 5-10% (NEON-friendly, but minimal impact)
- **Conclusion**: Scalar-only optimal for most operations

**Validates Existing Tools**:
- minimap2's scalar minimizer design (empirically confirmed)
- Identifies 2.2× speedup opportunity for DNABert k-mer extraction

---

## [1.0.0] - 2025-11-05 (Phase 4: Sequence Manipulation Primitives)

### 🎉 First Stable Release

Production-ready ARM-native bioinformatics library with streaming architecture and evidence-based optimization.

**Grade**: A (rust-code-quality-reviewer, post-Phase 4)
**Tests**: 279 passing (209 unit/integration + 70 doc)
**Date**: November 6, 2025

#### Core Sequence Operations (src/operations/sequence.rs)
- **Reverse complement**: DNA/RNA sequence transformation
- **Complement**: Base complementation with lookup table
- **Reverse**: Sequence reversal
- **In-place variants**: Zero-allocation transformations
- **Validation**: `is_valid_dna()`, `is_valid_rna()`, `count_invalid_bases()`
- **IUPAC support**: Full support for ambiguous codes (N, R, Y, S, W, K, M, B, D, H, V)

#### Record-Level Operations (src/operations/record_ops.rs)
- **Extract region**: Subsequence extraction with quality preservation
- **Reverse complement record**: RC with quality reversal
- **Length filtering**: `meets_length_requirement()` with min/max bounds
- **Sequence length**: `sequence_length()` helper

#### Trimming Operations (src/operations/trimming.rs)
- **Fixed-position**: `trim_start()`, `trim_end()`, `trim_both()`
- **Quality-based**: `trim_quality_end()`, `trim_quality_start()`, `trim_quality_both()`
- **Sliding window**: `trim_quality_window()` (Trimmomatic-style)
- **Phred+33 encoding**: Illumina 1.8+ standard

#### Masking Operations (src/operations/masking.rs)
- **Quality-based masking**: Replace low-quality bases with 'N'
- **Length preservation**: Unlike trimming, maintains read length
- **In-place + copy variants**: `mask_low_quality()`, `mask_low_quality_copy()`
- **Counting**: `count_masked_bases()` for QC metrics

#### API Improvements
- **FastqRecord::is_empty()**: Check for empty records after trimming
- **Empty record documentation**: Clear guidance on handling empty results
- **Validation helper**: DRY principle for sequence-quality alignment checks

#### Performance (Scalar Implementation)
Benchmarked on Mac M3 Max, evidence-based NEON deferral decision:
- reverse_complement: 3.7 GiB/s @ 150bp (Illumina standard)
- complement: 3.8 GiB/s @ 150bp
- reverse: 11.1 GiB/s @ 150bp (3× faster, proves table lookup is bottleneck)
- trim_quality_both: Optimized to single-pass (eliminated intermediate allocation)

**NEON Decision**: Deferred (evidence-based)
- Current scalar performance: 3-5 GiB/s (excellent)
- Operations are memory-bound (table lookup bottleneck)
- Estimated NEON speedup: <2× (fails ≥5× threshold)
- Better investment: Additional features, user documentation

#### Testing
- **Unit tests**: 209 passing (up from 121)
- **Property-based tests**: 32 tests with proptest
  - Involutive properties (RC(RC(x)) = x, C(C(x)) = x, R(R(x)) = x)
  - Decomposition (RC(x) = R(C(x)) = C(R(x)))
  - In-place ≡ allocating variants
  - Length preservation
  - Quality threshold validation
  - Window trimming edge cases
- **Doc tests**: 70 passing (comprehensive usage examples)

#### Documentation
- **Module-level docs**: Design principles, evidence citations, RNA limitations
- **Function-level docs**: Every public API with examples
- **Examples**: `examples/sequence_operations.rs` (279 LOC, 5 demo sections)
- **Benchmarks**: `benches/sequence_operations.rs` (249 LOC, 9 benchmark groups)
- **Analysis**: `PHASE4_BENCHMARK_ANALYSIS.md` (comprehensive NEON decision rationale)

#### Code Quality Improvements
Addressed all rust-code-quality-reviewer recommendations:
1. ✅ Optimized `trim_quality_both` (single-pass, no intermediate allocation)
2. ✅ Added `FastqRecord::is_empty()` helper method
3. ✅ Documented empty record handling in all trimming functions
4. ✅ Enhanced RNA documentation (module-level warning + function notes)
5. ✅ Added window trimming edge case property test
6. ✅ Refactored validation logic (DRY principle, 8 instances → 2 helpers)

### Changed

- **Test count**: 121 → 279 tests (158 new tests, +130% increase)
- **Memory efficiency**: `trim_quality_both` now single-pass

### Performance Analysis

**Evidence-Based Decision** (Category 2: Benchmark First):
- Scalar baseline measured: 3.7 GiB/s @ 150bp
- Throughput sufficient: 24.9M reads/second (100M reads in 4 seconds)
- Memory-bound workload: Table lookups cannot parallelize effectively
- NEON threshold not met: Estimated <2× speedup (need ≥5×)
- Documentation: Full analysis in PHASE4_BENCHMARK_ANALYSIS.md

**Comparison to Rule 1 Operations**:
- Base counting: 16.7× NEON speedup (element-wise comparison)
- GC content: 20.3× NEON speedup (element-wise comparison)
- Reverse complement: <2× estimated (table lookup anti-pattern)

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
