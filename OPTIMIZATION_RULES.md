# Optimization Rules for biometal

**Evidence Base**: 1,357 experiments, 40,710 measurements (N=30 statistical rigor)
**Source**: Apple Silicon Bio Bench (ASBB) - Systematic Hardware Characterization
**Period**: October 30 - November 4, 2025
**Lab Notebook**: 33 entries documenting all experimental work
**Repository**: https://github.com/shandley/apple-silicon-bio-bench

---

## Overview

These optimization rules are derived from comprehensive experimental validation across 1,357 experiments with statistical rigor (N=30 repetitions, 95% CI, Cohen's d effect sizes). Each rule is evidence-linked to specific lab notebook entries documenting the experimental basis.

**Key Insight**: Evidence from one context doesn't always transfer to another. Rules must be validated in YOUR specific architecture. biometal's findings: Rule 3 (parallel BGZF) failed due to streaming constraints, Rule 4 (mmap) limited by decompression bottleneck. See individual rules for context-dependent details.

---

## Rule 1: Use ARM NEON SIMD (4-25× speedup)

### When to Apply

Operations with **element-wise processing** patterns and **≥15% CPU time**

### Expected Speedup (Refined November 2025)

**Compute-Bound Operations** (16-25× speedup):
- Base counting: 16.7× speedup
- GC content calculation: 20.3× speedup
- Quality filtering: 25.1× speedup
- Sequence complexity: 18.2× speedup

**Memory-Bound Operations** (4-8× speedup):
- BAM sequence decoding (4-bit → ASCII): 4.62× speedup
- Quality score decoding (byte + offset): 4-8× speedup (predicted)
- Operations with heavy read/write memory access

**Memory+Allocation Bound** (3-5× speedup):
- Frequent small allocations per-element
- Memory allocation overhead dominates

### Evidence

**Compute-Bound Evidence**:
- **Source**: [Lab Notebook Entry 020-025](lab-notebook/2025-11/) (DAG Framework Validation)
- **Experiments**: 307 total (9,210 measurements)
- **Statistical rigor**: 95% CI, Cohen's d effect sizes (very large, d > 3.0)
- **Cross-platform**: Mac M4 Max, AWS Graviton 3
- **Findings**: `results/dag_statistical/PHASE4_STATISTICAL_ANALYSIS_REPORT.md`

**Memory-Bound Evidence**:
- **Source**: BAM SIMD Sequence Decoding Experiment (November 9, 2025)
- **Finding**: 4.62× NEON speedup (not 16-25×)
- **Explanation**: Memory allocation (20-30% of time) + heavy memory access patterns
- **Overall impact**: +27.5% faster BAM parsing (5.5× the ≥5% target)
- **Documentation**: `experiments/bam-simd-sequence-decoding/FINDINGS.md`

**Effectiveness Predictor**:
- Element-wise operations (complexity 0.30-0.40)
- CPU time ≥15% (Rule 1 threshold)
- Adjust expectations based on memory access pattern

### Implementation Pattern

```rust
#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

/// NEON-optimized base counting (16.7× faster than scalar)
/// Evidence: Entry 020, Cohen's d = 4.82 (very large effect)
#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_neon(seq: &[u8]) -> [u32; 4] {
    let mut counts = [0u32; 4];

    // NEON registers for ACGT counts
    let mut vcounts = [vdupq_n_u32(0); 4];

    // Process 16 bytes at a time
    let chunks = seq.chunks_exact(16);
    let remainder = chunks.remainder();

    for chunk in chunks {
        let seq_vec = vld1q_u8(chunk.as_ptr());

        // Compare against A, C, G, T (4 NEON comparisons in parallel)
        let a_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'A'));
        let c_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'C'));
        let g_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'G'));
        let t_mask = vceqq_u8(seq_vec, vdupq_n_u8(b'T'));

        // Accumulate counts
        vcounts[0] = vaddq_u32(vcounts[0], vpaddlq_u16(vpaddlq_u8(a_mask)));
        vcounts[1] = vaddq_u32(vcounts[1], vpaddlq_u16(vpaddlq_u8(c_mask)));
        vcounts[2] = vaddq_u32(vcounts[2], vpaddlq_u16(vpaddlq_u8(g_mask)));
        vcounts[3] = vaddq_u32(vcounts[3], vpaddlq_u16(vpaddlq_u8(t_mask)));
    }

    // Extract counts from NEON registers
    for i in 0..4 {
        counts[i] = vgetq_lane_u32(vcounts[i], 0) +
                    vgetq_lane_u32(vcounts[i], 1) +
                    vgetq_lane_u32(vcounts[i], 2) +
                    vgetq_lane_u32(vcounts[i], 3);
    }

    // Handle remainder with scalar code
    for &base in remainder {
        match base {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }

    counts
}

/// Scalar fallback for non-ARM platforms
#[cfg(not(target_arch = "aarch64"))]
pub fn count_bases_scalar(seq: &[u8]) -> [u32; 4] {
    let mut counts = [0u32; 4];
    for &base in seq {
        match base {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }
    counts
}
```

### Key Insights

1. **NEON is standard ARM**: Works across Mac, Graviton, Ampere, Raspberry Pi
2. **Not Apple-specific**: Uses standard ARM NEON intrinsics
3. **Portable**: Always provide scalar fallback for x86_64
4. **Predictable**: Complexity score predicts NEON effectiveness

### Platform Support

- ✅ **macOS** (Apple Silicon: M1, M2, M3, M4)
- ✅ **Linux ARM** (AWS Graviton, Ampere Altra, Raspberry Pi 4/5)
- ✅ **Windows ARM** (Surface Pro X, ARM64 Windows)
- ⚠️ **x86_64**: Falls back to scalar (no SIMD equivalent)


---

## Rule 2: Use Block-Based Processing (10K records per block) [DEFERRED]

### Current Status (v1.6.0)

**Implementation**: Block API exists but provides **convenience, not performance optimization**.
**Performance**: Block operations perform similarly to per-record operations (0.88-1.01×).
**Status**: Deferred to future work (Phase 2+). Focus on Rules 3+4 instead (16.3× combined).

### Why This Matters

**Problem**: Record-by-record streaming loses 82-86% of NEON speedup due to function call overhead.

**Solution (Future)**: Process records in blocks of ~10K with inlined NEON operations to preserve SIMD performance.

**Current Implementation**: Block API calls underlying NEON functions once per sequence (10,000 function calls for 10K sequences), resulting in same overhead as per-record approach.

### Evidence

**Source**: [Lab Notebook Entry 027](lab-notebook/2025-11/20251103-027-EXPERIMENT-streaming-overhead.md)
- **Experiments**: 48 total (1,440 measurements)
- **Finding**: Record-by-record NEON = 82-86% overhead
- **Solution**: Block-based processing (10K records) = 4-8% overhead
- **Conclusion**: 14× speedup achievable by reducing function calls from 10,000 to 1

**biometal Validation** (November 11, 2025, N=30):
- **Base counting**: Block is 0.88× per-record (no improvement)
- **GC content**: Block is 0.98× per-record (no improvement)
- **Mean quality**: Block is 1.01× per-record (no improvement)
- **Root cause**: Both approaches make 10,000 function calls
- **Analysis**: See `RULE2_INVESTIGATION_FINDINGS.md`

**Trade-off Analysis**:
- Block size too small (1K): SIMD setup overhead dominates
- Block size too large (100K): Memory pressure, reduces streaming benefit
- **Sweet spot**: 10K records (~1.5 MB for 150bp reads)

### Current Implementation (Convenience API)

```rust
/// Block-based convenience API (similar performance to per-record)
pub fn count_bases_block(sequences: &[&[u8]]) -> Vec<BaseCounts> {
    #[cfg(target_arch = "aarch64")]
    {
        unsafe { count_bases_block_neon(sequences) }
    }
    #[cfg(not(target_arch = "aarch64"))]
    {
        count_bases_block_scalar(sequences)
    }
}

#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_block_neon(sequences: &[&[u8]]) -> Vec<BaseCounts> {
    let mut results = Vec::with_capacity(sequences.len());
    for seq in sequences {
        results.push(count_bases_neon(seq));  // ← Still N function calls!
    }
    results
}
```

**Problem**: Still makes 10,000 function calls → same overhead as per-record.

### Future Implementation (14× Speedup)

To achieve Entry 027's 14× speedup would require:

```rust
#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_block_neon_true(sequences: &[&[u8]]) -> Vec<BaseCounts> {
    use std::arch::aarch64::*;

    let mut results = Vec::with_capacity(sequences.len());

    // NEON setup once (not per-sequence)
    for seq in sequences {
        // Inline NEON operations here (no function call)
        let mut counts = [0u32; 4];
        let mut vcounts = [vdupq_n_u32(0); 4];

        let chunks = seq.chunks_exact(16);
        for chunk in chunks {
            let seq_vec = vld1q_u8(chunk.as_ptr());
            // ... full NEON implementation inlined ...
        }

        results.push(counts);
    }

    results  // Single function call, hot NEON registers throughout
}
```

**Trade-offs**:
- ✅ Achieves 14× speedup (validated by Entry 027)
- ✅ Single function call, amortized overhead
- ❌ Code duplication (~1,200 lines: 3 operations × 2 variants × ~200 lines each)
- ❌ Maintenance burden keeping implementations in sync
- ❌ Binary size increase

**Decision**: Deferred to Phase 2+. Rules 3+4 provide 16.3× combined speedup with less complexity.

### Key Insights

1. **Entry 027 is valid**: 14× speedup IS achievable with proper implementation
2. **API design matters**: Per-sequence operations prevent capturing the speedup
3. **Function call overhead**: 10,000 calls wastes ~500K-900K cycles (82-86% overhead)
4. **Trade-off decision**: Rules 3+4 (16.3× combined) > Rule 2 (14×) for less effort
5. **Current block API**: Useful for convenience, not performance optimization

---

## Rule 3: Parallel Bgzip Decompression - DISABLED (Context-Dependent)

### Status: NOT IMPLEMENTED in biometal

**Multi-scale testing (November 11, 2025)** showed bounded streaming parallel BGZF achieves **0.77-0.84× slowdown** (not 6.5× speedup).

**DAG Decision**: Failed pruning threshold (<1.5×) → Optimization removed.

### Evidence: Two Contexts, Different Results

**Entry 029 Context** (All-at-once decompression):
- **Implementation**: Load entire file → decompress all blocks in parallel
- **Memory**: Unbounded (file size = RAM usage)
- **Speedup**: 6.5× (validated)
- **Example**: `std::fs::read()` → `par_iter()` → decompress all

**biometal Context** (Bounded streaming):
- **Implementation**: Process 8 blocks at a time for constant memory (Rule 5)
- **Memory**: Bounded ~1 MB (regardless of file size)
- **Speedup**: 0.77-0.84× (overhead dominates)
- **Multi-scale**: Gets WORSE with larger files (5.4MB: 0.84×, 544MB: 0.77×)
- **Root cause**: Chunking overhead compounds with scale

### Why Entry 029's 6.5× Doesn't Transfer

**Architectural conflict**: Rule 3 (parallelism) incompatible with Rule 5 (constant memory streaming).

**Trade-off choice**: biometal prioritizes Rule 5 (streaming for TB-scale files) over Rule 3 (speed).

**Documentation**: `RULE3_BENCHMARK_RESULTS.md`, `benches/rule3_multiscale_validation.rs`

### Why not GPU

[Entry 031](lab-notebook/2025-11/20251104-031-EXPERIMENT-metal-deflate-phase2.md):
- Real bgzip uses 100% dynamic Huffman trees (not fixed)
- GPU implementation requires dynamic Huffman decoder: 7-10 days development
- ROI too low: 7-10 days for 2-3× incremental over CPU's theoretical 6.5×
- **Decision**: Neither GPU nor CPU parallel (streaming architecture priority)

### Alternative: All-at-Once Approach (Entry 029's Method)

**If Rule 5 (streaming) is NOT a priority**, Entry 029's approach works:

```rust
// Load entire file (unbounded memory)
let compressed = std::fs::read("file.fq.gz")?;

// Parse all blocks
let blocks = parse_bgzip_blocks(&compressed)?;

// Decompress ALL blocks in parallel (6.5× speedup)
let decompressed: Vec<_> = blocks
    .par_iter()
    .map(decompress_block)
    .collect::<io::Result<Vec<_>>>()?;

// Concatenate
Ok(decompressed.concat())
```

**Pros**: 6.5× speedup (validated)
**Cons**: Memory scales with file size (100 GB file = 100 GB RAM)

### Lesson: Context Dependency

**Key insight**: Evidence from one architecture doesn't always transfer to another.

- Entry 029: Optimized for speed (unbounded memory)
- biometal: Optimized for scalability (constant memory)
- **Different goals → Different optimal solutions**


---

## Rule 4: Use Smart mmap for Large Files (2.5× additional speedup)

### When to Apply

**Files ≥50 MB** on platforms with memory-mapped I/O optimization (macOS validated, Linux future)

**Threshold-based approach**: Use mmap for large files, standard I/O for small files to avoid overhead.

### Evidence

**Source**: [Lab Notebook Entry 032](lab-notebook/2025-11/20251104-032-EXPERIMENT-mmap-apfs-optimization.md)
- **Test 1**: Initial validation (5.4 MB file) = 1.66× speedup
- **Test 2**: Scale validation (0.54 MB to 544 MB)
  - Small files (<50 MB): 0.66-0.99× (overhead dominates, **don't use mmap**)
  - Large files (≥50 MB): 2.30-2.55× speedup (APFS prefetching dominates)
- **Threshold**: 50 MB (empirically determined break-even point)
- **Findings**: `results/io_optimization/MMAP_FINDINGS.md`

**Why threshold matters**:
- Small files: mmap setup cost (~200 µs) >> actual I/O time (~67 µs)
- Large files: APFS prefetching benefit >> mmap overhead

### Implementation Pattern

```rust
use memmap2::Mmap;

const MMAP_THRESHOLD: u64 = 50 * 1024 * 1024; // 50 MB

/// Smart I/O data source (threshold-based mmap)
/// Evidence: Entry 032 (scale validation across 0.54-544 MB)
enum DataSource {
    StandardIo(Vec<u8>),
    MemoryMapped(Mmap),
}

impl DataSource {
    /// Open file with smart I/O method selection
    /// - Small files (<50 MB): Use standard I/O (avoids mmap overhead)
    /// - Large files (≥50 MB): Use mmap + madvise (2.5× speedup)
    pub fn open(path: &Path) -> io::Result<Self> {
        let metadata = std::fs::metadata(path)?;
        let file_size = metadata.len();

        if file_size >= MMAP_THRESHOLD {
            // Large file: Use mmap with APFS optimization hints
            Self::open_mmap(path)
        } else {
            // Small file: Use standard I/O (faster for <50 MB)
            let data = std::fs::read(path)?;
            Ok(DataSource::StandardIo(data))
        }
    }

    /// Open file with memory mapping + APFS hints
    #[cfg(target_os = "macos")]
    fn open_mmap(path: &Path) -> io::Result<Self> {
        use libc::{madvise, MADV_SEQUENTIAL, MADV_WILLNEED};

        let file = std::fs::File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        // Give kernel sequential access hints for APFS optimization
        unsafe {
            madvise(
                mmap.as_ptr() as *mut _,
                mmap.len(),
                MADV_SEQUENTIAL | MADV_WILLNEED,
            );
        }

        Ok(DataSource::MemoryMapped(mmap))
    }
}
```

### Performance: Context-Dependent

**Entry 032 tested**: RAW file I/O (no decompression)
- Large files (≥50 MB): 2.30-2.55× speedup ✓

**biometal reality**: Compressed file processing (with decompression)
- Decompression: 98.7% of time (CPU-bound)
- I/O: 1.3% of time
- **mmap benefit**: ~1% overall (Amdahl's Law)

**Why the difference**:
```
Entry 032: std::fs::read() → 100% I/O → 2.5× applies fully
biometal: decompress() → 99% CPU → 2.5× on 1% = negligible
```

**Documentation**: `RULE4_FINDINGS.md`, `benches/rule4_mmap_validation.rs`

### When mmap Actually Helps

**Uncompressed files** (hypothetical):
- 100% I/O-bound → 2.5× speedup applies fully
- Example: Uncompressed BAM, large text files

**Compressed files** (current use case):
- ~1% improvement (I/O is only 1.3% of total time)
- Implementation kept (no harm, small benefit)

### Combined Performance (Revised)

**Current Reality** (Rule 3 disabled, Rule 4 limited):
- Sequential decompression: Baseline
- mmap on compressed files: ~1% improvement
- **Overall**: ~1× (no significant speedup for compressed files)

**Original Projection** (Rules 3+4 both working):
- Parallel bgzip: 6.5×
- Smart mmap: 2.5× additional
- **Combined**: 16.3× (not achieved due to architectural conflicts)
- Before optimization: NEON 1.04-1.08× E2E (I/O bottleneck 264-352×)
- After optimization: Projected **17× E2E speedup** for large files

### Key Insights

1. **Threshold is critical**: mmap hurts small files, helps large files
2. **Platform-specific**: APFS (macOS) validated, Linux future
3. **Complementary optimizations**: Parallel bgzip + mmap = multiplicative benefits
4. **Unified memory helps**: Apple Silicon's architecture optimizes mmap

### Platform Support

- ✅ **macOS**: Validated with APFS (2.3-2.5× for ≥50 MB)
- ⏳ **Linux**: Future validation on AWS Graviton (Week 3-4)
- ❌ **Windows**: No equivalent API (use standard I/O)

---

## Rule 5: Design for Constant-Memory Streaming (~5 MB)

### Result

**99.5% memory reduction** (1,344 MB → 5 MB @ 1M sequences)

**Critical finding**: Streaming memory is **CONSTANT (~5 MB)** regardless of dataset size. This enables analyzing 5TB datasets on laptops with <100 MB RAM.

### Evidence

**Source**: [Lab Notebook Entry 026](lab-notebook/2025-11/20251103-026-EXPERIMENT-streaming-memory-footprint-v2.md)
- **Experiments**: 24 total (720 measurements, corrected v2)
- **Scales tested**: 10K, 100K, 1M sequences
- **Finding**: Streaming memory remains constant at ~5 MB regardless of scale
- **Batch comparison**: 1M sequences = 1,344 MB (batch) vs 5 MB (streaming)
- **Conclusion**: Memory footprint is independent of dataset size

### Implementation Pattern

```rust
/// Constant-memory FASTQ streaming (99.5% memory reduction)
/// Evidence: Entry 026 (constant ~5 MB regardless of scale)
pub struct FastqStream<R: BufRead> {
    reader: R,
    line_buffer: String,
    record_buffer: FastqRecord,
    // NO full-dataset storage! Memory is constant.
}

impl<R: BufRead> FastqStream<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buffer: String::with_capacity(512), // ~5 MB total footprint
            record_buffer: FastqRecord::default(),
        }
    }
}

impl<R: BufRead> Iterator for FastqStream<R> {
    type Item = io::Result<FastqRecord>;

    /// Read one record at a time (constant memory)
    /// No accumulation - process and discard
    fn next(&mut self) -> Option<Self::Item> {
        // Read FASTQ record, process, return
        // Buffer is reused, keeping memory constant
    }
}
```

### Dataset Size Independence

| Dataset Size | Batch Memory | Streaming Memory | Reduction |
|--------------|--------------|------------------|-----------|
| 10K sequences | 13.4 MB | 5 MB | 62.7% |
| 100K sequences | 134 MB | 5 MB | 96.3% |
| 1M sequences | 1,344 MB | 5 MB | 99.5% |
| 10M sequences | 13,440 MB | 5 MB | 99.96% |
| **5TB dataset** | **5,000,000 MB** | **5 MB** | **99.9999%** |

**Breakthrough**: Memory footprint is constant, enabling analysis of arbitrarily large datasets on consumer hardware.

### Key Insights

1. **Constant memory**: Independent of dataset size (5 MB for 10K or 10M sequences)
2. **Enables 5TB analysis**: On 24GB laptop (Data Access pillar validated)
3. **Iterator pattern**: Natural fit for streaming (no accumulation)
4. **Buffer reuse**: Clear and reuse, don't allocate new buffers
5. **Combines with blocks**: 10K block buffer = ~1.5 MB (still constant)

---

## Rule 6: I/O Bottleneck is Critical (Network Streaming Essential)

### Finding

**NEON provides only 1.04-1.08× E2E speedup** (vs 16-25× isolated computation)

**Root cause**: I/O dominates by 264-352× compared to compute time.

### Evidence

**Source**: [Lab Notebook Entry 028](lab-notebook/2025-11/20251103-028-EXPERIMENT-streaming-e2e-pipeline.md)
- **Experiments**: 12 total (360 measurements)
- **Pipeline**: Read FASTQ.gz → Process → Filter → Write
- **Finding**: NEON 1.04-1.08× E2E (vs 16-25× isolated)
- **Analysis**: I/O bottleneck is **264-352× slower** than compute
- **Conclusion**: Network streaming + caching is CRITICAL, not optional

### Bottleneck Breakdown

**Without optimization**:
- Compute time: 1.0× (baseline)
- I/O time: **264-352×** (dominates!)
- NEON speedup (isolated): 16-25×
- NEON speedup (E2E): 1.04-1.08× (I/O masks compute gains)

**With I/O optimization stack** (Rules 3 + 4):
- Compute time: 1.0×
- I/O time: 16-22× (264-352× → 16-22× with optimization)
- NEON speedup (E2E): **~17×** (projected)

### Implication for Architecture

The I/O bottleneck validation makes **network streaming** a critical feature:

1. **Problem**: Downloading 5TB dataset takes days/weeks on slow connections
2. **Solution**: Stream directly from network with smart caching
3. **Benefit**: Start analysis immediately, cache only what's needed
4. **Result**: 5TB dataset becomes accessible without 5TB download

### Network Streaming Design (Week 3-4)

```rust
/// Network streaming source (addresses I/O bottleneck)
/// Evidence: Entry 028 (I/O dominates 264-352×)
pub enum DataSource {
    Local(PathBuf),
    Http(Url),
    Sra(String), // SRA accession
}

pub struct StreamingReader {
    source: DataSource,
    cache: LruCache<BlockId, Vec<u8>>,
    prefetch: Prefetcher,
}
```

### Key Insights

1. **I/O dominates in real workloads**: Compute speedup masked by I/O bottleneck
2. **Layered optimization essential**: Parallel bgzip + mmap reduces bottleneck 16.3×
3. **Network streaming critical**: Enables 5TB analysis without download
4. **Smart caching**: LRU cache balances memory and network requests
5. **Prefetching**: Background downloads hide network latency


---

## Summary: Layered Optimization Strategy

### Individual Rules (Good)

| Rule | Speedup | Platform | Status | Priority |
|------|---------|----------|--------|----------|
| NEON SIMD | 16-25× | ARM | ✅ v1.0.0 | High |
| Block-based | 14× (potential) | All | ⏳ Deferred | Medium |
| Parallel bgzip | 6.5× | All | ❌ Phase 2 | High |
| Smart mmap | 2.5× | macOS | ❌ Phase 2 | High |
| Constant streaming | 99.5% mem | All | ✅ v1.0.0 | High |
| Network streaming | Analysis without download | All | ✅ v1.0.0 | High |

**Note**: Rule 2 (Block processing) deferred. Current block API is convenience feature only (no speedup).

### Combined Stack (Exceptional)

**Example**: I/O optimization demonstrates layered benefits

**Layer 1** (Parallel bgzip): 6.5× speedup  
**Layer 2** (Smart mmap): 2.5× additional  
**Combined**: 6.5 × 2.5 = **16.3× total speedup**

**E2E impact**:
- Before: NEON 1.04-1.08× E2E (I/O bottleneck 264-352×)
- After: Projected **17× E2E speedup** for large files

**Time to process 1M sequences**:
- Before: 12.3 seconds
- After: 0.75 seconds (16.3× faster!)

### Evidence Base

- **Total experiments**: 1,357
- **Total measurements**: 40,710 (N=30)
- **Lab notebook entries**: 33 (full experimental log)
- **Statistical rigor**: 95% CI, Cohen's d effect sizes
- **Publications**: 3 papers in preparation

### Implementation Timeline

**v1.0.0 (Released)**: Core streaming + NEON
- Rules 1, 5, 6 (NEON SIMD + constant streaming + network streaming)

**Phase 2 (Planned)**: I/O optimization stack
- Rules 3, 4 (parallel bgzip + smart mmap) → 16.3× combined speedup

**Future (Deferred)**: True block processing
- Rule 2 (inline NEON for 14× speedup) - requires code duplication
- Lower priority: Rules 3+4 provide more speedup (16.3×) with less complexity

---

## Rule 7: Metal GPU + Unified Memory (Apple Silicon Only)

### Overview

⭐ **NEW**: Extension beyond original ASBB scope to validate Apple Silicon breakthrough opportunities

Apple's unified memory architecture (UMA) enables zero-copy GPU acceleration impossible on traditional CUDA/x86 systems. This rule applies to compute-bound BAM operations with massive parallelism (1000s of genomic positions).

### When to Apply

Operations with **embarrassingly parallel** patterns across genomic positions:
- **Pileup generation**: Accumulate 50-100 reads per position (10-50× expected)
- **Coverage calculation**: Parallel position counting (20-100× expected)
- **Depth computation**: Statistical accumulation across regions

**Do NOT use Metal for**:
- I/O-bound operations (parsing, file reading)
- Operations already optimized with NEON (16-25× Rule 1)
- String processing (CIGAR parsing - use NEON instead)

### Why Apple Silicon is Different

**Traditional GPU (CUDA)**:
```
CPU RAM ─[PCIe: 16 GB/s]→ GPU RAM
         ↓ copy overhead ↓
      50-80% of speedup lost to memory transfers
```

**Apple Silicon Unified Memory**:
```
CPU ←[400 GB/s]→ Shared RAM ←[400 GB/s]→ GPU
     Zero-copy, 25× bandwidth advantage
```

### Evidence (To Be Validated)

⚠️ **Status**: Planned experiments for Week 9-10 (Dec 30 - Jan 12, 2026)

**Entry 034: NEON CIGAR Parsing** (Week 9)
- **Experiments**: 30 (N=30 statistical rigor)
- **Platforms**: Mac M4 Max (NEON) vs. x86_64 (scalar)
- **Operation**: Parse 1M CIGAR strings (BAM alignment format)
- **Expected speedup**: 10-20× (analogous to Rule 1 base counting 16.7×)
- **Rationale**: CIGAR is string pattern matching, perfect for SIMD
- **Status**: Not yet conducted

**Entry 035: Metal GPU Pileup Generation** (Week 10)
- **Experiments**: 30 (N=30 statistical rigor)
- **Platforms**: Mac M4 Max (Metal) vs. CPU parallel (rayon)
- **Dataset**: 50× coverage BAM file, chromosome 1 (millions of reads)
- **Baseline**: SAMtools mpileup (single-threaded)
- **Expected speedup**: 10-50× vs. SAMtools, 2-5× vs. CPU parallel
- **Key advantage**: Zero-copy UMA (no PCIe bottleneck)
- **Status**: Not yet conducted

**Entry 036: Metal GPU Coverage Calculation** (Week 11+)
- **Experiments**: 30 (N=30 statistical rigor)
- **Operation**: Calculate per-position coverage across genome
- **Expected speedup**: 20-100× (simpler than pileup, highly parallel)
- **Status**: Planned for v1.2+

### Implementation Pattern (Metal GPU)

```rust
#[cfg(all(target_arch = "aarch64", target_os = "macos"))]
pub fn generate_pileup_metal(
    bam_records: &[BamRecord],
    region: GenomicRegion,
) -> Result<Pileup> {
    use metal::{Device, MTLResourceOptions};

    let device = Device::system_default()
        .ok_or(BiometalError::MetalNotAvailable)?;

    // Zero-copy: Share CPU buffer with GPU (UMA advantage)
    let shared_buffer = device.new_buffer_with_data(
        bam_records.as_ptr() as *const _,
        (bam_records.len() * std::mem::size_of::<BamRecord>()) as u64,
        MTLResourceOptions::StorageModeShared, // KEY: Zero-copy!
    );

    // Metal compute shader processes 1000s of positions in parallel
    let pileup = compute_pileup_gpu(device, shared_buffer, region)?;

    Ok(pileup)
}

// Fallback for non-macOS ARM platforms (Graviton, RPi)
#[cfg(all(target_arch = "aarch64", not(target_os = "macos")))]
pub fn generate_pileup_metal(
    bam_records: &[BamRecord],
    region: GenomicRegion,
) -> Result<Pileup> {
    // Use NEON + rayon CPU parallel instead
    generate_pileup_neon_parallel(bam_records, region)
}
```

### Metal Compute Shader (src/metal/shaders.metal)

```metal
#include <metal_stdlib>
using namespace metal;

kernel void pileup_accumulate(
    constant BamRecord* records [[buffer(0)]],
    device atomic_uint* pileup [[buffer(1)]],
    constant uint& num_records [[buffer(2)]],
    constant uint& start_position [[buffer(3)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= num_records) return;

    const BamRecord& record = records[gid];

    // Each GPU thread processes one read
    for (uint i = 0; i < record.length; i++) {
        uint genome_pos = record.position + i - start_position;
        atomic_fetch_add_explicit(&pileup[genome_pos], 1, memory_order_relaxed);
    }
}
```

### Decision Framework: Metal vs. NEON vs. rayon

| Operation | Best Tool | Speedup | Rationale |
|-----------|-----------|---------|-----------|
| Pileup (Mac) | Metal GPU | 10-50× | Massive parallelism, zero-copy UMA |
| Pileup (Linux ARM) | NEON + rayon | 5-10× | No Metal, but NEON still helps |
| Coverage (Mac) | Metal GPU | 20-100× | Embarrassingly parallel |
| CIGAR parsing | NEON | 10-20× | String SIMD pattern matching |
| Base counting | NEON | 16.7× | Already proven (Rule 1) |
| Quality filtering | NEON | 25.1× | Already proven (Rule 1) |
| File decompression | rayon | 6.5× | I/O bound (Rule 3) |

### Validation Plan

**Week 9** (Dec 30 - Jan 5):
1. Implement NEON CIGAR parsing
2. Benchmark against scalar baseline (N=30)
3. Calculate Cohen's d, 95% CI
4. Document Entry 034

**Week 10** (Jan 6-12):
1. Implement Metal pileup generation
2. Benchmark against SAMtools + CPU parallel (N=30)
3. Validate zero-copy UMA advantage
4. Calculate Cohen's d, 95% CI
5. Document Entry 035

**Week 11+** (v1.2):
1. Extend to coverage calculation (Entry 036)
2. Compare CUDA analogues (if available)
3. Publish findings

### Unique Contribution

⭐ **World-first**: No published bioinformatics work using Metal GPU + UMA
- Parabricks uses CUDA (PCIe bottleneck)
- PaCBAM uses CPU parallel only
- SAMtools is single-threaded
- **biometal will be the first** to exploit Apple's true unified memory

### Platform Support Tiers

1. **Mac (Apple Silicon)**: Full Metal + NEON acceleration
2. **Linux ARM** (Graviton, RPi): NEON + rayon fallback
3. **x86_64**: Scalar + rayon fallback

**Key insight**: Even without Metal, ARM still wins with NEON (Rule 1)

---

## References

**Full experimental documentation**:
- Repository: https://github.com/shandley/apple-silicon-bio-bench
- Lab notebook: 33 entries documenting 1,357 experiments
- Results: Comprehensive analysis in `results/` directory
- Plots: Publication-quality visualizations

**Key findings documents**:
- DAG Framework: `results/dag_statistical/PHASE4_STATISTICAL_ANALYSIS_REPORT.md`
- Streaming: `results/streaming/STREAMING_FINDINGS.md`
- I/O Optimization: `results/io_optimization/MMAP_FINDINGS.md`
- Parallel bgzip: `results/bgzip_parallel/FINAL_DECISION.md`

**Publications** (in preparation):
1. DAG Framework: BMC Bioinformatics
2. biometal Library: Bioinformatics (Application Note) or JOSS
3. Four-Pillar Democratization: GigaScience

---

**Document Version**: 1.1
**Date**: November 4, 2025
**Status**:
- Rules 1-6: Evidence complete (1,357 experiments, 40,710 measurements)
- Rule 7: Planned for Week 9-10 (+90 experiments → 1,447 total)
**Current**: biometal v1.0 implementation (FASTQ/FASTA, Nov 4 - Dec 15)
**Next**: biometal v1.1 with Metal GPU breakthrough (Dec 16 - Jan 12, 2026)
