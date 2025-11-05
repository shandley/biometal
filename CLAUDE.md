# biometal: Claude Development Guide

**Project**: biometal - ARM-native bioinformatics library
**Status**: v1.0.0 - Production Release (November 5, 2025)
**Grade**: A+ (rust-code-quality-reviewer)

---

## Mission

Democratize bioinformatics by enabling 5TB dataset analysis on consumer hardware through:
- **Streaming architecture**: Constant ~5 MB memory (not load-all)
- **ARM-native performance**: 16-25x NEON speedup
- **Network streaming**: Analyze without downloading
- **Evidence-based optimization**: Every rule validated experimentally

**Target**: LMIC researchers, small labs, students, field researchers, ML practitioners

---

## Core Principles

### 1. Evidence-Based Design

Every optimization comes from validated experimental results (ASBB project):
- **Follow OPTIMIZATION_RULES.md**: 6 rules from 1,357 experiments (N=30)
- **Don't guess**: Reference ASBB evidence for optimization decisions
- **Document rationale**: Link implementations to specific rules/entries

Example:
```rust
// Rule 2: Block size from Entry 027 (1,440 measurements)
const BLOCK_SIZE: usize = 10_000; // Evidence-based, not arbitrary
```

### 2. Streaming-First Architecture

Always design for constant memory:
- Bad: `Vec<FastqRecord>` (accumulates in memory)
- Good: `Iterator<Item = FastqRecord>` (constant memory)
- Target: ~5 MB regardless of dataset size (Rule 5)

### 3. ARM-Native with Portable Fallback

Always provide both ARM and fallback implementations:
```rust
#[cfg(target_arch = "aarch64")]
pub fn operation_neon(input: &[u8]) -> Result { /* 16-25x faster */ }

#[cfg(not(target_arch = "aarch64"))]
pub fn operation_scalar(input: &[u8]) -> Result { /* x86_64 fallback */ }

pub fn operation(input: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { operation_neon(input) }
    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(input) }
}
```

Platform priority: Mac → Linux ARM (Graviton) → x86_64 fallback

### 4. Production Quality

- Use `Result<T, BiometalError>` (no panics in library)
- Document every public API with examples
- Property-based testing (proptest)
- Benchmarks (criterion, N=30)
- No `unwrap()` or `expect()` in library code

---

## Project Status

**v1.0.0 Released** (November 5, 2025)
- Production-ready ARM-native bioinformatics library
- 121 tests passing (87 unit + 7 integration + 27 doc)
- Python bindings (PyO3 0.27, Python 3.9-3.14)
- Cross-platform validated (Mac ARM, AWS Graviton, x86_64)
- Grade A+ (rust-code-quality-reviewer)

See [CHANGELOG.md](CHANGELOG.md) for complete release history.

### Recent Experiments

**sra-decoder** (Nov 5, 2025) - Status: NO-GO
- Hypothesis: Native ARM SRA decoder achieves ≥10x speedup
- Outcome: NO-GO (Day 2, evidence-based termination)
- Finding: VDB format complexity (5,000-10,000 LOC), projected 2-3x speedup (fails ≥10x threshold)
- Alternative: SRA Toolkit wrapper (500-1,000 LOC, NCBI maintains format)
- Documentation: experiments/sra-decoder/{PROPOSAL,RESEARCH_LOG,FINDINGS}.md
- Learning: Time-boxed experiments catch dead-ends early (saved 12 days)

See: `experiments/.experiments.toml` and `experiments/sra-decoder/`

---

## Project Structure

```
biometal/
├── src/
│   ├── lib.rs              # Public API, re-exports
│   ├── io/                 # Streaming parsers (Rules 3-5)
│   │   ├── fastq.rs        # FASTQ streaming (DONE)
│   │   ├── fasta.rs        # FASTA streaming (DONE)
│   │   ├── compression.rs  # Parallel bgzip + mmap (DONE)
│   │   ├── network.rs      # HTTP streaming (DONE)
│   │   └── sra.rs          # SRA (toolkit wrapper recommended)
│   ├── operations/         # NEON-optimized operations (Rule 1)
│   │   ├── base_counting.rs    # DONE
│   │   ├── gc_content.rs       # DONE
│   │   └── quality_filter.rs   # DONE
│   ├── optimization/       # Platform detection
│   │   ├── platform.rs     # DONE
│   │   └── thresholds.rs   # DONE
│   ├── error.rs            # Error types (DONE)
│   └── types.rs            # Common types
├── benches/                # Criterion benchmarks
├── examples/               # Usage examples
├── docs/                   # Documentation
│   ├── ARCHITECTURE.md           # Network streaming architecture
│   ├── PERFORMANCE_TUNING.md     # Configuration guide
│   └── CODE_QUALITY_IMPROVEMENTS.md  # Pending improvements
├── experiments/            # Research experiments
│   ├── .experiments.toml   # Experiment registry
│   ├── TEMPLATE/           # Experiment template
│   └── sra-decoder/        # Completed experiment (NO-GO)
├── OPTIMIZATION_RULES.md   # Evidence base (1,357 experiments)
├── README.md               # User documentation
├── CLAUDE.md               # This file
└── Cargo.toml
```

---

## Implementation Guidelines

### Rule 1: ARM NEON SIMD (16-25x speedup)

When: Element-wise operations (complexity 0.30-0.40)

```rust
#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

#[cfg(target_arch = "aarch64")]
pub unsafe fn count_bases_neon(seq: &[u8]) -> [u32; 4] {
    // Process 16 bytes at a time with NEON
    // See OPTIMIZATION_RULES.md Rule 1 for full example
}

#[cfg(not(target_arch = "aarch64"))]
pub fn count_bases_scalar(seq: &[u8]) -> [u32; 4] {
    // Scalar fallback
}
```

Evidence: Entry 020-025 (307 experiments, 9,210 measurements)

### Rule 2: Block-Based Processing (10K records)

Why: Preserves NEON speedup (avoids 82-86% overhead)

```rust
const BLOCK_SIZE: usize = 10_000; // From Entry 027

pub struct FastqStream<R: BufRead> {
    reader: R,
    block_buffer: Vec<FastqRecord>,
}

impl<R: BufRead> FastqStream<R> {
    fn process_block(&mut self) -> Result<ProcessedBlock> {
        self.block_buffer.clear();

        while self.block_buffer.len() < BLOCK_SIZE {
            match self.read_record()? {
                Some(record) => self.block_buffer.push(record),
                None => break,
            }
        }

        let results = unsafe { process_block_neon(&self.block_buffer) };
        Ok(ProcessedBlock::new(results))
    }
}
```

Evidence: Entry 027 (1,440 measurements)

### Rule 3: Parallel Bgzip (6.5x speedup)

When: All bgzip-compressed files

```rust
use rayon::prelude::*;

pub fn decompress_bgzip_parallel(compressed: &[u8]) -> io::Result<Vec<u8>> {
    let blocks = parse_bgzip_blocks(compressed)?;

    let decompressed: Vec<_> = blocks
        .par_iter()
        .map(|block| decompress_block(block))
        .collect::<io::Result<Vec<_>>>()?;

    Ok(decompressed.concat())
}
```

Evidence: Entry 029 (CPU parallel prototype)

### Rule 4: Smart mmap (2.5x additional, threshold-based)

When: Files ≥50 MB on macOS

```rust
const MMAP_THRESHOLD: u64 = 50 * 1024 * 1024; // 50 MB

enum DataSource {
    StandardIo(Vec<u8>),
    MemoryMapped(Mmap),
}

impl DataSource {
    pub fn open(path: &Path) -> io::Result<Self> {
        let size = std::fs::metadata(path)?.len();

        if size >= MMAP_THRESHOLD {
            Self::open_mmap(path)
        } else {
            Ok(Self::StandardIo(std::fs::read(path)?))
        }
    }
}
```

Evidence: Entry 032 (scale validation, 0.54-544 MB)

### Rule 5: Constant-Memory Streaming (~5 MB)

Always: Design for streaming, not batch

```rust
pub struct FastqStream<R: BufRead> {
    reader: R,
    line_buffer: String,
    // NO Vec<FastqRecord> accumulation
}

impl<R: BufRead> Iterator for FastqStream<R> {
    type Item = io::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        // Read one record, return, discard
        // Memory stays constant
    }
}
```

Evidence: Entry 026 (720 measurements, 99.5% reduction)

### Rule 6: Network Streaming

Why: I/O dominates 264-352x (makes network streaming critical)

```rust
pub enum DataSource {
    Local(PathBuf),
    Http(Url),
    Sra(String),
}

pub struct StreamingReader {
    source: DataSource,
    cache: LruCache<BlockId, Vec<u8>>,
    prefetch: Prefetcher,
}
```

Evidence: Entry 028 (360 measurements)

---

## Error Handling

Always use `Result` types:

```rust
#[derive(Debug, thiserror::Error)]
pub enum BiometalError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid FASTQ format at line {line}: {msg}")]
    InvalidFormat { line: usize, msg: String },

    #[error("Network error: {0}")]
    Network(String),

    #[error("Compression error: {0}")]
    Compression(String),
}

pub type Result<T> = std::result::Result<T, BiometalError>;
```

---

## Testing Strategy

### Property-Based Testing

```rust
use proptest::prelude::*;

proptest! {
    #[test]
    fn test_base_counting_matches_naive(seq in "[ACGT]{1,1000}") {
        let neon_result = count_bases_neon(seq.as_bytes());
        let naive_result = count_bases_naive(seq.as_bytes());
        prop_assert_eq!(neon_result, naive_result);
    }
}
```

### Benchmarking

```rust
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_base_counting(c: &mut Criterion) {
    let seq = generate_sequence(100_000);

    c.bench_function("base_counting_neon", |b| {
        b.iter(|| count_bases_neon(&seq))
    });
}

criterion_group!(benches, bench_base_counting);
criterion_main!(benches);
```

---

## Documentation Requirements

Every public API must have:
1. Doc comment explaining functionality
2. Example showing usage
3. Link to evidence (when implementing optimizations)

Example:
```rust
/// Stream FASTQ records from a file with constant memory.
///
/// Uses block-based processing (10K records) to preserve ARM NEON speedup
/// while maintaining streaming benefits. Memory footprint remains constant
/// at ~5 MB regardless of file size.
///
/// # Evidence
///
/// - Rule 2 (Block-based): Entry 027, 1,440 measurements
/// - Rule 5 (Streaming): Entry 026, 99.5% memory reduction
///
/// # Example
///
/// ```
/// use biometal::FastqStream;
///
/// let stream = FastqStream::from_path("large.fq.gz")?;
/// for record in stream {
///     let record = record?;
///     // Process one record at a time
/// }
/// # Ok::<(), biometal::Error>(())
/// ```
pub struct FastqStream<R: BufRead> { /* ... */ }
```

---

## Common Pitfalls

### 1. Accumulating Records in Memory

Bad:
```rust
let mut records = Vec::new();
for line in reader.lines() {
    records.push(parse_record(line)?);
}
```

Good:
```rust
for record in FastqStream::from_path(path)? {
    let record = record?;
    // Process immediately, no accumulation
}
```

### 2. Using unwrap() in Library Code

Bad:
```rust
pub fn operation(input: &[u8]) -> Output {
    parse(input).unwrap() // Panics on invalid input
}
```

Good:
```rust
pub fn operation(input: &[u8]) -> Result<Output> {
    parse(input).map_err(|e| BiometalError::InvalidFormat {
        line: 0,
        msg: e.to_string(),
    })
}
```

### 3. Implementing Optimizations Without Evidence

Bad:
```rust
const BLOCK_SIZE: usize = 8_192; // Arbitrary choice
```

Good:
```rust
// Rule 2: Block size from Entry 027 (1,440 measurements)
const BLOCK_SIZE: usize = 10_000;
```

### 4. Platform-Specific Code Without Fallback

Bad:
```rust
pub fn operation(input: &[u8]) -> Result {
    unsafe { operation_neon(input) } // Only works on ARM
}
```

Good:
```rust
pub fn operation(input: &[u8]) -> Result {
    #[cfg(target_arch = "aarch64")]
    { unsafe { operation_neon(input) } }

    #[cfg(not(target_arch = "aarch64"))]
    { operation_scalar(input) }
}
```

---

## For Claude: Session Guidelines

### What to Emphasize

- Evidence-based design (follow OPTIMIZATION_RULES.md)
- Streaming-first architecture (constant memory)
- ARM-native with portable fallback
- Production quality (error handling, docs, tests)

### What NOT to Do

- Don't make up optimization parameters (refer to evidence)
- Don't accumulate records in memory (streaming only)
- Don't panic in library code (use Result)
- Don't implement ARM-only code without scalar fallback

### Decision Framework

When user asks "how should I implement X?":
1. Check OPTIMIZATION_RULES.md for relevant rule
2. Follow the implementation pattern for that rule
3. Link to evidence (lab notebook entry)

When user proposes optimization:
1. Is this validated in ASBB experiments?
2. If yes: Which rule/entry documents it?
3. If no: Suggest validating first or using proven approach

### Experiments Management

When research/innovation ideas arise:
1. Use `experiments/` directory for time-boxed validation
2. Follow TEMPLATE structure (PROPOSAL, RESEARCH_LOG, FINDINGS)
3. Clear go/no-go criteria (quantitative thresholds)
4. Time-box research (typically 2 weeks max)
5. Document negative results (valuable for community)
6. Update `.experiments.toml` registry

See: `experiments/README.md` for full process

---

## Quick Reference

### Evidence Base
- 1,357 experiments, 40,710 measurements (N=30)
- Source: apple-silicon-bio-bench
- Rules: 6 optimization rules (OPTIMIZATION_RULES.md)

### Current Status
- **Version**: v1.0.0 (Released November 5, 2025)
- **Tests**: 121 passing (87 unit + 7 integration + 27 doc)
- **Grade**: A+ (rust-code-quality-reviewer)
- **Python**: PyO3 0.27 (Python 3.9-3.14)
- **Platforms**: Mac ARM (optimized), Graviton/x86_64 (portable)

### Platform Support
1. **Mac ARM** (M1/M2/M3/M4): 16-25× NEON speedup (optimized)
2. **Linux ARM** (Graviton): 6-10× NEON speedup (portable)
3. **x86_64**: 1× scalar fallback (portable)

---

**Last Updated**: November 5, 2025 (v1.0.0 Release)
