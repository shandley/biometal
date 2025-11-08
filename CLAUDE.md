# biometal: Claude Development Guide

**Project**: biometal - ARM-native bioinformatics library
**Latest Release**: v1.2.0 (November 6, 2025)
**Current Work**: BAM/SAM parser implementation (experiments/native-bam-implementation/)

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

Every optimization comes from validated experimental results (apple-silicon-bio-bench):
- Follow OPTIMIZATION_RULES.md: 6 rules from 1,357 experiments (N=30)
- Don't guess: Reference ASBB evidence for optimization decisions
- Document rationale: Link implementations to specific rules/entries

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

Platform priority: Mac ARM → Linux ARM (Graviton) → x86_64 fallback

### 4. Production Quality

- Use `Result<T, BiometalError>` (no panics in library)
- Document every public API with examples
- Property-based testing (proptest)
- Benchmarks (criterion, N=30)
- No `unwrap()` or `expect()` in library code

---

## Project Status

### Released (v1.2.0)
- FASTQ/FASTA streaming parsers (constant memory)
- ARM NEON operations (base counting, GC content, quality filtering)
- Sequence manipulation (reverse_complement, trimming, masking)
- K-mer operations (extraction, minimizers, spectrum)
- Network streaming (HTTP, SRA)
- Python bindings (PyO3 0.27, 40+ functions)
- 347 tests passing (260 library + 87 doc)

### In Progress
- **BAM/SAM Parser**: Native implementation in experiments/native-bam-implementation/
  - Phase 1-3 complete (record parsing, error types, robustness)
  - Currently evaluating for production integration
  - See: experiments/native-bam-implementation/NOODLES_LESSONS.md

### Distribution
- **PyPI**: biometal-rs (pip install biometal-rs)
- **crates.io**: biometal (cargo add biometal)

---

## Project Structure

```
biometal/
├── src/
│   ├── lib.rs              # Public API
│   ├── io/                 # Streaming parsers
│   │   ├── fastq.rs        # FASTQ streaming
│   │   ├── fasta.rs        # FASTA streaming
│   │   ├── bam/            # BAM/SAM (in development)
│   │   ├── compression.rs  # Parallel bgzip + mmap
│   │   ├── network.rs      # HTTP streaming
│   │   ├── paired.rs       # Paired-end reads
│   │   ├── sink.rs         # Output writers
│   │   └── sra.rs          # SRA toolkit wrapper
│   ├── operations/         # Analysis operations
│   │   ├── base_counting.rs    # NEON-optimized
│   │   ├── gc_content.rs       # NEON-optimized
│   │   ├── quality_filter.rs   # NEON-optimized
│   │   ├── sequence.rs         # Sequence manipulation
│   │   ├── record_ops.rs       # Record operations
│   │   ├── trimming.rs         # Quality/fixed trimming
│   │   ├── masking.rs          # Quality masking
│   │   ├── kmer.rs             # K-mer operations
│   │   └── complexity.rs       # Shannon entropy
│   ├── python/             # Python bindings (PyO3)
│   ├── optimization/       # Platform detection
│   ├── error.rs            # Error types
│   └── types.rs            # Common types
├── benches/                # Criterion benchmarks
├── examples/               # Usage examples
├── experiments/            # Research experiments
│   ├── .experiments.toml   # Experiment registry
│   ├── TEMPLATE/           # Experiment template
│   ├── sra-decoder/        # Completed (NO-GO)
│   └── native-bam-implementation/  # In progress
├── OPTIMIZATION_RULES.md   # Evidence base (1,357 experiments)
├── CLAUDE.md               # This file
└── CHANGELOG.md            # Version history
```

---

## The 6 Optimization Rules

### Rule 1: ARM NEON SIMD (16-25x speedup)
**When**: Element-wise operations (complexity 0.30-0.40)
**Evidence**: Entry 020-025 (307 experiments, 9,210 measurements)
**Example**: Base counting, GC content calculation

### Rule 2: Block-Based Processing (10K records)
**When**: Preserving NEON speedup in streaming contexts
**Evidence**: Entry 027 (1,440 measurements)
**Value**: Avoids 82-86% overhead from single-record processing

### Rule 3: Parallel Bgzip (6.5x speedup)
**When**: All bgzip-compressed files
**Evidence**: Entry 029 (CPU parallel prototype)
**Implementation**: Rayon-based parallel block decompression

### Rule 4: Smart mmap (2.5x additional)
**When**: Files ≥50 MB on macOS
**Evidence**: Entry 032 (scale validation, 0.54-544 MB)
**Threshold-based**: Only activate for large files

### Rule 5: Constant-Memory Streaming (~5 MB)
**Always**: Design for streaming, not batch
**Evidence**: Entry 026 (720 measurements, 99.5% reduction)
**Pattern**: Iterator-based APIs, no accumulation

### Rule 6: Network Streaming
**Why**: I/O dominates 264-352x
**Evidence**: Entry 028 (360 measurements)
**Critical**: Makes network streaming viable

See OPTIMIZATION_RULES.md for complete implementation details.

---

## Error Handling

Always use `Result` types with structured errors:

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

Never use `unwrap()` or `expect()` in library code.

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

## Session Guidelines

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

When implementing features:
1. Check OPTIMIZATION_RULES.md for relevant rule
2. Follow the implementation pattern for that rule
3. Link to evidence (lab notebook entry)

When evaluating optimizations:
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

See: experiments/README.md for full process

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

## Quick Reference

### Evidence Base
- 1,357 experiments, 40,710 measurements (N=30)
- Source: apple-silicon-bio-bench
- Rules: 6 optimization rules (OPTIMIZATION_RULES.md)

### Platform Support
1. **Mac ARM** (M1/M2/M3/M4): 16-25x NEON speedup (optimized)
2. **Linux ARM** (Graviton): 6-10x NEON speedup (portable)
3. **x86_64**: 1x scalar fallback (portable)

### Key Files
- `OPTIMIZATION_RULES.md` - 6 evidence-based rules
- `CHANGELOG.md` - Version history
- `README.md` - User documentation
- `experiments/.experiments.toml` - Experiment registry

### Current Focus Areas
- BAM/SAM parser (experiments/native-bam-implementation/)
- Community feedback on v1.2.0 Python bindings
- Performance benchmarking vs existing tools

---

**Last Updated**: November 8, 2025 (Post BAM Phase 3)
