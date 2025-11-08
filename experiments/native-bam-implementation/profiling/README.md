# BAM Profiling Tool

**Purpose**: Measure noodles BAM parsing performance on ARM to validate NEON optimization hypothesis.

**Critical Question**: Is sequence decoding a significant bottleneck? (Need ≥15% CPU time for GO decision)

---

## Quick Start

### 1. Build the Profiling Tool

```bash
cd profiling
cargo build --release
```

### 2. Get Test Data

Option A: Download public BAM file
```bash
# Small test file from 1000 Genomes (chromosome 22, ~100MB)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom22.ILLUMINA.bwa.GBR.low_coverage.20120522.bam

# Or use the helper script
./get_test_data.sh
```

Option B: Generate synthetic BAM
```bash
# Install samtools if needed: brew install samtools
./generate_test_bam.sh 10000    # 10K records
./generate_test_bam.sh 100000   # 100K records
```

### 3. Profile with Instruments (macOS)

```bash
# Build with debug symbols
cargo build --release

# Run with Instruments Time Profiler
instruments -t "Time Profiler" \
    ./target/release/bam-profiling \
    ../test-data/HG00096.chrom22.ILLUMINA.bwa.GBR.low_coverage.20120522.bam

# Open the trace file
open instrumentscli*.trace
```

### 4. Profile with cargo-flamegraph

```bash
# Install flamegraph
cargo install flamegraph

# Generate flamegraph
cargo flamegraph --bin bam-profiling -- ../test-data/test.bam

# Open flamegraph.svg in browser
open flamegraph.svg
```

---

## What to Look For

### Primary NEON Opportunity: Sequence Decoding

Look for CPU time in:
- `decode_base` (noodles-bam/src/record/codec/decoder/sequence.rs)
- `read_sequence` (4-bit unpacking logic)
- Any sequence-related functions

### GO/NO-GO Criteria

| CPU % in Sequence Decoding | Decision | Confidence |
|----------------------------|----------|------------|
| ≥30% | **STRONG GO** | 3-4× overall speedup likely |
| 15-29% | **GO** | 2× overall speedup likely |
| 10-14% | **BORDERLINE** | Re-evaluate with BGZF |
| <10% | **NO-GO** | Unlikely to reach 2× threshold |

### Other Hotspots

Also document:
- **BGZF decompression**: Expected ~30-40% (we have 6.5× parallel solution)
- **I/O operations**: File reading, buffering
- **CIGAR parsing**: TAG parsing, record construction
- **Quality score decoding**: ASCII to Phred conversion

---

## Example Output

```
Profiling BAM parsing with noodles
File: ../test-data/HG00096.chrom22.bam

Header parsing: 1.2ms
Reference sequences parsing: 0.8ms

Parsed 1000000 records (87245.3 rec/s)...

Parsing complete!

Results:
  Records parsed: 1000000
  Total sequence length: 100000000 bp
  Total quality scores: 100000000
  Elapsed time: 11.456 s

Throughput:
  87304.2 records/sec
  0.09 Mrec/sec
  8.7 Mbp/sec

Profiling Notes:
  - Run with: instruments -t 'Time Profiler' ./target/release/bam-profiling <file.bam>
  - Or use: cargo flamegraph --bin bam-profiling -- <file.bam>
  - Look for CPU time in:
    * decode_base (sequence decoding)
    * read_sequence (4-bit unpacking)
    * BGZF decompression
  - GO Criteria: Sequence decoding ≥15% CPU time
  - NO-GO Risk: Sequence decoding <10% CPU time
```

---

## Profiling Workflow

### Step 1: Baseline Measurement

```bash
# Run once to get baseline throughput
./target/release/bam-profiling ../test-data/test.bam
```

Note the throughput (Mrec/sec, Mbp/sec).

### Step 2: Time Profiler Analysis

```bash
# Run with Instruments
instruments -t "Time Profiler" ./target/release/bam-profiling ../test-data/test.bam

# Open trace
open instrumentscli*.trace
```

In Instruments:
1. Look at **Call Tree** view
2. Sort by **Self Time** (descending)
3. Find functions with highest CPU %
4. Drill down into noodles-bam crate
5. Locate sequence decoding functions

### Step 3: Flamegraph Visualization

```bash
# Generate flamegraph
cargo flamegraph --bin bam-profiling -- ../test-data/test.bam

# Open in browser
open flamegraph.svg
```

Look for wide bars (high CPU time) in:
- `decode_base`
- `read_sequence`
- BGZF-related functions

### Step 4: Document Findings

Update `PHASE_0_RESEARCH.md` with:
```markdown
## ARM Profiling Results (Days 5-6)

**Test Data**: HG00096.chrom22.bam (1M records, 100 Mbp)
**Platform**: M1/M2/M3 Mac (specify model)

### CPU Time Distribution

| Function | CPU % | Notes |
|----------|-------|-------|
| decode_base | X% | Sequence 4-bit → ASCII |
| BGZF decompress | Y% | Block decompression |
| read_record | Z% | Overall record parsing |
| Other | W% | I/O, overhead |

### Throughput

- Baseline: X Mrec/sec
- Throughput: Y Mbp/sec

### NEON Opportunity Validation

**Sequence decoding**: X% of CPU time

**Decision**: GO / NO-GO / BORDERLINE

**Rationale**: [Based on % CPU time in sequence decoding]
```

---

## Test Data Sizes

Recommended test files:

### Small (Quick Iteration)
- **10K records** (~1 MB compressed)
- Purpose: Fast iteration, development
- Generation time: <1 second

### Medium (Profiling)
- **100K records** (~10 MB compressed)
- Purpose: Profiling, validation
- Generation time: ~10 seconds

### Large (Scaling Validation)
- **1M records** (~100 MB compressed)
- Purpose: Real-world validation
- Source: Public data (1000 Genomes)

---

## Troubleshooting

### Error: "No such file or directory"

Ensure test data exists:
```bash
ls -lh ../test-data/
```

If empty, download or generate test data (see Step 2 above).

### Error: "Invalid BAM file"

Verify file is valid:
```bash
samtools view -H ../test-data/test.bam
```

### Instruments Not Found

macOS only. Install Xcode Command Line Tools:
```bash
xcode-select --install
```

### Flamegraph Not Found

Install cargo-flamegraph:
```bash
cargo install flamegraph
```

---

## Next Steps After Profiling

### If GO (Sequence Decoding ≥15%)

1. Document findings in `PHASE_0_RESEARCH.md`
2. Create `DECISION.md` with GO decision
3. Proceed to **Phase 1: Minimal BAM Reader**

### If NO-GO (Sequence Decoding <10%)

1. Document findings in `PHASE_0_RESEARCH.md`
2. Create `DECISION.md` with NO-GO decision
3. Pivot to noodles integration (document learnings)
4. Consider: Is BGZF parallelization alone worth it? (6.5× is substantial)

### If BORDERLINE (10-14%)

1. Calculate combined speedup: NEON sequence + parallel BGZF
2. If combined ≥2×: **GO** (but document risk)
3. If combined <2×: **NO-GO** (not worth maintenance burden)

---

**Status**: Profiling tool ready
**Next**: Download/generate test data, run profiling
**Decision**: Day 7 (after profiling validation)
