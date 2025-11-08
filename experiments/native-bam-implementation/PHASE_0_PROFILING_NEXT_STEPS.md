# Phase 0 Profiling - Next Steps

**Date**: November 8, 2025
**Status**: Profiling infrastructure ready, test data needed
**Phase**: Day 5-6 (Profiling validation)

---

## Current Status

### ✅ Completed

1. **Profiling tool built** (`profiling/target/release/bam-profiling`)
   - Compiles successfully with noodles 0.68
   - Measures throughput (records/sec, Mbp/sec)
   - Accesses sequences and quality scores to ensure decoding

2. **Helper scripts created**:
   - `get_test_data.sh` - Download 1000 Genomes BAM
   - `generate_test_bam.sh` - Generate synthetic BAM (requires samtools)

3. **Documentation complete**:
   - `README.md` - Usage guide
   - `PROFILING_STATUS.md` - Status and workflow

### ⏳ Pending

**Test data acquisition** - Need valid BAM file for profiling

---

## Options to Get Test Data

### Option 1: Install samtools and Generate Synthetic BAM

**Best for**: Quick testing, controlled data size

```bash
# Install samtools
brew install samtools

# Generate test BAM (100K records)
cd profiling
./generate_test_bam.sh 100000

# Run profiling
./target/release/bam-profiling ../test-data/synthetic_100000.bam
```

**Pros**: Fast, controlled, no large download
**Cons**: Requires samtools installation (~50 MB)

### Option 2: Download Public BAM

**Best for**: Real-world data, no dependencies

```bash
# Download small BAM from GIAB (Genome in a Bottle)
cd test-data

# HG002 (GIAB reference) - chr22 only (~200-500 MB)
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_BH8JCAADXX/Project_RM8391_L1/Sample_1192A1/alignment/1192A1.bam

# OR: 1000 Genomes Project (various sizes available)
# See: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/
```

**Pros**: Real-world data, representative workload
**Cons**: Large download (200 MB+), requires network

### Option 3: Use Existing Biomet al Test Data (if available)

```bash
# Check if biometal has any alignment files
find ../../.. -name "*.bam" -o -name "*.sam" 2>/dev/null

# If SAM files exist, convert to BAM with samtools
samtools view -b input.sam > output.bam
```

---

## Profiling Workflow (Once Data is Ready)

### Step 1: Run Baseline Profiling

```bash
cd profiling
./target/release/bam-profiling ../test-data/test.bam
```

**Expected output**:
```
Profiling BAM parsing with noodles
File: ../test-data/test.bam

Header parsing: 1.2ms

Parsed 100000 records (87245.3 rec/s)...

Results:
  Records parsed: 100000
  Total sequence length: 10000000 bp
  Elapsed time: 1.145 s

Throughput:
  87304.2 records/sec
  0.09 Mrec/sec
  8.7 Mbp/sec
```

### Step 2: Profile with macOS Instruments

```bash
# Run with Time Profiler
instruments -t "Time Profiler" \
    ./target/release/bam-profiling \
    ../test-data/test.bam

# Open trace file
open instrumentscli*.trace
```

**In Instruments**:
1. Switch to **Call Tree** view
2. Sort by **Self Time** (descending)
3. Find `decode_base` function
4. Find `read_sequence` function
5. Document CPU time %

### Step 3: Alternative - cargo-flamegraph

```bash
# Install flamegraph
cargo install flamegraph

# Generate flamegraph
cargo flamegraph --bin bam-profiling -- ../test-data/test.bam

# Open in browser
open flamegraph.svg
```

Look for wide bars in:
- `decode_base` (4-bit → ASCII conversion)
- `read_sequence` (sequence decoding)
- BGZF decompression

---

## Decision Criteria

### CPU Time in Sequence Decoding

| % CPU Time | Decision | Expected Speedup | Confidence |
|------------|----------|------------------|------------|
| ≥30% | **STRONG GO** | 3-4× overall | High |
| 15-29% | **GO** | 2× overall | Medium |
| 10-14% | **BORDERLINE** | 1.5-2× | Low |
| <10% | **NO-GO** | <1.5× | Risk |

### Combined with BGZF Parallel (6.5×)

Even if sequence decoding is modest (10-15%), combined with parallel BGZF decompression:
- Sequence NEON: 1.5-2× on decoding
- BGZF parallel: 6.5× on decompression
- **Overall**: Likely ≥2× total speedup

---

## Documentation Template

After profiling, update `PHASE_0_RESEARCH.md`:

```markdown
## ARM Profiling Results (Days 5-6)

**Date**: [Date]
**Test Data**: [File name] ([size], [num records])
**Platform**: [Mac model, M1/M2/M3/M4]

### CPU Time Distribution

| Function | CPU % | Notes |
|----------|-------|-------|
| decode_base | X% | 4-bit → ASCII conversion (NEON target) |
| read_sequence | Y% | Sequence decoding overhead |
| BGZF decompress | Z% | Block decompression (parallel target) |
| read_record | W% | Overall parsing overhead |
| Other | V% | I/O, allocation, etc. |

### Baseline Throughput

- **Throughput**: X Mrec/sec
- **Bandwidth**: Y Mbp/sec
- **Total time**: Z seconds

### NEON Opportunity Validation

**Sequence decoding**: X% of CPU time

**Decision**: [GO / NO-GO / BORDERLINE]

**Rationale**:
[Based on CPU time %, expected speedup calculation, risk assessment]

If GO:
- Estimated sequence NEON speedup: A×
- Estimated parallel BGZF speedup: 6.5× (proven)
- Combined estimate: B× overall

If NO-GO:
- Sequence decoding <10% → Not bottleneck
- Alternative: Use noodles directly, integrate parallel BGZF only
```

---

## Estimated Time

**Option 1** (samtools + synthetic):
- Install samtools: 5 min
- Generate BAM: 2 min
- Run profiling: 5-10 min
- Analyze: 10-15 min
- **Total**: ~25-35 min

**Option 2** (public download):
- Download BAM: 10-20 min (network dependent)
- Run profiling: 10-15 min
- Analyze: 10-15 min
- **Total**: ~30-50 min

---

## Recommendation

**For immediate testing**: Option 1 (samtools + synthetic)
- Fastest path to validation
- Controlled data size
- Representative of real workload

**For publication**: Option 2 (public data)
- Real-world validation
- Reproducible (public dataset)
- Stronger evidence for paper

**Best approach**: Do both!
1. Quick validation with synthetic (Option 1)
2. Validation with public data (Option 2)
3. Document both in PHASE_0_RESEARCH.md

---

## Current Blocker

**Test data acquisition** - Choose option above and execute

Once test data is available, profiling can be completed in ~30 minutes.

---

**Status**: Infrastructure ready ✅, awaiting test data ⏳
**Next**: Install samtools OR download public BAM
**Decision**: Day 7 (after profiling results available)
