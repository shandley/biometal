# BAM Profiling Infrastructure - Status

**Date**: November 8, 2025
**Status**: ✅ Ready to use
**Phase**: Phase 0 Day 5-6 (profiling preparation)

---

## What's Ready

### 1. Profiling Tool ✅

**Location**: `profiling/`

**Components**:
- `Cargo.toml` - Dependencies (noodles-bam 0.68)
- `src/main.rs` - Profiling binary (compiles successfully)
- `README.md` - Complete usage guide
- `get_test_data.sh` - Download public BAM files
- `generate_test_bam.sh` - Generate synthetic BAM files

**Build**:
```bash
cd profiling
cargo build --release
```

**Binary**: `./target/release/bam-profiling`

### 2. Documentation ✅

**README.md** contains:
- Quick start guide
- Profiling workflows (Instruments, flamegraph)
- GO/NO-GO criteria (≥15% CPU time in sequence decoding)
- Troubleshooting
- Next steps based on results

### 3. Helper Scripts ✅

**get_test_data.sh**:
- Downloads 1000 Genomes BAM file (~100-200 MB)
- Chromosome 22, Illumina low coverage
- ~500K-1M records

**generate_test_bam.sh**:
- Generates synthetic BAM files
- Configurable record count (10K, 100K, 1M)
- Requires: samtools (install via `brew install samtools`)

---

## Next Steps (To Complete Profiling)

### Step 1: Install samtools (if needed)

```bash
brew install samtools
```

**Note**: Only needed for synthetic BAM generation. Public data download doesn't require samtools.

### Step 2: Get Test Data

**Option A**: Download public BAM
```bash
cd profiling
./get_test_data.sh
```

**Option B**: Generate synthetic BAM
```bash
cd profiling
./generate_test_bam.sh 100000  # 100K records
```

### Step 3: Run Profiling Tool

```bash
cd profiling
./target/release/bam-profiling ../test-data/HG00096.chrom22.bam
```

Expected output:
```
Profiling BAM parsing with noodles
File: ../test-data/HG00096.chrom22.bam

Header parsing: 1.2ms

Parsed 1000000 records (87245.3 rec/s)...

Parsing complete!

Results:
  Records parsed: 1000000
  Total sequence length: 100000000 bp
  ...
  Throughput: 87304.2 records/sec
```

### Step 4: Profile with Instruments (macOS)

```bash
# Run with Time Profiler
instruments -t "Time Profiler" \
    ./target/release/bam-profiling \
    ../test-data/HG00096.chrom22.bam

# Open trace file
open instrumentscli*.trace
```

**What to look for**:
- CPU time % in `decode_base` function
- CPU time % in `read_sequence` function
- CPU time % in BGZF decompression

### Step 5: Analyze Results

**GO Criteria** (sequence decoding ≥15% CPU time):
- ✅ **STRONG GO**: ≥30% → 3-4× overall speedup likely
- ✅ **GO**: 15-29% → 2× overall speedup likely
- ⚠️ **BORDERLINE**: 10-14% → Re-evaluate with BGZF
- ❌ **NO-GO**: <10% → Unlikely to reach 2× threshold

### Step 6: Document Findings

Update `PHASE_0_RESEARCH.md` with profiling results:

```markdown
## ARM Profiling Results (Days 5-6)

**Test Data**: HG00096.chrom22.bam (1M records, 100 Mbp)
**Platform**: [Mac model]

### CPU Time Distribution

| Function | CPU % | Notes |
|----------|-------|-------|
| decode_base | X% | Sequence 4-bit → ASCII |
| BGZF decompress | Y% | Block decompression |
| read_record | Z% | Overall record parsing |
| Other | W% | I/O, overhead |

### Decision: GO / NO-GO / BORDERLINE

**Rationale**: [Based on evidence]
```

---

## Files Created

```
profiling/
├── Cargo.toml                  # Dependencies (noodles 0.68)
├── Cargo.lock                  # Lockfile
├── src/
│   └── main.rs                 # Profiling binary
├── target/
│   └── release/
│       └── bam-profiling       # Compiled binary ✅
├── README.md                   # Complete usage guide
├── get_test_data.sh            # Download script (executable)
├── generate_test_bam.sh        # Generation script (executable)
└── PROFILING_STATUS.md         # This file
```

---

## Current State

**Profiling tool**: ✅ Built successfully
**Test data**: ❌ Not yet downloaded (requires Step 2)
**Profiling**: ❌ Not yet run (requires Steps 3-5)
**Decision**: ⏳ Pending (Day 7 after profiling)

---

## Time Estimate

**Step 1** (Install samtools): 2 minutes
**Step 2** (Get test data): 5-10 minutes (download) or 1 minute (generate)
**Step 3** (Run profiling tool): 1-2 minutes
**Step 4** (Instruments profiling): 5-10 minutes
**Step 5** (Analyze results): 10-15 minutes
**Step 6** (Document findings): 10-15 minutes

**Total**: ~30-45 minutes to complete profiling validation

---

## Expected Outcome

### If GO (≥15% CPU time in sequence decoding)

1. Create `DECISION.md` with GO decision
2. Document findings in `PHASE_0_RESEARCH.md`
3. Proceed to **Phase 1: Minimal BAM Reader** (Week 2-3)

### If NO-GO (<10% CPU time)

1. Create `DECISION.md` with NO-GO decision
2. Document findings and learnings
3. Pivot to noodles integration (recommend direct usage)
4. Consider: Is BGZF parallelization alone worth it? (6.5× is substantial)

### If BORDERLINE (10-14%)

1. Calculate combined speedup estimate: NEON + parallel BGZF
2. If combined ≥2×: **GO** (document risk)
3. If combined <2×: **NO-GO** (not worth maintenance burden)

---

## Critical Question

**Will sequence decoding be a bottleneck?**

This profiling will definitively answer this question. The NEON optimization hypothesis depends on sequence decoding consuming significant CPU time (≥15%).

If validated, we proceed with confidence. If not, we've saved weeks of implementation effort by validating early.

---

**Status**: Infrastructure ready, awaiting test data and profiling run
**Next**: Download/generate test data, run profiling (Steps 2-5)
**Decision Point**: Day 7 (November 15, 2025)
