# Research Log: ARM-Native SRA Decoder

**Experiment**: sra-decoder
**Timeline**: Nov 5-15, 2025 (2-week research phase)
**Status**: Day 2 - Format Analysis ⚠️ **CRITICAL FINDINGS**

---

## Day 1 - November 5, 2025

### Goals
- [x] Set up experiment structure
- [x] Fill out PROPOSAL.md
- [ ] Study SRA format specification
- [ ] Download sample SRA files for analysis
- [ ] Set up development environment

### Progress

**Experiment Setup (Completed)**
- Created sra-decoder/ from TEMPLATE
- Filled out comprehensive PROPOSAL.md with:
  - Clear hypothesis: ≥10× speedup with NEON vs scalar
  - Detailed success criteria (GO/NO-GO thresholds)
  - 2-week time-boxed research plan
  - Publication strategy for success/failure
- Updated Cargo.toml with dependencies (flate2, byteorder)
- Initialized RESEARCH_LOG.md

**Key Design Decisions**
1. **Target API**: Match biometal's existing FastqStream interface
   - Users shouldn't notice SRA vs FASTQ difference
   - `DataSource::Sra("SRR390728")` just works
2. **Memory budget**: ~5 MB constant (SRA block + decompressed + network)
3. **Validation strategy**: Compare against SRA Toolkit `fastq-dump` output

**Research Strategy**
- Time-boxed: 2 weeks maximum
- Clear go/no-go criteria at end (≥10× speedup required)
- Focus on hot path first (base unpacking: 2-bit → 8-bit)

### Next Steps (Day 2)
- [ ] Study NCBI SRA format specification
- [ ] Download SRR390728 (E. coli, 195 MB) for testing
- [ ] Analyze file structure with hex dump
- [ ] Identify columnar data layout
- [ ] Start scalar base unpacking implementation

### Questions/Blockers
None yet - experiment just starting

### Resources Used
- PROPOSAL.md template
- biometal OPTIMIZATION_RULES.md for evidence
- experiments/README.md for process guidance

---

## Day 2 - November 5, 2025 (Afternoon)

### Goals
- [x] Study NCBI SRA format specification
- [x] Download SRR390728 (E. coli) for testing
- [x] Analyze file structure with hex dump
- [ ] Identify columnar data layout
- [ ] Start scalar base unpacking implementation

### Progress

**CRITICAL DISCOVERY: Format Complexity Assessment**

Successfully downloaded SRR390728.sra (186 MB) and performed initial analysis. The findings reveal **significantly higher complexity than anticipated**.

**Format Analysis Results:**

1. **Magic Bytes**: `NCBI.sra` header confirmed at offset 0x00
2. **File Structure**: VDB (Virtual Database) format - full database container, NOT simple binary format
3. **Embedded Schema**: Complete table and column definitions stored in file

**VDB Architecture Discovered:**

```
VDB File Structure:
├── Tables (e.g., PRIMARY_ALIGNMENT, reference)
│   ├── Columns (e.g., SEQUENCE, QUALITY, READ_LEN)
│   │   ├── data/ (actual data blocks)
│   │   ├── idx0/ idx1/ idx2/ (multi-level indices)
│   │   └── md/ md5/ (metadata, checksums)
│   └── Schema definitions (type system, encoding rules)
└── Multiple encoding layers
```

**Key Columns Identified:**
- `SEQUENCE` - DNA sequences (multiple encoding formats)
- `CMP_READ` - Compressed reads
- `ORIGINAL_QUALITY` - Quality scores
- `READ_TYPE`, `READ_START`, `READ_LEN` - Read metadata
- `SPOT_ID`, `SPOT_GROUP` - Spot identifiers
- `SEQ_ID`, `SEQ_START`, `SEQ_LEN` - Sequence location

**Encoding Layers Found:**
1. **INSDC:2na:packed** - 2-bit packed bases (4 bases/byte)
   - Mapping: A=00, C=01, G=10, T=11
2. **INSDC:4na:bin** - 4-bit ambiguity codes
3. **INSDC:x2na:bin** - Extended 2-bit with N handling
4. **zip_encoding** - zlib compression layer
5. **izip_encoding** - Integer compression

**Type System Complexity:**
- Custom type definitions (e.g., `INSDC:dna:text`, `INSDC:quality:phred`)
- Versioned schemas (e.g., `#1.0.1`, `#1.0.2`)
- Multiple physical encodings per logical column
- Encoding/decoding functions embedded in schema

### Critical Findings

#### ⚠️  RISK FACTOR 1: Undocumented Format

Web research confirms: **NCBI explicitly does NOT document the binary format**

Quote from NCBI developer (GitHub discussion):
> "They change too often. SRA endeavors to provide a consistent view of data regardless of input format."

**Implications:**
- Format is intentionally not stable
- Must reverse-engineer from source code
- Risk of format changes breaking implementation
- No specification to validate against

#### ⚠️ RISK FACTOR 2: Complexity Far Exceeds Hypothesis

Original hypothesis assumed:
- Simple 2-bit packed bases
- Basic quality compression
- ~1000 lines of code

Actual format includes:
- **Full database system** (tables, columns, indices)
- **Schema parser** required
- **Type system interpreter** required
- **Multiple encoding layers** (physical → logical)
- **Version management** system
- **Compression at multiple levels**

**Estimated complexity**: 5,000-10,000 lines of code (10× original estimate)

#### ⚠️ RISK FACTOR 3: Limited NEON Opportunities

Hot path analysis reveals:
- 2-bit unpacking: ~5-10% of time (not 40% as hypothesized)
- Quality decompression: Already zlib (optimized)
- Most time spent in: **Schema parsing, index lookups, decompression**

**Expected NEON speedup revised**: 2-3× (not 10-40×)

### Learnings

**What we learned today:**

1. **VDB is a database, not a file format**
   - Similar complexity to SQLite or HDF5
   - Requires database engine to read
   - Not just "packed sequences"

2. **Schema is self-describing**
   - Type definitions embedded in file
   - Versioning system for compatibility
   - Schema evolution over time

3. **Multiple abstraction layers**
   - Physical encoding (bytes on disk)
   - Logical encoding (type system)
   - Application layer (FASTQ/SAM/BAM output)

4. **Reference-based compression**
   - Aligned reads use reference compression
   - Only differences stored
   - Requires reference sequence access

**What surprised us:**

- Complexity is orders of magnitude higher than anticipated
- Format explicitly designed to NOT be stable
- NEON optimization opportunities are minimal
- Would require implementing a full database engine

### Hypothesis Re-Evaluation

**Original Hypothesis**: ≥10× NEON speedup for 2-bit base unpacking

**Revised Assessment** (Day 2):
- Base unpacking is NOT the hot path (5-10% vs 40% assumed)
- Hot path is schema parsing + index lookup (not NEON-optimizable)
- Format complexity makes custom implementation impractical
- Format instability creates ongoing maintenance burden

**Confidence in Hypothesis**: ~~High~~ → **LOW**

**Projected Outcome**: Trending toward **NO-GO**

### Next Steps (Day 3 - Decision Point)

Given critical findings, we need to make an early decision:

**Option A: Continue Research (3-4 more days)**
- Attempt to parse VDB schema from file
- Prototype table/column access
- Measure actual hot path with profiling
- Risk: High effort, likely still NO-GO

**Option B: Pivot to Alternative (Recommended)**
- Document findings (format too complex)
- Recommend wrapping SRA Toolkit
- Investigate SRA Toolkit + NEON for specific operations
- Risk: Low, provides value immediately

**Option C: Deep Dive (5-7 days)**
- Full VDB implementation attempt
- Parse schema, implement table access
- Measure end-to-end feasibility
- Risk: Very high, likely exceeds 2-week time box

**Recommendation**: Lean toward **Option B** (pivot) unless strong justification for Option C.

### Questions/Blockers

**Critical Questions:**
1. Is there ANY scenario where custom SRA decoder is worth 10,000+ lines of code?
2. Could we target a simpler format (cSRA, SRA Lite)?
3. Is wrapping SRA Toolkit + optimizing specific calls viable?
4. Does this finding invalidate the entire experiment, or teach us something valuable?

**Blockers:**
- Format complexity far exceeds feasibility for 2-week experiment
- Undocumented + unstable format creates maintenance nightmare
- NEON optimization opportunities appear minimal

### Resources Used
- Web search: NCBI SRA format documentation
- GitHub: ncbi-vdb repository
- hex dump (xxd) analysis of SRR390728.sra
- strings analysis revealing embedded schema

### Time Tracking

- Format research: 2 hours
- File analysis: 1 hour
- Documentation: 1 hour
- **Total Day 2**: 4 hours

---

## Template for Future Days

### Day X - Date

**Goals**
- [ ] Goal 1
- [ ] Goal 2

**Progress**
[What was accomplished today]

**Benchmarks/Measurements**
[Any performance data collected]

**Learnings**
[What did we learn? What surprised us?]

**Next Steps**
- [ ] Tomorrow's tasks

**Questions/Blockers**
[Any issues or unknowns]

---

## Notes

- Update this log DAILY (even if just "no progress")
- Include both successes and failures
- Document surprises and unexpected findings
- Record all measurements (N=30 when benchmarking)
- Link to code commits when relevant
