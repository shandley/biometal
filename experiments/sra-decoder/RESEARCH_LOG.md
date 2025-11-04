# Research Log: ARM-Native SRA Decoder

**Experiment**: sra-decoder
**Timeline**: Nov 5-15, 2025 (2-week research phase)
**Status**: Day 1 - Format Analysis

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
