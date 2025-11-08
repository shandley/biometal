# biometal Alignment Formats: Strategic Roadmap

**Date**: November 8, 2025
**Status**: Dual-format strategy approved
**Timeline**: 15-18 weeks total

---

## Vision

Build **two complementary alignment formats** optimized for modern ARM hardware:

1. **BAM** (compatibility first): ARM-optimized, standards-compliant
2. **CAF** (performance first): Columnar, NEON-optimized, 5-10× faster

**Philosophy**: Users choose trade-offs (compatibility vs performance)

**Publication Goal**: "Dual-Format Strategy for ARM-Native Bioinformatics: BAM and CAF"
- Target: Bioinformatics (Oxford) or BMC Bioinformatics
- Timeline: Submit Q2 2026 (~6 months from now)

---

## Strategic Rationale

### Why Both Formats?

**BAM (Native Implementation)**:
- ✅ Compatibility with ecosystem (samtools, IGV, etc.)
- ✅ Standards-compliant (SAM/BAM spec v1.6)
- ✅ Production-ready for clinical/sharing use
- ✅ Foundation for CAF (format understanding)
- 🎯 Target: **2-3× ARM speedup** via NEON + parallel BGZF

**CAF (Innovation)**:
- ✅ Optimized for modern hardware (NEON, GPU, Neural Engine)
- ✅ Analytical performance (5-10× faster for operations)
- ✅ Research tool (experimental, okay if low adoption)
- ✅ Pioneer format design (inspire ecosystem)
- 🎯 Target: **5-10× speedup** over BAM

### The 2009 vs 2025 Problem

**BAM was designed for 2009**:
- Disk expensive ($1/GB) → 4-bit compression
- CPU cheap → Burn cycles unpacking 4-bit
- gzip standard → No modern codecs

**2025 reality**:
- Disk cheap ($0.01/GB) → Storage not bottleneck
- ARM SIMD powerful → CPU cycles valuable
- zstd/lz4 available → Faster compression

**CAF addresses**: Modern hardware, modern trade-offs

---

## Timeline Overview

```
┌──────────────────────────────────────────────────────────────┐
│ Week 1-12: Native BAM Implementation (9-12 weeks)            │
├──────────────────────────────────────────────────────────────┤
│   Week 1:     Phase 0 - Research & Architecture              │
│   Week 2-3:   Phase 1 - Minimal BAM Reader                   │
│   Week 4:     Phase 2 - Streaming Architecture               │
│   Week 5-6:   Phase 3 - ARM NEON Optimization ◄─┐           │
│   Week 7:     Phase 4 - SAM Support              │           │
│   Week 8-10:  Phase 5 - CRAM (optional)          │           │
│   Week 11-12: Phase 6 - Production Polish        │           │
│                                                   │           │
├───────────────────────────────────────────────────┼───────────┤
│                                                   │           │
│ Week 13-18: Columnar Alignment Format (6 weeks)  │           │
├───────────────────────────────────────────────────┘  Builds   │
│   Week 13-14: Phase 1 - Core Format                   on      │
│   Week 15-16: Phase 2 - NEON Optimization            BAM      │
│   Week 17:    Phase 3 - Indexing & Queries          foundation│
│   Week 18:    Phase 4 - Production Polish                     │
│                                                                │
├────────────────────────────────────────────────────────────────┤
│ Week 19-20: Comprehensive Benchmarking (PAPER)                │
├────────────────────────────────────────────────────────────────┤
│   Week 19:    Data collection (N=30, diverse datasets)        │
│   Week 20:    Analysis, figures, tables generation            │
│                                                                │
├────────────────────────────────────────────────────────────────┤
│ Week 21-24: Paper Writing & Submission                        │
├────────────────────────────────────────────────────────────────┤
│   Week 21-22: Draft paper, create figures                     │
│   Week 23:    Internal review, revisions                      │
│   Week 24:    Final polish, submit to journal                 │
└────────────────────────────────────────────────────────────────┘

Total: 24 weeks (~6 months) to submission
Implementation: 18 weeks
Paper: 6 weeks (benchmarking + writing)
```

---

## Experiment 1: Native BAM Implementation

### Status
**Phase**: Phase 0 - Research & Architecture (Week 1, Day 1)
**Timeline**: 9-12 weeks
**Started**: November 8, 2025

### Key Deliverables (So Far)

✅ **Phase 0 Documents**:
- `PROPOSAL.md` (800+ lines) - Complete implementation plan
- `PHASE_0_RESEARCH.md` (600+ lines) - Format & noodles analysis
- `NEON_OPPORTUNITIES.md` (500+ lines) - Optimization analysis
- noodles source cloned and analyzed

### Primary NEON Opportunity

**4-bit Sequence Decoding**:
```rust
// Current noodles: Sequential match-based lookup
fn decode_base(n: u8) -> u8 {
    match n & 0x0f {
        0 => b'=', 1 => b'A', 2 => b'C', ...
    }
}

// NEON optimization: Parallel table lookup
unsafe fn decode_sequence_neon(packed: &[u8]) {
    let lookup = vld1q_u8(LOOKUP_TABLE.as_ptr());
    // Process 16 packed bytes (32 bases) at once
    let bases = vqtbl1q_u8(lookup, packed);  // 16 lookups in parallel
}
```

**Estimated speedup**:
- Sequence decoding: **3-6× operation speedup**
- Overall BAM parsing: **2-3× total** (validated by profiling this week)

### Success Criteria

| Phase | Metric | Target |
|-------|--------|--------|
| Phase 0 | NEON opportunities | ≥2 identified |
| Phase 1 | Correctness | 100% vs noodles |
| Phase 3 | ARM speedup | **≥2×** |
| Phase 6 | Tests | 200+ |

### Current Status

**Day 1 Complete**:
- ✅ Specification analysis (SAM/BAM v1.6)
- ✅ noodles source analysis
- ✅ NEON opportunity identified (high confidence)
- ⏳ Profiling pending (Day 5-6)
- ⏳ GO/NO-GO decision (Day 7)

### Next Steps (Week 1)

**Days 2-4**: Complete specification deep dive
**Days 5-6**: ARM profiling validation
- Download diverse BAM files
- Profile noodles with Instruments
- Validate sequence decoding is bottleneck (need ≥15% CPU time)

**Day 7**: GO/NO-GO Decision
- If sequence decoding ≥15% CPU time → **GO to Phase 1**
- If <10% → Pivot to noodles integration

---

## Experiment 2: Columnar Alignment Format (CAF)

### Status
**Phase**: Design Phase (Pre-implementation)
**Timeline**: 6 weeks (after BAM Phase 3-4)
**Expected Start**: ~December 20, 2025

### Key Deliverables (So Far)

✅ **Design Documents**:
- `PROPOSAL.md` (1,000+ lines) - 6-week implementation plan
- `CAF_SPECIFICATION.md` (1,500+ lines) - Format design v1.0
- `README.md` - Overview and strategic positioning

### Key Design Features

**1. Columnar Block Layout**:
```
Block (10,000 records):
  positions:    [i32; 10000]      (zstd compressed)
  mapq:         [u8; 10000]       (raw or RLE)
  sequences:    [u8; total_len]   (ASCII, lz4)  ← Pre-decoded!
  qualities:    [u8; total_len]   (raw)
  ...
```

**2. No 4-bit Encoding** (Key Difference):
- BAM: 4-bit → unpack every time (CPU overhead)
- CAF: ASCII → ready to use (zero overhead)
- Trade-off: 2× storage for sequence data

**3. Modern Compression**:
- zstd level 3 (positions, metadata): 2-3× faster than gzip
- lz4 (sequences): Extremely fast (>GB/s decompression)
- Raw (qualities): High entropy, doesn't compress

### Expected Performance

| Operation | BAM | CAF | Speedup |
|-----------|-----|-----|---------|
| Parse 100K | 2.56s | 0.25s | **10×** |
| Quality filter | 2.00s | 0.08s | **25×** |
| Base counting | 1.50s | 0.06s | **25×** |
| **Overall** | **1.0×** | **5-10×** | Target |

**Storage**: 1.5-2× larger than BAM (acceptable trade-off)

### Success Criteria

| Phase | Metric | Target |
|-------|--------|--------|
| Phase 1 | Round-trip | 100% lossless |
| Phase 2 | Overall speedup | **≥5×** |
| Phase 4 | Tests | 200+ |

### Implementation Phases

**Phase 1** (2 weeks): Core format, BAM ↔ CAF conversion
**Phase 2** (2 weeks): NEON optimization, benchmarks
**Phase 3** (1 week): Indexing, region queries
**Phase 4** (1 week): Production polish, docs

### Use Cases

**✅ CAF is great for**:
- Batch processing pipelines
- ML training data prep
- Quality control
- Research workflows

**❌ CAF is NOT for**:
- Genome browsers (use BAM)
- Sharing/archival (use BAM)
- Clinical pipelines (use BAM)

---

## Strategic Positioning

### The Dual-Format Philosophy

**biometal offers BOTH**:

1. **BAM** - When you need:
   - Compatibility (samtools, IGV, etc.)
   - Standards compliance
   - Sharing with others
   - Clinical pipelines

2. **CAF** - When you want:
   - Maximum performance (5-10× faster)
   - Analytical queries
   - Research workflows
   - Modern hardware optimization

**Users choose**: Trade compatibility for speed, or vice versa

### Seamless Conversion

```bash
# Workflow 1: Compatibility
biometal process input.bam --output output.bam

# Workflow 2: Performance
biometal convert input.bam temp.caf      # Convert once
biometal process temp.caf --fast         # 5-10× faster
biometal convert temp.caf final.bam      # Back to BAM for sharing

# Workflow 3: Hybrid
biometal process input.bam --internal-caf  # Use CAF internally, BAM I/O
```

### Documentation Strategy

**Clear trade-offs**:
- "CAF is 5-10× faster but 1.5-2× larger and incompatible (initially)"
- "Use BAM for sharing, CAF for performance"
- "Lossless conversion both ways"

**Adoption**: Organic
- Document thoroughly
- Benchmark transparently
- Let community decide
- "Okay if nobody uses CAF" - experimental innovation

---

## Resource Requirements

### Time Investment

**Total**: 15-18 weeks (3.5-4.5 months)
- BAM: 9-12 weeks
- CAF: 6 weeks

**Part-time (50%)**: 30-36 weeks (~7-9 months)

### Infrastructure

**Have**:
- ✅ Mac ARM (M1/M2/M3/M4)
- ✅ Rust toolchain
- ✅ NEON expertise (proven in biometal)
- ✅ Streaming architecture (FASTQ/FASTA)

**Need**:
- [ ] Test BAM files (diverse sources)
- [ ] Benchmarking harness (have criterion)
- [ ] zstd/lz4 Rust crates (available)

---

## Risk Assessment

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| BAM NEON speedup <2× | Low | Medium | BGZF gives 1.6× floor already |
| CAF speedup <5× | Low | Medium | NEON ops proven 16-25× |
| Community confusion | Medium | Low | Clear documentation |
| Maintenance burden | Low | Low | Comprehensive tests |

### Strategic Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| CAF not adopted | High | Low | Experimental tool, okay |
| Format fragmentation | Low | Low | Lossless conversion |

---

## Success Definition

### Minimum Viable Success

**BAM** (Phase 6):
- ✅ 100% spec-compliant
- ✅ ≥2× ARM speedup
- ✅ Production-ready (200+ tests)
- ✅ Compatible with ecosystem

**CAF** (Phase 4):
- ✅ Lossless BAM ↔ CAF
- ✅ ≥5× analytical speedup
- ✅ Documented trade-offs
- ✅ CLI tools

**Both**:
- ✅ Users understand when to use each
- ✅ Seamless conversion workflow
- ✅ Published benchmarks

### Stretch Goals

- ⭐ BAM: ≥3× speedup (exceptional)
- ⭐ CAF: ≥10× speedup (exceptional)
- ⭐ CAF: GPU integration (Metal)
- ⭐ Community adoption (others implement CAF)
- ⭐ Published paper (modern format design)

---

## Current State (November 8, 2025)

### Project Structure

```
biometal/experiments/
├── native-bam-implementation/        ← Week 1 (Phase 0)
│   ├── PROPOSAL.md                   ✅ Complete
│   ├── PHASE_0_RESEARCH.md           ✅ Complete
│   ├── NEON_OPPORTUNITIES.md         ✅ Complete
│   ├── noodles/                      ✅ Cloned
│   └── [Phase 1-6 work...]           ⏳ Pending
│
├── columnar-alignment-format/        ← Design phase
│   ├── PROPOSAL.md                   ✅ Complete
│   ├── CAF_SPECIFICATION.md          ✅ Complete (v1.0 draft)
│   ├── README.md                     ✅ Complete
│   └── [Phase 1-4 work...]           ⏳ After BAM Phase 3-4
│
├── .experiments.toml                 ✅ Updated (both experiments)
└── ALIGNMENT_FORMATS_ROADMAP.md      ✅ This file
```

### Experiment Registry

**Updated**: `.experiments.toml` tracks both:
- `native-bam-implementation`: Phase 0 (research)
- `columnar-alignment-format`: Proposed (design)

---

## Next Steps

### Immediate (This Week)

**BAM Phase 0** (Days 2-7):
- [ ] Complete specification analysis
- [ ] Profile noodles on ARM (Instruments)
- [ ] Validate NEON opportunities
- [ ] GO/NO-GO decision (Day 7)

### Short-term (Weeks 2-12)

**BAM Phases 1-6**:
- [ ] Implement minimal BAM reader
- [ ] Streaming architecture
- [ ] NEON optimization
- [ ] SAM support
- [ ] Production polish

### Medium-term (Weeks 13-18)

**CAF Phases 1-4**:
- [ ] Core format + conversion
- [ ] NEON optimization
- [ ] Indexing
- [ ] Production polish

---

## Commitment

**We commit to**:
- Evidence-based design (profiling, benchmarks)
- Comprehensive documentation
- Honest assessment (publish negative results)
- Having fun innovating! :)

**Success means**:
- Two excellent alignment format options
- Users choose based on needs
- ARM-native bioinformatics platform
- Pioneer modern format design

---

**Status**: BAM Phase 0 (Day 1), CAF Design Phase
**Next Decision**: BAM GO/NO-GO (November 15, 2025)
**Expected Completion**: March-April 2026

Let's build the future of alignment formats! 🚀
