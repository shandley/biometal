# Native BAM/SAM/CRAM Implementation: Proposal

**Date**: November 7, 2025
**Status**: Phase 0 - Research & Architecture
**Type**: Strategic Implementation Experiment
**Expected Duration**: 9-12 weeks

---

## Vision

Build ARM-first, biometal-native alignment format support that:
- Owns the full stack (no external format dependencies)
- Optimizes for ARM architecture (Apple Silicon, Graviton)
- Pioneers streaming-first patterns
- Provides compelling ergonomics

**Extended Vision**: BAM implementation serves as foundation for **Columnar Alignment Format (CAF)**, a modern columnar format optimized for ARM NEON (5-10× faster than BAM for analytical operations). See `experiments/columnar-alignment-format/` for CAF design.

---

## Hypothesis

A native BAM/SAM implementation optimized for ARM can provide:
1. **Performance**: ≥2× speedup on ARM vs architecture-agnostic libraries
2. **Control**: Full-stack ownership enables innovation
3. **Integration**: Seamless biometal API (streaming + NEON operations)
4. **Learning**: Deep format understanding enables future innovations

---

## Strategic Rationale

### Why Native Implementation?

**1. ARM-First Philosophy**
- biometal targets ARM (Mac, Graviton) as first-class platform
- Existing libraries are architecture-agnostic (not ARM-optimized)
- Opportunity to pioneer ARM-native bioinformatics

**2. Full-Stack Control**
- Own entire pipeline: FASTQ → operations → BAM
- Freedom to innovate (novel compression, streaming patterns)
- No dependency on external maintenance

**3. Learning & Platform Value**
- Deep understanding of alignment formats
- Foundation for future formats (VCF, CRAM, etc.)
- Platform for bioinformatics research

**4. Strategic Positioning**
- "Complete ARM-native bioinformatics toolkit"
- Not just operations, but comprehensive format support
- Differentiation in Rust bio ecosystem

### Why Now?

**Favorable Conditions**:
- ✅ No rush to publish (time to do it right)
- ✅ noodles source available (reference implementation)
- ✅ Proven ARM NEON capability (16-25× on operations)
- ✅ Parallel BGZF already implemented (6.5× speedup)
- ✅ Streaming architecture established (Rule 5: constant memory)

**Risk Mitigation**:
- Phased approach with decision points (can pivot/pause)
- noodles integration remains fallback option
- Each phase delivers value independently

---

## Scope

### In Scope

**Phase 1-2: Core BAM Support** (Must Have)
- ✅ BAM reading (binary format)
- ✅ Record parsing (all mandatory fields)
- ✅ CIGAR operations
- ✅ Optional tags
- ✅ Streaming API

**Phase 3: ARM Optimization** (Must Have)
- ✅ NEON-optimized parsing
- ✅ Parallel record processing
- ✅ Integrated operations (quality, GC, etc.)

**Phase 4: SAM Support** (Should Have)
- ✅ Text SAM reading/writing
- ✅ SAM ↔ BAM conversion

**Phase 5: CRAM Support** (Nice to Have)
- ⚠️ Reference-based compression (complex)
- ⚠️ Decision point: Defer if ROI unclear

**Phase 6: Production Polish** (Must Have)
- ✅ Comprehensive testing
- ✅ Documentation
- ✅ Examples

### Out of Scope (Initially)

- ❌ BAM indexing (BAI/CSI) - defer to Phase 7 or later
- ❌ Random access queries - defer
- ❌ BAM writing (focus on reading first)
- ❌ VCF/BCF/GFF formats - separate experiments

---

## Success Criteria

### Quantitative (Phase-by-Phase)

| Phase | Metric | Target | Threshold |
|-------|--------|--------|-----------|
| **Phase 0** | NEON opportunities | ≥3 identified | GO if ≥2 |
| **Phase 1** | Correctness | 100% vs noodles | Must match |
| **Phase 2** | API quality | User feedback | Cleaner than noodles |
| **Phase 3** | ARM speedup | ≥2× overall | GO if ≥1.5× |
| **Phase 3** | NEON ops | ≥5× per operation | Individual ops |
| **Phase 4** | SAM parity | Feature complete | vs noodles |
| **Phase 6** | Test coverage | 200+ tests | ≥95% coverage |

### Qualitative

**Must Achieve**:
- ✅ Correct: 100% BAM/SAM spec compliance
- ✅ Safe: No unsafe operations without validation
- ✅ Fast: Measurably faster on ARM
- ✅ Ergonomic: Better API than noodles integration
- ✅ Documented: Complete guides + examples

**Success Even If**:
- ⚠️ Speedup is "only" 1.5-2× (strategic value justifies)
- ⚠️ CRAM deferred (BAM/SAM sufficient for v1)
- ⚠️ Some advanced features deferred (indexing, etc.)

---

## Timeline & Phases

### Overview (9-12 weeks)

```
Week 1:     Phase 0 - Research & Architecture
Weeks 2-3:  Phase 1 - Minimal BAM Reader
Week 4:     Phase 2 - Streaming Architecture
Weeks 5-6:  Phase 3 - ARM NEON Optimization
Week 7:     Phase 4 - SAM Support
Weeks 8-10: Phase 5 - CRAM (optional, decision point)
Weeks 11-12: Phase 6 - Production Polish
```

### Decision Points

| Week | Decision | GO Criteria | NO-GO Action |
|------|----------|-------------|--------------|
| 1 | Continue to Phase 1? | ≥2 NEON opportunities | Pivot to noodles integration |
| 3 | Continue to Phase 2? | Clean code, no blockers | Reassess scope |
| 4 | Continue to Phase 3? | API compelling | Consider hybrid approach |
| 6 | Continue to Phase 4? | ≥1.5× ARM speedup | Ship BAM-only, defer SAM |
| 7 | Attempt CRAM? | User need + capacity | Defer CRAM |
| 12 | Ship it? | All must-haves met | Extended polish |

---

## Phase 0 Objectives (Week 1)

### Primary Goals

**1. Format Deep Dive**
- Study SAM/BAM/CRAM specifications thoroughly
- Map all edge cases and complexity
- Understand design decisions in spec

**2. noodles Analysis**
- Study noodles-bam source architecture
- Extract edge case handling
- Identify optimization opportunities
- Document lessons learned

**3. ARM Profiling**
- Profile noodles on Mac ARM with real BAM files
- Identify actual hotspots (not theoretical)
- Measure ARM vs x86 performance gap
- Find NEON opportunities

**4. Architecture Design**
- Design ARM-first architecture
- Plan NEON optimization strategy
- Define streaming API patterns
- Integration with biometal operations

### Deliverables

- [ ] `PHASE_0_RESEARCH.md` - Format analysis + noodles study
- [ ] `PROFILING_RESULTS.md` - ARM profiling data
- [ ] `ARCHITECTURE.md` - biometal-native design
- [ ] `NEON_OPPORTUNITIES.md` - Optimization analysis
- [ ] `DECISION.md` - GO/NO-GO with evidence

### Success Criteria

**GO to Phase 1 if**:
- ✅ Identified ≥2 NEON opportunities with ≥3× potential each
- ✅ Profiling shows clear ARM optimization path
- ✅ Architecture design is clean and feasible
- ✅ No major technical blockers discovered

**NO-GO (pivot to noodles) if**:
- ❌ No clear NEON opportunities found
- ❌ noodles already ARM-optimized (unlikely but check)
- ❌ Complexity far exceeds estimates
- ❌ Fundamental architectural issues

---

## Resource Requirements

### Time

- **Full-time equivalent**: 9-12 weeks
- **Part-time (50%)**: 18-24 weeks
- **Decision points**: Every 1-2 weeks (can pause/pivot)

### Infrastructure

**Required**:
- ✅ Mac ARM (M1/M2/M3/M4) - you have
- ✅ Rust toolchain - you have
- ✅ noodles source (reference) - available

**Needed**:
- [ ] Test BAM files (diverse sources)
- [ ] Profiling tools (cargo-flamegraph, Instruments)
- [ ] CI for x86 validation (optional but recommended)

### Knowledge

**Have**:
- ✅ Rust expertise
- ✅ ARM NEON experience (from operations)
- ✅ Parallel programming (BGZF, rayon)
- ✅ Streaming architecture (FASTQ/FASTA)

**Acquiring**:
- [ ] BAM/SAM format expertise (Phase 0)
- [ ] BGZF internals (partially have)
- [ ] Alignment format edge cases

---

## Risk Assessment

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| NEON gains < 2× | Medium | Medium | Strategic value still justifies |
| Edge cases complex | Medium | Medium | Use noodles as reference |
| Timeline overruns | Medium | Low | Decision points allow pause |
| CRAM too complex | High | Low | Make optional (Phase 5) |
| Maintenance burden | Low | Medium | Comprehensive tests, docs |

### Strategic Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| noodles improves | Low | Low | We're pioneering ARM-first, different value prop |
| Format spec changes | Low | Medium | Monitor samtools/htslib development |
| User adoption low | Medium | Low | Even internal use has value |

---

## Comparison to Alternatives

### vs Noodles Integration (Previous Plan)

| Aspect | Native Implementation | Noodles Integration |
|--------|----------------------|---------------------|
| **Development Time** | 9-12 weeks | 0 weeks (done) |
| **Performance** | ≥2× on ARM (goal) | 1.0× (baseline) |
| **Control** | Full stack ownership | Dependency on noodles |
| **API** | biometal-native | Adapter required |
| **Learning** | Deep format knowledge | Surface-level |
| **Innovation** | Platform for new ideas | Limited by noodles API |
| **Risk** | Medium (can pivot) | Low (proven) |
| **Strategic Value** | High (differentiation) | Medium (integration) |

**Decision**: Native implementation chosen for strategic value + learning + ARM optimization potential

### vs rust-htslib

| Aspect | biometal-native | rust-htslib |
|--------|----------------|-------------|
| **Dependencies** | Pure Rust | C (HTSlib) |
| **ARM Optimization** | First-class | Via HTSlib (unknown) |
| **Philosophy Fit** | Perfect | Conflicts (C deps) |
| **Build Complexity** | Cargo only | C compiler required |

**Decision**: rust-htslib rejected due to C dependency (established in format-integration experiment)

---

## Success Definition

### Minimum Viable Success (MVP)

At end of Phase 6, we have:
- ✅ BAM reading that matches noodles correctness
- ✅ ≥1.5× speedup on ARM
- ✅ Clean streaming API
- ✅ Integration with biometal operations
- ✅ Comprehensive tests
- ✅ Documentation

**Even if**:
- SAM support is basic
- CRAM is deferred
- Some features deferred (indexing, writing)

### Stretch Goals

- ⭐ ≥3× ARM speedup (exceptional)
- ⭐ Novel streaming patterns adopted by others
- ⭐ CRAM support
- ⭐ Published benchmarks influencing ecosystem

### Acceptable Outcomes

**Scenario A**: Achieve 1.5-2× speedup
- Strategic value justifies (ownership, learning, platform)
- Ship as production feature

**Scenario B**: Achieve <1.5× speedup but excellent API
- Consider hybrid (noodles backend, biometal API)
- Learning value still high

**Scenario C**: Discover fundamental blocker (unlikely)
- Pivot to noodles integration
- Document learnings for community

---

## Learnings from format-integration Experiment

### What We Validated

✅ **noodles is excellent** - high bar to beat
✅ **Iterator pattern has overhead** - design for zero-copy
✅ **Benchmark early** - time-boxed experiments work
✅ **Evidence-based thresholds** - but apply contextually

### What Changed Our Mind

**New context**:
- No rush to publish (time to do it right)
- Strategic value beyond just performance
- Learning + innovation platform
- noodles source available (reference)

**Revised threshold**:
- Originally: ≥5× speedup required for format reimplementation
- Now: ≥2× speedup + strategic value justifies native implementation
- Rationale: Strategic goals (ARM-first, ownership, learning) add value beyond pure performance

---

## Experiment Methodology

### Research Approach (Phase 0)

1. **Specification Study** (2 days)
   - Read SAM/BAM spec v1.6 completely
   - Document complexity hotspots
   - Map edge cases

2. **noodles Source Analysis** (2 days)
   - Study architecture decisions
   - Extract edge case handling
   - Document optimization opportunities

3. **ARM Profiling** (2 days)
   - Profile noodles with real BAM files
   - Identify bottlenecks
   - Measure ARM vs x86 gap

4. **Architecture Design** (1 day)
   - Synthesize learnings
   - Design biometal approach
   - Plan NEON strategy

### Validation Approach (Phase 1+)

**Differential Testing**:
```rust
// Compare every output to noodles
#[test]
fn test_matches_noodles() {
    let biometal_output = biometal::bam::read("test.bam")?;
    let noodles_output = noodles::bam::read("test.bam")?;
    assert_eq!(biometal_output, noodles_output);
}
```

**Property Testing**:
```rust
// Round-trip property
proptest! {
    #[test]
    fn roundtrip_parse_serialize(record in arbitrary_record()) {
        let serialized = serialize(&record)?;
        let parsed = parse(&serialized)?;
        assert_eq!(record, parsed);
    }
}
```

---

## Documentation Plan

### Phase 0 Outputs

```
experiments/native-bam-implementation/
├── PROPOSAL.md                 # This file
├── PHASE_0_RESEARCH.md        # Format + noodles analysis
├── PROFILING_RESULTS.md       # ARM profiling data
├── ARCHITECTURE.md            # Design document
├── NEON_OPPORTUNITIES.md      # Optimization analysis
└── DECISION.md                # GO/NO-GO with evidence
```

### Future Phases

```
experiments/native-bam-implementation/
├── PHASE_1_IMPLEMENTATION.md  # Core BAM reader
├── PHASE_2_STREAMING.md       # API design
├── PHASE_3_NEON.md           # ARM optimization
├── BENCHMARK_RESULTS.md       # Performance data
└── LESSONS_LEARNED.md         # Retrospective
```

---

## Next Steps (Week 1)

### Immediate Actions

**Today** (Nov 7):
- [x] Create experiment directory structure
- [x] Write PROPOSAL.md
- [ ] Begin BAM specification study
- [ ] Download test BAM files

**This Week**:
- [ ] Complete specification analysis
- [ ] Study noodles source
- [ ] Profile noodles on ARM
- [ ] Design architecture
- [ ] Write Phase 0 deliverables
- [ ] Make GO/NO-GO decision

### Preparation

**Get Test Data**:
```bash
# Download diverse BAM files
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam

# PacBio long reads
# ONT nanopore
# Different aligners (bwa, minimap2, etc.)
```

**Set Up Profiling**:
```bash
cargo install cargo-flamegraph
# macOS Instruments (already have)
```

---

## Commitment

**I commit to**:
- Following evidence-based approach (profiling before optimizing)
- Decision points every 1-2 weeks (can pivot if needed)
- Comprehensive documentation (learning for community)
- Honest assessment (will NO-GO if evidence doesn't support)

**Success means**:
- Building something ARM-native and excellent
- Deep understanding of alignment formats
- Platform for future bioinformatics innovation
- Full-stack biometal ownership

Let's build something great! 🚀

---

**Status**: Phase 0 Starting
**Next Decision Point**: End of Week 1
**Expected Outcome**: GO to Phase 1 with validated NEON opportunities
