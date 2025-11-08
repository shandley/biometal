# Columnar Alignment Format (CAF): Proposal

**Date**: November 8, 2025
**Status**: Design Phase (Start after BAM Phase 3-4)
**Type**: Innovation Experiment
**Expected Duration**: 6 weeks

---

## Vision

Build the world's first **ARM-native, columnar alignment format** optimized for 2025+ hardware (NEON, GPU, Neural Engine).

**Strategic Goal**: Pioneer modern bioinformatics file formats that leverage ARM SIMD, not constrained by 2009-era design (BAM's 4-bit encoding, gzip compression).

---

## Hypothesis

A columnar alignment format can achieve:
1. **Performance**: 5-10× faster than BAM for analytical operations
2. **Modern Hardware**: Full utilization of ARM NEON, GPU potential
3. **Usability**: Seamless BAM ↔ CAF conversion (lossless)
4. **Innovation**: Inspire ecosystem to rethink format design

---

## Strategic Rationale

### Why CAF?

**1. BAM's Limitations (2009 Design)**:
- 4-bit sequence encoding: Burns CPU cycles on every read (unpacking overhead)
- gzip compression: Slower than modern codecs (zstd, lz4)
- Row-oriented: Poor SIMD utilization (variable-length records)
- Random access: Optimized for genome browsers, not batch processing

**2. Modern Hardware Opportunities (2025)**:
- ARM NEON: 128-bit SIMD (16 operations/cycle)
- Apple Silicon: Neural Engine (11-15 TOPS), GPU (thousands of threads)
- NVMe storage: 7 GB/s (storage cheap, CPU expensive)
- Modern codecs: zstd (2-3× faster decompression than gzip)

**3. biometal's Differentiation**:
- "First ARM-native bioinformatics toolkit"
- Full-stack ownership (FASTQ → operations → BAM → CAF)
- Platform for format innovation
- Research tool (not constrained by clinical compliance)

### Why Now?

**Favorable Conditions**:
- ✅ BAM implementation underway (reference for compatibility)
- ✅ Proven NEON expertise (16-25× speedups in biometal)
- ✅ No rush to publish (time to innovate)
- ✅ Research/experimental context (okay if community doesn't adopt immediately)

**Risk Mitigation**:
- BAM remains primary format (compatibility)
- CAF is additive (not replacement)
- Lossless conversion (BAM ↔ CAF)
- Clear documentation of trade-offs

---

## Scope

### In Scope

**Phase 1: Core Format** (2 weeks)
- ✅ Columnar block layout (10,000 records)
- ✅ BAM → CAF conversion (lossless)
- ✅ CAF → BAM conversion (lossless)
- ✅ zstd/lz4 compression integration
- ✅ Basic reader/writer API

**Phase 2: NEON Optimization** (2 weeks)
- ✅ NEON quality filtering (16-25× target)
- ✅ NEON base counting (16-25× target)
- ✅ NEON MAPQ filtering (16× target)
- ✅ Parallel block processing (rayon)

**Phase 3: Indexing** (1 week)
- ✅ Block-level index (.caf.idx)
- ✅ Region queries (chr:start-end)
- ✅ Streaming API

**Phase 4: Production Polish** (1 week)
- ✅ Comprehensive testing (200+ tests)
- ✅ Documentation and examples
- ✅ Benchmarks (CAF vs BAM)
- ✅ CLI tools (convert, query, stats)

### Out of Scope (Initially)

**Deferred to future versions**:
- ❌ GPU integration (Metal/CUDA) - Phase 5 (optional)
- ❌ Neural Engine quantization - Phase 5 (optional)
- ❌ Base-level random access - Not designed for genome browsers
- ❌ Writing tools (e.g., aligner output) - Conversion from BAM sufficient
- ❌ Clinical pipeline validation - Research tool

---

## Success Criteria

### Quantitative

| Phase | Metric | Target | Threshold |
|-------|--------|--------|-----------|
| **Phase 1** | Round-trip correctness | 100% vs BAM | Must match |
| **Phase 2** | Quality filter speedup | ≥20× | GO if ≥10× |
| **Phase 2** | Base count speedup | ≥20× | GO if ≥10× |
| **Phase 2** | Overall speedup | ≥5× | GO if ≥3× |
| **Phase 3** | Index query overhead | <10% | Index size/access |
| **Phase 4** | Test coverage | 200+ tests | ≥95% coverage |

### Qualitative

**Must Achieve**:
- ✅ Correct: Lossless BAM ↔ CAF conversion
- ✅ Fast: Measurably faster on ARM (5-10× target)
- ✅ Usable: Clean API, good docs
- ✅ Documented: Trade-offs clearly explained

**Success Even If**:
- ⚠️ Community doesn't adopt immediately (experimental tool)
- ⚠️ Storage overhead is 2× (acceptable trade-off for 5-10× speed)
- ⚠️ Some features deferred (GPU, Neural Engine)

---

## Timeline & Phases

### Overview (6 weeks, after BAM Phase 3-4)

```
Weeks 1-2:  Phase 1 - Core Format
Weeks 3-4:  Phase 2 - NEON Optimization
Week 5:     Phase 3 - Indexing and Queries
Week 6:     Phase 4 - Production Polish
```

### Decision Points

| Week | Decision | GO Criteria | NO-GO Action |
|------|----------|-------------|--------------|
| 2 | Continue to Phase 2? | Lossless conversion works | Debug or pivot |
| 4 | Continue to Phase 3? | ≥3× speedup achieved | Ship as-is, defer indexing |
| 5 | Continue to Phase 4? | Indexing works | Ship without index |
| 6 | Ship it? | All must-haves met | Extended polish |

---

## Phase Objectives

### Phase 1: Core Format (Weeks 1-2)

**Goals**:
1. Design and implement columnar block layout
2. BAM → CAF converter (lossless)
3. CAF → BAM converter (lossless)
4. Basic reader API

**Deliverables**:
- [ ] `src/io/caf/writer.rs` - CAF writer
- [ ] `src/io/caf/reader.rs` - CAF reader
- [ ] `src/io/caf/block.rs` - Columnar block structure
- [ ] `src/io/caf/compression.rs` - zstd/lz4 integration
- [ ] `examples/caf_convert.rs` - Conversion example
- [ ] Tests: Round-trip correctness (100+ files)

**Validation**:
```rust
// Differential testing
#[test]
fn test_roundtrip_lossless() {
    let original_bam = read_bam("test.bam")?;

    // BAM → CAF → BAM
    convert_bam_to_caf("test.bam", "temp.caf")?;
    convert_caf_to_bam("temp.caf", "roundtrip.bam")?;

    let roundtrip_bam = read_bam("roundtrip.bam")?;

    // Must match exactly
    assert_eq!(original_bam, roundtrip_bam);
}
```

### Phase 2: NEON Optimization (Weeks 3-4)

**Goals**:
1. Implement NEON quality filtering
2. Implement NEON base counting
3. Implement NEON MAPQ filtering
4. Benchmark CAF vs BAM

**Deliverables**:
- [ ] `src/operations/caf_quality.rs` - CAF quality ops (NEON)
- [ ] `src/operations/caf_base_counting.rs` - CAF base counting (NEON)
- [ ] `src/operations/caf_filtering.rs` - CAF filtering (NEON)
- [ ] `benches/caf_vs_bam.rs` - Comprehensive benchmark
- [ ] `experiments/columnar-alignment-format/BENCHMARKS.md` - Results

**Target Benchmark**:
```
Operation: Quality filter (Q30, 100K records)
  BAM (noodles):     2.00 sec  (baseline)
  CAF (NEON):        0.08 sec  (25× faster)

Operation: Base counting (100K records)
  BAM (noodles):     1.50 sec  (baseline)
  CAF (NEON):        0.06 sec  (25× faster)

Operation: MAPQ filter (>30, 100K records)
  BAM (noodles):     0.50 sec  (baseline)
  CAF (NEON):        0.03 sec  (16× faster)

Overall parse + filter:
  BAM:               2.56 sec  (baseline)
  CAF:               0.25 sec  (10× faster)
```

**GO Criteria**: ≥3× overall speedup (target: 5-10×)

### Phase 3: Indexing (Week 5)

**Goals**:
1. Block-level index (.caf.idx)
2. Region query support
3. Streaming query API

**Deliverables**:
- [ ] `src/io/caf/index.rs` - Indexing module
- [ ] `src/io/caf/query.rs` - Query API
- [ ] CLI: `biometal query --region chr1:1000-2000`
- [ ] Tests: Query correctness

**Index Structure**:
```rust
pub struct CafIndex {
    block_offsets: Vec<u64>,        // File offset per block
    block_metadata: Vec<BlockMeta>,  // Genomic range per block
}

pub struct BlockMeta {
    ref_id: i32,
    start_pos: i32,
    end_pos: i32,
    num_records: u32,
}
```

**Query Performance**:
- Target: <10% overhead vs full scan
- Granularity: Block-level (10,000 records)

### Phase 4: Production Polish (Week 6)

**Goals**:
1. Comprehensive testing (edge cases, property tests)
2. Documentation (guides, examples, benchmarks)
3. CLI tools (convert, query, stats, validate)
4. Published benchmark results

**Deliverables**:
- [ ] 200+ tests (unit, integration, property)
- [ ] `docs/CAF_USER_GUIDE.md` - User documentation
- [ ] `experiments/columnar-alignment-format/BENCHMARKS.md` - Published results
- [ ] CLI: `biometal convert`, `biometal query`, `biometal stats`
- [ ] Examples: Real-world workflows

**Documentation Requirements**:
- Clear trade-offs (storage vs speed)
- Compatibility notes (CAF ↔ BAM)
- Performance characteristics
- When to use CAF vs BAM

---

## Resource Requirements

### Time

- **Full-time equivalent**: 6 weeks (after BAM complete)
- **Part-time (50%)**: 12 weeks
- **Decision points**: Every 1-2 weeks

### Infrastructure

**Required**:
- ✅ Mac ARM (M1/M2/M3/M4) - you have
- ✅ Rust toolchain - you have
- ✅ BAM implementation - in progress

**Needed**:
- [ ] Test BAM files (diverse sources) - from BAM experiment
- [ ] Benchmarking harness (criterion) - have
- [ ] Large datasets (100K-1M records) - download

### Knowledge

**Have**:
- ✅ Rust expertise
- ✅ ARM NEON experience (proven in biometal)
- ✅ BAM format knowledge (from Phase 0-3)
- ✅ Columnar data structures (common pattern)

**Acquiring**:
- [ ] zstd/lz4 compression libraries
- [ ] Indexing strategies (B-trees, block metadata)

---

## Risk Assessment

### Technical Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| NEON speedup <5× | Low | Medium | BAM already gives 2× NEON floor |
| Round-trip not lossless | Low | High | Comprehensive differential testing |
| Storage overhead >2× | Low | Low | Acceptable trade-off (documented) |
| Indexing complex | Medium | Low | Defer to Phase 5 if needed |

### Strategic Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Community doesn't adopt | High | Low | Experimental tool, research use case |
| Users confused (CAF vs BAM) | Medium | Low | Clear documentation, trade-offs |
| Maintenance burden | Low | Low | Well-tested, clear architecture |

---

## Comparison to Alternatives

### vs BAM (Row-Oriented)

| Aspect | BAM | CAF |
|--------|-----|-----|
| **Storage** | 1.0× (baseline) | 1.5-2.0× |
| **Parse speed** | 1.0× (39 Melem/s) | 5-10× faster |
| **NEON optimization** | Limited (4-bit decode) | Extensive (columnar) |
| **Compression** | gzip (1992) | zstd/lz4 (modern) |
| **Random access** | Base-level (BAI) | Block-level (10K) |
| **Compatibility** | Universal | biometal only (initially) |
| **Use case** | General purpose | Analytical queries |

**Decision**: CAF complements BAM, not replaces

### vs Parquet (Columnar, Generic)

**Apache Parquet**: Generic columnar format
- Pros: Mature, tool support, Spark/Hadoop integration
- Cons: Not bioinformatics-specific, no NEON optimization, complex spec

**CAF**:
- Pros: Bioinformatics-specific, ARM NEON optimized, simple
- Cons: New format, no tool support (yet)

**Decision**: CAF is domain-specific, optimized for ARM bioinformatics

---

## Success Definition

### Minimum Viable Success (MVP)

At end of Phase 4, we have:
- ✅ Lossless BAM ↔ CAF conversion
- ✅ ≥3× speedup for analytical operations
- ✅ Comprehensive tests and documentation
- ✅ CLI tools (convert, query, stats)
- ✅ Published benchmarks

**Even if**:
- Community doesn't adopt (experimental research tool)
- Storage is 2× larger (acceptable for 5-10× speed)
- Some features deferred (GPU, Neural Engine)

### Stretch Goals

- ⭐ ≥10× speedup (exceptional performance)
- ⭐ GPU integration (Metal kernels)
- ⭐ Neural Engine quantization
- ⭐ Community adoption (others implement CAF readers)
- ⭐ Published paper (new format design for modern hardware)

---

## Experiment Methodology

### Validation Approach

**Correctness** (Phase 1):
```rust
// Differential testing vs BAM
proptest! {
    #[test]
    fn roundtrip_preserves_data(bam_file in arb_bam_file()) {
        let original = read_bam(&bam_file)?;
        let caf = convert_to_caf(&bam_file)?;
        let roundtrip = convert_to_bam(&caf)?;
        prop_assert_eq!(original, roundtrip);
    }
}
```

**Performance** (Phase 2):
```rust
// Benchmark: CAF vs BAM
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_quality_filter(c: &mut Criterion) {
    let bam = load_test_bam("100k_records.bam");
    let caf = convert_to_caf(&bam);

    c.bench_function("BAM quality filter", |b| {
        b.iter(|| filter_quality_bam(&bam, 30))
    });

    c.bench_function("CAF quality filter NEON", |b| {
        b.iter(|| filter_quality_caf_neon(&caf, 30))
    });
}
```

---

## Documentation Plan

### Phase 1-4 Outputs

```
experiments/columnar-alignment-format/
├── PROPOSAL.md                # This file
├── CAF_SPECIFICATION.md       # Format design (v1.0)
├── IMPLEMENTATION_LOG.md      # Development diary
├── BENCHMARKS.md              # CAF vs BAM results
└── LESSONS_LEARNED.md         # Retrospective

docs/
└── CAF_USER_GUIDE.md          # User-facing guide

examples/
├── caf_convert.rs             # BAM ↔ CAF conversion
├── caf_query.rs               # Region queries
└── caf_analysis.rs            # Analytical workflows
```

---

## Next Steps (After BAM Phase 3-4)

### Preparation (During BAM Development)

**While implementing BAM**:
- [ ] Finalize CAF specification (this document)
- [ ] Design columnar block structure
- [ ] Research zstd/lz4 Rust crates
- [ ] Plan NEON optimization strategy

### Week 1 (Phase 1 Start)

**Day 1-2**: Core data structures
- [ ] Define CafBlock struct (columnar layout)
- [ ] Implement basic serialization/deserialization
- [ ] zstd compression integration

**Day 3-5**: BAM → CAF converter
- [ ] Read BAM records (using biometal BAM reader)
- [ ] Convert to columnar layout
- [ ] Write CAF blocks

**Day 6-7**: CAF → BAM converter
- [ ] Read CAF blocks
- [ ] Convert to row-oriented BAM
- [ ] Validation: Round-trip testing

---

## Commitment

**I commit to**:
- Following evidence-based approach (benchmark early and often)
- Decision points every 1-2 weeks (can pause/pivot if needed)
- Comprehensive documentation (help community understand CAF)
- Honest assessment (will document if CAF doesn't meet targets)

**Success means**:
- Building truly modern, ARM-native format
- Proving columnar layout benefits for bioinformatics
- Inspiring ecosystem to rethink format design
- Having fun innovating! :)

---

**Status**: Design Phase (November 2025)
**Next Decision Point**: After BAM Phase 3-4 complete
**Expected Start**: ~Week 6-7 from now
**Expected Outcome**: Production-ready CAF v1.0 with 5-10× speedup over BAM
