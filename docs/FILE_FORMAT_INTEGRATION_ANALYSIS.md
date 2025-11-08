# File Format Integration: Ecosystem Analysis & Strategy

**Date**: November 7, 2025
**Status**: Strategic Planning
**Goal**: Integrate comprehensive file format support into biometal

---

## Executive Summary

**Problem**: Rust bioinformatics ecosystem is fragmented - excellent individual crates but no unified, streaming-first, ARM-optimized solution

**Opportunity**: biometal can become the **unified bioinformatics primitives library** by:
1. Wrapping existing mature crates (noodles, needletail)
2. Providing consistent streaming-first API
3. Adding ARM NEON optimization where evidence supports it
4. Maintaining our evidence-based optimization methodology

**Recommendation**: **Hybrid Strategy** - Depend on mature crates, optimize selectively based on evidence

---

## Current Rust Bioinformatics Ecosystem

### Tier 1: Mature, Production-Ready

| Crate | Formats | Quality | Pure Rust | Streaming | Maintained |
|-------|---------|---------|-----------|-----------|------------|
| **noodles** | BAM, SAM, CRAM, VCF, BCF, GFF, GTF, BED, FASTA, FASTQ | ⭐⭐⭐⭐⭐ | ✅ | ✅ | ✅ Active |
| **needletail** | FASTA, FASTQ (+ compression) | ⭐⭐⭐⭐⭐ | ✅ | ✅ | ✅ Active |
| **rust-htslib** | BAM, SAM, CRAM, VCF, BCF, tabix | ⭐⭐⭐⭐⭐ | ❌ (C bindings) | ✅ | ✅ Active |

### Tier 2: Specialized/Emerging

| Crate | Focus | Quality | Notes |
|-------|-------|---------|-------|
| **paraseq** | Parallel FASTX processing | ⭐⭐⭐⭐ | Map-reduce pattern, adaptive buffering |
| **seq_io** | Fast FASTA/FASTQ parsing | ⭐⭐⭐⭐ | Zero-copy parsing |
| **rust-bio** | Algorithms library | ⭐⭐⭐⭐⭐ | Not I/O focused, uses rust-htslib |
| **bigtools** | BigWig/BigBed | ⭐⭐⭐ | Specialized format |

### Analysis: What's Missing?

**Strengths of ecosystem**:
- ✅ Comprehensive format coverage (noodles)
- ✅ Production-quality parsers
- ✅ Active maintenance
- ✅ Good performance

**Gaps (biometal opportunity)**:
- ❌ **No unified API** across formats
- ❌ **No ARM NEON optimization** (all x86_64 focused)
- ❌ **No evidence-based optimization** (performance not systematically measured)
- ❌ **No streaming-first philosophy** (some support streaming, but not unified)
- ❌ **Fragmentation** (users must learn multiple APIs)

---

## Format-by-Format Assessment

### Sequence Formats

#### FASTA
- **Existing**: biometal ✅, noodles, needletail, seq_io
- **Status**: biometal has good implementation
- **Action**: **Keep current**, potentially benchmark vs needletail
- **Optimization potential**: Medium (text parsing could use NEON)

#### FASTQ
- **Existing**: biometal ✅, noodles, needletail, paraseq
- **Status**: biometal has streaming implementation
- **Action**: **Keep current**, add quality score operations
- **Optimization potential**: High (quality filtering is NEON-friendly)

#### Compressed FASTX (gzip/bgzip)
- **Existing**: biometal (bgzip ✅), needletail (gzip)
- **Status**: biometal has parallel bgzip (Rule 3: 6.5× speedup)
- **Action**: **Keep current**, document advantage
- **Optimization potential**: Low (already optimized)

### Alignment Formats

#### BAM/SAM/CRAM
- **Complexity**: ⭐⭐⭐⭐⭐ (Very High)
- **Existing**: noodles (pure Rust), rust-htslib (C bindings)
- **Format details**:
  - Binary (BAM/CRAM) and text (SAM)
  - Complex compression (BGZF for BAM, custom for CRAM)
  - Indexing (.bai, .crai)
  - Large specification (~50 pages)
- **Decision**: **Wrapper over noodles** (pure Rust, well-maintained)
- **Rationale**:
  - Too complex to reimplement (5,000-10,000 LOC)
  - noodles is production-quality
  - Focus our effort on value-add (streaming API, operations)
- **Value-add opportunities**:
  - ARM-optimized decompression (Rule 3)
  - Streaming filtering operations
  - Quality-based filtering with NEON
- **Action**: Phase 1 priority

#### PAF (Pairwise Alignment Format)
- **Complexity**: ⭐⭐ (Low)
- **Existing**: Limited Rust support
- **Format details**: Simple tab-delimited text (minimap2 output)
- **Decision**: **Reimplement** (simple format)
- **Optimization potential**: Medium (text parsing)
- **Action**: Phase 2 (if users need it)

### Variant Formats

#### VCF/BCF
- **Complexity**: ⭐⭐⭐⭐ (High)
- **Existing**: noodles, rust-htslib
- **Format details**:
  - Text (VCF) and binary (BCF)
  - Complex header and INFO fields
  - Moderate specification (~30 pages)
- **Decision**: **Wrapper over noodles**
- **Value-add opportunities**:
  - Streaming variant filtering
  - ARM-optimized parsing of INFO fields
- **Action**: Phase 1 priority

### Annotation Formats

#### GFF/GTF
- **Complexity**: ⭐⭐⭐ (Medium)
- **Existing**: noodles (GFF3, GTF 2.2)
- **Format details**: Tab-delimited text with attributes
- **Decision**: **Evaluate** (wrapper vs reimplementation)
- **Optimization potential**: **HIGH** (text parsing is NEON-friendly)
- **Evidence needed**: Benchmark noodles, measure bottlenecks
- **Action**: Phase 2 - potential reimplementation candidate

#### BED
- **Complexity**: ⭐ (Low)
- **Existing**: noodles, multiple others
- **Format details**: Simple tab-delimited text
- **Decision**: **Reimplementation candidate**
- **Optimization potential**: Medium (simple enough to optimize)
- **Action**: Phase 2

### Specialized Formats

#### GenBank/EMBL
- **Complexity**: ⭐⭐⭐ (Medium)
- **Existing**: Limited Rust support
- **Decision**: **Phase 3** (low priority, niche use case)

#### SRA
- **Status**: ✅ Already supported via sra-tools wrapper
- **Action**: Keep current

---

## Integration Strategy

### Strategy: Hybrid Approach

**Philosophy**: "Integrate mature crates, optimize where we add value"

```
┌─────────────────────────────────────────────────────────────┐
│                     biometal unified API                     │
├─────────────────────────────────────────────────────────────┤
│  Streaming-first  │  ARM NEON ops  │  Evidence-based opt   │
├─────────────────┬──────────────────┬───────────────────────┤
│  FASTA/FASTQ    │   BAM/SAM/CRAM   │      VCF/BCF         │
│  (current impl) │  (noodles wrap)  │   (noodles wrap)     │
├─────────────────┼──────────────────┼───────────────────────┤
│    GFF/GTF      │       BED        │     PAF (future)     │
│  (reimplement?) │  (reimplement?)  │   (reimplement?)     │
└─────────────────┴──────────────────┴───────────────────────┘
```

### Decision Framework

For each format, ask:

1. **Complexity**: Can we realistically maintain a parser?
   - Low (BED, PAF): Reimplementation viable
   - Medium (GFF/GTF): Evaluate bottlenecks
   - High (BAM, VCF): Wrapper preferred
   - Very High (CRAM): Wrapper only

2. **Existing quality**: Is there a mature Rust implementation?
   - Yes, pure Rust (noodles): Wrapper by default
   - Yes, C bindings (rust-htslib): Avoid if pure Rust alternative exists
   - No: Reimplementation opportunity

3. **Optimization potential**: Can we add measurable value?
   - Measure bottlenecks (≥5× improvement threshold)
   - ARM NEON applicable? (text parsing, compression, filtering)
   - Streaming architecture benefit? (memory reduction)

4. **User demand**: Do users need this format?
   - Priority 1 (common): BAM, VCF, FASTQ
   - Priority 2 (moderate): GFF, BED
   - Priority 3 (niche): GenBank, EMBL

### Implementation Phases

#### Phase 1: Unified API Layer (2-3 weeks)

**Goal**: Single API for all major formats

**Deliverables**:
1. Core traits:
   ```rust
   pub trait SequenceRecord {
       fn id(&self) -> &str;
       fn sequence(&self) -> &[u8];
       fn quality(&self) -> Option<&[u8]>;
   }

   pub trait AlignmentRecord {
       fn id(&self) -> &str;
       fn sequence(&self) -> &[u8];
       fn reference(&self) -> &str;
       fn position(&self) -> u64;
       fn cigar(&self) -> &str;
   }

   pub trait VariantRecord {
       fn chromosome(&self) -> &str;
       fn position(&self) -> u64;
       fn reference_allele(&self) -> &str;
       fn alternate_alleles(&self) -> Vec<&str>;
   }

   pub trait RecordReader: Iterator<Item = Result<Self::Record>> {
       type Record;

       fn from_path(path: impl AsRef<Path>) -> Result<Self>;
       fn from_reader(reader: impl Read) -> Result<Self>;
   }
   ```

2. Wrapper implementations:
   - `biometal::io::bam` (wraps noodles::bam)
   - `biometal::io::vcf` (wraps noodles::vcf)
   - Unified error handling
   - Consistent API patterns

3. Documentation:
   - Migration guide from noodles
   - Performance characteristics
   - When to use each format

4. Validation:
   - Correctness tests (vs noodles directly)
   - Performance benchmarks (overhead < 5%)
   - Memory usage (streaming maintained)

**Success criteria**:
- ✅ User can read BAM, VCF, FASTQ with same pattern
- ✅ Performance within 5% of direct noodles usage
- ✅ Streaming architecture maintained
- ✅ Documentation clear and complete

#### Phase 2: Evidence-Based Optimization (1-2 months)

**Goal**: Add ARM NEON value where measurable

**Process** (for each format):
1. **Profile**: Measure bottlenecks in noodles
2. **Hypothesis**: Identify NEON optimization opportunity
3. **Prototype**: Implement ARM-optimized version
4. **Benchmark**: Compare vs noodles (N=30, rigorous)
5. **Decide**: Ship if ≥5× improvement, else stay with wrapper

**Candidate optimizations**:
1. **GFF/GTF parsing**:
   - Hypothesis: Text parsing benefits from NEON
   - Estimated: 3-5× speedup (text operations)
   - Priority: HIGH (simple format, high usage)

2. **BAM decompression**:
   - Hypothesis: Parallel BGZF (Rule 3: 6.5×)
   - Estimated: 2-3× additional (already have parallel bgzip)
   - Priority: MEDIUM (check if noodles already parallel)

3. **VCF INFO field parsing**:
   - Hypothesis: NEON-accelerated parsing
   - Estimated: 2-4× speedup
   - Priority: MEDIUM (complex format)

4. **BED parsing**:
   - Hypothesis: NEON text parsing
   - Estimated: 5-10× speedup (very simple format)
   - Priority: MEDIUM (if users need it)

**Success criteria**:
- ✅ At least one format with validated ≥5× speedup
- ✅ ASBB-style entry documenting optimization
- ✅ Fallback to noodles if optimization doesn't meet threshold

#### Phase 3: Format-Specific Operations (ongoing)

**Goal**: Cross-format intelligence

**Examples**:
1. **Quality filtering**:
   ```rust
   // FASTQ → filtered FASTQ
   fastq.filter(|r| r.mean_quality() >= 30)

   // BAM → filtered BAM
   bam.filter(|r| r.mapping_quality() >= 20)
   ```

2. **Format conversion**:
   ```rust
   // FASTQ → FASTA
   fastq.to_fasta()

   // SAM → BAM (streaming)
   sam.to_bam(writer)

   // VCF → BED (variants to intervals)
   vcf.to_bed()
   ```

3. **Integration operations**:
   ```rust
   // Annotate FASTQ with GFF features
   fastq.annotate_with(gff, |seq_id| /* lookup */)

   // Filter BAM by BED regions
   bam.filter_by_regions(bed)

   // Extract sequences from reference
   fasta.extract_regions(bed)
   ```

**Success criteria**:
- ✅ User-driven priorities (based on real feedback)
- ✅ Each operation benchmarked and optimized
- ✅ Documentation with real-world examples

---

## Dependency Analysis

### Recommended Dependencies (Phase 1)

```toml
[dependencies]
# Core I/O (existing)
flate2 = "1.0"
memmap2 = "0.9"
crc32fast = "1.4"
thiserror = "1.0"
rayon = "1.8"

# Format support (new)
noodles = { version = "0.80", features = ["bam", "sam", "cram", "vcf", "bcf", "gff", "bed"] }
noodles-core = "0.15"
noodles-bgzf = "0.32"

# Optional: needletail for FASTA/FASTQ benchmarking
needletail = { version = "0.5", optional = true }

# Network streaming (existing)
reqwest = { version = "0.11", features = ["stream"], optional = true }
# ... rest of existing dependencies
```

### Dependency Evaluation

**noodles**:
- ✅ Pure Rust (no C dependencies)
- ✅ Comprehensive (BAM, VCF, GFF, BED, etc.)
- ✅ Well-maintained (active development)
- ✅ Good performance (profiled, optimized)
- ✅ MIT licensed (compatible)
- ⚠️ Not ARM-optimized (opportunity for us)

**rust-htslib** (not recommended):
- ❌ C dependency (htslib)
- ❌ Complicates cross-compilation
- ✅ Production-proven
- Decision: Use noodles instead (pure Rust)

**needletail** (optional, for benchmarking):
- ✅ Fast FASTA/FASTQ parsing
- ✅ MIT licensed
- ⚠️ Overlaps with our implementation
- Decision: Optional dependency for comparison benchmarks

---

## API Design Proposal

### Unified Streaming API

```rust
// Core abstraction: All formats implement RecordReader
pub trait RecordReader: Iterator<Item = Result<Self::Record>> {
    type Record;

    fn from_path(path: impl AsRef<Path>) -> Result<Self>;
    fn from_reader(reader: impl Read) -> Result<Self>;
}

// Example: Reading BAM (same pattern as FASTQ)
use biometal::io::bam::BamReader;
use biometal::RecordReader;

let bam = BamReader::from_path("aligned.bam")?;
for record in bam {
    let record = record?;
    println!("Read: {}, Pos: {}", record.id(), record.position());
}

// Example: Reading VCF (same pattern)
use biometal::io::vcf::VcfReader;

let vcf = VcfReader::from_path("variants.vcf")?;
for variant in vcf {
    let variant = variant?;
    println!("{}:{} {} -> {}",
        variant.chromosome(),
        variant.position(),
        variant.reference_allele(),
        variant.alternate_alleles().join(",")
    );
}

// Example: Network streaming (existing pattern extended)
use biometal::io::bam::BamReader;

let bam = BamReader::from_url("https://example.com/sample.bam")?;
// Same streaming API, network-aware
```

### Consistency with Existing API

**Current biometal API** (FASTQ):
```rust
use biometal::FastqStream;

let stream = FastqStream::from_path("reads.fq.gz")?;
for record in stream {
    let record = record?;
    // Process record
}
```

**Proposed unified API** (all formats):
```rust
use biometal::{FastqReader, BamReader, VcfReader};

// Same pattern for all formats
let fastq = FastqReader::from_path("reads.fq.gz")?;
let bam = BamReader::from_path("aligned.bam")?;
let vcf = VcfReader::from_path("variants.vcf")?;

// All are iterators with same error handling
for record in fastq {
    let record = record?;
    // Format-specific operations on record
}
```

---

## Performance Targets

### Baseline Performance (Phase 1)

**Goal**: Match existing implementations

| Format | Baseline (noodles) | Target (biometal wrapper) | Overhead |
|--------|-------------------|---------------------------|----------|
| BAM    | ~200 Mbp/s       | ≥190 Mbp/s               | <5%      |
| VCF    | ~150 Mbp/s       | ≥140 Mbp/s               | <5%      |
| GFF    | ~100 Mbp/s       | ≥95 Mbp/s                | <5%      |

**Validation**: N=30 benchmarks for each format

### Optimization Targets (Phase 2)

**Goal**: ≥5× improvement from ARM NEON

| Format | Current | Target (ARM NEON) | Evidence Required |
|--------|---------|-------------------|-------------------|
| GFF    | 100 Mbp/s | ≥500 Mbp/s | Validate text parsing speedup |
| BED    | 150 Mbp/s | ≥750 Mbp/s | Simple format, high potential |
| VCF    | 150 Mbp/s | ≥750 Mbp/s | INFO field parsing |

**Process**: ASBB-style entries for each optimization

---

## Risk Assessment

### Technical Risks

**Risk 1: Wrapper overhead > 5%**
- **Likelihood**: Low (minimal abstraction)
- **Impact**: Medium (user preference for noodles)
- **Mitigation**: Benchmark early, optimize hot paths
- **Fallback**: Document when to use noodles directly

**Risk 2: noodles API changes**
- **Likelihood**: Medium (active development)
- **Impact**: Medium (maintenance burden)
- **Mitigation**: Pin versions, test in CI
- **Fallback**: Multiple version support

**Risk 3: ARM optimization doesn't meet ≥5× threshold**
- **Likelihood**: Medium (text parsing varies)
- **Impact**: Low (stay with wrapper)
- **Mitigation**: Evidence-based approach (measure first)
- **Fallback**: Document why wrapper is optimal

**Risk 4: Scope creep**
- **Likelihood**: High (many formats exist)
- **Impact**: High (dilutes focus)
- **Mitigation**: Phased approach, user-driven priorities
- **Fallback**: Prioritize top 5 formats

### Strategic Risks

**Risk 5: Fragmentation continues despite our effort**
- **Likelihood**: Medium (ecosystem momentum)
- **Impact**: Medium (adoption slower)
- **Mitigation**: Excellent docs, clear value proposition
- **Fallback**: Still useful for our own projects

**Risk 6: Users prefer specialized crates**
- **Likelihood**: Medium (existing tools work)
- **Impact**: Low (we still benefit)
- **Mitigation**: Don't force migration, show value
- **Fallback**: Use for our own work

---

## Success Metrics

### Phase 1 (Unified API)
- ✅ ≥5 formats supported with consistent API
- ✅ Documentation with examples for each format
- ✅ Performance within 5% of direct noodles usage
- ✅ All tests passing (correctness vs noodles)
- ✅ PyPI/crates.io published

### Phase 2 (ARM Optimization)
- ✅ At least 1 format with ≥5× NEON speedup
- ✅ ASBB entry documenting optimization
- ✅ Cross-platform testing (Mac ARM, Graviton, x86_64)
- ✅ Publication-quality evidence

### Phase 3 (Format Operations)
- ✅ ≥3 cross-format operations implemented
- ✅ Real-world usage examples
- ✅ User feedback incorporated

### Long-term (6-12 months)
- ✅ 100+ PyPI downloads/month
- ✅ 5+ GitHub issues/PRs from community
- ✅ 1+ publication citing biometal
- ✅ Positive user testimonials

---

## Recommendation

### ✅ **PROCEED with Hybrid Strategy**

**Phase 1 priorities** (next 2-3 weeks):
1. Integrate noodles for BAM, VCF, GFF
2. Create unified `RecordReader` trait
3. Wrap existing FASTQ implementation
4. Documentation and examples
5. Benchmark vs direct noodles usage

**Evidence-based approach**:
- Measure before optimizing
- ≥5× threshold for ARM NEON work
- User feedback drives format priorities
- ASBB-style documentation for all optimizations

**Why this is right**:
- ✅ Leverages mature ecosystem (noodles)
- ✅ Focuses our effort on value-add (ARM, streaming, unification)
- ✅ Follows lessons from SIMD integration (integrate > reimplement)
- ✅ Maintains evidence-based methodology
- ✅ Phased approach manages risk

**What we bring to ecosystem**:
1. **Unified API** (one pattern for all formats)
2. **ARM NEON optimization** (where evidence supports it)
3. **Streaming-first** (constant memory, consistent design)
4. **Evidence-based** (rigorous benchmarks, documented trade-offs)
5. **Cross-format operations** (annotation, filtering, conversion)

---

## Next Steps

### Immediate (This Week)

1. **Create experiment directory**:
   ```bash
   mkdir -p experiments/format-integration
   ```

2. **Write proposal**:
   - PROPOSAL.md (hypothesis, timeline, success criteria)
   - Target formats (Phase 1)
   - Integration approach

3. **Prototype noodles wrapper**:
   - Single format (BAM)
   - Measure overhead
   - Validate approach

4. **Decision point**:
   - If overhead < 5%: Proceed with Phase 1
   - If overhead > 5%: Investigate optimization
   - Document findings

### Week 1-2 (Phase 1 Implementation)

- [ ] Add noodles dependency
- [ ] Create `RecordReader` trait
- [ ] Implement BAM wrapper
- [ ] Implement VCF wrapper
- [ ] Implement GFF wrapper
- [ ] Documentation and examples
- [ ] Benchmarks and tests

### Week 3-4 (Phase 1 Validation)

- [ ] Cross-platform testing
- [ ] Performance validation (N=30)
- [ ] User documentation
- [ ] Integration examples
- [ ] PyPI/crates.io release

### Month 2-3 (Phase 2 - if Phase 1 succeeds)

- [ ] Profile GFF parsing
- [ ] Prototype NEON optimization
- [ ] Benchmark and validate
- [ ] ASBB entry if ≥5× achieved

---

**Status**: Ready for experimental validation
**Timeline**: 2-3 weeks for Phase 1, 1-2 months for Phase 2
**Confidence**: HIGH (following proven hybrid approach from SIMD integration)
