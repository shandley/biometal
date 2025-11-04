# Experiment Proposal: ARM-Native SRA Decoder

**Status**: Research
**Owner**: Scott Handley
**Date**: 2025-11-05
**Timeline**: 2 weeks (research), +4 weeks (implementation if GO)

---

## 1. Problem Statement

**What problem are we trying to solve?**

NCBI's Sequence Read Archive (SRA) stores petabases of sequencing data in a proprietary binary format that requires the SRA Toolkit to decode. This creates several barriers:

1. **Dependency hell**: SRA Toolkit has complex dependencies and inconsistent behavior across platforms
2. **Performance bottleneck**: Toolkit is not optimized for ARM architecture (no NEON SIMD)
3. **Integration friction**: C++ toolkit difficult to integrate into Rust streaming pipeline
4. **Streaming limitations**: Toolkit often requires full download before processing

**Why is this important?**

SRA is the world's largest public repository of sequencing data. Enabling direct streaming analysis of SRA files without downloads would:

- **Democratize access**: ~5 MB memory instead of downloading 50-200 MB files
- **Speed research**: Immediate analysis without wait time
- **Reduce costs**: No local storage required for LMIC researchers
- **Enable new workflows**: Stream → analyze → discard for large-scale meta-analysis

**Current approaches and their limitations:**

1. **NCBI SRA Toolkit** (C++)
   - Not ARM-optimized (no NEON SIMD)
   - Complex build process
   - Inconsistent behavior across platforms
   - Often requires full download

2. **Wrap SRA Toolkit** (pragmatic but suboptimal)
   - Still has performance limitations
   - External dependency complexity
   - Loses biometal's ARM-native advantages

3. **Convert to FASTQ upstream** (workaround)
   - Defeats purpose of network streaming
   - Requires preprocessing infrastructure
   - Not viable for on-demand analysis

---

## 2. Hypothesis

**Primary Hypothesis:**

ARM NEON-optimized SRA decoder can achieve **≥10× speedup** compared to scalar implementation for critical operations (base unpacking, quality decompression), while maintaining:
- Streaming architecture (constant ~5 MB memory)
- 100% correctness vs SRA Toolkit reference
- Clean integration with biometal's existing infrastructure

**Supporting Evidence:**

From OPTIMIZATION_RULES.md (ASBB project):

1. **Rule 1 (NEON SIMD)**: 16-25× speedup for element-wise operations
   - Complexity 0.30-0.40 operations show best gains
   - Base unpacking (2-bit → 8-bit) fits this profile perfectly
   - Entry 020-025: 307 experiments, 9,210 measurements

2. **Rule 2 (Block-based)**: 10K record blocks preserve SIMD speedup
   - Entry 027: 1,440 measurements
   - Directly applicable to SRA record processing

3. **Rule 5 (Streaming)**: 99.5% memory reduction
   - Entry 026: 720 measurements
   - SRA's columnar format naturally fits streaming

**Literature evidence:**

- Published SRA format specification available
- 2-bit packed bases are embarrassingly parallel
- Quality score decompression is well-understood (deflate/gzip)
- Columnar storage enables selective decompression

---

## 3. Proposed Approach

### 3.1 Technical Design

**Phase 1: Format Understanding & Scalar Baseline**
1. Implement SRA format parser (header, index, data blocks)
2. Scalar implementations of critical operations:
   - 2-bit packed base decoding
   - Quality score decompression (deflate)
   - Read name reconstruction
3. Validate correctness against SRA Toolkit output

**Phase 2: ARM NEON Optimization**
1. NEON-optimized base unpacking (2-bit → 8-bit FASTQ)
2. Parallel quality decompression (if beneficial)
3. Block-based processing (10K records, Rule 2)

**Phase 3: Streaming Integration**
1. Iterator-based API matching biometal's FastqStream
2. Constant memory architecture (~5 MB target)
3. Integration with existing network streaming (HttpReader)

**High-level API:**
```rust
// Should feel identical to existing FASTQ streaming
let source = DataSource::Sra("SRR390728".to_string());
let stream = FastqStream::new(source)?;

for record in stream {
    let record = record?;
    // Native ARM-optimized decoding under the hood
}
```

### 3.2 ARM NEON Opportunities

**Hot Path 1: Base Unpacking (2-bit → 8-bit)**
- **Operation**: Convert packed bases (00=A, 01=C, 10=G, 11=T) to ASCII
- **NEON advantage**: Process 32 bases per instruction (128-bit / 4 bits)
- **Expected speedup**: 16-25× (Rule 1, complexity ~0.35)
- **Evidence**: Similar to GC content calculation (Entry 023)

**Hot Path 2: Quality Score Decompression**
- **Operation**: Inflate compressed quality strings
- **NEON advantage**: Parallel block decompression if custom
- **Expected speedup**: 6.5× (Rule 3, parallel approach)
- **Alternative**: Use existing flate2 (already optimized)

**Hot Path 3: Read Name Reconstruction**
- **Operation**: Build read identifiers from columnar data
- **NEON advantage**: String formatting may benefit
- **Priority**: LOWER (likely I/O dominated)

**Performance Profile (estimated):**
- Base unpacking: 40% of CPU time → 16× speedup → 38% overall improvement
- Quality decompression: 30% of CPU time → already optimized (flate2)
- Parsing/overhead: 30% of CPU time → limited NEON benefit

**Conservative estimate**: 20-40× total speedup vs naive scalar, 10-15× vs optimized scalar

### 3.3 Streaming Architecture

**Memory Architecture:**
- SRA block buffer: ~1 MB (columnar blocks)
- Decompressed buffer: ~2 MB (10K records × 200 bytes)
- Network buffer: ~2 MB (HttpReader existing)
- **Total: ~5 MB constant memory** ✅

**Streaming Properties:**
- ✅ Iterator-based: `impl Iterator<Item = Result<FastqRecord>>`
- ✅ Constant memory: Block-based processing (Rule 2)
- ✅ Zero-copy: Direct slice access where possible
- ✅ Lazy decompression: Only decompress requested data

**Integration with biometal:**
```rust
// src/io/sra_native.rs (new module)
pub struct SraNativeDecoder<R: BufRead> {
    reader: R,
    block_buffer: Vec<u8>,
    record_buffer: Vec<FastqRecord>,
    position: usize,
}

impl<R: BufRead> Iterator for SraNativeDecoder<R> {
    type Item = io::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        // Decode one record, return, maintain constant memory
    }
}
```

---

## 4. Success Criteria

### 4.1 Go/No-Go Thresholds

**✅ GO Decision** (proceed to full implementation):
- **Performance**: ≥10× speedup (NEON vs scalar) on critical operations
- **Memory**: Constant ~5 MB regardless of file size
- **Correctness**: 100% match with SRA Toolkit output (N=1000 random reads)
- **Integration**: Clean API matching existing FastqStream interface
- **Maintainability**: Code is documented and testable

**🤔 MAYBE Decision** (deeper investigation):
- Speedup: 5-10× (marginal benefit, need more analysis)
- Some edge cases fail (need to understand scope)
- Format complexity higher than expected (assess effort)

**❌ NO-GO Decision** (wrap existing solution):
- Speedup: <5× (not worth custom implementation)
- Correctness: Fails validation on common inputs
- Complexity: Format reverse-engineering infeasible
- Maintenance: Would require ongoing format tracking

### 4.2 Validation Methodology

**Benchmarking:**
- **Sample size**: N=30 per configuration (statistical validity)
- **Statistics**: Mean, 95% CI, Cohen's d effect size
- **Baseline comparison**:
  - Scalar implementation (same Rust code, no NEON)
  - SRA Toolkit `fastq-dump` (for context, not primary comparison)
- **Test files**:
  - Small: 1,000 reads (~1 MB)
  - Medium: 100,000 reads (~10 MB)
  - Large: 1,000,000 reads (~100 MB)

**Correctness:**
- **Test against**: SRA Toolkit `fastq-dump` output (gold standard)
- **Property-based testing**:
  - Base sequences use only ACGTN characters
  - Quality lengths match sequence lengths
  - All reads have valid identifiers
  - Packed → unpacked → packed round-trip
- **Test datasets**:
  - SRR390728 (E. coli, well-characterized)
  - SRR1553425 (Human, longer reads)
  - Edge cases (single read, empty file, corrupted data)

---

## 5. Timeline

### Week 1-2: Research Phase (Nov 5-15, 2025)

**Days 1-3**: Format Analysis & Scalar Implementation
- [ ] Study SRA format specification (NCBI docs)
- [ ] Analyze sample SRA files (hex dump, structure)
- [ ] Implement scalar base unpacking (2-bit → 8-bit)
- [ ] Implement quality score decompression
- [ ] Write correctness tests vs SRA Toolkit

**Days 4-8**: NEON Optimization & Benchmarking
- [ ] Implement NEON-optimized base unpacking
- [ ] Profile to identify remaining hot paths
- [ ] Optimize additional hot paths if found
- [ ] Write comprehensive benchmarks (N=30)
- [ ] Property-based testing (proptest)

**Days 9-10**: Validation & Analysis
- [ ] Run full benchmark suite (N=30, multiple file sizes)
- [ ] Statistical analysis (mean, 95% CI, Cohen's d)
- [ ] Correctness validation (1000+ random reads)
- [ ] Document findings in FINDINGS.md
- [ ] Update RESEARCH_LOG.md with daily progress

**Days 11-12**: Decision & Documentation
- [ ] Make Go/No-Go decision based on data
- [ ] Update experiments/.experiments.toml with outcome
- [ ] If GO: Draft integration plan for biometal
- [ ] If NO-GO: Document why and recommend alternatives
- [ ] Share findings with team

### Week 3-6: Implementation Phase (IF GO)

**Week 3**: Core Integration
- [ ] Integrate SraNativeDecoder into biometal
- [ ] Add to DataSource enum
- [ ] Update FastqStream to dispatch to native decoder
- [ ] Comprehensive testing

**Week 4**: Network Streaming
- [ ] Enable SRA streaming over HTTP
- [ ] Range request optimization for SRA blocks
- [ ] Caching strategy for SRA index
- [ ] Performance validation

**Week 5**: Polish & Documentation
- [ ] Examples and documentation
- [ ] Cross-platform testing (Mac, Graviton, RPi)
- [ ] Benchmark results publication
- [ ] Integration testing with real pipelines

**Week 6**: Publication Preparation
- [ ] Draft paper/blog post
- [ ] Prepare benchmark data for publication
- [ ] Code cleanup and release
- [ ] Community feedback

---

## 6. Resources Required

**Development Time:**
- Research: 2 weeks (time-boxed, Nov 5-15)
- Prototype: Included in research phase
- Implementation (if GO): 4 weeks
- Total: 2-6 weeks depending on decision

**Dependencies:**
- `flate2` - Quality score decompression (already in Rust ecosystem)
- `byteorder` - Binary format parsing (lightweight)
- SRA Toolkit - For correctness validation only (not runtime dependency)

**Hardware:**
- Mac M3 (primary development)
- Linux ARM (Graviton) for cross-platform validation
- Optional: Raspberry Pi for low-end ARM validation

**Documentation:**
- [SRA Format Specification](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration)
- "The Sequence Read Archive" (Leinonen et al., 2011)
- NCBI SRA file format documentation
- Existing sra-tools source code for reference

---

## 7. Risks & Mitigation

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Format complexity higher than expected | Medium | High | Time-box research (2 weeks), have NO-GO ready |
| NCBI changes format frequently | Low | Medium | Document format version, detect changes |
| Speedup doesn't meet 10× threshold | Low | Medium | Lower threshold to 5×, or wrap toolkit |
| Correctness validation fails | Low | Critical | Strict testing, property-based tests |
| Integration complexity with biometal | Low | Medium | Early API design, incremental integration |
| Limited SRA format documentation | Medium | Medium | Reverse-engineer from sra-tools source |

**Unknown unknowns:**

- Are there undocumented format variants?
- Does NCBI have additional compression methods?
- Are there edge cases in real-world files not covered by spec?

**Discovery strategy:**

- Analyze multiple SRA files from different sources
- Test against wide variety of organisms and read types
- Monitor SRA Toolkit source for format changes
- Engage with NCBI if format questions arise

---

## 8. Publication Potential

**If successful (≥10× speedup):**
- **Target venue**: Bioinformatics (journal) or ISMB/RECOMB (conference)
- **Angle**:
  - "ARM-Native SRA Decoder: Democratizing Access to Petabases of Sequencing Data"
  - First ARM NEON-optimized implementation
  - Streaming architecture for resource-constrained environments
  - Evidence-based optimization methodology
- **Timeline**: Submit Q1 2026 (after 3-6 months of validation)
- **Broader impact**: Could influence NCBI to provide ARM-optimized toolkit

**If moderately successful (5-10× speedup):**
- **Target**: Preprint (bioRxiv) + blog post
- **Value**: Still useful, documents approach for others
- **Integration**: May still integrate if other benefits (cleaner API, streaming)

**If unsuccessful (<5× speedup):**
- **Value**: Negative results prevent others from wasting time
- **Publication**: Technical report or blog post
- **Lessons learned**:
  - Document why approach didn't work
  - What aspects of SRA format limit optimization
  - Recommendations for future attempts
- **Still valuable**: Shows evidence-based process works (caught dead-end early)

---

## 9. References

**SRA Format Specification:**
- [NCBI SRA File Format Documentation](https://github.com/ncbi/sra-tools/wiki)
- [sra-tools source code](https://github.com/ncbi/sra-tools) (primary reference)
- [SRA Handbook](https://www.ncbi.nlm.nih.gov/books/NBK47540/)

**Related Publications:**
- Leinonen R, et al. (2011) "The Sequence Read Archive" *Nucleic Acids Research* 39:D19-21
- Kodama Y, et al. (2012) "The Sequence Read Archive: explosive growth of sequencing data" *Nucleic Acids Research* 40:D54-56

**ARM NEON Optimization:**
- OPTIMIZATION_RULES.md (biometal, from ASBB project)
- [ARM NEON Intrinsics Reference](https://developer.arm.com/architectures/instruction-sets/intrinsics/)
- Handley S (2025) "Evidence-Based Optimization for Bioinformatics on Apple Silicon" (ASBB lab notebook)

**Existing Implementations:**
- NCBI SRA Toolkit (C++) - gold standard for correctness
- pysradb (Python) - wraps SRA Toolkit
- ngs (Java/C++) - NCBI's next-gen API

---

## 10. Approval

**Proposed by**: Scott Handley (with Claude)
**Date**: 2025-11-05

**Status**: APPROVED for 2-week research phase

**Rationale for approval:**
- Aligns perfectly with biometal's mission (ARM-native, streaming, evidence-based)
- Clear go/no-go criteria prevent wasted effort
- Time-boxed research limits risk
- High potential impact if successful
- Valuable learning even if unsuccessful

**Experiment Management:**
- Registered in experiments/.experiments.toml
- Daily updates to RESEARCH_LOG.md required
- Go/No-Go decision by Nov 15, 2025
- Full findings documented regardless of outcome

---

**Remember**: This is research. Failure is an option. Document everything.

**Success = Learning, not just achieving hypothesis.**
