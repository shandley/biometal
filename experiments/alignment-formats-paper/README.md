# Dual-Format Strategy Paper

**Title**: "Dual-Format Strategy for ARM-Native Bioinformatics: BAM and CAF"
**Target**: Bioinformatics (Oxford) or BMC Bioinformatics
**Timeline**: Submit Q2 2026 (~6 months from now)
**Status**: Planning phase (Week 1 of 24)

---

## Overview

This directory contains all materials for the dual-format alignment paper, combining:
1. **ARM-native BAM** (2-3× speedup, spec-compliant)
2. **CAF** (5-10× speedup, columnar format)

**Key Innovation**: First comparison of compatibility-first vs performance-first alignment formats optimized for ARM SIMD.

---

## Timeline

```
Week 1-12:   BAM implementation
Week 13-18:  CAF implementation
Week 19-20:  Comprehensive benchmarking ← CRITICAL FOR PAPER
Week 21-24:  Paper writing & submission
```

**Submission target**: May-June 2026

---

## Directory Structure

```
alignment-formats-paper/
├── README.md                    # This file
├── PAPER_OUTLINE.md             # Complete paper structure
├── BENCHMARKING_PLAN.md         # Data collection strategy
├── benchmarks/                  # (Created during weeks 19-20)
│   ├── raw_data/               # Criterion output
│   ├── processed/              # CSV for figures/tables
│   ├── figures/                # Generated figures (PDF/PNG)
│   └── scripts/                # Analysis scripts (Python/R)
├── draft/                       # (Created during weeks 21-24)
│   ├── manuscript.tex          # LaTeX source
│   ├── figures/                # Final figures
│   ├── tables/                 # Final tables
│   └── supplementary/          # Supplementary material
└── submission/                  # (Week 24)
    ├── manuscript.pdf
    ├── cover_letter.txt
    └── supplementary.pdf
```

---

## Paper Structure (4,000 words)

### Sections

**Abstract** (250 words)
- Problem: BAM designed for 2009 constraints
- Solution: Dual-format strategy
- Results: 2-3× (BAM), 5-10× (CAF) speedups
- Impact: Validated ARM optimization for bioinformatics

**Introduction** (800 words)
- Historical context (BAM in 2009)
- Modern hardware (ARM SIMD, cheap storage)
- Our contribution (dual-format + decision framework)

**Methods** (1,200 words)
- ARM-native BAM (NEON sequence decoding, parallel BGZF)
- CAF design (columnar, pre-decoded, modern compression)
- Implementation (Rust, NEON intrinsics)
- Validation (differential testing, diverse datasets)

**Results** (1,000 words)
- Benchmark suite (parse, filter, count, etc.)
- Performance tables/figures
- Storage analysis
- Real-world workflows

**Discussion** (800 words)
- When to use each format
- Design principles
- Limitations
- Future work

**Conclusions** (200 words)
- Summary and impact

---

## Key Results (Projected)

### Performance Targets

| Metric | BAM (native) | CAF | Speedup |
|--------|--------------|-----|---------|
| Parse throughput | 87 Melem/s | 385 Melem/s | 2.2× / 9.8× |
| Quality filter | - | - | 22× (CAF vs noodles) |
| Base counting | - | - | 21× (CAF vs noodles) |
| Storage overhead | 1.0× | 1.5× | Acceptable |

### Datasets (Public)

1. **HG00096** (1000 Genomes, Illumina)
2. **HG002 PacBio** (GIAB, long reads)
3. **HG002 ONT** (GIAB, ultra-long)
4. **Synthetic** (edge cases)

---

## Benchmarking Plan (Weeks 19-20)

### Core Experiments

1. **Parse throughput** (noodles vs BAM vs CAF)
2. **Quality filtering** (Q20, Q30, Q40)
3. **Base counting** (A/C/G/T/N)
4. **MAPQ filtering** (>30)
5. **Region extraction** (chr:start-end)
6. **Real-world pipeline** (complete QC workflow)
7. **Scaling analysis** (10K → 10M records)
8. **Memory usage** (constant ~5 MB)
9. **Storage analysis** (size comparisons)
10. **Conversion performance** (BAM ↔ CAF)

### Statistical rigor

- **N=30 repetitions** (criterion default)
- **95% confidence intervals**
- **Significance testing** (p < 0.05)
- **Reproducibility** (complete scripts + data)

---

## Figure Plan (4-5 main figures)

**Figure 1**: Format architecture
- Panel A: BAM row-oriented
- Panel B: CAF columnar
- Panel C: NEON vectorization

**Figure 2**: Performance benchmarks
- Panel A: Parse throughput
- Panel B: Operation speedups (bar chart)
- Panel C: Scaling analysis

**Figure 3**: Storage vs performance trade-off
- Scatter plot: overhead vs speedup

**Figure 4**: Real-world workflow
- Timeline comparison (BAM vs CAF)

---

## Table Plan (3-4 main tables)

**Table 1**: Benchmark datasets
- Dataset, Source, Records, Technology, Size

**Table 2**: Performance comparison
- Operation, BAM, CAF, Speedup

**Table 3**: Storage analysis
- Format, Size, Ratio, Compression

---

## Writing Timeline (Weeks 21-24)

### Week 21-22: Drafting

**Day 1-2**: Introduction + Methods
**Day 3-4**: Results (tables + figures)
**Day 5-6**: Discussion + Conclusions
**Day 7**: Abstract + polish

### Week 23: Internal Review

**Day 1-3**: Self-review, revisions
**Day 4-5**: Peer review (lab colleagues)
**Day 6-7**: Incorporate feedback

### Week 24: Submission

**Day 1-2**: Final polish
**Day 3-4**: Format for journal (LaTeX)
**Day 5**: Cover letter, supplementary
**Day 6-7**: Submit!

---

## Submission Checklist

### Main manuscript

- [ ] Title, abstract, keywords
- [ ] Introduction (800 words, ~4 citations)
- [ ] Methods (1,200 words, implementation details)
- [ ] Results (1,000 words, 4 figures, 3 tables)
- [ ] Discussion (800 words)
- [ ] Conclusions (200 words)
- [ ] References (BibTeX)
- [ ] Author contributions
- [ ] Acknowledgments
- [ ] Data availability statement

### Figures

- [ ] Figure 1 (architecture diagram)
- [ ] Figure 2 (performance benchmarks)
- [ ] Figure 3 (trade-offs)
- [ ] Figure 4 (workflows)
- [ ] All high-resolution (300 DPI)
- [ ] Colorblind-friendly palettes

### Tables

- [ ] Table 1 (datasets)
- [ ] Table 2 (performance)
- [ ] Table 3 (storage)
- [ ] LaTeX formatted

### Supplementary

- [ ] Extended benchmarks
- [ ] Additional datasets
- [ ] Reproducibility guide
- [ ] Source code links
- [ ] Zenodo archive DOI

### Administrative

- [ ] Cover letter
- [ ] Suggested reviewers (3-5)
- [ ] Competing interests statement
- [ ] ORCID IDs
- [ ] Funding acknowledgment

---

## Target Journals

### Primary: Bioinformatics (Oxford)

**Pros**:
- High impact (IF: 5.8)
- Methods focus
- Computational biology audience
- Fast review (~6-8 weeks)

**Cons**:
- Competitive
- Application notes format (shorter)

**Fit**: Excellent (methods + benchmarks)

### Secondary: BMC Bioinformatics

**Pros**:
- Open access
- Software/methods focus
- Generous word count
- Reproducibility emphasis

**Cons**:
- Lower impact (IF: 3.2)
- APC (~$2,500)

**Fit**: Excellent (software + open source)

### Tertiary: GigaScience

**Pros**:
- Data-intensive research
- Reproducibility focus
- Open access

**Cons**:
- APC (~$1,600)

---

## Success Criteria

### For publication

**Must have**:
- ✅ 2-3× BAM speedup (ARM NEON)
- ✅ 5-10× CAF speedup (columnar)
- ✅ Lossless BAM ↔ CAF conversion
- ✅ Comprehensive benchmarks (N=30)
- ✅ Public datasets (reproducibility)
- ✅ Open-source implementation

**Nice to have**:
- ⭐ Community adoption (citations)
- ⭐ Tool integration (samtools, etc.)
- ⭐ Follow-up papers (GPU, etc.)

---

## Anticipated Reviewer Questions

**Q1**: "Why not just improve BAM?"
**A**: We did! That's our ARM-native BAM (2-3× speedup). CAF shows what's possible without backward compatibility constraints.

**Q2**: "Will anyone use CAF?"
**A**: CAF demonstrates design principles. Even if adoption is low, the evidence-based approach to format innovation has value.

**Q3**: "Is ARM speedup portable?"
**A**: Yes, we include x86 fallback and validate on multiple platforms (Mac, Graviton, Intel).

**Q4**: "What about samtools?"
**A**: We include samtools in baseline comparison (C implementation).

**Q5**: "Can we reproduce this?"
**A**: Yes! Complete code, datasets, scripts, and Zenodo archive with DOI.

---

## Related Files

**Implementation**:
- `../native-bam-implementation/` - BAM experiments
- `../columnar-alignment-format/` - CAF experiments
- `../ALIGNMENT_FORMATS_ROADMAP.md` - Strategic overview

**Code** (future):
- `../../src/io/bam/` - BAM implementation
- `../../src/io/caf/` - CAF implementation
- `../../benches/` - Benchmark suite

---

## Current Status (Week 1)

**Completed**:
- ✅ Paper outline (PAPER_OUTLINE.md)
- ✅ Benchmarking plan (BENCHMARKING_PLAN.md)
- ✅ Timeline defined (24 weeks)

**In progress**:
- ⏳ BAM Phase 0 (research)
- ⏳ Specification analysis
- ⏳ noodles profiling

**Next milestone**: BAM Phase 1 (minimal reader, Week 2-3)

---

**Date**: November 8, 2025
**Expected submission**: May-June 2026
**Status**: Planning complete, implementation phase beginning

Let's build and publish something great! 🚀
