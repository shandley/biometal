# Experiment Findings: [Experiment Name]

**Experiment**: [Name]
**Owner**: [Your Name]
**Duration**: YYYY-MM-DD to YYYY-MM-DD
**Status**: [Completed / Abandoned]
**Decision**: [GO / NO-GO / MAYBE]

---

## Executive Summary

[2-3 paragraph summary of the entire experiment: what you set out to do, what you found, and what you recommend]

**Key Result**: [One sentence summary]

**Recommendation**: [GO / NO-GO and why]

---

## 1. Original Hypothesis

**Stated Hypothesis**:
[As written in PROPOSAL.md]

**Success Criteria**:
- ✅ GO: [Criteria]
- ❌ NO-GO: [Criteria]

---

## 2. Methodology

### 2.1 Approach

[Describe what you actually did. May differ from proposal based on learnings.]

### 2.2 Benchmarking Setup

**Hardware**:
- CPU: [e.g., Apple M4 Max]
- RAM: [e.g., 128 GB]
- OS: [e.g., macOS 15.1]

**Software**:
- Rust: [version]
- Criterion: [version]
- Baseline tool: [version]

**Methodology**:
- Sample size: N=30
- Warmup: X iterations
- Measurement: Y iterations
- Statistics: Mean, 95% CI, Cohen's d

### 2.3 Correctness Validation

[How did you ensure the implementation was correct?]

---

## 3. Results

### 3.1 Performance Benchmarks

#### Benchmark 1: [Operation Name]

| Implementation | Mean (ms) | Std Dev | 95% CI | Speedup |
|----------------|-----------|---------|--------|---------|
| Baseline (Scalar) | X.XX | ±Y.YY | [A, B] | 1.0× |
| NEON Optimized | X.XX | ±Y.YY | [A, B] | **Z.Z×** |
| SRA Toolkit (comparison) | X.XX | ±Y.YY | [A, B] | W.W× |

**Statistical Significance**:
- Cohen's d: [value] ([interpretation: large/medium/small effect])
- P-value: [if calculated]

**Interpretation**:
[What does this tell us? Did we meet success criteria?]

#### Benchmark 2: [Operation Name]

[Same format]

### 3.2 Memory Usage

| Implementation | Peak Memory | Constant? |
|----------------|-------------|-----------|
| Baseline | X MB | Yes/No |
| NEON | Y MB | Yes/No |
| SRA Toolkit | Z MB | Yes/No |

**Analysis**:
[Does this meet biometal's constant-memory guarantee?]

### 3.3 Correctness

**Test Results**:
- Unit tests: [X/Y passing]
- Integration tests: [X/Y passing]
- Property-based tests: [X cases, Y failures]
- Comparison vs baseline: [Match / Differ]

**Issues Found**:
- [Issue 1] - [How resolved]
- [Issue 2] - [Still open / Won't fix]

---

## 4. Analysis

### 4.1 Success Criteria Met?

**GO Criteria**:
- [ ] Criterion 1: [Met / Not Met - by how much]
- [ ] Criterion 2: [Met / Not Met - by how much]
- [ ] Criterion 3: [Met / Not Met - by how much]

**Overall**: [Did we meet GO criteria?]

### 4.2 Unexpected Findings

[What surprised you? What didn't match expectations?]

1. [Finding 1]
   - Expected: [...]
   - Actual: [...]
   - Explanation: [...]

2. [Finding 2]
   ...

### 4.3 Limitations

[What are the limitations of this implementation/approach?]

1. [Limitation 1]
2. [Limitation 2]
3. [Limitation 3]

---

## 5. Decision Rationale

### 5.1 Quantitative Analysis

[Data-driven reasons for GO/NO-GO decision]

**Achieved vs Target**:
- Speedup: [Achieved: X.X×] vs [Target: Y.Y×]
- Memory: [Achieved: X MB constant] vs [Target: Y MB constant]
- Correctness: [100% match vs baseline]

### 5.2 Qualitative Factors

[Non-quantitative considerations]

**Pros**:
- [Pro 1]
- [Pro 2]

**Cons**:
- [Con 1]
- [Con 2]

**Risk Assessment**:
- Implementation complexity: [High / Medium / Low]
- Maintenance burden: [High / Medium / Low]
- External dependencies: [Many / Few / None]

### 5.3 Final Decision

**Decision**: [✅ GO / ❌ NO-GO / 🤔 MAYBE]

**Reasoning**:
[Detailed explanation of why this decision was made]

---

## 6. Recommendations

### If GO Decision

**Next Steps**:
1. [Step 1]
2. [Step 2]
3. [Step 3]

**Timeline**: [Estimated time for full implementation]

**Integration Plan**:
[How will this be integrated into biometal main codebase?]

**Risks to Monitor**:
- [Risk 1]
- [Risk 2]

### If NO-GO Decision

**Alternative Approaches**:
1. [Alternative 1] - [Why might this work better?]
2. [Alternative 2] - [Why might this work better?]

**Fallback Solution**:
[What should we do instead? Wrap existing tool?]

**What We Learned**:
- [Learning 1 that will help future experiments]
- [Learning 2 that will help future experiments]

---

## 7. Code Artifacts

**Location**: `experiments/[experiment-name]/`

**Key Files**:
- `src/lib.rs` - [Description]
- `src/neon.rs` - [Description]
- `benches/benchmark.rs` - [Description]

**Commits**: [Hash range or tag]

**Preservation**:
- [ ] Code archived in Git
- [ ] Benchmarks documented
- [ ] Can be reproduced from docs

---

## 8. Publication Plan

### If Successful (GO)

**Paper Draft**: [Link to draft if started]

**Target Venue**:
- Primary: [Journal/Conference name]
- Secondary: [Fallback option]

**Key Contributions**:
1. [Contribution 1]
2. [Contribution 2]
3. [Contribution 3]

**Timeline**:
- Draft: [Date]
- Submit: [Date]
- Expected publication: [Date]

### If Unsuccessful (NO-GO)

**Negative Results Report**: [Link]

**Where to Share**:
- [ ] Blog post
- [ ] GitHub discussion
- [ ] Conference poster (if relevant)
- [ ] Technical report

**Value of Negative Results**:
[What does the community learn from this failure?]

---

## 9. Lessons Learned

### What Worked Well

1. [Success 1]
2. [Success 2]
3. [Success 3]

### What Could Be Improved

1. [Improvement 1]
2. [Improvement 2]
3. [Improvement 3]

### Advice for Future Experiments

[Based on this experience, what would you tell someone starting a similar experiment?]

---

## 10. References

**Papers Cited**:
- [Citation 1]
- [Citation 2]

**Code References**:
- [Repository 1]
- [Repository 2]

**Datasets Used**:
- [Dataset 1] - [Where to find it]
- [Dataset 2] - [Where to find it]

---

## Appendices

### Appendix A: Detailed Benchmark Data

[Raw benchmark output, CSV files, statistical analysis]

### Appendix B: Code Listings

[Key code snippets if not obvious from repository]

### Appendix C: Design Documents

[Any architectural diagrams, format specifications studied, etc.]

---

**Document Version**: 1.0
**Last Updated**: YYYY-MM-DD
**Status**: [Draft / Final]

---

## Sign-Off

**Prepared by**: [Name]
**Reviewed by**: [Name(s)]
**Approved for**: [Next Action: Full Implementation / Archive / Pivot]
**Date**: YYYY-MM-DD
