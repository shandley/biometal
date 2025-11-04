# biometal Experiments

This directory contains **research experiments** exploring novel approaches to bioinformatics computing on ARM hardware. Each experiment follows the scientific method with hypothesis, validation, and clear go/no-go criteria.

## 🎯 Purpose

**Experiments are NOT production code**. They are:
- Research prototypes
- Performance explorations
- Architectural investigations
- Validation of optimization hypotheses

Experiments may **succeed** (promoted to `src/`) or **fail** (documented as negative results).

## 📊 Experiment Lifecycle

### 1. **Proposal** → Create experiment directory

```bash
# Copy template
cp -r experiments/TEMPLATE experiments/my-experiment

# Edit proposal
vim experiments/my-experiment/PROPOSAL.md
```

### 2. **Research** → Time-boxed investigation (typically 1-2 weeks)

- Study existing approaches
- Identify ARM NEON opportunities
- Design prototype
- Set success criteria

### 3. **Prototype** → Build minimal viable implementation

- Focus on critical path (hot loops)
- Implement both NEON and scalar
- Write benchmarks (N=30, 95% CI, Cohen's d)

### 4. **Validate** → Benchmark against hypothesis

**Success criteria example:**
- ✅ GO: ≥10× speedup vs baseline
- 🤔 MAYBE: 5-10× speedup → deeper investigation
- ❌ NO-GO: <5× speedup → document and archive

### 5. **Decision** → Go/No-Go

**GO Decision:**
- Move to full implementation
- Integrate with main codebase
- Add to production roadmap
- Publish findings

**NO-GO Decision:**
- Document why it didn't work
- Archive experiment
- Share negative results (equally valuable!)
- Move to next experiment

### 6. **Publication** → Share findings (success or failure)

Each experiment should produce:
- Technical report (in experiment's `/docs`)
- Blog post (for negative results)
- Paper (for significant positive results)
- Lab notebook entry (in main project)

## 📁 Experiment Structure

Each experiment is a **self-contained Rust project**:

```
experiments/my-experiment/
├── Cargo.toml              # Independent package
├── PROPOSAL.md             # Initial hypothesis
├── RESEARCH_LOG.md         # Day-by-day progress
├── FINDINGS.md             # Final results
├── src/
│   ├── lib.rs
│   └── neon.rs
├── benches/
│   └── benchmark.rs        # N=30, statistical analysis
├── tests/
│   └── correctness.rs      # Validate vs baseline
└── docs/
    └── methodology.md      # Reproducible research
```

## 🧬 Active Experiments

See [.experiments.toml](.experiments.toml) for current status.

### sra-decoder (Week 5-6, 2025)

**Hypothesis**: ARM NEON-optimized SRA decoder can achieve 20-40× speedup vs SRA toolkit.

**Status**: Research phase
**Owner**: Scott Handley
**Timeline**: 2 weeks research, 4-6 weeks implementation if validated

## 📚 Experiment Guidelines

### 1. Evidence-Based Approach

Every experiment must:
- State clear hypothesis
- Define success criteria (quantitative)
- Use statistical rigor (N=30, 95% CI, Cohen's d)
- Document methodology (reproducible)

### 2. Time-Boxing

- **Research phase**: 1-2 weeks MAX
- **Prototype phase**: 2-3 weeks MAX
- **Decision point**: After prototype benchmarks
- Total: 3-5 weeks from start to go/no-go

### 3. Independent Development

- Each experiment is a separate Cargo package
- No dependencies on biometal internals (initially)
- Can break things without affecting production
- Fast iteration without stability constraints

### 4. Comparative Benchmarking

Always benchmark against:
- Baseline (scalar implementation)
- State-of-the-art (existing tools)
- Theoretical maximum (hardware limits)

**Example benchmark:**
```rust
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

fn benchmark_operation(c: &mut Criterion) {
    let mut group = c.benchmark_group("sra_decode");

    // N=30 for statistical validity
    group.sample_size(30);

    group.bench_function("scalar", |b| {
        b.iter(|| operation_scalar(&data))
    });

    group.bench_function("neon", |b| {
        b.iter(|| operation_neon(&data))
    });

    group.finish();
}

criterion_group!(benches, benchmark_operation);
criterion_main!(benches);
```

### 5. Documentation Requirements

Each experiment must produce:

**PROPOSAL.md** - Before starting:
- Problem statement
- Hypothesis
- Success criteria
- Timeline
- Resources needed

**RESEARCH_LOG.md** - During experiment:
- Daily progress notes
- Challenges encountered
- Design decisions
- Benchmark results

**FINDINGS.md** - After completion:
- Summary of results
- Speedup achieved (with statistics)
- Go/No-Go recommendation
- Future work (if applicable)

## 🎓 Publication Strategy

### Successful Experiments → Academic Papers

- **Bioinformatics**: Methods and algorithms
- **BMC Bioinformatics**: Software and tools
- **GigaScience**: Large-scale data analysis
- **Genome Biology**: High-impact discoveries

### Failed Experiments → Negative Results

- Blog posts (transparency)
- GitHub discussions
- Conference presentations (poster sessions)
- "What we learned" documentation

**Why publish negative results?**
- Saves others from repeating failed approaches
- Shows rigorous scientific process
- Builds trust in positive results
- Contributes to knowledge

## 🔬 Example: SRA Decoder Experiment

This is our first experiment. Let's walk through the process:

### Week 1-2: Research Phase

**Days 1-3**: Format Analysis
- Study NCBI SRA specification
- Identify hot paths (base unpacking, quality decompression)
- Map to ARM NEON instructions

**Days 4-8**: Prototype Critical Path
- Implement base unpacking (2-bit → ASCII) with NEON
- Implement quality decompression with NEON
- Write benchmarks (N=30)

**Days 9-10**: Benchmark & Analysis
- Run statistical analysis (mean, 95% CI, Cohen's d)
- Compare vs SRA toolkit
- Document findings

**Days 11-12**: Decision
- **GO**: ≥10× speedup → proceed to full implementation
- **NO-GO**: <5× speedup → wrap toolkit, document why

### Week 3-6: Implementation Phase (if GO)

- Full SRA decoder with streaming architecture
- Integration with biometal
- Comprehensive testing
- Documentation
- Benchmarking paper draft

## 🚀 Future Experiment Ideas

Potential experiments to explore:

### 1. **BAM/CRAM Decoders**
- Hypothesis: ARM NEON can accelerate binary alignment format parsing
- Timeline: 4-6 weeks
- Impact: Enables streaming variant calling

### 2. **GPU Acceleration (Apple Metal)**
- Hypothesis: Apple GPU can accelerate massive parallel operations
- Timeline: 6-8 weeks
- Impact: Leverage full M-series chip capabilities

### 3. **Custom Compression**
- Hypothesis: Genomics-specific compression beats general-purpose (gzip)
- Timeline: 4-6 weeks
- Impact: Reduce network bandwidth, storage costs

### 4. **Async Network Layer**
- Hypothesis: Tokio async improves network streaming performance
- Timeline: 3-4 weeks
- Impact: Better concurrency, lower latency

### 5. **FASTQ Format Extensions**
- Hypothesis: Enhanced FASTQ format enables better streaming
- Timeline: 2-3 weeks
- Impact: Propose new standard to community

## 📝 Creating a New Experiment

```bash
# 1. Copy template
cd experiments
cp -r TEMPLATE my-experiment
cd my-experiment

# 2. Edit proposal
vim PROPOSAL.md

# 3. Initialize Cargo project
cargo init --lib

# 4. Add to registry
cd ..
echo "[[experiment]]" >> .experiments.toml
echo "name = \"my-experiment\"" >> .experiments.toml
echo "status = \"research\"" >> .experiments.toml
echo "start_date = \"2025-11-05\"" >> .experiments.toml

# 5. Start research!
cd my-experiment
cargo build
cargo test
cargo bench
```

## 🎯 Success Metrics

An experiment is considered successful if it:

1. **Technical Merit**:
   - Achieves stated performance goals
   - Maintains or improves memory efficiency
   - Compatible with biometal architecture

2. **Scientific Rigor**:
   - Reproducible methodology
   - Statistical validation (N=30, 95% CI)
   - Clear documentation

3. **Practical Value**:
   - Solves real user problem
   - Integrates with existing workflows
   - Maintainable long-term

## 🤝 Contributing Experiments

External contributors welcome! To propose an experiment:

1. Open GitHub issue with "Experiment:" prefix
2. Describe hypothesis and success criteria
3. Estimate timeline and resources
4. Wait for approval before starting
5. Follow process documented here

## 📖 References

- **Evidence Base**: [OPTIMIZATION_RULES.md](../OPTIMIZATION_RULES.md)
- **Lab Notebook**: [apple-silicon-bio-bench](https://github.com/shandley/apple-silicon-bio-bench)
- **Architecture**: [docs/ARCHITECTURE.md](../docs/ARCHITECTURE.md)

---

**Remember**: Experiments can fail. That's the point. We learn from failures and build on successes.

**Last Updated**: November 4, 2025
**Active Experiments**: 1 (sra-decoder)
