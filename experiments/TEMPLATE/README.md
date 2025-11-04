# Experiment Template

This is a template for new biometal experiments. Copy this entire directory to create a new experiment.

## Quick Start

```bash
# From experiments/ directory
cp -r TEMPLATE my-experiment
cd my-experiment

# Rename package in Cargo.toml
sed -i '' 's/experiment-name/my-experiment/g' Cargo.toml

# Fill out proposal
vim PROPOSAL.md

# Initialize Git tracking
git add .
git commit -m "feat(experiment): Start my-experiment"

# Start research!
cargo build
cargo test
cargo bench
```

## Directory Structure

```
my-experiment/
├── Cargo.toml          # Rust package manifest
├── PROPOSAL.md         # Fill this out FIRST
├── RESEARCH_LOG.md     # Update DAILY
├── FINDINGS.md         # Fill out at END
├── README.md           # This file - update with specifics
├── src/
│   ├── lib.rs         # Core implementation
│   └── neon.rs        # ARM NEON optimizations
├── benches/
│   └── benchmark.rs   # Performance benchmarks (N=30)
├── tests/
│   └── correctness.rs # Validate vs baseline
└── docs/
    └── methodology.md # Detailed technical docs
```

## Experiment Checklist

### Before Starting (Week 0)

- [ ] Fill out `PROPOSAL.md` completely
- [ ] Define clear hypothesis
- [ ] Set quantitative success criteria
- [ ] Get approval (if needed)
- [ ] Update `../.experiments.toml` with entry

### Research Phase (Week 1-2)

- [ ] Study existing implementations
- [ ] Identify ARM NEON opportunities
- [ ] Implement scalar baseline
- [ ] Implement NEON optimization
- [ ] Write benchmarks (N=30)
- [ ] Update `RESEARCH_LOG.md` daily
- [ ] Run statistical analysis

### Decision Point (End of Week 2)

- [ ] Review benchmark results
- [ ] Check against success criteria
- [ ] Make GO / NO-GO / MAYBE decision
- [ ] Document in `FINDINGS.md`
- [ ] Update `../.experiments.toml` with outcome

### Implementation Phase (Week 3-6, if GO)

- [ ] Full implementation
- [ ] Integration with biometal
- [ ] Comprehensive testing
- [ ] Documentation
- [ ] Paper draft (if publishing)

### Completion

- [ ] Archive code in Git
- [ ] Complete `FINDINGS.md`
- [ ] Update `../.experiments.toml` with final status
- [ ] Publish results (paper or blog post)
- [ ] Share learnings with team

## Benchmarking Guidelines

All benchmarks must use:
- **N=30**: Minimum sample size
- **Statistical rigor**: Mean, 95% CI, Cohen's d
- **Fair comparison**: Same input, same conditions
- **Reproducible**: Document setup completely

Example:

```rust
use criterion::{criterion_group, criterion_main, Criterion};

fn benchmark_operation(c: &mut Criterion) {
    let mut group = c.benchmark_group("my_operation");

    // N=30 for statistical validity
    group.sample_size(30);

    let data = setup_test_data();

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

## Testing Guidelines

Ensure correctness before benchmarking:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn test_matches_baseline() {
        let input = /* ... */;
        assert_eq!(operation_scalar(&input), operation_neon(&input));
    }

    proptest! {
        #[test]
        fn test_property(input in /* generator */) {
            // Property that must hold for all inputs
            prop_assert_eq!(
                operation_scalar(&input),
                operation_neon(&input)
            );
        }
    }
}
```

## Documentation

Each experiment must produce:

1. **PROPOSAL.md** - Hypothesis and plan
2. **RESEARCH_LOG.md** - Daily progress
3. **FINDINGS.md** - Results and decision
4. **Code** - Well-commented, tested
5. **Benchmarks** - Reproducible results

## Success Criteria Example

**GO Decision** (proceed to full implementation):
- Speedup: ≥10× (NEON vs scalar)
- Memory: Constant (≤5 MB)
- Correctness: 100% match with baseline
- Implementation: Clean, maintainable

**NO-GO Decision** (wrap existing or abandon):
- Speedup: <5×
- Memory: Unbounded growth
- Correctness: Fails validation
- Implementation: Too complex

## Getting Help

- Review [experiments/README.md](../README.md) for process
- Check [OPTIMIZATION_RULES.md](../../OPTIMIZATION_RULES.md) for evidence
- See [sra-decoder/](../sra-decoder/) for example (once created)

## Tips

1. **Start small**: Prototype hot path first
2. **Benchmark early**: Validate hypothesis quickly
3. **Document daily**: Future you will thank current you
4. **Fail fast**: If it's not working, pivot or stop
5. **Share learnings**: Negative results are valuable

---

**Remember**: Experiments can fail. That's okay. We learn from failures and build on successes.
