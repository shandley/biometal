# Cross-Platform Testing Plan

**Objective**: Validate biometal performance and correctness across Mac ARM, Linux ARM (Graviton), and Linux x86_64.

## Testing Strategy

Based on the successful approach used in `apple-silicon-bio-bench` for Graviton validation, we will:

1. **Launch AWS instances** (Graviton ARM and x86_64)
2. **Setup environment** (Rust, dependencies)
3. **Run test suite** (Rust tests + Python bindings)
4. **Run benchmarks** (criterion benchmarks)
5. **Download results** (metrics and logs)
6. **Terminate instances** (cost control)
7. **Analyze results** (cross-platform comparison)

## Test Matrix

| Platform | Instance Type | vCPU | RAM | Cost/hr | NEON Support |
|----------|---------------|------|-----|---------|--------------|
| Mac M3 Max | Local | 14 | 36GB | $0 | ✅ Yes |
| AWS Graviton 3 | c7g.xlarge | 4 | 8GB | $0.145 | ✅ Yes |
| AWS x86_64 | c7i.xlarge | 4 | 8GB | $0.176 | ❌ No (scalar) |

**Total cost per run**: ~$0.30-0.50 (assuming 1-2 hours per platform)

## Test Coverage

### 1. Rust Tests (cargo test)
- Unit tests: All modules
- Integration tests: Network streaming
- Doc tests: API examples
- Property-based tests: proptest

**Expected**: All 121 tests passing on all platforms

### 2. Python Bindings (test_python_bindings.py)
- Streaming tests: FASTQ/FASTA
- Operations tests: NEON-accelerated functions
- K-mer tests: ML preprocessing
- Integration tests: End-to-end workflows

**Expected**: All 8 tests passing on all platforms

### 3. Benchmarks (cargo bench)
- Operations: base_counting, gc_content, quality_filter
- Network streaming: HTTP range requests, LRU cache
- Comparison: ARM NEON vs x86_64 scalar

**Expected**:
- Graviton NEON: 16-25× speedup vs x86_64
- Mac NEON: Similar to Graviton (within 10%)
- x86_64: Scalar baseline (no NEON)

## Implementation

### Phase 1: Script Infrastructure (30 minutes)

Create automation scripts in `scripts/aws/`:

1. `launch_graviton.sh` - Launch c7g.xlarge (Graviton 3)
2. `launch_x86.sh` - Launch c7i.xlarge (Intel Xeon)
3. `setup_instance.sh` - Install Rust, deps, compile
4. `run_tests.sh` - Execute tests and benchmarks
5. `download_results.sh` - Fetch results from instance
6. `terminate_instance.sh` - Clean up instance
7. `full_automation.sh` - Orchestrate entire workflow

### Phase 2: Test Execution (1-2 hours per platform)

**Graviton Testing**:
```bash
./scripts/aws/full_automation.sh --platform graviton
```

**x86_64 Testing**:
```bash
./scripts/aws/full_automation.sh --platform x86
```

Each run:
1. Launch instance (~2 min)
2. Setup environment (~5 min)
3. Compile project (~3 min)
4. Run tests (~1 min)
5. Run benchmarks (~30 min, N=30 for statistical rigor)
6. Download results (~1 min)
7. Terminate instance (<1 min)

**Total time**: ~45-60 minutes per platform

### Phase 3: Analysis (30 minutes)

Compare results across platforms:
- Test pass rates (expect 100% on all)
- Benchmark performance (NEON speedup validation)
- Memory usage (streaming validation)
- Python bindings correctness

Generate report: `results/cross_platform/FINDINGS.md`

## AWS Configuration

### Instance Selection Rationale

**c7g.xlarge (Graviton 3)**:
- Latest Graviton generation
- ARM Neoverse V1 cores (NEON support)
- Same architecture as Mac M-series
- Cost-effective ($0.145/hr)
- Proven in ASBB project

**c7i.xlarge (Intel Xeon)**:
- Comparable specs to Graviton
- x86_64 architecture (scalar fallback)
- Similar price point ($0.176/hr)
- Validates portable fallback code

### Region Selection

**Primary**: `us-east-1` (N. Virginia)
- Lowest cost
- Highest availability
- Used in ASBB project

**Fallback**: `us-west-2` (Oregon)
- Alternative if us-east-1 capacity issues

### AMI Selection

**Graviton**: Amazon Linux 2023 ARM64
- AMI ID: `ami-0ba5fd8ce786c1351` (us-east-1)
- Rust-friendly
- Modern toolchain

**x86_64**: Amazon Linux 2023 x86_64
- AMI ID: `ami-0453ec754f44f9a4a` (us-east-1)
- Same base OS as Graviton

### Security

- Key pair: `biometal-testing-key` (auto-generated)
- Security group: SSH only (port 22)
- Restrict to developer IP (production)
- Auto-terminate after testing

## Cost Estimation

### Per-Run Costs

| Component | Graviton | x86_64 | Total |
|-----------|----------|--------|-------|
| Compute (1 hr) | $0.145 | $0.176 | $0.32 |
| Storage (20GB, 1 day) | $0.01 | $0.01 | $0.02 |
| Data transfer | $0.01 | $0.01 | $0.02 |
| **Total** | **$0.17** | **$0.19** | **$0.36** |

### Development Phase (10 runs)

Testing during development: **~$3.60**

### CI/CD (monthly)

If running on every PR (assume 20 PRs/month):
- 20 runs × $0.36 = **$7.20/month**

**Recommendation**: Manual trigger for cross-platform tests (not every commit)

## Automation Workflow

### Manual Trigger (Recommended)

```bash
# Test both platforms sequentially
./scripts/aws/full_automation.sh --platforms graviton,x86

# Test single platform
./scripts/aws/full_automation.sh --platform graviton
```

### GitHub Actions (Future)

```yaml
name: Cross-Platform Tests
on:
  workflow_dispatch:  # Manual trigger only
  push:
    tags:
      - 'v*'  # Trigger on version tags

jobs:
  aws-testing:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Configure AWS credentials
        uses: aws-actions/configure-aws-credentials@v2
      - name: Run cross-platform tests
        run: ./scripts/aws/full_automation.sh --platforms graviton,x86
```

## Expected Results

### Test Suite

| Platform | Tests | Status | Notes |
|----------|-------|--------|-------|
| Mac M3 | 121 Rust + 8 Python | ✅ Pass | Baseline |
| Graviton 3 | 121 Rust + 8 Python | ✅ Pass | NEON parity |
| x86_64 | 121 Rust + 8 Python | ✅ Pass | Scalar fallback |

### Benchmarks

| Operation | Mac M3 | Graviton 3 | x86_64 | Speedup (ARM/x86) |
|-----------|--------|------------|--------|-------------------|
| GC content | ~5,954 Kseq/s | ~5,500 Kseq/s | ~294 Kseq/s | 19-20× |
| Base counting | ~5,254 Kseq/s | ~4,900 Kseq/s | ~315 Kseq/s | 15-17× |
| Quality filter | ~6,143 Kseq/s | ~5,700 Kseq/s | ~245 Kseq/s | 23-25× |

**Expected variance**: ±10% between Mac and Graviton (clock speed differences)

### Python Bindings

All 8 Python tests should pass on all platforms:
- Version test
- FASTQ streaming
- FASTA streaming
- GC content
- Base counting
- Mean quality
- K-mer extraction
- Integration workflow

## Success Criteria

1. **Correctness**: 100% test pass rate on all platforms
2. **Performance**: ARM NEON 16-25× faster than x86_64 scalar
3. **Portability**: x86_64 fallback works correctly (no crashes)
4. **Python bindings**: All tests pass on all platforms
5. **Memory**: Streaming maintains ~5 MB constant memory

## Timeline

| Phase | Duration | Cost |
|-------|----------|------|
| Script creation | 30 min | $0 |
| Graviton testing | 1 hr | $0.17 |
| x86_64 testing | 1 hr | $0.19 |
| Analysis | 30 min | $0 |
| **Total** | **3 hrs** | **$0.36** |

## Next Steps

1. Create `scripts/aws/` directory
2. Adapt ASBB scripts for biometal
3. Run initial Graviton test
4. Run x86_64 test
5. Compare results
6. Document findings
7. Update CLAUDE.md with cross-platform status

## References

- ASBB Graviton scripts: `/Users/scotthandley/Code/apple-silicon-bio-bench/scripts/`
- ASBB Graviton findings: Entry 021 in lab notebook
- AWS Graviton: https://aws.amazon.com/ec2/graviton/
- Criterion benchmarking: https://github.com/bheisler/criterion.rs

## Security Notes

- Never commit AWS credentials
- Use IAM roles with minimal permissions
- Restrict security groups to known IPs
- Auto-terminate instances after testing
- Use spot instances for cost savings (optional)

## Troubleshooting

### Issue: SSH connection timeout

**Solution**: Wait longer (up to 90 seconds), or check security group rules

### Issue: Instance launch failure

**Solution**: Try different region (us-west-2), or different instance type

### Issue: Compilation errors

**Solution**: Check Rust version (need 1.70+), verify dependencies installed

### Issue: Benchmark variance

**Solution**: Increase N=30 to N=50, or run multiple times and average

## Alternative: Docker Testing

If AWS access is limited, can use Docker for x86_64 testing:

```bash
# Build x86_64 image
docker build --platform linux/amd64 -t biometal-test .

# Run tests
docker run --rm biometal-test cargo test

# Run benchmarks
docker run --rm biometal-test cargo bench
```

**Limitation**: Cannot test Graviton without AWS access
