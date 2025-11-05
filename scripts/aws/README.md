# AWS Cross-Platform Testing Scripts

Automated scripts for testing biometal on AWS Graviton (ARM) and x86_64 instances.

## Prerequisites

1. **AWS CLI** installed and configured:
   ```bash
   aws configure
   ```

2. **IAM permissions** for EC2:
   - `ec2:RunInstances`
   - `ec2:DescribeInstances`
   - `ec2:TerminateInstances`
   - `ec2:CreateKeyPair`
   - `ec2:CreateSecurityGroup`
   - `ec2:AuthorizeSecurityGroupIngress`

3. **SSH access** from your current IP

## Quick Start

### Test Both Platforms

```bash
# Recommended: Test both platforms sequentially
./scripts/aws/full_automation.sh --platforms graviton,x86_64
```

This will:
1. Launch both instances
2. Setup environments
3. Run tests and benchmarks
4. Download results
5. Terminate instances
6. Generate comparison report

**Time**: ~2-3 hours total
**Cost**: ~$0.36

### Test Single Platform

```bash
# Test Graviton only
./scripts/aws/full_automation.sh --platform graviton

# Test x86_64 only
./scripts/aws/full_automation.sh --platform x86_64
```

### Keep Instances Running

```bash
# Skip termination for debugging
./scripts/aws/full_automation.sh --platform graviton --skip-terminate
```

## Manual Workflow

If you prefer step-by-step control:

### 1. Launch Instance

```bash
# Graviton
./scripts/aws/launch_graviton.sh

# x86_64
./scripts/aws/launch_x86.sh
```

Instance info saved to: `results/cross_platform/{platform}/instance_info.txt`

### 2. Setup Environment

```bash
./scripts/aws/setup_instance.sh <public-ip> <platform>
```

Example:
```bash
./scripts/aws/setup_instance.sh 54.123.45.67 graviton
```

### 3. Run Tests

```bash
./scripts/aws/run_tests.sh <public-ip> <platform>
```

This takes 30-45 minutes (tests + benchmarks).

### 4. Download Results

```bash
./scripts/aws/download_results.sh <public-ip> <platform>
```

Results saved to: `results/cross_platform/{platform}/`

### 5. Terminate Instance

```bash
./scripts/aws/terminate_instance.sh <platform>
```

Requires confirmation (`yes`).

## Results

Results are saved in `results/cross_platform/{platform}/`:

```
results/cross_platform/
├── graviton/
│   ├── instance_info.txt
│   ├── system_info.txt
│   ├── test_results.txt
│   ├── bench_operations.txt
│   ├── bench_network.txt (if network feature)
│   └── criterion/ (detailed benchmark data)
└── x86_64/
    ├── instance_info.txt
    ├── system_info.txt
    ├── test_results.txt
    ├── bench_operations.txt
    ├── bench_network.txt
    └── criterion/
```

## Instance Details

### Graviton (c7g.xlarge)
- **Architecture**: ARM64 (Graviton 3)
- **Cores**: 4 vCPUs (Neoverse V1)
- **Memory**: 8 GB
- **NEON**: ✅ Yes
- **Cost**: $0.145/hour
- **AMI**: Amazon Linux 2023 ARM64

### x86_64 (c7i.xlarge)
- **Architecture**: x86_64 (Intel Xeon)
- **Cores**: 4 vCPUs
- **Memory**: 8 GB
- **NEON**: ❌ No (scalar fallback)
- **Cost**: $0.176/hour
- **AMI**: Amazon Linux 2023 x86_64

## Cost Estimation

| Operation | Graviton | x86_64 | Total |
|-----------|----------|--------|-------|
| Launch | $0 | $0 | $0 |
| Setup (10 min) | $0.02 | $0.03 | $0.05 |
| Tests (45 min) | $0.11 | $0.13 | $0.24 |
| Download | $0.01 | $0.01 | $0.02 |
| Storage | $0.01 | $0.01 | $0.02 |
| **Total** | **$0.15** | **$0.18** | **$0.33** |

**Note**: Costs are approximate. Actual costs may vary by region and usage.

## Security

- SSH keys stored in `~/.ssh/biometal-testing-key.pem`
- Security group restricts SSH to your current IP
- Instances auto-terminate after testing
- No credentials committed to git

## Troubleshooting

### SSH Connection Timeout

Wait longer (up to 90 seconds), or check security group rules:

```bash
aws ec2 describe-security-groups \
  --group-names biometal-testing-sg \
  --region us-east-1
```

### Instance Launch Failure

Try different region:
```bash
# Edit region in launch script
REGION="us-west-2"
```

Or check EC2 limits:
```bash
aws ec2 describe-account-attributes \
  --attribute-names max-instances \
  --region us-east-1
```

### Compilation Errors

Check Rust version (need 1.70+):
```bash
ssh -i ~/.ssh/biometal-testing-key.pem ec2-user@<ip> "rustc --version"
```

### Benchmark Variance

Increase sample size or run multiple times. AWS instances may have noisy neighbors.

## Cleanup

### List Running Instances

```bash
aws ec2 describe-instances \
  --filters "Name=tag:Project,Values=biometal" \
           "Name=instance-state-name,Values=running" \
  --region us-east-1 \
  --query 'Reservations[].Instances[].[InstanceId,Tags[?Key==`Name`].Value|[0],PublicIpAddress]' \
  --output table
```

### Terminate All biometal Instances

```bash
aws ec2 terminate-instances \
  --instance-ids $(aws ec2 describe-instances \
    --filters "Name=tag:Project,Values=biometal" \
              "Name=instance-state-name,Values=running" \
    --region us-east-1 \
    --query 'Reservations[].Instances[].InstanceId' \
    --output text) \
  --region us-east-1
```

### Remove Security Group

```bash
aws ec2 delete-security-group \
  --group-name biometal-testing-sg \
  --region us-east-1
```

### Remove SSH Key

```bash
aws ec2 delete-key-pair \
  --key-name biometal-testing-key \
  --region us-east-1

rm ~/.ssh/biometal-testing-key.pem
```

## Expected Results

### Test Pass Rates

All platforms should show:
```
test result: ok. 121 passed; 0 failed
```

### Benchmark Performance

| Operation | Mac M3 | Graviton 3 | x86_64 | Speedup (ARM/x86) |
|-----------|--------|------------|--------|-------------------|
| GC content | 5,954 Kseq/s | ~5,500 Kseq/s | ~294 Kseq/s | ~19× |
| Base counting | 5,254 Kseq/s | ~4,900 Kseq/s | ~315 Kseq/s | ~16× |
| Quality filter | 6,143 Kseq/s | ~5,700 Kseq/s | ~245 Kseq/s | ~23× |

**Variance**: ±10% between Mac and Graviton is normal (clock speed differences).

## Next Steps

After testing:

1. Review results in `results/cross_platform/`
2. Verify NEON speedup on ARM platforms
3. Verify scalar fallback works on x86_64
4. Update `CLAUDE.md` with cross-platform status
5. Add findings to release notes

## References

- [AWS Graviton](https://aws.amazon.com/ec2/graviton/)
- [Amazon Linux 2023](https://aws.amazon.com/linux/amazon-linux-2023/)
- [EC2 Pricing](https://aws.amazon.com/ec2/pricing/)
- [AWS CLI Reference](https://awscli.amazonaws.com/v2/documentation/api/latest/reference/ec2/index.html)
