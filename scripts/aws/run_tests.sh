#!/bin/bash
#
# Run tests and benchmarks on AWS instance
#
# Usage: ./scripts/aws/run_tests.sh <public-ip> <platform>
#        platform: graviton | x86_64
#

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <public-ip> <platform>"
    echo "  platform: graviton | x86_64"
    exit 1
fi

PUBLIC_IP=$1
PLATFORM=$2
KEY_NAME="biometal-testing-key"
SSH_KEY=~/.ssh/${KEY_NAME}.pem

echo "=== Running Tests ($PLATFORM) ==="
echo "Public IP: $PUBLIC_IP"
echo

# Create test script to run on instance
cat > /tmp/biometal_run_tests.sh <<'EOF'
#!/bin/bash
set -e

source ~/.cargo/env
cd ~/biometal

echo "=== biometal Cross-Platform Testing ==="
echo "Platform: $(uname -m)"
echo "Date: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
echo

# 1. Run Rust tests
echo "=== Phase 1: Rust Tests ==="
echo "Running cargo test..."
cargo test --release 2>&1 | tee test_results.txt
echo "✅ Tests complete"
echo

# 2. Run benchmarks (reduced iterations for faster testing)
echo "=== Phase 2: Benchmarks ==="
echo "Running cargo bench (this will take 20-30 minutes)..."

# Run operations benchmarks
if [ -f benches/operations.rs ]; then
    echo "Running operations benchmarks..."
    cargo bench --bench operations 2>&1 | tee bench_operations.txt
fi

# Run network streaming benchmarks (if network feature enabled)
if cargo metadata --format-version 1 | grep -q '"network"'; then
    if [ -f benches/network_streaming.rs ]; then
        echo "Running network streaming benchmarks..."
        cargo bench --bench network_streaming 2>&1 | tee bench_network.txt
    fi
fi

echo "✅ Benchmarks complete"
echo

# 3. Collect system information
echo "=== Phase 3: System Information ==="
cat > system_info.txt <<SYSINFO
Platform: $(uname -m)
Kernel: $(uname -r)
OS: $(cat /etc/os-release | grep PRETTY_NAME | cut -d'"' -f2)
CPU: $(lscpu | grep "Model name" | cut -d':' -f2 | xargs)
Cores: $(nproc)
Memory: $(free -h | grep Mem | awk '{print $2}')
Rust: $(rustc --version)
Cargo: $(cargo --version)
Test Date: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
SYSINFO

cat system_info.txt
echo

echo "=== Testing Complete ==="
echo "Results saved to:"
echo "  - test_results.txt"
echo "  - bench_operations.txt"
echo "  - system_info.txt"
EOF

# Transfer and execute test script
echo "Transferring test script to instance..."
scp -i "$SSH_KEY" -o StrictHostKeyChecking=no /tmp/biometal_run_tests.sh ec2-user@${PUBLIC_IP}:~/run_tests.sh
echo

echo "Executing tests on instance..."
echo "⏱️  This will take 30-45 minutes (tests + benchmarks)..."
echo

# Run tests and show output in real-time
ssh -i "$SSH_KEY" -o StrictHostKeyChecking=no ec2-user@${PUBLIC_IP} "chmod +x ~/run_tests.sh && ~/run_tests.sh"

echo
echo "✅ Tests completed successfully"
echo
echo "Next step: ./scripts/aws/download_results.sh $PUBLIC_IP $PLATFORM"
