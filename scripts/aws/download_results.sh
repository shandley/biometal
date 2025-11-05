#!/bin/bash
#
# Download test results from AWS instance
#
# Usage: ./scripts/aws/download_results.sh <public-ip> <platform>
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

echo "=== Downloading Results ($PLATFORM) ==="
echo "Public IP: $PUBLIC_IP"
echo

# Create results directory
RESULTS_DIR="results/cross_platform/$PLATFORM"
mkdir -p "$RESULTS_DIR"

# Download test results
echo "Downloading test results..."
scp -i "$SSH_KEY" -o StrictHostKeyChecking=no \
    ec2-user@${PUBLIC_IP}:~/biometal/test_results.txt \
    "$RESULTS_DIR/" 2>/dev/null || echo "⚠️  test_results.txt not found"

# Download benchmark results
echo "Downloading benchmark results..."
scp -i "$SSH_KEY" -o StrictHostKeyChecking=no \
    ec2-user@${PUBLIC_IP}:~/biometal/bench_operations.txt \
    "$RESULTS_DIR/" 2>/dev/null || echo "⚠️  bench_operations.txt not found"

scp -i "$SSH_KEY" -o StrictHostKeyChecking=no \
    ec2-user@${PUBLIC_IP}:~/biometal/bench_network.txt \
    "$RESULTS_DIR/" 2>/dev/null || echo "⚠️  bench_network.txt not found"

# Download system info
echo "Downloading system information..."
scp -i "$SSH_KEY" -o StrictHostKeyChecking=no \
    ec2-user@${PUBLIC_IP}:~/biometal/system_info.txt \
    "$RESULTS_DIR/"

# Download criterion benchmark data (if exists)
echo "Downloading criterion data..."
scp -i "$SSH_KEY" -o StrictHostKeyChecking=no -r \
    ec2-user@${PUBLIC_IP}:~/biometal/target/criterion \
    "$RESULTS_DIR/" 2>/dev/null || echo "⚠️  criterion data not found"

echo
echo "✅ Results downloaded successfully"
echo "Results saved to: $RESULTS_DIR/"
echo
ls -lh "$RESULTS_DIR/"
echo
echo "Next step: ./scripts/aws/terminate_instance.sh $PLATFORM"
