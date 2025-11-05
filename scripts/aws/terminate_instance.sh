#!/bin/bash
#
# Terminate AWS instance
#
# Usage: ./scripts/aws/terminate_instance.sh <platform>
#        platform: graviton | x86_64
#

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <platform>"
    echo "  platform: graviton | x86_64"
    exit 1
fi

PLATFORM=$1
REGION="us-east-1"

echo "=== Terminating Instance ($PLATFORM) ==="
echo

# Load instance info
INSTANCE_INFO="results/cross_platform/$PLATFORM/instance_info.txt"

if [ ! -f "$INSTANCE_INFO" ]; then
    echo "❌ Instance info not found: $INSTANCE_INFO"
    echo "Cannot determine instance ID to terminate"
    exit 1
fi

source "$INSTANCE_INFO"

echo "Instance ID: $INSTANCE_ID"
echo "Region: $REGION"
echo

# Confirm termination
read -p "Are you sure you want to terminate this instance? (yes/no): " CONFIRM

if [ "$CONFIRM" != "yes" ]; then
    echo "Termination cancelled"
    exit 0
fi

# Terminate instance
echo "Terminating instance..."
aws ec2 terminate-instances \
    --instance-ids "$INSTANCE_ID" \
    --region "$REGION" \
    --output text

echo "✅ Termination initiated"
echo

# Wait for termination
echo "Waiting for instance to terminate..."
aws ec2 wait instance-terminated \
    --instance-ids "$INSTANCE_ID" \
    --region "$REGION"

echo "✅ Instance terminated"

# Archive instance info
TERMINATED_AT=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
echo "TERMINATED_AT=$TERMINATED_AT" >> "$INSTANCE_INFO"

echo
echo "=== Termination Complete ==="
echo "Instance info archived in: $INSTANCE_INFO"
