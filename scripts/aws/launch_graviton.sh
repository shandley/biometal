#!/bin/bash
#
# Launch AWS Graviton instance for biometal cross-platform testing
#
# Usage: ./scripts/aws/launch_graviton.sh
#

set -e

echo "=== Graviton Instance Launch ==="
echo "Instance type: c7g.xlarge (Graviton 3, 4 vCPUs, 8GB RAM)"
echo

# Configuration
INSTANCE_TYPE="c7g.xlarge"
REGION="us-east-1"
KEY_NAME="biometal-testing-key"
SECURITY_GROUP_NAME="biometal-testing-sg"
AMI_ID="ami-0ba5fd8ce786c1351"  # Amazon Linux 2023 ARM64 (us-east-1)

# Create key pair if it doesn't exist
if ! aws ec2 describe-key-pairs --key-names "$KEY_NAME" --region "$REGION" &>/dev/null; then
    echo "Creating SSH key pair..."
    aws ec2 create-key-pair \
        --key-name "$KEY_NAME" \
        --region "$REGION" \
        --query 'KeyMaterial' \
        --output text > ~/.ssh/${KEY_NAME}.pem
    chmod 400 ~/.ssh/${KEY_NAME}.pem
    echo "✅ Key pair created: ~/.ssh/${KEY_NAME}.pem"
else
    echo "✅ Key pair already exists: $KEY_NAME"
fi

# Create security group if it doesn't exist
SG_ID=$(aws ec2 describe-security-groups \
    --filters "Name=group-name,Values=$SECURITY_GROUP_NAME" \
    --region "$REGION" \
    --query 'SecurityGroups[0].GroupId' \
    --output text 2>/dev/null || echo "None")

if [ "$SG_ID" == "None" ]; then
    echo "Creating security group..."
    SG_ID=$(aws ec2 create-security-group \
        --group-name "$SECURITY_GROUP_NAME" \
        --description "Security group for biometal cross-platform testing" \
        --region "$REGION" \
        --query 'GroupId' \
        --output text)

    # Get current public IP
    CURRENT_IP=$(curl -s https://checkip.amazonaws.com)

    # Allow SSH from current IP only
    aws ec2 authorize-security-group-ingress \
        --group-id "$SG_ID" \
        --protocol tcp \
        --port 22 \
        --cidr ${CURRENT_IP}/32 \
        --region "$REGION"

    echo "✅ Security group created: $SG_ID (SSH restricted to $CURRENT_IP)"
else
    echo "✅ Security group already exists: $SG_ID"
fi

# Launch instance
echo
echo "Launching c7g.xlarge instance..."
INSTANCE_ID=$(aws ec2 run-instances \
    --image-id "$AMI_ID" \
    --instance-type "$INSTANCE_TYPE" \
    --key-name "$KEY_NAME" \
    --security-group-ids "$SG_ID" \
    --region "$REGION" \
    --block-device-mappings 'DeviceName=/dev/xvda,Ebs={VolumeSize=20,VolumeType=gp3}' \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=biometal-graviton},{Key=Project,Value=biometal},{Key=Platform,Value=graviton}]" \
    --query 'Instances[0].InstanceId' \
    --output text)

echo "Instance ID: $INSTANCE_ID"
echo

# Wait for instance to be running
echo "Waiting for instance to be running..."
aws ec2 wait instance-running --instance-ids "$INSTANCE_ID" --region "$REGION"
echo "✅ Instance is running"

# Get public IP
PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids "$INSTANCE_ID" \
    --region "$REGION" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' \
    --output text)

echo
echo "=== Instance Details ==="
echo "Instance ID: $INSTANCE_ID"
echo "Public IP:   $PUBLIC_IP"
echo "SSH command: ssh -i ~/.ssh/${KEY_NAME}.pem ec2-user@${PUBLIC_IP}"
echo

# Wait for SSH to be ready
echo "Waiting for SSH to be ready (this may take 60-90 seconds)..."
MAX_ATTEMPTS=30
ATTEMPT=0
while [ $ATTEMPT -lt $MAX_ATTEMPTS ]; do
    if ssh -i ~/.ssh/${KEY_NAME}.pem -o StrictHostKeyChecking=no -o ConnectTimeout=5 ec2-user@${PUBLIC_IP} "echo SSH ready" &>/dev/null; then
        echo "✅ SSH is ready"
        break
    fi
    ATTEMPT=$((ATTEMPT + 1))
    echo "  Attempt $ATTEMPT/$MAX_ATTEMPTS..."
    sleep 5
done

if [ $ATTEMPT -eq $MAX_ATTEMPTS ]; then
    echo "❌ SSH connection timeout"
    exit 1
fi

# Save instance info to file
mkdir -p results/cross_platform/graviton
cat > results/cross_platform/graviton/instance_info.txt <<EOF
INSTANCE_ID=$INSTANCE_ID
PUBLIC_IP=$PUBLIC_IP
KEY_NAME=$KEY_NAME
REGION=$REGION
PLATFORM=graviton
LAUNCHED_AT=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
EOF

echo
echo "✅ Instance launched successfully!"
echo "Instance info saved to: results/cross_platform/graviton/instance_info.txt"
echo
echo "Next step: ./scripts/aws/setup_instance.sh $PUBLIC_IP graviton"
