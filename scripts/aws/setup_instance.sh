#!/bin/bash
#
# Setup AWS instance: Install dependencies and compile biometal
#
# Usage: ./scripts/aws/setup_instance.sh <public-ip> <platform>
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

echo "=== Instance Setup ($PLATFORM) ==="
echo "Public IP: $PUBLIC_IP"
echo

# Create setup script to run on instance
cat > /tmp/biometal_setup_remote.sh <<'EOF'
#!/bin/bash
set -e

echo "=== Remote Setup Started ==="
echo

# Update system
echo "1. Updating system packages..."
sudo yum update -y
echo "✅ System updated"
echo

# Install development tools
echo "2. Installing development tools..."
sudo yum install -y gcc git cmake
echo "✅ Development tools installed"
echo

# Install Rust
echo "3. Installing Rust toolchain..."
if [ ! -f ~/.cargo/bin/rustc ]; then
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source ~/.cargo/env
    echo "✅ Rust installed"
else
    source ~/.cargo/env
    echo "✅ Rust already installed"
fi
echo

# Verify Rust installation
echo "4. Verifying Rust installation..."
rustc --version
cargo --version
echo "✅ Rust verified"
echo

# Detect architecture
echo "5. Detecting architecture..."
ARCH=$(uname -m)
echo "Architecture: $ARCH"
if [ "$ARCH" == "aarch64" ]; then
    echo "✅ ARM64 detected (NEON support)"
elif [ "$ARCH" == "x86_64" ]; then
    echo "✅ x86_64 detected (scalar fallback)"
else
    echo "⚠️  Unknown architecture: $ARCH"
fi
echo

# Create biometal directory for code transfer
echo "6. Creating biometal directory..."
mkdir -p ~/biometal
echo "✅ Directory created"

echo "=== Remote Setup Complete ==="
EOF

# Transfer and execute setup script
echo "Transferring setup script to instance..."
scp -i "$SSH_KEY" -o StrictHostKeyChecking=no /tmp/biometal_setup_remote.sh ec2-user@${PUBLIC_IP}:~/setup.sh
echo

echo "Executing setup script on instance (this may take 5-10 minutes)..."
ssh -i "$SSH_KEY" -o StrictHostKeyChecking=no ec2-user@${PUBLIC_IP} "chmod +x ~/setup.sh && ~/setup.sh"
echo

# Transfer biometal code
echo "Transferring biometal code to instance..."

# Create tarball of local code
echo "Creating code tarball..."
tar -czf /tmp/biometal-code.tar.gz \
    --exclude=target \
    --exclude=.git \
    --exclude=results \
    --exclude=.venv \
    --exclude='*.csv' \
    --exclude='*.txt' \
    -C "$(dirname "$(dirname "$(pwd)")")" \
    biometal

# Transfer tarball
echo "Uploading tarball to instance..."
scp -i "$SSH_KEY" -o StrictHostKeyChecking=no /tmp/biometal-code.tar.gz ec2-user@${PUBLIC_IP}:~/

# Extract on instance
echo "Extracting code on instance..."
ssh -i "$SSH_KEY" -o StrictHostKeyChecking=no ec2-user@${PUBLIC_IP} "tar -xzf biometal-code.tar.gz && rm biometal-code.tar.gz"

# Cleanup local tarball
rm /tmp/biometal-code.tar.gz
echo "✅ Code transferred successfully"
echo

# Compile project (default features, no Python for now)
echo "Compiling biometal (Rust only)..."
ssh -i "$SSH_KEY" -o StrictHostKeyChecking=no ec2-user@${PUBLIC_IP} "source ~/.cargo/env && cd ~/biometal && cargo build --release"
echo

# Verify compilation
echo "Verifying compilation..."
BUILD_SUCCESS=$(ssh -i "$SSH_KEY" -o StrictHostKeyChecking=no ec2-user@${PUBLIC_IP} "[ -d ~/biometal/target/release ] && echo 'yes' || echo 'no'")

if [ "$BUILD_SUCCESS" == "yes" ]; then
    echo "✅ Build successful"
else
    echo "❌ Build failed"
    exit 1
fi

echo
echo "=== Setup Complete ==="
echo "Next step: ./scripts/aws/run_tests.sh $PUBLIC_IP $PLATFORM"
