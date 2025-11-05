#!/bin/bash
#
# Full automation: Launch, test, and analyze cross-platform results
#
# Usage: ./scripts/aws/full_automation.sh --platform <platform>
#        ./scripts/aws/full_automation.sh --platforms graviton,x86_64
#
# Options:
#   --platform    Single platform to test (graviton | x86_64)
#   --platforms   Multiple platforms, comma-separated
#   --skip-terminate  Don't terminate instances after testing
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

# Parse arguments
PLATFORMS=""
SKIP_TERMINATE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --platform)
            PLATFORMS="$2"
            shift 2
            ;;
        --platforms)
            PLATFORMS="$2"
            shift 2
            ;;
        --skip-terminate)
            SKIP_TERMINATE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 --platform <platform> [--skip-terminate]"
            echo "       $0 --platforms graviton,x86_64 [--skip-terminate]"
            exit 1
            ;;
    esac
done

if [ -z "$PLATFORMS" ]; then
    echo "Error: No platform specified"
    echo "Usage: $0 --platform <platform> [--skip-terminate]"
    echo "       $0 --platforms graviton,x86_64 [--skip-terminate]"
    exit 1
fi

echo "========================================"
echo "  biometal Cross-Platform Testing"
echo "  Platforms: $PLATFORMS"
echo "========================================"
echo

# Convert comma-separated list to array
IFS=',' read -ra PLATFORM_ARRAY <<< "$PLATFORMS"

# Process each platform
for PLATFORM in "${PLATFORM_ARRAY[@]}"; do
    PLATFORM=$(echo "$PLATFORM" | xargs)  # Trim whitespace

    echo
    echo "========================================"
    echo "  Testing Platform: $PLATFORM"
    echo "========================================"
    echo

    # Validate platform
    if [ "$PLATFORM" != "graviton" ] && [ "$PLATFORM" != "x86_64" ]; then
        echo "❌ Invalid platform: $PLATFORM"
        echo "Valid platforms: graviton, x86_64"
        exit 1
    fi

    # Phase 1: Launch instance
    echo "=== Phase 1: Launch Instance ==="
    if [ "$PLATFORM" == "graviton" ]; then
        ./scripts/aws/launch_graviton.sh
    else
        ./scripts/aws/launch_x86.sh
    fi

    # Load instance info
    source "results/cross_platform/$PLATFORM/instance_info.txt"
    echo

    # Phase 2: Setup instance
    echo "=== Phase 2: Setup Instance ==="
    ./scripts/aws/setup_instance.sh "$PUBLIC_IP" "$PLATFORM"
    echo

    # Phase 3: Run tests
    echo "=== Phase 3: Run Tests (30-45 minutes) ==="
    ./scripts/aws/run_tests.sh "$PUBLIC_IP" "$PLATFORM"
    echo

    # Phase 4: Download results
    echo "=== Phase 4: Download Results ==="
    ./scripts/aws/download_results.sh "$PUBLIC_IP" "$PLATFORM"
    echo

    # Phase 5: Terminate instance (optional)
    if [ "$SKIP_TERMINATE" = false ]; then
        echo "=== Phase 5: Terminate Instance ==="
        echo "yes" | ./scripts/aws/terminate_instance.sh "$PLATFORM"
        echo
    else
        echo "=== Skipping Termination ==="
        echo "Instance still running: $PUBLIC_IP"
        echo "To terminate manually: ./scripts/aws/terminate_instance.sh $PLATFORM"
        echo
    fi

    echo "✅ $PLATFORM testing complete"
    echo
done

# Generate comparison report if multiple platforms tested
if [ ${#PLATFORM_ARRAY[@]} -gt 1 ]; then
    echo "========================================"
    echo "  Cross-Platform Comparison"
    echo "========================================"
    echo

    echo "Test results:"
    for PLATFORM in "${PLATFORM_ARRAY[@]}"; do
        PLATFORM=$(echo "$PLATFORM" | xargs)
        echo
        echo "=== $PLATFORM ==="
        if [ -f "results/cross_platform/$PLATFORM/system_info.txt" ]; then
            cat "results/cross_platform/$PLATFORM/system_info.txt"
        fi
        echo
        if [ -f "results/cross_platform/$PLATFORM/test_results.txt" ]; then
            echo "Test summary:"
            grep -i "test result" "results/cross_platform/$PLATFORM/test_results.txt" | tail -1
        fi
    done

    echo
    echo "Detailed results in: results/cross_platform/"
    echo
fi

# Done!
echo "========================================"
echo "  ✅ CROSS-PLATFORM TESTING COMPLETE!"
echo "========================================"
echo
echo "Results saved to: results/cross_platform/"
echo
echo "Next steps:"
echo "  1. Review test results in results/cross_platform/"
echo "  2. Compare benchmark performance across platforms"
echo "  3. Verify NEON speedup on ARM platforms"
echo "  4. Verify scalar fallback works correctly on x86_64"
echo "  5. Update CLAUDE.md with cross-platform status"
echo
