#!/bin/bash

# generate_test_bam.sh - Generate synthetic BAM files for profiling
#
# Usage: ./generate_test_bam.sh <num_records>
#
# Requires: samtools (install via: brew install samtools)

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <num_records>"
    echo
    echo "Example:"
    echo "  $0 10000     # Generate 10K records (~1 MB)"
    echo "  $0 100000    # Generate 100K records (~10 MB)"
    echo "  $0 1000000   # Generate 1M records (~100 MB)"
    exit 1
fi

NUM_RECORDS=$1
TEST_DATA_DIR="../test-data"
mkdir -p "$TEST_DATA_DIR"

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found"
    echo "Install via: brew install samtools"
    exit 1
fi

OUTPUT_SAM="$TEST_DATA_DIR/synthetic_${NUM_RECORDS}.sam"
OUTPUT_BAM="$TEST_DATA_DIR/synthetic_${NUM_RECORDS}.bam"

echo "Generating synthetic BAM file"
echo "============================="
echo "Records: $NUM_RECORDS"
echo "Output: $OUTPUT_BAM"
echo

# Generate SAM header
cat > "$OUTPUT_SAM" << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@SQ	SN:chr2	LN:242193529
@SQ	SN:chr22	LN:50818468
@PG	ID:generate_test_bam	PN:generate_test_bam	VN:1.0	CL:generate_test_bam.sh
EOF

echo "Generating $NUM_RECORDS synthetic records..."

# Generate synthetic reads
# Format: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
for i in $(seq 1 $NUM_RECORDS); do
    # Random position on chr1 (0-248M)
    POS=$((RANDOM % 248956422))

    # Random MAPQ (0-60)
    MAPQ=$((RANDOM % 61))

    # Simple 100bp read
    SEQ="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    QUAL="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"

    # CIGAR: 100M (100 matches)
    CIGAR="100M"

    # Read name
    QNAME="read_${i}"

    # FLAG: 0 (unmapped pair) or 163 (properly paired, first in pair)
    if [ $((i % 2)) -eq 0 ]; then
        FLAG=163
    else
        FLAG=83
    fi

    echo -e "${QNAME}\t${FLAG}\tchr1\t${POS}\t${MAPQ}\t${CIGAR}\t*\t0\t0\t${SEQ}\t${QUAL}" >> "$OUTPUT_SAM"

    # Progress indicator
    if [ $((i % 10000)) -eq 0 ]; then
        echo "  Generated $i / $NUM_RECORDS records..."
    fi
done

echo "Converting SAM to BAM..."
samtools view -b -o "$OUTPUT_BAM" "$OUTPUT_SAM"

echo "Sorting BAM..."
samtools sort -o "${OUTPUT_BAM}.sorted" "$OUTPUT_BAM"
mv "${OUTPUT_BAM}.sorted" "$OUTPUT_BAM"

echo "Indexing BAM..."
samtools index "$OUTPUT_BAM"

# Clean up SAM
rm "$OUTPUT_SAM"

# Report file sizes
BAM_SIZE=$(du -h "$OUTPUT_BAM" | cut -f1)
echo
echo "Complete!"
echo "  BAM file: $OUTPUT_BAM ($BAM_SIZE)"
echo "  Index: ${OUTPUT_BAM}.bai"
echo

# Validate
echo "Validating BAM..."
samtools quickcheck "$OUTPUT_BAM" && echo "  ✓ BAM file is valid"

echo
echo "Stats:"
samtools flagstat "$OUTPUT_BAM"

echo
echo "Ready to profile!"
echo "Run: cd profiling && cargo build --release && ./target/release/bam-profiling $OUTPUT_BAM"
