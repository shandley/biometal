#!/bin/bash

# get_test_data.sh - Download public BAM test files
#
# This script downloads small-to-medium BAM files from public sources
# for profiling noodles performance on ARM.

set -euo pipefail

TEST_DATA_DIR="../test-data"
mkdir -p "$TEST_DATA_DIR"

echo "BAM Test Data Downloader"
echo "========================"
echo

# Option 1: Small BAM from 1000 Genomes (chromosome 22, ~100-200MB)
echo "Option 1: 1000 Genomes Project - HG00096 chromosome 22"
echo "  Source: 1000 Genomes Phase 3"
echo "  Technology: Illumina (low coverage)"
echo "  Size: ~100-200 MB compressed"
echo "  Records: ~500K-1M"
echo

read -p "Download this file? [y/N] " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Downloading HG00096.chrom22.ILLUMINA.bwa.GBR.low_coverage.20120522.bam..."

    URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom22.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"

    if command -v wget &> /dev/null; then
        wget -O "$TEST_DATA_DIR/HG00096.chrom22.bam" "$URL"
    elif command -v curl &> /dev/null; then
        curl -o "$TEST_DATA_DIR/HG00096.chrom22.bam" "$URL"
    else
        echo "Error: Neither wget nor curl found. Install one and try again."
        exit 1
    fi

    echo "Download complete: $TEST_DATA_DIR/HG00096.chrom22.bam"

    # Download index too
    echo "Downloading BAM index..."
    INDEX_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom22.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai"

    if command -v wget &> /dev/null; then
        wget -O "$TEST_DATA_DIR/HG00096.chrom22.bam.bai" "$INDEX_URL"
    elif command -v curl &> /dev/null; then
        curl -o "$TEST_DATA_DIR/HG00096.chrom22.bam.bai" "$INDEX_URL"
    fi

    echo "Index downloaded: $TEST_DATA_DIR/HG00096.chrom22.bam.bai"
fi

echo
echo "Alternative: Generate synthetic BAM files"
echo "  Run: ./generate_test_bam.sh <num_records>"
echo "  Example: ./generate_test_bam.sh 10000"
echo

echo "Test data directory: $TEST_DATA_DIR"
ls -lh "$TEST_DATA_DIR" || true
