#!/usr/bin/env python3
"""Test script for BAM statistics functions."""

import biometal
import sys

def test_statistics():
    """Test BAM statistics functions."""
    bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam"
    print(f"Opening {bam_path}...")

    # Test 1: count_by_flag
    print("\n1. Testing count_by_flag()...")
    stats = biometal.count_by_flag(bam_path)

    print(f"\n   Flag Statistics:")
    print(f"   Total records: {stats['total']:,}")
    print(f"   Mapped: {stats['mapped']:,} ({100*stats['mapped']/stats['total']:.1f}%)")
    print(f"   Unmapped: {stats['unmapped']:,} ({100*stats['unmapped']/stats['total']:.1f}%)")
    print(f"   Primary: {stats['primary']:,} ({100*stats['primary']/stats['total']:.1f}%)")
    print(f"   Secondary: {stats['secondary']:,}")
    print(f"   Supplementary: {stats['supplementary']:,}")
    print(f"   Paired: {stats['paired']:,} ({100*stats['paired']/stats['total']:.1f}%)")
    print(f"   Proper pair: {stats['proper_pair']:,} ({100*stats['proper_pair']/stats['total']:.1f}%)")
    print(f"   Forward: {stats['forward']:,} ({100*stats['forward']/stats['total']:.1f}%)")
    print(f"   Reverse: {stats['reverse']:,} ({100*stats['reverse']/stats['total']:.1f}%)")
    print(f"   QC fail: {stats['qc_fail']:,}")
    print(f"   Duplicate: {stats['duplicate']:,}")

    # Test 2: mapq_distribution
    print("\n\n2. Testing mapq_distribution()...")
    mapq_dist = biometal.mapq_distribution(bam_path)

    print(f"\n   MAPQ Distribution (showing MAPQ scores with >1000 reads):")
    for mapq in sorted(mapq_dist.keys()):
        count = mapq_dist[mapq]
        if count > 1000:
            print(f"   MAPQ {mapq:3d}: {count:7,} reads ({100*count/stats['total']:5.2f}%)")

    # Test 3: mapq_distribution with reference filter
    print("\n\n3. Testing mapq_distribution() with reference_id filter...")
    mapq_dist_chr1 = biometal.mapq_distribution(bam_path, reference_id=0)

    total_chr1 = sum(mapq_dist_chr1.values())
    print(f"\n   MAPQ Distribution for chr1 (reference_id=0):")
    print(f"   Total chr1 reads: {total_chr1:,}")

    # Show top 5 MAPQ scores
    top_mapq = sorted(mapq_dist_chr1.items(), key=lambda x: x[1], reverse=True)[:5]
    for mapq, count in top_mapq:
        print(f"   MAPQ {mapq:3d}: {count:7,} reads ({100*count/total_chr1:5.2f}%)")

    # Test 4: calculate_coverage
    print("\n\n4. Testing calculate_coverage()...")
    print(f"   Calculating coverage for chr1:0-1000...")

    coverage = biometal.calculate_coverage(bam_path, reference_id=0, start=0, end=1000)

    total_positions = len(coverage)
    positions_covered = sum(1 for depth in coverage.values() if depth > 0)
    max_depth = max(coverage.values()) if coverage else 0
    avg_depth = sum(coverage.values()) / len(coverage) if coverage else 0

    print(f"\n   Coverage Statistics (chr1:0-1000):")
    print(f"   Positions with coverage: {positions_covered:,} / 1,000")
    print(f"   Maximum depth: {max_depth}")
    print(f"   Average depth: {avg_depth:.2f}")

    # Show first few positions
    print(f"\n   Sample positions (first 5 with coverage):")
    count = 0
    for pos in sorted(coverage.keys())[:1000]:
        depth = coverage[pos]
        if depth > 0:
            print(f"   Position {pos}: {depth}× coverage")
            count += 1
            if count >= 5:
                break

    # Test 5: calculate_coverage for different regions
    print("\n\n5. Testing calculate_coverage() for multiple regions...")

    regions = [
        (0, 500),
        (500, 1000),
        (1000, 1500),
    ]

    for start, end in regions:
        cov = biometal.calculate_coverage(bam_path, reference_id=0, start=start, end=end)
        positions = len(cov)
        avg = sum(cov.values()) / len(cov) if cov else 0
        max_d = max(cov.values()) if cov else 0
        print(f"   Region {start:4d}-{end:4d}: {positions:4d} positions, avg depth {avg:5.2f}, max depth {max_d}")

    print(f"\n\n✅ All statistics tests passed!")

if __name__ == "__main__":
    try:
        test_statistics()
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
