#!/usr/bin/env python3
"""Test script for BAM filtering helper properties."""

import biometal
import sys

def test_filtering():
    """Test filtering helper properties."""
    bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam"
    print(f"Opening {bam_path}...")
    bam = biometal.BamReader.from_path(bam_path)

    # Collect statistics
    total = 0
    mapped = 0
    unmapped = 0
    primary = 0
    secondary = 0
    supplementary = 0
    forward = 0
    reverse = 0
    paired = 0
    proper_pair = 0
    first = 0
    second = 0
    qc_fail = 0
    duplicate = 0

    mapq_30_plus = 0

    print("\nScanning records...")
    for record in bam:
        total += 1

        # Count by mapping status
        if record.is_mapped:
            mapped += 1
        else:
            unmapped += 1

        # Count by alignment type
        if record.is_primary:
            primary += 1
        if record.is_secondary:
            secondary += 1
        if record.is_supplementary:
            supplementary += 1

        # Count by strand
        if record.is_forward:
            forward += 1
        if record.is_reverse:
            reverse += 1

        # Count by pairing
        if record.is_paired:
            paired += 1
        if record.is_proper_pair:
            proper_pair += 1
        if record.is_first:
            first += 1
        if record.is_second:
            second += 1

        # Count quality issues
        if record.is_qc_fail:
            qc_fail += 1
        if record.is_duplicate:
            duplicate += 1

        # Count high-quality alignments
        if record.mapq and record.mapq >= 30:
            mapq_30_plus += 1

        # Stop after processing enough for testing
        if total >= 10000:
            break

    # Print statistics
    print(f"\n📊 Statistics (first {total:,} records):\n")

    print(f"Mapping status:")
    print(f"  Mapped:    {mapped:6,} ({100*mapped/total:.1f}%)")
    print(f"  Unmapped:  {unmapped:6,} ({100*unmapped/total:.1f}%)")

    print(f"\nAlignment type:")
    print(f"  Primary:       {primary:6,} ({100*primary/total:.1f}%)")
    print(f"  Secondary:     {secondary:6,} ({100*secondary/total:.1f}%)")
    print(f"  Supplementary: {supplementary:6,} ({100*supplementary/total:.1f}%)")

    print(f"\nStrand:")
    print(f"  Forward: {forward:6,} ({100*forward/total:.1f}%)")
    print(f"  Reverse: {reverse:6,} ({100*reverse/total:.1f}%)")

    print(f"\nPairing:")
    print(f"  Paired:       {paired:6,} ({100*paired/total:.1f}%)")
    print(f"  Proper pair:  {proper_pair:6,} ({100*proper_pair/total:.1f}%)")
    print(f"  First in pair: {first:6,} ({100*first/total:.1f}%)")
    print(f"  Second in pair: {second:6,} ({100*second/total:.1f}%)")

    print(f"\nQuality:")
    print(f"  QC fail:    {qc_fail:6,} ({100*qc_fail/total:.1f}%)")
    print(f"  Duplicate:  {duplicate:6,} ({100*duplicate/total:.1f}%)")
    print(f"  MAPQ ≥ 30:  {mapq_30_plus:6,} ({100*mapq_30_plus/total:.1f}%)")

    # Test 2: Complex filtering example
    print(f"\n\n🔍 Complex Filter Example:")
    print(f"   Filter: high-quality primary forward-strand alignments")
    print(f"   Criteria: is_mapped AND is_primary AND is_forward AND mapq ≥ 30\n")

    bam2 = biometal.BamReader.from_path(bam_path)
    count = 0
    for record in bam2:
        if (record.is_mapped and
            record.is_primary and
            record.is_forward and
            record.mapq and record.mapq >= 30):
            if count < 3:
                print(f"   ✓ {record.name}: MAPQ={record.mapq}, forward strand")
            count += 1
        if count >= 1000:  # Limit for speed
            break

    print(f"\n   Found {count:,} matching records")

    # Test 3: Proper pair filtering
    print(f"\n\n👫 Proper Pair Filter Example:")
    bam3 = biometal.BamReader.from_path(bam_path)
    count = 0
    for record in bam3:
        if record.is_paired and record.is_proper_pair:
            if count < 3:
                print(f"   ✓ {record.name}: properly paired")
            count += 1
        if count >= 100:
            break

    print(f"\n   Found {count:,} properly paired reads")

    # Test 4: Sequence length filtering
    print(f"\n\n📏 Sequence Length Filter Example:")
    print(f"   Filter: reads with length 50-150 bp\n")

    bam4 = biometal.BamReader.from_path(bam_path)
    count = 0
    for record in bam4:
        if 50 <= record.sequence_length <= 150:
            if count < 3:
                print(f"   ✓ {record.name}: {record.sequence_length} bp")
            count += 1
        if count >= 100:
            break

    print(f"\n   Found {count:,} records with length 50-150 bp")

    print(f"\n\n✅ All filtering tests passed!")

if __name__ == "__main__":
    try:
        test_filtering()
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
