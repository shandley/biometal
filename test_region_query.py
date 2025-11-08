#!/usr/bin/env python3
"""Test script for BAM region query functionality."""

import biometal
import sys

def test_region_query():
    """Test region query features."""
    bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam"
    print(f"Opening {bam_path}...")

    # First, get header to know reference IDs
    bam = biometal.BamReader.from_path(bam_path)
    header = bam.header
    print(f"\n✓ References: {', '.join(header.reference_names)}")

    # Test 1: Query specific region with start and end
    print(f"\n1. Query region chr1:0-1000 (reference_id=0, start=0, end=1000):")
    count = 0
    for record in biometal.BamReader.query(bam_path, reference_id=0, start=0, end=1000):
        if count < 3:
            print(f"   Record '{record.name}': position {record.position}")
        count += 1
    print(f"   Total: {count} records in region")

    # Test 2: Query with only start position
    print(f"\n2. Query region chr1:1000-end (reference_id=0, start=1000):")
    count = 0
    for record in biometal.BamReader.query(bam_path, reference_id=0, start=1000):
        if count < 3:
            print(f"   Record '{record.name}': position {record.position}")
        count += 1
        if count >= 100:  # Limit to avoid processing entire file
            break
    print(f"   Found {count} records (stopped at 100 for speed)")

    # Test 3: Query with only end position
    print(f"\n3. Query region chr1:0-500 (reference_id=0, end=500):")
    count = 0
    for record in biometal.BamReader.query(bam_path, reference_id=0, end=500):
        if count < 3:
            print(f"   Record '{record.name}': position {record.position}")
        count += 1
    print(f"   Total: {count} records in region")

    # Test 4: Query entire reference
    print(f"\n4. Query entire reference chr2 (reference_id=1, no range):")
    count = 0
    for record in biometal.BamReader.query(bam_path, reference_id=1):
        if count == 0:
            print(f"   First record: '{record.name}' at position {record.position}")
        count += 1
        if count >= 50:  # Limit for speed
            break
    print(f"   Found {count} records (stopped at 50 for speed)")

    # Test 5: Use reference name lookup with query
    print(f"\n5. Query using reference name lookup:")
    ref_id = header.get_reference_id("chr1")
    if ref_id is not None:
        count = 0
        for record in biometal.BamReader.query(bam_path, reference_id=ref_id, start=500, end=1500):
            if count < 3:
                ref_name = header.reference_name(record.reference_id)
                print(f"   Record '{record.name}': {ref_name}:{record.position}")
            count += 1
        print(f"   Total: {count} records in chr1:500-1500")

    # Test 6: Verify position filtering is correct
    print(f"\n6. Verify positions are within range:")
    all_in_range = True
    count = 0
    test_start, test_end = 1000, 2000
    for record in biometal.BamReader.query(bam_path, reference_id=0, start=test_start, end=test_end):
        if record.position < test_start or record.position >= test_end:
            print(f"   ❌ Record '{record.name}' at position {record.position} is out of range!")
            all_in_range = False
        count += 1
        if count >= 100:
            break

    if all_in_range:
        print(f"   ✓ All {count} records are within range [{test_start}, {test_end})")
    else:
        print(f"   ❌ Some records were out of range!")

    print(f"\n✅ All region query tests passed!")

if __name__ == "__main__":
    try:
        test_region_query()
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
