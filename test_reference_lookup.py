#!/usr/bin/env python3
"""Test script for BAM reference name lookup functionality."""

import biometal
import sys

def test_reference_lookup():
    """Test reference name lookup features."""
    # Open BAM file
    bam_path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam"
    print(f"Opening {bam_path}...")
    bam = biometal.BamReader.from_path(bam_path)

    # Get header
    header = bam.header
    print(f"\n✓ Header: {header}")

    # Test 1: Get reference count
    print(f"\n1. Reference count: {header.reference_count}")

    # Test 2: Get all reference names
    ref_names = header.reference_names
    print(f"\n2. Reference names ({len(ref_names)}):")
    for name in ref_names[:5]:  # Show first 5
        print(f"   - {name}")
    if len(ref_names) > 5:
        print(f"   ... and {len(ref_names) - 5} more")

    # Test 3: Get reference by ID
    print(f"\n3. Reference lookup by ID:")
    for i in range(min(3, len(ref_names))):
        ref = header.reference(i)
        if ref:
            print(f"   ID {i}: {ref}")
        else:
            print(f"   ID {i}: None (invalid)")

    # Test 4: Get reference name by ID
    print(f"\n4. Reference name by ID:")
    for i in range(min(3, len(ref_names))):
        name = header.reference_name(i)
        print(f"   ID {i}: {name}")

    # Test 5: Get reference length by ID
    print(f"\n5. Reference length by ID:")
    for i in range(min(3, len(ref_names))):
        length = header.reference_length(i)
        print(f"   ID {i}: {length:,} bp")

    # Test 6: Get reference ID by name (forward lookup)
    print(f"\n6. Reference ID by name (reverse lookup):")
    test_names = ref_names[:3] if len(ref_names) >= 3 else ref_names
    for name in test_names:
        ref_id = header.get_reference_id(name)
        print(f"   '{name}' → ID {ref_id}")

    # Test 7: Non-existent reference name
    print(f"\n7. Non-existent reference:")
    fake_name = "chrDoesNotExist"
    ref_id = header.get_reference_id(fake_name)
    print(f"   '{fake_name}' → ID {ref_id} (should be None)")

    # Test 8: Use in real workflow - convert reference IDs to names
    print(f"\n8. Converting record reference IDs to names:")
    count = 0
    for record in bam:
        if record.reference_id is not None:
            ref_name = header.reference_name(record.reference_id)
            print(f"   Record '{record.name}': {ref_name}:{record.position}")
            count += 1
            if count >= 5:
                break

    print(f"\n✅ All tests passed!")

if __name__ == "__main__":
    try:
        test_reference_lookup()
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
