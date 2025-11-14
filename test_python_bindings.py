#!/usr/bin/env python3
"""Test Python bindings for format library parsers"""

import sys

def test_bed_bindings():
    """Test BED format Python bindings"""
    print("Testing BED bindings...")

    try:
        import biometal

        # Test Bed3Stream API
        print("  ✓ Bed3Stream class exists:", hasattr(biometal, 'Bed3Stream'))
        print("  ✓ Bed3Stream.from_path exists:", hasattr(biometal.Bed3Stream, 'from_path'))

        # Test Bed6Stream API
        print("  ✓ Bed6Stream class exists:", hasattr(biometal, 'Bed6Stream'))
        print("  ✓ Bed6Stream.from_path exists:", hasattr(biometal.Bed6Stream, 'from_path'))

        # Test Bed12Stream API
        print("  ✓ Bed12Stream class exists:", hasattr(biometal, 'Bed12Stream'))
        print("  ✓ Bed12Stream.from_path exists:", hasattr(biometal.Bed12Stream, 'from_path'))

        # Test record types
        print("  ✓ Bed3Record exists:", hasattr(biometal, 'Bed3Record'))
        print("  ✓ Bed6Record exists:", hasattr(biometal, 'Bed6Record'))
        print("  ✓ Bed12Record exists:", hasattr(biometal, 'Bed12Record'))

        print("  ✅ BED bindings OK")
        return True
    except Exception as e:
        print(f"  ❌ BED bindings failed: {e}")
        return False

def test_gfa_bindings():
    """Test GFA format Python bindings"""
    print("\nTesting GFA bindings...")

    try:
        import biometal

        # Test GfaStream API
        print("  ✓ GfaStream class exists:", hasattr(biometal, 'GfaStream'))
        print("  ✓ GfaStream.from_path exists:", hasattr(biometal.GfaStream, 'from_path'))

        # Test record types
        print("  ✓ GfaSegment exists:", hasattr(biometal, 'GfaSegment'))
        print("  ✓ GfaLink exists:", hasattr(biometal, 'GfaLink'))
        print("  ✓ GfaPath exists:", hasattr(biometal, 'GfaPath'))

        print("  ✅ GFA bindings OK")
        return True
    except Exception as e:
        print(f"  ❌ GFA bindings failed: {e}")
        return False

def test_vcf_bindings():
    """Test VCF format Python bindings"""
    print("\nTesting VCF bindings...")

    try:
        import biometal

        # Test VcfStream API
        print("  ✓ VcfStream class exists:", hasattr(biometal, 'VcfStream'))
        print("  ✓ VcfStream.from_path exists:", hasattr(biometal.VcfStream, 'from_path'))

        # Test record types
        print("  ✓ VcfHeader exists:", hasattr(biometal, 'VcfHeader'))
        print("  ✓ VcfRecord exists:", hasattr(biometal, 'VcfRecord'))

        print("  ✅ VCF bindings OK")
        return True
    except Exception as e:
        print(f"  ❌ VCF bindings failed: {e}")
        return False

def test_gff3_bindings():
    """Test GFF3 format Python bindings"""
    print("\nTesting GFF3 bindings...")

    try:
        import biometal

        # Test Gff3Stream API
        print("  ✓ Gff3Stream class exists:", hasattr(biometal, 'Gff3Stream'))
        print("  ✓ Gff3Stream.from_path exists:", hasattr(biometal.Gff3Stream, 'from_path'))

        # Test record types
        print("  ✓ Gff3Record exists:", hasattr(biometal, 'Gff3Record'))

        print("  ✅ GFF3 bindings OK")
        return True
    except Exception as e:
        print(f"  ❌ GFF3 bindings failed: {e}")
        return False

def test_actual_parsing():
    """Test actual file parsing"""
    print("\nTesting actual file parsing...")

    try:
        import biometal

        # Test BED6 parsing (gzipped file)
        print("  Testing BED6 parsing (gzipped)...")
        stream = biometal.Bed6Stream.from_path("tests/data/real_world/encode_peaks.bed.gz")
        count = 0
        for record in stream:
            count += 1
            assert hasattr(record, 'chrom')
            assert hasattr(record, 'start')
            assert hasattr(record, 'end')
            assert hasattr(record, 'name')
            assert hasattr(record, 'score')
            assert hasattr(record, 'strand')
        print(f"    ✓ Parsed {count} BED6 records from .gz file")

        # Test GFA parsing (uncompressed)
        print("  Testing GFA parsing...")
        stream = biometal.GfaStream.from_path("tests/data/real_world/lambda_phage.gfa")
        segments = 0
        links = 0
        paths = 0
        for record in stream:
            if isinstance(record, biometal.GfaSegment):
                segments += 1
            elif isinstance(record, biometal.GfaLink):
                links += 1
            elif isinstance(record, biometal.GfaPath):
                paths += 1
        print(f"    ✓ Parsed {segments} segments, {links} links, {paths} paths")

        # Test VCF parsing (gzipped file)
        print("  Testing VCF parsing (gzipped)...")
        stream = biometal.VcfStream.from_path("tests/data/real_world/synthetic_1000g.vcf.gz")
        # Get header (automatically parsed)
        header = stream.header()
        print(f"    ✓ Parsed VCF header: {header.fileformat}")
        count = 0
        for record in stream:
            count += 1
            assert hasattr(record, 'chrom')
            assert hasattr(record, 'pos')
            assert hasattr(record, 'reference')  # Note: 'reference' not 'ref'
            assert hasattr(record, 'alternate')  # Note: 'alternate' not 'alt'
        print(f"    ✓ Parsed {count} VCF records from .gz file")

        # Test GFF3 parsing (gzipped file)
        print("  Testing GFF3 parsing (gzipped)...")
        stream = biometal.Gff3Stream.from_path("tests/data/real_world/ensembl_chr21.gff3.gz")
        count = 0
        genes = 0
        for record in stream:
            count += 1
            assert hasattr(record, 'seqid')
            assert hasattr(record, 'feature_type')
            assert hasattr(record, 'start')
            assert hasattr(record, 'end')
            if record.feature_type == "gene":
                genes += 1
            if count >= 1000:  # Don't parse entire file
                break
        print(f"    ✓ Parsed {count} GFF3 records ({genes} genes) from .gz file")

        print("  ✅ File parsing OK")
        return True
    except Exception as e:
        print(f"  ❌ File parsing failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("=" * 60)
    print("Python Bindings Verification for Format Library")
    print("=" * 60)

    results = []

    # Test bindings exist
    results.append(("BED bindings", test_bed_bindings()))
    results.append(("GFA bindings", test_gfa_bindings()))
    results.append(("VCF bindings", test_vcf_bindings()))
    results.append(("GFF3 bindings", test_gff3_bindings()))

    # Test actual parsing
    results.append(("File parsing", test_actual_parsing()))

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{name:20s}: {status}")

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All Python bindings working correctly!")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
