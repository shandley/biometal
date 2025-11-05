#!/usr/bin/env python3
"""Test script for biometal Python bindings."""

import sys
import biometal

def test_version():
    """Test module version."""
    print(f"biometal version: {biometal.__version__}")
    assert biometal.__version__ == "0.2.3"
    print("✓ Version test passed")

def test_fastq_stream():
    """Test FASTQ streaming."""
    print("\nTesting FASTQ streaming...")
    stream = biometal.FastqStream.from_path("tests/test_data.fq.gz")

    records = []
    for record in stream:
        records.append(record)
        print(f"  Read: {record.id}")
        print(f"    Sequence: {record.sequence_str}")
        print(f"    Quality: {record.quality_str}")

    assert len(records) == 3
    assert records[0].id == "read1"
    assert records[0].sequence_str == "ATGCATGC"
    assert records[1].id == "read2"
    assert records[2].id == "read3"
    assert len(records[2].sequence) == 16

    print("✓ FASTQ streaming test passed")

def test_fasta_stream():
    """Test FASTA streaming."""
    print("\nTesting FASTA streaming...")
    stream = biometal.FastaStream.from_path("tests/test_data.fa.gz")

    records = []
    for record in stream:
        records.append(record)
        print(f"  Sequence: {record.id}")
        print(f"    Length: {len(record.sequence)} bp")

    assert len(records) == 3
    assert records[0].id == "seq1"
    assert records[0].sequence_str == "ATGCATGCATGC"
    assert records[1].id == "seq2"
    assert records[2].id == "seq3"

    print("✓ FASTA streaming test passed")

def test_gc_content():
    """Test GC content calculation."""
    print("\nTesting GC content...")

    # 50% GC (4 G/C out of 8)
    seq1 = b"ATGCATGC"
    gc1 = biometal.gc_content(seq1)
    print(f"  {seq1.decode()}: {gc1:.2%}")
    assert abs(gc1 - 0.5) < 0.001

    # 100% GC
    seq2 = b"GCGCGCGC"
    gc2 = biometal.gc_content(seq2)
    print(f"  {seq2.decode()}: {gc2:.2%}")
    assert abs(gc2 - 1.0) < 0.001

    # 0% GC
    seq3 = b"ATATATAT"
    gc3 = biometal.gc_content(seq3)
    print(f"  {seq3.decode()}: {gc3:.2%}")
    assert abs(gc3 - 0.0) < 0.001

    print("✓ GC content test passed")

def test_count_bases():
    """Test base counting."""
    print("\nTesting base counting...")

    seq = b"ATGCATGC"
    counts = biometal.count_bases(seq)
    print(f"  Sequence: {seq.decode()}")
    print(f"  Counts: A={counts['A']}, C={counts['C']}, G={counts['G']}, T={counts['T']}")

    assert counts['A'] == 2
    assert counts['T'] == 2
    assert counts['G'] == 2
    assert counts['C'] == 2

    print("✓ Base counting test passed")

def test_mean_quality():
    """Test mean quality calculation."""
    print("\nTesting mean quality...")

    # Phred+33 encoding
    # I = ASCII 73, Phred = 73 - 33 = 40
    qual = b"IIIIIIII"
    mean_q = biometal.mean_quality(qual)
    print(f"  Quality string: {qual.decode()}")
    print(f"  Mean quality: {mean_q:.1f}")
    assert abs(mean_q - 40.0) < 0.1

    print("✓ Mean quality test passed")

def test_extract_kmers():
    """Test k-mer extraction."""
    print("\nTesting k-mer extraction...")

    seq = b"ATGCATGC"

    # 3-mers (overlapping)
    kmers_3 = biometal.extract_kmers(seq, 3)
    print(f"  Sequence: {seq.decode()}")
    print(f"  3-mers: {kmers_3}")
    assert len(kmers_3) == 6  # 8 - 3 + 1
    assert kmers_3[0] == "ATG"
    assert kmers_3[-1] == "TGC"

    # 4-mers (overlapping)
    kmers_4 = biometal.extract_kmers(seq, 4)
    print(f"  4-mers: {kmers_4}")
    assert len(kmers_4) == 5  # 8 - 4 + 1

    # Non-overlapping 4-mers
    kmers_non_overlap = biometal.extract_kmers_non_overlapping(seq, 4)
    print(f"  4-mers (non-overlapping): {kmers_non_overlap}")
    assert len(kmers_non_overlap) == 2  # 8 / 4
    assert kmers_non_overlap[0] == "ATGC"
    assert kmers_non_overlap[1] == "ATGC"

    print("✓ K-mer extraction test passed")

def test_integration():
    """Test integrated workflow."""
    print("\nTesting integrated workflow...")

    stream = biometal.FastqStream.from_path("tests/test_data.fq.gz")

    total_bases = 0
    total_gc = 0.0

    for record in stream:
        # Convert to bytes for operations
        seq_bytes = bytes(record.sequence)
        qual_bytes = bytes(record.quality)

        # Calculate GC content
        gc = biometal.gc_content(seq_bytes)

        # Count bases
        counts = biometal.count_bases(seq_bytes)
        base_count = sum(counts.values())

        # Mean quality
        mean_q = biometal.mean_quality(qual_bytes)

        # Extract k-mers
        kmers = biometal.extract_kmers(seq_bytes, 3)

        print(f"  {record.id}:")
        print(f"    Bases: {base_count}")
        print(f"    GC: {gc:.2%}")
        print(f"    Mean Q: {mean_q:.1f}")
        print(f"    3-mers: {len(kmers)}")

        total_bases += base_count
        total_gc += gc * base_count

    avg_gc = total_gc / total_bases if total_bases > 0 else 0
    print(f"\n  Total bases: {total_bases}")
    print(f"  Average GC: {avg_gc:.2%}")

    print("✓ Integration test passed")

def main():
    """Run all tests."""
    print("=" * 60)
    print("biometal Python Bindings Test Suite")
    print("=" * 60)

    try:
        test_version()
        test_fastq_stream()
        test_fasta_stream()
        test_gc_content()
        test_count_bases()
        test_mean_quality()
        test_extract_kmers()
        test_integration()

        print("\n" + "=" * 60)
        print("✓ ALL TESTS PASSED")
        print("=" * 60)
        return 0

    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
