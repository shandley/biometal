#!/usr/bin/env python3
"""
Basic usage examples for biometal Python bindings.

This script demonstrates the core functionality:
- FASTQ/FASTA streaming
- NEON-accelerated operations
- K-mer extraction
"""

import biometal

def example_fastq_streaming():
    """Stream and analyze a FASTQ file."""
    print("=" * 60)
    print("Example 1: FASTQ Streaming")
    print("=" * 60)

    stream = biometal.FastqStream.from_path("../tests/test_data.fq.gz")

    for record in stream:
        seq_bytes = bytes(record.sequence)
        qual_bytes = bytes(record.quality)

        # Calculate metrics
        gc = biometal.gc_content(seq_bytes)
        mean_q = biometal.mean_quality(qual_bytes)
        counts = biometal.count_bases(seq_bytes)

        print(f"\nRead: {record.id}")
        print(f"  Length: {len(record.sequence)} bp")
        print(f"  GC content: {gc:.2%}")
        print(f"  Mean quality: {mean_q:.1f}")
        print(f"  Base counts: {counts}")


def example_fasta_streaming():
    """Stream and analyze a FASTA file."""
    print("\n" + "=" * 60)
    print("Example 2: FASTA Streaming")
    print("=" * 60)

    stream = biometal.FastaStream.from_path("../tests/test_data.fa.gz")

    for record in stream:
        seq_bytes = bytes(record.sequence)
        gc = biometal.gc_content(seq_bytes)

        print(f"\nSequence: {record.id}")
        print(f"  Length: {len(record.sequence)} bp")
        print(f"  GC content: {gc:.2%}")
        print(f"  Sequence: {record.sequence_str}")


def example_quality_filter():
    """Filter reads by quality score."""
    print("\n" + "=" * 60)
    print("Example 3: Quality Filtering")
    print("=" * 60)

    stream = biometal.FastqStream.from_path("../tests/test_data.fq.gz")

    high_quality = []
    low_quality = []

    for record in stream:
        qual_bytes = bytes(record.quality)
        mean_q = biometal.mean_quality(qual_bytes)

        if mean_q >= 30.0:
            high_quality.append(record.id)
        else:
            low_quality.append(record.id)

    print(f"\nHigh quality (Q≥30): {len(high_quality)} reads")
    for read_id in high_quality:
        print(f"  ✓ {read_id}")

    print(f"\nLow quality (Q<30): {len(low_quality)} reads")
    for read_id in low_quality:
        print(f"  ✗ {read_id}")


def example_kmer_extraction():
    """Extract k-mers for ML preprocessing."""
    print("\n" + "=" * 60)
    print("Example 4: K-mer Extraction for ML")
    print("=" * 60)

    stream = biometal.FastqStream.from_path("../tests/test_data.fq.gz")

    print("\nExtracting 3-mers from high-quality reads...")

    all_kmers = []
    for record in stream:
        seq_bytes = bytes(record.sequence)
        qual_bytes = bytes(record.quality)

        # Only use high-quality reads
        mean_q = biometal.mean_quality(qual_bytes)
        if mean_q >= 30.0:
            kmers = biometal.extract_kmers(seq_bytes, 3)
            all_kmers.extend(kmers)
            print(f"\n{record.id}: {len(kmers)} k-mers")
            print(f"  First 5: {kmers[:5]}")

    # Count frequencies
    from collections import Counter
    kmer_counts = Counter(all_kmers)

    print(f"\nTotal unique k-mers: {len(kmer_counts)}")
    print("Most common k-mers:")
    for kmer, count in kmer_counts.most_common(5):
        print(f"  {kmer}: {count}")


def example_gc_distribution():
    """Calculate GC content distribution."""
    print("\n" + "=" * 60)
    print("Example 5: GC Content Distribution")
    print("=" * 60)

    stream = biometal.FastqStream.from_path("../tests/test_data.fq.gz")

    gc_values = []
    for record in stream:
        seq_bytes = bytes(record.sequence)
        gc = biometal.gc_content(seq_bytes)
        gc_values.append(gc)

    avg_gc = sum(gc_values) / len(gc_values)
    min_gc = min(gc_values)
    max_gc = max(gc_values)

    print(f"\nGC Content Statistics:")
    print(f"  Average: {avg_gc:.2%}")
    print(f"  Min: {min_gc:.2%}")
    print(f"  Max: {max_gc:.2%}")
    print(f"  Range: {(max_gc - min_gc):.2%}")


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("biometal Python Bindings - Basic Usage Examples")
    print(f"Version: {biometal.__version__}")
    print("=" * 60)

    example_fastq_streaming()
    example_fasta_streaming()
    example_quality_filter()
    example_kmer_extraction()
    example_gc_distribution()

    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
