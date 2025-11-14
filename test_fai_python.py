#!/usr/bin/env python3
"""Test FAI Python bindings"""

import biometal
import tempfile
import os

def test_fai_build_and_query():
    """Test building FAI index and querying sequences"""

    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        fasta_path = f.name
        f.write('>chr1 first chromosome\n')
        f.write('ACGTACGTACGTACGTACGTACGTACGT\n')
        f.write('TGCATGCATGCATGCA\n')
        f.write('>chr2 second chromosome\n')
        f.write('GGGGCCCCAAAATTTT\n')
        f.write('>chr3 third chromosome\n')
        f.write('ATCGATCGATCGATCGATCGATCGATCG\n')
        f.write('ATCGATCGATCGATCGATCGATCGATCG\n')
        f.write('ATCGATCG\n')

    try:
        # Build index
        print("Building FAI index...")
        index = biometal.FaiIndex.build(fasta_path)
        print(f"✓ Built index with {len(index)} sequences")

        # Check sequences
        names = index.sequence_names()
        assert names == ['chr1', 'chr2', 'chr3'], f"Unexpected sequence names: {names}"
        print(f"✓ Sequence names: {names}")

        # Get sequence info
        chr1_info = index.get_info('chr1')
        print(f"✓ chr1 info: length={chr1_info['length']} bp")
        assert chr1_info['length'] == 44

        # Fetch entire sequence
        chr1_seq = index.fetch('chr1', fasta_path)
        assert len(chr1_seq) == 44, f"Expected 44 bp, got {len(chr1_seq)}"
        print(f"✓ Fetched chr1: {len(chr1_seq)} bp")

        # Fetch chr2
        chr2_seq = index.fetch('chr2', fasta_path)
        assert chr2_seq == 'GGGGCCCCAAAATTTT', f"Unexpected chr2 sequence: {chr2_seq}"
        print(f"✓ Fetched chr2: {chr2_seq}")

        # Fetch region
        region = index.fetch_region('chr1', 0, 10, fasta_path)
        assert region == 'ACGTACGTAC', f"Unexpected region: {region}"
        print(f"✓ Fetched chr1:0-10: {region}")

        # Fetch region crossing line boundary
        region = index.fetch_region('chr1', 28, 38, fasta_path)
        assert region == 'TGCATGCATG', f"Unexpected region: {region}"
        print(f"✓ Fetched chr1:28-38: {region}")

        # Check contains
        assert index.contains('chr1'), "chr1 should exist"
        assert not index.contains('chr99'), "chr99 should not exist"
        print("✓ contains() works correctly")

        # Write and reload index
        fai_path = fasta_path + '.fai'
        index.write(fai_path)
        print(f"✓ Wrote index to {fai_path}")

        # Load from file
        loaded_index = biometal.FaiIndex.from_path(fai_path)
        assert len(loaded_index) == 3, f"Expected 3 sequences, got {len(loaded_index)}"
        print(f"✓ Loaded index from file: {len(loaded_index)} sequences")

        # Verify loaded index works
        chr2_reload = loaded_index.fetch('chr2', fasta_path)
        assert chr2_reload == chr2_seq, "Reloaded index gives different results"
        print("✓ Reloaded index works correctly")

        # Cleanup
        os.unlink(fai_path)

        print("\n✅ All FAI Python binding tests passed!")

    finally:
        # Cleanup
        os.unlink(fasta_path)

if __name__ == '__main__':
    test_fai_build_and_query()
