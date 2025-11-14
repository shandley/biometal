#!/usr/bin/env python3
"""Test TBI Python bindings"""

import sys
sys.path.insert(0, '/Users/scotthandley/Code/biometal/.venv/lib/python3.14/site-packages')

import biometal
import tempfile
import os
import struct

def create_minimal_tbi():
    """Create a minimal valid TBI file for testing"""
    data = bytearray()

    # Magic string "TBI\1"
    data.extend(b"TBI\x01")

    # n_ref = 2 (int32)
    data.extend(struct.pack('<i', 2))

    # format = 2 (VCF, int32)
    data.extend(struct.pack('<i', 2))

    # col_seq = 0 (int32)
    data.extend(struct.pack('<i', 0))

    # col_beg = 1 (int32)
    data.extend(struct.pack('<i', 1))

    # col_end = 0 (int32, same as beg for VCF)
    data.extend(struct.pack('<i', 0))

    # meta = '#' (int32)
    data.extend(struct.pack('<i', ord('#')))

    # skip = 0 (int32)
    data.extend(struct.pack('<i', 0))

    # l_nm = 10 (length of "chr1\0chr2\0", int32)
    data.extend(struct.pack('<i', 10))

    # names = "chr1\0chr2\0"
    data.extend(b"chr1\0chr2\0")

    # Index for chr1:
    # n_bin = 2 (int32)
    data.extend(struct.pack('<i', 2))

    # Bin 0:
    # bin_id = 0 (uint32)
    data.extend(struct.pack('<I', 0))

    # n_chunk = 1 (int32)
    data.extend(struct.pack('<i', 1))

    # Chunk: start = 0x1000 (uint64)
    data.extend(struct.pack('<Q', 0x1000))

    # Chunk: end = 0x2000 (uint64)
    data.extend(struct.pack('<Q', 0x2000))

    # Bin 4681 (smallest bin, level 5):
    # bin_id = 4681 (uint32)
    data.extend(struct.pack('<I', 4681))

    # n_chunk = 1 (int32)
    data.extend(struct.pack('<i', 1))

    # Chunk: start = 0x1500 (uint64)
    data.extend(struct.pack('<Q', 0x1500))

    # Chunk: end = 0x1800 (uint64)
    data.extend(struct.pack('<Q', 0x1800))

    # Linear index:
    # n_intv = 2 (int32)
    data.extend(struct.pack('<i', 2))

    # Interval 0: offset = 0x1000
    data.extend(struct.pack('<Q', 0x1000))

    # Interval 1: offset = 0x1500
    data.extend(struct.pack('<Q', 0x1500))

    # Index for chr2:
    # n_bin = 1 (int32)
    data.extend(struct.pack('<i', 1))

    # Bin 0:
    # bin_id = 0 (uint32)
    data.extend(struct.pack('<I', 0))

    # n_chunk = 1 (int32)
    data.extend(struct.pack('<i', 1))

    # Chunk: start = 0x3000 (uint64)
    data.extend(struct.pack('<Q', 0x3000))

    # Chunk: end = 0x4000 (uint64)
    data.extend(struct.pack('<Q', 0x4000))

    # Linear index:
    # n_intv = 1 (int32)
    data.extend(struct.pack('<i', 1))

    # Interval 0: offset = 0x3000
    data.extend(struct.pack('<Q', 0x3000))

    return bytes(data)

def test_tbi_load_and_query():
    """Test loading TBI index and querying regions"""

    # Create temporary TBI file
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.tbi', delete=False) as f:
        tbi_path = f.name
        f.write(create_minimal_tbi())

    try:
        # Load index
        print("Loading TBI index...")
        index = biometal.TbiIndex.from_path(tbi_path)
        print(f"✓ Loaded TBI index: {index}")

        # Check metadata
        assert index.format() == "Vcf", f"Expected VCF format, got {index.format()}"
        print(f"✓ Format: {index.format()}")

        assert index.col_seq() == 0, f"Expected col_seq=0, got {index.col_seq()}"
        print(f"✓ Sequence column: {index.col_seq()}")

        assert index.col_beg() == 1, f"Expected col_beg=1, got {index.col_beg()}"
        print(f"✓ Start column: {index.col_beg()}")

        assert index.meta_char() == '#', f"Expected meta_char='#', got {index.meta_char()}"
        print(f"✓ Comment character: '{index.meta_char()}'")

        # Check references
        refs = index.references()
        assert len(refs) == 2, f"Expected 2 references, got {len(refs)}"
        assert refs == ['chr1', 'chr2'], f"Unexpected references: {refs}"
        print(f"✓ References: {refs}")

        # Check __len__
        assert len(index) == 2, f"Expected len=2, got {len(index)}"
        print(f"✓ Index length: {len(index)}")

        # Check contains
        assert index.contains('chr1'), "chr1 should exist"
        assert index.contains('chr2'), "chr2 should exist"
        assert not index.contains('chr99'), "chr99 should not exist"
        print("✓ contains() works correctly")

        # Get reference info
        chr1_info = index.get_info('chr1')
        assert chr1_info is not None, "chr1 info should exist"
        assert chr1_info['name'] == 'chr1', f"Expected name='chr1', got {chr1_info['name']}"
        assert chr1_info['n_bins'] == 2, f"Expected 2 bins, got {chr1_info['n_bins']}"
        print(f"✓ chr1 info: {chr1_info}")

        chr2_info = index.get_info('chr2')
        assert chr2_info is not None, "chr2 info should exist"
        assert chr2_info['name'] == 'chr2'
        print(f"✓ chr2 info: {chr2_info}")

        # Query region on chr1
        chunks = index.query('chr1', 0, 100000)
        assert len(chunks) > 0, "Should return at least one chunk"
        print(f"✓ Query chr1:0-100000 returned {len(chunks)} chunks")

        for i, (start, end) in enumerate(chunks):
            assert isinstance(start, int), "Chunk start should be int"
            assert isinstance(end, int), "Chunk end should be int"
            assert start < end, "Chunk start should be < end"
            print(f"  Chunk {i}: {start:016x} - {end:016x}")

        # Query region on chr2
        chunks_chr2 = index.query('chr2', 0, 50000)
        assert len(chunks_chr2) > 0, "Should return chunks for chr2"
        print(f"✓ Query chr2:0-50000 returned {len(chunks_chr2)} chunks")

        # Test invalid queries
        try:
            index.query('chr99', 0, 100)
            assert False, "Should raise error for non-existent reference"
        except ValueError as e:
            assert "not found" in str(e).lower()
            print(f"✓ Correctly raises error for non-existent reference")

        try:
            index.query('chr1', 100, 50)
            assert False, "Should raise error for invalid range"
        except ValueError as e:
            assert "invalid range" in str(e).lower()
            print(f"✓ Correctly raises error for invalid range")

        print("\n✅ All TBI Python binding tests passed!")

    finally:
        # Cleanup
        os.unlink(tbi_path)

if __name__ == '__main__':
    test_tbi_load_and_query()
