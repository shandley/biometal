#!/usr/bin/env python3
"""
Comprehensive test suite for BAM/SAM Python bindings

Tests all BAM functionality including:
- BamReader streaming
- Header and reference access
- Record fields and properties
- CIGAR operations
- Tag parsing
- Region queries
- Statistics functions
- Flag properties
"""
import pytest
from pathlib import Path
import biometal


# Test data path
BAM_FILE = Path("experiments/native-bam-implementation/test-data/synthetic_100000.bam")


@pytest.fixture
def bam_path():
    """Fixture providing path to test BAM file"""
    if not BAM_FILE.exists():
        pytest.skip(f"Test BAM file not found: {BAM_FILE}")
    return str(BAM_FILE)


@pytest.fixture
def bam_reader(bam_path):
    """Fixture providing BamReader instance"""
    return biometal.BamReader.from_path(bam_path)


# ============================================================================
# BamReader Tests
# ============================================================================

def test_bam_reader_creation(bam_path):
    """Test BamReader can be created from path"""
    reader = biometal.BamReader.from_path(bam_path)
    assert reader is not None


def test_bam_reader_invalid_path():
    """Test BamReader raises error for invalid path"""
    with pytest.raises(Exception):  # Should raise IOError or similar
        biometal.BamReader.from_path("/nonexistent/file.bam")


def test_bam_reader_iteration(bam_reader):
    """Test BamReader can be iterated"""
    count = 0
    for record in bam_reader:
        count += 1
        if count >= 10:
            break
    assert count == 10


def test_bam_reader_repr(bam_reader):
    """Test BamReader has string representation"""
    repr_str = repr(bam_reader)
    assert "BamReader" in repr_str
    assert "references" in repr_str


# ============================================================================
# Header Tests
# ============================================================================

def test_header_access(bam_reader):
    """Test header can be accessed"""
    header = bam_reader.header
    assert header is not None


def test_header_reference_count(bam_reader):
    """Test header has reference count"""
    header = bam_reader.header
    assert header.reference_count >= 0


def test_header_reference_lookup(bam_reader):
    """Test reference lookup by index"""
    header = bam_reader.header
    if header.reference_count > 0:
        ref = header.reference(0)
        assert ref is not None
        assert isinstance(ref.name, str)
        assert ref.length > 0


def test_header_reference_by_name(bam_reader):
    """Test reference lookup by name - currently tests reference_name() method"""
    header = bam_reader.header
    if header.reference_count > 0:
        name = header.reference_name(0)
        assert name is not None
        assert isinstance(name, str)


def test_header_reference_names(bam_reader):
    """Test reference_names property"""
    header = bam_reader.header
    names = header.reference_names
    assert isinstance(names, list)
    assert len(names) == header.reference_count
    if len(names) > 0:
        assert all(isinstance(name, str) for name in names)


def test_reference_repr():
    """Test Reference has string representation"""
    # We can't create Reference directly, but we can test the format
    # This will be tested via header.reference()
    pass


# ============================================================================
# Record Field Tests
# ============================================================================

def test_record_fields(bam_reader):
    """Test record has all expected fields"""
    for record in bam_reader:
        # Required fields
        assert isinstance(record.name, str)
        assert isinstance(record.flags, int)
        assert isinstance(record.sequence, bytes)
        assert isinstance(record.quality, bytes)
        assert isinstance(record.template_length, int)

        # Optional fields (may be None)
        assert record.reference_id is None or isinstance(record.reference_id, int)
        assert record.position is None or isinstance(record.position, int)
        assert record.mapq is None or isinstance(record.mapq, int)
        assert record.mate_reference_id is None or isinstance(record.mate_reference_id, int)
        assert record.mate_position is None or isinstance(record.mate_position, int)
        break


def test_record_sequence_str(bam_reader):
    """Test sequence can be converted to string"""
    for record in bam_reader:
        seq_str = record.sequence_str  # Property, not method
        assert isinstance(seq_str, str)
        assert len(seq_str) == len(record.sequence)
        # Check all characters are valid DNA bases
        valid_bases = set('ACGTNacgtn')
        assert all(c in valid_bases for c in seq_str)
        break


def test_record_quality_str(bam_reader):
    """Test quality can be converted to string"""
    for record in bam_reader:
        qual_str = record.quality_str  # Property, not method
        assert isinstance(qual_str, str)
        assert len(qual_str) == len(record.quality)
        break


def test_record_repr(bam_reader):
    """Test record has string representation"""
    for record in bam_reader:
        repr_str = repr(record)
        assert "BamRecord" in repr_str
        assert record.name in repr_str
        break


# ============================================================================
# CIGAR Tests
# ============================================================================

def test_cigar_access(bam_reader):
    """Test CIGAR operations can be accessed"""
    for record in bam_reader:
        cigar = record.cigar
        assert cigar is not None
        assert isinstance(cigar, list)
        if len(cigar) > 0:
            assert hasattr(cigar[0], 'length')
            assert hasattr(cigar[0], 'op_char')
        break


def test_cigar_string(bam_reader):
    """Test CIGAR string generation"""
    for record in bam_reader:
        cigar_str = record.cigar_string()
        assert isinstance(cigar_str, str)
        # Should contain numbers and operation characters
        if len(cigar_str) > 0:
            assert any(c.isdigit() for c in cigar_str)
            assert any(c in 'MIDNSHP=X' for c in cigar_str)
        break


def test_cigar_op_properties(bam_reader):
    """Test CigarOp properties"""
    for record in bam_reader:
        cigar = record.cigar
        if len(cigar) > 0:
            op = cigar[0]

            # Test length
            assert isinstance(op.length, int)
            assert op.length > 0

            # Test op_char
            assert isinstance(op.op_char, str)
            assert len(op.op_char) == 1
            assert op.op_char in 'MIDNSHP=X'

            # Test string representation
            op_str = str(op)
            assert str(op.length) in op_str
            assert op.op_char in op_str
            break


def test_cigar_op_type_checking(bam_reader):
    """Test CigarOp type checking methods"""
    for record in bam_reader:
        cigar = record.cigar
        if len(cigar) > 0:
            op = cigar[0]

            # At least one should be True
            type_checks = [
                op.is_match(),
                op.is_insertion(),
                op.is_deletion(),
                op.is_ref_skip(),
                op.is_soft_clip(),
                op.is_hard_clip(),
                op.is_padding(),
                op.is_seq_match(),
                op.is_seq_mismatch()
            ]
            assert any(type_checks)
            break


def test_cigar_consumption_methods(bam_reader):
    """Test CigarOp consumption methods"""
    for record in bam_reader:
        cigar = record.cigar
        if len(cigar) > 0:
            op = cigar[0]

            # Both should return boolean
            assert isinstance(op.consumes_reference(), bool)
            assert isinstance(op.consumes_query(), bool)

            # Match operations should consume both
            if op.is_match():
                assert op.consumes_reference()
                assert op.consumes_query()

            # Insertions consume query but not reference
            if op.is_insertion():
                assert not op.consumes_reference()
                assert op.consumes_query()

            # Deletions consume reference but not query
            if op.is_deletion():
                assert op.consumes_reference()
                assert not op.consumes_query()
            break


def test_reference_length(bam_reader):
    """Test reference_length calculation"""
    for record in bam_reader:
        ref_len = record.reference_length()
        assert isinstance(ref_len, int)
        assert ref_len >= 0

        # Should equal sum of reference-consuming operations
        cigar = record.cigar
        expected = sum(
            op.length for op in cigar
            if op.consumes_reference()
        )
        assert ref_len == expected
        break


def test_query_length(bam_reader):
    """Test query_length calculation"""
    for record in bam_reader:
        query_len = record.query_length()
        assert isinstance(query_len, int)
        assert query_len >= 0

        # Should equal sum of query-consuming operations
        cigar = record.cigar
        expected = sum(
            op.length for op in cigar
            if op.consumes_query()
        )
        assert query_len == expected
        break


def test_reference_end(bam_reader):
    """Test reference_end calculation"""
    for record in bam_reader:
        ref_end = record.reference_end()

        if record.position is not None:
            assert ref_end is not None
            assert isinstance(ref_end, int)
            assert ref_end >= record.position

            # Should equal position + reference_length
            expected = record.position + record.reference_length()
            assert ref_end == expected
        else:
            # Unmapped records should return None
            assert ref_end is None
        break


# ============================================================================
# Flag Property Tests
# ============================================================================

def test_flag_properties(bam_reader):
    """Test all flag property methods"""
    for record in bam_reader:
        # All flag properties should return boolean
        assert isinstance(record.is_paired, bool)
        assert isinstance(record.is_proper_pair, bool)
        assert isinstance(record.is_mapped, bool)
        assert isinstance(record.is_mate_mapped, bool)
        assert isinstance(record.is_reverse, bool)
        assert isinstance(record.is_mate_reverse, bool)
        assert isinstance(record.is_first, bool)
        assert isinstance(record.is_second, bool)
        assert isinstance(record.is_secondary, bool)
        assert isinstance(record.is_qc_fail, bool)
        assert isinstance(record.is_duplicate, bool)
        assert isinstance(record.is_supplementary, bool)

        # Test is_forward (should be opposite of is_reverse)
        assert record.is_forward == (not record.is_reverse)

        # Test is_primary (should be opposite of is_secondary and is_supplementary)
        assert record.is_primary == (not record.is_secondary and not record.is_supplementary)
        break


def test_flag_consistency(bam_reader):
    """Test flag properties are consistent with flags field"""
    for record in bam_reader:
        # is_paired should match bit 0x1
        assert record.is_paired == bool(record.flags & 0x1)

        # is_mapped should be opposite of bit 0x4
        assert record.is_mapped == (not bool(record.flags & 0x4))

        # is_reverse should match bit 0x10
        assert record.is_reverse == bool(record.flags & 0x10)
        break


# ============================================================================
# Tag Tests
# ============================================================================

def test_tag_methods_exist(bam_reader):
    """Test tag access methods exist"""
    for record in bam_reader:
        # Methods should exist
        assert hasattr(record, 'get_tag')
        assert hasattr(record, 'has_tag')
        assert hasattr(record, 'tags')
        break


def test_tags_access(bam_reader):
    """Test tags() method returns list"""
    for record in bam_reader:
        tags = record.tags()
        assert isinstance(tags, list)
        # Synthetic data may not have tags
        break


def test_has_tag(bam_reader):
    """Test has_tag() method"""
    for record in bam_reader:
        # Should return boolean
        has_nm = record.has_tag("NM")
        assert isinstance(has_nm, bool)
        break


def test_get_tag(bam_reader):
    """Test get_tag() method"""
    for record in bam_reader:
        tag = record.get_tag("NM")
        # May be None if tag doesn't exist
        if tag is not None:
            assert hasattr(tag, 'name')
            assert hasattr(tag, 'value')
        break


# ============================================================================
# Region Query Tests
# ============================================================================

def test_region_query(bam_path):
    """Test region query functionality"""
    # Query first 1000 bases of first reference
    region_iter = biometal.BamReader.query(bam_path, 0, 0, 1000)
    assert region_iter is not None

    count = 0
    for record in region_iter:
        count += 1
        # Records should be from reference 0
        if record.reference_id is not None:
            assert record.reference_id == 0
        # Records should overlap the query region
        if record.position is not None:
            # Position should be before region end
            assert record.position < 1000

    # Should have found some records
    assert count > 0


def test_region_query_full_reference(bam_path):
    """Test region query with no start/end (full reference)"""
    region_iter = biometal.BamReader.query(bam_path, 0, None, None)

    count = 0
    for record in region_iter:
        count += 1
        if count >= 10:
            break

    assert count > 0


# ============================================================================
# Statistics Function Tests
# ============================================================================

def test_calculate_coverage(bam_path):
    """Test calculate_coverage function"""
    coverage = biometal.calculate_coverage(bam_path, 0, 0, 1000)

    assert isinstance(coverage, dict)
    # Should have integer keys (positions) and integer values (counts)
    for pos, count in coverage.items():
        assert isinstance(pos, int)
        assert isinstance(count, int)
        assert count > 0


def test_mapq_distribution(bam_path):
    """Test mapq_distribution function"""
    mapq_dist = biometal.mapq_distribution(bam_path)

    assert isinstance(mapq_dist, dict)
    # Should have integer keys (MAPQ) and integer values (counts)
    for mapq, count in mapq_dist.items():
        assert isinstance(mapq, int)
        assert 0 <= mapq <= 255
        assert isinstance(count, int)
        assert count > 0


def test_mapq_distribution_filtered(bam_path):
    """Test mapq_distribution with reference filter"""
    mapq_dist = biometal.mapq_distribution(bam_path, reference_id=0)

    assert isinstance(mapq_dist, dict)
    # Should have fewer or equal records compared to unfiltered
    total_filtered = sum(mapq_dist.values())
    assert total_filtered > 0


def test_count_by_flag(bam_path):
    """Test count_by_flag function"""
    flag_counts = biometal.count_by_flag(bam_path)

    assert isinstance(flag_counts, dict)

    # Should have expected keys (based on actual implementation)
    expected_keys = [
        'total', 'mapped', 'unmapped', 'primary', 'secondary',
        'supplementary', 'paired', 'proper_pair', 'forward',
        'reverse', 'qc_fail', 'duplicate'
    ]

    for key in expected_keys:
        assert key in flag_counts, f"Missing key: {key}"
        assert isinstance(flag_counts[key], int)
        assert flag_counts[key] >= 0

    # Total should be largest
    assert flag_counts['total'] >= flag_counts['paired']
    assert flag_counts['total'] >= flag_counts['mapped']
    # Mapped + unmapped should equal total
    assert flag_counts['mapped'] + flag_counts['unmapped'] == flag_counts['total']


# ============================================================================
# Integration Tests
# ============================================================================

def test_full_workflow(bam_reader):
    """Test complete workflow: read, filter, analyze"""
    high_quality_count = 0
    total_ref_length = 0

    for record in bam_reader:
        # Filter high-quality mapped reads
        if record.is_mapped and record.mapq is not None and record.mapq >= 30:
            high_quality_count += 1
            total_ref_length += record.reference_length()

        if high_quality_count >= 10:
            break

    assert high_quality_count > 0
    assert total_ref_length > 0


def test_cigar_based_coverage(bam_reader):
    """Test coverage calculation using CIGAR"""
    coverage = {}

    for record in bam_reader:
        if record.is_mapped and record.position is not None:
            pos = record.position

            # Calculate coverage from CIGAR
            for op in record.cigar:
                if op.consumes_reference():
                    for i in range(op.length):
                        coverage[pos] = coverage.get(pos, 0) + 1
                        pos += 1

        if len(coverage) >= 1000:
            break

    assert len(coverage) > 0
    assert all(count > 0 for count in coverage.values())


# ============================================================================
# Tag convenience methods (v1.4.0+)
# ============================================================================

def test_get_int(bam_reader):
    """Test get_int() convenience method"""
    for record in bam_reader:
        # Try to get edit distance (NM tag)
        nm = record.get_int("NM")

        # If NM exists, verify it matches get_tag()
        if nm is not None:
            tag = record.get_tag("NM")
            assert tag is not None
            assert tag.value.is_int()
            assert tag.value.as_int() == nm
            assert nm >= 0  # Edit distance should be non-negative
            break


def test_get_string(bam_reader):
    """Test get_string() convenience method"""
    for record in bam_reader:
        # Try to get MD string
        md = record.get_string("MD")

        # If MD exists, verify it matches get_tag()
        if md is not None:
            tag = record.get_tag("MD")
            assert tag is not None
            assert tag.value.is_string()
            assert tag.value.as_string() == md
            assert len(md) > 0  # MD string should not be empty
            break

    # Test with read group (RG) as well
    for record in bam_reader:
        rg = record.get_string("RG")
        if rg is not None:
            assert isinstance(rg, str)
            assert len(rg) > 0
            break


def test_get_int_invalid_name(bam_reader):
    """Test get_int() with invalid tag name"""
    for record in bam_reader:
        # Too short
        with pytest.raises(ValueError, match="exactly 2 characters"):
            record.get_int("N")

        # Too long
        with pytest.raises(ValueError, match="exactly 2 characters"):
            record.get_int("NMX")

        break


def test_get_string_invalid_name(bam_reader):
    """Test get_string() with invalid tag name"""
    for record in bam_reader:
        # Too short
        with pytest.raises(ValueError, match="exactly 2 characters"):
            record.get_string("M")

        # Too long
        with pytest.raises(ValueError, match="exactly 2 characters"):
            record.get_string("MDZ")

        break


def test_edit_distance(bam_reader):
    """Test edit_distance() convenience method"""
    found_nm = False

    for record in bam_reader:
        edit_dist = record.edit_distance()

        if edit_dist is not None:
            found_nm = True
            # Should match get_int("NM")
            assert edit_dist == record.get_int("NM")
            assert edit_dist >= 0

            # Should be <= read length (can't have more mismatches than bases)
            assert edit_dist <= len(record.sequence)

        if found_nm:
            break

    # Note: NM tag is optional, so finding it is best-effort


def test_alignment_score(bam_reader):
    """Test alignment_score() convenience method"""
    found_as = False

    for record in bam_reader:
        score = record.alignment_score()

        if score is not None:
            found_as = True
            # Should match get_int("AS")
            assert score == record.get_int("AS")
            # Alignment score can be negative for poor alignments
            assert isinstance(score, int)

        if found_as:
            break

    # Note: AS tag is optional, so finding it is best-effort


def test_read_group(bam_reader):
    """Test read_group() convenience method"""
    found_rg = False

    for record in bam_reader:
        rg = record.read_group()

        if rg is not None:
            found_rg = True
            # Should match get_string("RG")
            assert rg == record.get_string("RG")
            assert isinstance(rg, str)
            assert len(rg) > 0

        if found_rg:
            break

    # Note: RG tag is optional, so finding it is best-effort


def test_md_string(bam_reader):
    """Test md_string() convenience method"""
    found_md = False

    for record in bam_reader:
        md = record.md_string()

        if md is not None:
            found_md = True
            # Should match get_string("MD")
            assert md == record.get_string("MD")
            assert isinstance(md, str)
            assert len(md) > 0
            # MD string should contain digits and/or letters
            assert any(c.isalnum() for c in md)

        if found_md:
            break

    # Note: MD tag is optional, so finding it is best-effort


# ============================================================================
# Statistics functions (v1.4.0+)
# ============================================================================

def test_insert_size_distribution(bam_path):
    """Test insert_size_distribution() for paired-end QC"""
    dist = biometal.insert_size_distribution(bam_path)

    assert isinstance(dist, dict)

    # If we have insert sizes, verify properties
    if dist:
        # All keys should be positive integers
        for insert_size, count in dist.items():
            assert isinstance(insert_size, int)
            assert insert_size > 0
            assert isinstance(count, int)
            assert count > 0

        # Typical insert sizes are 100-1000bp for most libraries
        sizes = list(dist.keys())
        assert min(sizes) >= 0
        assert max(sizes) <= 100000  # Sanity check


def test_insert_size_distribution_filtered(bam_path):
    """Test insert_size_distribution() with reference filter"""
    # Test with reference 0
    dist_ref0 = biometal.insert_size_distribution(bam_path, reference_id=0)
    assert isinstance(dist_ref0, dict)

    # Test with all references
    dist_all = biometal.insert_size_distribution(bam_path)

    # Reference-filtered should be subset of all
    if dist_ref0 and dist_all:
        assert sum(dist_ref0.values()) <= sum(dist_all.values())


def test_edit_distance_stats(bam_path):
    """Test edit_distance_stats() for alignment quality"""
    stats = biometal.edit_distance_stats(bam_path)

    assert isinstance(stats, dict)
    assert 'total_records' in stats
    assert 'with_nm_tag' in stats
    assert 'mean' in stats
    assert 'median' in stats
    assert 'min' in stats
    assert 'max' in stats
    assert 'distribution' in stats

    assert stats['total_records'] >= 0
    assert stats['with_nm_tag'] >= 0
    assert stats['with_nm_tag'] <= stats['total_records']

    # If we have NM tags, verify statistics
    if stats['with_nm_tag'] > 0:
        assert stats['mean'] is not None
        assert stats['median'] is not None
        assert stats['min'] is not None
        assert stats['max'] is not None
        assert isinstance(stats['distribution'], dict)

        # Mean should be reasonable
        assert stats['mean'] >= 0
        assert stats['min'] <= stats['median'] <= stats['max']

        # Distribution should match count
        assert sum(stats['distribution'].values()) == stats['with_nm_tag']


def test_strand_bias(bam_path):
    """Test strand_bias() for variant calling QC"""
    # Test at position 1000 on reference 0
    bias = biometal.strand_bias(bam_path, reference_id=0, position=1000, window_size=100)

    assert isinstance(bias, dict)
    assert 'forward' in bias
    assert 'reverse' in bias
    assert 'total' in bias
    assert 'ratio' in bias

    assert bias['forward'] >= 0
    assert bias['reverse'] >= 0
    assert bias['total'] == bias['forward'] + bias['reverse']

    # Ratio should be 0-1
    assert 0.0 <= bias['ratio'] <= 1.0

    # If we have reads, ratio should match calculation
    if bias['total'] > 0:
        expected_ratio = bias['forward'] / bias['total']
        assert abs(bias['ratio'] - expected_ratio) < 0.001


def test_alignment_length_distribution(bam_path):
    """Test alignment_length_distribution() for RNA-seq QC"""
    dist = biometal.alignment_length_distribution(bam_path)

    assert isinstance(dist, dict)

    # If we have alignments, verify properties
    if dist:
        # All keys and values should be positive
        for length, count in dist.items():
            assert isinstance(length, int)
            assert length > 0
            assert isinstance(count, int)
            assert count > 0

        # Typical alignment lengths are 50-500bp
        lengths = list(dist.keys())
        assert min(lengths) > 0
        assert max(lengths) <= 100000  # Sanity check


def test_statistics_empty_filters(bam_path):
    """Test statistics functions with filters that return no records"""
    # Test with non-existent reference ID (likely doesn't exist)
    dist = biometal.insert_size_distribution(bam_path, reference_id=9999)
    assert isinstance(dist, dict)
    # May be empty or have records depending on BAM file

    stats = biometal.edit_distance_stats(bam_path, reference_id=9999)
    assert isinstance(stats, dict)
    assert stats['total_records'] >= 0


# ============================================================================
# Run tests
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
