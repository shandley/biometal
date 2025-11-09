#!/usr/bin/env python3
"""
Advanced BAM Filtering Examples with biometal

Demonstrates production-ready filtering techniques for BAM alignment files:
- Region queries (chr:start-end)
- Complex multi-criteria filters
- Common bioinformatics workflows
- Performance-optimized patterns

Performance:
- 4.54 million records/sec throughput
- 43.0 MiB/s compressed file processing
- 4× speedup via parallel BGZF decompression
- Constant ~5 MB memory footprint

Requirements:
- biometal-rs >= 1.3.0 (pip install biometal-rs)
- Python 3.9+

New in v1.3.0:
- CIGAR operations (alignment details)
- SAM writing (BAM → SAM conversion)
- CIGAR-aware coverage calculation

New in v1.4.0:
- Tag convenience methods (edit_distance, alignment_score, read_group, md_string)
- get_int/get_string accessors for direct tag value access
- Tag analysis functions (edit distance distribution, read group filtering)
"""

import biometal
from typing import Iterator, Dict, List, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict


# ============================================================================
# Region Query Helpers
# ============================================================================

@dataclass
class Region:
    """Genomic region for filtering."""
    reference_id: int
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive

    def contains(self, position: int) -> bool:
        """Check if position falls within region."""
        return self.start <= position < self.end

    def overlaps(self, start: int, end: int) -> bool:
        """Check if interval [start, end) overlaps this region."""
        return not (end <= self.start or start >= self.end)

    @classmethod
    def from_string(cls, region_str: str, ref_map: Dict[str, int]) -> 'Region':
        """
        Parse region string like 'chr1:1000-2000' or '0:1000-2000'.

        Args:
            region_str: Region string (e.g., 'chr1:1000-2000' or '0:1000-2000')
            ref_map: Dictionary mapping reference names to IDs

        Returns:
            Region object

        Examples:
            >>> ref_map = {'chr1': 0, 'chr2': 1}
            >>> region = Region.from_string('chr1:1000-2000', ref_map)
            >>> region.reference_id
            0
        """
        parts = region_str.split(':')
        if len(parts) != 2:
            raise ValueError(f"Invalid region format: {region_str}")

        ref_name = parts[0]

        # Handle numeric reference IDs
        if ref_name.isdigit():
            reference_id = int(ref_name)
        elif ref_name in ref_map:
            reference_id = ref_map[ref_name]
        else:
            raise ValueError(f"Unknown reference: {ref_name}")

        # Parse coordinates
        coords = parts[1].split('-')
        if len(coords) != 2:
            raise ValueError(f"Invalid coordinates: {parts[1]}")

        start = int(coords[0])
        end = int(coords[1])

        if start >= end:
            raise ValueError(f"Invalid range: start ({start}) >= end ({end})")

        return cls(reference_id, start, end)


def filter_by_region(
    bam_path: str,
    region: Region,
    require_primary: bool = True,
    min_mapq: int = 0
) -> Iterator[biometal.BamRecord]:
    """
    Filter BAM records by genomic region.

    Args:
        bam_path: Path to BAM file
        region: Region to query
        require_primary: Only return primary alignments
        min_mapq: Minimum mapping quality

    Yields:
        BamRecord objects within the specified region

    Examples:
        >>> region = Region(reference_id=0, start=1000, end=2000)
        >>> for record in filter_by_region("alignments.bam", region, min_mapq=30):
        ...     print(f"{record.name}: {record.position}")
    """
    bam = biometal.BamReader.from_path(bam_path)

    for record in bam:
        # Check reference matches
        if record.reference_id != region.reference_id:
            continue

        # Check if unmapped
        if record.position is None:
            continue

        # Check if overlaps region
        if not region.contains(record.position):
            continue

        # Apply filters
        if require_primary and not record.is_primary:
            continue

        if record.mapq is not None and record.mapq < min_mapq:
            continue

        yield record


# ============================================================================
# Complex Filters
# ============================================================================

@dataclass
class FilterCriteria:
    """Comprehensive filtering criteria for BAM records."""

    # Alignment quality
    min_mapq: Optional[int] = None
    max_mapq: Optional[int] = None

    # Alignment type
    require_mapped: bool = False
    require_primary: bool = False
    require_proper_pair: bool = False

    # Strand
    require_forward: bool = False
    require_reverse: bool = False

    # Read properties
    min_length: Optional[int] = None
    max_length: Optional[int] = None

    # Reference
    allowed_references: Optional[List[int]] = None
    excluded_references: Optional[List[int]] = None

    # Region
    region: Optional[Region] = None

    def passes(self, record: biometal.BamRecord) -> bool:
        """
        Check if record passes all criteria.

        Args:
            record: BamRecord to evaluate

        Returns:
            True if record passes all filters
        """
        # MAPQ filters
        if self.min_mapq is not None:
            if record.mapq is None or record.mapq < self.min_mapq:
                return False

        if self.max_mapq is not None:
            if record.mapq is None or record.mapq > self.max_mapq:
                return False

        # Alignment type filters
        if self.require_mapped and not record.is_mapped:
            return False

        if self.require_primary and not record.is_primary:
            return False

        if self.require_proper_pair and not record.is_paired:
            return False

        # Strand filters
        if self.require_forward and record.is_reverse:
            return False

        if self.require_reverse and not record.is_reverse:
            return False

        # Length filters
        seq_len = len(record.sequence)
        if self.min_length is not None and seq_len < self.min_length:
            return False

        if self.max_length is not None and seq_len > self.max_length:
            return False

        # Reference filters
        if self.allowed_references is not None:
            if record.reference_id not in self.allowed_references:
                return False

        if self.excluded_references is not None:
            if record.reference_id in self.excluded_references:
                return False

        # Region filter
        if self.region is not None:
            if record.reference_id != self.region.reference_id:
                return False
            if record.position is None:
                return False
            if not self.region.contains(record.position):
                return False

        return True


def filter_bam(
    bam_path: str,
    criteria: FilterCriteria,
    limit: Optional[int] = None
) -> Iterator[biometal.BamRecord]:
    """
    Filter BAM records with complex criteria.

    Args:
        bam_path: Path to BAM file
        criteria: FilterCriteria object with filter settings
        limit: Optional maximum number of records to return

    Yields:
        BamRecord objects that pass all criteria

    Examples:
        >>> # High-quality, primary, forward-strand alignments
        >>> criteria = FilterCriteria(
        ...     min_mapq=30,
        ...     require_primary=True,
        ...     require_forward=True,
        ...     min_length=50
        ... )
        >>> for record in filter_bam("alignments.bam", criteria, limit=1000):
        ...     process(record)
    """
    bam = biometal.BamReader.from_path(bam_path)
    count = 0

    for record in bam:
        if criteria.passes(record):
            yield record
            count += 1
            if limit is not None and count >= limit:
                break


# ============================================================================
# Common Workflows
# ============================================================================

def variant_calling_filter(bam_path: str) -> Iterator[biometal.BamRecord]:
    """
    Filter for variant calling: high quality, primary, mapped.

    Standard filters:
    - MAPQ >= 30 (high confidence)
    - Primary alignments only
    - Mapped to reference

    Args:
        bam_path: Path to BAM file

    Yields:
        High-quality alignments suitable for variant calling
    """
    criteria = FilterCriteria(
        min_mapq=30,
        require_mapped=True,
        require_primary=True
    )
    return filter_bam(bam_path, criteria)


def strand_specific_rna_seq(
    bam_path: str,
    forward_only: bool = True
) -> Iterator[biometal.BamRecord]:
    """
    Filter for strand-specific RNA-seq analysis.

    Args:
        bam_path: Path to BAM file
        forward_only: If True, return only forward strand; else reverse

    Yields:
        Strand-specific alignments
    """
    criteria = FilterCriteria(
        require_mapped=True,
        require_primary=True,
        require_forward=forward_only,
        require_reverse=not forward_only
    )
    return filter_bam(bam_path, criteria)


def paired_end_properly_paired(bam_path: str) -> Iterator[biometal.BamRecord]:
    """
    Filter for properly paired reads (quality control).

    Args:
        bam_path: Path to BAM file

    Yields:
        Properly paired alignments
    """
    criteria = FilterCriteria(
        require_proper_pair=True,
        require_primary=True
    )
    return filter_bam(bam_path, criteria)


def extract_region(
    bam_path: str,
    reference_id: int,
    start: int,
    end: int,
    min_mapq: int = 20
) -> List[biometal.BamRecord]:
    """
    Extract all alignments from a specific genomic region.

    Args:
        bam_path: Path to BAM file
        reference_id: Reference sequence ID
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        min_mapq: Minimum mapping quality

    Returns:
        List of alignments in the region

    Examples:
        >>> # Extract chr1:1000-2000
        >>> alignments = extract_region("alignments.bam", 0, 1000, 2000)
        >>> print(f"Found {len(alignments)} alignments in region")
    """
    region = Region(reference_id, start, end)
    return list(filter_by_region(bam_path, region, min_mapq=min_mapq))


# ============================================================================
# Coverage Analysis
# ============================================================================

def calculate_coverage(
    bam_path: str,
    reference_id: int,
    start: int,
    end: int,
    min_mapq: int = 20
) -> Dict[int, int]:
    """
    Calculate per-base coverage for a region using CIGAR operations.

    Uses CIGAR strings to accurately determine which reference positions
    are covered by each alignment (accounting for insertions, deletions, etc.).

    Args:
        bam_path: Path to BAM file
        reference_id: Reference sequence ID
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        min_mapq: Minimum mapping quality

    Returns:
        Dictionary mapping position -> coverage depth

    Examples:
        >>> coverage = calculate_coverage("alignments.bam", 0, 1000, 2000)
        >>> mean_cov = sum(coverage.values()) / len(coverage)
        >>> print(f"Mean coverage: {mean_cov:.1f}×")
    """
    coverage = defaultdict(int)
    region = Region(reference_id, start, end)

    for record in filter_by_region(bam_path, region, min_mapq=min_mapq):
        # Use CIGAR operations for accurate coverage calculation
        pos = record.position

        for op in record.cigar:
            # Only count operations that consume reference sequence
            if op.consumes_reference():
                for i in range(op.length):
                    ref_pos = pos + i
                    if start <= ref_pos < end:
                        coverage[ref_pos] += 1
                pos += op.length
            # Insertions and clips don't advance reference position

    return dict(coverage)


def find_high_coverage_regions(
    bam_path: str,
    reference_id: int,
    min_coverage: int = 20,
    min_length: int = 100
) -> List[Tuple[int, int, float]]:
    """
    Identify high-coverage regions (e.g., for CNV detection).

    Args:
        bam_path: Path to BAM file
        reference_id: Reference sequence ID
        min_coverage: Minimum coverage depth
        min_length: Minimum region length

    Returns:
        List of (start, end, mean_coverage) tuples

    Examples:
        >>> regions = find_high_coverage_regions("alignments.bam", 0)
        >>> for start, end, cov in regions:
        ...     print(f"Region {start}-{end}: {cov:.1f}× coverage")
    """
    # This is a simplified example - production would use sliding windows
    high_cov_regions = []

    # Count alignments starting at each position
    position_counts = defaultdict(int)

    criteria = FilterCriteria(
        require_mapped=True,
        require_primary=True,
        allowed_references=[reference_id]
    )

    for record in filter_bam(bam_path, criteria):
        if record.position is not None:
            position_counts[record.position] += 1

    # Find contiguous high-coverage regions
    in_region = False
    region_start = None
    region_positions = []

    for pos in sorted(position_counts.keys()):
        cov = position_counts[pos]

        if cov >= min_coverage:
            if not in_region:
                in_region = True
                region_start = pos
            region_positions.append((pos, cov))
        else:
            if in_region and len(region_positions) >= min_length:
                mean_cov = sum(c for _, c in region_positions) / len(region_positions)
                high_cov_regions.append((region_start, pos, mean_cov))
            in_region = False
            region_positions = []

    return high_cov_regions


# ============================================================================
# CIGAR Operations (v1.3.0)
# ============================================================================

def analyze_cigar_operations(bam_path: str, limit: int = 1000) -> Dict[str, int]:
    """
    Analyze CIGAR operation distribution in alignments.

    CIGAR operations describe how reads align to the reference:
    - M: Match/mismatch
    - I: Insertion
    - D: Deletion
    - N: Skipped reference (introns in RNA-seq)
    - S: Soft clipping
    - H: Hard clipping
    - P: Padding
    - =: Sequence match
    - X: Sequence mismatch

    Args:
        bam_path: Path to BAM file
        limit: Number of records to analyze

    Returns:
        Dictionary with operation counts

    Examples:
        >>> ops = analyze_cigar_operations("alignments.bam")
        >>> print(f"Matches: {ops['M']}, Insertions: {ops['I']}, Deletions: {ops['D']}")
    """
    op_counts = defaultdict(int)
    bam = biometal.BamReader.from_path(bam_path)

    count = 0
    for record in bam:
        for op in record.cigar:
            op_counts[op.op_char] += op.length

        count += 1
        if count >= limit:
            break

    return dict(op_counts)


def find_indels(
    bam_path: str,
    reference_id: int,
    min_indel_length: int = 3
) -> List[Tuple[str, int, int, str]]:
    """
    Find insertions and deletions using CIGAR operations.

    Args:
        bam_path: Path to BAM file
        reference_id: Reference sequence ID
        min_indel_length: Minimum indel length to report

    Returns:
        List of (type, position, length, read_name) tuples

    Examples:
        >>> indels = find_indels("alignments.bam", 0, min_indel_length=5)
        >>> for indel_type, pos, length, name in indels:
        ...     print(f"{indel_type} at {pos}: {length}bp in {name}")
    """
    indels = []

    criteria = FilterCriteria(
        require_mapped=True,
        require_primary=True,
        allowed_references=[reference_id]
    )

    for record in filter_bam(bam_path, criteria):
        pos = record.position

        for op in record.cigar:
            if op.is_insertion() and op.length >= min_indel_length:
                indels.append(("INS", pos, op.length, record.name))
            elif op.is_deletion() and op.length >= min_indel_length:
                indels.append(("DEL", pos, op.length, record.name))

            # Advance position for reference-consuming operations
            if op.consumes_reference():
                pos += op.length

    return indels


def calculate_alignment_metrics(record: biometal.BamRecord) -> Dict[str, int]:
    """
    Calculate alignment quality metrics from CIGAR operations.

    Args:
        record: BamRecord to analyze

    Returns:
        Dictionary with alignment metrics

    Examples:
        >>> reader = biometal.BamReader.from_path("alignments.bam")
        >>> for record in reader:
        ...     metrics = calculate_alignment_metrics(record)
        ...     print(f"Matches: {metrics['matches']}, Indels: {metrics['indels']}")
        ...     break
    """
    metrics = {
        'matches': 0,
        'insertions': 0,
        'deletions': 0,
        'soft_clips': 0,
        'hard_clips': 0,
        'indels': 0,
        'query_length': record.query_length(),
        'reference_length': record.reference_length(),
    }

    for op in record.cigar:
        if op.is_match() or op.is_seq_match() or op.is_seq_mismatch():
            metrics['matches'] += op.length
        elif op.is_insertion():
            metrics['insertions'] += op.length
            metrics['indels'] += 1
        elif op.is_deletion():
            metrics['deletions'] += op.length
            metrics['indels'] += 1
        elif op.is_soft_clip():
            metrics['soft_clips'] += op.length
        elif op.is_hard_clip():
            metrics['hard_clips'] += op.length

    return metrics


# ============================================================================
# SAM Writing (v1.3.0)
# ============================================================================

def convert_bam_to_sam(
    bam_path: str,
    sam_path: str,
    criteria: Optional[FilterCriteria] = None,
    limit: Optional[int] = None
) -> int:
    """
    Convert BAM to SAM format with optional filtering.

    Args:
        bam_path: Input BAM file path
        sam_path: Output SAM file path
        criteria: Optional filter criteria
        limit: Optional limit on records to convert

    Returns:
        Number of records written

    Examples:
        >>> # Convert entire file
        >>> count = convert_bam_to_sam("input.bam", "output.sam")
        >>> print(f"Converted {count} records")

        >>> # Convert only high-quality alignments
        >>> criteria = FilterCriteria(min_mapq=30, require_primary=True)
        >>> count = convert_bam_to_sam("input.bam", "filtered.sam", criteria)
    """
    reader = biometal.BamReader.from_path(bam_path)
    writer = biometal.SamWriter.create(sam_path)

    # Write header
    writer.write_header(reader.header)

    # Write records
    count = 0
    for record in reader:
        # Apply filter if specified
        if criteria is not None and not criteria.passes(record):
            continue

        writer.write_record(record)
        count += 1

        if limit is not None and count >= limit:
            break

    writer.close()
    return count


def extract_region_to_sam(
    bam_path: str,
    sam_path: str,
    reference_id: int,
    start: int,
    end: int,
    min_mapq: int = 20
) -> int:
    """
    Extract a genomic region to SAM format.

    Args:
        bam_path: Input BAM file path
        sam_path: Output SAM file path
        reference_id: Reference sequence ID
        start: Start position (0-based, inclusive)
        end: End position (0-based, exclusive)
        min_mapq: Minimum mapping quality

    Returns:
        Number of records written

    Examples:
        >>> # Extract chr1:1000-2000 to SAM
        >>> count = extract_region_to_sam("alignments.bam", "region.sam", 0, 1000, 2000)
        >>> print(f"Extracted {count} alignments to region.sam")
    """
    reader = biometal.BamReader.from_path(bam_path)
    writer = biometal.SamWriter.create(sam_path)

    # Write header
    writer.write_header(reader.header)

    # Write filtered records
    region = Region(reference_id, start, end)
    count = 0

    for record in filter_by_region(bam_path, region, min_mapq=min_mapq):
        writer.write_record(record)
        count += 1

    writer.close()
    return count


# ============================================================================
# Statistics and Summary
# ============================================================================

def count_by_flag(
    bam_path: str,
    limit: Optional[int] = None
) -> Dict[str, int]:
    """
    Count alignments by various flag properties.

    Args:
        bam_path: Path to BAM file
        limit: Optional limit on records to process

    Returns:
        Dictionary with counts for each property

    Examples:
        >>> stats = count_by_flag("alignments.bam", limit=10000)
        >>> print(f"Mapped: {stats['mapped']} / {stats['total']}")
    """
    counts = {
        'total': 0,
        'mapped': 0,
        'unmapped': 0,
        'primary': 0,
        'secondary': 0,
        'paired': 0,
        'forward': 0,
        'reverse': 0,
        'first': 0,
        'second': 0,
    }

    bam = biometal.BamReader.from_path(bam_path)

    for record in bam:
        counts['total'] += 1

        if record.is_mapped:
            counts['mapped'] += 1
        else:
            counts['unmapped'] += 1

        if record.is_primary:
            counts['primary'] += 1
        else:
            counts['secondary'] += 1

        if record.is_paired:
            counts['paired'] += 1

        if record.is_reverse:
            counts['reverse'] += 1
        else:
            counts['forward'] += 1

        if record.is_first:
            counts['first'] += 1

        if record.is_second:
            counts['second'] += 1

        if limit is not None and counts['total'] >= limit:
            break

    return counts


def mapq_distribution(
    bam_path: str,
    limit: Optional[int] = None
) -> Dict[int, int]:
    """
    Calculate MAPQ score distribution.

    Args:
        bam_path: Path to BAM file
        limit: Optional limit on records to process

    Returns:
        Dictionary mapping MAPQ -> count

    Examples:
        >>> dist = mapq_distribution("alignments.bam", limit=10000)
        >>> for mapq in sorted(dist.keys()):
        ...     print(f"MAPQ {mapq}: {dist[mapq]} alignments")
    """
    distribution = defaultdict(int)

    bam = biometal.BamReader.from_path(bam_path)
    count = 0

    for record in bam:
        if record.mapq is not None:
            distribution[record.mapq] += 1
        count += 1

        if limit is not None and count >= limit:
            break

    return dict(distribution)


# ============================================================================
# Tag Analysis Functions (v1.4.0+)
# ============================================================================

def analyze_edit_distance(
    bam_path: str,
    reference_id: Optional[int] = None,
    limit: Optional[int] = None
) -> Dict[str, any]:
    """
    Analyze edit distance (NM tag) distribution for alignment quality assessment.

    Args:
        bam_path: Path to BAM file
        reference_id: Optional reference ID to filter by
        limit: Optional limit on number of records to analyze

    Returns:
        Dictionary containing:
        - distribution: Dict[int, int] mapping edit distance → count
        - mean: Average edit distance
        - median: Median edit distance
        - max: Maximum edit distance observed
        - with_nm_tag: Number of records with NM tag
        - total_records: Total records analyzed

    Example:
        >>> stats = analyze_edit_distance("alignments.bam", limit=100000)
        >>> print(f"Mean edit distance: {stats['mean']:.2f}")
        >>> print(f"Records with NM tag: {stats['with_nm_tag']:,} / {stats['total_records']:,}")
    """
    reader = biometal.BamReader.from_path(bam_path)
    distribution = defaultdict(int)
    edit_distances = []
    count = 0
    with_nm = 0

    for record in reader:
        # Filter by reference if specified
        if reference_id is not None and record.reference_id != reference_id:
            continue

        # Get edit distance using convenience method
        edit_dist = record.edit_distance()

        if edit_dist is not None:
            distribution[edit_dist] += 1
            edit_distances.append(edit_dist)
            with_nm += 1

        count += 1

        if limit is not None and count >= limit:
            break

    # Calculate statistics
    stats = {
        "distribution": dict(distribution),
        "with_nm_tag": with_nm,
        "total_records": count,
    }

    if edit_distances:
        edit_distances.sort()
        stats["mean"] = sum(edit_distances) / len(edit_distances)
        stats["median"] = edit_distances[len(edit_distances) // 2]
        stats["max"] = max(edit_distances)
        stats["min"] = min(edit_distances)
    else:
        stats["mean"] = None
        stats["median"] = None
        stats["max"] = None
        stats["min"] = None

    return stats


def filter_by_read_group(
    bam_path: str,
    read_groups: List[str],
    output_path: Optional[str] = None
) -> Iterator[biometal.BamRecord]:
    """
    Filter BAM records by read group (RG tag).

    Useful for:
    - Sample-specific analysis in multiplexed sequencing
    - Lane-specific QC
    - Flowcell-specific filtering

    Args:
        bam_path: Path to BAM file
        read_groups: List of read group IDs to include
        output_path: Optional output SAM path for filtered records

    Yields:
        BamRecord: Records matching specified read groups

    Example:
        >>> # Filter records from specific read groups
        >>> for record in filter_by_read_group("alignments.bam", ["RG001", "RG002"]):
        ...     print(f"{record.name}: RG={record.read_group()}")
        >>>
        >>> # Export filtered records to SAM
        >>> count = 0
        >>> for record in filter_by_read_group("alignments.bam", ["RG001"], "filtered.sam"):
        ...     count += 1
        >>> print(f"Exported {count:,} records")
    """
    reader = biometal.BamReader.from_path(bam_path)
    writer = None

    if output_path:
        writer = biometal.SamWriter.create(output_path)
        writer.write_header(reader.header)

    read_group_set = set(read_groups)

    for record in reader:
        rg = record.read_group()

        if rg in read_group_set:
            if writer:
                writer.write_record(record)
            yield record

    if writer:
        writer.close()


def analyze_tag_coverage(
    bam_path: str,
    tag_names: List[str],
    limit: Optional[int] = None
) -> Dict[str, Dict[str, any]]:
    """
    Analyze coverage of optional tags across BAM records.

    Useful for:
    - Understanding which tags are present in your BAM files
    - Quality control for alignment pipelines
    - Debugging missing tag issues

    Args:
        bam_path: Path to BAM file
        tag_names: List of 2-character tag names to analyze (e.g., ["NM", "AS", "RG"])
        limit: Optional limit on number of records to analyze

    Returns:
        Dictionary mapping tag name → statistics:
        - present: Number of records with this tag
        - missing: Number of records without this tag
        - coverage: Percentage of records with this tag

    Example:
        >>> stats = analyze_tag_coverage("alignments.bam", ["NM", "AS", "RG", "MD"])
        >>> for tag, info in stats.items():
        ...     print(f"{tag}: {info['coverage']:.1f}% coverage ({info['present']:,} / {info['present'] + info['missing']:,})")
    """
    reader = biometal.BamReader.from_path(bam_path)
    stats = {tag: {"present": 0, "missing": 0} for tag in tag_names}
    count = 0

    for record in reader:
        for tag in tag_names:
            if record.has_tag(tag):
                stats[tag]["present"] += 1
            else:
                stats[tag]["missing"] += 1

        count += 1

        if limit is not None and count >= limit:
            break

    # Calculate coverage percentages
    for tag in tag_names:
        total = stats[tag]["present"] + stats[tag]["missing"]
        stats[tag]["coverage"] = (stats[tag]["present"] / total * 100) if total > 0 else 0.0

    return stats


# ============================================================================
# Example Usage
# ============================================================================

def main():
    """Demonstrate advanced filtering techniques."""
    import sys

    if len(sys.argv) < 2:
        print("Usage: python bam_advanced_filtering.py <bam_file>")
        print("\nExamples:")
        print("  python bam_advanced_filtering.py alignments.bam")
        sys.exit(1)

    bam_path = sys.argv[1]

    print(f"📊 Advanced BAM Filtering Examples")
    print(f"{'=' * 60}\n")
    print(f"File: {bam_path}\n")

    # Example 1: Variant calling filter
    print("1. Variant Calling Filter (MAPQ ≥ 30, primary, mapped):")
    count = sum(1 for _ in variant_calling_filter(bam_path))
    print(f"   Found {count:,} high-quality alignments\n")

    # Example 2: Region query
    print("2. Region Query (chr0:0-1000):")
    region = Region(reference_id=0, start=0, end=1000)
    region_records = list(filter_by_region(bam_path, region, min_mapq=20))
    print(f"   Found {len(region_records):,} alignments in region\n")

    # Example 3: Complex filter
    print("3. Complex Filter (high-quality, forward strand, 50-150 bp):")
    criteria = FilterCriteria(
        min_mapq=30,
        require_primary=True,
        require_forward=True,
        min_length=50,
        max_length=150
    )
    complex_count = sum(1 for _ in filter_bam(bam_path, criteria, limit=10000))
    print(f"   Found {complex_count:,} matching alignments (first 10K)\n")

    # Example 4: Flag statistics
    print("4. Flag Statistics (first 10K records):")
    stats = count_by_flag(bam_path, limit=10000)
    print(f"   Total: {stats['total']:,}")
    print(f"   Mapped: {stats['mapped']:,} ({100*stats['mapped']/stats['total']:.1f}%)")
    print(f"   Primary: {stats['primary']:,} ({100*stats['primary']/stats['total']:.1f}%)")
    print(f"   Strand: {stats['forward']:,} forward, {stats['reverse']:,} reverse\n")

    # Example 5: MAPQ distribution
    print("5. MAPQ Distribution (first 10K records):")
    dist = mapq_distribution(bam_path, limit=10000)
    for mapq in sorted(dist.keys())[:10]:  # Show first 10
        count = dist[mapq]
        print(f"   MAPQ {mapq}: {count:,} alignments")
    print()

    # Example 6: CIGAR operation analysis (v1.3.0)
    print("6. CIGAR Operation Analysis (first 1K records):")
    ops = analyze_cigar_operations(bam_path, limit=1000)
    for op_char in sorted(ops.keys()):
        count = ops[op_char]
        print(f"   {op_char}: {count:,} bases")
    print()

    # Example 7: Indel detection (v1.3.0)
    print("7. Indel Detection (≥5bp, reference 0, first 1K records):")
    indels = find_indels(bam_path, reference_id=0, min_indel_length=5)[:10]
    if indels:
        for indel_type, pos, length, name in indels[:5]:  # Show first 5
            print(f"   {indel_type} at {pos}: {length}bp ({name})")
    else:
        print("   No indels ≥5bp found")
    print()

    # Example 8: SAM writing (v1.3.0)
    print("8. SAM Writing (extract region to SAM):")
    sam_path = "region_sample.sam"
    sam_count = extract_region_to_sam(bam_path, sam_path, 0, 0, 1000, min_mapq=20)
    print(f"   Wrote {sam_count:,} alignments to {sam_path}")
    print()

    # Example 9: Edit distance analysis (v1.4.0)
    print("9. Edit Distance Analysis (first 10K records):")
    ed_stats = analyze_edit_distance(bam_path, limit=10000)
    if ed_stats['mean'] is not None:
        print(f"   Records with NM tag: {ed_stats['with_nm_tag']:,} / {ed_stats['total_records']:,}")
        print(f"   Mean edit distance: {ed_stats['mean']:.2f}")
        print(f"   Median edit distance: {ed_stats['median']}")
        print(f"   Range: {ed_stats['min']} - {ed_stats['max']}")
    else:
        print("   No NM tags found")
    print()

    # Example 10: Tag coverage analysis (v1.4.0)
    print("10. Tag Coverage Analysis (first 10K records):")
    tag_stats = analyze_tag_coverage(bam_path, ["NM", "AS", "RG", "MD"], limit=10000)
    for tag, info in sorted(tag_stats.items()):
        print(f"    {tag}: {info['coverage']:.1f}% ({info['present']:,} / {info['present'] + info['missing']:,})")
    print()

    # Example 11: Insert size distribution (v1.4.0)
    print("11. Insert Size Distribution (paired-end QC, first 1K records):")
    insert_dist = biometal.insert_size_distribution(bam_path)
    if insert_dist:
        sizes = list(insert_dist.keys())
        total_pairs = sum(insert_dist.values())
        mean_insert = sum(s * insert_dist[s] for s in sizes) / total_pairs
        print(f"    Total pairs: {total_pairs:,}")
        print(f"    Mean insert size: {mean_insert:.0f}bp")
        print(f"    Range: {min(sizes)}-{max(sizes)}bp")
        # Show most common sizes
        top_sizes = sorted(insert_dist.items(), key=lambda x: x[1], reverse=True)[:5]
        for size, count in top_sizes:
            print(f"    {size}bp: {count:,} pairs")
    else:
        print("    No properly paired reads found")
    print()

    # Example 12: Edit distance statistics (v1.4.0)
    print("12. Edit Distance Statistics (alignment quality):")
    ed_stats = biometal.edit_distance_stats(bam_path)
    print(f"    Total records: {ed_stats['total_records']:,}")
    print(f"    With NM tag: {ed_stats['with_nm_tag']:,} ({ed_stats['with_nm_tag']/ed_stats['total_records']*100:.1f}%)")
    if ed_stats['mean'] is not None:
        print(f"    Mean edit distance: {ed_stats['mean']:.2f}")
        print(f"    Median: {ed_stats['median']}")
        print(f"    Range: {ed_stats['min']}-{ed_stats['max']}")
    print()

    # Example 13: Strand bias analysis (v1.4.0)
    print("13. Strand Bias Analysis (variant QC at chr0:1000):")
    bias = biometal.strand_bias(bam_path, reference_id=0, position=1000, window_size=100)
    print(f"    Forward strand: {bias['forward']:,} reads")
    print(f"    Reverse strand: {bias['reverse']:,} reads")
    print(f"    Total: {bias['total']:,} reads")
    if bias['total'] > 0:
        print(f"    Strand ratio: {bias['ratio']:.3f} (0.5 = no bias)")
        if abs(bias['ratio'] - 0.5) > 0.2:
            print(f"    ⚠️  Warning: Potential strand bias detected")
    print()

    # Example 14: Alignment length distribution (v1.4.0)
    print("14. Alignment Length Distribution (RNA-seq QC):")
    len_dist = biometal.alignment_length_distribution(bam_path)
    if len_dist:
        lengths = list(len_dist.keys())
        total_alns = sum(len_dist.values())
        mean_len = sum(l * len_dist[l] for l in lengths) / total_alns
        print(f"    Total alignments: {total_alns:,}")
        print(f"    Mean length: {mean_len:.0f}bp")
        print(f"    Range: {min(lengths)}-{max(lengths)}bp")
        # Show most common lengths
        top_lens = sorted(len_dist.items(), key=lambda x: x[1], reverse=True)[:5]
        for length, count in top_lens:
            print(f"    {length}bp: {count:,} alignments")
    else:
        print("    No mapped reads found")
    print()

    print(f"{'=' * 60}")
    print("✅ Analysis complete")
    print(f"\n💡 New in v1.3.0:")
    print(f"   - CIGAR operations for detailed alignment analysis")
    print(f"   - SAM writing for format conversion and region extraction")
    print(f"   - CIGAR-aware coverage calculation")
    print(f"\n💡 New in v1.4.0:")
    print(f"   - Tag convenience methods (edit_distance, alignment_score, etc.)")
    print(f"   - Insert size distribution for paired-end QC")
    print(f"   - Edit distance statistics for alignment quality")
    print(f"   - Strand bias analysis for variant calling")
    print(f"   - Alignment length distribution for RNA-seq")
    print(f"\n💡 Tip: Modify criteria objects to customize filters")
    print(f"   Performance: 4.54M records/sec, constant ~5 MB memory")


if __name__ == '__main__':
    main()
