# BAM/SAM API Reference

Complete reference for biometal's Python BAM/SAM parser with streaming architecture and production-grade performance.

**Version**: 1.4.0+
**Performance**: 4.54M records/sec, 43.0 MiB/s throughput
**Memory**: Constant ~5 MB footprint (streams terabyte-scale files)

---

## Quick Start

```python
import biometal

# Open BAM file (streaming, constant memory)
reader = biometal.BamReader.from_path("alignments.bam")

# Iterate through records
for record in reader:
    # Access basic fields
    print(f"{record.name}: chr{record.reference_id}:{record.position}")
    print(f"  MAPQ: {record.mapq}, Flags: {record.flags}")

    # Check alignment quality
    if record.is_mapped and record.is_primary and record.mapq >= 30:
        # Analyze CIGAR operations
        for op in record.cigar:
            if op.is_insertion() and op.length >= 5:
                print(f"  Large insertion: {op.length}bp")

        # Access tags
        edit_dist = record.edit_distance()
        if edit_dist is not None:
            print(f"  Mismatches: {edit_dist}")
```

---

## Table of Contents

1. [BamReader Class](#bamreader-class)
2. [BamRecord Class](#bamrecord-class)
3. [BamHeader Class](#bamheader-class)
4. [CigarOp Class](#cigarop-class)
5. [Tag Classes](#tag-classes)
6. [SamWriter Class](#samwriter-class)
7. [Statistics Functions](#statistics-functions)
8. [Error Handling](#error-handling)
9. [Performance Tips](#performance-tips)

---

## BamReader Class

### Constructor

#### `BamReader.from_path(path: str) -> BamReader`

Open a BAM file for streaming reading.

**Parameters**:
- `path` (str): Path to BAM file (local or HTTP URL)

**Returns**: BamReader instance

**Raises**:
- `IOError`: If file cannot be opened or is not valid BAM format

**Example**:
```python
# Local file
reader = biometal.BamReader.from_path("data.bam")

# HTTP streaming (no download needed!)
reader = biometal.BamReader.from_path("https://example.com/data.bam")
```

---

### Properties

#### `reader.header -> BamHeader`

Access BAM file header containing reference sequences and metadata.

**Example**:
```python
header = reader.header
print(f"References: {header.reference_count}")
for i, ref in enumerate(header.references):
    print(f"  {i}: {ref.name} ({ref.length:,}bp)")
```

---

### Methods

#### `reader.__iter__() -> Iterator[BamRecord]`

Iterate through all records in the BAM file (streaming, constant memory).

**Example**:
```python
for record in reader:
    # Process one record at a time
    # Memory usage: constant ~5 MB
    process(record)
```

#### `reader.query(reference_name: str, start: int = None, end: int = None) -> BamRegionIter`

Query records in a genomic region.

**Parameters**:
- `reference_name` (str): Reference sequence name (e.g., "chr1")
- `start` (int, optional): Start position (0-based, inclusive)
- `end` (int, optional): End position (0-based, exclusive)

**Returns**: Iterator yielding BamRecord objects in the region

**Note**: Currently performs full file scan. BAI/CSI index support planned for v1.5.0.

**Example**:
```python
# Query specific region
for record in reader.query("chr1", start=1000, end=2000):
    print(f"{record.name} at {record.position}")

# Query entire chromosome
for record in reader.query("chr1"):
    process(record)
```

---

## BamRecord Class

### Core Properties

#### `record.name -> str`
Read name/identifier

#### `record.reference_id -> int | None`
Reference sequence ID (0-based), None for unmapped

#### `record.position -> int | None`
0-based leftmost mapping position, None for unmapped

#### `record.mapq -> int | None`
Mapping quality (0-255), None if unavailable

#### `record.flags -> int`
SAM flags (bitwise), use flag properties below

#### `record.mate_reference_id -> int | None`
Mate's reference sequence ID

#### `record.mate_position -> int | None`
Mate's mapping position

#### `record.template_length -> int`
Template length (insert size for paired-end)

#### `record.sequence -> bytes`
Read sequence (raw bytes)

#### `record.quality -> bytes`
Phred quality scores (raw bytes)

**Example**:
```python
print(f"Name: {record.name}")
print(f"Position: chr{record.reference_id}:{record.position}")
print(f"MAPQ: {record.mapq}")
print(f"Insert size: {record.template_length}bp")
```

---

### Sequence Methods

#### `record.sequence_str() -> str`
Get sequence as string (decodes from 4-bit encoding)

#### `record.quality_str() -> str`
Get quality scores as string

**Example**:
```python
seq = record.sequence_str()
qual = record.quality_str()
print(f"Sequence: {seq}")
print(f"Quality:  {qual}")
```

---

### Flag Properties

Boolean properties for SAM flags:

#### `record.is_paired -> bool`
Read is paired-end

#### `record.is_proper_pair -> bool`
Read and mate are properly paired

#### `record.is_mapped -> bool`
Read is mapped

#### `record.is_unmapped -> bool`
Read is unmapped (convenience: `not is_mapped`)

#### `record.is_mate_mapped -> bool`
Mate is mapped

#### `record.is_reverse -> bool`
Read is on reverse strand

#### `record.is_forward -> bool`
Read is on forward strand (convenience: `not is_reverse`)

#### `record.is_mate_reverse -> bool`
Mate is on reverse strand

#### `record.is_first -> bool`
Read is first in pair

#### `record.is_second -> bool`
Read is second in pair

#### `record.is_primary -> bool`
Primary alignment (not secondary/supplementary)

#### `record.is_secondary -> bool`
Secondary alignment

#### `record.is_supplementary -> bool`
Supplementary alignment

#### `record.is_qc_fail -> bool`
Failed quality checks

#### `record.is_duplicate -> bool`
PCR or optical duplicate

**Example**:
```python
if record.is_mapped and record.is_primary:
    if record.is_paired and record.is_proper_pair:
        print("High-quality paired read")

    if record.is_duplicate:
        print("  (duplicate, skip for counting)")
```

---

### CIGAR Methods (v1.3.0+)

#### `record.cigar -> list[CigarOp]`
Get CIGAR operations as list

#### `record.cigar_string() -> str`
Get human-readable CIGAR string (e.g., "100M2I50M")

#### `record.reference_length() -> int`
Calculate reference span from CIGAR

#### `record.query_length() -> int`
Calculate query span from CIGAR

#### `record.reference_end() -> int | None`
Calculate alignment end position (position + reference_length)

**Example**:
```python
print(f"CIGAR: {record.cigar_string()}")
print(f"Reference span: {record.reference_length()}bp")
print(f"Query span: {record.query_length()}bp")
print(f"Alignment: {record.position}-{record.reference_end()}")

# Detailed CIGAR analysis
for op in record.cigar:
    if op.is_insertion():
        print(f"  Insertion: {op.length}bp at reference")
    elif op.is_deletion():
        print(f"  Deletion: {op.length}bp from reference")
```

---

### Tag Methods (v1.4.0+)

#### `record.get_tag(name: str) -> Tag | None`
Get tag by name (2-character string)

#### `record.has_tag(name: str) -> bool`
Check if tag exists

#### `record.tags() -> list[Tag]`
Get all tags

#### `record.get_int(name: str) -> int | None`
Get integer tag value directly

#### `record.get_string(name: str) -> str | None`
Get string tag value directly

#### `record.edit_distance() -> int | None`
Get NM tag (number of mismatches)

#### `record.alignment_score() -> int | None`
Get AS tag (alignment score)

#### `record.read_group() -> str | None`
Get RG tag (read group identifier)

#### `record.md_string() -> str | None`
Get MD tag (mismatch/deletion string)

**Example**:
```python
# Check for specific tag
if record.has_tag("NM"):
    nm = record.edit_distance()
    print(f"Edit distance: {nm}")

# Get tag value directly
as_score = record.get_int("AS")
rg = record.get_string("RG")

# Iterate all tags
for tag in record.tags():
    print(f"{tag.name}: {tag.value}")
```

---

## BamHeader Class

### Properties

#### `header.reference_count -> int`
Number of reference sequences

#### `header.references -> list[Reference]`
List of reference sequences

#### `header.reference_names -> list[str]`
List of reference names (convenience)

### Methods

#### `header.reference_name(id: int) -> str`
Get reference name by ID

**Example**:
```python
header = reader.header
print(f"References: {header.reference_count}")

# Get reference name
ref_name = header.reference_name(0)  # "chr1"

# Iterate references
for ref in header.references:
    print(f"{ref.name}: {ref.length:,}bp")

# Get all names
names = header.reference_names  # ["chr1", "chr2", ...]
```

---

## Reference Class

### Properties

#### `reference.name -> str`
Reference sequence name (e.g., "chr1", "chrM")

#### `reference.length -> int`
Reference sequence length in bases

---

## CigarOp Class

### Properties

#### `op.op_char -> str`
Operation character: M, I, D, N, S, H, P, =, X

#### `op.length -> int`
Operation length in bases

### Type Checking Methods

#### `op.is_match() -> bool`
Match/mismatch (M)

#### `op.is_insertion() -> bool`
Insertion to reference (I)

#### `op.is_deletion() -> bool`
Deletion from reference (D)

#### `op.is_ref_skip() -> bool`
Skipped region (N, RNA-seq introns)

#### `op.is_soft_clip() -> bool`
Soft clipping (S, present in SEQ)

#### `op.is_hard_clip() -> bool`
Hard clipping (H, absent in SEQ)

#### `op.is_padding() -> bool`
Padding (P, silent deletion)

#### `op.is_seq_match() -> bool`
Sequence match (=)

#### `op.is_seq_mismatch() -> bool`
Sequence mismatch (X)

### Consumption Methods

#### `op.consumes_reference() -> bool`
True if operation advances reference position (M, D, N, =, X)

#### `op.consumes_query() -> bool`
True if operation consumes query sequence (M, I, S, =, X)

**Example**:
```python
for op in record.cigar:
    if op.is_insertion() and op.length >= 5:
        print(f"Large insertion: {op.length}bp")

    if op.is_deletion() and op.length >= 10:
        print(f"Large deletion: {op.length}bp")

    if op.is_ref_skip():
        print(f"Intron: {op.length}bp (RNA-seq)")
```

---

## Tag Classes

### Tag Class

#### Properties

#### `tag.name -> str`
Two-character tag name (e.g., "NM", "AS")

#### `tag.value -> TagValue`
Tag value (typed)

### TagValue Class

#### Type Checking Methods

#### `value.is_int() -> bool`
#### `value.is_float() -> bool`
#### `value.is_string() -> bool`
#### `value.is_char() -> bool`
#### `value.is_array() -> bool`

#### Conversion Methods

#### `value.as_int() -> int | None`
#### `value.as_float() -> float | None`
#### `value.as_string() -> str | None`
#### `value.as_char() -> str | None`
#### `value.as_array() -> list | None`

**Example**:
```python
tag = record.get_tag("NM")
if tag:
    if tag.value.is_int():
        nm = tag.value.as_int()
        print(f"Edit distance: {nm}")
```

---

## SamWriter Class

### Constructor

#### `SamWriter.create(path: str) -> SamWriter`

Create SAM writer for BAM → SAM conversion.

**Example**:
```python
writer = biometal.SamWriter.create("output.sam")
```

### Methods

#### `writer.write_header(header: BamHeader)`
Write SAM header from BamHeader

#### `writer.write_record(record: BamRecord)`
Write BamRecord in SAM format

#### `writer.close()`
Close writer and flush buffers

**Example**:
```python
# BAM → SAM conversion
reader = biometal.BamReader.from_path("input.bam")
writer = biometal.SamWriter.create("output.sam")

# Write header
writer.write_header(reader.header)

# Write filtered records
for record in reader:
    if record.is_primary and record.mapq >= 30:
        writer.write_record(record)

writer.close()
```

---

## Statistics Functions

### Coverage Analysis

#### `calculate_coverage(path: str, reference_id: int, start: int = None, end: int = None) -> dict[int, int]`

Calculate per-position coverage using CIGAR operations.

**Returns**: Dict mapping position → coverage depth

**Example**:
```python
coverage = biometal.calculate_coverage("alignments.bam", reference_id=0, start=1000, end=2000)
print(f"Position 1500: {coverage.get(1500, 0)}× coverage")
```

### MAPQ Distribution

#### `mapq_distribution(path: str, reference_id: int = None) -> dict[int, int]`

Calculate MAPQ score distribution.

**Returns**: Dict mapping MAPQ → count

**Example**:
```python
dist = biometal.mapq_distribution("alignments.bam")
for mapq in sorted(dist.keys()):
    print(f"MAPQ {mapq}: {dist[mapq]:,} reads")
```

### Flag Statistics

#### `count_by_flag(path: str) -> dict[str, int]`

Count reads by flag categories.

**Returns**: Dict with keys: total, mapped, unmapped, primary, secondary, supplementary, paired, proper_pair, forward, reverse, qc_fail, duplicate

**Example**:
```python
stats = biometal.count_by_flag("alignments.bam")
print(f"Mapped: {stats['mapped']:,} / {stats['total']:,}")
print(f"Primary: {stats['primary']:,}")
print(f"Duplicates: {stats['duplicate']:,}")
```

### Paired-End QC (v1.4.0+)

#### `insert_size_distribution(path: str, reference_id: int = None) -> dict[int, int]`

Calculate insert size distribution for properly paired reads.

**Returns**: Dict mapping insert size → count

**Example**:
```python
dist = biometal.insert_size_distribution("alignments.bam")
sizes = list(dist.keys())
mean_insert = sum(s * dist[s] for s in sizes) / sum(dist.values())
print(f"Mean insert size: {mean_insert:.0f}bp")
```

### Alignment Quality (v1.4.0+)

#### `edit_distance_stats(path: str, reference_id: int = None) -> dict`

Analyze NM tag (edit distance) distribution.

**Returns**: Dict with keys: mean, median, min, max, distribution, with_nm_tag, total_records

**Example**:
```python
stats = biometal.edit_distance_stats("alignments.bam")
print(f"Mean mismatches: {stats['mean']:.2f}")
print(f"Median: {stats['median']}")
print(f"Coverage: {stats['with_nm_tag'] / stats['total_records'] * 100:.1f}%")
```

### Variant QC (v1.4.0+)

#### `strand_bias(path: str, reference_id: int, position: int, window_size: int = 1) -> dict`

Calculate strand bias at a genomic position.

**Returns**: Dict with keys: forward, reverse, total, ratio (0.5 = no bias)

**Example**:
```python
bias = biometal.strand_bias("alignments.bam", reference_id=0, position=1000, window_size=100)
print(f"Forward: {bias['forward']}, Reverse: {bias['reverse']}")
print(f"Ratio: {bias['ratio']:.3f}")
if abs(bias['ratio'] - 0.5) > 0.2:
    print("⚠️  Potential strand bias!")
```

### RNA-seq QC (v1.4.0+)

#### `alignment_length_distribution(path: str, reference_id: int = None) -> dict[int, int]`

Calculate alignment length distribution (CIGAR-based).

**Returns**: Dict mapping alignment length → count

**Example**:
```python
dist = biometal.alignment_length_distribution("alignments.bam")
lengths = list(dist.keys())
mean_len = sum(l * dist[l] for l in lengths) / sum(dist.values())
print(f"Mean alignment length: {mean_len:.0f}bp")
```

---

## Error Handling

All functions raise standard Python exceptions:

- `IOError`: File I/O errors, invalid BAM format
- `ValueError`: Invalid parameters (e.g., tag name not 2 characters)
- `RuntimeError`: Internal errors (rare, please report)

**Example**:
```python
try:
    reader = biometal.BamReader.from_path("missing.bam")
except IOError as e:
    print(f"Cannot open file: {e}")

try:
    tag = record.get_tag("TOOLONG")
except ValueError as e:
    print(f"Invalid tag name: {e}")
```

---

## Performance Tips

### Memory Efficiency

✅ **DO**: Stream records one at a time
```python
for record in reader:
    process(record)  # Constant ~5 MB memory
```

❌ **DON'T**: Accumulate all records in memory
```python
records = list(reader)  # ⚠️ Loads entire file into RAM
```

### Filtering Performance

✅ **DO**: Filter early in the loop
```python
for record in reader:
    if not (record.is_mapped and record.is_primary):
        continue  # Skip early
    # Process only relevant records
```

✅ **DO**: Use region queries when possible
```python
for record in reader.query("chr1", 1000, 2000):
    # Only processes region of interest
```

### Tag Access

✅ **DO**: Use convenience methods
```python
nm = record.edit_distance()  # Direct access
```

❌ **AVOID**: Repeated tag lookups
```python
# In a tight loop
for i in range(1000):
    tag = record.get_tag("NM")  # Slow: parses tags each time
```

### Batch Processing

✅ **DO**: Process in batches for statistics
```python
# Collect data in batches
batch_size = 10000
batch = []
for record in reader:
    batch.append(record.mapq)
    if len(batch) >= batch_size:
        process_batch(batch)
        batch = []
```

---

## Complete Example

```python
import biometal

# Open BAM file
reader = biometal.BamReader.from_path("alignments.bam")

# Show header info
print(f"References: {reader.header.reference_count}")
for ref in reader.header.references:
    print(f"  {ref.name}: {ref.length:,}bp")

# Quality control pipeline
high_quality = 0
total = 0
indels_found = []

for record in reader:
    total += 1

    # Filter for high-quality reads
    if record.is_mapped and record.is_primary and record.mapq and record.mapq >= 30:
        high_quality += 1

        # Check edit distance
        edit_dist = record.edit_distance()
        if edit_dist is not None and edit_dist <= 2:
            # Analyze CIGAR for variants
            for op in record.cigar:
                if (op.is_insertion() or op.is_deletion()) and op.length >= 5:
                    indels_found.append({
                        'read': record.name,
                        'type': 'INS' if op.is_insertion() else 'DEL',
                        'length': op.length,
                        'position': record.position
                    })

    # Progress
    if total % 100000 == 0:
        print(f"Processed {total:,} reads...")

# Summary
print(f"\n=== Summary ===")
print(f"Total reads: {total:,}")
print(f"High-quality: {high_quality:,} ({high_quality/total*100:.1f}%)")
print(f"Large indels found: {len(indels_found)}")

# Export filtered reads to SAM
writer = biometal.SamWriter.create("filtered.sam")
writer.write_header(reader.header)

reader2 = biometal.BamReader.from_path("alignments.bam")
for record in reader2:
    if record.is_mapped and record.is_primary and record.mapq >= 30:
        writer.write_record(record)
writer.close()

print("✅ Exported to filtered.sam")
```

---

## See Also

- [Benchmarks](BENCHMARKS.md) - Performance analysis
- [Examples](../examples/) - Production workflows
- [Notebooks](../notebooks/) - Interactive tutorials
- [README](../README.md) - Quick start guide

---

**biometal v1.4.0** | ARM-native streaming bioinformatics | [GitHub](https://github.com/shandley/biometal)
