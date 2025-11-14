# Python Bindings Verification Report

**Date**: November 13, 2025
**Version**: 1.8.0
**Status**: ✅ Bindings Exist, ⚠️ Documentation Needs Updates

---

## Executive Summary

All format library Python bindings (BED, GFA, VCF, GFF3) **ARE** registered and exposed in the Python module. However, the documentation examples use an **incorrect API** that doesn't match the actual implementation.

### Findings

✅ **GOOD**: All parsers registered in `src/python/mod.rs` (lines 73-94)
⚠️ **ISSUE**: Documentation uses non-existent API (`TabDelimitedParser.new()`, `GfaParser.new()`, etc.)
✅ **SOLUTION**: Update documentation to use correct APIs (`Bed6Stream.from_path()`, `GfaStream.from_path()`, etc.)

---

## Python Bindings Status

### Registered Classes (src/python/mod.rs)

All format parsers are registered in the `biometal` Python module:

```rust
// BED format types (lines 73-79)
m.add_class::<PyBed3Record>()?;
m.add_class::<PyBed6Record>()?;
m.add_class::<PyBed12Record>()?;
m.add_class::<PyBed3Stream>()?;
m.add_class::<PyBed6Stream>()?;
m.add_class::<PyBed12Stream>()?;

// GFA format types (lines 81-85)
m.add_class::<PyGfaSegment>()?;
m.add_class::<PyGfaLink>()?;
m.add_class::<PyGfaPath>()?;
m.add_class::<PyGfaStream>()?;

// VCF format types (lines 87-90)
m.add_class::<PyVcfHeader>()?;
m.add_class::<PyVcfRecord>()?;
m.add_class::<PyVcfStream>()?;

// GFF3 format types (lines 92-94)
m.add_class::<PyGff3Record>()?;
m.add_class::<PyGff3Stream>()?;
```

---

## Correct Python API

### BED Format

**Exposed Classes**:
- `Bed3Record`: BED3 record type
- `Bed6Record`: BED6 record type
- `Bed12Record`: BED12 record type
- `Bed3Stream`: Iterator for BED3 files
- `Bed6Stream`: Iterator for BED6 files
- `Bed12Stream`: Iterator for BED12 files

**Correct API** (from src/python/bed.rs):

```python
import biometal

# BED3: Basic intervals
stream = biometal.Bed3Stream.from_path("intervals.bed")
for record in stream:
    print(f"{record.chrom}:{record.start}-{record.end}")
    length = record.length()  # Method available

# BED6: Extended intervals with name, score, strand
stream = biometal.Bed6Stream.from_path("peaks.bed.gz")
for record in stream:
    print(f"{record.chrom}:{record.start}-{record.end}")
    print(f"  Name: {record.name}, Score: {record.score}, Strand: {record.strand}")
    length = record.length()  # Method available

# BED12: Full format with blocks
stream = biometal.Bed12Stream.from_path("transcripts.bed")
for record in stream:
    print(f"{record.name}: {record.block_count} exons")
    length = record.length()  # Method available
```

**Record Attributes**:
- `Bed3Record`: `chrom`, `start`, `end`
- `Bed6Record`: `chrom`, `start`, `end`, `name`, `score`, `strand`
- `Bed12Record`: `chrom`, `start`, `end`, `name`, `score`, `strand`, `thick_start`, `thick_end`, `item_rgb`, `block_count`, `block_sizes`, `block_starts`

**Methods**:
- All records: `length()`, `to_line()`, `from_line(line)`

---

### GFA Format

**Exposed Classes**:
- `GfaSegment`: Segment (contig/node)
- `GfaLink`: Link (edge)
- `GfaPath`: Path (ordered traversal)
- `GfaStream`: Iterator for GFA files

**Correct API** (from src/python/gfa.rs):

```python
import biometal

# Parse GFA assembly graph
stream = biometal.GfaStream.from_path("assembly.gfa")

segments = []
links = []
paths = []

for record in stream:
    # Record type is determined dynamically
    if isinstance(record, biometal.GfaSegment):
        segments.append(record)
        print(f"Segment {record.name}: {len(record.sequence)} bp")
        print(f"  Tags: {record.tags}")
        print(f"  Length: {record.length()}")
        print(f"  Coverage: {record.coverage()}")

    elif isinstance(record, biometal.GfaLink):
        links.append(record)
        print(f"Link: {record.from_segment}{record.from_orient} -> {record.to_segment}{record.to_orient}")
        print(f"  Overlap: {record.overlap}")

    elif isinstance(record, biometal.GfaPath):
        paths.append(record)
        print(f"Path {record.name}: {len(record.segments)} segments")

print(f"Graph: {len(segments)} segments, {len(links)} links, {len(paths)} paths")
```

**Record Attributes**:
- `GfaSegment`: `name`, `sequence`, `tags` (dict)
- `GfaLink`: `from_segment`, `from_orient`, `to_segment`, `to_orient`, `overlap`, `tags`
- `GfaPath`: `name`, `segments` (list), `overlaps` (list), `tags`

**Methods**:
- `GfaSegment`: `length()`, `coverage()`, `from_line(line)`
- `GfaLink`: `from_line(line)`
- `GfaPath`: `from_line(line)`

---

### VCF Format

**Exposed Classes**:
- `VcfHeader`: VCF header with metadata
- `VcfRecord`: Variant record
- `VcfStream`: Iterator for VCF files

**Correct API** (from src/python/vcf.rs):

```python
import biometal

# Parse VCF file
# Note: Python bindings don't support gzip yet - use uncompressed files
stream = biometal.VcfStream.from_path("variants.vcf")

# Get header (automatically parsed on file open)
header = stream.header()  # Note: header() not parse_header()
print(f"VCF version: {header.fileformat}")
print(f"Samples: {', '.join(header.samples)}")
print(f"Contigs: {len(header.contigs)}")
print(f"INFO fields: {len(header.info_fields)}")

# Then iterate variants
for variant in stream:
    print(f"{variant.chrom}:{variant.pos} {variant.reference}→{variant.alternate}")

    # Access optional fields
    if variant.id:
        print(f"  ID: {variant.id}")
    if variant.quality is not None:
        print(f"  Quality: {variant.quality}")
    if variant.filter:
        print(f"  Filter: {variant.filter}")

    # Access INFO fields
    for key, value in variant.info.items():
        print(f"  {key}={value}")

    # Classify variants
    if variant.is_snp():
        print("  Type: SNP")
    elif variant.is_indel():
        print("  Type: Indel")

    if variant.is_multi_allelic():
        print(f"  Multi-allelic: {len(variant.alternate)} alleles")
```

**Header Attributes**:
- `fileformat`: VCF version string
- `samples`: List of sample names
- `contigs`: Dictionary of contig definitions
- `info_fields`: Dictionary of INFO field definitions
- `format_fields`: Dictionary of FORMAT field definitions
- `filters`: Dictionary of FILTER definitions

**Record Attributes**:
- `chrom`: Chromosome name
- `pos`: Position (1-based)
- `id`: Variant ID (optional)
- `reference`: Reference allele (Note: `reference` not `ref`)
- `alternate`: List of alternate alleles (Note: `alternate` not `alt`)
- `quality`: Quality score (optional float)
- `filter`: Filter status (optional string)
- `info`: Dictionary of INFO fields
- `format`: FORMAT column (optional)
- `samples`: List of sample values

**Methods**:
- `is_snp()`: Check if variant is a SNP
- `is_indel()`: Check if variant is an indel
- `is_multi_allelic()`: Check if variant has multiple alternates
- `from_line(line)`: Parse from VCF line

---

### GFF3 Format

**Exposed Classes**:
- `Gff3Record`: GFF3 feature record
- `Gff3Stream`: Iterator for GFF3 files

**Correct API** (from src/python/gff.rs):

```python
import biometal

# Parse GFF3 annotations
stream = biometal.Gff3Stream.from_path("annotations.gff3.gz")

genes = {}
exons = []
cds_features = []

for feature in stream:
    # Access record fields
    print(f"{feature.seqid}:{feature.start}-{feature.end} ({feature.feature_type})")
    print(f"  Source: {feature.source}")
    print(f"  Strand: {feature.strand}")

    # Optional fields
    if feature.score is not None:
        print(f"  Score: {feature.score}")
    if feature.phase is not None:
        print(f"  Phase: {feature.phase}")

    # Access attributes (ID, Parent, Name, etc.)
    for key, value in feature.attributes.items():
        print(f"  {key}={value}")

    # Convenience methods
    feature_id = feature.get_id()  # Get ID attribute
    parent_id = feature.get_parent()  # Get Parent attribute
    name = feature.get_name()  # Get Name attribute

    # Calculate length (1-based inclusive)
    length = feature.length()  # end - start + 1

    # Process by type
    if feature.feature_type == "gene":
        if feature_id:
            genes[feature_id] = feature

    elif feature.feature_type == "exon":
        if parent_id:
            exons.append((parent_id, feature))

    elif feature.feature_type == "CDS":
        if parent_id:
            cds_features.append((parent_id, feature))

print(f"Parsed {len(genes)} genes, {len(exons)} exons, {len(cds_features)} CDS")
```

**Record Attributes**:
- `seqid`: Sequence ID (chromosome)
- `source`: Annotation source
- `feature_type`: Feature type (gene, exon, CDS, etc.)
- `start`: Start position (1-based, inclusive)
- `end`: End position (1-based, inclusive)
- `score`: Score (optional float)
- `strand`: Strand (+, -, .)
- `phase`: Reading frame phase (optional 0/1/2)
- `attributes`: Dictionary of attributes

**Methods**:
- `get_id()`: Get ID attribute (returns string or None)
- `get_parent()`: Get Parent attribute (returns string or None)
- `get_name()`: Get Name attribute (returns string or None)
- `length()`: Calculate feature length (end - start + 1)
- `from_line(line)`: Parse from GFF3 line

---

## Documentation Issues Found

### Issue 1: Incorrect API in README.md

**Location**: README.md lines 424-462

**Current (WRONG)**:
```python
# BED: Parse genomic intervals (ChIP-seq peaks, gene annotations)
bed_parser = biometal.TabDelimitedParser.new("peaks.bed.gz", record_type="Bed6")
for record in bed_parser:
    print(f"{record.chrom}:{record.start}-{record.end} score={record.score}")

# GFA: Parse assembly graphs (genome assembly, pangenomes)
gfa_parser = biometal.GfaParser.new("assembly.gfa")

# VCF: Parse genetic variants (SNPs, indels)
vcf_parser = biometal.VcfParser.new("variants.vcf.gz")

# GFF3: Parse hierarchical gene annotations
gff_parser = biometal.Gff3Parser.new("annotations.gff3.gz")
```

**Should Be (CORRECT)**:
```python
# BED: Parse genomic intervals (ChIP-seq peaks, gene annotations)
stream = biometal.Bed6Stream.from_path("peaks.bed.gz")
for record in stream:
    print(f"{record.chrom}:{record.start}-{record.end} score={record.score}")

# GFA: Parse assembly graphs (genome assembly, pangenomes)
stream = biometal.GfaStream.from_path("assembly.gfa")

# VCF: Parse genetic variants (SNPs, indels)
stream = biometal.VcfStream.from_path("variants.vcf")
header = stream.header()  # Get header

# GFF3: Parse hierarchical gene annotations
stream = biometal.Gff3Stream.from_path("annotations.gff3.gz")
```

---

### Issue 2: Incorrect API in User Guide

**Location**: docs/USER_GUIDE.md lines 660-670, 697-715, 746-777, 836-877, 940-988

All examples use non-existent API:
- `TabDelimitedParser.new()` → Should be `Bed3Stream.from_path()`/`Bed6Stream.from_path()`/`Bed12Stream.from_path()`
- `GfaParser.new()` → Should be `GfaStream.from_path()`
- `VcfParser.new()` → Should be `VcfStream.from_path()`
- `Gff3Parser.new()` → Should be `Gff3Stream.from_path()`

---

## Required Documentation Updates

### Files to Update

1. **README.md**
   - Lines 424-462: Format Library usage examples
   - Replace all parser instances with correct Stream API

2. **docs/USER_GUIDE.md**
   - Lines 660-670: BED3 Python example
   - Lines 697-715: BED6 Python example
   - Lines 746-777: GFA Python example
   - Lines 836-877: VCF Python example
   - Lines 940-988: GFF3 Python example
   - Lines 1122-1130: BED best practices example
   - Lines 1150-1160: GFA best practices example
   - Lines 1181-1189: VCF best practices example
   - Lines 1211-1222: GFF3 best practices example

3. **examples/*.py** (if they exist)
   - Update all Python examples with correct API

### Correction Pattern

**Find**:
```python
parser = biometal.TabDelimitedParser.new("file.bed", record_type="Bed6")
```

**Replace with**:
```python
stream = biometal.Bed6Stream.from_path("file.bed")
```

**Find**:
```python
parser = biometal.GfaParser.new("file.gfa")
```

**Replace with**:
```python
stream = biometal.GfaStream.from_path("file.gfa")
```

**Find**:
```python
parser = biometal.VcfParser.new("file.vcf.gz")
```

**Replace with**:
```python
stream = biometal.VcfStream.from_path("file.vcf")
header = stream.header()  # Get header
```

**Find**:
```python
parser = biometal.Gff3Parser.new("file.gff3.gz")
```

**Replace with**:
```python
stream = biometal.Gff3Stream.from_path("file.gff3.gz")
```

---

## Verification Checklist

### Bindings Exist ✅

- [x] BED: Bed3Record, Bed6Record, Bed12Record, Bed3Stream, Bed6Stream, Bed12Stream
- [x] GFA: GfaSegment, GfaLink, GfaPath, GfaStream
- [x] VCF: VcfHeader, VcfRecord, VcfStream
- [x] GFF3: Gff3Record, Gff3Stream

### Bindings Registered ✅

- [x] All classes added to PyModule in src/python/mod.rs (lines 73-94)
- [x] Python docstrings present in source files
- [x] Example usage in docstrings

### Documentation Needs Updates ⚠️

- [ ] README.md: Update format library examples
- [ ] User Guide: Update BED examples (3 sections)
- [ ] User Guide: Update GFA examples (2 sections)
- [ ] User Guide: Update VCF examples (2 sections)
- [ ] User Guide: Update GFF3 examples (3 sections)
- [ ] User Guide: Update best practices examples (4 sections)

---

## Recommended Actions

### High Priority (Must Do Before Release)

1. **Update README.md**
   - Correct all format library examples (lines 424-462)
   - Ensure examples are copy-paste ready

2. **Update User Guide**
   - Correct all Python examples (~10 code blocks)
   - Verify all examples use correct API

3. **Add Python Test**
   - Create `tests/python/test_format_library.py`
   - Test BED, GFA, VCF, GFF3 parsing with real files
   - Verify all documented examples actually work

### Medium Priority (Nice to Have)

4. **Example Scripts**
   - Create `examples/bed_parser.py` (correct API)
   - Create `examples/gfa_parser.py` (correct API)
   - Create `examples/vcf_parser.py` (correct API)
   - Create `examples/gff3_parser.py` (correct API)

5. **API Documentation**
   - Document in docs/FORMAT_LIBRARY_API.md
   - Python API reference section
   - All classes, attributes, methods

---

## Conclusion

**Status**: ⚠️ **Bindings Exist but Documentation is Incorrect**

### Summary

✅ **GOOD**: All format library parsers (BED, GFA, VCF, GFF3) are fully implemented and registered in Python bindings

⚠️ **ISSUE**: Documentation uses non-existent API that doesn't match actual implementation

✅ **SOLUTION**: Update documentation to use correct API (`Bed6Stream.from_path()`, etc.)

### Impact

**Before Fix**:
- Users will get `AttributeError` when trying documented examples
- Examples are completely broken
- No user can actually use the format library from Python

**After Fix**:
- All examples will work as documented
- Users can copy-paste examples successfully
- Format library becomes usable from Python

### Estimated Effort

- **README.md updates**: 10-15 minutes
- **User Guide updates**: 30-45 minutes (10 code blocks)
- **Testing**: 15-20 minutes
- **Total**: ~1-1.5 hours

---

**Last Updated**: November 13, 2025
**Next Step**: Update documentation with correct API
