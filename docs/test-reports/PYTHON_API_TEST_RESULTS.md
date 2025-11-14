# Python API Test Results

**Date**: November 13, 2025
**Version**: 1.8.0
**Status**: ✅ All Format Parsers Verified Working

---

## Executive Summary

All Python bindings for the format library (BED, GFA, VCF, GFF3) have been successfully tested and verified working. However, testing revealed **additional API inconsistencies** between documentation and actual implementation that require fixes.

### Test Results
- ✅ **BED bindings**: PASS
- ✅ **GFA bindings**: PASS
- ✅ **VCF bindings**: PASS
- ✅ **GFF3 bindings**: PASS
- ✅ **File parsing**: PASS (all 4 formats)

**Total**: 5/5 tests passed 🎉

---

## Key Findings

### 1. ✅ All Bindings Exist and Work

All classes are correctly registered and functional:
- **BED**: `Bed3Stream`, `Bed6Stream`, `Bed12Stream` + record types
- **GFA**: `GfaStream`, `GfaSegment`, `GfaLink`, `GfaPath`
- **VCF**: `VcfStream`, `VcfHeader`, `VcfRecord`
- **GFF3**: `Gff3Stream`, `Gff3Record`

### 2. ⚠️ Additional API Inconsistencies Found

Testing revealed **additional differences** between documented API and actual implementation:

#### VCF Header Access
**Documented (WRONG)**:
```python
stream = biometal.VcfStream.from_path("variants.vcf.gz")
header = stream.parse_header()  # Method doesn't exist!
```

**Actual (CORRECT)**:
```python
stream = biometal.VcfStream.from_path("variants.vcf")
header = stream.header()  # Method (not property)
```

**Issue**: Documentation says `parse_header()` but actual method is `header()`

#### VCF Record Attributes
**Documented (WRONG)**:
```python
print(variant.ref)  # Attribute doesn't exist!
print(variant.alt)  # Attribute doesn't exist!
```

**Actual (CORRECT)**:
```python
print(variant.reference)  # Correct attribute name
print(variant.alternate)  # Correct attribute name
```

**Issue**: Documentation uses `ref`/`alt` but actual attributes are `reference`/`alternate`

### 3. ⚠️ Gzip Support Not Available in Python

**Finding**: Python bindings do **not** support gzipped files, only uncompressed files work.

**Test**:
```python
# ❌ FAILS with "I/O error: stream did not contain valid UTF-8"
stream = biometal.Bed6Stream.from_path("peaks.bed.gz")

# ✅ WORKS
stream = biometal.Bed6Stream.from_path("peaks.bed")
```

**Impact**: All documentation examples using `.gz` files will fail for Python users.

**Workaround**: Users must decompress files first:
```bash
gunzip -c file.bed.gz > file.bed
python script.py  # Use uncompressed file
```

---

## Verified Working APIs

### BED Format ✅
```python
import biometal

# BED6 parsing (uncompressed only)
stream = biometal.Bed6Stream.from_path("peaks.bed")
for record in stream:
    print(f"{record.chrom}:{record.start}-{record.end}")
    print(f"  Score: {record.score}, Strand: {record.strand}")
    length = record.length()  # Method works
```

**Attributes**: `chrom`, `start`, `end`, `name`, `score`, `strand`
**Methods**: `length()`, `to_line()`, `from_line()`

### GFA Format ✅
```python
import biometal

# GFA parsing
stream = biometal.GfaStream.from_path("assembly.gfa")
for record in stream:
    if isinstance(record, biometal.GfaSegment):
        print(f"Segment: {record.name}, {len(record.sequence)} bp")
    elif isinstance(record, biometal.GfaLink):
        print(f"Link: {record.from_segment} -> {record.to_segment}")
    elif isinstance(record, biometal.GfaPath):
        print(f"Path: {record.name}")
```

**Record Types**: Dynamic (use `isinstance()` to check type)
**Segment Attributes**: `name`, `sequence`, `tags`
**Link Attributes**: `from_segment`, `from_orient`, `to_segment`, `to_orient`, `overlap`

### VCF Format ✅
```python
import biometal

# VCF parsing (uncompressed only)
stream = biometal.VcfStream.from_path("variants.vcf")

# Get header (method, not property)
header = stream.header()
print(f"VCF version: {header.fileformat}")
print(f"Samples: {header.samples}")

# Iterate variants
for variant in stream:
    print(f"{variant.chrom}:{variant.pos}")
    print(f"  REF: {variant.reference}")  # Note: 'reference' not 'ref'
    print(f"  ALT: {variant.alternate}")  # Note: 'alternate' not 'alt'

    # Variant classification
    if variant.is_snp():
        print("  Type: SNP")
    elif variant.is_indel():
        print("  Type: Indel")
```

**Header Attributes**: `fileformat`, `samples`, `contigs`, `info_fields`, `format_fields`, `filters`
**Record Attributes**: `chrom`, `pos`, `reference`, `alternate`, `id`, `quality`, `filter`, `info`, `format`, `samples`
**Methods**: `is_snp()`, `is_indel()`, `is_insertion()`, `is_deletion()`, `from_line()`, `to_line()`

### GFF3 Format ✅
```python
import biometal

# GFF3 parsing (uncompressed only)
stream = biometal.Gff3Stream.from_path("annotations.gff3")
for feature in stream:
    print(f"{feature.seqid}:{feature.start}-{feature.end}")
    print(f"  Type: {feature.feature_type}")

    # Convenience methods
    gene_id = feature.get_id()
    parent_id = feature.get_parent()
    name = feature.get_name()
    length = feature.length()  # 1-based inclusive
```

**Attributes**: `seqid`, `source`, `feature_type`, `start`, `end`, `score`, `strand`, `phase`, `attributes`
**Methods**: `get_id()`, `get_parent()`, `get_name()`, `length()`, `from_line()`

---

## Required Documentation Fixes

Based on actual API testing, the following documentation updates are required:

### 1. VCF Header Access

**Files to Fix**: README.md, docs/USER_GUIDE.md (all VCF examples)

**Find**:
```python
header = stream.parse_header()
```

**Replace with**:
```python
header = stream.header()  # Note: header() is a method, not parse_header()
```

**Affected Lines**:
- README.md: Line ~443
- USER_GUIDE.md: Lines 842, 866, 1188 (3 locations)

### 2. VCF Record Attributes

**Files to Fix**: README.md, docs/USER_GUIDE.md, PYTHON_BINDINGS_VERIFICATION.md

**Find**:
```python
variant.ref
variant.alt
```

**Replace with**:
```python
variant.reference  # Note: 'reference' not 'ref'
variant.alternate  # Note: 'alternate' not 'alt'
```

**Affected Lines**:
- README.md: Line ~445
- USER_GUIDE.md: Lines 852-869 (VCF parsing examples)
- USER_GUIDE.md: Lines 1188-1194 (VCF best practices)
- PYTHON_BINDINGS_VERIFICATION.md: Lines 179, 196, 199

### 3. Gzip Support Warning

**Action**: Add warning to all Python examples that gzipped files are not yet supported.

**Recommended Note**:
```python
# Note: Python bindings do not support gzipped files yet.
# Decompress first: gunzip -c file.bed.gz > file.bed
stream = biometal.Bed6Stream.from_path("peaks.bed")  # Use uncompressed file
```

**Add to**:
- README.md: After Python examples section
- USER_GUIDE.md: At start of "Working with Genomic Formats" section
- PYTHON_BINDINGS_VERIFICATION.md: In each format's "Correct API" section

---

## Test Environment

**Setup**:
```bash
# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install maturin and build bindings
pip install maturin
maturin develop --release

# Run tests
python test_python_bindings.py
```

**Python Version**: 3.14.0
**Maturin Version**: 1.10.1
**Platform**: macOS 14.6.0 (Darwin 24.6.0), ARM64

---

## Test Files Used

All tests used real-world data files:
- **BED6**: `tests/data/real_world/encode_peaks.bed.gz` (10 peaks, decompressed)
- **GFA**: `tests/data/real_world/lambda_phage.gfa` (5 segments, 5 links, 2 paths)
- **VCF**: `tests/data/real_world/synthetic_1000g.vcf.gz` (10 variants, decompressed)
- **GFF3**: `tests/data/real_world/ensembl_chr21.gff3.gz` (61,547 features, first 1,000 tested, decompressed)

---

## Recommendations

### High Priority (Must Fix)

1. **Update all VCF examples**:
   - Replace `parse_header()` → `header()`
   - Replace `ref`/`alt` → `reference`/`alternate`

2. **Add gzip warning to Python documentation**:
   - Explain limitation
   - Provide workaround (decompress first)

3. **Update PYTHON_BINDINGS_VERIFICATION.md**:
   - Correct VCF API examples
   - Add gzip limitation note

### Medium Priority (Nice to Have)

4. **Add gzip support to Python bindings**:
   - Currently missing in PyO3 wrapper
   - Rust code has gzip support (via flate2/libdeflate)
   - Need to expose decompression to Python layer

5. **Create comprehensive Python examples**:
   - `examples/python/bed_parser.py` (correct, working example)
   - `examples/python/gfa_parser.py`
   - `examples/python/vcf_parser.py`
   - `examples/python/gff3_parser.py`

6. **Add Python-specific test suite**:
   - Move `test_python_bindings.py` to `tests/python/`
   - Expand test coverage
   - Add to CI/CD pipeline

---

## Summary Statistics

### API Verification
- **Classes tested**: 13
- **Methods tested**: 20+
- **Attributes tested**: 30+
- **File formats tested**: 4 (BED, GFA, VCF, GFF3)
- **Test files parsed**: 4 (1 per format)
- **Records parsed**: 1,030 total
  - BED6: 10 records
  - GFA: 12 records (5 segments + 5 links + 2 paths)
  - VCF: 10 variants
  - GFF3: 1,000 features (22 genes)

### Documentation Issues Found
- **Original issues** (from PYTHON_BINDINGS_VERIFICATION.md): 2
  - TabDelimitedParser/GfaParser/VcfParser/Gff3Parser → Stream classes
  - GFA type checking (is_segment() → isinstance())
- **New issues** (from actual testing): 3
  - VCF: parse_header() → header()
  - VCF: ref/alt → reference/alternate
  - Gzip support missing

**Total documentation fixes needed**: 5 API corrections + gzip warning

---

## Conclusion

**Status**: ✅ **Python Bindings Fully Functional**

All format library Python bindings work correctly once the proper API is used. Testing revealed:

1. ✅ **Good**: All classes registered and functional
2. ✅ **Good**: Streaming architecture works (constant memory)
3. ✅ **Good**: All parsing logic correct
4. ⚠️ **Issue**: VCF API differs from documentation (header(), reference/alternate)
5. ⚠️ **Limitation**: Gzip support not available in Python (workaround: decompress first)

**Next Steps**:
1. Update VCF documentation examples (3-5 locations)
2. Add gzip limitation warning to Python docs
3. Consider adding gzip support to Python bindings (future enhancement)

**Impact on Users**:
- **After doc fixes**: All Python examples will work correctly with uncompressed files
- **Current workaround**: Users can decompress .gz files before using Python bindings
- **Future**: Add gzip support to Python layer for seamless experience

---

**Last Updated**: November 13, 2025
**Test Status**: All 5/5 tests passing
**Documentation Status**: Additional fixes needed (VCF API corrections)
