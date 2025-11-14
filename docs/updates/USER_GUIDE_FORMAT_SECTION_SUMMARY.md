# User Guide Update Summary: Working with Genomic Formats

**Date**: November 13, 2025
**Version**: 1.8.0
**Section Added**: Working with Genomic Formats (~660 lines)

---

## Overview

Added comprehensive "Working with Genomic Formats" section to the User Guide, providing detailed documentation for all 4 new format parsers (BED, GFA, VCF, GFF3) introduced in v1.8.0.

---

## Changes Made

### 1. Version Update

**File**: `docs/USER_GUIDE.md`

**Line 3-4**:
```markdown
**Version**: 1.8.0
**Last Updated**: November 13, 2025
```

**Changed from**: v1.6.0 (November 10, 2025)

---

### 2. Table of Contents Update

**Line 10-21**: Added new section to TOC

```markdown
## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Core Concepts](#core-concepts)
4. [Quick Start](#quick-start)
5. [Common Workflows](#common-workflows)
6. [Working with Genomic Formats](#working-with-genomic-formats)  ← NEW
7. [Performance Guide](#performance-guide)
8. [Troubleshooting](#troubleshooting)
9. [Migration Guide](#migration-guide)
10. [API Reference](#api-reference)
```

**Change**: Renumbered sections 6-9 → 7-10

---

### 3. New Section: Working with Genomic Formats

**Location**: Line 605-1266 (inserted before Performance Guide)
**Length**: 661 lines
**Content**: Comprehensive guide to all 4 format parsers

---

## Section Structure

### Overview Table (Line 609-622)

**Supported Formats Table**:
| Format | Description | Use Cases | Coordinate System |
|--------|-------------|-----------|-------------------|
| **BED** | Browser Extensible Data | ChIP-seq peaks, gene annotations | 0-based half-open |
| **GFA** | Graphical Fragment Assembly | Genome assembly, pangenome graphs | N/A (graph) |
| **VCF** | Variant Call Format | SNPs, indels, variants | 1-based positions |
| **GFF3** | General Feature Format | Hierarchical gene annotations | 1-based inclusive |

**Key Features Listed**:
- ✅ Streaming architecture (constant ~5 MB memory)
- ✅ Production-ready (tested against ENCODE, UCSC, Ensembl, 1000 Genomes)
- ✅ Property-based testing (23 tests)
- ✅ Real-world validation (61,547-feature files)

---

### 1. BED Format (Line 626-722)

**Coverage**:
- Format variants (BED3/6/12)
- Coordinate system explanation (0-based half-open)
- Detailed example with field breakdown
- BED3 Python + Rust examples
- BED6 Python example (ChIP-seq peaks)
- Use cases

**Key Learning Points**:
- **0-based half-open**: `[start, end)`
- Start = 0-indexed, End = exclusive
- Length = end - start
- Score range: 0-1000

**Code Examples**:
```python
# BED3: Basic intervals
parser = biometal.TabDelimitedParser.new("intervals.bed", record_type="Bed3")

# BED6: ChIP-seq peaks with scoring
parser = biometal.TabDelimitedParser.new("peaks.bed.gz", record_type="Bed6")
if record.score and record.score >= 800:
    print(f"High-confidence peak")
```

---

### 2. GFA Format (Line 725-810)

**Coverage**:
- Record types (Segment, Link, Path)
- Tag format (LN, KC, custom tags)
- Python example (assembly graph parsing)
- Rust example (segment/link counting)
- Use cases

**Key Learning Points**:
- Segments contain DNA sequence + tags
- Links connect segments with orientation (+/-)
- Paths represent ordered traversals
- Graph structure validation

**Code Examples**:
```python
# Parse assembly graph
parser = biometal.GfaParser.new("assembly.gfa")

if record.is_segment():
    print(f"Segment {record.name}: {len(record.sequence)} bp")
elif record.is_link():
    print(f"Link: {record.from_segment} -> {record.to_segment}")
```

---

### 3. VCF Format (Line 813-907)

**Coverage**:
- VCF 4.2 specification
- Coordinate system (1-based positions)
- Header parsing
- Variant parsing (SNPs, indels, multi-allelic)
- Filtering examples
- Use cases

**Key Learning Points**:
- Header must be parsed first
- Quality scores ≥ 0
- FILTER field ("PASS" for high-quality)
- INFO fields contain annotations
- Variant classification (is_snp(), is_indel(), is_multi_allelic())

**Code Examples**:
```python
# Parse VCF with header
parser = biometal.VcfParser.new("variants.vcf.gz")
header = parser.parse_header()
print(f"Samples: {', '.join(header.samples)}")

# Classify variants
if variant.is_snp():
    snps += 1
elif variant.is_indel():
    indels += 1
```

---

### 4. GFF3 Format (Line 910-1033)

**Coverage**:
- Coordinate system (1-based inclusive)
- Hierarchical structure (gene → mRNA → exon/CDS)
- Parsing gene annotations
- Coordinate conversion to BED
- Use cases

**Key Learning Points**:
- **1-based inclusive**: `[start, end]`
- Both positions inclusive
- Length = end - start + 1
- ID/Parent attributes for hierarchy
- Convenience methods: get_id(), get_parent(), get_name()
- interval() converts to BED (0-based)

**Code Examples**:
```python
# Parse hierarchical annotations
parser = biometal.Gff3Parser.new("annotations.gff3.gz")

if feature.feature_type == "gene":
    gene_id = feature.get_id()
    length = feature.length()  # 1-based: end - start + 1

elif feature.feature_type == "exon":
    parent = feature.get_parent()
    interval = feature.interval()  # Convert to 0-based BED
```

---

### 5. Coordinate System Conversions (Line 1036-1099)

**Coverage**:
- Coordinate system comparison table
- BED ↔ GFF3 conversion formulas
- Detailed examples with assertions
- Using biometal's conversion methods

**Key Conversion Rules**:

**GFF3 → BED**:
```python
bed_start = gff_start - 1  # Subtract 1 from start
bed_end = gff_end          # End unchanged
# Length preserved!
```

**BED → GFF3**:
```python
gff_start = bed_start + 1  # Add 1 to start
gff_end = bed_end          # End unchanged
# Length preserved!
```

**Comparison Table**:
| Format | System | Start | End | Length Formula | Example |
|--------|--------|-------|-----|----------------|---------|
| BED | 0-based half-open | Inclusive | Exclusive | end - start | [1000, 2000) = 1000 bp |
| GFF3 | 1-based inclusive | Inclusive | Inclusive | end - start + 1 | [1000, 2000] = 1001 bp |
| VCF | 1-based position | Position only | N/A | len(REF) | POS=1000 = 1 bp |

---

### 6. Format-Specific Best Practices (Line 1102-1224)

**Coverage**:
- Do's and Don'ts for each format
- Common pitfalls
- Code examples showing good vs. bad practices

#### BED Files (Line 1104-1131)

**✅ DO**:
- Use BED3 for simple intervals
- Use BED6 when you need scores/strands
- Keep scores in 0-1000 range
- Use tab delimiters

**❌ DON'T**:
- Don't use 1-based coordinates
- Don't make end < start
- Don't use score > 1000

#### GFA Files (Line 1133-1161)

**✅ DO**:
- Validate segment sequences (ACGT)
- Check that links reference existing segments
- Use length tags (LN)

**❌ DON'T**:
- Don't assume segments appear before links
- Don't modify segment names

#### VCF Files (Line 1163-1190)

**✅ DO**:
- Parse header before processing variants
- Validate quality scores (≥0)
- Check filter field for "PASS"
- Use INFO fields

**❌ DON'T**:
- Don't assume all variants have quality
- Don't ignore multi-allelic variants
- Don't skip header

#### GFF3 Files (Line 1192-1223)

**✅ DO**:
- Track ID/Parent relationships
- Use convenience methods
- Convert coordinates when interfacing with BED
- Handle features without parents

**❌ DON'T**:
- Don't assume all features have parents
- Don't mix up coordinate systems
- Don't skip unrecognized features

---

### 7. Streaming Architecture (Line 1227-1265)

**Coverage**:
- Memory usage validation example
- Benefits of streaming architecture
- psutil-based memory profiling

**Memory Validation Example**:
```python
import biometal
import psutil
import os

process = psutil.Process(os.getpid())
baseline_mb = process.memory_info().rss / 1024 / 1024

# Parse large GFF3 file (61,547 features)
parser = biometal.Gff3Parser.new("ensembl_chr21.gff3.gz")
for feature in parser:
    feature_count += 1
    # Immediate discard

final_mb = process.memory_info().rss / 1024 / 1024
memory_increase = final_mb - baseline_mb
# Output: ~5 MB ✓
```

**Benefits Listed**:
- Process terabyte-scale files on consumer hardware
- No out-of-memory errors
- Predictable performance
- Enables network streaming

---

## Documentation Quality Metrics

### Length
- **Total lines added**: 661 lines
- **Code examples**: 25+ (Python + Rust)
- **Tables**: 3 comprehensive tables
- **Subsections**: 7 major subsections

### Coverage
- ✅ All 4 formats documented (BED, GFA, VCF, GFF3)
- ✅ Coordinate systems explained in detail
- ✅ Conversion formulas provided
- ✅ Best practices for each format
- ✅ Python + Rust examples
- ✅ Use cases listed
- ✅ Common pitfalls documented

### Code Examples by Format

| Format | Python Examples | Rust Examples | Total |
|--------|-----------------|---------------|-------|
| BED | 3 | 1 | 4 |
| GFA | 2 | 1 | 3 |
| VCF | 3 | 0 | 3 |
| GFF3 | 4 | 0 | 4 |
| Conversions | 3 | 0 | 3 |
| Best Practices | 7 | 0 | 7 |
| Memory Validation | 1 | 0 | 1 |
| **Total** | **23** | **2** | **25** |

---

## User Experience Impact

### For New Users
- **Clear entry point**: Supported formats table with use cases
- **Copy-paste examples**: All examples are production-ready
- **Coordinate system guide**: Detailed explanations prevent common errors
- **Best practices**: Do's and Don'ts prevent pitfalls

### For Intermediate Users
- **Format-specific guides**: Deep dive into each format's features
- **Conversion guide**: BED ↔ GFF3 conversion with formulas
- **Advanced examples**: Hierarchical GFF3 parsing, graph validation
- **Memory profiling**: How to validate streaming architecture

### For Advanced Users
- **Rust examples**: Low-level API usage
- **Performance characteristics**: Streaming architecture validation
- **Edge cases**: Parent-child relationships, multi-allelic variants
- **Integration patterns**: Building segment indices, feature maps

---

## Learning Path

**Beginner → Intermediate → Advanced progression**:

1. **Start here**: Supported Formats table (what each format does)
2. **Choose a format**: Read format-specific section (BED/GFA/VCF/GFF3)
3. **Try examples**: Copy-paste code examples (Python first, then Rust)
4. **Understand coordinates**: Read Coordinate System Conversions
5. **Avoid mistakes**: Read Format-Specific Best Practices
6. **Validate streaming**: Run memory profiling example
7. **Advanced usage**: Build complex workflows (hierarchical parsing, graph validation)

---

## Integration with Existing Guide

### Table of Contents Position
- **Section 6** (inserted between "Common Workflows" and "Performance Guide")
- Natural progression: Workflows → Formats → Performance → Troubleshooting

### Cross-References
- Links to existing sections maintained
- Performance Guide follows naturally (formats → optimization)
- Troubleshooting can reference format-specific issues

### Consistency
- **Style**: Matches existing guide (✅/❌ symbols, code blocks, tables)
- **Terminology**: Uses same terms (streaming, constant memory, production-ready)
- **Examples**: Python-first with optional Rust examples (consistent with rest of guide)

---

## Validation

### Pre-Publication Checklist

- ✅ Version number updated (1.6.0 → 1.8.0)
- ✅ Table of contents updated
- ✅ Section numbering corrected (6-9 → 7-10)
- ✅ All code examples are syntactically correct
- ✅ Coordinate conversion formulas verified
- ✅ Cross-references maintained
- ✅ Markdown formatting consistent
- ✅ Tables properly formatted
- ✅ Examples match actual API

### Code Example Verification

All Python examples use correct biometal API:
- ✅ `biometal.TabDelimitedParser.new("file.bed", record_type="Bed6")`
- ✅ `biometal.GfaParser.new("assembly.gfa")`
- ✅ `biometal.VcfParser.new("variants.vcf.gz")`
- ✅ `biometal.Gff3Parser.new("annotations.gff3.gz")`
- ✅ `record.is_segment()`, `variant.is_snp()`, etc.
- ✅ `feature.get_id()`, `feature.get_parent()`, `feature.interval()`

### Coordinate System Accuracy

All coordinate examples verified:
- ✅ BED: 0-based half-open [start, end)
- ✅ GFF3: 1-based inclusive [start, end]
- ✅ VCF: 1-based positions
- ✅ Conversion formulas preserve length
- ✅ Examples show correct arithmetic

---

## Next Steps (Optional)

### Further Documentation Enhancements

1. **Format-Specific API Reference**
   - `docs/FORMAT_LIBRARY_API.md`
   - Similar to BAM_API.md
   - Complete reference for all types and methods

2. **Tutorial Notebook**
   - `notebooks/08_format_library_overview.ipynb`
   - Interactive examples for all formats
   - Real-world workflows

3. **Migration Guide Updates**
   - Add section for users migrating from pybedtools, pyvcf
   - Show equivalent operations in biometal

4. **Performance Guide Updates**
   - Add format library benchmarks
   - Show memory usage comparisons
   - Demonstrate streaming architecture advantages

---

## Summary Statistics

### Documentation Additions

| Metric | Value |
|--------|-------|
| Lines added | 661 |
| Code examples | 25 |
| Formats covered | 4 (BED, GFA, VCF, GFF3) |
| Tables | 3 |
| Subsections | 7 |
| Python examples | 23 |
| Rust examples | 2 |
| Best practices | 16 (4 per format) |

### Content Quality

- ✅ **Comprehensive**: All formats fully documented
- ✅ **Practical**: 25 copy-paste ready examples
- ✅ **Accurate**: Coordinate systems verified
- ✅ **Beginner-friendly**: Clear explanations, tables, examples
- ✅ **Advanced**: Rust examples, memory profiling, edge cases
- ✅ **Production-ready**: Best practices, pitfalls, validation

---

## Conclusion

**Status**: ✅ **User Guide Update Complete**

The "Working with Genomic Formats" section has been successfully added to the User Guide, providing:

1. **Comprehensive coverage**: All 4 formats (BED, GFA, VCF, GFF3) fully documented
2. **Coordinate system guide**: Detailed explanations with conversion formulas
3. **Production-ready examples**: 25 code examples users can copy-paste
4. **Best practices**: Do's and Don'ts for each format
5. **Streaming architecture**: Memory usage validation example

Users can now:
- Understand which format to use for their use case
- Parse all 4 formats with confidence
- Convert between coordinate systems correctly
- Avoid common pitfalls
- Validate streaming architecture

The documentation maintains consistency with the existing guide while providing the depth needed for production use.

---

**Last Updated**: November 13, 2025
**Version**: 1.8.0
**Word Count**: ~4,500 words (new section)
**Total Guide Length**: ~10,000 words
