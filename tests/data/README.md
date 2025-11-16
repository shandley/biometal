# biometal Test Data Directory

This directory contains test datasets for validating biometal's file format parsers and operations. Files are organized by source (real-world vs synthetic) and format type.

**Last Updated**: November 16, 2025
**Version**: v1.12.0

---

## 📋 Directory Structure

```
tests/data/
├── README.md                    # This file
├── real_world/                  # Real data from public datasets
│   ├── alignment/              # BAM, CRAM, SAM files
│   ├── variants/               # VCF, BCF files
│   ├── annotation/             # BED, GFF3, GTF files
│   ├── assembly/               # GFA, PAF files
│   ├── sequence/               # FASTQ, FASTA files
│   └── indices/                # BAI, CSI, TBI, FAI files
├── synthetic/                   # Generated test files
│   ├── alignment/
│   ├── variants/
│   └── sequence/
└── benchmark/                   # Performance testing data
```

---

## 🌍 Real-World Data Files

### Alignment Formats

#### `real_world/1000g_sample.cram` (196 bytes)
- **Source**: 1000 Genomes Project
- **Created**: November 15, 2025
- **Purpose**: CRAM format validation with real 1000 Genomes data
- **Features Tested**: CRAM 3.0/3.1 reader, reference-based compression
- **Integration Tests**: `tests/cram_real_file_test.rs`
- **Expected**: Small sample from chr22 region
- **Note**: Requires `mini_reference.fa` for decompression

### Variant Formats

#### `real_world/synthetic_1000g.vcf.gz` (708 bytes)
- **Source**: 1000 Genomes Project (Phase 3)
- **Created**: November 13, 2025
- **Purpose**: VCF format validation with real variant data
- **Features Tested**: VCF 4.2 parser, gzip decompression, INFO/FORMAT fields
- **Integration Tests**: `tests/real_world_data_integration.rs::test_1000genomes_vcf`
- **Expected**: ~5 variant records, multiple samples
- **Chromosome**: chr21
- **Format**: VCF 4.2, gzip compressed

### Annotation Formats

#### `real_world/encode_peaks.bed.gz` (157 bytes)
- **Source**: ENCODE Project (ChIP-seq peaks)
- **Created**: November 13, 2025
- **Purpose**: BED narrowPeak format validation
- **Features Tested**: BED parser, narrowPeak fields, gzip decompression
- **Integration Tests**: `tests/real_world_data_integration.rs::test_encode_peaks_bed`
- **Expected**: ~10 peak regions
- **Format**: narrowPeak (BED6+4)

#### `real_world/ucsc_genes.bed.gz` (18 MB)
- **Source**: UCSC Genome Browser (knownGene table)
- **Created**: November 13, 2025
- **Purpose**: Large BED12 format validation
- **Features Tested**: BED12 parser, exon blocks, gene structures
- **Integration Tests**: `tests/real_world_data_integration.rs::test_ucsc_genes_bed12`
- **Expected**: ~82,000 gene isoforms
- **Format**: BED12, gzip compressed
- **Genome**: hg38

#### `real_world/annotation/ensembl_chr21.gff3.gz` (533 KB)
- **Source**: Ensembl Release 110 (Human chr21)
- **Created**: November 13, 2025
- **Purpose**: GFF3 format validation with hierarchical features
- **Features Tested**: GFF3 parser, Parent/ID relationships, attributes
- **Integration Tests**: `tests/real_world_data_integration.rs::test_ensembl_gff3`
- **Expected**: ~5,000 features (genes, transcripts, exons)
- **Chromosome**: chr21
- **Format**: GFF3, gzip compressed
- **Genome**: GRCh38

#### `real_world/annotation/ensembl_chr21.gtf.gz` (87 KB)
- **Source**: Ensembl Release 110 (Human chr21)
- **Created**: November 16, 2025
- **Command**: `curl http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz | gunzip | grep "^21\t" | head -5000 | gzip`
- **Purpose**: GTF format validation with real gene annotations
- **Features Tested**: GTF parser, gene/transcript/exon features, attributes
- **Integration Tests**: `tests/gtf_real_world_test.rs`
- **Expected**: 5,000 GTF records from chr21
- **Chromosome**: chr21
- **Format**: GTF (Gene Transfer Format), gzip compressed
- **Genome**: GRCh38

### Assembly Formats

#### `real_world/lambda_phage.gfa` (464 bytes)
- **Source**: Lambda phage genome assembly
- **Created**: November 13, 2025
- **Purpose**: GFA format validation with assembly graph
- **Features Tested**: GFA parser, segments, links, paths
- **Integration Tests**: Referenced in `tests/gfa_integration.rs`
- **Expected**: ~6 segments, multiple links
- **Format**: GFA v1

### Sequence Formats

#### `real_world/sequence/small_sample.fastq.gz` (57 KB)
- **Source**: SRA dataset ERR000589 (1000 Genomes Project)
- **Created**: November 16, 2025
- **Command**: `curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR000/ERR000589/ERR000589_1.fastq.gz | gunzip | head -4000 | gzip`
- **Purpose**: Real-world FASTQ format validation
- **Features Tested**: FASTQ parser, quality scores, sequence IDs
- **Integration Tests**: `tests/fastq_real_world_test.rs`
- **Expected**: 1,000 reads from real Illumina sequencing run
- **Format**: FASTQ, gzip compressed
- **Sample**: Human genome sequencing data

### Alignment Formats

#### `real_world/alignments/minimap2_alignment.paf` (1.7 KB)
- **Source**: Generated by minimap2 v2.30
- **Created**: November 16, 2025
- **Command**: `minimap2 -x map-ont mini_reference.fa query_sequences.fa > minimap2_alignment.paf`
- **Purpose**: PAF format validation with real minimap2 output
- **Features Tested**: PAF parser, alignment coordinates, CIGAR strings, mapping quality
- **Integration Tests**: `tests/paf_real_world_test.rs`
- **Expected**: 18 alignment records
- **Query**: 3 sequences extracted from reference (overlapping regions)
- **Target**: mini_reference.fa (chr1)
- **Format**: PAF (Pairwise mApping Format)

---

## 🧪 Synthetic Test Files

### Alignment Formats

#### `synthetic_100k.bam` (969 KB)
- **Source**: Generated by `generate_test_bam.py`
- **Created**: November 9, 2025
- **Purpose**: BAM format validation, indexing, queries
- **Features Tested**: BAM parser, BGZF decompression, CIGAR, tags
- **Integration Tests**:
  - `tests/bam_writer_integration.rs` (multiple tests)
  - `tests/csi_real_data_test.rs` (CSI index validation)
- **Expected**: 100,000 synthetic reads
- **Reference**: chr22, chr21, chrX
- **Indices**: `synthetic_100k.bam.bai`, `synthetic_100k.bam.csi`

#### `synthetic_100k.bam.bai` (192 bytes)
- **Source**: Generated by samtools index
- **Created**: November 10, 2025
- **Purpose**: BAI index format validation
- **Features Tested**: BAI parser, chunk queries, linear index
- **Integration Tests**: BAM indexed query tests
- **Index Type**: BAI (standard BAM index)

#### `synthetic_100k.bam.csi` (139 bytes)
- **Source**: Generated by `samtools index -c`
- **Created**: November 16, 2025
- **Purpose**: CSI index format validation
- **Features Tested**: CSI parser, configurable binning (14-bit, depth 5)
- **Integration Tests**: `tests/csi_real_data_test.rs`
- **Index Type**: CSI (coordinate-sorted index)
- **Parameters**: min_shift=14, depth=5

#### `test_mini.cram` (105 KB)
- **Source**: Created from synthetic_100k.bam using samtools
- **Created**: November 15, 2025
- **Command**: `samtools view -b synthetic_100k.bam chr22:1-10000 | samtools view -C -T mini_reference.fa -o test_mini.cram -`
- **Purpose**: Small CRAM file for unit testing
- **Features Tested**: CRAM decoder, reference compression, slice parsing
- **Integration Tests**: `tests/cram_real_file_test.rs`
- **Expected**: ~100-200 reads from chr22:1-10000
- **Reference**: `mini_reference.fa`

#### `test.bam` (826 bytes)
- **Source**: Hand-crafted minimal BAM
- **Created**: November 9, 2025
- **Purpose**: Basic BAM parsing, header validation
- **Features Tested**: BAM magic, header parsing, single record
- **Expected**: Minimal valid BAM with 1-2 reads

#### `test.sam` (2.6 KB)
- **Source**: Hand-crafted SAM file
- **Created**: November 9, 2025
- **Purpose**: SAM text format parsing
- **Features Tested**: SAM header, alignment records
- **Expected**: Multiple alignment records with various CIGAR operations

### Variant Formats

#### `test_variants.bcf` (782 bytes)
- **Source**: Generated by bcftools from test_variants.vcf
- **Created**: November 16, 2025
- **Command**: `bcftools view -O b test_variants.vcf -o test_variants.bcf`
- **Purpose**: BCF format validation with typed values
- **Features Tested**: BCF reader, typed values, IDX field parsing, FORMAT fields
- **Integration Tests**: `tests/bcf_integration_test.rs` (5 tests, all passing)
- **Expected**: 5 variant records, 3 samples
- **Chromosomes**: chr1, chr2, chr3
- **Features**:
  - Multi-allelic variants
  - INFO fields: DP, AF, AC
  - FORMAT fields: GT, GQ, DP
  - FILTER: PASS, LowQual

#### `test_variants.vcf` (1.1 KB)
- **Source**: Hand-crafted VCF 4.2
- **Created**: November 16, 2025
- **Purpose**: Source for BCF generation, VCF format testing
- **Features Tested**: VCF 4.2 spec compliance
- **Expected**: 5 variants, 3 samples
- **Format**: VCF 4.2 (text)

### Sequence Formats

#### `mini_reference.fa` (9.9 KB)
- **Source**: Synthetic reference for CRAM testing
- **Created**: November 15, 2025
- **Purpose**: Reference sequences for CRAM decompression
- **Sequences**: chr22, chr21, chrX (small regions)
- **Format**: FASTA
- **Index**: `mini_reference.fa.fai`

#### `mini_reference.fa.fai` (19 bytes)
- **Source**: Generated by samtools faidx
- **Created**: November 15, 2025
- **Purpose**: FAI index format validation
- **Features Tested**: FAI parser, random access
- **Index Type**: FAI (FASTA index)

#### `fasta/test.fa` (219 bytes)
- **Source**: Hand-crafted FASTA
- **Created**: November 14, 2025
- **Purpose**: Basic FASTA parsing
- **Features Tested**: FASTA reader, multi-sequence files
- **Expected**: 2-3 sequences

---

## 📊 Test Coverage by Format

| Format | Synthetic Files | Real-World Files | Integration Tests |
|--------|----------------|------------------|-------------------|
| **BAM** | ✅ synthetic_100k.bam, test.bam | ⚠️ Use real data | ✅ bam_writer_integration.rs |
| **SAM** | ✅ test.sam | N/A (text format) | ✅ Via BAM tests |
| **CRAM** | ✅ test_mini.cram | ✅ 1000g_sample.cram | ✅ cram_real_file_test.rs |
| **BCF** | ✅ test_variants.bcf | ❌ Missing | ✅ bcf_integration_test.rs |
| **VCF** | ✅ test_variants.vcf | ✅ synthetic_1000g.vcf.gz | ✅ real_world_data_integration.rs |
| **BED** | ⚠️ Inline in tests | ✅ ENCODE, UCSC | ✅ real_world_data_integration.rs |
| **GFF3** | ⚠️ Inline in tests | ✅ ensembl_chr21.gff3.gz | ✅ real_world_data_integration.rs |
| **GTF** | ⚠️ Inline in tests | ✅ ensembl_chr21.gtf.gz | ✅ gtf_real_world_test.rs |
| **GFA** | ⚠️ Inline in tests | ✅ lambda_phage.gfa | ✅ gfa_integration.rs |
| **PAF** | ⚠️ Inline in tests | ✅ minimap2_alignment.paf | ✅ paf_real_world_test.rs |
| **FASTQ** | ⚠️ Inline in tests | ✅ small_sample.fastq.gz | ✅ fastq_real_world_test.rs |
| **FASTA** | ✅ test.fa, mini_reference.fa | ❌ Missing | ⚠️ Library tests only |
| **BAI** | ✅ synthetic_100k.bam.bai | N/A (index) | ✅ Via BAM tests |
| **CSI** | ✅ synthetic_100k.bam.csi | N/A (index) | ✅ csi_real_data_test.rs |
| **FAI** | ✅ mini_reference.fa.fai | ❌ Missing | ⚠️ Library tests only |
| **TBI** | ❌ Missing | ❌ Missing | ❌ No tests |

---

## 🚨 Missing Real-World Test Data

### Medium Priority:
4. **BCF** - Real 1000 Genomes BCF file
5. **TBI** - Tabix-indexed VCF/BED file
6. **FASTA** - Real reference sequence (e.g., chr21 from hg38)

---

## 📥 Adding New Test Files

When adding new test files, please:

1. **Document the source**:
   - URL or command used to generate
   - Date created
   - Dataset version/release

2. **Add to appropriate subdirectory**:
   - `real_world/` for public datasets
   - `synthetic/` for generated files
   - Further organize by format type

3. **Create integration test**:
   - Add test to `tests/` directory
   - Validate expected characteristics
   - Test edge cases specific to that dataset

4. **Update this README**:
   - Add file details to appropriate section
   - Document expected content
   - Link to integration tests

5. **Keep files small**:
   - Real-world: Use small regions or subsets
   - Target: < 1 MB per file
   - Exception: Large files for performance testing

---

## 🔗 External Data Sources

### Public Datasets:
- **1000 Genomes**: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/
- **ENCODE**: https://www.encodeproject.org/
- **Ensembl**: http://ftp.ensembl.org/pub/
- **UCSC Genome Browser**: https://hgdownload.soe.ucsc.edu/
- **SRA**: https://www.ncbi.nlm.nih.gov/sra

### Tools for Creating Test Data:
- **samtools**: BAM/CRAM/index generation
- **bcftools**: BCF/VCF manipulation
- **bedtools**: BED file operations
- **minimap2**: PAF alignment generation
- **seqtk**: FASTQ/FASTA sampling

---

## 📝 Notes

- All real-world files should be small (< 1 MB) to keep repository size manageable
- Use `.gz` compression for text formats (VCF, BED, GFF3)
- BGZF compression (.gz) allows random access with TBI/CSI indices
- Reference files (FASTA) needed for CRAM decompression
- Index files (BAI, CSI, FAI, TBI) must match their data files

---

## 🔄 Maintenance

**Last Audit**: November 16, 2025
**Next Audit**: After adding FASTQ/GTF/PAF real-world data

**Cleanup Opportunities**:
- Move inline test data to files
- Add TBI index examples
- Add FASTQ/FASTA real-world samples
- Consolidate duplicate test data
