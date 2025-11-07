# Jupyter Notebooks Plan - biometal Python Tutorials

**Purpose**: Progressive tutorial series teaching biometal's Python API through real bioinformatics workflows
**Target Audience**: Bioinformaticians, students, ML practitioners using Python
**Approach**: Learn by doing with real data and complete examples

---

## Notebook Series Structure

### Beginner Level (Getting Started)

#### 01_getting_started.ipynb
**Goal**: Introduce basic concepts and streaming architecture
**Topics**:
- Installing biometal
- Streaming FASTQ files (constant memory)
- Basic operations: GC content, base counting, quality scores
- Understanding ARM NEON speedup
**Data**: Small test FASTQ file (included)
**Duration**: 15-20 minutes

#### 02_quality_control_pipeline.ipynb
**Goal**: Build complete QC pipeline
**Topics**:
- Quality-based trimming (Trimmomatic-style)
- Length filtering
- Quality-based masking
- QC metrics and visualization
- Complete trim → filter → mask workflow
**Data**: Real E. coli dataset (small, ~1000 reads)
**Duration**: 30-40 minutes

---

### Intermediate Level (Real Use Cases)

#### 03_kmer_analysis.ipynb
**Goal**: K-mer extraction for machine learning
**Topics**:
- K-mer extraction for DNABert/DNABERT-2
- Minimizers for indexing
- K-mer spectrum analysis
- Parallel extraction for large datasets
- Feeding to ML models
**Data**: Bacterial genome sequences
**Duration**: 30-40 minutes
**ML Context**: Preprocessing for transformer models

#### 04_sra_streaming.ipynb
**Goal**: Analyze without downloading
**Topics**:
- Streaming from NCBI SRA
- Network streaming architecture
- Memory efficiency demonstration
- Real-world E. coli analysis (SRR390728)
- Performance tuning
**Data**: SRA accession (streamed, not downloaded)
**Duration**: 30-40 minutes
**Impact**: Analyze TB-scale data on laptop

---

### Advanced Level (Production Workflows)

#### 05_complete_pipeline.ipynb
**Goal**: End-to-end analysis combining everything
**Topics**:
- Stream from SRA
- QC pipeline
- K-mer extraction
- Complexity filtering
- Output to formats (FASTA/FASTQ)
- Performance metrics
**Data**: Complete SRA run
**Duration**: 45-60 minutes
**Output**: Production-ready pipeline

#### 06_performance_comparison.ipynb (Optional - Future)
**Goal**: Benchmark against existing tools
**Topics**:
- biometal vs pure Python
- biometal vs cutadapt (trimming)
- biometal vs Trimmomatic
- Memory usage comparison
- Speed benchmarks
**Data**: Various dataset sizes
**Duration**: 30 minutes
**Purpose**: Validate performance claims

---

## Notebook Features

### Every Notebook Includes:
1. **Learning Objectives** - What you'll learn
2. **Prerequisites** - What you need to know
3. **Setup Section** - Installation and imports
4. **Real Data** - No toy examples
5. **Code Explanations** - Why, not just what
6. **Visualizations** - Charts showing results
7. **Best Practices** - Production tips
8. **Exercise Section** - Try it yourself
9. **What's Next** - Link to next notebook

### Common Elements:
- Markdown cells with context
- Code cells with extensive comments
- Output cells showing real results
- Plots using matplotlib/seaborn
- Memory usage demonstrations
- Timing comparisons
- Error handling examples

---

## Data Strategy

### Test Data (Small, Included in Repo)
- `data/test_reads.fq.gz` - 100 reads for quick testing
- `data/ecoli_sample.fq.gz` - 1000 reads for QC pipeline
- `data/sequences.fasta` - Bacterial sequences for k-mer analysis

### Real Data (Streamed from SRA)
- SRR390728 - E. coli K-12 (~250K reads, 40 MB)
- Streamed directly, not downloaded
- Used in notebooks 04 and 05

### Generate Test Data:
```python
# Script to generate test data
# Will create synthetic FASTQ files with known properties
# Stored in notebooks/generate_test_data.py
```

---

## Technical Approach

### Notebook Format
- **Format**: Jupyter Notebook (.ipynb)
- **Kernel**: Python 3.9+
- **Dependencies**: biometal-rs, matplotlib, seaborn, pandas (for visualization)
- **Style**: Educational, not just code dumps

### Code Style
```python
# Good: Explanatory
# Calculate GC content for quality assessment
# High GC regions may indicate contamination
gc = biometal.gc_content(record.sequence)
if gc > 0.65:  # Typical threshold for bacterial contamination
    print(f"Warning: High GC content {gc:.2%}")

# Bad: Just code
gc = biometal.gc_content(record.sequence)
if gc > 0.65:
    print(f"Warning: {gc}")
```

### Progressive Complexity
- Start simple (single operations)
- Build up (combine operations)
- End complex (full pipelines)

---

## Visualization Strategy

### Plots to Include:
1. **Quality Score Distributions** (before/after trimming)
2. **Read Length Histograms** (before/after filtering)
3. **GC Content Distribution**
4. **K-mer Frequency Spectra**
5. **Memory Usage Over Time** (constant line)
6. **Processing Speed** (reads/second)
7. **QC Metrics Dashboard** (pass/fail rates)

### Visualization Tools:
- matplotlib for basic plots
- seaborn for statistical plots
- pandas for data manipulation
- Custom plotting functions for common patterns

---

## Learning Path

```
01_getting_started.ipynb
    ↓ (Learn basics)
02_quality_control_pipeline.ipynb
    ↓ (Build QC skills)
    ├→ 03_kmer_analysis.ipynb (ML path)
    └→ 04_sra_streaming.ipynb (Network streaming path)
         ↓
05_complete_pipeline.ipynb (Everything together)
```

---

## Distribution Strategy

### In Repository:
- `notebooks/` directory at repo root
- Each notebook is self-contained
- README.md with learning path
- Test data in `notebooks/data/`

### On GitHub:
- Notebooks render natively
- Users can browse without running
- Easy to download and run locally

### Future (Optional):
- **Binder** - Run in browser without installation
- **Google Colab** - Run with GPU access
- **Blog Posts** - Convert to articles
- **Video Walkthroughs** - Screen recordings

---

## Success Metrics

### For Each Notebook:
- [ ] Runs without errors on fresh install
- [ ] Produces expected output
- [ ] All plots render correctly
- [ ] Takes stated duration to complete
- [ ] Clear learning progression
- [ ] Real data (no toy examples)
- [ ] Production best practices shown

### Overall Series:
- [ ] Covers all major biometal features
- [ ] Progressive difficulty
- [ ] Real bioinformatics use cases
- [ ] Complete workflows (not fragments)
- [ ] Copy-paste ready for production

---

## Implementation Checklist

### Phase 1: Core Notebooks (High Priority)
- [ ] Create notebooks directory structure
- [ ] Generate test data
- [ ] 01_getting_started.ipynb
- [ ] 02_quality_control_pipeline.ipynb
- [ ] 03_kmer_analysis.ipynb
- [ ] 04_sra_streaming.ipynb
- [ ] notebooks/README.md

### Phase 2: Advanced Content (Medium Priority)
- [ ] 05_complete_pipeline.ipynb
- [ ] Add more visualizations
- [ ] Add exercise sections

### Phase 3: Polish & Distribution (Future)
- [ ] 06_performance_comparison.ipynb (optional)
- [ ] Set up Binder
- [ ] Convert to blog posts
- [ ] Create video walkthroughs

---

## Estimated Timeline

- **Phase 1**: 4-6 hours (5 notebooks + data)
- **Phase 2**: 2-3 hours (advanced content)
- **Phase 3**: TBD (distribution)

**Total for Phase 1**: 4-6 hours of development
**User Value**: 2-3 hours of learning content

---

## Next Steps

1. Create notebooks directory structure ✅
2. Create test data generation script
3. Implement notebooks 01-04 (Phase 1)
4. Test all notebooks end-to-end
5. Create notebooks/README.md
6. Update main README.md with links
7. Commit and push

---

**Status**: Planning Complete
**Ready to Implement**: Yes
**Priority**: High (post-v1.2.0 strategic priority)
