# biometal Python Tutorials - Jupyter Notebooks

**Learn biometal through real bioinformatics workflows**

Interactive Jupyter notebooks teaching biometal's Python API from basics to production pipelines.

---

## 📚 Learning Path

### Beginner Level (Start Here!)

#### [01_getting_started.ipynb](01_getting_started.ipynb)
**Duration**: 15-20 minutes | **Status**: ✅ Complete

Learn the basics of biometal:
- Installing and importing biometal
- Streaming FASTQ files with constant memory
- Calculating GC content and base composition
- Analyzing quality scores
- Understanding ARM NEON performance benefits

**Perfect for**: First-time users, understanding streaming architecture

---

### Intermediate Level (Coming Soon)

#### 02_quality_control_pipeline.ipynb
**Duration**: 30-40 minutes | **Status**: 🚧 In Development

Build complete QC pipelines:
- Quality-based trimming (Trimmomatic-style sliding window)
- Length filtering
- Quality-based masking
- QC metrics and visualization
- Complete trim → filter → mask workflow

**Perfect for**: Pre-processing reads for downstream analysis

#### 03_kmer_analysis.ipynb
**Duration**: 30-40 minutes | **Status**: 🚧 Planned

K-mer extraction for machine learning:
- K-mer extraction for DNABert/DNABERT-2
- Minimizers for indexing (minimap2-style)
- K-mer spectrum analysis
- Parallel extraction for large datasets
- Feeding to ML models

**Perfect for**: ML practitioners, transformer model preprocessing

#### 04_sra_streaming.ipynb
**Duration**: 30-40 minutes | **Status**: 🚧 Planned

Analyze without downloading:
- Streaming directly from NCBI SRA
- Network streaming architecture
- Memory efficiency demonstration
- Real-world E. coli analysis (SRR390728)
- Performance tuning

**Perfect for**: Working with public datasets, cloud analysis

---

### Advanced Level (Future)

#### 05_complete_pipeline.ipynb
**Status**: 📋 Planned

End-to-end production workflow:
- Stream from SRA
- Complete QC pipeline
- K-mer extraction
- Complexity filtering
- Output to formats
- Performance metrics

**Perfect for**: Production bioinformatics pipelines

---

## 🚀 Quick Start

### Installation

```bash
# Install biometal
pip install biometal-rs

# Install Jupyter
pip install jupyter matplotlib seaborn pandas

# Clone repo (for notebooks)
git clone https://github.com/shandley/biometal
cd biometal/notebooks

# Start Jupyter
jupyter notebook
```

### Run in Browser

Open `01_getting_started.ipynb` and run all cells (Cell → Run All).

---

## 📊 What You'll Learn

### Core Concepts
✅ **Streaming Architecture**: Analyze 5TB datasets with 5 MB memory
✅ **ARM NEON Performance**: 16-25× speedup on Apple Silicon
✅ **Network Streaming**: Analyze without downloading
✅ **Production Quality**: Real workflows, not toy examples

### Practical Skills
✅ Quality control pipelines (trim, filter, mask)
✅ K-mer extraction for ML (DNABert preprocessing)
✅ SRA streaming (access public data instantly)
✅ Memory-efficient processing (constant memory)

### biom

etal Features (v1.2.0)
✅ **Core Operations**: GC content, base counting, quality scores (v1.0.0)
✅ **K-mer Operations**: Extraction, minimizers, spectrum (v1.1.0)
✅ **Phase 4 Operations**: Trimming, masking, sequence manipulation (v1.2.0)

---

## 💻 Requirements

### Software
- **Python**: 3.9+ (tested on 3.9-3.14)
- **biometal-rs**: 1.2.0+ (`pip install biometal-rs`)
- **Jupyter**: `pip install jupyter`
- **Visualization** (optional): `pip install matplotlib seaborn pandas`

### Hardware
- **Minimum**: Any modern laptop (x86_64 or ARM)
- **Recommended**: Apple Silicon (M1/M2/M3/M4) for NEON acceleration
- **Memory**: 8 GB RAM (biometal uses ~5 MB, but Jupyter needs more)

### Data
- Test data included in notebooks (small synthetic files)
- Real data streamed from NCBI SRA (no download needed!)

---

## 🎯 Learning Outcomes

After completing these notebooks, you will be able to:

1. **Stream and analyze** large FASTQ/FASTA files with constant memory
2. **Build QC pipelines** for read preprocessing (trim → filter → mask)
3. **Extract k-mers** for machine learning (DNABert, transformers)
4. **Stream from SRA** to analyze public data without downloading
5. **Write production pipelines** using biometal's streaming API
6. **Optimize performance** on Apple Silicon (ARM NEON)

---

## 📖 Tutorial Format

Each notebook follows this structure:

1. **Learning Objectives** - What you'll learn
2. **Prerequisites** - What you need to know
3. **Concepts** - Background and motivation
4. **Hands-On Code** - Real examples with explanations
5. **Visualizations** - Charts and plots showing results
6. **Best Practices** - Production tips
7. **Exercises** - Try it yourself
8. **What's Next** - Link to next tutorial

---

## 🔬 Real Data Examples

### Included Test Data
- Synthetic FASTQ files (generated in notebooks)
- Small E. coli samples (~1000 reads)
- Bacterial sequences for k-mer analysis

### Streamed Real Data
- **SRR390728**: E. coli K-12 (~250K reads, 40 MB)
  - Streamed directly from NCBI SRA
  - No download required!
  - Used in notebooks 04 and 05

---

## 🆘 Getting Help

### Common Issues

**Q: Import error `ModuleNotFoundError: No module named 'biometal'`**
A: Install with `pip install biometal-rs` (note the `-rs`)

**Q: Notebook says "Kernel starting, please wait"**
A: Give it 10-30 seconds on first start

**Q: Can't find test data files**
A: Test data is generated in the notebooks (run all cells)

**Q: Performance seems slow**
A: Check your platform (run platform check cell). NEON acceleration only on ARM.

### Support

- **GitHub Issues**: https://github.com/shandley/biometal/issues
- **Discussions**: https://github.com/shandley/biometal/discussions
- **Documentation**: https://docs.rs/biometal

---

## 🤝 Contributing

Found an issue or want to improve a notebook?

1. Open an issue: https://github.com/shandley/biometal/issues
2. Suggest improvements
3. Share your own examples!

---

## 📜 License

These tutorials are part of biometal and licensed under MIT OR Apache-2.0 (your choice).

---

## 🙏 Acknowledgments

- **Evidence Base**: 1,357 experiments from apple-silicon-bio-bench
- **Community**: Feedback from bioinformaticians worldwide
- **Platform**: Optimized for Apple Silicon (M1/M2/M3/M4)

---

**biometal v1.2.0** - ARM-native bioinformatics with streaming architecture

**Mission**: Democratizing bioinformatics compute for LMIC researchers, small labs, students, and field researchers.
