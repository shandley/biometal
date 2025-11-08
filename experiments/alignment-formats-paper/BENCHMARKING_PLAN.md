# Comprehensive Benchmarking Plan for Paper

**Purpose**: Generate publication-quality performance data for dual-format paper
**Timeline**: Weeks 19-20 (after BAM + CAF implementation complete)
**Deliverables**: Figures 1-4, Tables 1-3, Supplementary data

---

## Benchmark Infrastructure

### Hardware Platforms

**Primary (ARM)**:
- Mac M2 (Apple Silicon): Personal development machine
- AWS Graviton3: Cloud ARM validation
- Specifications to report: CPU model, cores, RAM, storage type

**Comparison (x86)**:
- Intel Xeon: Baseline comparison
- Validate portable fallback performance

### Software Environment

**Record for reproducibility**:
```toml
[benchmark-environment]
rust_version = "1.74+"
biometal_version = "X.X.X"
noodles_version = "0.80+"
criterion_version = "0.5+"
rayon_version = "1.8+"

[hardware.mac-m2]
cpu = "Apple M2"
cores = 8
ram = "16 GB"
storage = "NVMe SSD"

[hardware.graviton3]
cpu = "AWS Graviton3"
cores = 64
ram = "128 GB"
storage = "EBS gp3"
```

---

## Dataset Selection

### Diversity Requirements

**Technologies**:
- Illumina (short reads, 100-150 bp)
- PacBio HiFi (long reads, 10-20 Kbp)
- ONT Nanopore (ultra-long, 5-50 Kbp)

**Sizes**:
- Small: 10K records (~1 MB BAM)
- Medium: 100K records (~50 MB BAM)
- Large: 1M records (~500 MB BAM)
- Very large: 10M records (~5 GB BAM)

### Specific Datasets (Public)

**Dataset 1: 1000 Genomes (Illumina)**
```
Name: HG00096 (GBR population)
Source: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/
Technology: Illumina HiSeq
Read length: 100 bp paired-end
Coverage: Low coverage (4-6×)
Use: Baseline benchmark, most common use case
```

**Dataset 2: GIAB HG002 (PacBio)**
```
Name: HG002 (Ashkenazi trio)
Source: https://www.nist.gov/programs-projects/genome-bottle
Technology: PacBio HiFi
Read length: 10-20 Kbp
Coverage: 30×
Use: Long read performance
```

**Dataset 3: GIAB HG002 (ONT)**
```
Name: HG002 (same sample)
Source: GIAB consortium
Technology: ONT PromethION
Read length: 5-50 Kbp (ultra-long)
Coverage: 50×
Use: Ultra-long read stress test
```

**Dataset 4: Synthetic (Edge Cases)**
```
Name: Generated test data
Properties:
  - Extreme quality variation
  - Many unmapped reads
  - High duplicate rate
  - Unusual CIGAR operations
Use: Edge case validation
```

---

## Benchmark Experiments

### Experiment 1: Parse Throughput

**Metric**: Records/second, MB/second

**Variants**:
- noodles BAM (baseline)
- biometal BAM (ARM-optimized)
- biometal CAF

**Operations**:
```rust
// Minimal parsing (just iterate)
let mut count = 0;
let mut record = Record::default();
while reader.read_record(&mut record)? != 0 {
    count += 1;
}
```

**Parameters**:
- Dataset sizes: 10K, 100K, 1M, 10M records
- Repetitions: N=30 (criterion default)
- Warmup: 5 iterations

**Expected results** (Table 2, Figure 2A):
```
noodles:      39.1 Melem/s (baseline)
biometal BAM:  87.3 Melem/s (2.2× faster)
biometal CAF: 385.2 Melem/s (9.8× faster)
```

### Experiment 2: Quality Filtering

**Metric**: Time to filter records by mean quality ≥ threshold

**Operation**:
```rust
// Filter records by quality score
let passing: Vec<Record> = records
    .filter(|r| mean_quality(&r.quality) >= 30.0)
    .collect();
```

**Parameters**:
- Thresholds: Q20, Q30, Q40
- Dataset: 100K Illumina records
- Repetitions: N=30

**Expected results** (Table 2, Figure 2B):
```
BAM:  2.01 seconds
CAF:  0.09 seconds (22.3× faster)
```

### Experiment 3: Base Counting

**Metric**: Time to count A/C/G/T/N across all sequences

**Operation**:
```rust
// Count base composition
let mut counts = [0u32; 5];  // A, C, G, T, N
for record in records {
    for &base in record.sequence() {
        match base {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            // ...
        }
    }
}
```

**Expected results**:
```
BAM:  1.48 seconds
CAF:  0.07 seconds (21.1× faster, NEON parallel count)
```

### Experiment 4: MAPQ Filtering

**Metric**: Time to filter by mapping quality

**Operation**:
```rust
let high_quality: Vec<Record> = records
    .filter(|r| r.mapq >= 30)
    .collect();
```

**Expected results**:
```
BAM:  0.51 seconds
CAF:  0.04 seconds (12.8× faster, NEON parallel compare)
```

### Experiment 5: Region Extraction

**Metric**: Time to extract records in genomic region

**Operation**:
```rust
// Extract chr1:1000000-2000000
let region_records: Vec<Record> = records
    .filter(|r| {
        r.ref_id == chr1_id &&
        r.position >= 1_000_000 &&
        r.position < 2_000_000
    })
    .collect();
```

**Expected results**:
```
BAM:  0.89 seconds (sequential scan)
CAF:  0.12 seconds (block-level index + NEON filter)
```

### Experiment 6: Real-World Pipeline

**Quality Control Workflow**:
```rust
// Realistic bioinformatics workflow
fn qc_pipeline(input: &Path) -> QcReport {
    let records = read_alignment(input)?;

    // Filter Q30
    let high_quality = records.filter(|r| mean_quality(r) >= 30.0);

    // Count bases
    let base_comp = count_bases(&high_quality);

    // Calculate GC content
    let gc_content = (base_comp.G + base_comp.C) as f64 / total as f64;

    // Extract high-quality regions
    let regions = extract_high_quality_regions(&high_quality);

    QcReport { ... }
}
```

**Benchmark**: End-to-end time for complete workflow

**Expected results** (Figure 4):
```
Input: 1M Illumina records
BAM workflow:  18.2 seconds
CAF workflow:   2.1 seconds (8.7× faster)
```

### Experiment 7: Scaling Analysis

**Metric**: Throughput vs dataset size

**Datasets**: 10K, 50K, 100K, 500K, 1M, 5M, 10M records

**Operations**: Parse, filter Q30, count bases

**Expected** (Figure 2C):
- CAF: Linear scaling (constant per-record cost)
- BAM: Sublinear (decompression overhead)

### Experiment 8: Memory Usage

**Metric**: Peak RSS during processing

**Operation**: Process entire file with streaming

**Expected** (Supplementary):
```
BAM:  ~5 MB (constant, streaming)
CAF:  ~5 MB (constant, streaming)
noodles (iterator):  ~8 MB (slight overhead)
```

### Experiment 9: Storage Analysis

**Metric**: Compressed file size

**Datasets**: All datasets (Illumina, PacBio, ONT)

**Expected results** (Table 3):
```
Format          | Illumina | PacBio | ONT
BAM (gzip)      | 47.3 MB  | 89.2 MB | 156.3 MB
CAF (zstd/lz4)  | 71.8 MB  | 132.1 MB | 221.7 MB
Ratio           | 1.52×    | 1.48×   | 1.42×
```

### Experiment 10: Conversion Performance

**Metric**: Time to convert BAM ↔ CAF

**Operation**:
```rust
// BAM → CAF
convert_bam_to_caf(input_bam, output_caf)?;

// CAF → BAM
convert_caf_to_bam(input_caf, output_bam)?;

// Round-trip validation
assert_eq!(original_bam, roundtrip_bam);
```

**Expected**:
```
BAM → CAF: ~1.5× parse time (write columnar layout)
CAF → BAM: ~1.2× parse time (write row-oriented)
```

---

## Statistical Analysis

### Methodology

**Criterion.rs configuration**:
```rust
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

fn benchmark_parse(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse");
    group.sample_size(30);  // N=30
    group.measurement_time(Duration::from_secs(60));  // 1 min per benchmark

    for size in [10_000, 100_000, 1_000_000].iter() {
        group.bench_with_input(BenchmarkId::new("noodles", size), size, |b, &size| {
            b.iter(|| parse_noodles(size))
        });

        group.bench_with_input(BenchmarkId::new("biometal_bam", size), size, |b, &size| {
            b.iter(|| parse_biometal_bam(size))
        });

        group.bench_with_input(BenchmarkId::new("biometal_caf", size), size, |b, &size| {
            b.iter(|| parse_biometal_caf(size))
        });
    }

    group.finish();
}
```

**Statistical reporting**:
- Mean ± standard deviation
- Median (robust to outliers)
- 95% confidence intervals
- Significance testing (t-test, p < 0.05)

### Reproducibility

**Documentation**:
```markdown
# Reproducing Benchmarks

## Prerequisites
- Rust 1.74+
- ARM64 hardware (Mac M1+ or AWS Graviton)
- 16 GB RAM minimum
- 50 GB free disk space

## Setup
```bash
git clone https://github.com/[user]/biometal
cd biometal
cargo build --release

## Download datasets
bash scripts/download_benchmark_data.sh

## Run benchmarks
cargo bench --bench parse_throughput
cargo bench --bench quality_filter
cargo bench --bench base_counting
# ... etc

## Generate figures
python scripts/plot_results.py
```

---

## Data Collection & Organization

### Directory Structure

```
experiments/alignment-formats-paper/
├── benchmarks/
│   ├── raw_data/               # Criterion output (JSON)
│   │   ├── parse_throughput/
│   │   ├── quality_filter/
│   │   └── ...
│   ├── processed/              # Processed for figures/tables
│   │   ├── table_2.csv
│   │   ├── figure_2a.csv
│   │   └── ...
│   ├── figures/                # Generated figures
│   │   ├── figure_1.pdf
│   │   ├── figure_2.pdf
│   │   └── ...
│   └── scripts/                # Analysis scripts
│       ├── plot_figure_2.py
│       ├── generate_table_2.R
│       └── ...
```

### Data Formats

**Raw benchmark output** (Criterion JSON):
```json
{
  "group_id": "parse_throughput",
  "function_id": "biometal_bam",
  "value": {
    "throughput": 87300000.0,
    "unit": "records/sec"
  },
  "mean": {
    "point_estimate": 87326415.2,
    "confidence_interval": {
      "lower_bound": 86894321.1,
      "upper_bound": 87758509.3
    }
  }
}
```

**Processed for tables** (CSV):
```csv
operation,format,mean_time,std_dev,speedup
parse,noodles,2.56,0.08,1.0
parse,biometal_bam,1.15,0.04,2.2
parse,biometal_caf,0.26,0.01,9.8
```

---

## Figure Generation

### Figure 1: Format Architecture

**Tool**: Manual diagram (Inkscape, OmniGraffle, or similar)

**Panels**:
- A: BAM row-oriented layout (schema diagram)
- B: CAF columnar layout (schema diagram)
- C: NEON vectorization (before/after illustration)

### Figure 2: Performance Benchmarks

**Tool**: Python (matplotlib/seaborn) or R (ggplot2)

**Script**: `scripts/plot_figure_2.py`

```python
import matplotlib.pyplot as plt
import pandas as pd

# Panel A: Parse throughput
df = pd.read_csv('processed/parse_throughput.csv')
ax = df.plot(x='dataset_size', y='throughput', kind='line')
ax.set_xlabel('Dataset Size (records)')
ax.set_ylabel('Throughput (Melem/s)')

# Panel B: Operation speedups (bar chart)
df = pd.read_csv('processed/operation_speedups.csv')
df.plot(x='operation', y='speedup', kind='bar')

# Panel C: Scaling analysis
# ...
```

### Figure 3: Storage vs Performance

**Tool**: Python scatter plot

**Data**: Processed CSV with columns: storage_ratio, speedup, operation

### Figure 4: Real-World Workflow

**Tool**: Timeline diagram (matplotlib or manual)

**Data**: End-to-end pipeline times

---

## Table Generation

### Table 1: Benchmark Datasets

**Auto-generated** from metadata:

```python
# scripts/generate_table_1.py
datasets = [
    {"name": "HG00096", "source": "1000G", "records": "100K", "tech": "Illumina"},
    {"name": "HG002 PacBio", "source": "GIAB", "records": "50K", "tech": "PacBio"},
    # ...
]

df = pd.DataFrame(datasets)
df.to_latex('tables/table_1.tex', index=False)
```

### Table 2: Performance Comparison

**Auto-generated** from benchmark results:

```python
# scripts/generate_table_2.py
results = load_benchmark_results('raw_data/')
summary = compute_summary_stats(results)
summary.to_latex('tables/table_2.tex', float_format='%.2f')
```

---

## Quality Control

### Validation Checklist

- [ ] All benchmarks run N=30 repetitions
- [ ] Statistical significance tested (p < 0.05)
- [ ] Outliers identified and documented
- [ ] Reproducibility verified (re-run subset)
- [ ] Data backed up (multiple locations)
- [ ] Figures meet journal requirements (resolution, format)
- [ ] Tables formatted correctly (LaTeX)
- [ ] Supplementary data organized

### Peer Review Preparation

**Anticipate questions**:
1. "Why not compare to samtools?" → Include samtools in baseline
2. "Is ARM speedup portable?" → Include x86 comparison
3. "What about other operations?" → Supplementary comprehensive tests
4. "Can we reproduce this?" → Provide complete scripts + data

---

## Timeline (Weeks 19-20)

### Week 19: Data Collection

**Days 1-2**: Infrastructure setup
- Download datasets
- Configure benchmark harness
- Validate environment

**Days 3-5**: Core benchmarks
- Experiments 1-6 (parse, filter, count, etc.)
- N=30 repetitions each
- Preliminary analysis

**Days 6-7**: Extended benchmarks
- Experiments 7-10 (scaling, memory, storage, conversion)
- Validate results
- Identify any issues

### Week 20: Analysis & Figures

**Days 1-2**: Data processing
- Process raw criterion output
- Statistical analysis
- Generate summary tables

**Days 3-4**: Figure generation
- Create all figures (1-4)
- Iterate on design
- Ensure publication quality

**Days 5**: Supplementary material
- Extended benchmarks
- Additional datasets
- Reproducibility documentation

**Days 6-7**: Quality control
- Re-run subset for validation
- Peer review (internal)
- Finalize data package

---

## Deliverables

**For paper**:
- [ ] Figure 1 (format architecture)
- [ ] Figure 2 (performance benchmarks)
- [ ] Figure 3 (storage vs performance)
- [ ] Figure 4 (real-world workflow)
- [ ] Table 1 (datasets)
- [ ] Table 2 (performance comparison)
- [ ] Table 3 (storage analysis)

**Supplementary**:
- [ ] All raw data (Criterion JSON)
- [ ] Processed data (CSV)
- [ ] Analysis scripts (Python/R)
- [ ] Reproducibility guide
- [ ] Extended benchmarks

**Archive** (Zenodo):
- [ ] Complete dataset
- [ ] DOI for citation
- [ ] README with instructions

---

**Status**: Plan complete, ready for execution after BAM + CAF implementation
**Next**: Complete implementations (Weeks 1-18), then execute benchmark plan (Weeks 19-20)
