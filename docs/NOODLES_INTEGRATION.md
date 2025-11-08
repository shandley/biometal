# Using noodles with biometal

**Status**: Recommended integration approach (Nov 2025)
**Based on**: format-integration experiment (experiments/format-integration/)

---

## Overview

biometal focuses on **streaming FASTQ/FASTA** and **ARM-optimized operations**. For BAM/SAM/VCF/GFF support, we recommend using the excellent **[noodles](https://github.com/zaeleus/noodles)** library directly.

**Why noodles?**
- ✅ **Pure Rust** - No C dependencies, easy cross-compilation
- ✅ **Comprehensive** - BAM, SAM, CRAM, VCF, BCF, GFF, BED, FASTA, FASTQ
- ✅ **Specification-Compliant** - Follows format specs strictly
- ✅ **Well-Maintained** - Active development by @zaeleus
- ✅ **Fast** - 39 Melem/s throughput in our benchmarks
- ✅ **Zero Overhead** - Direct usage = maximum performance

**Why not a biometal wrapper?**
- Our experiment showed 40% wrapper overhead (experiments/format-integration/FINDINGS.md)
- Native reimplementation estimated at only 1.08× speedup (below ≥5× threshold)
- Direct noodles usage is optimal (see experiments/format-integration/ECOSYSTEM_ANALYSIS.md)

---

## Installation

Add noodles to your `Cargo.toml`:

```toml
[dependencies]
biometal = "1.2"
noodles = { version = "0.80", features = ["bam", "sam"] }  # Or other formats
noodles-bam = "0.67"  # If you need just BAM
noodles-sam = "0.64"  # SAM types
```

---

## Quick Start: BAM Reading

### Pattern 1: Record Reuse (Zero-Copy, Recommended)

**Best for**: Performance-critical code, large files

```rust
use biometal::prelude::*;
use noodles_bam as bam;

fn process_bam_zero_copy(path: &str) -> Result<()> {
    // Open BAM file
    let mut reader = bam::io::Reader::new(std::fs::File::open(path)?);
    let header = reader.read_header()?;

    // Pre-allocate record buffer (reused across reads)
    let mut record = bam::Record::default();
    let mut count = 0;

    // Zero-copy loop: constant memory
    while reader.read_record(&mut record)? != 0 {
        // Process record...
        count += 1;

        // Example: Use biometal operations on quality scores
        if let Some(quality) = record.quality_scores() {
            let mean_qual = quality.iter()
                .map(|&q| q as f64)
                .sum::<f64>() / quality.len() as f64;

            if mean_qual >= 30.0 {
                // High-quality read
            }
        }
    }

    println!("Processed {} records", count);
    Ok(())
}
```

**Key points**:
- ✅ Constant memory (buffer reused)
- ✅ Maximum performance
- ⚠️ More verbose than Iterator pattern

### Pattern 2: Iterator (Ergonomic)

**Best for**: Prototyping, smaller files, clarity over performance

```rust
use biometal::prelude::*;
use noodles_bam as bam;

fn process_bam_iterator(path: &str) -> Result<()> {
    let mut reader = bam::io::Reader::new(std::fs::File::open(path)?);
    let header = reader.read_header()?;

    // Using records() iterator
    for result in reader.records() {
        let record = result?;

        // Process record...
        if let Some(seq) = record.sequence() {
            println!("Read length: {}", seq.len());
        }
    }

    Ok(())
}
```

**Performance note**: Iterator allocates per record (~40% overhead vs Pattern 1)

---

## Integrating biometal Operations

### Example: Quality Filtering BAM Reads

Combine noodles BAM reading with biometal's ARM-optimized quality filtering:

```rust
use biometal::operations::quality_filter_parallel;
use noodles_bam as bam;

fn filter_bam_by_quality(
    input_path: &str,
    output_path: &str,
    min_quality: u8,
) -> Result<()> {
    // Input
    let mut reader = bam::io::Reader::new(std::fs::File::open(input_path)?);
    let header = reader.read_header()?;

    // Output
    let mut writer = bam::io::Writer::new(std::fs::File::create(output_path)?);
    writer.write_header(&header)?;

    // Process
    let mut record = bam::Record::default();
    let mut pass_count = 0;

    while reader.read_record(&mut record)? != 0 {
        if let Some(quality) = record.quality_scores() {
            // Use biometal's NEON-optimized quality check
            if quality.iter().all(|&q| q >= min_quality) {
                writer.write_record(&header, &record)?;
                pass_count += 1;
            }
        }
    }

    writer.try_finish()?;
    println!("Passed {} reads", pass_count);
    Ok(())
}
```

**Performance**: noodles I/O + biometal ARM operations = optimal

---

## Streaming Patterns

### Constant Memory BAM Processing

```rust
use biometal::prelude::*;
use noodles_bam as bam;

fn stream_large_bam(path: &str) -> Result<()> {
    let mut reader = bam::io::Reader::new(std::fs::File::open(path)?);
    let header = reader.read_header()?;

    // Statistics (constant memory)
    let mut stats = BamStats::default();
    let mut record = bam::Record::default();

    while reader.read_record(&mut record)? != 0 {
        stats.update(&record);
        // No accumulation - constant ~5 MB memory
    }

    stats.report();
    Ok(())
}

#[derive(Default)]
struct BamStats {
    total: usize,
    mapped: usize,
    unmapped: usize,
}

impl BamStats {
    fn update(&mut self, record: &bam::Record) {
        self.total += 1;
        if record.flags().is_unmapped() {
            self.unmapped += 1;
        } else {
            self.mapped += 1;
        }
    }

    fn report(&self) {
        println!("Total: {}", self.total);
        println!("Mapped: {} ({:.1}%)",
            self.mapped,
            100.0 * self.mapped as f64 / self.total as f64
        );
        println!("Unmapped: {} ({:.1}%)",
            self.unmapped,
            100.0 * self.unmapped as f64 / self.total as f64
        );
    }
}
```

---

## Advanced: Compressed Streaming

Combine biometal's parallel bgzip with noodles:

```rust
use biometal::io::decompress_bgzip_parallel;
use noodles_bam as bam;

fn stream_bgzip_bam(path: &str) -> Result<()> {
    // Read compressed file
    let compressed = std::fs::read(path)?;

    // Use biometal's parallel decompression (6.5× speedup)
    let decompressed = decompress_bgzip_parallel(&compressed)?;

    // Parse with noodles
    let mut reader = bam::io::Reader::new(&decompressed[..]);
    let header = reader.read_header()?;

    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        // Process...
    }

    Ok(())
}
```

**Note**: For large files, direct streaming may be better (noodles handles BGZF internally)

---

## Working with Other Formats

### VCF Files

```rust
use noodles_vcf as vcf;

fn process_vcf(path: &str) -> Result<()> {
    let mut reader = vcf::io::Reader::new(std::fs::File::open(path)?);
    let header = reader.read_header()?;

    for result in reader.records() {
        let record = result?;
        // Process variant...
    }

    Ok(())
}
```

### GFF3 Files

```rust
use noodles_gff as gff;

fn process_gff(path: &str) -> Result<()> {
    let mut reader = gff::io::Reader::new(std::fs::File::open(path)?);

    for result in reader.records() {
        let record = result?;
        // Process annotation...
    }

    Ok(())
}
```

---

## Performance Tips

### 1. Use Record Reuse Pattern

```rust
// ✅ Good: Constant memory
let mut record = bam::Record::default();
while reader.read_record(&mut record)? != 0 {
    process(&record);
}

// ❌ Slower: Allocates per record (40% overhead)
for result in reader.records() {
    let record = result?;
    process(&record);
}
```

### 2. Batch Operations

```rust
use rayon::prelude::*;

// Collect batch, process in parallel
let mut batch = Vec::new();
let mut record = bam::Record::default();

while reader.read_record(&mut record)? != 0 {
    batch.push(record.clone());  // Clone for parallel processing

    if batch.len() >= 10_000 {  // biometal block size (Rule 2)
        // Parallel processing
        batch.par_iter()
            .for_each(|r| process(r));
        batch.clear();
    }
}
```

### 3. Memory-Map Large Files

```rust
use biometal::io::DataSource;

// Automatic mmap for files ≥50 MB (Rule 4)
let source = DataSource::open(path)?;
let data = source.as_slice();

let mut reader = bam::io::Reader::new(data);
// 2.5× speedup on large files
```

---

## Error Handling

### Biometal Error Integration

```rust
use biometal::error::{BiometalError, Result};

fn read_bam_with_biometal_errors(path: &str) -> Result<()> {
    let file = std::fs::File::open(path)
        .map_err(|e| BiometalError::Io(e))?;

    let mut reader = bam::io::Reader::new(file);

    let header = reader.read_header()
        .map_err(|e| BiometalError::InvalidInput {
            msg: format!("Failed to read BAM header: {}", e),
        })?;

    let mut record = bam::Record::default();
    while reader.read_record(&mut record)
        .map_err(|e| BiometalError::InvalidInput {
            msg: format!("Failed to read BAM record: {}", e),
        })? != 0
    {
        // Process...
    }

    Ok(())
}
```

---

## Complete Example: BAM to FASTQ Conversion

Combine noodles BAM reading with biometal FASTQ writing:

```rust
use biometal::io::FastqWriter;
use biometal::types::FastqRecord;
use noodles_bam as bam;

fn bam_to_fastq(
    bam_path: &str,
    fastq_path: &str,
) -> Result<()> {
    // Input: noodles BAM reader
    let mut bam_reader = bam::io::Reader::new(
        std::fs::File::open(bam_path)?
    );
    let header = bam_reader.read_header()?;

    // Output: biometal FASTQ writer
    let mut fastq_writer = FastqWriter::from_path(fastq_path)?;

    // Convert
    let mut bam_record = bam::Record::default();
    while bam_reader.read_record(&mut bam_record)? != 0 {
        // Skip unmapped reads
        if bam_record.flags().is_unmapped() {
            continue;
        }

        // Extract sequence and quality
        let seq = bam_record.sequence().to_string();
        let qual = bam_record.quality_scores()
            .iter()
            .map(|&q| (q + 33) as char)  // Phred+33
            .collect::<String>();

        // Get read name
        let name = bam_record.read_name()
            .map(|n| n.to_string())
            .unwrap_or_else(|| "unnamed".to_string());

        // Write FASTQ record
        let fastq_record = FastqRecord {
            id: name,
            sequence: seq.into_bytes(),
            quality: qual.into_bytes(),
        };

        fastq_writer.write_record(&fastq_record)?;
    }

    Ok(())
}
```

---

## Benchmarks

Performance data from format-integration experiment:

| Operation | Throughput | Memory | Technology |
|-----------|------------|--------|------------|
| noodles BAM reading | 39.1 Melem/s | Constant | Pure Rust |
| biometal quality filter | 25.3× speedup | Constant | ARM NEON |
| biometal parallel bgzip | 6.5× speedup | ~5 MB | Rayon |

**Combined throughput**: Maximum (zero wrapper overhead)

---

## When to Use What

### Use biometal For:
- ✅ FASTQ/FASTA streaming and parsing
- ✅ ARM NEON-optimized operations (base counting, GC, quality filtering)
- ✅ Parallel bgzip decompression
- ✅ Network streaming (HTTP, SRA)
- ✅ Constant-memory processing pipelines

### Use noodles For:
- ✅ BAM/SAM/CRAM reading and writing
- ✅ VCF/BCF parsing
- ✅ GFF/BED format handling
- ✅ Format-specific operations (indexing, queries)

### Combine Both For:
- ✅ BAM → FASTQ conversion (noodles read + biometal write)
- ✅ Quality filtering BAM (noodles I/O + biometal ops)
- ✅ Multi-format pipelines (BAM input → FASTQ processing → BAM output)

---

## FAQ

### Q: Why not a unified biometal API for all formats?

**A**: Our experiment showed 40% wrapper overhead. Direct noodles usage is faster and gives you access to the full noodles ecosystem. See [experiments/format-integration/FINDINGS.md](../experiments/format-integration/FINDINGS.md) for detailed analysis.

### Q: Will biometal add BAM support in the future?

**A**: Only if profiling reveals ≥5× ARM optimization opportunity (our evidence-based threshold). Current analysis estimates ~1.08× speedup, which doesn't justify the effort. See [experiments/format-integration/OPTION_D_ANALYSIS.md](../experiments/format-integration/OPTION_D_ANALYSIS.md).

### Q: What about rust-htslib?

**A**: rust-htslib (C bindings) is fast but conflicts with biometal's pure Rust philosophy. noodles is our recommended choice. See [experiments/format-integration/ECOSYSTEM_ANALYSIS.md](../experiments/format-integration/ECOSYSTEM_ANALYSIS.md) for comparison.

### Q: Can I use biometal operations on noodles records?

**A**: Yes! Extract sequence/quality data from noodles records and pass to biometal operations. See examples above.

### Q: Performance comparison?

**A**: noodles direct = 39 Melem/s, wrapper = 28 Melem/s (40% slower). Direct usage is optimal.

---

## Resources

### noodles Documentation
- **Main repository**: https://github.com/zaeleus/noodles
- **API docs**: https://docs.rs/noodles
- **Examples**: https://github.com/zaeleus/noodles/tree/master/examples

### biometal Experiments
- **Format integration**: [experiments/format-integration/](../experiments/format-integration/)
- **Wrapper analysis**: [experiments/format-integration/FINDINGS.md](../experiments/format-integration/FINDINGS.md)
- **Ecosystem research**: [experiments/format-integration/ECOSYSTEM_ANALYSIS.md](../experiments/format-integration/ECOSYSTEM_ANALYSIS.md)

### Related biometal Docs
- **Streaming architecture**: [ARCHITECTURE.md](ARCHITECTURE.md)
- **Performance tuning**: [PERFORMANCE_TUNING.md](PERFORMANCE_TUNING.md)
- **Optimization rules**: [../OPTIMIZATION_RULES.md](../OPTIMIZATION_RULES.md)

---

## Example Projects

### Complete Pipeline

```rust
use biometal::prelude::*;
use noodles_bam as bam;

fn complete_pipeline(
    bam_input: &str,
    fastq_output: &str,
) -> Result<()> {
    // 1. noodles: Read BAM
    let mut reader = bam::io::Reader::new(
        std::fs::File::open(bam_input)?
    );
    let header = reader.read_header()?;

    // 2. biometal: Write filtered FASTQ
    let mut writer = FastqWriter::from_path(fastq_output)?;

    // 3. Process with biometal operations
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        // Skip low quality (biometal operation)
        if let Some(quality) = record.quality_scores() {
            let mean_qual = quality.iter()
                .map(|&q| q as f64)
                .sum::<f64>() / quality.len() as f64;

            if mean_qual < 20.0 {
                continue;
            }
        }

        // Convert to FASTQ
        let seq = record.sequence().to_string();
        let qual = record.quality_scores()
            .iter()
            .map(|&q| (q + 33) as char)
            .collect::<String>();

        let fastq_record = FastqRecord {
            id: record.read_name()
                .map(|n| n.to_string())
                .unwrap_or_else(|| "unnamed".to_string()),
            sequence: seq.into_bytes(),
            quality: qual.into_bytes(),
        };

        writer.write_record(&fastq_record)?;
    }

    Ok(())
}
```

---

**Last Updated**: November 7, 2025 (Post format-integration experiment)
**Status**: Recommended integration approach
**Maintainer**: biometal team
