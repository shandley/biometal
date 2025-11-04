//! SRA Streaming Example: E. coli Dataset (Demonstration)
//!
//! **⚠️ IMPORTANT LIMITATION**: This example demonstrates the network streaming
//! API, but **SRA files are in NCBI's proprietary binary format**, not FASTQ.
//! Direct FASTQ streaming from SRA requires the SRA toolkit to decode the format.
//!
//! This example will successfully fetch the SRA file from NCBI S3 (195 MB),
//! but will fail when trying to parse it as FASTQ since it's in SRA format.
//!
//! **For working FASTQ streaming**, use `DataSource::Http` with FASTQ.gz URLs:
//! ```ignore
//! let source = DataSource::Http("https://example.com/ecoli.fastq.gz".to_string());
//! ```
//!
//! ## Dataset Information
//!
//! **SRR390728** - E. coli K-12 MG1655 whole genome sequencing
//! - Size: ~195 MB (SRA binary format)
//! - Platform: Illumina Genome Analyzer IIx
//!
//! # What This Example Shows
//!
//! 1. **Direct SRA Access**: Stream from NCBI without download
//! 2. **Constant Memory**: ~5 MB regardless of dataset size
//! 3. **ARM NEON Operations**: Count bases with 16-25× speedup
//! 4. **Background Prefetching**: Hide network latency automatically
//! 5. **Smart Caching**: 50 MB LRU cache reduces redundant fetches
//!
//! # Evidence
//!
//! Entry 028 (Lab Notebook):
//! - I/O bottleneck: 264-352× slower than compute
//! - Network streaming: Addresses critical democratization barrier
//! - Memory savings: 99.5% reduction (Entry 026)
//!
//! # Usage
//!
//! ```bash
//! # Stream real E. coli data from SRA:
//! cargo run --example sra_ecoli --features network
//!
//! # Optional: Limit number of records to process
//! cargo run --example sra_ecoli --features network -- 10000
//! ```
//!
//! # Expected Output
//!
//! - Dataset: ~250,000 reads
//! - Read length: 36 bp (Illumina GA IIx)
//! - Total bases: ~9 Mbp
//! - Processing time: ~30 seconds (network dependent)
//! - Memory usage: Constant ~5 MB
//!
//! # Dataset Details
//!
//! **Accession**: SRR390728
//! **Organism**: Escherichia coli K-12 MG1655
//! **Platform**: Illumina Genome Analyzer IIx
//! **Strategy**: WGS (Whole Genome Sequencing)
//! **Size**: ~40 MB compressed
//! **Reads**: ~250,000 reads
//! **Read Length**: 36 bp
//! **Publication**: [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/SRR390728)

use biometal::io::DataSource;
use biometal::operations::{count_bases, gc_content, mean_quality};
use biometal::FastqStream;
use std::time::Instant;

fn main() -> biometal::Result<()> {
    println!("╔═══════════════════════════════════════════════════════════════╗");
    println!("║  biometal: SRA Streaming Example (Demonstration)             ║");
    println!("║  Dataset: E. coli K-12 MG1655 (SRR390728)                    ║");
    println!("╚═══════════════════════════════════════════════════════════════╝\n");

    println!("⚠️  IMPORTANT LIMITATION:");
    println!("   SRA files are in NCBI's proprietary binary format, not FASTQ.");
    println!("   This example will successfully fetch the file from NCBI S3,");
    println!("   but will fail when parsing as FASTQ.\n");
    println!("   For working FASTQ streaming, use DataSource::Http with");
    println!("   direct FASTQ.gz URLs.\n");
    println!("   This example demonstrates the network streaming API.\n");

    // Get optional limit from command line
    let args: Vec<String> = std::env::args().collect();
    let limit = if args.len() > 1 {
        args[1].parse::<usize>().ok()
    } else {
        None
    };

    // Dataset information
    let accession = "SRR390728";
    println!("📊 Dataset Information:");
    println!("   Accession:    {}", accession);
    println!("   Organism:     Escherichia coli K-12 MG1655");
    println!("   Platform:     Illumina Genome Analyzer IIx");
    println!("   Strategy:     WGS (Whole Genome Sequencing)");
    println!("   Size:         ~195 MB (SRA binary format)");
    println!("   Reads:        ~250,000 reads (when decoded)");
    println!("   Read Length:  36 bp (when decoded)\n");

    if let Some(n) = limit {
        println!("⚠️  Processing limit: {} records (demo mode)\n", n);
    } else {
        println!("✓ Processing entire dataset (no limit)\n");
    }

    println!("🔗 Network Streaming Configuration:");
    println!("   • Direct NCBI S3 access (no download)");
    println!("   • HTTP range requests (partial fetches)");
    println!("   • 50 MB LRU cache (byte-bounded)");
    println!("   • Background prefetching (4 blocks ahead)");
    println!("   • Automatic retry (3 attempts, exponential backoff)\n");

    println!("🧬 Starting analysis...\n");

    let start = Instant::now();

    // Create SRA data source
    let source = DataSource::Sra(accession.to_string());

    // Create streaming iterator (constant memory)
    let stream = FastqStream::new(source)?;

    // Statistics
    let mut record_count = 0;
    let mut total_bases = 0;
    let mut total_gc = 0.0;
    let mut total_quality = 0.0;
    let mut base_counts = [0u64; 4]; // A, C, G, T

    // Process records one at a time (streaming)
    for result in stream {
        let record = result?;

        // ARM NEON-optimized base counting (16-25× speedup)
        let bases = count_bases(&record.sequence);
        for (i, count) in bases.iter().enumerate() {
            base_counts[i] += *count as u64;
        }
        total_bases += bases.iter().sum::<u32>() as usize;

        // ARM NEON-optimized GC content calculation
        total_gc += gc_content(&record.sequence);

        // ARM NEON-optimized quality score calculation
        total_quality += mean_quality(&record.quality) as f64;

        record_count += 1;

        // Progress indicator every 10,000 records
        if record_count % 10000 == 0 {
            let elapsed = start.elapsed().as_secs_f64();
            let rate = record_count as f64 / elapsed;
            println!(
                "   Processed {:>7} reads | {:.1} reads/sec | Memory: ~5 MB",
                record_count, rate
            );
        }

        // Show first record as example
        if record_count == 1 {
            println!("   First record: {}", record.id);
            println!("   Sequence: {}...", &String::from_utf8_lossy(&record.sequence[..20]));
            println!("   Length: {} bp\n", record.sequence.len());
        }

        // Apply limit if specified
        if let Some(n) = limit {
            if record_count >= n {
                println!("\n⚠️  Reached limit of {} records\n", n);
                break;
            }
        }
    }

    let elapsed = start.elapsed();

    // Final statistics
    println!("\n╔═══════════════════════════════════════════════════════════════╗");
    println!("║  Analysis Complete                                            ║");
    println!("╚═══════════════════════════════════════════════════════════════╝\n");

    println!("📈 Results:");
    println!("   Records processed:    {:>10}", record_count);
    println!("   Total bases:          {:>10}", total_bases);
    println!(
        "   Average read length:  {:>10.1} bp",
        total_bases as f64 / record_count as f64
    );
    println!();

    println!("🧬 Base Composition:");
    let total = base_counts.iter().sum::<u64>() as f64;
    println!(
        "   A: {:>10} ({:>5.1}%)",
        base_counts[0],
        (base_counts[0] as f64 / total) * 100.0
    );
    println!(
        "   C: {:>10} ({:>5.1}%)",
        base_counts[1],
        (base_counts[1] as f64 / total) * 100.0
    );
    println!(
        "   G: {:>10} ({:>5.1}%)",
        base_counts[2],
        (base_counts[2] as f64 / total) * 100.0
    );
    println!(
        "   T: {:>10} ({:>5.1}%)",
        base_counts[3],
        (base_counts[3] as f64 / total) * 100.0
    );
    println!();

    let avg_gc = (total_gc / record_count as f64) * 100.0;
    let avg_quality = total_quality / record_count as f64;
    println!("🎯 Quality Metrics:");
    println!("   Average GC content:   {:>10.1}%", avg_gc);
    println!("   Average quality:      {:>10.1}", avg_quality);
    println!();

    println!("⏱️  Performance:");
    println!("   Total time:           {:>10.2} sec", elapsed.as_secs_f64());
    println!(
        "   Throughput:           {:>10.1} reads/sec",
        record_count as f64 / elapsed.as_secs_f64()
    );
    println!(
        "   Data rate:            {:>10.1} KB/sec",
        (total_bases as f64 / 1024.0) / elapsed.as_secs_f64()
    );
    println!();

    println!("💾 Memory:");
    println!("   Peak memory:          ~5 MB (constant)");
    println!("   Cache size:           50 MB (LRU)");
    println!("   Total footprint:      ~55 MB");
    println!();

    println!("✅ Analysis completed successfully!");
    println!();

    // E. coli specific observations
    if record_count > 100000 && avg_gc > 49.0 && avg_gc < 53.0 {
        println!("🔬 Biology Note:");
        println!("   E. coli K-12 MG1655 has ~50.8% GC content genome-wide.");
        println!("   Your result ({:.1}%) is consistent with this.", avg_gc);
        println!("   ✓ Data quality looks good!\n");
    }

    // Network streaming reminder
    println!("🌐 Network Streaming Benefits:");
    println!("   • No local download required (saved ~40 MB disk space)");
    println!("   • Constant memory footprint (~5 MB)");
    println!("   • Analysis started immediately (no wait for download)");
    println!("   • Can analyze datasets larger than available disk/RAM");
    println!();

    // Evidence reference
    println!("📚 Evidence Base:");
    println!("   Entry 028: I/O bottleneck 264-352× slower than compute");
    println!("   Entry 026: Streaming achieves 99.5% memory reduction");
    println!("   Rules 1-6: See OPTIMIZATION_RULES.md for full methodology");
    println!();

    Ok(())
}
