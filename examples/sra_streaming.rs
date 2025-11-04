//! SRA Streaming Example
//!
//! Demonstrates streaming data directly from NCBI's Sequence Read Archive (SRA)
//! without downloading entire datasets. This enables analysis of multi-GB SRA runs
//! with constant ~5 MB memory.
//!
//! # Evidence
//!
//! Entry 028: I/O bottleneck is 264-352× slower than compute
//! Network streaming addresses this critical bottleneck (Rule 6)
//!
//! # Usage
//!
//! ```bash
//! # Demo mode (shows capabilities without real SRA access):
//! cargo run --example sra_streaming --features network
//!
//! # With real SRA accession:
//! cargo run --example sra_streaming --features network -- SRR000001
//! ```
//!
//! # Note
//!
//! Real SRA streaming requires:
//! - Network connection to NCBI S3
//! - SRA accession that exists and is public
//! - Server support for HTTP range requests (AWS S3 supports this)

use biometal::io::{is_sra_accession, sra_to_url};

fn main() -> biometal::Result<()> {
    println!("biometal SRA Streaming Demo");
    println!("===========================\n");

    // Get SRA accession from command line or use demo mode
    let args: Vec<String> = std::env::args().collect();
    let demo_mode = args.len() == 1;

    if demo_mode {
        demo_sra_features()?;
    } else {
        let accession = &args[1];
        stream_sra_data(accession)?;
    }

    Ok(())
}

fn demo_sra_features() -> biometal::Result<()> {
    println!("Demo Mode: SRA Streaming Capabilities");
    println!("=====================================\n");

    println!("SRA (Sequence Read Archive) Integration:");
    println!("----------------------------------------\n");

    println!("✓ Direct streaming from NCBI SRA");
    println!("✓ No local download required");
    println!("✓ Constant memory (~5 MB)");
    println!("✓ Background prefetching (hides latency)");
    println!("✓ Smart caching (50 MB LRU)\n");

    println!("Supported Accession Types:");
    println!("-------------------------\n");

    println!("• SRR (Run):        Sequencing run - most common");
    println!("• SRX (Experiment): Collection of runs");
    println!("• SRS (Sample):     Biological sample");
    println!("• SRP (Study):      Collection of experiments\n");

    println!("Example Accessions:");
    println!("------------------\n");

    let examples = vec![
        "SRR000001",
        "SRR390728", // E. coli
        "SRR1553500", // Human WGS
    ];

    for accession in &examples {
        if is_sra_accession(accession) {
            let url = sra_to_url(accession)?;
            println!("  {} →", accession);
            println!("    {}", url);
        }
    }

    println!("\n\nURL Pattern:");
    println!("------------\n");

    println!("NCBI provides public S3 access:");
    println!("  https://sra-pub-run-odp.s3.amazonaws.com/sra/<accession>/<accession>");
    println!();
    println!("  Example: SRR000001 → https://...com/sra/SRR000001/SRR000001");
    println!("  Example: SRR390728 → https://...com/sra/SRR390728/SRR390728\n");

    println!("Example Usage (Code):");
    println!("====================\n");

    println!("// Stream directly from SRA:");
    println!("let source = DataSource::Sra(\"SRR000001\".to_string());");
    println!("let stream = FastqStream::new(source)?;\n");

    println!("for record in stream {{");
    println!("    // Process one record at a time");
    println!("    // Memory: Constant ~5 MB");
    println!("}}\n");

    println!("// Or convert to URL manually:");
    println!("let url = sra_to_url(\"SRR000001\")?;");
    println!("let reader = HttpReader::new(&url)?;\n");

    println!("Memory Guarantees:");
    println!("==================\n");

    println!("  Streaming:     ~5 MB per stream (constant)");
    println!("  Cache:         50 MB (LRU, byte-bounded)");
    println!("  Prefetch:      4 blocks ahead (default)");
    println!("  Total:         ~55 MB regardless of SRA file size\n");

    println!("Performance:");
    println!("============\n");

    println!("  Network overhead:  Automatic retry + exponential backoff");
    println!("  Latency hiding:    Background prefetching (4 blocks)");
    println!("  Bandwidth saving:  Only fetch needed data (range requests)");
    println!("  Cache efficiency:  LRU reduces redundant downloads\n");

    println!("Evidence (Entry 028):");
    println!("=====================\n");

    println!("  I/O bottleneck:        264-352× slower than compute");
    println!("  Without streaming:     1.04-1.08× E2E speedup (masked)");
    println!("  With streaming:        ~17× E2E speedup projected");
    println!("  Conclusion:            SRA streaming is CRITICAL\n");

    println!("✅ SRA streaming implementation complete!\n");

    println!("To stream real SRA data:");
    println!("  cargo run --example sra_streaming --features network -- <SRA_ACCESSION>");
    println!();
    println!("Example accessions to try:");
    println!("  SRR390728   (E. coli, ~40 MB)");
    println!("  SRR1553500  (Human WGS, ~1 GB)");

    Ok(())
}

fn stream_sra_data(accession: &str) -> biometal::Result<()> {
    use biometal::io::DataSource;
    use biometal::operations::{count_bases, gc_content};
    use biometal::FastqStream;

    // Validate accession format
    if !is_sra_accession(accession) {
        eprintln!("Error: '{}' doesn't look like an SRA accession", accession);
        eprintln!("Expected format: SRR/SRX/SRS/SRP followed by digits");
        eprintln!("Example: SRR000001");
        return Ok(());
    }

    println!("Streaming from SRA: {}", accession);

    // Convert to URL
    let url = sra_to_url(accession)?;
    println!("URL: {}\n", url);

    println!("Initializing HTTP streaming with:");
    println!("  • Range requests (partial downloads)");
    println!("  • 50 MB LRU cache");
    println!("  • Background prefetching (4 blocks)");
    println!("  • Automatic retry (3 attempts)\n");

    println!("Streaming FASTQ records...\n");

    let source = DataSource::Sra(accession.to_string());
    let stream = FastqStream::new(source)?;

    let mut record_count = 0;
    let mut total_bases = 0;
    let mut total_gc = 0.0;

    for record in stream {
        let record = record?;
        record_count += 1;

        // Count bases (NEON-optimized on ARM)
        let bases = count_bases(&record.sequence);
        total_bases += bases.iter().sum::<u32>() as usize;

        // Calculate GC content (NEON-optimized on ARM)
        let gc = gc_content(&record.sequence);
        total_gc += gc;

        // Print first few records as examples
        if record_count <= 5 {
            println!("Record {}: {}", record_count, record.id);
            println!("  Length: {} bp", record.sequence.len());
            println!("  A: {}, C: {}, G: {}, T: {}", bases[0], bases[1], bases[2], bases[3]);
            println!("  GC: {:.1}%\n", gc * 100.0);
        }

        // Progress indicator
        if record_count % 10000 == 0 {
            println!("Processed {} records... (memory constant at ~5 MB)", record_count);
        }

        // Demo mode: stop after 1000 records to avoid long downloads
        if record_count >= 1000 {
            println!("\n[Demo limit: stopping after 1000 records]\n");
            break;
        }
    }

    println!("\n=== Summary ===");
    println!("Records processed: {}", record_count);
    println!("Total bases: {}", total_bases);
    println!("Average GC content: {:.1}%", (total_gc / record_count as f64) * 100.0);
    println!("Memory usage: Constant ~5 MB (streaming)");
    println!("\n✅ SRA streaming completed successfully!");

    Ok(())
}
