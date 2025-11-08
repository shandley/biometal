use std::env;
use std::fs::File;
use std::io::{self, BufReader};
use std::time::Instant;

use noodles_bam as bam;

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <input.bam>", args[0]);
        eprintln!();
        eprintln!("This profiling tool parses a BAM file using noodles and reports");
        eprintln!("throughput metrics. Use with macOS Instruments Time Profiler to");
        eprintln!("identify CPU hotspots (especially sequence decoding).");
        std::process::exit(1);
    }

    let path = &args[1];

    println!("Profiling BAM parsing with noodles");
    println!("File: {}", path);
    println!();

    profile_parse(path)?;

    Ok(())
}

fn profile_parse(path: &str) -> io::Result<()> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut bam_reader = bam::io::Reader::new(reader);

    // Read header
    let header_start = Instant::now();
    let _header = bam_reader.read_header()?;
    let header_time = header_start.elapsed();

    println!("Header parsing: {:?}", header_time);
    println!();

    // Parse all records
    let mut record = bam::Record::default();
    let mut count = 0u64;
    let mut total_seq_len = 0u64;
    let mut total_qual_len = 0u64;

    let parse_start = Instant::now();

    loop {
        match bam_reader.read_record(&mut record) {
            Ok(0) => break, // EOF
            Ok(_) => {
                count += 1;

                // Access sequence to ensure it's decoded
                let seq = record.sequence();
                total_seq_len += seq.len() as u64;

                // Access quality scores
                let qual = record.quality_scores();
                total_qual_len += qual.len() as u64;

                // Print progress every 10K records
                if count % 10_000 == 0 {
                    let elapsed = parse_start.elapsed();
                    let records_per_sec = count as f64 / elapsed.as_secs_f64();
                    print!("\rParsed {} records ({:.1} rec/s)...", count, records_per_sec);
                    io::Write::flush(&mut io::stdout())?;
                }
            }
            Err(e) => {
                eprintln!("\nError reading record {}: {}", count + 1, e);
                return Err(e);
            }
        }
    }

    let total_time = parse_start.elapsed();

    println!("\r                                                              ");
    println!("Parsing complete!");
    println!();
    println!("Results:");
    println!("  Records parsed: {}", count);
    println!("  Total sequence length: {} bp", total_seq_len);
    println!("  Total quality scores: {}", total_qual_len);
    println!("  Elapsed time: {:.3} s", total_time.as_secs_f64());
    println!();
    println!("Throughput:");
    println!("  {:.1} records/sec", count as f64 / total_time.as_secs_f64());
    println!("  {:.1} Mrec/sec", count as f64 / total_time.as_secs_f64() / 1e6);
    println!("  {:.1} Mbp/sec", total_seq_len as f64 / total_time.as_secs_f64() / 1e6);
    println!();
    println!("Profiling Notes:");
    println!("  - Run with: cargo build --release && instruments -t 'Time Profiler' ./target/release/bam-profiling <file.bam>");
    println!("  - Or use: cargo flamegraph --bin bam-profiling -- <file.bam>");
    println!("  - Look for CPU time in:");
    println!("    * decode_base (sequence decoding)");
    println!("    * read_sequence (4-bit unpacking)");
    println!("    * BGZF decompression");
    println!("  - GO Criteria: Sequence decoding ≥15% CPU time");
    println!("  - NO-GO Risk: Sequence decoding <10% CPU time");

    Ok(())
}
