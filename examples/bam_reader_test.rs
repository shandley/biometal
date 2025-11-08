//! Simple test of BAM reader with synthetic test data.
//!
//! This example demonstrates basic BAM reading functionality.

use biometal::io::bam::BamReader;
use std::io::BufReader;
use std::fs::File;

fn main() -> std::io::Result<()> {
    // Path to test BAM file
    let path = "experiments/native-bam-implementation/test-data/synthetic_100000.bam";

    println!("Opening BAM file: {}", path);

    // Open BAM file
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut bam = BamReader::new(reader)?;

    // Print header info
    let header = bam.header();
    println!("\nHeader:");
    println!("  References: {}", header.reference_count());
    for (i, reference) in header.references.iter().enumerate() {
        println!("    {}: {} ({} bp)", i, reference.name, reference.length);
    }

    // Read and count records
    let mut count = 0;
    let mut total_bases = 0u64;

    for result in bam.records() {
        let record = result?;
        count += 1;
        total_bases += record.sequence.len() as u64;

        // Print first 5 records for verification
        if count <= 5 {
            println!("\nRecord {}:", count);
            println!("  Name: {}", record.name);
            println!("  Reference: {:?}", record.reference_id);
            println!("  Position: {:?}", record.position);
            println!("  MAPQ: {:?}", record.mapq);
            println!("  Flags: 0x{:04x}", record.flags);
            println!("  Sequence length: {}", record.sequence.len());
            println!("  CIGAR ops: {}", record.cigar.len());
            if record.sequence.len() <= 20 {
                println!("  Sequence: {}", String::from_utf8_lossy(&record.sequence));
            }
        }
    }

    println!("\n\nSummary:");
    println!("  Total records: {}", count);
    println!("  Total bases: {}", total_bases);
    println!("  Average read length: {:.1}", total_bases as f64 / count as f64);

    Ok(())
}
