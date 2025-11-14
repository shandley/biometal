//! Demonstration of format primitives.
//!
//! This example shows how to use the generic tab-delimited parser
//! to create a simple custom format parser.

use biometal::formats::primitives::{
    FormatError, GenomicInterval, Result, Strand, TabDelimitedParser, TabDelimitedRecord,
};
use std::str::FromStr;

/// A simple BED-like record (chromosome, start, end, name, strand).
#[derive(Debug, PartialEq)]
struct SimpleBedRecord {
    interval: GenomicInterval,
    name: String,
    strand: Strand,
}

impl TabDelimitedRecord for SimpleBedRecord {
    fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 6 {
            return Err(FormatError::FieldCount {
                expected: 6,
                actual: fields.len(),
                line: 0,
            });
        }

        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().map_err(|e| FormatError::InvalidField {
            field: "start".to_string(),
            line: 0,
            reason: format!("{}", e),
        })?;
        let end: u64 = fields[2].parse().map_err(|e| FormatError::InvalidField {
            field: "end".to_string(),
            line: 0,
            reason: format!("{}", e),
        })?;

        let interval = GenomicInterval::new(chrom, start, end)?;
        let name = fields[3].to_string();
        let strand = Strand::from_str(fields[5])?;

        Ok(SimpleBedRecord {
            interval,
            name,
            strand,
        })
    }

    fn to_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t0\t{}",
            self.interval.chrom,
            self.interval.start,
            self.interval.end,
            self.name,
            self.strand
        )
    }

    fn expected_fields() -> Option<usize> {
        Some(6)
    }
}

fn main() -> Result<()> {
    println!("=== Format Primitives Demo ===\n");

    // Example 1: Genomic Intervals
    println!("1. Genomic Intervals:");
    let int1 = GenomicInterval::new("chr1".to_string(), 100, 200)?;
    let int2 = GenomicInterval::new("chr1".to_string(), 150, 250)?;

    println!("  Interval 1: {}", int1);
    println!("  Interval 2: {}", int2);
    println!("  Length: {}", int1.length());
    println!("  Overlaps: {}", int1.overlaps(&int2));
    println!("  Contains: {}\n", int1.contains(&int2));

    // Example 2: Strand parsing
    println!("2. Strand Parsing:");
    let forward = Strand::from_str("+")?;
    let reverse = Strand::from_str("-")?;
    let unknown = Strand::from_str(".")?;

    println!("  Forward: {}", forward);
    println!("  Reverse: {}", reverse);
    println!("  Unknown: {}\n", unknown);

    // Example 3: Tab-delimited parsing
    println!("3. Tab-Delimited Parsing:");

    let bed_data = "\
# This is a comment
chr1\t1000\t2000\tgene1\t0\t+
chr1\t3000\t4000\tgene2\t0\t-
# Another comment
chr2\t5000\t6000\tgene3\t0\t.
";

    let parser = TabDelimitedParser::<_, SimpleBedRecord>::new(bed_data.as_bytes());

    println!("  Parsing BED-like data:");
    for (i, record) in parser.enumerate() {
        let record = record?;
        println!(
            "  Record {}: {} ({}:{}-{}, strand: {})",
            i + 1,
            record.name,
            record.interval.chrom,
            record.interval.start,
            record.interval.end,
            record.strand
        );
    }

    // Example 4: Round-trip serialization
    println!("\n4. Round-Trip Serialization:");
    let original = SimpleBedRecord {
        interval: GenomicInterval::new("chr1".to_string(), 1000, 2000)?,
        name: "test_gene".to_string(),
        strand: Strand::Forward,
    };

    let serialized = original.to_line();
    println!("  Serialized: {}", serialized);

    let parsed = SimpleBedRecord::from_line(&serialized)?;
    println!("  Parsed back: {} ({})", parsed.name, parsed.interval);
    println!("  Round-trip success: {}", original == parsed);

    println!("\n=== Demo Complete ===");

    Ok(())
}
