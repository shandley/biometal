//! Example: Using CSI (Coordinate-Sorted Index) for fast region queries
//!
//! This example demonstrates how to:
//! 1. Load a CSI index file
//! 2. Query regions by reference name or index
//! 3. Use the chunks for targeted file access
//!
//! CSI is the successor to BAI/TBI with support for larger chromosomes
//! and configurable binning parameters.

use biometal::formats::index::CsiIndex;
use biometal::Result;

fn main() -> Result<()> {
    // Example 1: Load CSI index and inspect metadata
    println!("=== Example 1: Load and Inspect CSI Index ===\n");

    // Load the index (in real usage, use actual .csi file)
    // let index = CsiIndex::from_path("alignments.bam.csi")?;

    // For this example, create a minimal in-memory index
    let index = create_example_index()?;

    println!("CSI Index Parameters:");
    println!("  Min shift: {} (bin size: {} kb)",
        index.min_shift(),
        (1 << index.min_shift()) / 1024);
    println!("  Depth: {} levels", index.depth());
    println!("  References: {}", index.references().len());
    println!();

    // Example 2: Query region by index
    println!("=== Example 2: Query by Reference Index ===\n");

    let start = 1_000_000;
    let end = 2_000_000;

    match index.query_by_index(0, start, end)? {
        Some(chunks) => {
            println!("Query: ref_idx=0, region {}-{}", start, end);
            println!("Found {} chunks:", chunks.len());
            for (i, chunk) in chunks.iter().enumerate() {
                println!(
                    "  Chunk {}: offset {:#x} - {:#x} ({} bytes)",
                    i + 1,
                    chunk.start.as_raw(),
                    chunk.end.as_raw(),
                    chunk.end.as_raw() - chunk.start.as_raw()
                );
            }
        }
        None => println!("Reference index 0 not found"),
    }
    println!();

    // Example 3: Query region by name (if names available in aux data)
    println!("=== Example 3: Query by Reference Name ===\n");

    if let Some(chunks) = index.query("chr1", start, end)? {
        println!("Query: chr1:{}-{}", start, end);
        println!("Found {} chunks", chunks.len());
    } else {
        println!("Reference 'chr1' not found in index");
        println!("(This is expected for this example - aux data would need reference names)");
    }
    println!();

    // Example 4: Inspect reference information
    println!("=== Example 4: Reference Information ===\n");

    for (idx, reference) in index.references().iter().enumerate() {
        println!("Reference {}:", idx);
        if let Some(name) = &reference.name {
            println!("  Name: {}", name);
        } else {
            println!("  Name: <not available>");
        }
        println!("  Bins: {} (Note: CSI uses bins only, no linear intervals)", reference.bins.len());

        if !reference.bins.is_empty() {
            println!("  First bin:");
            let bin = &reference.bins[0];
            println!("    ID: {}", bin.bin_id);
            println!("    Left offset: {:#x}", bin.loffset.as_raw());
            println!("    Chunks: {}", bin.chunks.len());
        }
    }
    println!();

    // Example 5: Compare CSI parameters for different use cases
    println!("=== Example 5: CSI Parameter Examples ===\n");

    println!("Common CSI configurations:");
    println!("  Default (samtools):  min_shift=14 (16kb), depth=5");
    println!("    → Max chromosome: 2^(14+3*5) = 512 Mbp (same as BAI)");
    println!();
    println!("  Large genomes:       min_shift=18 (256kb), depth=5");
    println!("    → Max chromosome: 2^(18+3*5) = 8 Gbp");
    println!();
    println!("  Very large:          min_shift=20 (1Mbp), depth=6");
    println!("    → Max chromosome: 2^(20+3*6) = 64 Gbp");
    println!();

    println!("Current index parameters:");
    println!("  Min shift: {} → bin size: {} kb",
        index.min_shift(),
        (1 << index.min_shift()) / 1024);
    println!("  Depth: {} → max chromosome: {} Gbp",
        index.depth(),
        (1u64 << (index.min_shift() + 3 * index.depth())) / 1_000_000_000);

    Ok(())
}

/// Helper function to create an example CSI index for demonstration
fn create_example_index() -> Result<CsiIndex> {
    use std::io::Cursor;

    // Create a minimal valid CSI index
    let mut data = Vec::new();

    // Magic
    data.extend_from_slice(b"CSI\x01");

    // min_shift (14 = 16kb)
    data.extend_from_slice(&14i32.to_le_bytes());

    // depth (5 levels)
    data.extend_from_slice(&5i32.to_le_bytes());

    // aux_size (0 = no aux data)
    data.extend_from_slice(&0i32.to_le_bytes());

    // n_ref (1 reference)
    data.extend_from_slice(&1i32.to_le_bytes());

    // Reference 0:
    // n_bin (1 bin)
    data.extend_from_slice(&1i32.to_le_bytes());

    // Bin 0 (entire sequence):
    data.extend_from_slice(&0u32.to_le_bytes()); // bin_id
    data.extend_from_slice(&0u64.to_le_bytes()); // loffset
    data.extend_from_slice(&2i32.to_le_bytes()); // n_chunk

    // Chunk 1:
    data.extend_from_slice(&0x1000u64.to_le_bytes()); // chunk_beg
    data.extend_from_slice(&0x2000u64.to_le_bytes()); // chunk_end

    // Chunk 2:
    data.extend_from_slice(&0x2000u64.to_le_bytes()); // chunk_beg
    data.extend_from_slice(&0x3000u64.to_le_bytes()); // chunk_end

    // n_intv (0 intervals)
    data.extend_from_slice(&0i32.to_le_bytes());

    let mut cursor = Cursor::new(data);
    CsiIndex::parse(&mut cursor)
}
