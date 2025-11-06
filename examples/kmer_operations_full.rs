//! Comprehensive k-mer operations example
//!
//! Demonstrates all k-mer functionality in biometal with evidence-based design.
//! This example shows:
//! - Simple k-mer extraction (scalar)
//! - Parallel k-mer extraction (2.2× for large datasets)
//! - Minimizer extraction (minimap2-style)
//! - K-mer spectrum analysis (frequency counting)
//! - Streaming integration (constant memory)
//!
//! Evidence: ASBB Entry 034 - K-mer operations are data-structure-bound

use biometal::operations::kmer::{
    extract_kmers, extract_minimizers, kmer_spectrum, KmerExtractor,
};
use biometal::FastqStream;

fn main() -> biometal::Result<()> {
    println!("=== biometal K-mer Operations ===\n");

    // ===== 1. Simple K-mer Extraction =====
    println!("1. Simple K-mer Extraction (Scalar)");
    println!("-----------------------------------");

    let sequence = b"ATGCATGCATGCATGC";
    let kmers = extract_kmers(sequence, 6);

    println!("Sequence: {}", String::from_utf8_lossy(sequence));
    println!("K-mers (k=6): {} total", kmers.len());
    for (i, kmer) in kmers.iter().take(5).enumerate() {
        println!("  [{}] {}", i, String::from_utf8_lossy(kmer));
    }
    println!();

    // ===== 2. Minimizer Extraction (minimap2-style) =====
    println!("2. Minimizer Extraction (minimap2-style)");
    println!("----------------------------------------");

    let sequence = b"ATGCATGCATGCATGCATGCATGCATGC";
    let minimizers = extract_minimizers(sequence, 6, 10);

    println!("Sequence: {} ({}bp)", String::from_utf8_lossy(sequence), sequence.len());
    println!("Minimizers (k=6, w=10): {} total", minimizers.len());
    for minimizer in &minimizers {
        println!("  Position {}: {} (hash={})",
                 minimizer.position,
                 String::from_utf8_lossy(&minimizer.kmer),
                 minimizer.hash);
    }
    println!("\nNote: Entry 034 validated minimap2's scalar design (1.26× parallel, below threshold)");
    println!();

    // ===== 3. K-mer Spectrum (Frequency Counting) =====
    println!("3. K-mer Spectrum (Frequency Counting)");
    println!("--------------------------------------");

    let sequences = vec![
        b"ATGCATGC".as_ref(),
        b"GCATGCAT".as_ref(),
        b"CATGCATG".as_ref(),
    ];

    let spectrum = kmer_spectrum(&sequences, 4);

    println!("Sequences: {} total", sequences.len());
    println!("K-mer spectrum (k=4): {} unique k-mers", spectrum.len());

    // Show top 5 most frequent k-mers
    let mut spectrum_vec: Vec<_> = spectrum.iter().collect();
    spectrum_vec.sort_by(|a, b| b.1.cmp(a.1));

    println!("Top 5 most frequent:");
    for (kmer, count) in spectrum_vec.iter().take(5) {
        println!("  {} : {}", String::from_utf8_lossy(kmer), count);
    }
    println!("\nNote: Entry 034 found parallel makes this SLOWER (HashMap contention)");
    println!();

    // ===== 4. Parallel K-mer Extraction (2.2× for large datasets) =====
    println!("4. Parallel K-mer Extraction (Opt-in)");
    println!("-------------------------------------");

    // Simulate large dataset (1000+ sequences)
    let large_dataset: Vec<&[u8]> = (0..1000)
        .map(|_| b"ATGCATGCATGCATGCATGCATGC".as_ref())
        .collect();

    // Scalar (default)
    let scalar_extractor = KmerExtractor::new();
    let scalar_kmers = scalar_extractor.extract(&large_dataset, 6);
    println!("Scalar extraction: {} k-mers from {} sequences",
             scalar_kmers.len(), large_dataset.len());

    // Parallel (opt-in, 2.2× speedup for large datasets)
    let parallel_extractor = KmerExtractor::with_parallel(4);
    let parallel_kmers = parallel_extractor.extract(&large_dataset, 6);
    println!("Parallel extraction (4 threads): {} k-mers from {} sequences",
             parallel_kmers.len(), large_dataset.len());

    println!("\nNote: Entry 034 measured 2.19-2.38× speedup with Parallel-4t");
    println!("Threads capped at 4 (optimal, no benefit beyond)");
    println!();

    // ===== 5. Streaming Integration (Constant Memory) =====
    println!("5. Streaming Integration (Constant Memory)");
    println!("------------------------------------------");

    println!("K-mer extraction integrates with streaming for constant memory:");
    println!();
    println!("```rust");
    println!("let stream = FastqStream::from_path(\"large_dataset.fq.gz\")?;");
    println!();
    println!("for record in stream {{");
    println!("    // Extract k-mers from each record (constant memory)");
    println!("    let kmers = extract_kmers(&record.sequence, 6);");
    println!("    ");
    println!("    // Feed to DNABert/ML model immediately");
    println!("    // No accumulation - memory stays at ~5 MB");
    println!("}}");
    println!("```");
    println!();
    println!("This enables processing 5TB datasets on consumer hardware!");
    println!();

    // ===== 6. Evidence Summary =====
    println!("=== Evidence Summary (ASBB Entry 034) ===");
    println!("-----------------------------------------");
    println!();
    println!("K-mer operations are DATA-STRUCTURE-BOUND (not compute-bound):");
    println!();
    println!("Runtime breakdown:");
    println!("  - Hash computation: 50-60% (sequential, can't vectorize)");
    println!("  - HashMap operations: 30-40% (sequential, thread contention)");
    println!("  - Base validation: 5-10% (only NEON-friendly part)");
    println!();
    println!("Results:");
    println!("  - Minimizers:   1.02-1.26× (NEON/Parallel) → Scalar-only");
    println!("  - Spectrum:     0.95-1.88× (sometimes SLOWER!) → Scalar-only");
    println!("  - Extraction:   2.19-2.38× (Parallel-4t) → Opt-in parallel");
    println!();
    println!("Key insight: NEON/GPU provide NO benefit for k-mer operations");
    println!("(unlike base_counting 16.7×, gc_content 20.3×, quality_filter 25.1×)");
    println!();
    println!("Validated existing tools:");
    println!("  - minimap2: Scalar minimizers confirmed optimal");
    println!("  - DNABert: 2.2× k-mer extraction speedup opportunity");

    Ok(())
}
