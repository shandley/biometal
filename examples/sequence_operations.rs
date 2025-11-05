//! Sequence Manipulation Operations Demo (Phase 4)
//!
//! Demonstrates all sequence manipulation operations implemented in Phase 4:
//! - Core sequence transformations (reverse complement, complement, reverse)
//! - Sequence validation (DNA/RNA checking, invalid base counting)
//! - Record-level operations (extract regions, reverse complement records)
//! - Quality-based trimming (end/start/both, sliding window)
//! - Fixed-position trimming
//! - Quality-based masking (replace low-quality bases with 'N')

use biometal::operations::{
    complement, count_invalid_bases, count_masked_bases, extract_region,
    is_valid_dna, is_valid_rna, mask_low_quality, meets_length_requirement,
    reverse, reverse_complement, reverse_complement_record, sequence_length,
    trim_both, trim_quality_both, trim_quality_end, trim_quality_window,
};
use biometal::{FastqRecord, FastqStream};
use std::io::Cursor;

fn main() -> biometal::Result<()> {
    println!("biometal Sequence Operations Demo (Phase 4)");
    println!("============================================\n");

    demo_sequence_transformations();
    println!();
    demo_sequence_validation();
    println!();
    demo_record_operations()?;
    println!();
    demo_trimming_operations()?;
    println!();
    demo_masking_operations()?;

    Ok(())
}

/// Demo 1: Core sequence transformations
fn demo_sequence_transformations() {
    println!("1. Sequence Transformations");
    println!("---------------------------");

    let seq = b"GATTACA";

    println!("Original sequence:      {}", String::from_utf8_lossy(seq));
    println!("Reverse:                {}", String::from_utf8_lossy(&reverse(seq)));
    println!("Complement:             {}", String::from_utf8_lossy(&complement(seq)));
    println!("Reverse complement:     {}", String::from_utf8_lossy(&reverse_complement(seq)));

    println!("\nAmbiguous bases (IUPAC codes):");
    let ambig = b"ACGTRYSWKM";
    println!("Original:               {}", String::from_utf8_lossy(ambig));
    println!("Complement:             {}", String::from_utf8_lossy(&complement(ambig)));

    println!("\nRNA handling (produces DNA-style output):");
    let rna = b"AUGC";
    println!("RNA input:              {}", String::from_utf8_lossy(rna));
    println!("Complement:             {} (note: T not U)", String::from_utf8_lossy(&complement(rna)));
}

/// Demo 2: Sequence validation
fn demo_sequence_validation() {
    println!("2. Sequence Validation");
    println!("----------------------");

    let sequences = vec![
        (b"ACGT" as &[u8], "Pure DNA"),
        (b"ACGTRYSWKM" as &[u8], "DNA with ambiguous"),
        (b"AUGC" as &[u8], "RNA sequence"),
        (b"ACGTX" as &[u8], "Invalid base (X)"),
        (b"ACG123" as &[u8], "Non-nucleotide"),
    ];

    for (seq, label) in sequences {
        let is_dna = is_valid_dna(seq);
        let is_rna = is_valid_rna(seq);
        let invalid_count = count_invalid_bases(seq);

        println!("{:25} seq={:12} dna={:5} rna={:5} invalid={}",
                 label,
                 String::from_utf8_lossy(seq),
                 is_dna,
                 is_rna,
                 invalid_count);
    }
}

/// Demo 3: Record-level operations
fn demo_record_operations() -> biometal::Result<()> {
    println!("3. Record Operations");
    println!("--------------------");

    let record = FastqRecord::new(
        "forward_read".to_string(),
        b"GATTACAGATTACA".to_vec(),
        b"IIIIIIIIIIIIII".to_vec(),
    );

    println!("Original record:");
    println!("  ID:       {}", record.id);
    println!("  Sequence: {}", String::from_utf8_lossy(&record.sequence));
    println!("  Length:   {}", sequence_length(&record));

    // Extract region
    println!("\nExtract region [5..10]:");
    let region = extract_region(&record, 5, 10)?;
    println!("  Sequence: {}", String::from_utf8_lossy(&region.sequence));
    println!("  Quality:  {}", String::from_utf8_lossy(&region.quality));

    // Reverse complement
    println!("\nReverse complement:");
    let rc = reverse_complement_record(&record);
    println!("  Sequence: {}", String::from_utf8_lossy(&rc.sequence));

    // Length filtering
    println!("\nLength filtering:");
    let length_tests = vec![(10, 20, "10-20 bp"), (14, 14, "exactly 14 bp"), (15, 100, "15-100 bp")];
    for (min_len, max_len, label) in length_tests {
        let meets = meets_length_requirement(&record, min_len, max_len);
        println!("  {}: {}", label, if meets { "PASS ✓" } else { "FAIL ✗" });
    }

    Ok(())
}

/// Demo 4: Trimming operations (fixed + quality-based)
fn demo_trimming_operations() -> biometal::Result<()> {
    println!("4. Trimming Operations");
    println!("----------------------");

    // Demo 4a: Fixed-position trimming
    println!("4a. Fixed-position trimming:");
    let record = FastqRecord::new(
        "adapter_contamination".to_string(),
        b"NNNGATTACAGATTACANNN".to_vec(),
        b"!!!IIIIIIIIIIIIII!!!".to_vec(),
    );

    println!("  Original: {}", String::from_utf8_lossy(&record.sequence));
    println!("  Quality:  {}", String::from_utf8_lossy(&record.quality));

    let trimmed = trim_both(&record, 3, 3)?;
    println!("  Trim 3bp each end: {}", String::from_utf8_lossy(&trimmed.sequence));

    // Demo 4b: Quality-based trimming
    println!("\n4b. Quality-based trimming (threshold Q20):");

    let qual_record = FastqRecord::new(
        "degraded_end".to_string(),
        b"GATTACAGATTACA".to_vec(),
        b"IIIIIIIIII!!!!".to_vec(),  // High quality start, low quality end
    );

    println!("  Original: {}", String::from_utf8_lossy(&qual_record.sequence));
    println!("  Quality:  {}", String::from_utf8_lossy(&qual_record.quality));
    println!("            {}{}",
             "^^^^^^^^^^".to_string(),  // High quality (I = Q40)
             "!!!!".to_string());       // Low quality (! = Q0)

    let trimmed = trim_quality_end(&qual_record, 20)?;
    println!("  After trim_quality_end(Q20): {}", String::from_utf8_lossy(&trimmed.sequence));

    // Demo 4c: Sliding window trimming
    println!("\n4c. Sliding window trimming (Q20, window=4):");

    let window_record = FastqRecord::new(
        "variable_quality".to_string(),
        b"GATTACAGATTACAGATTACA".to_vec(),
        b"IIIII55555IIIII555555".to_vec(),  // Variable quality
    );

    println!("  Original: {}", String::from_utf8_lossy(&window_record.sequence));
    println!("  Quality:  {}", String::from_utf8_lossy(&window_record.quality));

    let trimmed = trim_quality_window(&window_record, 20, 4)?;
    println!("  After sliding window: {}", String::from_utf8_lossy(&trimmed.sequence));
    println!("  (Stops when 4bp window drops below Q20)");

    Ok(())
}

/// Demo 5: Quality-based masking
fn demo_masking_operations() -> biometal::Result<()> {
    println!("5. Masking Operations");
    println!("---------------------");
    println!("(Replace low-quality bases with 'N' instead of removing)");

    let mut record = FastqRecord::new(
        "variant_calling".to_string(),
        b"GATTACAGATTACA".to_vec(),
        b"IIII!!IIII!!II".to_vec(),  // Mixed quality
    );

    println!("\nOriginal:");
    println!("  Sequence: {}", String::from_utf8_lossy(&record.sequence));
    println!("  Quality:  {}", String::from_utf8_lossy(&record.quality));
    println!("            ^^^^  ^^^^  ^^  (high quality)");
    println!("                ^^    ^^    (low quality, will be masked)");

    // Apply masking
    mask_low_quality(&mut record, 20)?;
    println!("\nAfter masking (Q20):");
    println!("  Sequence: {}", String::from_utf8_lossy(&record.sequence));
    println!("  Quality:  {}", String::from_utf8_lossy(&record.quality));

    // Count masked bases
    let mask_count = count_masked_bases(&record);
    println!("\n  Masked bases: {}", mask_count);
    println!("  Note: Length preserved (important for alignment)");

    // Use case explanation
    println!("\nWhy masking instead of trimming?");
    println!("  - Preserves read length (required by some aligners)");
    println!("  - Maintains positional information");
    println!("  - Useful for variant calling (exclude low-quality calls)");
    println!("  - Better for paired-end read alignment");

    Ok(())
}

/// Demo 6: Streaming processing pipeline
#[allow(dead_code)]
fn demo_streaming_pipeline() -> biometal::Result<()> {
    println!("6. Streaming Pipeline");
    println!("---------------------");
    println!("Processing FASTQ with quality trimming and length filtering:\n");

    let fastq_data = b"@read1 good quality
GATTACAGATTACAGATTACA
+
IIIIIIIIIIIIIIIIIIIII
@read2 low quality ends
GATTACAGATTACAGATTACA
+
!!!IIIIIIIIIIIII!!!!!
@read3 very short
GATTACA
+
IIIIIII
";

    let stream = FastqStream::from_reader(Cursor::new(fastq_data));
    let min_length = 15;
    let min_quality = 20;

    let mut passed = 0;
    let mut too_short = 0;
    let mut trimmed = 0;

    for record in stream {
        let record = record?;
        println!("Processing: {}", record.id);

        // Step 1: Quality trim both ends
        let record = trim_quality_both(&record, min_quality)?;
        let trimmed_len = sequence_length(&record);

        println!("  After quality trim: {} bp", trimmed_len);

        // Step 2: Check length requirement
        if meets_length_requirement(&record, min_length, 1000) {
            println!("  Status: PASS ✓ ({} bp in range)", trimmed_len);
            passed += 1;
        } else {
            println!("  Status: FILTERED ✗ ({} bp too short)", trimmed_len);
            too_short += 1;
        }

        if trimmed_len < sequence_length(&record) {
            trimmed += 1;
        }
    }

    println!("\nSummary:");
    println!("  Passed:      {}", passed);
    println!("  Too short:   {}", too_short);
    println!("  Trimmed:     {}", trimmed);

    Ok(())
}
