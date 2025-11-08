use std::fs::File;
use std::io::{self, BufWriter};

use noodles_bam as bam;
use noodles_sam as sam;

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: {} <output.bam> <num_records>", args[0]);
        std::process::exit(1);
    }

    let output_path = &args[1];
    let num_records: usize = args[2].parse().expect("Invalid number");

    println!("Generating {} records to {}...", num_records, output_path);

    // Create header
    let header = sam::Header::builder()
        .set_header(sam::header::header::Header::default())
        .add_reference_sequence(
            "chr1".parse().unwrap(),
            sam::header::reference_sequence::ReferenceSequence::new("chr1".into(), 248956422)
                .map(|b| b.build())
                .unwrap(),
        )
        .add_reference_sequence(
            "chr22".parse().unwrap(),
            sam::header::reference_sequence::ReferenceSequence::new("chr22".into(), 50818468)
                .map(|b| b.build())
                .unwrap(),
        )
        .add_program("generate_bam".parse().unwrap(), sam::header::Program::default())
        .build();

    // Open output file
    let file = File::create(output_path)?;
    let writer = BufWriter::new(file);
    let mut bam_writer = bam::io::Writer::new(writer);

    // Write header
    bam_writer.write_header(&header)?;

    // Generate records
    for i in 0..num_records {
        let mut record = bam::Record::default();

        // Read name
        let name = format!("read_{}", i);
        *record.name_mut() = name.as_bytes().into();

        // Reference sequence ID (chr1 = 0)
        *record.reference_sequence_id_mut() = Some(0);

        // Position (random within chr1)
        let pos = ((i * 12345) % 248956422) as i32 + 1;
        *record.alignment_start_mut() = Some(pos.try_into().unwrap());

        // MAPQ (random 0-60)
        *record.mapping_quality_mut() = Some(((i * 7) % 61) as u8);

        // Flags (properly paired, first in pair)
        *record.flags_mut() = sam::alignment::record::Flags::PAIRED
            | sam::alignment::record::Flags::PROPER_PAIR
            | sam::alignment::record::Flags::READ_1;

        // Sequence (100bp)
        let seq_bytes = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        *record.sequence_mut() = seq_bytes.to_vec().into();

        // Quality scores (all 'I' = Q40)
        let qual_bytes = vec![b'I'; 100];
        *record.quality_scores_mut() = qual_bytes.into();

        // CIGAR (100M)
        use sam::alignment::record::cigar::op::Kind;
        let cigar = vec![sam::alignment::record::cigar::Op::new(Kind::Match, 100)];
        *record.cigar_mut() = cigar.try_into().unwrap();

        bam_writer.write_record(&header, &record)?;

        if (i + 1) % 10000 == 0 {
            print!("\rGenerated {} / {} records...", i + 1, num_records);
            io::Write::flush(&mut io::stdout())?;
        }
    }

    println!("\rGenerated {} records successfully!", num_records);
    println!("Output: {}", output_path);

    Ok(())
}
