/// Minimal test to debug SeekableBgzfReader seeking bug

use std::fs::File;
use std::io::Read;

#[test]
fn test_seekable_bgzf_basic_seek() {
    use biometal::io::compression::SeekableBgzfReader;

    // Open the BAM file
    let file = File::open("tests/data/synthetic/alignment/synthetic_100k.bam")
        .expect("Failed to open BAM file");

    let mut reader = SeekableBgzfReader::new(file)
        .expect("Failed to create SeekableBgzfReader");

    // Read first few bytes to verify reader works
    let mut buf = vec![0u8; 100];
    let n = reader.read(&mut buf).expect("Failed to read initial data");
    println!("Read {} bytes from start", n);
    assert!(n > 0, "Should read some data");

    // Get the virtual offset we want to seek to
    // From the BAI index test, we know chunk 0 starts at 0x1160000
    // Compressed offset: 0x116 = 278 bytes
    // Uncompressed offset: 0x0000 = 0 bytes
    let virtual_offset: u64 = 0x1160000;

    println!("Attempting to seek to virtual offset: {:#x}", virtual_offset);
    println!("  Compressed offset: {:#x}", virtual_offset >> 16);
    println!("  Uncompressed offset: {:#x}", virtual_offset & 0xFFFF);

    // Try to seek - this is where it fails
    match reader.seek_to_virtual_offset(virtual_offset) {
        Ok(()) => {
            println!("✓ Seek succeeded!");

            // Try to read after seek
            let mut buf2 = vec![0u8; 100];
            match reader.read(&mut buf2) {
                Ok(n2) => println!("✓ Read {} bytes after seek", n2),
                Err(e) => println!("✗ Failed to read after seek: {}", e),
            }
        }
        Err(e) => {
            println!("✗ Seek failed: {}", e);
            panic!("Seek failed: {}", e);
        }
    }
}

#[test]
fn test_bgzf_block_at_offset_278() {
    // Let's manually check what's at offset 278 in the file
    let mut file = File::open("tests/data/synthetic/alignment/synthetic_100k.bam")
        .expect("Failed to open BAM file");

    use std::io::{Seek, SeekFrom};

    // Seek to byte 278
    file.seek(SeekFrom::Start(278)).expect("Failed to seek");

    // Read potential BGZF header (18 bytes)
    let mut header = [0u8; 18];
    file.read_exact(&mut header).expect("Failed to read header");

    println!("Bytes at offset 278:");
    println!("  Magic: {:02x} {:02x} (should be 1f 8b)", header[0], header[1]);
    println!("  Method: {:02x} (should be 08)", header[2]);
    println!("  Flags: {:02x} (should have 0x04 for FEXTRA)", header[3]);
    println!("  XLEN: {:02x} {:02x}", header[10], header[11]);

    let xlen = u16::from_le_bytes([header[10], header[11]]) as usize;
    println!("  XLEN value: {}", xlen);

    // Read extra field
    let mut extra = vec![0u8; xlen];
    file.read_exact(&mut extra).expect("Failed to read extra field");

    println!("  Extra field ({} bytes):", xlen);
    for (i, byte) in extra.iter().enumerate() {
        print!("{:02x} ", byte);
        if (i + 1) % 16 == 0 {
            println!();
        }
    }
    println!();

    // Check for BC subfield
    let mut pos = 0;
    while pos + 4 <= xlen {
        let si1 = extra[pos];
        let si2 = extra[pos + 1];
        let slen = u16::from_le_bytes([extra[pos + 2], extra[pos + 3]]) as usize;

        println!("  Subfield at pos {}: SI1={:02x} SI2={:02x} SLEN={}",
                 pos, si1, si2, slen);

        if si1 == 66 && si2 == 67 {
            println!("    ✓ Found BC subfield!");
            if pos + 6 <= xlen {
                let bsize = u16::from_le_bytes([extra[pos + 4], extra[pos + 5]]);
                println!("    BSIZE = {} (block size = {})", bsize, bsize as usize + 1);
            }
        }

        pos += 4 + slen;
    }
}
