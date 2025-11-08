# Phase 1: Minimal BAM Reader - Kickoff

**Date**: November 8, 2025
**Phase**: Phase 1 (Weeks 2-3)
**Status**: Ready to begin
**Decision**: GO ✅ (Phase 0 complete)

---

## Phase 0 Recap

✅ **Profiling complete**: BGZF decompression is bottleneck (66-80% CPU time)
✅ **Decision made**: GO for native implementation (full-stack ownership)
✅ **Strategy defined**: Parallel BGZF first (captures main bottleneck), native parsing second

---

## Phase 1 Goals (Weeks 2-3)

**Primary Goal**: Build minimal, correct BAM parser

**Success Criteria**:
- ✅ Parse BAM header (SAM header + reference sequences)
- ✅ Parse BAM records (all fields: name, position, MAPQ, FLAGS, sequence, quality, CIGAR, tags)
- ✅ Streaming iterator interface
- ✅ Differential testing passes (vs noodles)

**Non-Goals (defer to later phases)**:
- ❌ Performance optimization (Phase 2-3)
- ❌ NEON intrinsics (Phase 4, if needed)
- ❌ SAM format support (Phase 5)
- ❌ Complete error handling (Phase 6)

**Focus**: **Correctness over speed**

---

## Architecture Overview

### High-Level Structure

```
biometal/
└── src/
    └── io/
        └── bam/
            ├── mod.rs              # Public API
            ├── header.rs           # Header parsing
            ├── record.rs           # Record structure
            ├── reader.rs           # Streaming reader
            ├── sequence.rs         # 4-bit sequence decoding
            ├── quality.rs          # Quality score parsing
            ├── cigar.rs            # CIGAR string parsing
            ├── tags.rs             # Optional tags
            └── compression.rs      # BGZF wrapper (use existing biometal)
```

### Key Components

**1. BamReader** (public API):
```rust
pub struct BamReader<R> {
    bgzf: BgzfReader<R>,
    header: Header,
    buffer: Vec<u8>,
}

impl<R: BufRead> BamReader<R> {
    pub fn new(reader: R) -> io::Result<Self>;
    pub fn read_header(&mut self) -> io::Result<Header>;
    pub fn records(&mut self) -> Records<'_, R>;
}

pub struct Records<'a, R> {
    reader: &'a mut BamReader<R>,
}

impl<R: BufRead> Iterator for Records<'_, R> {
    type Item = io::Result<Record>;
    fn next(&mut self) -> Option<Self::Item>;
}
```

**2. Header** (parsed from BAM):
```rust
pub struct Header {
    pub text: String,           // SAM header text
    pub references: Vec<Reference>,
}

pub struct Reference {
    pub name: String,
    pub length: u32,
}
```

**3. Record** (BAM alignment):
```rust
pub struct Record {
    pub name: String,
    pub reference_id: Option<usize>,
    pub position: Option<i32>,
    pub mapq: Option<u8>,
    pub flags: u16,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub cigar: Vec<CigarOp>,
    pub tags: Tags,
    // ... mate info, template length, etc.
}
```

---

## Implementation Plan

### Week 2: Basic Structure

**Day 1-2: Project Setup**
```bash
# Create BAM module structure
mkdir -p src/io/bam
touch src/io/bam/{mod.rs,header.rs,record.rs,reader.rs}

# Update src/io/mod.rs
# pub mod bam;

# Add dependencies to Cargo.toml
# (none yet - pure Rust, use existing biometal BGZF later)
```

**Day 3-4: Header Parsing**
```rust
// src/io/bam/header.rs

/// Parse BAM magic bytes ("BAM\1")
fn read_magic(reader: &mut impl Read) -> io::Result<()> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != b"BAM\1" {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid BAM magic"));
    }
    Ok(())
}

/// Parse SAM header text (null-terminated string)
fn read_header_text(reader: &mut impl Read) -> io::Result<String> {
    let mut len = [0u8; 4];
    reader.read_exact(&mut len)?;
    let len = i32::from_le_bytes(len) as usize;

    let mut text = vec![0u8; len];
    reader.read_exact(&mut text)?;

    String::from_utf8(text)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Parse reference sequences
fn read_references(reader: &mut impl Read) -> io::Result<Vec<Reference>> {
    let mut count = [0u8; 4];
    reader.read_exact(&mut count)?;
    let count = i32::from_le_bytes(count) as usize;

    let mut refs = Vec::with_capacity(count);
    for _ in 0..count {
        // Read name length, name, sequence length
        let name = read_reference_name(reader)?;
        let length = read_reference_length(reader)?;
        refs.push(Reference { name, length });
    }
    Ok(refs)
}
```

**Day 5-7: Record Parsing**
```rust
// src/io/bam/record.rs

/// Parse a single BAM record (binary format)
pub fn parse_record(data: &[u8]) -> io::Result<Record> {
    let mut cursor = 0;

    // Block size (skip, already know from read)
    cursor += 4;

    // Reference ID
    let ref_id = read_i32(&data[cursor..])?;
    cursor += 4;

    // Position
    let pos = read_i32(&data[cursor..])?;
    cursor += 4;

    // Read name length
    let l_read_name = data[cursor] as usize;
    cursor += 1;

    // MAPQ
    let mapq = data[cursor];
    cursor += 1;

    // BAI bin (skip for now)
    cursor += 2;

    // CIGAR length
    let n_cigar_op = read_u16(&data[cursor..])? as usize;
    cursor += 2;

    // FLAGS
    let flags = read_u16(&data[cursor..])?;
    cursor += 2;

    // Sequence length
    let l_seq = read_i32(&data[cursor..])? as usize;
    cursor += 4;

    // Next reference ID (mate)
    let next_ref_id = read_i32(&data[cursor..])?;
    cursor += 4;

    // Next position (mate)
    let next_pos = read_i32(&data[cursor..])?;
    cursor += 4;

    // Template length
    let tlen = read_i32(&data[cursor..])?;
    cursor += 4;

    // Read name (null-terminated)
    let name = read_cstring(&data[cursor..], l_read_name)?;
    cursor += l_read_name;

    // CIGAR operations
    let cigar = read_cigar(&data[cursor..], n_cigar_op)?;
    cursor += n_cigar_op * 4;

    // Sequence (4-bit encoded, 2 bases per byte)
    let sequence = decode_sequence(&data[cursor..], l_seq)?;
    cursor += (l_seq + 1) / 2;

    // Quality scores
    let quality = read_quality(&data[cursor..], l_seq)?;
    cursor += l_seq;

    // Optional tags (remaining bytes)
    let tags = parse_tags(&data[cursor..])?;

    Ok(Record {
        name,
        reference_id: if ref_id >= 0 { Some(ref_id as usize) } else { None },
        position: if pos >= 0 { Some(pos) } else { None },
        mapq: if mapq != 255 { Some(mapq) } else { None },
        flags,
        sequence,
        quality,
        cigar,
        tags,
        // ... mate info, tlen
    })
}
```

### Week 3: Integration & Testing

**Day 1-2: Streaming Iterator**
```rust
// src/io/bam/reader.rs

impl<R: BufRead> Iterator for Records<'_, R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        // Read record block size
        let block_size = match self.reader.read_block_size() {
            Ok(0) => return None,  // EOF
            Ok(size) => size,
            Err(e) => return Some(Err(e)),
        };

        // Read record data
        self.reader.buffer.resize(block_size, 0);
        if let Err(e) = self.reader.bgzf.read_exact(&mut self.reader.buffer) {
            return Some(Err(e));
        }

        // Parse record
        Some(parse_record(&self.reader.buffer))
    }
}
```

**Day 3-4: Differential Testing**
```rust
// tests/differential_bam.rs

#[test]
fn test_against_noodles() {
    let bam_path = "test-data/synthetic_100000.bam";

    // Parse with biometal
    let mut biometal_reader = biometal::io::bam::BamReader::from_path(bam_path).unwrap();
    let biometal_records: Vec<_> = biometal_reader.records()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();

    // Parse with noodles
    let mut noodles_reader = noodles_bam::io::Reader::from_path(bam_path).unwrap();
    let noodles_records: Vec<_> = /* ... collect records ... */;

    // Compare record-by-record
    assert_eq!(biometal_records.len(), noodles_records.len());
    for (i, (br, nr)) in biometal_records.iter().zip(noodles_records.iter()).enumerate() {
        assert_eq!(br.name, nr.name, "Record {}: name mismatch", i);
        assert_eq!(br.position, nr.position, "Record {}: position mismatch", i);
        assert_eq!(br.sequence, nr.sequence, "Record {}: sequence mismatch", i);
        // ... compare all fields
    }
}
```

**Day 5-7: Bug Fixes & Polish**
- Fix any differential testing failures
- Handle edge cases (unmapped reads, secondary alignments)
- Basic error handling
- Documentation

---

## Critical Implementation Details

### 1. 4-Bit Sequence Decoding

**Format**: 2 bases per byte, high nibble first

```rust
// src/io/bam/sequence.rs

const SEQ_LOOKUP: [u8; 16] = [
    b'=', b'A', b'C', b'M',
    b'G', b'R', b'S', b'V',
    b'T', b'W', b'Y', b'H',
    b'K', b'D', b'B', b'N',
];

pub fn decode_sequence(data: &[u8], length: usize) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(length);

    for i in 0..length {
        let byte_idx = i / 2;
        let nibble = if i % 2 == 0 {
            data[byte_idx] >> 4  // High nibble
        } else {
            data[byte_idx] & 0x0F  // Low nibble
        };
        sequence.push(SEQ_LOOKUP[nibble as usize]);
    }

    sequence
}
```

**Note**: No NEON optimization yet (Phase 4 if needed)

### 2. CIGAR Parsing

**Format**: 32-bit integer per operation (4 bits op, 28 bits length)

```rust
// src/io/bam/cigar.rs

pub enum CigarOp {
    Match(u32),       // M
    Insertion(u32),   // I
    Deletion(u32),    // D
    RefSkip(u32),     // N
    SoftClip(u32),    // S
    HardClip(u32),    // H
    Padding(u32),     // P
    SeqMatch(u32),    // =
    SeqMismatch(u32), // X
}

pub fn parse_cigar(data: &[u8], n_ops: usize) -> io::Result<Vec<CigarOp>> {
    let mut ops = Vec::with_capacity(n_ops);

    for i in 0..n_ops {
        let cigar_int = u32::from_le_bytes([
            data[i*4], data[i*4+1], data[i*4+2], data[i*4+3]
        ]);

        let length = cigar_int >> 4;
        let op = cigar_int & 0x0F;

        let cigar_op = match op {
            0 => CigarOp::Match(length),
            1 => CigarOp::Insertion(length),
            2 => CigarOp::Deletion(length),
            3 => CigarOp::RefSkip(length),
            4 => CigarOp::SoftClip(length),
            5 => CigarOp::HardClip(length),
            6 => CigarOp::Padding(length),
            7 => CigarOp::SeqMatch(length),
            8 => CigarOp::SeqMismatch(length),
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid CIGAR op")),
        };

        ops.push(cigar_op);
    }

    Ok(ops)
}
```

### 3. Optional Tags

**Format**: TAG:TYPE:VALUE (variable length)

```rust
// src/io/bam/tags.rs

pub struct Tags {
    data: Vec<u8>,  // Raw tag data for now
}

pub fn parse_tags(data: &[u8]) -> io::Result<Tags> {
    // Phase 1: Just store raw bytes
    // Phase 6: Parse individual tags
    Ok(Tags {
        data: data.to_vec(),
    })
}
```

**Note**: Full tag parsing deferred to Phase 6

---

## Testing Strategy

### Unit Tests

**Per-module tests**:
- `header.rs`: Test header parsing with known inputs
- `sequence.rs`: Test 4-bit decoding with all 16 bases
- `cigar.rs`: Test all CIGAR operations
- `tags.rs`: Test tag parsing (later)

### Differential Tests

**Against noodles**:
- Parse same BAM file with both biometal and noodles
- Compare every field of every record
- Use test data: `test-data/synthetic_100000.bam`

**Property tests** (later):
- Round-trip: write → read → compare
- Edge cases: empty sequences, unmapped reads, etc.

### Test Data

**Existing**:
- ✅ `synthetic_100000.bam` (100K records, generated in Phase 0)

**Additional** (generate in Phase 1):
- Small: 100 records (quick iteration)
- Edge cases: unmapped reads, secondary alignments, etc.

---

## Success Criteria

**Phase 1 Complete when**:
- ✅ BAM header parsing works
- ✅ BAM record parsing works (all fields)
- ✅ Streaming iterator works (constant memory)
- ✅ Differential tests pass (100% match vs noodles)
- ✅ Code compiles, tests pass
- ✅ Basic documentation exists

**Quality Bar**:
- Correctness: 100% differential test pass
- Memory: Streaming (no accumulation)
- Performance: Not measured yet (Phase 2)

---

## Dependencies

**Cargo.toml additions**:
```toml
# For differential testing only
[dev-dependencies]
noodles-bam = "0.68"
noodles-sam = "0.68"

# No additional runtime dependencies
# (will use existing biometal BGZF in Phase 2)
```

---

## Timeline

**Week 2**: Basic parsing
- Days 1-2: Project setup
- Days 3-4: Header parsing
- Days 5-7: Record parsing

**Week 3**: Integration & testing
- Days 1-2: Streaming iterator
- Days 3-4: Differential testing
- Days 5-7: Bug fixes, edge cases

**Deliverable**: Minimal, correct BAM parser (no optimization yet)

---

## Next Phase Preview

**Phase 2** (Weeks 3-4): Parallel BGZF Integration
- Replace sequential BGZF with biometal parallel version
- Benchmark: baseline vs parallel
- Expected: ~6.5× BGZF speedup = ~4-5× overall
- This is where the BIG performance win happens!

---

## Getting Started

**Immediate next steps**:

1. **Create module structure**:
   ```bash
   mkdir -p src/io/bam
   touch src/io/bam/{mod.rs,header.rs,record.rs,reader.rs,sequence.rs,cigar.rs,tags.rs}
   ```

2. **Start with header parsing** (simplest):
   - Read magic bytes ("BAM\1")
   - Parse SAM header text
   - Parse reference sequences

3. **Write unit tests as you go**:
   - Test each function with known inputs
   - Build confidence incrementally

4. **Iterate quickly**:
   - Focus on correctness
   - Don't worry about performance yet
   - Use differential testing to validate

---

**Phase 1**: Minimal BAM Reader
**Focus**: Correctness over speed
**Timeline**: Weeks 2-3
**Status**: Ready to begin! 🚀

Let's build a native, ARM-optimized BAM reader for biometal!
