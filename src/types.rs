//! Common types used throughout biometal

/// A FASTQ record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    /// Sequence identifier (without '@' prefix)
    pub id: String,
    /// DNA/RNA sequence
    pub sequence: Vec<u8>,
    /// Quality scores (Phred+33)
    pub quality: Vec<u8>,
}

impl FastqRecord {
    /// Create a new FASTQ record
    pub fn new(id: String, sequence: Vec<u8>, quality: Vec<u8>) -> Self {
        Self { id, sequence, quality }
    }

    /// Check if the record has an empty sequence
    ///
    /// Returns `true` if the sequence length is zero. This can occur when
    /// quality-based trimming removes all bases (all below threshold).
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::FastqRecord;
    ///
    /// let empty = FastqRecord::new("read1".to_string(), Vec::new(), Vec::new());
    /// assert!(empty.is_empty());
    ///
    /// let non_empty = FastqRecord::new("read2".to_string(), b"ACGT".to_vec(), b"IIII".to_vec());
    /// assert!(!non_empty.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// A FASTA record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// Sequence identifier (without '>' prefix)
    pub id: String,
    /// DNA/RNA/protein sequence
    pub sequence: Vec<u8>,
}

impl FastaRecord {
    /// Create a new FASTA record
    pub fn new(id: String, sequence: Vec<u8>) -> Self {
        Self { id, sequence }
    }
}
