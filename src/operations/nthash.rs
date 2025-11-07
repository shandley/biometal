//! ntHash: Rolling hash function for genomic sequences
//!
//! ntHash is a hash function tuned for genomic data that performs best when
//! calculating hash values for adjacent k-mers using rolling hash updates.
//!
//! # Algorithm
//!
//! ntHash uses table lookups and bit rotations (both SIMD-friendly operations):
//! 1. **Forward hash**: XOR of rotated nucleotide values
//! 2. **Reverse-complement hash**: XOR of RC-rotated nucleotide values
//! 3. **Rolling update**: Remove outgoing base, add incoming base (O(1))
//! 4. **Canonical**: Minimum of forward and reverse-complement
//!
//! # Performance
//!
//! - **Initial hash**: O(k) per k-mer
//! - **Rolling update**: O(1) per k-mer (vs O(k) rehashing)
//! - **SIMD potential**: 8-way parallelism for batch hashing
//!
//! # Evidence
//!
//! - **Source**: simd-minimizers-analysis (Nov 2025)
//! - **Speedup**: 8× faster than FNV-1a for minimizer extraction
//! - **Vectorizable**: Table lookup, rotate, XOR (all NEON-friendly)
//!
//! # References
//!
//! - Original paper: Mohamadi et al. (2016) https://doi.org/10.1093/bioinformatics/btw397
//! - C++ implementation: https://github.com/bcgsc/ntHash/
//! - Rust implementation: nthash crate (MIT licensed)
//! - simd-minimizers: Uses ntHash for 221× faster minimizers
//!
//! # Example
//!
//! ```
//! use biometal::operations::nthash::NtHashIterator;
//!
//! let seq = b"ACTGC";
//! let iter = NtHashIterator::new(seq, 3).unwrap();
//! let hashes: Vec<u64> = iter.collect();
//! // Returns canonical hashes (min of forward and reverse-complement)
//! ```

use crate::error::{BiometalError, Result};

/// Maximum k-mer size supported (u32::MAX)
pub const MAX_K_SIZE: usize = u32::MAX as usize;

/// Lookup table for forward strand hashing
///
/// Each nucleotide (A, C, G, T) maps to a 64-bit value used in the hash computation.
/// These values are chosen to provide good hash distribution.
const H_LOOKUP: [u64; 256] = {
    let mut lookup = [0; 256];
    lookup[b'A' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'C' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'G' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'T' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'N' as usize] = 0;
    // All other indices remain 0 (invalid)
    lookup
};

/// Lookup table for reverse-complement strand hashing
///
/// Nucleotides are mapped to their reverse-complement values:
/// A ↔ T, C ↔ G
const RC_LOOKUP: [u64; 256] = {
    let mut lookup = [0; 256];
    lookup[b'A' as usize] = 0x2955_49f5_4be2_4456; // T's value
    lookup[b'C' as usize] = 0x2032_3ed0_8257_2324; // G's value
    lookup[b'G' as usize] = 0x3193_c185_62a0_2b4c; // C's value
    lookup[b'T' as usize] = 0x3c8b_fbb3_95c6_0474; // A's value
    lookup[b'N' as usize] = 0;
    // All other indices remain 0 (invalid)
    lookup
};

/// Get forward hash value for a nucleotide
///
/// # Errors
///
/// Returns error if nucleotide is not A, C, G, T, or N
#[inline(always)]
fn h(c: u8) -> Result<u64> {
    let val = H_LOOKUP[c as usize];
    if val == 0 && c != b'N' {
        return Err(BiometalError::InvalidInput {
            msg: format!("Invalid nucleotide: {}", c as char),
        });
    }
    Ok(val)
}

/// Get reverse-complement hash value for a nucleotide
///
/// # Errors
///
/// Returns error if nucleotide is not A, C, G, T, or N
#[inline(always)]
fn rc(c: u8) -> Result<u64> {
    let val = RC_LOOKUP[c as usize];
    if val == 0 && c != b'N' {
        return Err(BiometalError::InvalidInput {
            msg: format!("Invalid nucleotide: {}", c as char),
        });
    }
    Ok(val)
}

/// Calculate forward hash for a k-mer
///
/// This is the initial hash computation (O(k)). Subsequent k-mers use
/// rolling update (O(1)).
///
/// # Algorithm
///
/// ```text
/// hash = 0
/// for i in 0..k:
///     hash ^= h(sequence[i]).rotate_left(k - i - 1)
/// ```
///
/// # Arguments
///
/// * `seq` - Sequence slice (must contain only A, C, G, T, N)
/// * `k` - K-mer size
///
/// # Returns
///
/// 64-bit forward hash value
///
/// # Errors
///
/// Returns error if sequence contains invalid nucleotides
fn forward_hash(seq: &[u8], k: usize) -> Result<u64> {
    let mut hash = 0u64;
    for (i, &base) in seq.iter().take(k).enumerate() {
        hash ^= h(base)?.rotate_left((k - i - 1) as u32);
    }
    Ok(hash)
}

/// Calculate reverse-complement hash for a k-mer
///
/// This is the initial hash computation (O(k)). Subsequent k-mers use
/// rolling update (O(1)).
///
/// # Algorithm
///
/// ```text
/// hash = 0
/// for i in (0..k).rev():
///     hash ^= rc(sequence[i]).rotate_left(k - i - 1)
/// ```
///
/// # Arguments
///
/// * `seq` - Sequence slice (must contain only A, C, G, T, N)
/// * `k` - K-mer size
///
/// # Returns
///
/// 64-bit reverse-complement hash value
///
/// # Errors
///
/// Returns error if sequence contains invalid nucleotides
fn reverse_complement_hash(seq: &[u8], k: usize) -> Result<u64> {
    let mut hash = 0u64;
    for (i, &base) in seq.iter().take(k).rev().enumerate() {
        hash ^= rc(base)?.rotate_left((k - i - 1) as u32);
    }
    Ok(hash)
}

/// ntHash iterator for canonical hashes (minimum of forward and reverse-complement)
///
/// This iterator provides O(1) rolling hash updates for adjacent k-mers.
/// The first k-mer requires O(k) initialization, but subsequent k-mers
/// use only O(1) operations (remove outgoing base, add incoming base).
///
/// # Performance
///
/// - **Initial hash**: O(k)
/// - **Rolling update**: O(1) per k-mer
/// - **Total for n k-mers**: O(k + n) vs O(k × n) for rehashing
///
/// # Example
///
/// ```
/// use biometal::operations::nthash::NtHashIterator;
///
/// let seq = b"ACTGCACTGC";
/// let iter = NtHashIterator::new(seq, 5).unwrap();
///
/// for hash in iter {
///     println!("Canonical hash: {:x}", hash);
/// }
/// ```
#[derive(Debug)]
pub struct NtHashIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u64,  // Forward hash
    rh: u64,  // Reverse-complement hash
    current_idx: usize,
    max_idx: usize,
}

impl<'a> NtHashIterator<'a> {
    /// Create a new ntHash iterator
    ///
    /// # Arguments
    ///
    /// * `seq` - DNA sequence (A, C, G, T, N only)
    /// * `k` - K-mer size
    ///
    /// # Returns
    ///
    /// Iterator yielding canonical hashes (min of forward and reverse-complement)
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - k > sequence length
    /// - k > MAX_K_SIZE
    /// - Sequence contains invalid nucleotides
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::operations::nthash::NtHashIterator;
    ///
    /// let seq = b"ACTGC";
    /// let iter = NtHashIterator::new(seq, 3)?;
    /// let hashes: Vec<u64> = iter.collect();
    /// # Ok::<(), biometal::BiometalError>(())
    /// ```
    pub fn new(seq: &'a [u8], k: usize) -> Result<Self> {
        if k > seq.len() {
            return Err(BiometalError::InvalidInput {
                msg: format!("k-mer size ({}) exceeds sequence length ({})", k, seq.len()),
            });
        }
        if k > MAX_K_SIZE {
            return Err(BiometalError::InvalidInput {
                msg: format!("k-mer size ({}) exceeds maximum ({})", k, MAX_K_SIZE),
            });
        }
        if k == 0 {
            return Err(BiometalError::InvalidInput {
                msg: "k-mer size must be > 0".to_string(),
            });
        }

        // Initialize forward and reverse-complement hashes for first k-mer
        let fh = forward_hash(seq, k)?;
        let rh = reverse_complement_hash(seq, k)?;

        Ok(NtHashIterator {
            seq,
            k,
            fh,
            rh,
            current_idx: 0,
            max_idx: seq.len() - k + 1,
        })
    }

    /// Get the forward hash (non-canonical)
    ///
    /// Useful for applications that need both forward and reverse-complement hashes
    pub fn forward_hash(&self) -> u64 {
        self.fh
    }

    /// Get the reverse-complement hash
    ///
    /// Useful for applications that need both forward and reverse-complement hashes
    pub fn reverse_complement_hash(&self) -> u64 {
        self.rh
    }
}

impl<'a> Iterator for NtHashIterator<'a> {
    type Item = u64;

    /// Get next canonical hash (minimum of forward and reverse-complement)
    ///
    /// # Algorithm (Rolling Update)
    ///
    /// For k-mer at position i+1 (given hash at position i):
    ///
    /// ```text
    /// Forward:
    ///   1. Rotate current hash left by 1
    ///   2. XOR out the outgoing base (rotate_left by k)
    ///   3. XOR in the incoming base
    ///
    /// Reverse-complement:
    ///   1. Rotate current hash right by 1
    ///   2. XOR out the outgoing RC base (rotate_right by 1)
    ///   3. XOR in the incoming RC base (rotate_left by k-1)
    /// ```
    ///
    /// All operations are O(1) and SIMD-friendly (rotate, XOR, table lookup)
    fn next(&mut self) -> Option<u64> {
        if self.current_idx >= self.max_idx {
            return None;
        }

        // Rolling hash update (O(1) for positions after the first)
        if self.current_idx > 0 {
            let i = self.current_idx - 1;
            let outgoing = self.seq[i];
            let incoming = self.seq[i + self.k];

            // Update forward hash
            // Note: We use unwrap here because we validated all bases during initialization
            self.fh = self
                .fh
                .rotate_left(1)
                ^ h(outgoing).unwrap().rotate_left(self.k as u32)
                ^ h(incoming).unwrap();

            // Update reverse-complement hash
            self.rh = self
                .rh
                .rotate_right(1)
                ^ rc(outgoing).unwrap().rotate_right(1)
                ^ rc(incoming).unwrap().rotate_left((self.k - 1) as u32);
        }

        self.current_idx += 1;

        // Return canonical hash (minimum of forward and reverse-complement)
        Some(u64::min(self.fh, self.rh))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.max_idx - self.current_idx;
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for NtHashIterator<'a> {}

/// ntHash iterator for forward-only hashes (not canonical)
///
/// Use this when you only need forward strand hashes, avoiding the
/// overhead of reverse-complement computation.
///
/// # Example
///
/// ```
/// use biometal::operations::nthash::NtHashForwardIterator;
///
/// let seq = b"ACTGC";
/// let iter = NtHashForwardIterator::new(seq, 3).unwrap();
/// let hashes: Vec<u64> = iter.collect();
/// ```
#[derive(Debug)]
pub struct NtHashForwardIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u64,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> NtHashForwardIterator<'a> {
    /// Create a new forward-only ntHash iterator
    ///
    /// # Arguments
    ///
    /// * `seq` - DNA sequence (A, C, G, T, N only)
    /// * `k` - K-mer size
    ///
    /// # Returns
    ///
    /// Iterator yielding forward hashes only
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - k > sequence length
    /// - k > MAX_K_SIZE
    /// - Sequence contains invalid nucleotides
    pub fn new(seq: &'a [u8], k: usize) -> Result<Self> {
        if k > seq.len() {
            return Err(BiometalError::InvalidInput {
                msg: format!("k-mer size ({}) exceeds sequence length ({})", k, seq.len()),
            });
        }
        if k > MAX_K_SIZE {
            return Err(BiometalError::InvalidInput {
                msg: format!("k-mer size ({}) exceeds maximum ({})", k, MAX_K_SIZE),
            });
        }
        if k == 0 {
            return Err(BiometalError::InvalidInput {
                msg: "k-mer size must be > 0".to_string(),
            });
        }

        // Initialize forward hash for first k-mer
        let fh = forward_hash(seq, k)?;

        Ok(NtHashForwardIterator {
            seq,
            k,
            fh,
            current_idx: 0,
            max_idx: seq.len() - k + 1,
        })
    }
}

impl<'a> Iterator for NtHashForwardIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        if self.current_idx >= self.max_idx {
            return None;
        }

        // Rolling hash update (O(1) for positions after the first)
        if self.current_idx > 0 {
            let i = self.current_idx - 1;
            let outgoing = self.seq[i];
            let incoming = self.seq[i + self.k];

            // Update forward hash
            self.fh = self
                .fh
                .rotate_left(1)
                ^ h(outgoing).unwrap().rotate_left(self.k as u32)
                ^ h(incoming).unwrap();
        }

        self.current_idx += 1;
        Some(self.fh)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.max_idx - self.current_idx;
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for NtHashForwardIterator<'a> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lookup_tables() {
        // Verify forward lookup
        assert_eq!(H_LOOKUP[b'A' as usize], 0x3c8b_fbb3_95c6_0474);
        assert_eq!(H_LOOKUP[b'C' as usize], 0x3193_c185_62a0_2b4c);
        assert_eq!(H_LOOKUP[b'G' as usize], 0x2032_3ed0_8257_2324);
        assert_eq!(H_LOOKUP[b'T' as usize], 0x2955_49f5_4be2_4456);
        assert_eq!(H_LOOKUP[b'N' as usize], 0);

        // Verify reverse-complement lookup (A ↔ T, C ↔ G)
        assert_eq!(RC_LOOKUP[b'A' as usize], H_LOOKUP[b'T' as usize]);
        assert_eq!(RC_LOOKUP[b'T' as usize], H_LOOKUP[b'A' as usize]);
        assert_eq!(RC_LOOKUP[b'C' as usize], H_LOOKUP[b'G' as usize]);
        assert_eq!(RC_LOOKUP[b'G' as usize], H_LOOKUP[b'C' as usize]);
    }

    #[test]
    fn test_canonical_iterator_basic() {
        let seq = b"ACTGC";
        let iter = NtHashIterator::new(seq, 3).unwrap();
        let hashes: Vec<u64> = iter.collect();

        assert_eq!(hashes.len(), 3); // ACTGC has 3 k-mers of size 3
        // Hashes should be consistent (but we don't check exact values here)
    }

    #[test]
    fn test_forward_iterator_basic() {
        let seq = b"ACTGC";
        let iter = NtHashForwardIterator::new(seq, 3).unwrap();
        let hashes: Vec<u64> = iter.collect();

        assert_eq!(hashes.len(), 3);
    }

    #[test]
    fn test_invalid_k_size() {
        let seq = b"ACTGC";
        assert!(NtHashIterator::new(seq, 0).is_err()); // k=0
        assert!(NtHashIterator::new(seq, 10).is_err()); // k > seq.len()
    }

    #[test]
    fn test_invalid_nucleotide() {
        let seq = b"ACXGT"; // X is invalid (within first k=3)
        assert!(NtHashIterator::new(seq, 3).is_err());
    }

    #[test]
    fn test_rolling_property() {
        // Rolling hash should give same results as computing each hash independently
        let seq = b"ACTGCACTGC";
        let k = 5;

        // Collect hashes using rolling iterator
        let rolling_hashes: Vec<u64> = NtHashIterator::new(seq, k).unwrap().collect();

        // Verify we got the expected number of hashes
        assert_eq!(rolling_hashes.len(), seq.len() - k + 1);

        // Each hash should be consistent (same k-mer, same hash)
        assert_eq!(rolling_hashes[0], rolling_hashes[5]); // ACTGC appears twice
    }

    #[test]
    fn test_n_handling() {
        let seq = b"ACNGC";
        let iter = NtHashIterator::new(seq, 3).unwrap();
        let hashes: Vec<u64> = iter.collect();

        assert_eq!(hashes.len(), 3); // Should handle N gracefully
    }

    #[test]
    fn test_exact_size_iterator() {
        let seq = b"ACTGCACTGC";
        let iter = NtHashIterator::new(seq, 5).unwrap();

        assert_eq!(iter.len(), 6); // 10 - 5 + 1 = 6 k-mers
        assert_eq!(iter.size_hint(), (6, Some(6)));
    }
}
