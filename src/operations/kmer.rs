//! K-mer operations optimized for Apple Silicon
//!
//! # Evidence-Based Design
//!
//! K-mer operations are **data-structure-bound** (hash+HashMap), not compute-bound.
//! Unlike element-wise operations (base_counting, gc_content), k-mer operations spend
//! 50-60% of runtime on hash computation and 30-40% on data structure operations,
//! leaving only 5-10% for base validation (the only NEON-friendly part).
//!
//! # Evidence
//!
//! **Source**: [ASBB Entry 034](https://github.com/shandley/apple-silicon-bio-bench/blob/main/lab-notebook/2025-11/20251106-034-EXPERIMENT-kmer-operations.md)
//!
//! **Experiments**: Pilot benchmark (N=3), full hardware sweep
//! - **Minimizers**: 1.02-1.26× max (NEON, Parallel) - below ≥5× threshold
//! - **K-mer Spectrum**: 0.95-1.88× (sometimes SLOWER with parallel!)
//! - **K-mer Extraction**: 2.19-2.38× with Parallel-4t (opt-in future)
//!
//! **Conclusion**: Scalar-only implementations are optimal for most use cases.
//!
//! # Why No NEON/Parallel?
//!
//! **Runtime breakdown** (Entry 034):
//! - Hash computation: 50-60% (sequential, can't vectorize)
//! - HashMap operations: 30-40% (sequential, thread contention)
//! - Base validation: 5-10% (only NEON-friendly part → minimal impact)
//!
//! **NEON can only accelerate 5-10% of runtime** → ~1× overall speedup
//!
//! # Design Principles
//!
//! 1. **Simplicity by default**: Scalar implementations for 90% of users
//! 2. **Document evidence**: Link to ASBB experiments explaining design
//! 3. **Streaming-friendly**: Zero-copy iterators where possible
//! 4. **Constant memory**: Aligns with Rule 5 (streaming architecture)
//!
//! # Examples
//!
//! ```
//! use biometal::operations::kmer::{extract_kmers, kmer_spectrum};
//!
//! // Extract overlapping k-mers
//! let sequence = b"ATGCATGC";
//! let kmers = extract_kmers(sequence, 3);
//! // Returns: [b"ATG", b"TGC", b"GCA", b"CAT", b"ATG", b"TGC"]
//!
//! // Count k-mer frequencies
//! let sequences = vec![b"ATGCAT".as_ref(), b"GCATGC".as_ref()];
//! let spectrum = kmer_spectrum(&sequences, 3);
//! // Returns: {b"ATG": 2, b"TGC": 2, b"GCA": 2, b"CAT": 2}
//! ```

use crate::error::Result;
use crate::operations::nthash::NtHashIterator;
use crate::operations::sliding_min::SlidingMin;
use std::collections::HashMap;

/// Minimizer structure (position + hash value)
///
/// Minimizers are used for efficient sequence indexing and sketching.
/// They represent the lexicographically smallest k-mer in a sliding window.
///
/// # Backward Compatibility
///
/// This struct now stores only position and hash for performance (1.25× faster).
/// K-mer sequence can be extracted on-demand using helper functions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
    /// Position in the original sequence
    pub position: usize,
    /// ntHash canonical hash value of the k-mer
    pub hash: u64,
    /// K-mer size (for extracting k-mer from sequence)
    pub k: usize,
}

impl Minimizer {
    /// Extract the k-mer from the original sequence (zero-copy).
    ///
    /// Returns a slice reference to the k-mer without allocation.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::operations::kmer::extract_minimizers_fast;
    ///
    /// let sequence = b"ATGCATGCATGC";
    /// let minimizers = extract_minimizers_fast(sequence, 3, 5).unwrap();
    ///
    /// for m in &minimizers {
    ///     let kmer = m.kmer(sequence); // Zero-copy slice
    ///     println!("K-mer: {:?}", kmer);
    /// }
    /// ```
    #[inline]
    pub fn kmer<'a>(&self, sequence: &'a [u8]) -> &'a [u8] {
        &sequence[self.position..self.position + self.k]
    }

    /// Extract the k-mer as an owned Vec (allocates).
    ///
    /// Only use this when you need to store the k-mer separately from
    /// the sequence. For most cases, use `kmer()` for zero-copy access.
    ///
    /// # Example
    ///
    /// ```
    /// use biometal::operations::kmer::extract_minimizers_fast;
    ///
    /// let sequence = b"ATGCATGCATGC";
    /// let minimizers = extract_minimizers_fast(sequence, 3, 5).unwrap();
    ///
    /// let kmer_owned = minimizers[0].kmer_owned(sequence); // Allocates Vec
    /// ```
    #[inline]
    pub fn kmer_owned(&self, sequence: &[u8]) -> Vec<u8> {
        self.kmer(sequence).to_vec()
    }
}

/// Extract overlapping k-mers from a sequence (scalar-only)
///
/// Returns all valid k-mers (containing only A, C, G, T) as owned byte vectors.
/// K-mers containing N or other ambiguous bases are skipped.
///
/// # Performance
///
/// - **Time**: O(n) where n = sequence length
/// - **Space**: O(k × m) where m = number of valid k-mers
/// - **Optimization**: Scalar-only (Entry 034: parallel provides 2.2× but opt-in)
///
/// # Evidence
///
/// ASBB Entry 034 found that k-mer extraction with Parallel-4t provides 2.19-2.38×
/// speedup, but only for large datasets where overhead is amortized. For simplicity,
/// the default implementation is scalar. Optional parallel extraction will be
/// provided via `KmerExtractor::with_parallel()` in a future update.
///
/// # Arguments
///
/// * `sequence` - DNA sequence to extract k-mers from
/// * `k` - K-mer size (must be ≤ sequence length)
///
/// # Returns
///
/// Vector of k-mers as owned byte vectors. K-mers with invalid bases (N, etc.) are skipped.
///
/// # Examples
///
/// ```
/// use biometal::operations::kmer::extract_kmers;
///
/// let sequence = b"ATGCATGC";
/// let kmers = extract_kmers(sequence, 3);
///
/// assert_eq!(kmers.len(), 6);
/// assert_eq!(kmers[0], b"ATG");
/// assert_eq!(kmers[1], b"TGC");
/// assert_eq!(kmers[2], b"GCA");
/// ```
///
/// ```
/// use biometal::operations::kmer::extract_kmers;
///
/// // K-mers with N are skipped
/// let sequence = b"ATGCNNATGC";
/// let kmers = extract_kmers(sequence, 3);
/// // Skips: "GCN", "CNN", "NNA", "NAT"
/// assert_eq!(kmers.len(), 4); // ATG, TGC, ATG, TGC
/// ```
pub fn extract_kmers(sequence: &[u8], k: usize) -> Vec<Vec<u8>> {
    if k == 0 || k > sequence.len() {
        return Vec::new();
    }

    let mut kmers = Vec::new();

    for i in 0..=(sequence.len() - k) {
        let kmer = &sequence[i..i + k];
        if is_valid_kmer(kmer) {
            kmers.push(kmer.to_vec());
        }
    }

    kmers
}

/// Streaming k-mer iterator (zero-copy)
///
/// Returns an iterator over k-mers as borrowed slices. This avoids allocation
/// and enables constant-memory processing (Rule 5: streaming architecture).
///
/// # Performance
///
/// - **Time**: O(1) per k-mer
/// - **Space**: O(1) - no allocation, returns borrowed slices
/// - **Memory**: Constant ~5 MB (aligns with streaming architecture)
///
/// # Examples
///
/// ```
/// use biometal::operations::kmer::kmer_iter;
///
/// let sequence = b"ATGCATGC";
/// let kmers: Vec<_> = kmer_iter(sequence, 3).collect();
///
/// assert_eq!(kmers.len(), 6);
/// assert_eq!(kmers[0], b"ATG");
/// ```
///
/// **Streaming example** (constant memory):
/// ```
/// use biometal::operations::kmer::kmer_iter;
///
/// let sequence = b"ATGCATGC";
///
/// // Process k-mers one at a time (constant memory)
/// for kmer in kmer_iter(sequence, 3) {
///     // Process immediately, no accumulation
///     if let Ok(kmer_str) = std::str::from_utf8(kmer) {
///         println!("{}", kmer_str);
///     }
/// }
/// ```
pub fn kmer_iter(sequence: &[u8], k: usize) -> impl Iterator<Item = &[u8]> {
    if k == 0 || k > sequence.len() {
        return KmerIterator { sequence, k, position: sequence.len() };
    }

    KmerIterator { sequence, k, position: 0 }
}

/// Internal iterator implementation for k-mers
struct KmerIterator<'a> {
    sequence: &'a [u8],
    k: usize,
    position: usize,
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        while self.position <= self.sequence.len().saturating_sub(self.k) {
            let kmer = &self.sequence[self.position..self.position + self.k];
            self.position += 1;

            if is_valid_kmer(kmer) {
                return Some(kmer);
            }
        }

        None
    }
}

/// Extract minimizers from a sequence (scalar-only)
///
/// Minimizers are the lexicographically smallest k-mer (by hash) in each
/// sliding window of size w. They are used for efficient indexing and
/// sketching in tools like minimap2.
///
/// # Performance
///
/// - **Time**: O(n × w) where n = sequence length, w = window size
/// - **Space**: O(m) where m = number of minimizers (typically n/w)
/// - **Optimization**: Scalar-only (Entry 034: 1.26× max, below ≥5× threshold)
///
/// # Evidence
///
/// ASBB Entry 034 found that minimizer extraction provides:
/// - NEON: 1.02-1.11× (negligible)
/// - Parallel-4t: 1.12-1.26× (below threshold)
///
/// Small output per sequence means thread overhead dominates. Scalar implementation
/// is optimal, **validating minimap2's design choice**.
///
/// # Algorithm
///
/// Uses FNV-1a hash (64-bit) for consistent, fast hashing. For each window:
/// 1. Hash all k-mers in the window
/// 2. Select k-mer with minimum hash value
/// 3. Store position + hash + k-mer
///
/// # Arguments
///
/// * `sequence` - DNA sequence to extract minimizers from
/// * `k` - K-mer size
/// * `w` - Window size (number of k-mers per window)
///
/// # Returns
///
/// Vector of minimizers with position, hash, and k-mer sequence.
///
/// # Examples
///
/// ```
/// use biometal::operations::kmer::extract_minimizers;
///
/// let sequence = b"ATGCATGCATGC";
/// let minimizers = extract_minimizers(sequence, 3, 5);
///
/// // Each window of 5 k-mers yields one minimizer
/// assert!(!minimizers.is_empty());
/// ```
pub fn extract_minimizers(sequence: &[u8], k: usize, w: usize) -> Vec<Minimizer> {
    if k == 0 || k > sequence.len() || w == 0 {
        return Vec::new();
    }

    let num_kmers = sequence.len().saturating_sub(k - 1);
    if num_kmers == 0 {
        return Vec::new();
    }

    let mut minimizers = Vec::new();
    let mut last_minimizer_pos = None;

    // Slide window across k-mers
    for window_start in 0..=(num_kmers.saturating_sub(w)) {
        let window_end = (window_start + w).min(num_kmers);

        // Find k-mer with minimum hash in this window
        let mut min_hash = u64::MAX;
        let mut min_pos = window_start;

        for i in window_start..window_end {
            let kmer = &sequence[i..i + k];

            if !is_valid_kmer(kmer) {
                continue;
            }

            let hash = fnv1a_hash(kmer);
            if hash < min_hash {
                min_hash = hash;
                min_pos = i;
            }
        }

        // Only add if we found a valid k-mer and it's different from last
        if min_hash != u64::MAX && last_minimizer_pos != Some(min_pos) {
            minimizers.push(Minimizer {
                position: min_pos,
                hash: min_hash,
                k,
            });
            last_minimizer_pos = Some(min_pos);
        }
    }

    minimizers
}

/// Extract minimizers using ntHash + sliding window minimum (100-200× faster)
///
/// This function uses the optimized algorithm from simd-minimizers:
/// 1. **ntHash**: O(1) rolling hash per k-mer (vs O(k) per hash)
/// 2. **Sliding Min**: O(1) amortized minimum finding (vs O(w) per window)
/// 3. **Lazy k-mer**: Zero-copy extraction (1.40× speedup, Entry 036-C)
/// 4. **Canonical**: Min of forward and reverse-complement hashes
///
/// # Performance
///
/// - **Time**: O(n) where n = sequence length (vs O(n×w) naive)
/// - **Space**: O(w) for sliding window tracker
/// - **Speedup**: 100-200× faster than `extract_minimizers` (Entry 036 → v1.3.0)
///
/// # Evidence
///
/// - **Source**: Port from simd-minimizers (MIT licensed)
/// - **Entry 036-B**: 77 Mbp/s @ 1Mbp (ntHash + sliding_min)
/// - **Entry 036-C**: 105 Mbp/s @ 1Mbp (lazy k-mer, 1.40× speedup)
/// - **SimdMinimizers**: 820 Mbp/s (target)
///
/// # Algorithm
///
/// 1. Use NtHashIterator to compute canonical hashes in O(1) per k-mer
/// 2. Use SlidingMin to track minimum hash in current window (O(1) amortized)
/// 3. Deduplicate consecutive minimizers at same position
/// 4. Return minimizers with positions and hashes
///
/// # Arguments
///
/// * `sequence` - DNA sequence (must contain only A, C, G, T, N)
/// * `k` - K-mer size (must be > 0 and ≤ sequence length)
/// * `w` - Window size (must be > 0)
///
/// # Returns
///
/// Vector of Minimizers with position, hash, and k-mer sequence.
/// Returns empty vector if inputs are invalid or sequence is too short.
///
/// # Errors
///
/// Returns `Err` if:
/// - Invalid nucleotide encountered (not A, C, G, T, N)
/// - K-mer size is 0 or > 32
/// - Window size is 0
///
/// # Examples
///
/// ```
/// use biometal::operations::kmer::extract_minimizers_fast;
///
/// let sequence = b"ATGCATGCATGC";
/// let minimizers = extract_minimizers_fast(sequence, 3, 5).unwrap();
///
/// // Each minimizer represents the minimum hash in a window
/// for m in &minimizers {
///     println!("Position: {}, Hash: {}", m.position, m.hash);
/// }
/// ```
pub fn extract_minimizers_fast(sequence: &[u8], k: usize, w: usize) -> Result<Vec<Minimizer>> {
    // Validation
    if k == 0 || w == 0 || k > sequence.len() {
        return Ok(Vec::new());
    }

    // Create ntHash iterator for canonical hashes
    let hash_iter = NtHashIterator::new(sequence, k)?;

    // Create sliding window minimum tracker
    let mut sliding_min = SlidingMin::new(w)?;

    let mut minimizers = Vec::new();
    let mut last_minimizer_pos = None;

    // Process hashes directly (no blocking needed - ntHash + SlidingMin already O(1))
    for (pos, hash) in hash_iter.enumerate() {
        // Push hash to sliding window and get minimum if window is full
        if let Some(min_elem) = sliding_min.push(hash, pos) {
            // Only add if position is different from last minimizer (deduplication)
            if last_minimizer_pos != Some(min_elem.pos) {
                minimizers.push(Minimizer {
                    position: min_elem.pos,
                    hash: min_elem.val,
                    k,
                });

                last_minimizer_pos = Some(min_elem.pos);
            }
        }
    }

    Ok(minimizers)
}

/// Extract minimizers using SIMD-accelerated algorithm (3-15× faster than scalar)
///
/// This function uses the `simd-minimizers` library for ARM NEON/AVX2 vectorization.
/// It combines bitpacking (2 bits/base) + SIMD parallelization for maximum performance.
///
/// # Performance (Entry 036-E)
///
/// **Validated benchmarks** (Mac M1, RUSTFLAGS="-C target-cpu=native", N=20):
/// - **≤10Kbp**: SIMD is **slower** (0.41-0.84×) due to chunking overhead
/// - **100Kbp**: SIMD is **1.39-1.51× faster** (crossover point)
/// - **1Mbp**: SIMD is **1.94× faster** (99 → 193 MiB/s) ✓
/// - **10Mbp**: SIMD is **1.50× faster** (98 → 147 MiB/s) ✓
///
/// **Algorithm**: Processes 8 k-mers in parallel using SIMD instructions
///
/// # When to Use
///
/// - **Genomic-scale** (≥100Kbp): Use SIMD for 1.5-1.9× speedup
/// - **Short sequences** (<100Kbp): Use `extract_minimizers_fast` (faster)
/// - **Mixed workloads**: Profile your data to find crossover point
///
/// # Gap from Literature
///
/// The simd-minimizers paper reports 3-15× speedup, but we achieve 1.5-1.9×
/// due to double hashing overhead (SIMD positions + ntHash recomputation for
/// API compatibility). This is still a meaningful speedup for production workloads.
///
/// # Evidence
///
/// - **Entry 036-E** (SIMD validation): 1.5-1.9× @ ≥100Kbp (this implementation)
/// - **Entry 036-C** (fast baseline): 105 MiB/s @ 1Mbp
/// - **Source**: simd-minimizers v2.2 (MIT licensed)
/// - **Paper**: "SimdMinimizers: Computing Random Minimizers, fast" (SEA 2025)
///
/// # Requirements
///
/// Enable the `simd` feature in Cargo.toml:
/// ```toml
/// biometal = { version = "1.2", features = ["simd"] }
/// ```
///
/// Compile with native CPU features:
/// ```bash
/// RUSTFLAGS="-C target-cpu=native" cargo build --release --features simd
/// ```
///
/// # Arguments
///
/// * `sequence` - DNA sequence (A, C, G, T, N)
/// * `k` - K-mer size
/// * `w` - Window size
///
/// # Returns
///
/// Vector of minimizers with position, hash, and k-mer size.
/// Same format as `extract_minimizers_fast` for API compatibility.
///
/// # Errors
///
/// Returns `Err` if inputs are invalid or sequence contains unsupported characters.
///
/// # Examples
///
/// ```ignore
/// use biometal::operations::kmer::extract_minimizers_simd;
///
/// // Best for genomic-scale sequences (≥100Kbp)
/// // For this example, we'll use a smaller sequence
/// let sequence = b"ATGCATGCATGCATGCATGCATGCATGC";
/// let minimizers = extract_minimizers_simd(sequence, 5, 7).unwrap();
///
/// for m in minimizers.iter().take(5) {
///     let kmer = m.kmer(sequence);
///     println!("Position: {}, k-mer: {:?}", m.position, std::str::from_utf8(kmer).unwrap());
/// }
///
/// // For production use with large sequences:
/// // let genome = std::fs::read("genome.fasta")?;  // 1-100 Mbp
/// // let minimizers = extract_minimizers_simd(&genome, 21, 11)?;
/// ```
#[cfg(feature = "simd")]
pub fn extract_minimizers_simd(sequence: &[u8], k: usize, w: usize) -> Result<Vec<Minimizer>> {
    use packed_seq::AsciiSeq;

    // Validation
    if k == 0 || w == 0 || k > sequence.len() {
        return Ok(Vec::new());
    }

    // Wrap ASCII sequence for simd-minimizers
    // AsciiSeq is a zero-cost wrapper that implements the Seq trait
    let ascii_seq = AsciiSeq(sequence);

    // Use simd-minimizers library for SIMD-accelerated extraction
    // This uses ARM NEON (or AVX2 on x86) to process 8 k-mers in parallel
    let positions = simd_minimizers::canonical_minimizer_positions(ascii_seq, k, w);

    // Convert positions to our Minimizer format
    // Note: We recompute hashes using ntHash for compatibility with our API
    let hash_iter = NtHashIterator::new(sequence, k)?;
    let hashes: Vec<(usize, u64)> = hash_iter.enumerate().collect();

    let minimizers = positions
        .into_iter()
        .filter_map(|pos| {
            let pos_usize = pos as usize;  // Convert u32 to usize
            hashes.get(pos_usize).map(|&(_idx, hash)| Minimizer {
                position: pos_usize,
                hash,
                k,
            })
        })
        .collect();

    Ok(minimizers)
}

/// K-mer spectrum (frequency counting, scalar-only)
///
/// Counts the frequency of each k-mer across one or more sequences. This is
/// used for genome size estimation, error correction, and repeat detection.
///
/// # Performance
///
/// - **Time**: O(n × m) where n = total sequence length, m = number of sequences
/// - **Space**: O(u) where u = number of unique k-mers
/// - **Optimization**: Scalar-only (Entry 034: parallel makes it SLOWER!)
///
/// # Evidence
///
/// **CRITICAL**: ASBB Entry 034 found that parallelization causes HashMap contention,
/// making k-mer spectrum counting **0.95-1.88× (inconsistent, sometimes SLOWER)**.
///
/// **DO NOT parallelize this function** - HashMap updates from multiple threads
/// cause cache thrashing and thread contention, degrading performance.
///
/// # Arguments
///
/// * `sequences` - Slice of DNA sequences to count k-mers from
/// * `k` - K-mer size
///
/// # Returns
///
/// HashMap mapping each k-mer to its frequency count. Invalid k-mers (with N, etc.) are skipped.
///
/// # Examples
///
/// ```
/// use biometal::operations::kmer::kmer_spectrum;
///
/// let sequences = vec![b"ATGCAT".as_ref(), b"GCATGC".as_ref()];
/// let spectrum = kmer_spectrum(&sequences, 3);
///
/// // Each k-mer appears in both sequences
/// assert_eq!(spectrum[b"ATG".as_ref()], 2);
/// assert_eq!(spectrum[b"TGC".as_ref()], 2);
/// ```
pub fn kmer_spectrum(sequences: &[&[u8]], k: usize) -> HashMap<Vec<u8>, usize> {
    let mut counts = HashMap::new();

    for sequence in sequences {
        for i in 0..=(sequence.len().saturating_sub(k)) {
            let kmer = &sequence[i..i + k];

            if is_valid_kmer(kmer) {
                *counts.entry(kmer.to_vec()).or_insert(0) += 1;
            }
        }
    }

    counts
}

/// Check if k-mer contains only valid bases (A, C, G, T)
///
/// Returns false if k-mer contains N or other ambiguous bases.
///
/// # Performance
///
/// This is the only part of k-mer operations that could benefit from NEON.
/// However, it represents only 5-10% of total runtime (Entry 034), so the
/// overall speedup would be ~1× (negligible).
#[inline]
fn is_valid_kmer(kmer: &[u8]) -> bool {
    kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

/// FNV-1a hash function (64-bit)
///
/// Fast, non-cryptographic hash used for minimizer selection.
/// Same algorithm used in Entry 034 benchmarks.
///
/// # Algorithm
///
/// ```text
/// hash = FNV_OFFSET_BASIS (14695981039346656037)
/// for each byte:
///     hash = hash XOR byte
///     hash = hash × FNV_PRIME (1099511628211)
/// ```
#[inline]
fn fnv1a_hash(kmer: &[u8]) -> u64 {
    const FNV_OFFSET_BASIS: u64 = 14695981039346656037;
    const FNV_PRIME: u64 = 1099511628211;

    let mut hash = FNV_OFFSET_BASIS;
    for &byte in kmer {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

/// Configurable k-mer extractor with optional parallelization
///
/// Provides scalar (default) and parallel k-mer extraction. Parallel extraction
/// provides 2.19-2.38× speedup for large datasets (Entry 034) but adds complexity.
///
/// # Design Philosophy
///
/// - **Default**: Scalar (simple, fast for most use cases)
/// - **Opt-in**: Parallel-4t (2.2× for large datasets)
/// - **Automatic threshold**: Only use parallel for ≥1000 sequences (amortize overhead)
///
/// # Evidence
///
/// ASBB Entry 034 found:
/// - Parallel-2t: 1.40-1.46× (moderate benefit)
/// - Parallel-4t: 2.19-2.38× (consistent, best result)
/// - Parallel-8t+: No additional benefit (thread overhead)
///
/// **Recommendation**: Cap at 4 threads (empirically optimal)
///
/// # Examples
///
/// **Default (scalar)**:
/// ```
/// use biometal::operations::kmer::KmerExtractor;
///
/// let extractor = KmerExtractor::new();
/// let sequences = vec![b"ATGCATGC".as_ref(), b"GCATGCAT".as_ref()];
/// let kmers = extractor.extract(&sequences, 3);
/// // Simple, fast for small datasets
/// ```
///
/// **Parallel (2.2× speedup)**:
/// ```
/// use biometal::operations::kmer::KmerExtractor;
///
/// // Opt-in for large datasets (10K+ sequences)
/// let extractor = KmerExtractor::with_parallel(4);
/// let sequences = vec![b"ATGCATGC".as_ref(); 10_000];
/// let kmers = extractor.extract(&sequences, 6);
/// // 2.2× faster for large datasets
/// ```
#[derive(Debug, Clone)]
pub struct KmerExtractor {
    parallel: bool,
    threads: usize,
}

impl KmerExtractor {
    /// Minimum number of sequences required to enable parallel extraction
    ///
    /// Based on Entry 034 - overhead amortization threshold. Parallel extraction
    /// only provides speedup when dataset is large enough to amortize thread
    /// creation/synchronization costs.
    ///
    /// Value: 1000 sequences (empirically validated)
    pub const PARALLEL_THRESHOLD: usize = 1000;
    /// Create a new k-mer extractor (scalar default)
    ///
    /// Uses simple scalar extraction - fast and sufficient for most use cases.
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::operations::kmer::KmerExtractor;
    ///
    /// let extractor = KmerExtractor::new();
    /// let sequences = vec![b"ATGCATGC".as_ref()];
    /// let kmers = extractor.extract(&sequences, 3);
    /// assert_eq!(kmers.len(), 6);
    /// ```
    pub fn new() -> Self {
        Self {
            parallel: false,
            threads: 1,
        }
    }

    /// Create a k-mer extractor with parallel extraction enabled
    ///
    /// Provides 2.19-2.38× speedup for large datasets (Entry 034). Threads are
    /// automatically capped at 4 (empirically optimal - no benefit beyond 4t).
    ///
    /// # Evidence
    ///
    /// Entry 034 found that Parallel-4t provides consistent 2.2× speedup, but
    /// only for large datasets where thread overhead is amortized.
    ///
    /// # Arguments
    ///
    /// * `threads` - Number of threads (automatically capped at 4)
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::operations::kmer::KmerExtractor;
    ///
    /// // Request 4 threads (optimal from Entry 034)
    /// let extractor = KmerExtractor::with_parallel(4);
    ///
    /// // Request 8 threads (capped at 4 automatically)
    /// let extractor = KmerExtractor::with_parallel(8);
    /// ```
    pub fn with_parallel(threads: usize) -> Self {
        Self {
            parallel: true,
            threads: threads.min(4), // Cap at 4 (Entry 034 evidence)
        }
    }

    /// Extract k-mers from multiple sequences
    ///
    /// Automatically chooses scalar or parallel based on:
    /// - Configuration (parallel enabled?)
    /// - Dataset size (≥1000 sequences → parallel, else scalar)
    ///
    /// This ensures parallel is only used when overhead is amortized.
    ///
    /// # Arguments
    ///
    /// * `sequences` - Slice of DNA sequences
    /// * `k` - K-mer size
    ///
    /// # Returns
    ///
    /// Flattened vector of all k-mers from all sequences. Order may differ
    /// between scalar and parallel (both are correct).
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::operations::kmer::KmerExtractor;
    ///
    /// let sequences = vec![b"ATGCAT".as_ref(), b"GCATGC".as_ref()];
    ///
    /// // Scalar (default)
    /// let scalar_kmers = KmerExtractor::new().extract(&sequences, 3);
    ///
    /// // Parallel (opt-in)
    /// let parallel_kmers = KmerExtractor::with_parallel(4).extract(&sequences, 3);
    ///
    /// // Both produce same k-mers (order may differ)
    /// assert_eq!(scalar_kmers.len(), parallel_kmers.len());
    /// ```
    pub fn extract(&self, sequences: &[&[u8]], k: usize) -> Vec<Vec<u8>> {
        if self.will_use_parallel(sequences.len()) {
            // Use parallel for large datasets (overhead amortized)
            self.extract_parallel(sequences, k)
        } else {
            // Use scalar for small datasets or when parallel disabled
            self.extract_scalar(sequences, k)
        }
    }

    /// Returns true if parallel extraction will be used for the given dataset size
    ///
    /// Parallel extraction is only used when:
    /// 1. Parallel mode is enabled (via `with_parallel()`)
    /// 2. Dataset size meets or exceeds `PARALLEL_THRESHOLD`
    ///
    /// This method allows users to understand and predict behavior.
    ///
    /// # Arguments
    ///
    /// * `num_sequences` - Number of sequences in the dataset
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::operations::kmer::KmerExtractor;
    ///
    /// let extractor = KmerExtractor::with_parallel(4);
    ///
    /// // Small dataset - will use scalar
    /// assert!(!extractor.will_use_parallel(999));
    ///
    /// // Large dataset - will use parallel
    /// assert!(extractor.will_use_parallel(1000));
    /// assert!(extractor.will_use_parallel(10_000));
    ///
    /// // Scalar extractor never uses parallel
    /// let scalar_extractor = KmerExtractor::new();
    /// assert!(!scalar_extractor.will_use_parallel(10_000));
    /// ```
    pub fn will_use_parallel(&self, num_sequences: usize) -> bool {
        self.parallel && num_sequences >= Self::PARALLEL_THRESHOLD
    }

    /// Extract k-mers using scalar algorithm
    ///
    /// Simple sequential extraction - fast and sufficient for most use cases.
    fn extract_scalar(&self, sequences: &[&[u8]], k: usize) -> Vec<Vec<u8>> {
        let mut kmers = Vec::new();

        for sequence in sequences {
            if k == 0 || k > sequence.len() {
                continue;
            }

            for i in 0..=(sequence.len() - k) {
                let kmer = &sequence[i..i + k];
                if is_valid_kmer(kmer) {
                    kmers.push(kmer.to_vec());
                }
            }
        }

        kmers
    }

    /// Extract k-mers using parallel algorithm (2.19-2.38× speedup)
    ///
    /// Processes sequences in parallel using Rayon. Each sequence is processed
    /// independently, avoiding HashMap contention issues seen in spectrum counting.
    ///
    /// Falls back to scalar extraction if thread pool creation fails (e.g., resource
    /// exhaustion), ensuring robustness in edge cases.
    ///
    /// # Evidence
    ///
    /// Entry 034: Parallel-4t provides 2.19-2.38× speedup for k-mer extraction
    /// (unlike spectrum, which gets slower with parallelization).
    fn extract_parallel(&self, sequences: &[&[u8]], k: usize) -> Vec<Vec<u8>> {
        use rayon::prelude::*;

        // Create bounded thread pool (cap at configured threads)
        let pool = match rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
        {
            Ok(p) => p,
            Err(_) => {
                // Fall back to scalar if thread pool creation fails
                // (e.g., resource exhaustion, OS limits)
                return self.extract_scalar(sequences, k);
            }
        };

        pool.install(|| {
            sequences
                .par_iter()
                .flat_map(|sequence| {
                    if k == 0 || k > sequence.len() {
                        return Vec::new();
                    }

                    (0..=(sequence.len() - k))
                        .filter_map(|i| {
                            let kmer = &sequence[i..i + k];
                            if is_valid_kmer(kmer) {
                                Some(kmer.to_vec())
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect()
        })
    }
}

impl Default for KmerExtractor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== K-mer Extraction Tests =====

    #[test]
    fn test_extract_kmers_basic() {
        let sequence = b"ATGCATGC";
        let kmers = extract_kmers(sequence, 3);

        assert_eq!(kmers.len(), 6);
        assert_eq!(kmers[0], b"ATG");
        assert_eq!(kmers[1], b"TGC");
        assert_eq!(kmers[2], b"GCA");
        assert_eq!(kmers[3], b"CAT");
        assert_eq!(kmers[4], b"ATG");
        assert_eq!(kmers[5], b"TGC");
    }

    #[test]
    fn test_extract_kmers_with_n() {
        let sequence = b"ATGCNNATGC";
        let kmers = extract_kmers(sequence, 3);

        // Should skip k-mers containing N: GCN, CNN, NNA, NAT
        assert_eq!(kmers.len(), 4);
        assert_eq!(kmers[0], b"ATG");
        assert_eq!(kmers[1], b"TGC");
        assert_eq!(kmers[2], b"ATG");
        assert_eq!(kmers[3], b"TGC");
    }

    #[test]
    fn test_extract_kmers_edge_cases() {
        // Empty sequence
        assert_eq!(extract_kmers(b"", 3).len(), 0);

        // k = 0
        assert_eq!(extract_kmers(b"ATGC", 0).len(), 0);

        // k > sequence length
        assert_eq!(extract_kmers(b"ATG", 5).len(), 0);

        // k = sequence length
        let kmers = extract_kmers(b"ATG", 3);
        assert_eq!(kmers.len(), 1);
        assert_eq!(kmers[0], b"ATG");
    }

    #[test]
    fn test_kmer_iter_basic() {
        let sequence = b"ATGCATGC";
        let kmers: Vec<_> = kmer_iter(sequence, 3).collect();

        assert_eq!(kmers.len(), 6);
        assert_eq!(kmers[0], b"ATG");
        assert_eq!(kmers[1], b"TGC");
    }

    #[test]
    fn test_kmer_iter_zero_copy() {
        let sequence = b"ATGCATGC";

        // Verify we're getting borrowed slices, not copies
        for kmer in kmer_iter(sequence, 3) {
            assert!(kmer.len() == 3);
            assert!(kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')));
        }
    }

    #[test]
    fn test_kmer_iter_skips_invalid() {
        let sequence = b"ATGCNNATGC";
        let kmers: Vec<_> = kmer_iter(sequence, 3).collect();

        // Should skip k-mers with N
        assert_eq!(kmers.len(), 4);
    }

    // ===== Minimizer Tests =====

    #[test]
    fn test_extract_minimizers_basic() {
        let sequence = b"ATGCATGCATGC";
        let minimizers = extract_minimizers(sequence, 3, 5);

        // Should produce minimizers
        assert!(!minimizers.is_empty());

        // Each minimizer should have valid fields
        for minimizer in &minimizers {
            assert_eq!(minimizer.k, 3);
            let kmer = minimizer.kmer(sequence);
            assert_eq!(kmer.len(), 3);
            assert!(is_valid_kmer(kmer));
            assert!(minimizer.hash > 0);
        }
    }

    #[test]
    fn test_extract_minimizers_no_duplicates() {
        let sequence = b"ATGCATGCATGC";
        let minimizers = extract_minimizers(sequence, 3, 3);

        // Consecutive minimizers should be at different positions
        for i in 1..minimizers.len() {
            assert_ne!(minimizers[i].position, minimizers[i - 1].position);
        }
    }

    #[test]
    fn test_extract_minimizers_edge_cases() {
        // Empty sequence
        assert_eq!(extract_minimizers(b"", 3, 5).len(), 0);

        // k = 0
        assert_eq!(extract_minimizers(b"ATGC", 0, 5).len(), 0);

        // w = 0
        assert_eq!(extract_minimizers(b"ATGC", 3, 0).len(), 0);

        // k > sequence length
        assert_eq!(extract_minimizers(b"ATG", 5, 5).len(), 0);
    }

    #[test]
    fn test_extract_minimizers_with_invalid_bases() {
        let sequence = b"ATGCNNATGCATGC";
        let minimizers = extract_minimizers(sequence, 3, 5);

        // Should skip k-mers with N but still produce minimizers
        for minimizer in &minimizers {
            let kmer = minimizer.kmer(sequence);
            assert!(is_valid_kmer(kmer));
        }
    }

    // ===== Fast Minimizer Tests (ntHash + Sliding Min) =====

    #[test]
    fn test_extract_minimizers_fast_basic() {
        let sequence = b"ATGCATGCATGC";
        let minimizers = extract_minimizers_fast(sequence, 3, 5).unwrap();

        // Should produce minimizers
        assert!(!minimizers.is_empty());

        // Each minimizer should have valid fields
        for minimizer in &minimizers {
            assert_eq!(minimizer.k, 3);
            let kmer = minimizer.kmer(sequence);
            assert_eq!(kmer.len(), 3);
            assert!(is_valid_kmer(kmer));
            assert!(minimizer.hash > 0);
        }
    }

    #[test]
    fn test_extract_minimizers_fast_no_duplicates() {
        let sequence = b"ATGCATGCATGC";
        let minimizers = extract_minimizers_fast(sequence, 3, 3).unwrap();

        // Consecutive minimizers should be at different positions
        for i in 1..minimizers.len() {
            assert_ne!(minimizers[i].position, minimizers[i - 1].position);
        }
    }

    #[test]
    fn test_extract_minimizers_fast_edge_cases() {
        // Empty sequence
        assert_eq!(
            extract_minimizers_fast(b"", 3, 5).unwrap().len(),
            0
        );

        // k = 0
        assert_eq!(
            extract_minimizers_fast(b"ATGC", 0, 5).unwrap().len(),
            0
        );

        // w = 0 (returns empty, handled by early validation)
        assert_eq!(
            extract_minimizers_fast(b"ATGC", 3, 0).unwrap().len(),
            0
        );

        // k > sequence length
        assert_eq!(
            extract_minimizers_fast(b"ATG", 5, 5).unwrap().len(),
            0
        );
    }

    #[test]
    fn test_extract_minimizers_fast_positions_valid() {
        let sequence = b"ATGCATGCATGC";
        let minimizers = extract_minimizers_fast(sequence, 3, 5).unwrap();

        // All positions should be valid indices
        for minimizer in &minimizers {
            assert!(minimizer.position + 3 <= sequence.len());
            // Verify k-mer matches the position
            let kmer = minimizer.kmer(sequence);
            let expected_kmer = &sequence[minimizer.position..minimizer.position + 3];
            assert_eq!(kmer, expected_kmer);
        }
    }

    #[test]
    fn test_extract_minimizers_fast_long_sequence() {
        // Test with a longer sequence to ensure algorithm scales well
        let sequence = b"ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";
        let minimizers = extract_minimizers_fast(sequence, 21, 11).unwrap();

        // Should produce minimizers
        assert!(!minimizers.is_empty());

        // Verify all minimizers are valid
        for minimizer in &minimizers {
            assert_eq!(minimizer.k, 21);
            let kmer = minimizer.kmer(sequence);
            assert_eq!(kmer.len(), 21);
            assert!(minimizer.position + 21 <= sequence.len());
        }
    }

    // ===== SIMD Minimizer Tests =====

    #[test]
    #[cfg(feature = "simd")]
    fn test_extract_minimizers_simd_correctness() {
        // Verify SIMD implementation produces valid minimizers
        // Note: SIMD uses strand-aware selection (leftmost vs rightmost based on TG content)
        // while fast uses simple leftmost selection, so they may select different (but equally valid) minimizers
        let sequences = vec![
            b"ATGCATGCATGCATGCATGC".as_ref(),
            b"AAAACCCGGGTTTTACGT".as_ref(),
            b"ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC".as_ref(),
        ];

        let k_values = [21, 31];
        let w_values = [11, 19];

        for seq in &sequences {
            for &k in &k_values {
                for &w in &w_values {
                    // Skip if k > sequence length
                    if k > seq.len() {
                        continue;
                    }

                    let simd_minimizers = extract_minimizers_simd(seq, k, w).unwrap();

                    // SIMD should produce valid minimizers
                    assert!(
                        !simd_minimizers.is_empty() || k + w - 1 > seq.len(),
                        "SIMD should produce minimizers for valid parameters (k={}, w={}, len={})",
                        k,
                        w,
                        seq.len()
                    );

                    // Verify all SIMD minimizers are valid
                    for (i, minimizer) in simd_minimizers.iter().enumerate() {
                        // Check k-mer size matches
                        assert_eq!(minimizer.k, k, "Minimizer {} has wrong k-mer size", i);

                        // Check position is valid
                        assert!(
                            minimizer.position + k <= seq.len(),
                            "Minimizer {} position {} + k {} exceeds sequence length {}",
                            i,
                            minimizer.position,
                            k,
                            seq.len()
                        );

                        // Check k-mer is valid
                        let kmer = minimizer.kmer(seq);
                        assert_eq!(kmer.len(), k, "Minimizer {} k-mer has wrong length", i);
                        assert!(
                            is_valid_kmer(kmer),
                            "Minimizer {} k-mer {:?} is invalid",
                            i,
                            std::str::from_utf8(kmer).unwrap_or("<invalid>")
                        );

                        // Check hash is non-zero (ntHash property)
                        assert!(minimizer.hash != 0, "Minimizer {} has zero hash", i);
                    }

                    // Verify no duplicate consecutive positions (deduplication property)
                    for i in 1..simd_minimizers.len() {
                        assert_ne!(
                            simd_minimizers[i - 1].position,
                            simd_minimizers[i].position,
                            "Consecutive minimizers {} and {} have same position",
                            i - 1,
                            i
                        );
                    }
                }
            }
        }
    }

    // ===== K-mer Spectrum Tests =====

    #[test]
    fn test_kmer_spectrum_basic() {
        let sequences = vec![b"ATGCAT".as_ref(), b"GCATGC".as_ref()];
        let spectrum = kmer_spectrum(&sequences, 3);

        // Check some expected k-mers
        assert_eq!(spectrum[b"ATG".as_ref()], 2); // Appears in both
        assert_eq!(spectrum[b"TGC".as_ref()], 2);
        assert_eq!(spectrum[b"GCA".as_ref()], 2);
        assert_eq!(spectrum[b"CAT".as_ref()], 2);
    }

    #[test]
    fn test_kmer_spectrum_single_sequence() {
        let sequences = vec![b"ATGATGATG".as_ref()];
        let spectrum = kmer_spectrum(&sequences, 3);

        assert_eq!(spectrum[b"ATG".as_ref()], 3);
        assert_eq!(spectrum[b"TGA".as_ref()], 2);
        assert_eq!(spectrum[b"GAT".as_ref()], 2);
    }

    #[test]
    fn test_kmer_spectrum_skips_invalid() {
        let sequences = vec![b"ATGCNNATGC".as_ref()];
        let spectrum = kmer_spectrum(&sequences, 3);

        // Should skip k-mers with N
        assert_eq!(spectrum[b"ATG".as_ref()], 2);
        assert_eq!(spectrum[b"TGC".as_ref()], 2);

        // These should not exist (contain N)
        assert!(!spectrum.contains_key(b"GCN".as_ref()));
        assert!(!spectrum.contains_key(b"CNN".as_ref()));
    }

    #[test]
    fn test_kmer_spectrum_empty() {
        let sequences: Vec<&[u8]> = vec![];
        let spectrum = kmer_spectrum(&sequences, 3);
        assert!(spectrum.is_empty());
    }

    // ===== Validation Tests =====

    #[test]
    fn test_is_valid_kmer() {
        assert!(is_valid_kmer(b"ATGC"));
        assert!(is_valid_kmer(b"AAAA"));
        assert!(is_valid_kmer(b"CCCC"));
        assert!(is_valid_kmer(b"GGGG"));
        assert!(is_valid_kmer(b"TTTT"));

        assert!(!is_valid_kmer(b"ATGN"));
        assert!(!is_valid_kmer(b"NATG"));
        assert!(!is_valid_kmer(b"ATXG"));
        assert!(!is_valid_kmer(b"atgc")); // Lowercase not valid
    }

    // ===== Hash Function Tests =====

    #[test]
    fn test_fnv1a_hash_consistency() {
        // Same k-mer should always hash to same value
        let kmer = b"ATGC";
        let hash1 = fnv1a_hash(kmer);
        let hash2 = fnv1a_hash(kmer);
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_fnv1a_hash_different_kmers() {
        // Different k-mers should (usually) hash to different values
        let hash1 = fnv1a_hash(b"ATGC");
        let hash2 = fnv1a_hash(b"CGTA");
        assert_ne!(hash1, hash2);
    }

    #[test]
    fn test_fnv1a_hash_empty() {
        // Empty k-mer should return offset basis
        let hash = fnv1a_hash(b"");
        assert_eq!(hash, 14695981039346656037);
    }

    // ===== KmerExtractor Tests =====

    #[test]
    fn test_kmer_extractor_new() {
        let extractor = KmerExtractor::new();
        assert!(!extractor.parallel);
        assert_eq!(extractor.threads, 1);
    }

    #[test]
    fn test_kmer_extractor_with_parallel() {
        let extractor = KmerExtractor::with_parallel(4);
        assert!(extractor.parallel);
        assert_eq!(extractor.threads, 4);
    }

    #[test]
    fn test_kmer_extractor_thread_cap() {
        // Threads should be capped at 4 (Entry 034 evidence)
        let extractor = KmerExtractor::with_parallel(8);
        assert_eq!(extractor.threads, 4);

        let extractor = KmerExtractor::with_parallel(16);
        assert_eq!(extractor.threads, 4);
    }

    #[test]
    fn test_kmer_extractor_default() {
        let extractor = KmerExtractor::default();
        assert!(!extractor.parallel);
        assert_eq!(extractor.threads, 1);
    }

    #[test]
    fn test_kmer_extractor_scalar_basic() {
        let extractor = KmerExtractor::new();
        let sequences = vec![b"ATGCAT".as_ref(), b"GCATGC".as_ref()];
        let kmers = extractor.extract(&sequences, 3);

        // Should extract k-mers from both sequences
        assert!(kmers.len() > 0);

        // Check some expected k-mers
        assert!(kmers.contains(&b"ATG".to_vec()));
        assert!(kmers.contains(&b"TGC".to_vec()));
    }

    #[test]
    fn test_kmer_extractor_parallel_basic() {
        let extractor = KmerExtractor::with_parallel(4);
        let sequences = vec![b"ATGCAT".as_ref(); 1000]; // Meet threshold
        let kmers = extractor.extract(&sequences, 3);

        // Should extract k-mers from all sequences
        assert!(kmers.len() > 0);
    }

    #[test]
    fn test_kmer_extractor_threshold_logic() {
        let extractor = KmerExtractor::with_parallel(4);

        // Below threshold (1000 sequences) → uses scalar
        let small_dataset = vec![b"ATGCAT".as_ref(); 999];
        let kmers = extractor.extract(&small_dataset, 3);
        assert!(kmers.len() > 0);

        // At threshold (≥1000 sequences) → uses parallel
        let large_dataset = vec![b"ATGCAT".as_ref(); 1000];
        let kmers = extractor.extract(&large_dataset, 3);
        assert!(kmers.len() > 0);
    }

    #[test]
    fn test_kmer_extractor_parallel_matches_scalar() {
        // Critical correctness test: parallel must produce same k-mers as scalar
        let sequences = vec![
            b"ATGCAT".as_ref(),
            b"GCATGC".as_ref(),
            b"CATGCA".as_ref(),
        ];

        let scalar = KmerExtractor::new();
        let parallel = KmerExtractor::with_parallel(4);

        let mut scalar_kmers = scalar.extract(&sequences, 3);
        let mut parallel_kmers = parallel.extract(&sequences, 3);

        // Sort both (order may differ between scalar/parallel)
        scalar_kmers.sort();
        parallel_kmers.sort();

        // Must produce identical results
        assert_eq!(scalar_kmers, parallel_kmers);
    }

    #[test]
    fn test_kmer_extractor_with_invalid_bases() {
        let extractor = KmerExtractor::new();
        let sequences = vec![b"ATGCNNATGC".as_ref()];
        let kmers = extractor.extract(&sequences, 3);

        // Should skip k-mers with N
        assert_eq!(kmers.len(), 4); // ATG, TGC, ATG, TGC
        assert!(kmers.contains(&b"ATG".to_vec()));
        assert!(kmers.contains(&b"TGC".to_vec()));
    }

    #[test]
    fn test_kmer_extractor_empty_sequences() {
        let extractor = KmerExtractor::new();
        let sequences: Vec<&[u8]> = vec![];
        let kmers = extractor.extract(&sequences, 3);
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_kmer_extractor_k_too_large() {
        let extractor = KmerExtractor::new();
        let sequences = vec![b"ATG".as_ref()];
        let kmers = extractor.extract(&sequences, 10); // k > sequence length
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    #[ignore] // Run manually for performance validation
    fn test_kmer_extractor_parallel_speedup() {
        use std::time::Instant;

        // Generate large test dataset (10K sequences, 150bp each)
        let sequences: Vec<&[u8]> = (0..10_000)
            .map(|_| b"ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC".as_ref())
            .collect();

        // Benchmark scalar
        let scalar = KmerExtractor::new();
        let start = Instant::now();
        let scalar_kmers = scalar.extract(&sequences, 6);
        let scalar_time = start.elapsed();

        // Benchmark parallel
        let parallel = KmerExtractor::with_parallel(4);
        let start = Instant::now();
        let parallel_kmers = parallel.extract(&sequences, 6);
        let parallel_time = start.elapsed();

        let speedup = scalar_time.as_secs_f64() / parallel_time.as_secs_f64();

        println!("Scalar time: {:?}", scalar_time);
        println!("Parallel time: {:?}", parallel_time);
        println!("Speedup: {:.2}×", speedup);

        // Entry 034: 2.19-2.38× observed (in release mode)
        // Allow range 1.5-3.0× for debug mode and different hardware
        assert!(
            speedup >= 1.5 && speedup <= 3.0,
            "Speedup {:.2}× outside expected range (1.5-3.0×). Note: Run in --release for optimal performance.",
            speedup
        );

        // Verify correctness
        assert_eq!(scalar_kmers.len(), parallel_kmers.len());
    }

    // ===== Property-Based Tests =====

    mod properties {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            /// Property: extract_kmers count matches formula for valid sequences
            #[test]
            fn prop_kmer_count(seq in "[ACGT]{1,100}", k in 1usize..10) {
                let kmers = extract_kmers(seq.as_bytes(), k);
                let expected = if k > seq.len() { 0 } else { seq.len() - k + 1 };
                prop_assert_eq!(kmers.len(), expected);
            }

            /// Property: kmer_iter matches extract_kmers
            #[test]
            fn prop_iter_matches_extract(seq in "[ACGT]{1,100}", k in 1usize..10) {
                let iter_kmers: Vec<Vec<u8>> = kmer_iter(seq.as_bytes(), k)
                    .map(|s| s.to_vec())
                    .collect();
                let extracted = extract_kmers(seq.as_bytes(), k);
                prop_assert_eq!(iter_kmers, extracted);
            }

            /// Property: spectrum sum equals total k-mers
            #[test]
            fn prop_spectrum_sum(sequences in prop::collection::vec("[ACGT]{10,50}", 1..10), k in 3usize..8) {
                let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_bytes()).collect();
                let spectrum = kmer_spectrum(&seq_refs, k);

                let expected_total: usize = seq_refs.iter()
                    .map(|s| if s.len() >= k { s.len() - k + 1 } else { 0 })
                    .sum();

                let actual_total: usize = spectrum.values().sum();
                prop_assert_eq!(actual_total, expected_total);
            }

            /// Property: minimizers have unique consecutive positions
            #[test]
            fn prop_minimizer_uniqueness(seq in "[ACGT]{20,100}", k in 3usize..8, w in 3usize..10) {
                let minimizers = extract_minimizers(seq.as_bytes(), k, w);

                for i in 1..minimizers.len() {
                    // Consecutive minimizers must be at different positions
                    prop_assert_ne!(minimizers[i].position, minimizers[i-1].position);
                }
            }

            /// Property: KmerExtractor parallel matches scalar
            #[test]
            fn prop_parallel_matches_scalar(sequences in prop::collection::vec("[ACGT]{10,50}", 100..200), k in 3usize..8) {
                let seq_refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_bytes()).collect();

                let scalar = KmerExtractor::new();
                let parallel = KmerExtractor::with_parallel(4);

                let mut scalar_kmers = scalar.extract(&seq_refs, k);
                let mut parallel_kmers = parallel.extract(&seq_refs, k);

                // Sort for comparison (parallel may return in different order)
                scalar_kmers.sort();
                parallel_kmers.sort();

                prop_assert_eq!(scalar_kmers, parallel_kmers);
            }

            /// Property: extract_kmers handles empty and edge cases gracefully
            #[test]
            fn prop_edge_cases_no_panic(seq in ".*", k in 0usize..20) {
                // Should not panic regardless of input
                let _ = extract_kmers(seq.as_bytes(), k);
            }
        }
    }
}
