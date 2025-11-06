//! Sequence complexity calculation with ARM NEON SIMD optimization
//!
//! # Evidence
//!
//! Entry 020-025 (Lab Notebook) - Base counting pattern:
//! - **Expected speedup**: 15-18× (NEON-optimized counting dominates)
//! - **Algorithm**: Shannon entropy H = -Σ(p_i × log₂(p_i))
//! - **Category 1**: Similar to base_counting (16.7×), gc_content (20.3×)
//!
//! # Algorithm
//!
//! Shannon entropy measures sequence randomness:
//! - **0.0**: Homopolymer (AAAA...) - minimal complexity
//! - **0.5**: Dinucleotide repeat (ATATATAT...) - low complexity
//! - **1.0**: Random sequence - maximal complexity
//!
//! # Architecture
//!
//! 1. **Count bases** (NEON-optimized, dominates runtime)
//! 2. **Calculate entropy** (scalar, negligible cost)
//! 3. **Normalize** to 0.0-1.0 range
//!
//! # Use Cases
//!
//! - Filter low-complexity regions (repeats, homopolymers)
//! - QC for metagenomics (avoid biased sequences)
//! - Validate random sequence generation
//!
//! # Examples
//!
//! ```
//! use biometal::operations::complexity_score;
//!
//! // High complexity (random)
//! let random = b"ATGCTAGCTAGCTAGC";
//! let score = complexity_score(random);
//! assert!(score > 0.9);
//!
//! // Low complexity (homopolymer)
//! let homopolymer = b"AAAAAAAAAAAAAAAA";
//! let score = complexity_score(homopolymer);
//! assert!(score < 0.1);
//!
//! // Medium complexity (dinucleotide repeat)
//! let repeat = b"ATATATATATATATAT";
//! let score = complexity_score(repeat);
//! assert!(score > 0.3 && score < 0.7);
//! ```

use crate::operations::base_counting::{count_bases, BaseCounts};

/// Calculate sequence complexity score (0.0-1.0)
///
/// Uses Shannon entropy to measure sequence randomness. Higher scores indicate
/// more complex (random) sequences, while lower scores indicate repetitive or
/// low-complexity regions.
///
/// # Algorithm
///
/// 1. Count base frequencies (NEON-optimized)
/// 2. Calculate Shannon entropy: H = -Σ(p_i × log₂(p_i))
/// 3. Normalize to 0.0-1.0 range (H / log₂(4))
///
/// # Platform-Specific Optimization
///
/// - **ARM (aarch64)**: NEON SIMD for base counting (15-18× speedup)
/// - **x86_64**: Scalar fallback (portable)
///
/// # Returns
///
/// - **0.0**: Homopolymer (minimal complexity)
/// - **0.5**: Low complexity (e.g., dinucleotide repeats)
/// - **1.0**: Random sequence (maximal complexity)
///
/// # Evidence
///
/// Entry 020-025: Base counting achieves 16.7× NEON speedup (Cohen's d = 4.82)
/// - Complexity calculation reuses this pattern (counting dominates runtime)
/// - Expected speedup: 15-18× (entropy calculation is negligible)
///
/// # Examples
///
/// ```
/// use biometal::operations::complexity_score;
///
/// // High complexity (random sequence)
/// let random = b"ATGCTAGCTAGCTAGC";
/// let score = complexity_score(random);
/// assert!(score > 0.9);
///
/// // Low complexity (homopolymer)
/// let homopolymer = b"AAAAAAAAAAAAAAAA";
/// let score = complexity_score(homopolymer);
/// assert!(score < 0.1);
///
/// // Filter low-complexity regions (metagenomics example from README)
/// let sequence = b"ATGCATGCATGCATGC";
/// if complexity_score(sequence) > 0.5 {
///     // High enough complexity - process this sequence
///     println!("Sequence passed complexity filter");
/// }
/// ```
pub fn complexity_score(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }

    // Count bases using NEON-optimized implementation (16.7× speedup)
    let counts = count_bases(seq);

    // Calculate Shannon entropy from base counts
    calculate_entropy_from_counts(&counts, seq.len())
}

/// Calculate Shannon entropy from base counts
///
/// H = -Σ(p_i × log₂(p_i)), normalized to 0.0-1.0
///
/// # Arguments
///
/// * `counts` - Base counts [A, C, G, T]
/// * `total` - Total sequence length
///
/// # Returns
///
/// Normalized Shannon entropy (0.0-1.0)
///
/// # Performance
///
/// This function is negligible cost compared to base counting:
/// - Base counting: ~300-500 CPU cycles (NEON) or ~5000 cycles (scalar)
/// - Entropy calculation: ~50 CPU cycles (4 additions, 4 divisions, 4 logs)
/// - **Speedup is dominated by NEON base counting** (15-18× expected)
fn calculate_entropy_from_counts(counts: &BaseCounts, total: usize) -> f64 {
    if total == 0 {
        return 0.0;
    }

    let total_f64 = total as f64;
    let mut entropy = 0.0;

    for &count in counts {
        if count > 0 {
            let p = count as f64 / total_f64;
            // Shannon entropy: -Σ(p_i × log₂(p_i))
            entropy -= p * p.log2();
        }
    }

    // Normalize to 0.0-1.0 range
    // Maximum entropy for 4-base alphabet is log₂(4) = 2.0
    let max_entropy = 2.0;
    entropy / max_entropy
}

/// Calculate complexity score using scalar base counting
///
/// Exposed for testing and benchmarking. Users should call `complexity_score()`
/// which automatically selects the best implementation.
pub fn complexity_score_scalar(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }

    // Use scalar base counting
    use crate::operations::base_counting::count_bases_scalar;
    let counts = count_bases_scalar(seq);

    calculate_entropy_from_counts(&counts, seq.len())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complexity_homopolymer() {
        // All A's - minimal complexity
        let seq = b"AAAAAAAAAAAAAAAA";
        let score = complexity_score(seq);
        assert!(score < 0.1, "Homopolymer should have low complexity: {}", score);
        assert!(score >= 0.0, "Score should be non-negative");
    }

    #[test]
    fn test_complexity_random() {
        // Random-looking sequence - high complexity
        let seq = b"ATGCTAGCTAGCGATC";
        let score = complexity_score(seq);
        assert!(score > 0.8, "Random sequence should have high complexity: {}", score);
        assert!(score <= 1.0, "Score should not exceed 1.0");
    }

    #[test]
    fn test_complexity_dinucleotide_repeat() {
        // ATATATAT... - medium complexity
        let seq = b"ATATATATATATATAT";
        let score = complexity_score(seq);
        assert!(score > 0.3 && score < 0.7, "Dinucleotide repeat should have medium complexity: {}", score);
    }

    #[test]
    fn test_complexity_equal_bases() {
        // Equal distribution of all bases - maximal complexity
        let seq = b"ACGTACGTACGTACGT";
        let score = complexity_score(seq);
        assert!(score > 0.95, "Equal base distribution should have maximal complexity: {}", score);
    }

    #[test]
    fn test_complexity_empty() {
        let seq = b"";
        let score = complexity_score(seq);
        assert_eq!(score, 0.0, "Empty sequence should have zero complexity");
    }

    #[test]
    fn test_complexity_single_base() {
        let seq = b"A";
        let score = complexity_score(seq);
        assert_eq!(score, 0.0, "Single base should have zero complexity");
    }

    #[test]
    fn test_complexity_two_bases() {
        let seq = b"AT";
        let score = complexity_score(seq);
        assert!(score > 0.4 && score < 0.6, "Two different bases should have medium complexity");
    }

    #[test]
    fn test_complexity_with_n() {
        // N's are ignored by base_counting, so they don't contribute to entropy
        let seq = b"ACGTNNNACGT";
        let score = complexity_score(seq);
        // Only ACGT bases count (2 of each = perfect distribution), should have high complexity
        // Score is slightly below 1.0 due to N's not being counted
        assert!(score > 0.85, "High complexity bases with N's should still score high: {}", score);
    }

    #[test]
    fn test_complexity_long_homopolymer() {
        // Long homopolymer (tests NEON chunking)
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 40 A's
        let score = complexity_score(seq);
        assert!(score < 0.1, "Long homopolymer should have low complexity");
    }

    #[test]
    fn test_complexity_long_random() {
        // Long random-looking sequence (tests NEON chunking)
        let seq = b"ATGCTAGCTAGCGATCTAGCGATCGATCGATCGATCGATC"; // 40 bases
        let score = complexity_score(seq);
        assert!(score > 0.8, "Long random sequence should have high complexity");
    }

    #[test]
    fn test_metagenomics_filter_example() {
        // Example from README.md line 513
        let low_complexity = b"AAAAAAAAAAAAAAAA";
        let high_complexity = b"ATGCTAGCTAGCGATC";

        // Low complexity should be filtered out
        assert!(complexity_score(low_complexity) <= 0.5);

        // High complexity should pass
        assert!(complexity_score(high_complexity) > 0.5);
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn test_complexity_score_uses_neon_optimized_counting() {
        // Verify that complexity_score() produces correct results
        // (NEON optimization is tested in base_counting module)

        // Homopolymer: 0 bits entropy
        assert!((complexity_score(b"AAAA") - 0.0).abs() < 0.01);

        // Two bases: 1 bit entropy, normalized to 0.5
        assert!((complexity_score(b"ATATATAT") - 0.5).abs() < 0.01);

        // Four bases with equal distribution: ~2 bits entropy, normalized to 1.0
        assert!(complexity_score(b"ATGCTAGCTAGCGATC") > 0.95);
        assert!(complexity_score(b"ACGTACGTACGTACGT") > 0.95);
        assert!(complexity_score(b"ACGTACGTACGTACGTACGT") > 0.95);
    }
}
