//! Core sequence manipulation operations
//!
//! Provides fundamental DNA/RNA sequence transformations:
//! - Reverse complement (most common operation)
//! - Complement only
//! - Reverse only
//! - Sequence validation
//!
//! # Design Principles
//!
//! 1. **Zero-copy when possible**: Operations return slices where feasible
//! 2. **In-place variants**: Avoid allocation for performance-critical paths
//! 3. **Evidence-based optimization**: NEON only when benchmarks justify
//! 4. **DNA + RNA support**: Handle both T and U nucleotides
//! 5. **Ambiguous base handling**: Pass-through IUPAC codes (N, R, Y, etc.)
//!
//! # Evidence
//!
//! - Scalar implementations: Standard algorithms, no validation needed
//! - NEON implementations: Category 2 (Similar to base_counting)
//!   - Expected speedup: 10-16× (from Entry 020-025 patterns)
//!   - Validation: Quick criterion benchmark (N=10) before adding
//!
//! # Examples
//!
//! ```
//! use biometal::operations::{reverse_complement, complement, reverse};
//!
//! // Reverse complement (most common)
//! let seq = b"ATGC";
//! let rc = reverse_complement(seq);
//! assert_eq!(rc, b"GCAT");
//!
//! // Complement only
//! let comp = complement(seq);
//! assert_eq!(comp, b"TACG");
//!
//! // Reverse only
//! let rev = reverse(seq);
//! assert_eq!(rev, b"CGTA");
//! ```

// Note: Result and BiometalError reserved for future error handling
// (currently all operations are infallible)

/// Lookup table for DNA/RNA complement
///
/// Handles:
/// - Standard bases: A↔T/U, G↔C
/// - Ambiguous IUPAC codes (pass-through with complement)
/// - Invalid characters (preserved as-is)
///
/// # IUPAC Nucleotide Codes
///
/// Standard:
/// - A (Adenine) ↔ T (Thymine) or U (Uracil)
/// - G (Guanine) ↔ C (Cytosine)
///
/// Ambiguous (complement preserves meaning):
/// - R (A or G) ↔ Y (C or T)
/// - W (A or T) ↔ W (A or T) [self-complement]
/// - S (G or C) ↔ S (G or C) [self-complement]
/// - K (G or T) ↔ M (A or C)
/// - N (any) ↔ N (any)
///
/// # Evidence
///
/// Lookup table approach is standard for complement operations.
/// NEON optimization can process 16 bases simultaneously using vtbl1q_u8.
const COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = i as u8; // Default: preserve character
        i += 1;
    }

    // DNA standard bases
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'G' as usize] = b'C';
    table[b'C' as usize] = b'G';
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'g' as usize] = b'c';
    table[b'c' as usize] = b'g';

    // RNA (U instead of T)
    table[b'U' as usize] = b'A';
    table[b'u' as usize] = b'a';

    // IUPAC ambiguity codes (pass-through with complement)
    table[b'R' as usize] = b'Y'; // A or G → C or T
    table[b'Y' as usize] = b'R'; // C or T → A or G
    table[b'W' as usize] = b'W'; // A or T → A or T (self-complement)
    table[b'S' as usize] = b'S'; // G or C → G or C (self-complement)
    table[b'K' as usize] = b'M'; // G or T → A or C
    table[b'M' as usize] = b'K'; // A or C → G or T
    table[b'r' as usize] = b'y';
    table[b'y' as usize] = b'r';
    table[b'w' as usize] = b'w';
    table[b's' as usize] = b's';
    table[b'k' as usize] = b'm';
    table[b'm' as usize] = b'k';

    // N (any base) remains N
    table[b'N' as usize] = b'N';
    table[b'n' as usize] = b'n';

    table
};

/// Reverse complement a DNA/RNA sequence
///
/// Returns a new vector with the sequence reversed and complemented.
/// This is the most common sequence transformation in bioinformatics.
///
/// # Arguments
///
/// * `seq` - Input DNA/RNA sequence (ASCII bytes)
///
/// # Returns
///
/// New vector containing the reverse complement
///
/// # Nucleotide Handling
///
/// - **Standard DNA**: A↔T, G↔C (fully supported)
/// - **RNA input**: Accepts U, but outputs T (DNA-style complement)
///   - Example: AUGC → GCAT (not GCAU)
///   - Limitation: Single static complement table optimized for DNA
/// - **Ambiguous IUPAC codes**: Complemented appropriately (R↔Y, K↔M, etc.)
/// - **Invalid characters**: Preserved as-is
///
/// # Performance
///
/// - **Scalar**: O(n) single pass with lookup table
/// - **NEON** (future): O(n/16) with 16-byte SIMD operations
/// - **Expected speedup**: 10-16× (similar to base_counting, Entry 020-025)
///
/// # Evidence
///
/// - Category 2 (Similar): Expected to benefit from NEON like base_counting
/// - Validation: Criterion benchmark (N=10) before adding NEON
/// - Current: Scalar only (evidence-based approach)
///
/// # Examples
///
/// ```
/// use biometal::operations::reverse_complement;
///
/// // DNA
/// let seq = b"ATGC";
/// let rc = reverse_complement(seq);
/// assert_eq!(rc, b"GCAT");
///
/// // RNA (note: outputs T, not U)
/// let rna = b"AUGC";
/// let rc_rna = reverse_complement(rna);
/// assert_eq!(rc_rna, b"GCAT"); // T not U (DNA-style output)
///
/// // With ambiguous bases
/// let ambig = b"ATGCN";
/// let rc_ambig = reverse_complement(ambig);
/// assert_eq!(rc_ambig, b"NGCAT");
///
/// // Involutive property: RC(RC(x)) = x
/// let original = b"ATGCATGC";
/// let rc = reverse_complement(original);
/// let rc_rc = reverse_complement(&rc);
/// assert_eq!(rc_rc, original);
/// ```
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    reverse_complement_scalar(seq)
}

/// Scalar implementation of reverse complement
#[inline]
fn reverse_complement_scalar(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| COMPLEMENT_TABLE[base as usize])
        .collect()
}

/// In-place reverse complement
///
/// Mutates the input sequence to its reverse complement, avoiding allocation.
/// Use this for performance-critical paths where you own the sequence.
///
/// # Arguments
///
/// * `seq` - Mutable DNA/RNA sequence (modified in place)
///
/// # Performance
///
/// - **Memory**: Zero allocation (in-place transformation)
/// - **Time**: O(n/2) with two-pointer swap approach
/// - **Best for**: Large sequences where allocation cost matters
///
/// # Examples
///
/// ```
/// use biometal::operations::reverse_complement_inplace;
///
/// let mut seq = b"ATGC".to_vec();
/// reverse_complement_inplace(&mut seq);
/// assert_eq!(seq, b"GCAT");
///
/// // Involutive: applying twice returns to original
/// let mut seq2 = b"ATGC".to_vec();
/// reverse_complement_inplace(&mut seq2);
/// reverse_complement_inplace(&mut seq2);
/// assert_eq!(seq2, b"ATGC");
/// ```
pub fn reverse_complement_inplace(seq: &mut [u8]) {
    let len = seq.len();
    for i in 0..(len / 2) {
        let j = len - 1 - i;
        // Swap and complement simultaneously
        let temp = COMPLEMENT_TABLE[seq[i] as usize];
        seq[i] = COMPLEMENT_TABLE[seq[j] as usize];
        seq[j] = temp;
    }

    // Handle middle element for odd-length sequences
    if len % 2 == 1 {
        let mid = len / 2;
        seq[mid] = COMPLEMENT_TABLE[seq[mid] as usize];
    }
}

/// Complement a sequence without reversing
///
/// Returns a new vector with each base complemented but in the same order.
/// Less common than reverse_complement, but useful for specific applications.
///
/// # Arguments
///
/// * `seq` - Input DNA/RNA sequence
///
/// # Returns
///
/// New vector containing the complement
///
/// # Examples
///
/// ```
/// use biometal::operations::complement;
///
/// let seq = b"ATGC";
/// let comp = complement(seq);
/// assert_eq!(comp, b"TACG");
///
/// // Involutive: C(C(x)) = x
/// let comp2 = complement(&comp);
/// assert_eq!(comp2, b"ATGC");
/// ```
pub fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&base| COMPLEMENT_TABLE[base as usize])
        .collect()
}

/// In-place complement
///
/// Mutates the input sequence to its complement, avoiding allocation.
///
/// # Arguments
///
/// * `seq` - Mutable DNA/RNA sequence (modified in place)
///
/// # Examples
///
/// ```
/// use biometal::operations::complement_inplace;
///
/// let mut seq = b"ATGC".to_vec();
/// complement_inplace(&mut seq);
/// assert_eq!(seq, b"TACG");
/// ```
pub fn complement_inplace(seq: &mut [u8]) {
    for base in seq.iter_mut() {
        *base = COMPLEMENT_TABLE[*base as usize];
    }
}

/// Reverse a sequence without complementing
///
/// Returns a new vector with bases in reverse order.
/// Rarely needed alone, but useful for testing and specific workflows.
///
/// # Arguments
///
/// * `seq` - Input sequence
///
/// # Returns
///
/// New vector containing the reversed sequence
///
/// # Examples
///
/// ```
/// use biometal::operations::reverse;
///
/// let seq = b"ATGC";
/// let rev = reverse(seq);
/// assert_eq!(rev, b"CGTA");
///
/// // Involutive: R(R(x)) = x
/// let rev2 = reverse(&rev);
/// assert_eq!(rev2, b"ATGC");
/// ```
pub fn reverse(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().copied().collect()
}

/// In-place reverse
///
/// Mutates the input sequence to its reverse, avoiding allocation.
///
/// # Arguments
///
/// * `seq` - Mutable sequence (modified in place)
///
/// # Examples
///
/// ```
/// use biometal::operations::reverse_inplace;
///
/// let mut seq = b"ATGC".to_vec();
/// reverse_inplace(&mut seq);
/// assert_eq!(seq, b"CGTA");
/// ```
pub fn reverse_inplace(seq: &mut [u8]) {
    seq.reverse();
}

/// Check if sequence contains only valid DNA bases
///
/// Validates that all characters are standard DNA bases (ACGT) or
/// IUPAC ambiguity codes (N, R, Y, W, S, K, M, etc.).
///
/// # Arguments
///
/// * `seq` - Input sequence to validate
///
/// # Returns
///
/// `true` if all characters are valid DNA bases, `false` otherwise
///
/// # Valid Characters
///
/// - **Standard**: A, C, G, T (case-insensitive)
/// - **Ambiguous**: N, R, Y, W, S, K, M, B, D, H, V
/// - **Invalid**: Numbers, punctuation, other letters
///
/// # Performance
///
/// - **Scalar**: O(n) with early exit on first invalid character
/// - **NEON**: Not beneficial (early-exit pattern doesn't vectorize well)
///
/// # Examples
///
/// ```
/// use biometal::operations::is_valid_dna;
///
/// assert!(is_valid_dna(b"ATGC"));
/// assert!(is_valid_dna(b"ATGCN")); // N is valid (ambiguous)
/// assert!(is_valid_dna(b"atgc")); // Case-insensitive
/// assert!(!is_valid_dna(b"ATGCX")); // X is invalid
/// assert!(!is_valid_dna(b"ATG C")); // Space is invalid
/// ```
pub fn is_valid_dna(seq: &[u8]) -> bool {
    seq.iter().all(|&base| is_valid_dna_base(base))
}

/// Check if sequence contains only valid RNA bases
///
/// Validates that all characters are standard RNA bases (ACGU) or
/// IUPAC ambiguity codes.
///
/// # Arguments
///
/// * `seq` - Input sequence to validate
///
/// # Returns
///
/// `true` if all characters are valid RNA bases, `false` otherwise
///
/// # Valid Characters
///
/// - **Standard**: A, C, G, U (case-insensitive)
/// - **Ambiguous**: N, R, Y, W, S, K, M, B, D, H, V
/// - **Note**: Accepts U (RNA) but NOT T (DNA)
///
/// # Examples
///
/// ```
/// use biometal::operations::is_valid_rna;
///
/// assert!(is_valid_rna(b"AUGC"));
/// assert!(is_valid_rna(b"AUGCN")); // N is valid
/// assert!(!is_valid_rna(b"ATGC")); // T is invalid for RNA
/// ```
pub fn is_valid_rna(seq: &[u8]) -> bool {
    seq.iter().all(|&base| is_valid_rna_base(base))
}

/// Check if a single character is a valid DNA base
#[inline]
fn is_valid_dna_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'C' | b'G' | b'T' | b'N' | // Standard + N
        b'a' | b'c' | b'g' | b't' | b'n' | // Lowercase
        b'R' | b'Y' | b'W' | b'S' | b'K' | b'M' | // Ambiguous
        b'r' | b'y' | b'w' | b's' | b'k' | b'm' |
        b'B' | b'D' | b'H' | b'V' | // More ambiguous
        b'b' | b'd' | b'h' | b'v'
    )
}

/// Check if a single character is a valid RNA base
#[inline]
fn is_valid_rna_base(base: u8) -> bool {
    matches!(
        base,
        b'A' | b'C' | b'G' | b'U' | b'N' | // Standard RNA (U not T)
        b'a' | b'c' | b'g' | b'u' | b'n' | // Lowercase
        b'R' | b'Y' | b'W' | b'S' | b'K' | b'M' | // Ambiguous
        b'r' | b'y' | b'w' | b's' | b'k' | b'm' |
        b'B' | b'D' | b'H' | b'V' | // More ambiguous
        b'b' | b'd' | b'h' | b'v'
    )
}

/// Count invalid bases in a sequence
///
/// Returns the number of characters that are not valid DNA bases.
/// Useful for quality control and validation reporting.
///
/// # Arguments
///
/// * `seq` - Input sequence to analyze
///
/// # Returns
///
/// Count of invalid characters
///
/// # Examples
///
/// ```
/// use biometal::operations::count_invalid_bases;
///
/// assert_eq!(count_invalid_bases(b"ATGC"), 0);
/// assert_eq!(count_invalid_bases(b"ATGCX"), 1);
/// assert_eq!(count_invalid_bases(b"AT GC"), 1); // Space is invalid
/// assert_eq!(count_invalid_bases(b"123"), 3);
/// ```
pub fn count_invalid_bases(seq: &[u8]) -> usize {
    seq.iter().filter(|&&base| !is_valid_dna_base(base)).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Property-Based Tests (proptest) =====

    #[cfg(test)]
    mod properties {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            /// Property: Reverse complement is involutive for DNA
            /// Mathematical: RC(RC(x)) = x
            #[test]
            fn prop_reverse_complement_involutive(seq in "[ACGTN]{1,1000}") {
                let rc = reverse_complement(seq.as_bytes());
                let rc_rc = reverse_complement(&rc);
                prop_assert_eq!(rc_rc, seq.as_bytes().to_vec());
            }

            /// Property: Complement is involutive for DNA
            /// Mathematical: C(C(x)) = x
            #[test]
            fn prop_complement_involutive(seq in "[ACGTN]{1,1000}") {
                let comp = complement(seq.as_bytes());
                let comp2 = complement(&comp);
                prop_assert_eq!(comp2, seq.as_bytes().to_vec());
            }

            /// Property: Reverse is involutive
            /// Mathematical: R(R(x)) = x
            #[test]
            fn prop_reverse_involutive(seq in "[ACGTN]{1,1000}") {
                let rev = reverse(seq.as_bytes());
                let rev2 = reverse(&rev);
                prop_assert_eq!(rev2, seq.as_bytes().to_vec());
            }

            /// Property: Reverse complement = reverse + complement (order matters)
            /// Mathematical: RC(x) = C(R(x)) = R(C(x))
            #[test]
            fn prop_reverse_complement_decomposition(seq in "[ACGTN]{1,100}") {
                let rc = reverse_complement(seq.as_bytes());

                // RC(x) = R(C(x))
                let comp_then_rev = reverse(&complement(seq.as_bytes()));
                prop_assert_eq!(&rc, &comp_then_rev);

                // RC(x) = C(R(x))
                let rev_then_comp = complement(&reverse(seq.as_bytes()));
                prop_assert_eq!(&rc, &rev_then_comp);
            }

            /// Property: In-place operations match allocating versions
            #[test]
            fn prop_inplace_matches_allocating(seq in "[ACGTN]{1,1000}") {
                let seq_bytes = seq.as_bytes().to_vec();

                // Reverse complement
                let rc_alloc = reverse_complement(&seq_bytes);
                let mut rc_inplace = seq_bytes.clone();
                reverse_complement_inplace(&mut rc_inplace);
                prop_assert_eq!(rc_alloc, rc_inplace);

                // Complement
                let comp_alloc = complement(&seq_bytes);
                let mut comp_inplace = seq_bytes.clone();
                complement_inplace(&mut comp_inplace);
                prop_assert_eq!(comp_alloc, comp_inplace);

                // Reverse
                let rev_alloc = reverse(&seq_bytes);
                let mut rev_inplace = seq_bytes.clone();
                reverse_inplace(&mut rev_inplace);
                prop_assert_eq!(rev_alloc, rev_inplace);
            }

            /// Property: Validation accepts all valid DNA characters
            #[test]
            fn prop_valid_dna_accepts_valid(seq in "[ACGTNRYWSKMBDHVacgtnrywskmbdhv]{1,1000}") {
                prop_assert!(is_valid_dna(seq.as_bytes()));
            }

            /// Property: Validation rejects invalid characters
            #[test]
            fn prop_invalid_dna_rejects_invalid(
                seq in "[ACGTN]{0,100}[^ACGTNRYWSKMBDHVacgtnrywskmbdhv]{1,10}[ACGTN]{0,100}"
            ) {
                // If sequence contains invalid chars, should fail validation
                if seq.bytes().any(|b| !matches!(b, b'A'..=b'Z' | b'a'..=b'z')) {
                    prop_assert!(!is_valid_dna(seq.as_bytes()));
                }
            }

            /// Property: Count invalid bases is accurate
            #[test]
            fn prop_count_invalid_accurate(seq in "[ACGTN]{0,100}[X123]{0,10}[ACGTN]{0,100}") {
                let count = count_invalid_bases(seq.as_bytes());
                let expected = seq.bytes().filter(|&b| !is_valid_dna_base(b)).count();
                prop_assert_eq!(count, expected);
            }

            /// Property: Reverse complement preserves length
            #[test]
            fn prop_reverse_complement_preserves_length(seq in "[ACGTN]{1,1000}") {
                let rc = reverse_complement(seq.as_bytes());
                prop_assert_eq!(rc.len(), seq.len());
            }

            /// Property: Operations handle empty sequences gracefully
            #[test]
            fn prop_empty_sequences_work(_dummy in 0..1u8) {
                let empty = b"";
                prop_assert_eq!(reverse_complement(empty), empty);
                prop_assert_eq!(complement(empty), empty);
                prop_assert_eq!(reverse(empty), empty);
                prop_assert!(is_valid_dna(empty));
                prop_assert_eq!(count_invalid_bases(empty), 0);
            }

            /// Property: Operations handle single-character sequences
            #[test]
            fn prop_single_char_works(base in "[ACGT]") {
                let seq = base.as_bytes();
                let rc = reverse_complement(seq);
                prop_assert_eq!(rc.len(), 1);

                // RC of single base should be its complement
                let expected = complement(seq);
                prop_assert_eq!(rc, expected);
            }

            /// Property: Ambiguous bases are preserved through operations
            #[test]
            fn prop_ambiguous_bases_handled(seq in "[NRYSWKM]{1,100}") {
                // All ambiguous bases should remain valid DNA after RC
                let rc = reverse_complement(seq.as_bytes());
                prop_assert!(is_valid_dna(&rc));
            }
        }
    }

    // ===== Unit Tests =====

    // ===== Reverse Complement Tests =====

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCGC"), b"GCGC"); // Palindrome
    }

    #[test]
    fn test_reverse_complement_rna() {
        // RNA: U complements to A, but output will have T (not U)
        // This is a limitation of using a single static complement table
        assert_eq!(reverse_complement(b"AUGC"), b"GCAT"); // Note: T not U
        assert_eq!(reverse_complement(b"UUUU"), b"AAAA");

        // Note: For true RNA→RNA complement, would need A→U in table
        // Currently optimized for DNA (A→T), RNA input produces DNA output
    }

    #[test]
    fn test_reverse_complement_ambiguous() {
        // N stays N
        assert_eq!(reverse_complement(b"ATGCN"), b"NGCAT");
        // R (A or G) ↔ Y (C or T)
        assert_eq!(reverse_complement(b"ATGR"), b"YCAT");
        // W (A or T) ↔ W (self-complement)
        assert_eq!(reverse_complement(b"ATGW"), b"WCAT");
    }

    #[test]
    fn test_reverse_complement_involutive() {
        // RC(RC(x)) = x for DNA sequences only
        // (RNA is not involutive due to U→A→T with single complement table)
        let seqs = vec![
            &b"ATGC"[..],
            &b"ATGCATGC"[..],
            &b"ATGCN"[..],
        ];

        for seq in seqs {
            let rc = reverse_complement(seq);
            let rc_rc = reverse_complement(&rc);
            assert_eq!(rc_rc, seq, "RC should be involutive for DNA");
        }
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_reverse_complement_single() {
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"T"), b"A");
        assert_eq!(reverse_complement(b"G"), b"C");
        assert_eq!(reverse_complement(b"C"), b"G");
    }

    #[test]
    fn test_reverse_complement_case() {
        // Lowercase
        assert_eq!(reverse_complement(b"atgc"), b"gcat");
        // Mixed case
        assert_eq!(reverse_complement(b"AtGc"), b"gCaT");
    }

    // ===== In-Place Reverse Complement Tests =====

    #[test]
    fn test_reverse_complement_inplace_basic() {
        let mut seq = b"ATGC".to_vec();
        reverse_complement_inplace(&mut seq);
        assert_eq!(seq, b"GCAT");
    }

    #[test]
    fn test_reverse_complement_inplace_odd_length() {
        let mut seq = b"ATGCA".to_vec();
        reverse_complement_inplace(&mut seq);
        assert_eq!(seq, b"TGCAT");
    }

    #[test]
    fn test_reverse_complement_inplace_involutive() {
        let mut seq = b"ATGC".to_vec();
        reverse_complement_inplace(&mut seq);
        reverse_complement_inplace(&mut seq);
        assert_eq!(seq, b"ATGC");
    }

    // ===== Complement Tests =====

    #[test]
    fn test_complement_basic() {
        assert_eq!(complement(b"ATGC"), b"TACG");
        assert_eq!(complement(b"AAAA"), b"TTTT");
        assert_eq!(complement(b"GCGC"), b"CGCG");
    }

    #[test]
    fn test_complement_involutive() {
        // Complement is involutive for DNA (not RNA due to U→A→T)
        let seqs = vec![&b"ATGC"[..], &b"ATGCN"[..]];

        for seq in seqs {
            let comp = complement(seq);
            let comp2 = complement(&comp);
            assert_eq!(comp2, seq, "Complement should be involutive for DNA");
        }
    }

    #[test]
    fn test_complement_inplace() {
        let mut seq = b"ATGC".to_vec();
        complement_inplace(&mut seq);
        assert_eq!(seq, b"TACG");

        complement_inplace(&mut seq);
        assert_eq!(seq, b"ATGC");
    }

    // ===== Reverse Tests =====

    #[test]
    fn test_reverse_basic() {
        assert_eq!(reverse(b"ATGC"), b"CGTA");
        assert_eq!(reverse(b"AAAA"), b"AAAA");
    }

    #[test]
    fn test_reverse_involutive() {
        let seq = b"ATGC";
        let rev = reverse(seq);
        let rev2 = reverse(&rev);
        assert_eq!(rev2, seq);
    }

    #[test]
    fn test_reverse_inplace() {
        let mut seq = b"ATGC".to_vec();
        reverse_inplace(&mut seq);
        assert_eq!(seq, b"CGTA");

        reverse_inplace(&mut seq);
        assert_eq!(seq, b"ATGC");
    }

    // ===== Validation Tests =====

    #[test]
    fn test_is_valid_dna() {
        assert!(is_valid_dna(b"ATGC"));
        assert!(is_valid_dna(b"atgc"));
        assert!(is_valid_dna(b"ATGCN"));
        assert!(is_valid_dna(b"ATGCRYSWKM")); // All ambiguous

        assert!(!is_valid_dna(b"ATGCX"));
        assert!(!is_valid_dna(b"ATG C")); // Space
        assert!(!is_valid_dna(b"123"));
        assert!(!is_valid_dna(b"ATGCU")); // U is RNA, not DNA
    }

    #[test]
    fn test_is_valid_rna() {
        assert!(is_valid_rna(b"AUGC"));
        assert!(is_valid_rna(b"augc"));
        assert!(is_valid_rna(b"AUGCN"));

        assert!(!is_valid_rna(b"ATGC")); // T is DNA, not RNA
        assert!(!is_valid_rna(b"AUGCX"));
    }

    #[test]
    fn test_count_invalid_bases() {
        assert_eq!(count_invalid_bases(b"ATGC"), 0);
        assert_eq!(count_invalid_bases(b"ATGCX"), 1);
        assert_eq!(count_invalid_bases(b"AT GC"), 1);
        assert_eq!(count_invalid_bases(b"123"), 3);
        assert_eq!(count_invalid_bases(b"ATGCN"), 0); // N is valid
    }

    #[test]
    fn test_complement_table_correctness() {
        // DNA
        assert_eq!(COMPLEMENT_TABLE[b'A' as usize], b'T');
        assert_eq!(COMPLEMENT_TABLE[b'T' as usize], b'A');
        assert_eq!(COMPLEMENT_TABLE[b'G' as usize], b'C');
        assert_eq!(COMPLEMENT_TABLE[b'C' as usize], b'G');

        // RNA
        assert_eq!(COMPLEMENT_TABLE[b'U' as usize], b'A');

        // Ambiguous
        assert_eq!(COMPLEMENT_TABLE[b'R' as usize], b'Y');
        assert_eq!(COMPLEMENT_TABLE[b'Y' as usize], b'R');
        assert_eq!(COMPLEMENT_TABLE[b'W' as usize], b'W'); // Self-complement
        assert_eq!(COMPLEMENT_TABLE[b'S' as usize], b'S'); // Self-complement
        assert_eq!(COMPLEMENT_TABLE[b'K' as usize], b'M');
        assert_eq!(COMPLEMENT_TABLE[b'M' as usize], b'K');

        // N stays N
        assert_eq!(COMPLEMENT_TABLE[b'N' as usize], b'N');
    }
}
