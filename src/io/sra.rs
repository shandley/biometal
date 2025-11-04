//! SRA (Sequence Read Archive) integration
//!
//! This module provides utilities for streaming data directly from NCBI's
//! Sequence Read Archive without downloading entire datasets.
//!
//! # Evidence
//!
//! Entry 028 (Lab Notebook):
//! - Network streaming addresses I/O bottleneck (264-352× slower than compute)
//! - Enables analysis of large SRA datasets without local download
//! - Critical for democratizing bioinformatics access
//!
//! # Supported Accessions
//!
//! - **SRR** (Run): Most common, represents a sequencing run
//! - **SRX** (Experiment): Collection of runs
//! - **SRS** (Sample): Biological sample
//! - **SRP** (Project/Study): Collection of experiments
//!
//! # Important Limitation
//!
//! **SRA files are in NCBI's proprietary binary format**, not FASTQ. Direct streaming
//! from SRA requires the SRA toolkit to decode the format. This module provides
//! URL conversion utilities, but actual FASTQ streaming from SRA is not currently
//! supported without additional tooling.
//!
//! For direct HTTP FASTQ streaming, use `DataSource::Http` with FASTQ.gz URLs instead.
//!
//! # Example (Experimental)
//!
//! ```no_run
//! use biometal::io::sra_to_url;
//!
//! # fn main() -> biometal::Result<()> {
//! // Get SRA file URL (note: returns SRA format, not FASTQ)
//! let url = sra_to_url("SRR000001")?;
//! println!("SRA file URL: {}", url);
//!
//! // For actual FASTQ streaming, use direct FASTQ.gz URLs:
//! // let source = DataSource::Http("https://example.com/data.fastq.gz".to_string());
//! # Ok(())
//! # }
//! ```
//!
//! # URL Patterns
//!
//! NCBI provides public S3 access:
//! ```text
//! SRR000001 → https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000001/SRR000001
//! SRR390728 → https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390728/SRR390728
//! ```
//!
//! **Note**: As of 2024, NCBI simplified their S3 structure to use direct accession paths
//! without intermediate directory prefixes. This pattern is stable and documented in:
//! https://www.ncbi.nlm.nih.gov/sra/docs/sra-aws-download/

use crate::error::{BiometalError, Result};

/// Base URL for NCBI SRA public S3 access
const SRA_BASE_URL: &str = "https://sra-pub-run-odp.s3.amazonaws.com/sra";

/// Convert an SRA accession to an HTTP URL for streaming
///
/// Supports SRR (run), SRX (experiment), SRS (sample), and SRP (study) accessions.
/// Returns the URL for the SRA file which can be streamed via HTTP range requests.
///
/// # Arguments
///
/// * `accession` - SRA accession (e.g., "SRR000001", "SRX123456")
///
/// # Returns
///
/// HTTP URL for streaming the SRA data
///
/// # Errors
///
/// Returns error if accession format is invalid
///
/// # Example
///
/// ```
/// # use biometal::io::sra::sra_to_url;
/// let url = sra_to_url("SRR000001")?;
/// assert_eq!(url, "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000001/SRR000001");
///
/// let url = sra_to_url("SRR390728")?;
/// assert_eq!(url, "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390728/SRR390728");
/// # Ok::<(), biometal::BiometalError>(())
/// ```
pub fn sra_to_url(accession: &str) -> Result<String> {
    // Validate accession format
    if accession.len() < 7 {
        return Err(BiometalError::Network(format!(
            "Invalid SRA accession '{}': too short (expected at least 7 characters)",
            accession
        )));
    }

    // Check prefix (SRR, SRX, SRS, SRP)
    let prefix = &accession[..3];
    if !matches!(prefix, "SRR" | "SRX" | "SRS" | "SRP") {
        return Err(BiometalError::Network(format!(
            "Invalid SRA accession '{}': must start with SRR, SRX, SRS, or SRP",
            accession
        )));
    }

    // Check that remaining characters are digits
    if !accession[3..].chars().all(|c| c.is_ascii_digit()) {
        return Err(BiometalError::Network(format!(
            "Invalid SRA accession '{}': must be prefix + digits (e.g., SRR000001)",
            accession
        )));
    }

    // Build URL: base/accession/accession
    // Note: NCBI changed their S3 structure - no longer uses 6-char prefix directory
    let url = format!("{}/{}/{}", SRA_BASE_URL, accession, accession);

    Ok(url)
}

/// Check if a string looks like an SRA accession
///
/// Returns true if the string matches SRA accession patterns (SRR, SRX, SRS, SRP).
/// Does not validate that the accession exists in SRA database.
///
/// # Example
///
/// ```
/// # use biometal::io::sra::is_sra_accession;
/// assert!(is_sra_accession("SRR000001"));
/// assert!(is_sra_accession("SRX123456"));
/// assert!(!is_sra_accession("https://example.com/file.fq"));
/// assert!(!is_sra_accession("file.fastq"));
/// ```
pub fn is_sra_accession(s: &str) -> bool {
    if s.len() < 7 {
        return false;
    }

    let prefix = &s[..3];
    if !matches!(prefix, "SRR" | "SRX" | "SRS" | "SRP") {
        return false;
    }

    s[3..].chars().all(|c| c.is_ascii_digit())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sra_to_url_basic() {
        let url = sra_to_url("SRR000001").unwrap();
        assert_eq!(
            url,
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR000001/SRR000001"
        );
    }

    #[test]
    fn test_sra_to_url_six_digits() {
        let url = sra_to_url("SRR123456").unwrap();
        assert_eq!(
            url,
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR123456/SRR123456"
        );
    }

    #[test]
    fn test_sra_to_url_seven_digits() {
        let url = sra_to_url("SRR1234567").unwrap();
        assert_eq!(
            url,
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1234567/SRR1234567"
        );
    }

    #[test]
    fn test_sra_to_url_different_prefixes() {
        let url = sra_to_url("SRX000001").unwrap();
        assert!(url.contains("/SRX000001/SRX000001"));

        let url = sra_to_url("SRS000001").unwrap();
        assert!(url.contains("/SRS000001/SRS000001"));

        let url = sra_to_url("SRP000001").unwrap();
        assert!(url.contains("/SRP000001/SRP000001"));
    }

    #[test]
    fn test_sra_to_url_invalid_too_short() {
        let result = sra_to_url("SRR001");
        assert!(result.is_err());
        match result.unwrap_err() {
            BiometalError::Network(msg) => assert!(msg.contains("too short")),
            _ => panic!("Expected Network error"),
        }
    }

    #[test]
    fn test_sra_to_url_invalid_prefix() {
        let result = sra_to_url("ABC123456");
        assert!(result.is_err());
        match result.unwrap_err() {
            BiometalError::Network(msg) => assert!(msg.contains("must start with")),
            _ => panic!("Expected Network error"),
        }
    }

    #[test]
    fn test_sra_to_url_invalid_non_numeric() {
        let result = sra_to_url("SRR00000A");
        assert!(result.is_err());
        match result.unwrap_err() {
            BiometalError::Network(msg) => assert!(msg.contains("must be prefix + digits")),
            _ => panic!("Expected Network error"),
        }
    }

    #[test]
    fn test_is_sra_accession_valid() {
        assert!(is_sra_accession("SRR000001"));
        assert!(is_sra_accession("SRR123456"));
        assert!(is_sra_accession("SRX000001"));
        assert!(is_sra_accession("SRS000001"));
        assert!(is_sra_accession("SRP000001"));
        assert!(is_sra_accession("SRR1234567890")); // Long accession
    }

    #[test]
    fn test_is_sra_accession_invalid() {
        assert!(!is_sra_accession("SRR001")); // Too short
        assert!(!is_sra_accession("ABC123456")); // Wrong prefix
        assert!(!is_sra_accession("SRR00000A")); // Non-numeric
        assert!(!is_sra_accession("https://example.com")); // URL
        assert!(!is_sra_accession("file.fastq")); // Filename
        assert!(!is_sra_accession("")); // Empty
    }
}
