//! Python wrappers for VCF format (variant calling)

use pyo3::prelude::*;
use crate::formats::vcf::{VcfHeader, VcfRecord, VcfParser};
use crate::formats::TabDelimitedRecord;
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::io::Read;
use flate2::read::MultiGzDecoder;

/// Helper function to open a file, detecting and handling gzip compression
fn open_file(path: &Path) -> std::io::Result<Box<dyn Read>> {
    let file = File::open(path)?;

    // Check if file is gzipped by extension
    if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// VCF file header
///
/// Contains metadata from VCF header lines (##fileformat, ##INFO, etc.)
///
/// Attributes:
///     fileformat (str): VCF version (e.g., "VCFv4.2")
///     info_fields (dict[str, str]): INFO field definitions
///     format_fields (dict[str, str]): FORMAT field definitions
///     filters (dict[str, str]): FILTER definitions
///     contigs (dict[str, int | None]): Contig names and lengths
///     samples (list[str]): Sample IDs
///
/// Example:
///     >>> parser = VcfStream.from_path("variants.vcf")
///     >>> header = parser.header()
///     >>> print(f"VCF version: {header.fileformat}")
///     >>> print(f"Samples: {', '.join(header.samples)}")
#[pyclass(name = "VcfHeader")]
#[derive(Clone)]
pub struct PyVcfHeader {
    #[pyo3(get)]
    pub fileformat: String,
    #[pyo3(get)]
    pub info_fields: HashMap<String, String>,
    #[pyo3(get)]
    pub format_fields: HashMap<String, String>,
    #[pyo3(get)]
    pub filters: HashMap<String, String>,
    #[pyo3(get)]
    pub contigs: HashMap<String, Option<u64>>,
    #[pyo3(get)]
    pub samples: Vec<String>,
}

#[pymethods]
impl PyVcfHeader {
    fn __repr__(&self) -> String {
        format!(
            "VcfHeader(fileformat='{}', samples={})",
            self.fileformat,
            self.samples.len()
        )
    }

    fn __str__(&self) -> String {
        format!(
            "VCF {}: {} samples, {} contigs",
            self.fileformat,
            self.samples.len(),
            self.contigs.len()
        )
    }
}

impl From<VcfHeader> for PyVcfHeader {
    fn from(header: VcfHeader) -> Self {
        PyVcfHeader {
            fileformat: header.fileformat,
            info_fields: header.info_fields,
            format_fields: header.format_fields,
            filters: header.filters,
            contigs: header.contigs,
            samples: header.samples,
        }
    }
}

/// VCF variant record
///
/// Represents a single genetic variant.
///
/// Attributes:
///     chrom (str): Chromosome name
///     pos (int): Position (1-based)
///     id (str | None): Variant ID (e.g., rs12345)
///     reference (str): Reference allele
///     alternate (list[str]): Alternate alleles
///     quality (float | None): Variant quality score
///     filter (str | None): Filter status (e.g., "PASS")
///     info (dict[str, str]): INFO field key-value pairs
///     format (str | None): FORMAT field specification
///     samples (list[str]): Sample genotype data
///
/// Example:
///     >>> record = VcfRecord.from_line("chr1\t12345\trs123\tA\tT\t30\tPASS\tDP=100")
///     >>> print(f"{record.chrom}:{record.pos} {record.reference}>{record.alternate[0]}")
///     chr1:12345 A>T
///     >>> if record.filter == "PASS":
///     ...     depth = record.info.get("DP")
///     ...     print(f"Depth: {depth}")
///     Depth: 100
#[pyclass(name = "VcfRecord")]
#[derive(Clone)]
pub struct PyVcfRecord {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub pos: u64,
    #[pyo3(get)]
    pub id: Option<String>,
    #[pyo3(get)]
    pub reference: String,
    #[pyo3(get)]
    pub alternate: Vec<String>,
    #[pyo3(get)]
    pub quality: Option<f64>,
    #[pyo3(get)]
    pub filter: Option<String>,
    #[pyo3(get)]
    pub info: HashMap<String, String>,
    #[pyo3(get)]
    pub format: Option<String>,
    #[pyo3(get)]
    pub samples: Vec<String>,
}

#[pymethods]
impl PyVcfRecord {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = VcfRecord::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    fn to_line(&self) -> String {
        VcfRecord {
            chrom: self.chrom.clone(),
            pos: self.pos,
            id: self.id.clone(),
            reference: self.reference.clone(),
            alternate: self.alternate.clone(),
            quality: self.quality,
            filter: self.filter.clone(),
            info: self.info.clone(),
            format: self.format.clone(),
            samples: self.samples.clone(),
        }
        .to_line()
    }

    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternate.iter().all(|a| a.len() == 1)
    }

    fn is_insertion(&self) -> bool {
        self.alternate.iter().any(|a| a.len() > self.reference.len())
    }

    fn is_deletion(&self) -> bool {
        self.alternate.iter().any(|a| a.len() < self.reference.len())
    }

    fn is_indel(&self) -> bool {
        self.is_insertion() || self.is_deletion()
    }

    fn __repr__(&self) -> String {
        format!(
            "VcfRecord(chrom='{}', pos={}, ref='{}', alt={:?})",
            self.chrom, self.pos, self.reference, self.alternate
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{}:{} {}->{}",
            self.chrom,
            self.pos,
            self.reference,
            self.alternate.join(",")
        )
    }
}

impl From<VcfRecord> for PyVcfRecord {
    fn from(record: VcfRecord) -> Self {
        PyVcfRecord {
            chrom: record.chrom,
            pos: record.pos,
            id: record.id,
            reference: record.reference,
            alternate: record.alternate,
            quality: record.quality,
            filter: record.filter,
            info: record.info,
            format: record.format,
            samples: record.samples,
        }
    }
}

/// Stream VCF records with constant memory
///
/// Streaming iterator that processes VCF variant files one record at a time.
/// Header must be parsed before iterating records.
///
/// Args:
///     path (str): Path to VCF file
///
/// Example:
///     >>> stream = VcfStream.from_path("variants.vcf")
///     >>> header = stream.header()
///     >>> print(f"VCF version: {header.fileformat}")
///     >>>
///     >>> snps = 0
///     >>> indels = 0
///     >>> for record in stream:
///     ...     if record.is_snp():
///     ...         snps += 1
///     ...     elif record.is_indel():
///     ...         indels += 1
///     >>> print(f"SNPs: {snps}, Indels: {indels}")
#[pyclass(name = "VcfStream", unsendable)]
pub struct PyVcfStream {
    inner: Option<VcfParser<Box<dyn Read>>>,
    header: Option<PyVcfHeader>,
}

#[pymethods]
impl PyVcfStream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let mut parser = VcfParser::new(reader);

        // Parse header immediately
        let header = parser
            .parse_header()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;

        Ok(PyVcfStream {
            inner: Some(parser),
            header: Some(header.into()),
        })
    }

    fn header(&self) -> PyResult<PyVcfHeader> {
        self.header
            .clone()
            .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("Header not parsed"))
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyVcfRecord> {
        if let Some(ref mut parser) = slf.inner {
            match parser.next() {
                Some(Ok(record)) => Ok(record.into()),
                Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                None => Err(pyo3::exceptions::PyStopIteration::new_err("no more records")),
            }
        } else {
            Err(pyo3::exceptions::PyStopIteration::new_err("stream exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "VcfStream(...)".to_string()
    }
}
