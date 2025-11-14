//! Python wrappers for GFF3 format (gene annotations)

use pyo3::prelude::*;
use crate::formats::gff::{Gff3Record, Gff3Parser};
use crate::formats::primitives::Strand;
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

/// GFF3 feature record
///
/// Represents a genomic feature annotation (gene, mRNA, exon, CDS, etc.)
/// Coordinates are 1-based, inclusive [start, end].
///
/// Attributes:
///     seqid (str): Sequence/chromosome name
///     source (str): Annotation source (e.g., "Ensembl", "RefSeq")
///     feature_type (str): Feature type (e.g., "gene", "mRNA", "exon")
///     start (int): Start position (1-based, inclusive)
///     end (int): End position (1-based, inclusive)
///     score (float | None): Feature score
///     strand (str): Strand ('+', '-', or '.')
///     phase (int | None): Reading frame phase (0, 1, or 2)
///     attributes (dict[str, str]): Feature attributes (ID, Parent, Name, etc.)
///
/// Example:
///     >>> line = "chr1\tEnsembl\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=ABC1"
///     >>> record = Gff3Record.from_line(line)
///     >>> print(f"{record.get_name()}: {record.feature_type} on {record.seqid}")
///     ABC1: gene on chr1
///     >>> print(f"Location: {record.start}-{record.end} ({record.strand})")
///     Location: 1000-5000 (+)
#[pyclass(name = "Gff3Record")]
#[derive(Clone)]
pub struct PyGff3Record {
    #[pyo3(get)]
    pub seqid: String,
    #[pyo3(get)]
    pub source: String,
    #[pyo3(get)]
    pub feature_type: String,
    #[pyo3(get)]
    pub start: u64,
    #[pyo3(get)]
    pub end: u64,
    #[pyo3(get)]
    pub score: Option<f64>,
    #[pyo3(get)]
    pub strand: String,
    #[pyo3(get)]
    pub phase: Option<u8>,
    #[pyo3(get)]
    pub attributes: HashMap<String, String>,
}

#[pymethods]
impl PyGff3Record {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let record = Gff3Record::from_line(line)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(record.into())
    }

    fn to_line(&self) -> String {
        let strand = match self.strand.as_str() {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        };

        Gff3Record {
            seqid: self.seqid.clone(),
            source: self.source.clone(),
            feature_type: self.feature_type.clone(),
            start: self.start,
            end: self.end,
            score: self.score,
            strand,
            phase: self.phase,
            attributes: self.attributes.clone(),
        }
        .to_line()
    }

    fn length(&self) -> u64 {
        self.end - self.start + 1 // GFF3 is 1-based inclusive
    }

    fn get_id(&self) -> Option<String> {
        self.attributes.get("ID").cloned()
    }

    fn get_parent(&self) -> Option<String> {
        self.attributes.get("Parent").cloned()
    }

    fn get_name(&self) -> Option<String> {
        self.attributes.get("Name").cloned()
    }

    /// Convert to 0-based half-open interval [start, end)
    ///
    /// Returns:
    ///     tuple[int, int]: (start, end) in 0-based coordinates
    ///
    /// Example:
    ///     >>> record = Gff3Record.from_line("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1")
    ///     >>> start, end = record.to_0based()
    ///     >>> print(f"0-based: [{start}, {end})")
    ///     0-based: [999, 2000)
    fn to_0based(&self) -> (u64, u64) {
        (self.start - 1, self.end)
    }

    fn __repr__(&self) -> String {
        format!(
            "Gff3Record(seqid='{}', type='{}', start={}, end={})",
            self.seqid, self.feature_type, self.start, self.end
        )
    }

    fn __str__(&self) -> String {
        let name = self.get_name().unwrap_or_else(|| self.get_id().unwrap_or_else(|| "?".to_string()));
        format!(
            "{} ({}): {}:{}-{} [{}]",
            name, self.feature_type, self.seqid, self.start, self.end, self.strand
        )
    }
}

impl From<Gff3Record> for PyGff3Record {
    fn from(record: Gff3Record) -> Self {
        PyGff3Record {
            seqid: record.seqid,
            source: record.source,
            feature_type: record.feature_type,
            start: record.start,
            end: record.end,
            score: record.score,
            strand: record.strand.to_string(),
            phase: record.phase,
            attributes: record.attributes,
        }
    }
}

/// Stream GFF3 records with constant memory
///
/// Streaming iterator that processes GFF3 annotation files one record at a time.
/// Automatically skips header lines (##) and directive lines.
///
/// Args:
///     path (str): Path to GFF3 file
///
/// Example:
///     >>> stream = Gff3Stream.from_path("annotations.gff3")
///     >>> genes = []
///     >>> exons = []
///     >>> for record in stream:
///     ...     if record.feature_type == "gene":
///     ...         genes.append(record)
///     ...     elif record.feature_type == "exon":
///     ...         exons.append(record)
///     >>> print(f"Found {len(genes)} genes, {len(exons)} exons")
///
///     # Hierarchical analysis
///     >>> stream = Gff3Stream.from_path("annotations.gff3")
///     >>> features = {record.get_id(): record for record in stream if record.get_id()}
///     >>> for record in features.values():
///     ...     if record.feature_type == "mRNA":
///     ...         parent_id = record.get_parent()
///     ...         if parent_id and parent_id in features:
///     ...             gene = features[parent_id]
///     ...             print(f"Transcript {record.get_name()} belongs to gene {gene.get_name()}")
#[pyclass(name = "Gff3Stream", unsendable)]
pub struct PyGff3Stream {
    inner: Option<Gff3Parser<Box<dyn Read>>>,
}

#[pymethods]
impl PyGff3Stream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let parser = Gff3Parser::new(reader);

        Ok(PyGff3Stream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyGff3Record> {
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
        "Gff3Stream(...)".to_string()
    }
}
