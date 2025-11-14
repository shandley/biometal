//! Python wrappers for GFA format (assembly graphs)

use pyo3::prelude::*;
use crate::formats::gfa::{GfaSegment, GfaLink, GfaPath, GfaParser, Orientation};
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

/// GFA segment (contig/node in assembly graph)
///
/// Represents a sequence segment in an assembly graph.
///
/// Attributes:
///     name (str): Segment identifier
///     sequence (str): DNA sequence
///     tags (dict[str, str]): Optional tags (e.g., LN:i:1000 for length)
///
/// Example:
///     >>> seg = GfaSegment.from_line("S\tcontig1\tACGT\tLN:i:4")
///     >>> print(f"{seg.name}: {seg.length()} bp")
///     contig1: 4 bp
#[pyclass(name = "GfaSegment")]
#[derive(Clone)]
pub struct PyGfaSegment {
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub sequence: String,
    #[pyo3(get)]
    pub tags: HashMap<String, String>,
}

#[pymethods]
impl PyGfaSegment {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() || fields[0] != "S" {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Not a segment line"));
        }
        if fields.len() < 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid segment format"));
        }

        let name = fields[1].to_string();
        let sequence = fields[2].to_string();
        let mut tags = HashMap::new();

        for tag in fields.iter().skip(3) {
            if let Some((key, value)) = tag.split_once(':') {
                tags.insert(key.to_string(), value.to_string());
            }
        }

        Ok(PyGfaSegment { name, sequence, tags })
    }

    fn length(&self) -> Option<usize> {
        self.tags
            .get("LN")
            .and_then(|s| s.split(':').last())
            .and_then(|s| s.parse::<usize>().ok())
            .or_else(|| Some(self.sequence.len()))
    }

    fn coverage(&self) -> Option<f64> {
        self.tags
            .get("KC")
            .and_then(|s| s.split(':').last())
            .and_then(|s| s.parse::<f64>().ok())
    }

    fn __repr__(&self) -> String {
        format!("GfaSegment(name='{}', length={:?})", self.name, self.length())
    }

    fn __str__(&self) -> String {
        format!("{} ({} bp)", self.name, self.length().unwrap_or(0))
    }
}

impl From<GfaSegment> for PyGfaSegment {
    fn from(seg: GfaSegment) -> Self {
        PyGfaSegment {
            name: seg.name,
            sequence: seg.sequence,
            tags: seg.tags,
        }
    }
}

/// GFA link (edge in assembly graph)
///
/// Represents a connection between two segments.
///
/// Attributes:
///     from_segment (str): Source segment name
///     from_orient (str): Source orientation ('+' or '-')
///     to_segment (str): Target segment name
///     to_orient (str): Target orientation ('+' or '-')
///     overlap (str): Overlap CIGAR string
///     tags (dict[str, str]): Optional tags
///
/// Example:
///     >>> link = GfaLink.from_line("L\tcontig1\t+\tcontig2\t+\t10M")
///     >>> print(f"{link.from_segment} -> {link.to_segment}")
///     contig1 -> contig2
#[pyclass(name = "GfaLink")]
#[derive(Clone)]
pub struct PyGfaLink {
    #[pyo3(get)]
    pub from_segment: String,
    #[pyo3(get)]
    pub from_orient: String,
    #[pyo3(get)]
    pub to_segment: String,
    #[pyo3(get)]
    pub to_orient: String,
    #[pyo3(get)]
    pub overlap: String,
    #[pyo3(get)]
    pub tags: HashMap<String, String>,
}

#[pymethods]
impl PyGfaLink {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() || fields[0] != "L" {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Not a link line"));
        }
        if fields.len() < 6 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid link format"));
        }

        let from_segment = fields[1].to_string();
        let from_orient = fields[2].to_string();
        let to_segment = fields[3].to_string();
        let to_orient = fields[4].to_string();
        let overlap = fields[5].to_string();
        let mut tags = HashMap::new();

        for tag in fields.iter().skip(6) {
            if let Some((key, value)) = tag.split_once(':') {
                tags.insert(key.to_string(), value.to_string());
            }
        }

        Ok(PyGfaLink {
            from_segment,
            from_orient,
            to_segment,
            to_orient,
            overlap,
            tags,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "GfaLink(from='{}{}', to='{}{}')",
            self.from_segment, self.from_orient, self.to_segment, self.to_orient
        )
    }

    fn __str__(&self) -> String {
        format!("{}{} -> {}{}", self.from_segment, self.from_orient, self.to_segment, self.to_orient)
    }
}

impl From<GfaLink> for PyGfaLink {
    fn from(link: GfaLink) -> Self {
        PyGfaLink {
            from_segment: link.from_segment,
            from_orient: match link.from_orient {
                Orientation::Forward => "+".to_string(),
                Orientation::Reverse => "-".to_string(),
            },
            to_segment: link.to_segment,
            to_orient: match link.to_orient {
                Orientation::Forward => "+".to_string(),
                Orientation::Reverse => "-".to_string(),
            },
            overlap: link.overlap,
            tags: link.tags,
        }
    }
}

/// GFA path (walk through assembly graph)
///
/// Represents an ordered path through segments.
///
/// Attributes:
///     name (str): Path identifier
///     segments (list[str]): Ordered segment names
///     overlaps (list[str]): Overlap CIGARs between segments
///     tags (dict[str, str]): Optional tags
///
/// Example:
///     >>> path = GfaPath.from_line("P\tpath1\tcontig1+,contig2+\t10M,5M")
///     >>> print(f"Path {path.name}: {len(path.segments)} segments")
///     Path path1: 2 segments
#[pyclass(name = "GfaPath")]
#[derive(Clone)]
pub struct PyGfaPath {
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub segments: Vec<String>,
    #[pyo3(get)]
    pub overlaps: Vec<String>,
    #[pyo3(get)]
    pub tags: HashMap<String, String>,
}

#[pymethods]
impl PyGfaPath {
    #[staticmethod]
    fn from_line(line: &str) -> PyResult<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() || fields[0] != "P" {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Not a path line"));
        }
        if fields.len() < 4 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid path format"));
        }

        let name = fields[1].to_string();
        let segments: Vec<String> = fields[2].split(',').map(|s| s.to_string()).collect();
        let overlaps: Vec<String> = fields[3].split(',').map(|s| s.to_string()).collect();
        let mut tags = HashMap::new();

        for tag in fields.iter().skip(4) {
            if let Some((key, value)) = tag.split_once(':') {
                tags.insert(key.to_string(), value.to_string());
            }
        }

        Ok(PyGfaPath {
            name,
            segments,
            overlaps,
            tags,
        })
    }

    fn length(&self) -> usize {
        self.segments.len()
    }

    fn __repr__(&self) -> String {
        format!("GfaPath(name='{}', segments={})", self.name, self.segments.len())
    }

    fn __str__(&self) -> String {
        format!("{}: {}", self.name, self.segments.join(" -> "))
    }
}

impl From<GfaPath> for PyGfaPath {
    fn from(path: GfaPath) -> Self {
        PyGfaPath {
            name: path.name,
            segments: path.segments,
            overlaps: path.overlaps,
            tags: path.tags,
        }
    }
}

/// Stream GFA records with constant memory
///
/// Streaming iterator that processes GFA assembly graphs one record at a time.
/// Returns segments (S), links (L), and paths (P).
///
/// Args:
///     path (str): Path to GFA file
///
/// Example:
///     >>> stream = GfaStream.from_path("assembly.gfa")
///     >>> segments = []
///     >>> links = []
///     >>> paths = []
///     >>> for record in stream:
///     ...     if isinstance(record, GfaSegment):
///     ...         segments.append(record)
///     ...     elif isinstance(record, GfaLink):
///     ...         links.append(record)
///     ...     elif isinstance(record, GfaPath):
///     ...         paths.append(record)
#[pyclass(name = "GfaStream", unsendable)]
pub struct PyGfaStream {
    inner: Option<GfaParser<Box<dyn Read>>>,
}

#[pymethods]
impl PyGfaStream {
    #[staticmethod]
    fn from_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let reader = open_file(&path_buf)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        let parser = GfaParser::new(reader);

        Ok(PyGfaStream {
            inner: Some(parser),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<PyObject> {
        let py = slf.py();
        if let Some(ref mut parser) = slf.inner {
            // Loop to skip header records and comments
            loop {
                match parser.next() {
                    Some(Ok(record)) => {
                        match record {
                            crate::formats::gfa::GfaRecord::Segment(seg) => {
                                let py_seg: PyGfaSegment = seg.into();
                                return Ok(Bound::new(py, py_seg)?.into_any().unbind());
                            }
                            crate::formats::gfa::GfaRecord::Link(link) => {
                                let py_link: PyGfaLink = link.into();
                                return Ok(Bound::new(py, py_link)?.into_any().unbind());
                            }
                            crate::formats::gfa::GfaRecord::Path(path) => {
                                let py_path: PyGfaPath = path.into();
                                return Ok(Bound::new(py, py_path)?.into_any().unbind());
                            }
                            crate::formats::gfa::GfaRecord::Header(_) => {
                                // Skip headers, continue loop
                                continue;
                            }
                            crate::formats::gfa::GfaRecord::Comment(_) => {
                                // Skip comments, continue loop
                                continue;
                            }
                        }
                    }
                    Some(Err(e)) => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())),
                    None => return Err(pyo3::exceptions::PyStopIteration::new_err("no more records")),
                }
            }
        } else {
            Err(pyo3::exceptions::PyStopIteration::new_err("stream exhausted"))
        }
    }

    fn __repr__(&self) -> String {
        "GfaStream(...)".to_string()
    }
}
