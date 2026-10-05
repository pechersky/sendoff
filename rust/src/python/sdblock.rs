use std::fs::File;
use std::io::{ErrorKind, Read};
use std::path::PathBuf;

use pyo3::exceptions::{
    PyNotImplementedError, PyOSError, PyTypeError, PyUnicodeDecodeError, PyValueError,
};
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyIterator, PyList, PyString, PyTuple};

use crate::core::block::{self, Block, Record};
use crate::core::scan::{Scanner, translate_newlines};
use crate::python::ctable::{PyCTable, error_to_py};

/// Characters or bytes requested per `read` call.
const CHUNK: usize = 1 << 20;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Text,
    Bytes,
}

impl Kind {
    fn value<'py>(self, py: Python<'py>, bytes: &[u8]) -> PyResult<Bound<'py, PyAny>> {
        match self {
            Kind::Text => PyString::from_bytes(py, bytes).map(Bound::into_any),
            Kind::Bytes => Ok(PyBytes::new(py, bytes).into_any()),
        }
    }
}

/// Return the kind and UTF-8 or raw bytes of a str or bytes object.
fn bytes_of<'a>(value: &'a Bound<'_, PyAny>) -> PyResult<Option<(Kind, &'a [u8])>> {
    if let Ok(text) = value.cast::<PyString>() {
        Ok(Some((Kind::Text, text.to_str()?.as_bytes())))
    } else if let Ok(bytes) = value.cast::<PyBytes>() {
        Ok(Some((Kind::Bytes, bytes.as_bytes())))
    } else {
        Ok(None)
    }
}

fn set_kind(kind: &mut Option<Kind>, found: Kind) -> PyResult<()> {
    match kind {
        Some(existing) if *existing != found => Err(PyTypeError::new_err(
            "SDBlock.from_lines does not accept mixed str and bytes lines",
        )),
        _ => {
            *kind = Some(found);
            Ok(())
        }
    }
}

enum Source {
    /// A Python iterable of str or bytes lines.
    Lines(Py<PyIterator>),
    /// A file object's bound `read1` or `read` method.
    Reader(Py<PyAny>),
    /// A path to open on the first read; `text` translates newlines.
    Path {
        path: PathBuf,
        filename: Py<PyAny>,
        text: bool,
    },
    /// A file opened from `Path`, read without holding the GIL.
    File {
        file: File,
        filename: Py<PyAny>,
        text: bool,
        pending_cr: bool,
    },
    Closed,
}

/// Build the `OSError` subclass Python's `open` would raise.
fn os_error(py: Python<'_>, error: std::io::Error, filename: &Py<PyAny>) -> PyErr {
    let Some(errno) = error.raw_os_error() else {
        return error.into();
    };
    match PyModule::import(py, "os").and_then(|os| os.call_method1("strerror", (errno,))) {
        Ok(strerror) => PyOSError::new_err((errno, strerror.unbind(), filename.clone_ref(py))),
        Err(error) => error,
    }
}

impl Source {
    /// Append input to `buf`, returning `true` at end of input.
    fn fill(
        &mut self,
        py: Python<'_>,
        kind: &mut Option<Kind>,
        buf: &mut Vec<u8>,
    ) -> PyResult<bool> {
        match self {
            Source::Lines(lines) => {
                let mut lines = lines.bind(py).clone();
                loop {
                    let Some(line) = lines.next() else {
                        return Ok(true);
                    };
                    let line = line?;
                    let Some((found, bytes)) = bytes_of(&line)? else {
                        return Err(PyTypeError::new_err(
                            "SDBlock.from_lines lines must be str or bytes",
                        ));
                    };
                    set_kind(kind, found)?;
                    buf.extend_from_slice(bytes);
                    if !bytes.ends_with(b"\n") {
                        buf.push(b'\n');
                    }
                    // Stop at a delimiter so blocks are produced as lines arrive.
                    if bytes.starts_with(b"$$$$") {
                        return Ok(false);
                    }
                }
            }
            Source::Reader(read) => {
                let chunk = read.bind(py).call1((CHUNK,))?;
                let Some((found, bytes)) = bytes_of(&chunk)? else {
                    return Err(PyTypeError::new_err(
                        "SDBlock.from_lines file objects must read str or bytes",
                    ));
                };
                if bytes.is_empty() {
                    return Ok(true);
                }
                set_kind(kind, found)?;
                buf.extend_from_slice(bytes);
                Ok(false)
            }
            Source::Path {
                path,
                filename,
                text,
            } => {
                let file = File::open(&*path).map_err(|error| os_error(py, error, filename))?;
                *kind = Some(if *text { Kind::Text } else { Kind::Bytes });
                *self = Source::File {
                    file,
                    filename: filename.clone_ref(py),
                    text: *text,
                    pending_cr: false,
                };
                self.fill(py, kind, buf)
            }
            Source::File {
                file,
                filename,
                text,
                pending_cr,
            } => {
                let start = buf.len();
                buf.resize(start + CHUNK, 0);
                let read = loop {
                    match py.detach(|| file.read(&mut buf[start..])) {
                        Err(error) if error.kind() == ErrorKind::Interrupted => {
                            py.check_signals()?;
                        }
                        result => break result,
                    }
                };
                let read = read.map_err(|error| os_error(py, error, filename));
                buf.truncate(start + *read.as_ref().unwrap_or(&0));
                match read? {
                    0 => Ok(true),
                    _ if *text => {
                        translate_newlines(buf, start, pending_cr);
                        Ok(false)
                    }
                    _ => Ok(false),
                }
            }
            Source::Closed => Ok(true),
        }
    }
}

/// Lazily parse SDBlocks from an input source.
#[pyclass(name = "_SDBlockIterator", module = "sendoff.native")]
pub(crate) struct PySDBlockIterator {
    source: Source,
    kind: Option<Kind>,
    scanner: Scanner,
}

impl PySDBlockIterator {
    fn new(source: Source) -> Self {
        Self {
            source,
            kind: None,
            scanner: Scanner::default(),
        }
    }
}

#[pymethods]
impl PySDBlockIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PySDBlock>> {
        let Self {
            source,
            kind,
            scanner,
        } = self;
        match scanner.next_block(|buf| source.fill(py, kind, buf)) {
            Ok(Some(block)) => Ok(Some(PySDBlock {
                block,
                kind: kind.unwrap_or(Kind::Text),
            })),
            Ok(None) => {
                self.close();
                Ok(None)
            }
            Err(error) => {
                self.close();
                Err(error)
            }
        }
    }

    /// Stop iterating and release the input source, closing files that
    /// `parse_sdf` opened.
    fn close(&mut self) {
        self.source = Source::Closed;
        self.scanner = Scanner::default();
    }
}

/// Lazily iterate over one block's SD data records.
#[pyclass(name = "_SDRecordIterator", module = "sendoff.native")]
pub(crate) struct PySDRecordIterator {
    block: Py<PySDBlock>,
    pos: usize,
    finished: bool,
}

#[pymethods]
impl PySDRecordIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyTuple>>> {
        if self.finished {
            return Ok(None);
        }
        let block = self.block.borrow(py);
        let result = match block.block.next_record(&mut self.pos) {
            Ok(Some(record)) => block.record(py, record).map(Some),
            Ok(None) => Ok(None),
            Err(error) => Err(PyValueError::new_err(error.to_string())),
        };
        if !matches!(result, Ok(Some(_))) {
            self.finished = true;
        }
        result
    }
}

/// Handle a single molecule block from an SD file.
#[pyclass(name = "SDBlock", module = "sendoff.native")]
pub(crate) struct PySDBlock {
    block: Block,
    kind: Kind,
}

impl PySDBlock {
    fn record<'py>(&self, py: Python<'py>, record: Record<'_>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            [
                self.kind.value(py, record.name)?,
                self.kind.value(py, &record.value)?,
            ],
        )
    }

    fn require_text(&self) -> PyResult<()> {
        match self.kind {
            Kind::Text => Ok(()),
            Kind::Bytes => Err(PyNotImplementedError::new_err(
                "ctable operations are not supported for byte-backed SDBlocks",
            )),
        }
    }
}

fn decode<'a>(py: Python<'_>, bytes: &'a [u8]) -> PyResult<&'a str> {
    std::str::from_utf8(bytes)
        .map_err(|error| PyUnicodeDecodeError::new_err_from_utf8(py, bytes, error))
}

#[pymethods]
impl PySDBlock {
    /// Parse SDBlocks from homogeneous str or bytes lines or a file object.
    ///
    /// File objects are read in large chunks with `read1`, when available,
    /// or `read`, instead of line by line.
    ///
    /// Args:
    ///     lines: an Iterable of str or bytes lines separated by the $$$$
    ///         delimiter, or a text or binary file object
    ///
    /// Returns:
    ///     An iterator that produces owned SDBlocks as input is consumed
    #[staticmethod]
    fn from_lines(lines: &Bound<'_, PyAny>) -> PyResult<PySDBlockIterator> {
        let py = lines.py();
        for method in [intern!(py, "read1"), intern!(py, "read")] {
            if let Ok(read) = lines.getattr(method) {
                return Ok(PySDBlockIterator::new(Source::Reader(read.unbind())));
            }
        }
        Ok(PySDBlockIterator::new(Source::Lines(
            lines.try_iter()?.unbind(),
        )))
    }

    /// Return the stripped molecule title in the input type.
    #[getter]
    fn title<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind.value(py, self.block.title())
    }

    /// Return all SD data records.
    ///
    /// Returns:
    ///     Tuples of record name and value in the input type
    fn records<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty(py);
        let mut pos = 0;
        while let Some(record) = self
            .block
            .next_record(&mut pos)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
        {
            list.append(self.record(py, record)?)?;
        }
        Ok(list)
    }

    /// Lazily iterate over the SD data records.
    fn iter_records(slf: Bound<'_, Self>) -> PySDRecordIterator {
        PySDRecordIterator {
            block: slf.unbind(),
            pos: 0,
            finished: false,
        }
    }

    /// Write the block, ending with $$$$, to a file-like handle.
    ///
    /// Args:
    ///     outh: a text or binary file-like handle matching the input type
    fn write(&self, outh: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = outh.py();
        let mut output = Vec::with_capacity(self.block.serialized_len());
        self.block.serialize_into(&mut output);
        outh.call_method1(intern!(py, "write"), (self.kind.value(py, &output)?,))?;
        Ok(())
    }

    /// Append a field and value of the same type as the input block.
    ///
    /// Args:
    ///     record_name: field name
    ///     value: field value, with no terminating newline
    fn append_record(
        &mut self,
        record_name: &Bound<'_, PyAny>,
        value: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        match (bytes_of(record_name)?, bytes_of(value)?) {
            (Some((name_kind, name)), Some((value_kind, value)))
                if name_kind == self.kind && value_kind == self.kind =>
            {
                self.block.append_record(name, value);
                Ok(())
            }
            _ => Err(PyTypeError::new_err(match self.kind {
                Kind::Text => "text SDBlock values must be appended as str",
                Kind::Bytes => "byte-backed SDBlock values must be appended as bytes",
            })),
        }
    }

    /// Parse out the underlying CTable object from a text SDBlock.
    ///
    /// Byte-backed SDBlocks cannot be converted without decoding their title.
    fn ctable(&self, py: Python<'_>) -> PyResult<PyCTable> {
        self.require_text()?;
        let mut lines = vec![decode(py, self.block.title())?.to_owned()];
        for line in block::lines(self.block.mdl()) {
            lines.push(decode(py, line)?.to_owned());
        }
        PyCTable::from_lines(lines)
    }

    /// Get the number of atoms declared by the counts line.
    fn num_atoms(&self) -> PyResult<i64> {
        let (_, atoms, _) = self.block.counts().map_err(error_to_py)?;
        Ok(atoms)
    }

    /// Get the number of bonds declared by the counts line.
    fn num_bonds(&self) -> PyResult<i64> {
        let (_, _, bonds) = self.block.counts().map_err(error_to_py)?;
        Ok(bonds)
    }

    /// Renumber atom and bond indices to be 1-indexed and in order.
    ///
    /// Byte-backed blocks raise `NotImplementedError`.
    fn renumber_indices(&mut self, py: Python<'_>) -> PyResult<()> {
        let mut ctable = self.ctable(py)?;
        ctable.renumber_indices_for_block(py)?;
        let lines = ctable.lines_after_title();
        self.block.replace_mdl(lines.iter().map(String::as_bytes));
        Ok(())
    }
}

/// Lazily parse SDBlocks from an SDF path.
///
/// The file is opened on the first `next()` and read in chunks without
/// holding the GIL. str paths give str values decoded as UTF-8 on access,
/// with universal newlines; bytes paths give bytes split on b"\n" only.
///
/// Args:
///     path: a str, bytes, or os.PathLike path
///
/// Returns:
///     An iterator of SDBlocks; `close()` closes the file early
#[pyfunction]
pub(crate) fn parse_sdf(path: &Bound<'_, PyAny>) -> PyResult<PySDBlockIterator> {
    let os = PyModule::import(path.py(), "os")?;
    let filename = os.call_method1("fspath", (path,))?;
    let text = filename.is_instance_of::<PyString>();
    // fsdecode round-trips bytes paths losslessly through surrogateescape.
    let path = os.call_method1("fsdecode", (&filename,))?.extract()?;
    Ok(PySDBlockIterator::new(Source::Path {
        path,
        filename: filename.unbind(),
        text,
    }))
}
