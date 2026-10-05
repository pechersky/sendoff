use std::borrow::Cow;
use std::fs::File;
use std::io::{ErrorKind, Read};
use std::path::PathBuf;

use pyo3::exceptions::{
    PyNotImplementedError, PyOSError, PyTypeError, PyUnicodeDecodeError, PyValueError,
};
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyGenericAlias, PyIterator, PyList, PyString, PyTuple, PyType};

use crate::core::block::{self, Block, Record, Strip};
use crate::core::scan::{Scanner, TRANSLATE_SLACK, Tail, single_block, translate_newlines};
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

    /// The Python `strip()` that values of this kind follow.
    fn strip(self) -> Strip {
        match self {
            Kind::Text => Strip::Str,
            Kind::Bytes => Strip::Bytes,
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

/// Translate `\r\n` and lone `\r` to `\n`, as text files read in Python.
fn universal_newlines(bytes: &[u8]) -> Cow<'_, [u8]> {
    if memchr::memchr(b'\r', bytes).is_none() {
        return Cow::Borrowed(bytes);
    }
    let mut buffer = Vec::with_capacity(bytes.len() + TRANSLATE_SLACK);
    buffer.extend_from_slice(bytes);
    buffer.resize(bytes.len() + TRANSLATE_SLACK, 0);
    let len = translate_newlines(&mut buffer, bytes.len(), &mut false);
    buffer.truncate(len);
    Cow::Owned(buffer)
}

/// Collect `(name, value)` pairs from a mapping or an iterable of pairs.
fn record_pairs<'py>(
    records: &Bound<'py, PyAny>,
) -> PyResult<Vec<(Bound<'py, PyAny>, Bound<'py, PyAny>)>> {
    let py = records.py();
    let items = match records.getattr(intern!(py, "items")) {
        Ok(items) => items.call0()?,
        Err(_) => records.clone(),
    };
    items.try_iter()?.map(|item| item?.extract()).collect()
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
    /// An in-memory str or bytes, copied out a chunk at a time from `pos`;
    /// `text` translates newlines.
    Buffer {
        data: Py<PyAny>,
        pos: usize,
        text: bool,
        pending_cr: bool,
    },
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
    /// Append input to `tail`, returning `true` at end of input.
    fn fill(
        &mut self,
        py: Python<'_>,
        kind: &mut Option<Kind>,
        tail: &mut Tail<'_>,
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
                    tail.extend(bytes);
                    if !bytes.ends_with(b"\n") {
                        tail.extend(b"\n");
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
                tail.extend(bytes);
                Ok(false)
            }
            Source::Buffer {
                data,
                pos,
                text,
                pending_cr,
            } => {
                let Some((_, bytes)) = bytes_of(data.bind(py))? else {
                    unreachable!("Source::Buffer holds str or bytes");
                };
                let chunk = &bytes[*pos..bytes.len().min(*pos + CHUNK)];
                *pos += chunk.len();
                if *text {
                    let space = tail.space(chunk.len() + TRANSLATE_SLACK);
                    space[..chunk.len()].copy_from_slice(chunk);
                    let len = translate_newlines(space, chunk.len(), pending_cr);
                    tail.commit(len);
                } else {
                    tail.extend(chunk);
                }
                Ok(*pos == bytes.len())
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
                self.fill(py, kind, tail)
            }
            Source::File {
                file,
                filename,
                text,
                pending_cr,
            } => {
                let space = tail.space(CHUNK + TRANSLATE_SLACK);
                let read = loop {
                    match py.detach(|| file.read(&mut space[..CHUNK])) {
                        Err(error) if error.kind() == ErrorKind::Interrupted => {
                            py.check_signals()?;
                        }
                        result => break result,
                    }
                };
                let read = read.map_err(|error| os_error(py, error, filename))?;
                let len = if *text {
                    translate_newlines(space, read, pending_cr)
                } else {
                    read
                };
                tail.commit(len);
                Ok(read == 0)
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

    fn of_kind(source: Source, kind: Kind) -> Self {
        Self {
            kind: Some(kind),
            ..Self::new(source)
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
        match scanner.next_block(|tail| source.fill(py, kind, tail)) {
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
        let result = match block.block.next_record(&mut self.pos, block.kind.strip()) {
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

    /// The bytes of a str or bytes argument, which must match the block's
    /// type; `what` names the argument in the error.
    fn arg_bytes<'a>(&self, value: &'a Bound<'_, PyAny>, what: &str) -> PyResult<&'a [u8]> {
        match bytes_of(value)? {
            Some((kind, bytes)) if kind == self.kind => Ok(bytes),
            _ => Err(PyTypeError::new_err(match self.kind {
                Kind::Text => format!("text SDBlock {what} must be str"),
                Kind::Bytes => format!("byte-backed SDBlock {what} must be bytes"),
            })),
        }
    }
}

fn decode<'a>(py: Python<'_>, bytes: &'a [u8]) -> PyResult<&'a str> {
    std::str::from_utf8(bytes)
        .map_err(|error| PyUnicodeDecodeError::new_err_from_utf8(py, bytes, error))
}

#[pymethods]
impl PySDBlock {
    /// Allow `SDBlock[str]` and `SDBlock[bytes]` in runtime annotations.
    #[classmethod]
    fn __class_getitem__<'py>(
        cls: &Bound<'py, PyType>,
        item: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyGenericAlias>> {
        PyGenericAlias::new(cls.py(), cls.as_any(), item)
    }

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
        if lines.is_instance_of::<PyString>() || lines.is_instance_of::<PyBytes>() {
            return Err(PyTypeError::new_err(
                "SDBlock.from_lines takes lines or a file object; \
                 use SDBlock.from_sdf to parse SDF str or bytes",
            ));
        }
        for method in [intern!(py, "read1"), intern!(py, "read")] {
            if let Ok(read) = lines.getattr(method) {
                return Ok(PySDBlockIterator::new(Source::Reader(read.unbind())));
            }
        }
        Ok(PySDBlockIterator::new(Source::Lines(
            lines.try_iter()?.unbind(),
        )))
    }

    /// Parse one block from SDF str or bytes, such as a molblock.
    ///
    /// The block ends at an optional `$$$$` line, after which only
    /// whitespace may follow. str input has universal newlines.
    ///
    /// Args:
    ///     sdf: a molblock or one SDF record, as str or bytes
    ///     records: SD data to append after any in `sdf`, as a mapping or
    ///         an iterable of (name, value) tuples of the same type
    #[new]
    #[pyo3(signature = (sdf, records=None))]
    fn new(sdf: &Bound<'_, PyAny>, records: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let Some((kind, bytes)) = bytes_of(sdf)? else {
            return Err(PyTypeError::new_err(
                "SDBlock() argument must be str or bytes",
            ));
        };
        let data = match kind {
            Kind::Text => universal_newlines(bytes).into_owned(),
            Kind::Bytes => bytes.to_vec(),
        };
        let block = single_block(data).ok_or_else(|| {
            PyValueError::new_err(
                "SDBlock() takes one block; use SDBlock.from_sdf to parse several",
            )
        })?;
        let mut slf = Self { block, kind };
        if let Some(records) = records {
            for (name, value) in record_pairs(records)? {
                let name = slf.arg_bytes(&name, "record names")?;
                let value = slf.arg_bytes(&value, "values")?;
                slf.block.append_record(name, value);
            }
        }
        Ok(slf)
    }

    /// Lazily parse SDBlocks from SDF str or bytes held in memory.
    ///
    /// Equivalent to `from_lines(io.StringIO(sdf))` or
    /// `from_lines(io.BytesIO(sdf))`, without splitting lines: str input
    /// has universal newlines, and bytes are split on b"\n" only.
    ///
    /// Args:
    ///     sdf: SDF content as str or bytes
    ///
    /// Returns:
    ///     An iterator of SDBlocks of the input type
    #[staticmethod]
    fn from_sdf(sdf: &Bound<'_, PyAny>) -> PyResult<PySDBlockIterator> {
        let Some((kind, _)) = bytes_of(sdf)? else {
            return Err(PyTypeError::new_err(
                "SDBlock.from_sdf argument must be str or bytes",
            ));
        };
        Ok(PySDBlockIterator::of_kind(
            Source::Buffer {
                data: sdf.clone().unbind(),
                pos: 0,
                text: kind == Kind::Text,
                pending_cr: false,
            },
            kind,
        ))
    }

    /// Return the stripped molecule title in the input type.
    #[getter]
    fn title<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind.value(py, self.block.title(self.kind.strip()))
    }

    /// Replace the title line.
    #[setter(title)]
    fn set_title(&mut self, title: &Bound<'_, PyAny>) -> PyResult<()> {
        let title = self.arg_bytes(title, "titles")?;
        if title.iter().any(|&byte| byte == b'\n' || byte == b'\r') {
            return Err(PyValueError::new_err(
                "SDBlock titles must be a single line",
            ));
        }
        self.block.set_title(title);
        Ok(())
    }

    /// Return the title line and CTAB through `M  END` as one value.
    ///
    /// Suitable for RDKit's `Chem.MolFromMolBlock` or OpenEye's
    /// `OEReadMolFromBytes`.
    #[getter]
    fn molblock<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind.value(py, &self.block.molblock(self.kind.strip()))
    }

    /// Replace the title line and CTAB, keeping the SD data.
    #[setter(molblock)]
    fn set_molblock(&mut self, molblock: &Bound<'_, PyAny>) -> PyResult<()> {
        let molblock = self.arg_bytes(molblock, "molblocks")?;
        let molblock = match self.kind {
            Kind::Text => universal_newlines(molblock),
            Kind::Bytes => Cow::Borrowed(molblock),
        };
        self.block
            .set_molblock(&molblock)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Return whether `needle` occurs in the CTAB text after the title.
    ///
    /// Like `needle in "\\n".join(block.mdl)` without building the lines;
    /// text blocks are searched as UTF-8 without being decoded.
    ///
    /// Args:
    ///     needle: substring in the block's type
    fn mdl_contains(&self, needle: &Bound<'_, PyAny>) -> PyResult<bool> {
        let needle = self.arg_bytes(needle, "needles")?;
        let mdl = self.block.mdl();
        let mdl = mdl.strip_suffix(b"\n").unwrap_or(mdl);
        Ok(memchr::memmem::find(mdl, needle).is_some())
    }

    /// Return the CTAB lines after the title, through `M  END`, without
    /// line endings.
    #[getter]
    fn mdl<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let lines = block::lines(self.block.mdl())
            .map(|line| self.kind.value(py, line))
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, lines)
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
            .next_record(&mut pos, self.kind.strip())
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

    /// Return the value of the last record named `name`, or `default`.
    ///
    /// Equivalent to `dict(block.records()).get(name, default)`, without
    /// building the other values.
    ///
    /// Args:
    ///     name: record name in the block's type
    ///     default: returned when no record has that name
    ///
    /// Returns:
    ///     The record value in the block's type, or `default`
    #[pyo3(signature = (name, default=None))]
    fn get<'py>(
        &self,
        name: &Bound<'py, PyAny>,
        default: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Option<Bound<'py, PyAny>>> {
        let py = name.py();
        let name = self.arg_bytes(name, "record names")?;
        let strip = self.kind.strip();
        let mut pos = 0;
        let mut found = None;
        while let Some(record) = self
            .block
            .next_record(&mut pos, strip)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
        {
            if record.name == name {
                found = Some(record.value);
            }
        }
        match found {
            Some(value) => self.kind.value(py, &value).map(Some),
            None => Ok(default),
        }
    }

    /// Return the records as a dict, later duplicates winning.
    ///
    /// Equivalent to `dict(block.records())`.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        let strip = self.kind.strip();
        let mut pos = 0;
        while let Some(record) = self
            .block
            .next_record(&mut pos, strip)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
        {
            dict.set_item(
                self.kind.value(py, record.name)?,
                self.kind.value(py, &record.value)?,
            )?;
        }
        Ok(dict)
    }

    /// Return the block as SDF text ending with $$$$, in the input type.
    ///
    /// The inverse of `SDBlock(sdf)`, and what `write` writes.
    fn to_sdf<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind
            .value(py, &self.block.serialized(self.kind.strip()))
    }

    /// Write the block, ending with $$$$, to a file-like handle.
    ///
    /// Args:
    ///     outh: a text or binary file-like handle matching the input type
    fn write(&self, outh: &Bound<'_, PyAny>) -> PyResult<()> {
        let py = outh.py();
        let output = self.block.serialized(self.kind.strip());
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
        let name = self.arg_bytes(record_name, "record names")?;
        let value = self.arg_bytes(value, "values")?;
        self.block.append_record(name, value);
        Ok(())
    }

    /// Set a record, keeping the position of the first with that name.
    ///
    /// Later records with the same name are removed; if there is none, the
    /// record is appended.
    ///
    /// Args:
    ///     name: record name in the block's type
    ///     value: record value in the block's type
    fn set_record(&mut self, name: &Bound<'_, PyAny>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let name = self.arg_bytes(name, "record names")?;
        let value = self.arg_bytes(value, "values")?;
        self.block
            .set_record(name, value, self.kind.strip())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Remove every record named `name`.
    ///
    /// Args:
    ///     name: record name in the block's type
    ///
    /// Returns:
    ///     How many records were removed
    fn remove_record(&mut self, name: &Bound<'_, PyAny>) -> PyResult<usize> {
        let name = self.arg_bytes(name, "record names")?;
        self.block
            .remove_records(name, self.kind.strip())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Replace all SD data with the given records.
    ///
    /// Args:
    ///     records: a mapping, or an iterable of (name, value) tuples, in
    ///         the block's type
    fn replace_records(&mut self, records: &Bound<'_, PyAny>) -> PyResult<()> {
        let pairs = record_pairs(records)?;
        let records = pairs
            .iter()
            .map(|(name, value)| {
                Ok((
                    self.arg_bytes(name, "record names")?,
                    self.arg_bytes(value, "values")?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        self.block.replace_records(records);
        Ok(())
    }

    /// Parse out the underlying CTable object from a text SDBlock.
    ///
    /// Byte-backed SDBlocks cannot be converted without decoding their title.
    fn ctable(&self, py: Python<'_>) -> PyResult<PyCTable> {
        self.require_text()?;
        let mut lines = vec![decode(py, self.block.title(Strip::Str))?.to_owned()];
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
