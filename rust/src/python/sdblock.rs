use std::borrow::Cow;
use std::cell::{RefCell, RefMut};
use std::collections::HashMap;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{ErrorKind, Read, Write};
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

    /// Call `f` with a `Namer` that reuses name objects across blocks, or
    /// creates new ones in a re-entrant call that finds the cache in use.
    fn with_names<R>(self, f: impl FnOnce(&mut Namer<'_>) -> R) -> R {
        NAMES.with(|cache| {
            f(&mut Namer {
                kind: self,
                cache: cache.try_borrow_mut().ok(),
                prev: NO_NAME,
            })
        })
    }

    /// The Python `strip()` that values of this kind follow.
    fn strip(self) -> Strip {
        match self {
            Kind::Text => Strip::Str,
            Kind::Bytes => Strip::Bytes,
        }
    }
}

/// Names cached per kind before the cache is emptied, which bounds what
/// one-off names can hold on to.
const MAX_NAMES: usize = 4096;
/// Longer names are rare and are not cached.
const MAX_CACHED_NAME: usize = 64;

/// A fast hash for short record names, eight bytes at a time.
#[derive(Default)]
struct NameHasher(u64);

impl Hasher for NameHasher {
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut word = [0; 8];
            word[..chunk.len()].copy_from_slice(chunk);
            self.0 = (self.0.rotate_left(5) ^ u64::from_le_bytes(word))
                .wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }

    fn finish(&self) -> u64 {
        // Fold the well-mixed high bits into the low bits buckets use.
        self.0 ^ (self.0 >> 32)
    }
}

/// Marks no cached name: the start of a block, or an unknown successor.
const NO_NAME: u32 = u32::MAX;

/// A cached record-name object.
struct Name {
    bytes: Box<[u8]>,
    object: Py<PyAny>,
    /// The name that followed this one when last seen, or `NO_NAME`.
    next: u32,
}

/// Cached record names of one kind. Each links to the name that last
/// followed it, so a block repeating an earlier record order finds its names
/// by comparing bytes instead of hashing them, even when some records are
/// missing or new.
struct Names {
    index: HashMap<Box<[u8]>, u32, BuildHasherDefault<NameHasher>>,
    names: Vec<Name>,
    /// The first name of the last block.
    first: u32,
}

impl Default for Names {
    fn default() -> Self {
        Names {
            index: HashMap::default(),
            names: Vec::new(),
            first: NO_NAME,
        }
    }
}

impl Names {
    /// Return the index of the name `bytes`, caching a new object for it.
    /// A full cache is emptied first, which resets `prev`.
    fn find(&mut self, py: Python<'_>, kind: Kind, bytes: &[u8], prev: &mut u32) -> PyResult<u32> {
        if let Some(&found) = self.index.get(bytes) {
            return Ok(found);
        }
        let object = kind.value(py, bytes)?.unbind();
        if self.names.len() >= MAX_NAMES {
            // Releasing exact str and bytes objects runs no Python code, so
            // it cannot re-enter the cache.
            self.index.clear();
            self.names.clear();
            self.first = NO_NAME;
            *prev = NO_NAME;
        }
        let found = self.names.len() as u32;
        self.index.insert(bytes.into(), found);
        self.names.push(Name {
            bytes: bytes.into(),
            object,
            next: NO_NAME,
        });
        Ok(found)
    }
}

thread_local! {
    /// Record-name objects by their bytes, for text and then bytes names.
    static NAMES: RefCell<[Names; 2]> = RefCell::default();
}

/// Record-name objects of one kind, from the cache when it is available.
struct Namer<'a> {
    kind: Kind,
    cache: Option<RefMut<'a, [Names; 2]>>,
    /// The previous cached name in this block, or `NO_NAME` at its start.
    prev: u32,
}

impl Namer<'_> {
    /// Return the object for the next record name of a block, read in order.
    fn name<'py>(&mut self, py: Python<'py>, bytes: &[u8]) -> PyResult<Bound<'py, PyAny>> {
        let Some(cache) = self.cache.as_deref_mut() else {
            return self.kind.value(py, bytes);
        };
        if bytes.len() > MAX_CACHED_NAME {
            return self.kind.value(py, bytes);
        }
        let names = &mut cache[self.kind as usize];
        let predicted = match self.prev {
            NO_NAME => names.first,
            prev => names.names[prev as usize].next,
        };
        if let Some(name) = names.names.get(predicted as usize)
            && *name.bytes == *bytes
        {
            self.prev = predicted;
            return Ok(name.object.bind(py).clone());
        }
        let found = names.find(py, self.kind, bytes, &mut self.prev)?;
        match self.prev {
            NO_NAME => names.first = found,
            prev => names.names[prev as usize].next = found,
        }
        self.prev = found;
        Ok(names.names[found as usize].object.bind(py).clone())
    }

    /// Return the object for a record name read out of block order.
    fn lookup<'py>(&mut self, py: Python<'py>, bytes: &[u8]) -> PyResult<Bound<'py, PyAny>> {
        let Some(cache) = self.cache.as_deref_mut() else {
            return self.kind.value(py, bytes);
        };
        if bytes.len() > MAX_CACHED_NAME {
            return self.kind.value(py, bytes);
        }
        let names = &mut cache[self.kind as usize];
        let mut unlinked = NO_NAME;
        let found = names.find(py, self.kind, bytes, &mut unlinked)?;
        Ok(names.names[found as usize].object.bind(py).clone())
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
            Ok(Some(record)) => block
                .kind
                .with_names(|names| names.lookup(py, record.name))
                .and_then(|name| PyTuple::new(py, [name, block.kind.value(py, &record.value)?]))
                .map(Some),
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
///
/// `SDBlock[str]` blocks come from text input and `SDBlock[bytes]` blocks
/// from binary input; names, values, and lines have the block's type.
#[pyclass(name = "SDBlock", module = "sendoff.native", eq, skip_from_py_object)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySDBlock {
    block: Block,
    kind: Kind,
}

impl PySDBlock {
    fn record<'py>(
        &self,
        py: Python<'py>,
        names: &mut Namer<'_>,
        record: Record<'_>,
    ) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            [
                names.name(py, record.name)?,
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
    /// File objects are read in 1 MiB chunks with `read1`, when available,
    /// or `read`, so the handle's position runs ahead of the blocks
    /// produced. Line iterables are consumed only up to each block's `$$$$`
    /// line, and the iterator does not close them.
    ///
    /// Args:
    ///     lines: homogeneous str or bytes lines separated by the `$$$$`
    ///         delimiter, or a text or binary file object
    ///
    /// Returns:
    ///     An iterator that produces owned SDBlocks as input is consumed
    ///
    /// Raises:
    ///     TypeError: for a str or bytes argument; use `from_sdf` for those
    #[staticmethod]
    #[pyo3(signature = (
        lines: "Iterable[str] | Iterable[bytes] | _SupportsRead[str] | _SupportsRead[bytes]",
    ) -> "_SDBlockIterator[Any]")]
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
    ///
    /// Raises:
    ///     ValueError: if `sdf` holds more than one block
    ///     TypeError: if `sdf` is not str or bytes, or a record is not its type
    #[new]
    #[pyo3(signature = (
        sdf: "AnyStr",
        records: "Mapping[AnyStr, AnyStr] | Iterable[tuple[AnyStr, AnyStr]] | None" = None,
    ) -> "SDBlock[AnyStr]")]
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
    ///
    /// Raises:
    ///     TypeError: if `sdf` is not str or bytes
    #[staticmethod]
    #[pyo3(signature = (sdf: "str | bytes") -> "_SDBlockIterator[Any]")]
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

    /// Return the molecule title in the input type.
    ///
    /// Titles, record names, and record values are stripped as `str.strip()`
    /// does for str blocks and as `bytes.strip()` does for bytes blocks.
    #[getter]
    fn title<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind.value(py, self.block.title(self.kind.strip()))
    }

    /// Replace the title line.
    ///
    /// Raises:
    ///     ValueError: if the title contains a line break
    ///     TypeError: if the title is not the block's type
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
    ///
    /// For example, `block.molblock = Chem.MolToMolBlock(mol)`. Text
    /// molblocks have CRLF and CR line endings translated to LF.
    ///
    /// Raises:
    ///     ValueError: if the molblock has no `M  END` line, has more than
    ///         whitespace after it, or has a line starting with `$$$$`
    ///     TypeError: if the molblock is not the block's type
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

    /// Return the CTAB after the title, through `M  END`, as one value.
    ///
    /// Every line ends in `\n`, so `needle in block.mdl` searches the CTAB,
    /// `block.mdl.splitlines()` gives its lines, and
    /// `block.title + "\n" + block.mdl` is `block.molblock`.
    #[getter]
    fn mdl<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind.value(py, self.block.mdl())
    }

    /// Return all SD data records.
    ///
    /// Returns:
    ///     Tuples of record name and value in the input type; multiline
    ///     values include their inner newlines
    #[pyo3(signature = () -> "list[tuple[AnyStr, AnyStr]]")]
    fn records<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty(py);
        let mut pos = 0;
        self.kind.with_names(|names| {
            while let Some(record) = self
                .block
                .next_record(&mut pos, self.kind.strip())
                .map_err(|error| PyValueError::new_err(error.to_string()))?
            {
                list.append(self.record(py, names, record)?)?;
            }
            Ok(list)
        })
    }

    /// Lazily iterate over SD data records in input order.
    ///
    /// Names and values are owned and have the block's type.
    #[pyo3(signature = () -> "_SDRecordIterator[AnyStr]")]
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
    ///
    /// Raises:
    ///     TypeError: if `name` is not the block's type
    #[pyo3(signature = (name: "AnyStr", default: "_T | None" = None) -> "AnyStr | _T | None")]
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

    /// Return the record names in input order, duplicates included.
    ///
    /// Equivalent to `[name for name, _ in block.records()]`, without
    /// building the values. Blocks need not share names, so collect them
    /// across a file with `Counter(n for b in blocks for n in b.names())`.
    #[pyo3(signature = () -> "list[AnyStr]")]
    fn names<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty(py);
        let strip = self.kind.strip();
        let mut pos = 0;
        self.kind.with_names(|names| {
            while let Some(name) = self
                .block
                .next_record_name(&mut pos, strip)
                .map_err(|error| PyValueError::new_err(error.to_string()))?
            {
                list.append(names.name(py, name)?)?;
            }
            Ok(list)
        })
    }

    /// The number of records, duplicates included: `len(block.records())`.
    fn __len__(&self) -> PyResult<usize> {
        let strip = self.kind.strip();
        let (mut pos, mut count) = (0, 0);
        while self
            .block
            .next_record_name(&mut pos, strip)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
            .is_some()
        {
            count += 1;
        }
        Ok(count)
    }

    /// A block is true even without records, as before it had a length.
    fn __bool__(&self) -> bool {
        true
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!("<SDBlock title={}>", self.title(py)?.repr()?))
    }

    /// Return an independent copy; a block holds no other objects, so
    /// `copy.copy` and `copy.deepcopy` are the same.
    #[pyo3(signature = () -> "SDBlock[AnyStr]")]
    fn __copy__(&self) -> Self {
        self.clone()
    }

    #[pyo3(signature = (memo: "dict[int, Any]", /) -> "SDBlock[AnyStr]")]
    fn __deepcopy__(&self, memo: &Bound<'_, PyAny>) -> Self {
        let _ = memo;
        self.clone()
    }

    /// Pickle the block's exact state rather than its SDF text, which
    /// does not always parse back to the same block.
    #[pyo3(signature = () -> "tuple[Any, ...]")]
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let py = slf.py();
        let this = slf.borrow();
        let (data, title_end, mdl_end) = this.block.parts();
        let state = (
            PyBytes::new(py, data),
            title_end,
            mdl_end,
            this.kind == Kind::Text,
        );
        (slf.get_type().getattr(intern!(py, "_from_state"))?, state).into_pyobject(py)
    }

    /// Rebuild a block pickled by `__reduce__`.
    ///
    /// Raises:
    ///     ValueError: if the state cannot come from a block
    #[classmethod]
    #[pyo3(signature = (
        data: "bytes", title_end: "int", mdl_end: "int", text: "bool",
    ) -> "SDBlock[Any]")]
    fn _from_state(
        _cls: &Bound<'_, PyType>,
        data: &[u8],
        title_end: usize,
        mdl_end: usize,
        text: bool,
    ) -> PyResult<Self> {
        let kind = if text { Kind::Text } else { Kind::Bytes };
        Block::from_parts(data.to_vec(), title_end, mdl_end)
            .filter(|block| !text || std::str::from_utf8(block.parts().0).is_ok())
            .map(|block| Self { block, kind })
            .ok_or_else(|| PyValueError::new_err("invalid SDBlock state"))
    }

    /// Whether any record is named `name`: `name in block.names()`.
    ///
    /// Like `get`, every record header is checked, so a malformed one raises
    /// even after a match.
    ///
    /// Raises:
    ///     TypeError: if `name` is not the block's type
    ///     ValueError: if a record header has no `<name>`
    fn __contains__(&self, name: &Bound<'_, PyAny>) -> PyResult<bool> {
        let name = self.arg_bytes(name, "record names")?;
        let strip = self.kind.strip();
        let (mut pos, mut seen) = (0, false);
        while let Some(found) = self
            .block
            .next_record_name(&mut pos, strip)
            .map_err(|error| PyValueError::new_err(error.to_string()))?
        {
            seen |= found == name;
        }
        Ok(seen)
    }

    /// Return the records as a dict, later duplicates winning.
    ///
    /// Equivalent to `dict(block.records())`.
    #[pyo3(signature = () -> "dict[AnyStr, AnyStr]")]
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        let strip = self.kind.strip();
        let mut pos = 0;
        self.kind.with_names(|names| {
            while let Some(record) = self
                .block
                .next_record(&mut pos, strip)
                .map_err(|error| PyValueError::new_err(error.to_string()))?
            {
                dict.set_item(
                    names.name(py, record.name)?,
                    self.kind.value(py, &record.value)?,
                )?;
            }
            Ok(dict)
        })
    }

    /// Return the block as SDF text ending with $$$$, in the input type.
    ///
    /// The inverse of `SDBlock(sdf)`, and what `write` writes.
    #[pyo3(signature = () -> "AnyStr")]
    fn to_sdf<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind
            .value(py, &self.block.serialized(self.kind.strip()))
    }

    /// Write the block, ending with $$$$, to a file-like handle.
    ///
    /// Every line ends with "\n", and the block is written in one `write`
    /// call.
    ///
    /// Args:
    ///     outh: a text or binary file-like handle matching the input type
    #[pyo3(signature = (outh: "_SupportsWrite[AnyStr]"))]
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
    ///
    /// Raises:
    ///     TypeError: if `record_name` or `value` is not the block's type
    #[pyo3(signature = (record_name: "AnyStr", value: "AnyStr"))]
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
    ///
    /// Raises:
    ///     TypeError: if `name` or `value` is not the block's type
    #[pyo3(signature = (name: "AnyStr", value: "AnyStr"))]
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
    ///
    /// Raises:
    ///     TypeError: if `name` is not the block's type
    #[pyo3(signature = (name: "AnyStr"))]
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
    ///
    /// Raises:
    ///     TypeError: if a name or value is not the block's type
    #[pyo3(signature = (records: "Mapping[AnyStr, AnyStr] | Iterable[tuple[AnyStr, AnyStr]]"))]
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
    /// Returns:
    ///     The CTable titled with this block's title, parsed from its MDL
    ///     lines without validation
    ///
    /// Raises:
    ///     NotImplementedError: for byte-backed SDBlocks
    fn ctable(&self, py: Python<'_>) -> PyResult<PyCTable> {
        self.require_text()?;
        let mut lines = vec![decode(py, self.block.title(Strip::Str))?.to_owned()];
        for line in block::lines(self.block.mdl()) {
            lines.push(decode(py, line)?.to_owned());
        }
        PyCTable::from_lines(lines)
    }

    /// Get the number of atoms declared by the counts line.
    ///
    /// The count may differ from the number of atom lines.
    fn num_atoms(&self) -> PyResult<i64> {
        let (_, atoms, _) = self.block.counts().map_err(error_to_py)?;
        Ok(atoms)
    }

    /// Get the number of bonds declared by the counts line.
    ///
    /// The count may differ from the number of bond lines.
    fn num_bonds(&self) -> PyResult<i64> {
        let (_, _, bonds) = self.block.counts().map_err(error_to_py)?;
        Ok(bonds)
    }

    /// Renumber atom and bond indices to be 1-indexed and in order.
    ///
    /// Uses the CTable renumbering and replaces this block's MDL lines.
    ///
    /// Raises:
    ///     NotImplementedError: for byte-backed SDBlocks
    fn renumber_indices(&mut self, py: Python<'_>) -> PyResult<()> {
        let mut ctable = self.ctable(py)?;
        ctable.renumber_indices_for_block(py)?;
        self.block
            .replace_mdl(ctable.lines_after_title().iter().map(String::as_bytes));
        Ok(())
    }
}

/// Where `write_sdf` sends serialized blocks, a chunk at a time.
enum Sink<'py> {
    File {
        file: File,
        filename: Py<PyAny>,
    },
    Handle {
        write: Bound<'py, PyAny>,
        kind: Option<Kind>,
    },
}

impl Sink<'_> {
    /// Check a block's kind before its serialization joins the chunk, and
    /// return whether to write the block on its own.
    ///
    /// CPython stores and encodes ASCII str much faster than wider str, and
    /// one non-ASCII character widens a whole str, so non-ASCII text goes
    /// to a file object a block at a time, keeping ASCII chunks narrow.
    fn accept(&mut self, py: Python<'_>, kind: Kind, serialized: &[u8]) -> PyResult<bool> {
        match self {
            // A file object's write validates text chunks as it decodes them.
            Sink::Handle { kind: expected, .. } => {
                set_kind(expected, kind).map_err(|_| {
                    PyTypeError::new_err("write_sdf does not accept mixed str and bytes blocks")
                })?;
                Ok(kind == Kind::Text && !serialized.is_ascii())
            }
            Sink::File { .. } if kind == Kind::Text => decode(py, serialized).map(|_| false),
            Sink::File { .. } => Ok(false),
        }
    }

    fn write(&mut self, py: Python<'_>, chunk: &[u8]) -> PyResult<()> {
        match self {
            Sink::File { file, filename } => py
                .detach(|| file.write_all(chunk))
                .map_err(|error| os_error(py, error, filename)),
            Sink::Handle { write, kind } => {
                write.call1((kind.unwrap_or(Kind::Bytes).value(py, chunk)?,))?;
                Ok(())
            }
        }
    }
}

/// Write SDBlocks to a path or a file object, and return how many.
///
/// Blocks are written as `SDBlock.write` would, but a MiB or so at a
/// time: to a path without holding the GIL, or to a file object with one
/// `write` call per chunk, as str for text blocks or bytes for byte-backed
/// ones. Text blocks must be valid UTF-8.
///
/// Args:
///     dest: a str, bytes, or os.PathLike path to create or truncate, or a
///         text or binary file object
///     blocks: an iterable of SDBlocks, all str or all bytes for a file
///         object
///
/// Returns:
///     The number of blocks written
#[pyfunction]
#[pyo3(signature = (
    dest: "str | bytes | os.PathLike[str] | os.PathLike[bytes] | _SupportsWrite[Any]",
    blocks: "Iterable[SDBlock[Any]]",
))]
pub(crate) fn write_sdf(dest: &Bound<'_, PyAny>, blocks: &Bound<'_, PyAny>) -> PyResult<usize> {
    let py = dest.py();
    let mut sink = match dest.getattr(intern!(py, "write")) {
        Ok(write) => Sink::Handle { write, kind: None },
        Err(_) => {
            let os = PyModule::import(py, "os")?;
            let filename = os.call_method1("fspath", (dest,))?;
            let path: PathBuf = os.call_method1("fsdecode", (&filename,))?.extract()?;
            let filename = filename.unbind();
            let file = py
                .detach(|| File::create(&path))
                .map_err(|error| os_error(py, error, &filename))?;
            Sink::File { file, filename }
        }
    };
    // Chunks for a file object become new Python objects, which stay
    // cheaper below malloc's mmap threshold.
    let limit = match sink {
        Sink::File { .. } => CHUNK,
        Sink::Handle { .. } => 64 << 10,
    };
    let mut chunk = Vec::with_capacity(2 * limit);
    let mut count = 0;
    for block in blocks.try_iter()? {
        let block = block?;
        let block = block.cast::<PySDBlock>()?.borrow();
        let start = chunk.len();
        block.block.serialize_into(block.kind.strip(), &mut chunk);
        let alone = sink.accept(py, block.kind, &chunk[start..])?;
        count += 1;
        if alone {
            if start > 0 {
                sink.write(py, &chunk[..start])?;
            }
            sink.write(py, &chunk[start..])?;
            chunk.clear();
        } else if chunk.len() >= limit {
            sink.write(py, &chunk)?;
            chunk.clear();
        }
    }
    if !chunk.is_empty() {
        sink.write(py, &chunk)?;
    }
    Ok(count)
}

/// Lazily parse SDBlocks from an SDF path.
///
/// The file is opened on the first `next()` and read in chunks without
/// holding the GIL. str paths give str values decoded as UTF-8 on access,
/// with universal newlines; bytes paths give bytes split on b"\n" only.
/// Use a bytes path for feeds that are not valid UTF-8.
///
/// Args:
///     path: a str, bytes, or os.PathLike path
///
/// Returns:
///     An iterator of SDBlocks; `close()` closes the file early
#[pyfunction]
#[pyo3(signature = (path: "str | bytes | os.PathLike[str] | os.PathLike[bytes]") -> "_SDBlockIterator[Any]")]
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
