//! The `SDBlock` Python class and its record iterator.

use std::borrow::Cow;
use std::ffi::c_int;
use std::sync::atomic::{AtomicBool, Ordering};

use pyo3::exceptions::{
    PyKeyError, PyRuntimeError, PyTypeError, PyUnicodeDecodeError, PyValueError,
};
use pyo3::ffi;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyGenericAlias, PyList, PyString, PyTuple, PyType};

use crate::core::block::{Block, Parts};
use crate::core::ctable::CTable;
use crate::core::ctable::error::CTableError;
use crate::core::lines::lines;
use crate::core::newline::{translated, universal_newlines};
use crate::core::record::{HeaderStyle, Record};
use crate::core::split::single_block;
use crate::python::ctable::PyCTable;
use crate::python::exports::Exported;
use crate::python::kind::{Kind, bytes_of, decode};
use crate::python::names::Namer;
use crate::python::read::{PySDBlockIterator, Source, reader};

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

/// The style a `header` argument names.
fn header_style(header: &str) -> PyResult<HeaderStyle> {
    match header {
        "plain" => Ok(HeaderStyle::Plain),
        "padded" => Ok(HeaderStyle::Padded),
        "copy" => Ok(HeaderStyle::Copy),
        _ => Err(PyValueError::new_err(format!(
            "header must be 'plain', 'padded', or 'copy', not '{header}'"
        ))),
    }
}

/// Lazily iterate over one block's SD data records.
#[pyclass(name = "_SDRecordIterator", module = "sendoff.native")]
pub(crate) struct PySDRecordIterator {
    block: Py<PySDBlock>,
    edits: u64,
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
        let block = unchanged(self.block.bind(py), self.edits, &mut self.finished)?;
        let mut records = block.block.records_at(self.pos);
        let result = match records.next() {
            Some(Ok(record)) => Namer::with(block.kind, |names| names.lookup(py, record.name))
                .and_then(|name| PyTuple::new(py, [name, block.kind.value(py, &record.value)?]))
                .map_err(|error| block.noted(py, error, || reading_record(record.name)))
                .map(Some),
            None => Ok(None),
            Some(Err(error)) => Err(error.into()),
        };
        self.pos = records.pos();
        if !matches!(result, Ok(Some(_))) {
            self.finished = true;
        }
        result
    }
}

/// Lazily iterate over one block's record names.
#[pyclass(name = "_SDNameIterator", module = "sendoff.native")]
pub(crate) struct PySDNameIterator {
    block: Py<PySDBlock>,
    edits: u64,
    pos: usize,
    finished: bool,
}

#[pymethods]
impl PySDNameIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        if self.finished {
            return Ok(None);
        }
        let block = unchanged(self.block.bind(py), self.edits, &mut self.finished)?;
        let mut headers = block.block.records_at(self.pos).headers();
        let result = match headers.next() {
            Some(Ok(header)) => {
                self.pos = header.end;
                Namer::with(block.kind, |names| names.lookup(py, header.name))
                    .map_err(|error| block.noted(py, error, || reading_record(header.name)))
                    .map(Some)
            }
            None => Ok(None),
            Some(Err(error)) => Err(error.into()),
        };
        if !matches!(result, Ok(Some(_))) {
            self.finished = true;
        }
        result
    }
}

/// Borrow an iterator's block, or end the iteration with RuntimeError if
/// its records changed since the iterator was made.
fn unchanged<'py>(
    block: &Bound<'py, PySDBlock>,
    edits: u64,
    finished: &mut bool,
) -> PyResult<PyRef<'py, PySDBlock>> {
    let block = block.try_borrow()?;
    if block.edits != edits {
        *finished = true;
        return Err(PyRuntimeError::new_err(
            "SDBlock records changed during iteration",
        ));
    }
    Ok(block)
}

/// Handle a single molecule block from an SD file.
///
/// `SDBlock[str]` blocks come from text input and `SDBlock[bytes]` blocks
/// from binary input; names, values, and lines have the block's type.
///
/// Record names are read and written as given, unchecked. The CTfile
/// format reserves `>`, and its field names start with a letter followed by
/// letters, digits, and underscores: no `-`, `.`, `<`, `>`, `=`, `%`, or
/// spaces. Keeping names valid is up to you; other toolkits may read
/// invalid ones differently.
///
/// A record's value is its lines up to a blank one, as read but without
/// their line endings, joined with "\n". Lines of only ASCII whitespace,
/// which `bytes.strip()` removes, are blank, in str blocks too, as in RDKit
/// and OEChem. Methods that read
/// records raise sendoff.MalformedRecordError, a ValueError, at a header
/// line with no `<name>`. Values are written unchecked too: a value with a
/// blank line or a `$$$$` line reads back as more records or blocks.
#[pyclass(
    name = "SDBlock",
    module = "sendoff.native",
    eq,
    mapping,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PySDBlock {
    pub(super) block: Exported,
    pub(super) kind: Kind,
    /// Counts edits to the SD data, so record iterators can tell they are
    /// stale without copying the records.
    edits: u64,
    /// Whether a text block's data is known to be UTF-8, as it is when read
    /// from str or already decoded whole. Edits keep it: they replace whole
    /// lines with str arguments.
    utf8: KnownUtf8,
}

/// An atomic flag, since writes only borrow blocks.
#[derive(Default)]
struct KnownUtf8(AtomicBool);

impl Clone for KnownUtf8 {
    fn clone(&self) -> Self {
        Self(AtomicBool::new(self.0.load(Ordering::Relaxed)))
    }
}

impl PartialEq for PySDBlock {
    fn eq(&self, other: &Self) -> bool {
        self.kind == other.kind && self.block == other.block
    }
}

impl PySDBlock {
    /// Wrap `block`; `utf8` says its data is known to be UTF-8.
    pub(super) fn wrap(block: Block, kind: Kind, utf8: bool) -> Self {
        Self {
            block: Exported::new(block),
            kind,
            edits: 0,
            utf8: KnownUtf8(AtomicBool::new(utf8)),
        }
    }

    /// Check that a text block is UTF-8 unless that is already known.
    pub(super) fn check_utf8(&self, py: Python<'_>) -> PyResult<()> {
        if self.kind == Kind::Text && !self.utf8.0.load(Ordering::Relaxed) {
            decode(py, self.block.serialized())?;
            self.utf8.0.store(true, Ordering::Relaxed);
        }
        Ok(())
    }

    fn edited(&mut self) {
        self.edits = self.edits.wrapping_add(1);
    }

    /// The value of the last record named `name`, if any.
    fn last_value<'py>(&self, name: &Bound<'py, PyAny>) -> PyResult<Option<Bound<'py, PyAny>>> {
        let py = name.py();
        let name = self.arg_bytes(name, "record names")?;
        let value = match self.block.index()? {
            Some(index) => index.last(self.block.metadata(), name),
            None => {
                let mut found = None;
                for header in self.block.headers() {
                    let header = header?;
                    if header.name == name {
                        found = Some(header.start);
                    }
                }
                let Some(start) = found else {
                    return Ok(None);
                };
                let record = self.block.records_at(start).next();
                let record = record.and_then(Result::ok).expect("a record starts here");
                Some(record.value)
            }
        };
        let Some(value) = value else {
            return Ok(None);
        };
        self.kind
            .value(py, &value)
            .map(Some)
            .map_err(|error| self.noted(py, error, || reading_record(name)))
    }

    /// Set the record `name` as `set_record` does.
    fn put(
        &mut self,
        name: &Bound<'_, PyAny>,
        value: &Bound<'_, PyAny>,
        style: HeaderStyle,
    ) -> PyResult<()> {
        let name = self.arg_bytes(name, "record names")?;
        let value = self.arg_bytes(value, "values")?;
        self.block.get_mut()?.set_record(name, value, style)?;
        self.edited();
        Ok(())
    }

    /// The name and value objects of a record read in block order.
    ///
    /// Inlined into records() and to_dict(): called out of line once per
    /// record, it cost those 2-3%.
    #[inline(always)]
    fn record<'py>(
        &self,
        py: Python<'py>,
        names: &mut Namer<'_>,
        record: &Record<'_>,
    ) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
        let noted = |error| self.noted(py, error, || reading_record(record.name));
        Ok((
            names.name(py, record.name).map_err(noted)?,
            self.kind.value(py, &record.value).map_err(noted)?,
        ))
    }

    /// Parse the title and MDL lines as a CTable, decoding a bytes block's
    /// as Latin-1, so every byte round-trips.
    fn parse_ctable(&self, py: Python<'_>) -> PyResult<CTable> {
        let title = self.block.title();
        let decoded = match self.kind {
            Kind::Text => {
                let noted = |error| self.noted(py, error, || "Reading the CTAB of".into());
                let mut decoded = vec![decode(py, title).map_err(noted)?.to_owned()];
                for line in lines(self.block.mdl()) {
                    decoded.push(decode(py, line).map_err(noted)?.to_owned());
                }
                decoded
            }
            Kind::Bytes => std::iter::once(title)
                .chain(lines(self.block.mdl()))
                .map(|line| line.iter().map(|&byte| char::from(byte)).collect())
                .collect(),
        };
        Ok(CTable::parse(decoded)?)
    }

    /// Apply `edit` to the block's CTable, then replace the MDL lines with
    /// the edited table's.
    fn edit_ctable(
        &mut self,
        py: Python<'_>,
        edit: impl FnOnce(&mut CTable) -> Result<(), CTableError>,
    ) -> PyResult<()> {
        let mut ctable = self.parse_ctable(py)?;
        edit(&mut ctable)?;
        let lines = &ctable.lines[1..];
        match self.kind {
            Kind::Text => self
                .block
                .get_mut()?
                .replace_mdl(lines.iter().map(String::as_bytes)),
            Kind::Bytes => {
                // Edits write ASCII, so every char is a Latin-1 byte.
                let encoded: Vec<Vec<u8>> = lines
                    .iter()
                    .map(|line| line.chars().map(|c| c as u8).collect())
                    .collect();
                self.block
                    .get_mut()?
                    .replace_mdl(encoded.iter().map(Vec::as_slice));
            }
        }
        Ok(())
    }

    /// Add to a decode error a note naming this block, what of it `action`
    /// was doing, and how to read other encodings.
    pub(super) fn noted(
        &self,
        py: Python<'_>,
        error: PyErr,
        action: impl FnOnce() -> String,
    ) -> PyErr {
        if error.is_instance_of::<PyUnicodeDecodeError>(py) {
            let title = String::from_utf8_lossy(self.block.title());
            let note = format!("{} the SDBlock titled '{title}'. {DECODE_HINT}", action());
            if let Err(error) = error.add_note(py, note) {
                return error;
            }
        }
        error
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

/// How to read SDF that is not UTF-8, for notes on decode errors.
const DECODE_HINT: &str = "parse_sdf(os.fsencode(path)) reads such files as bytes, and \
    parse_sdf(open(path, encoding=...)) decodes another encoding.";

fn reading_record(name: &[u8]) -> String {
    format!("Reading record '{}' of", String::from_utf8_lossy(name))
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
    /// A file object, which has a `read1` or `read` method, is read in
    /// 128 KiB chunks with `read1`, when available, or `read`, so the
    /// handle's position runs ahead of the blocks produced. Any other
    /// iterable, such as a list or a generator, gives whole lines: each item
    /// ends a line and may hold several. An item that does not end in "\n"
    /// gets its block's title line ending ("\n" after a "\r"), so an empty
    /// item is an empty line; pass a stream's arbitrary chunks as a file
    /// object instead. A block is produced once an item starts with `$$$$`
    /// or the iterable ends, so line-by-line input is consumed only up to
    /// each block's `$$$$` line; the iterator does not close iterables.
    ///
    /// Lines end at "\n", and lines keep their endings, so CRLF input
    /// writes back as CRLF. `from_lines(f)` gives the blocks
    /// `from_lines(list(f))` does, except that a file opened with
    /// `newline=""` also ends its lines at a lone "\r", which chunks do not:
    /// open such files in default text mode, which translates "\r\n" and
    /// "\r" to "\n" as `parse_sdf` and `from_sdf` do for str input.
    ///
    /// Args:
    ///     lines: homogeneous str or bytes lines separated by the `$$$$`
    ///         delimiter, or a text or binary file object
    ///
    /// Returns:
    ///     An iterator that produces owned SDBlocks as input is consumed
    ///
    /// Raises:
    ///     TypeError: for a str or bytes argument (use `from_sdf` for those)
    ///         or a non-iterable; and, as blocks are produced, for items or
    ///         reads that are not str or bytes or that mix the two
    #[staticmethod]
    #[pyo3(signature = (
        lines: "Iterable[str] | Iterable[bytes] | _SupportsRead[str] | _SupportsRead[bytes]",
    ) -> "_SDBlockIterator[Any]")]
    fn from_lines(lines: &Bound<'_, PyAny>) -> PyResult<PySDBlockIterator> {
        if lines.is_instance_of::<PyString>() || lines.is_instance_of::<PyBytes>() {
            return Err(PyTypeError::new_err(
                "SDBlock.from_lines takes lines or a file object; \
                 use SDBlock.from_sdf to parse SDF str or bytes",
            ));
        }
        if let Some(read) = reader(lines) {
            return Ok(PySDBlockIterator::new(Source::Reader(read)));
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
    ///     header: how to write the records' header lines, as in
    ///         `append_record`
    ///
    /// Raises:
    ///     ValueError: if `sdf` holds more than one block, or for another
    ///         `header`
    ///     TypeError: if `sdf` is not str or bytes, or a record is not its type
    #[new]
    #[pyo3(signature = (
        sdf: "AnyStr",
        records: "Mapping[AnyStr, AnyStr] | Iterable[tuple[AnyStr, AnyStr]] | None" = None,
        *,
        header: "Literal['plain', 'padded', 'copy']" = "plain",
    ) -> "SDBlock[AnyStr]")]
    fn new(
        sdf: &Bound<'_, PyAny>,
        records: Option<&Bound<'_, PyAny>>,
        header: &str,
    ) -> PyResult<Self> {
        let style = header_style(header)?;
        let Some((kind, bytes)) = bytes_of(sdf)? else {
            return Err(PyTypeError::new_err(
                "SDBlock() argument must be str or bytes",
            ));
        };
        let data = match kind {
            Kind::Text => translated(bytes),
            Kind::Bytes => bytes.to_vec(),
        };
        let block = single_block(data).ok_or_else(|| {
            PyValueError::new_err(
                "SDBlock() takes one block; use SDBlock.from_sdf to parse several",
            )
        })?;
        let mut slf = Self::wrap(block, kind, true);
        if let Some(records) = records {
            for (name, value) in record_pairs(records)? {
                let name = slf.arg_bytes(&name, "record names")?;
                let value = slf.arg_bytes(&value, "values")?;
                slf.block.get_mut()?.append_record(name, value, style);
            }
        }
        Ok(slf)
    }

    /// Lazily parse SDBlocks from SDF str or bytes held in memory.
    ///
    /// Equivalent to `from_lines(io.StringIO(sdf, newline=None))` or
    /// `from_lines(io.BytesIO(sdf))`, without splitting lines: str input
    /// has universal newlines, translated to "\n", and bytes are split on
    /// b"\n" only.
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

    /// Return the title line without trailing whitespace, in the input type.
    ///
    /// Only ASCII whitespace is stripped, as `bytes.rstrip()` does, in str
    /// blocks too; `molblock` and writes keep the whole title line.
    #[getter]
    fn title<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind
            .value(py, self.block.title())
            .map_err(|error| self.noted(py, error, || "Reading the title of".into()))
    }

    /// Replace the title line.
    ///
    /// Raises:
    ///     ValueError: if the title contains a line break or starts with
    ///         `$$$$`, which would end the block
    ///     TypeError: if the title is not the block's type
    #[setter(title)]
    fn set_title(&mut self, title: &Bound<'_, PyAny>) -> PyResult<()> {
        let title = self.arg_bytes(title, "titles")?;
        self.block.get_mut()?.set_title(title).map_err(PyErr::from)
    }

    /// Return the title line and CTAB through `M  END` as one value.
    ///
    /// Suitable for RDKit's `Chem.MolFromMolBlock` or OpenEye's
    /// `OEReadMolFromBytes`.
    #[getter]
    fn molblock<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind
            .value(py, self.block.molblock())
            .map_err(|error| self.noted(py, error, || "Reading the molblock of".into()))
    }

    /// Replace the title line and CTAB, keeping the SD data.
    ///
    /// For example, `block.molblock = Chem.MolToMolBlock(mol)`. Text
    /// molblocks have CRLF and CR line endings translated to LF, and a
    /// molblock without a final newline gets its title line's ending. So
    /// `block.molblock = block.molblock` translates a str block's title line
    /// and CTAB to LF and changes nothing else.
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
            .get_mut()?
            .set_molblock(&molblock)
            .map_err(PyErr::from)
    }

    /// Return the CTAB after the title, through `M  END`, as one value.
    ///
    /// Every line ends in a newline, so `needle in block.mdl` searches the
    /// CTAB and `block.mdl.splitlines()` gives its lines. `block.molblock`
    /// is `block.title`, the title line's ending, then `block.mdl`.
    #[getter]
    fn mdl<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.kind
            .value(py, self.block.mdl())
            .map_err(|error| self.noted(py, error, || "Reading the MDL of".into()))
    }

    /// Return all SD data records.
    ///
    /// Returns:
    ///     Tuples of record name and value in the input type; multiline
    ///     values include their inner newlines
    #[pyo3(signature = () -> "list[tuple[AnyStr, AnyStr]]")]
    fn records<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let list = PyList::empty(py);
        Namer::with(self.kind, |names| {
            for record in self.block.records() {
                list.append(self.record(py, names, &record?)?)?;
            }
            Ok(list)
        })
    }

    /// Lazily iterate over SD data records in input order.
    ///
    /// Names and values are owned and have the block's type. Like a dict,
    /// the iterator raises RuntimeError if the records change; to edit them
    /// in a loop, loop over `records()`.
    #[pyo3(signature = () -> "_SDRecordIterator[AnyStr]")]
    fn iter_records(slf: Bound<'_, Self>) -> PySDRecordIterator {
        let edits = slf.borrow().edits;
        PySDRecordIterator {
            block: slf.unbind(),
            edits,
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
        default: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Option<Bound<'py, PyAny>>> {
        match self.last_value(name)? {
            Some(value) => Ok(Some(value)),
            None => Ok(default.cloned()),
        }
    }

    /// Return the value of the last record named `name`, like a dict.
    ///
    /// Raises:
    ///     KeyError: if no record has that name
    ///     TypeError: if `name` is not the block's type
    fn __getitem__<'py>(&self, name: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        self.last_value(name)?
            .ok_or_else(|| PyKeyError::new_err(name.clone().unbind()))
    }

    /// Set a record as `set_record(name, value)` does.
    ///
    /// Raises:
    ///     TypeError: if `name` or `value` is not the block's type
    fn __setitem__(&mut self, name: &Bound<'_, PyAny>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        self.put(name, value, HeaderStyle::Plain)
    }

    /// Remove every record named `name`, as `remove_record(name)` does.
    ///
    /// Raises:
    ///     KeyError: if no record has that name
    ///     TypeError: if `name` is not the block's type
    fn __delitem__(&mut self, name: &Bound<'_, PyAny>) -> PyResult<()> {
        if self.remove_record(name)? == 0 {
            return Err(PyKeyError::new_err(name.clone().unbind()));
        }
        Ok(())
    }

    /// Lazily iterate over the record names in input order, duplicates
    /// included, as `names()` lists them. Like a dict, the iterator raises
    /// RuntimeError if the records change.
    fn __iter__(slf: Bound<'_, Self>) -> PySDNameIterator {
        let edits = slf.borrow().edits;
        PySDNameIterator {
            block: slf.unbind(),
            edits,
            pos: 0,
            finished: false,
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
        Namer::with(self.kind, |names| {
            for header in self.block.headers() {
                let header = header?;
                let name = names.name(py, header.name);
                list.append(
                    name.map_err(|error| self.noted(py, error, || reading_record(header.name)))?,
                )?;
            }
            Ok(list)
        })
    }

    /// The number of records, duplicates included: `len(block.records())`.
    fn __len__(&self) -> PyResult<usize> {
        if let Some(index) = self.block.index()? {
            return Ok(index.len());
        }
        let mut count = 0;
        for header in self.block.headers() {
            header?;
            count += 1;
        }
        Ok(count)
    }

    /// A block is true even without records, as before it had a length.
    fn __bool__(&self) -> bool {
        true
    }

    /// Show the title; a str block's undecodable title bytes show as U+FFFD
    /// instead of raising.
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let title = self.block.title();
        let title = match self.kind {
            Kind::Text => PyString::new(py, &String::from_utf8_lossy(title)).into_any(),
            Kind::Bytes => PyBytes::new(py, title).into_any(),
        };
        Ok(format!("<SDBlock title={}>", title.repr()?))
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
        let this = slf.try_borrow()?;
        let Parts {
            data,
            title_end,
            mdl_end,
        } = this.block.parts();
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
        let kind = Kind::of(text);
        // Text blocks from str paths decode on access, so their state need
        // not be UTF-8.
        Block::from_parts(data.to_vec(), title_end, mdl_end)
            .map(|block| Self::wrap(block, kind, false))
            .ok_or_else(|| PyValueError::new_err("invalid SDBlock state"))
    }

    /// Whether any record is named `name`: `name in block.names()`.
    ///
    /// Like `get`, every record header is checked, so a malformed one raises
    /// even after a match.
    ///
    /// Raises:
    ///     TypeError: if `name` is not the block's type
    ///     MalformedRecordError: if a record header has no `<name>`
    fn __contains__(&self, name: &Bound<'_, PyAny>) -> PyResult<bool> {
        let name = self.arg_bytes(name, "record names")?;
        if let Some(index) = self.block.index()? {
            return Ok(index.contains(self.block.metadata(), name));
        }
        let mut seen = false;
        for header in self.block.headers() {
            seen |= header?.name == name;
        }
        Ok(seen)
    }

    /// Return the records as a dict, later duplicates winning.
    ///
    /// Equivalent to `dict(block.records())`.
    #[pyo3(signature = () -> "dict[AnyStr, AnyStr]")]
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        Namer::with(self.kind, |names| {
            for record in self.block.records() {
                let (name, value) = self.record(py, names, &record?)?;
                dict.set_item(name, value)?;
            }
            Ok(dict)
        })
    }

    /// A str block's SDF, ending with $$$$: the inverse of `SDBlock(sdf)`.
    ///
    /// `bytes(block)` is a bytes block's SDF; `str` gives its repr, as
    /// `str(b"")` does. Lines end as they were read, so blocks from a bytes
    /// path keep "\r\n"; str paths read universal newlines, so their lines
    /// end in "\n". Lines the block adds or ends (new records, a missing
    /// $$$$, a last line without "\n") end like the title line.
    ///
    /// Raises:
    ///     UnicodeDecodeError: for a block from a str path that is not UTF-8
    fn __str__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyString>> {
        if self.kind == Kind::Bytes {
            return Ok(PyString::new(py, &self.__repr__(py)?));
        }
        let text = PyString::from_bytes(py, self.block.serialized())
            .map_err(|error| self.noted(py, error, || "Writing".into()))?;
        self.utf8.0.store(true, Ordering::Relaxed);
        Ok(text)
    }

    /// Export a bytes block's SDF without copying.
    ///
    /// Edits raise BufferError until every export is released.
    ///
    /// Raises:
    ///     TypeError: for a str block
    ///     BufferError: for a writable buffer
    #[allow(unsafe_code)]
    unsafe fn __getbuffer__(
        slf: &Bound<'_, Self>,
        view: *mut ffi::Py_buffer,
        flags: c_int,
    ) -> PyResult<()> {
        let this = slf.try_borrow()?;
        if this.kind == Kind::Text {
            return Err(PyTypeError::new_err(
                "a str SDBlock is not bytes-like; use str(block)",
            ));
        }
        let data = this.block.serialized();
        // SAFETY: Python passes a view to fill, which takes a reference to
        // `slf`, keeping the block alive, and edits, which could move or
        // free `data`, raise until `__releasebuffer__` releases the export.
        let filled = unsafe {
            ffi::PyBuffer_FillInfo(
                view,
                slf.as_ptr(),
                data.as_ptr().cast_mut().cast(),
                data.len() as ffi::Py_ssize_t,
                1,
                flags,
            )
        };
        if filled == -1 {
            return Err(PyErr::fetch(slf.py()));
        }
        // SAFETY: the view was filled, so it is valid to write.
        unsafe { (*view).internal = this.block.export() };
        Ok(())
    }

    #[allow(unsafe_code)]
    unsafe fn __releasebuffer__(_slf: Bound<'_, Self>, view: *mut ffi::Py_buffer) {
        // SAFETY: Python releases each view `__getbuffer__` filled once.
        unsafe { Exported::release((*view).internal) };
    }

    /// Append a field and value of the same type as the input block.
    ///
    /// Args:
    ///     name: field name, written as given; see SDBlock
    ///     value: field value, with no terminating newline
    ///     header: how to write the header line: "plain" `> <C>`; "padded"
    ///         `>  <C>`, as RDKit and Open Babel write; or "copy" the text
    ///         around the name on the first named header line, as
    ///         `>  <C>  (1)` after RDKit's `>  <A>  (1)`, or `> <C>` without
    ///         one. "copy" reads that header first, so it is a bit slower.
    ///
    /// Raises:
    ///     TypeError: if `name` or `value` is not the block's type
    ///     ValueError: for another `header`
    #[pyo3(signature = (name: "AnyStr", value: "AnyStr", *, header: "Literal['plain', 'padded', 'copy']" = "plain"))]
    fn append_record(
        &mut self,
        name: &Bound<'_, PyAny>,
        value: &Bound<'_, PyAny>,
        header: &str,
    ) -> PyResult<()> {
        let style = header_style(header)?;
        let name = self.arg_bytes(name, "record names")?;
        let value = self.arg_bytes(value, "values")?;
        self.block.get_mut()?.append_record(name, value, style);
        self.edited();
        Ok(())
    }

    /// Set a record, keeping the header line and position of the first with
    /// that name.
    ///
    /// Later records with the same name are removed; if there is none, the
    /// record is appended.
    ///
    /// Args:
    ///     name: record name in the block's type, written as given; see
    ///         SDBlock
    ///     value: record value in the block's type
    ///     header: how to write an appended record's header line, as in
    ///         `append_record`
    ///
    /// Raises:
    ///     TypeError: if `name` or `value` is not the block's type
    ///     ValueError: for another `header`
    #[pyo3(signature = (name: "AnyStr", value: "AnyStr", *, header: "Literal['plain', 'padded', 'copy']" = "plain"))]
    fn set_record(
        &mut self,
        name: &Bound<'_, PyAny>,
        value: &Bound<'_, PyAny>,
        header: &str,
    ) -> PyResult<()> {
        self.put(name, value, header_style(header)?)
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
        let removed = self.block.get_mut()?.remove_records(name)?;
        if removed > 0 {
            self.edited();
        }
        Ok(removed)
    }

    /// Replace all SD data with the given records.
    ///
    /// The k-th record with a name keeps the k-th old header line with that
    /// name.
    ///
    /// Args:
    ///     records: a mapping, or an iterable of (name, value) tuples, in
    ///         the block's type; names are written as given, see SDBlock
    ///     header: how to write other header lines, as in `append_record`;
    ///         "copy" copies from the old SD data
    ///
    /// Raises:
    ///     TypeError: if a name or value is not the block's type
    ///     ValueError: for another `header`
    #[pyo3(signature = (
        records: "Mapping[AnyStr, AnyStr] | Iterable[tuple[AnyStr, AnyStr]]",
        *,
        header: "Literal['plain', 'padded', 'copy']" = "plain",
    ))]
    fn replace_records(
        slf: &Bound<'_, Self>,
        records: &Bound<'_, PyAny>,
        header: &str,
    ) -> PyResult<()> {
        let style = header_style(header)?;
        // Collect the records before borrowing the block to change it, as
        // they may come from the block itself.
        let pairs = record_pairs(records)?;
        let mut this = slf.try_borrow_mut()?;
        let records = pairs
            .iter()
            .map(|(name, value)| {
                Ok((
                    this.arg_bytes(name, "record names")?,
                    this.arg_bytes(value, "values")?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        this.block.get_mut()?.replace_records(records, style);
        this.edited();
        Ok(())
    }

    /// Parse out the underlying CTable object.
    ///
    /// A bytes block's lines are decoded as Latin-1, so
    /// `line.encode("latin-1")` gives the bytes back.
    ///
    /// Returns:
    ///     The CTable titled with this block's title, parsed from its MDL
    ///     lines without validation
    ///
    /// Raises:
    ///     MalformedCTableError: if the table is shorter than its header,
    ///         or its format or counts cannot be read
    fn ctable(&self, py: Python<'_>) -> PyResult<PyCTable> {
        Ok(self.parse_ctable(py)?.into())
    }

    /// The number of atoms declared by the counts line.
    ///
    /// The count may differ from the number of atom lines. Only the
    /// header is read, so unlike ctable a block that doesn't decode is not
    /// an error.
    ///
    /// Raises:
    ///     MalformedCTableError: if the format or counts cannot be read
    #[getter]
    fn num_atoms(&self) -> PyResult<i64> {
        Ok(self.block.counts()?.atoms)
    }

    /// The number of bonds declared by the counts line.
    ///
    /// The count may differ from the number of bond lines. Only the
    /// header is read, so unlike ctable a block that doesn't decode is not
    /// an error.
    ///
    /// Raises:
    ///     MalformedCTableError: if the format or counts cannot be read
    #[getter]
    fn num_bonds(&self) -> PyResult<i64> {
        Ok(self.block.counts()?.bonds)
    }

    /// Renumber atom and bond indices to be 1-indexed and in order.
    ///
    /// Uses the CTable renumbering and replaces this block's MDL lines. The
    /// counts line stays as it is; see fix_counts.
    ///
    /// Raises:
    ///     MalformedCTableError: as ctable and CTable.renumber_indices raise
    ///     IndicesDuplicateError: as CTable.renumber_indices raises
    ///     NotImplementedError: as CTable.renumber_indices raises
    fn renumber_indices(&mut self, py: Python<'_>) -> PyResult<()> {
        self.edit_ctable(py, CTable::renumber_indices)
    }

    /// Set the counts line to the number of atom and bond lines.
    ///
    /// Uses CTable.fix_counts and replaces this block's MDL lines.
    ///
    /// Raises:
    ///     MalformedCTableError: as ctable and CTable.fix_counts raise
    ///     NotImplementedError: as CTable.fix_counts raises
    fn fix_counts(&mut self, py: Python<'_>) -> PyResult<()> {
        self.edit_ctable(py, CTable::fix_counts)
    }
}
