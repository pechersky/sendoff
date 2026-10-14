//! Writing SDBlocks to paths and file objects a chunk at a time.

use std::fs::File;
use std::io::Write;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use pyo3::{PyTraverseError, PyVisit, intern};

use crate::python::kind::{Kind, same_kind};
use crate::python::os::{fs_path, os_error};
use crate::python::sdblock::PySDBlock;

/// Bytes buffered before each write to a path. 64 KiB to 1 MiB timed the
/// same, and 4 MiB was slower.
const WRITE_CHUNK: usize = 256 << 10;

/// Where serialized blocks go, a chunk at a time.
enum Sink {
    File {
        file: File,
        filename: Py<PyAny>,
    },
    Handle {
        write: Py<PyAny>,
        kind: Option<Kind>,
    },
}

impl Sink {
    /// Check a block's kind before its serialization joins the chunk, and
    /// return whether to write the block on its own.
    ///
    /// CPython stores and encodes ASCII str much faster than wider str, and
    /// one non-ASCII character widens a whole str, so non-ASCII text goes
    /// to a file object a block at a time, keeping ASCII chunks narrow.
    fn accept(&mut self, py: Python<'_>, block: &PySDBlock) -> PyResult<bool> {
        let kind = block.kind;
        match self {
            // A file object's write validates text chunks as it decodes them.
            Sink::Handle { kind: expected, .. } => {
                same_kind(
                    expected,
                    kind,
                    "cannot write mixed str and bytes blocks to a file object",
                )?;
                Ok(kind == Kind::Text && !block.block.serialized().is_ascii())
            }
            Sink::File { .. } => block.check_utf8(py).map(|()| false),
        }
    }

    fn write(&mut self, py: Python<'_>, chunk: &[u8]) -> PyResult<()> {
        match self {
            Sink::File { file, filename } => py
                .detach(|| file.write_all(chunk))
                .map_err(|error| os_error(py, error, filename)),
            Sink::Handle { write, kind } => {
                write.call1(py, (kind.unwrap_or(Kind::Bytes).value(py, chunk)?,))?;
                Ok(())
            }
        }
    }
}

/// Serialized blocks waiting to be written to a sink in one call.
struct Writer {
    sink: Sink,
    /// Only blocks not yet passed to a write, so after an error it holds
    /// the blocks before the failing one.
    chunk: Vec<u8>,
    limit: usize,
    count: usize,
}

impl Writer {
    /// Create or truncate a path, or take a file object's `write`.
    fn open(dest: &Bound<'_, PyAny>) -> PyResult<Self> {
        let py = dest.py();
        let sink = match dest.getattr(intern!(py, "write")) {
            Ok(write) => Sink::Handle {
                write: write.unbind(),
                kind: None,
            },
            Err(_) => {
                let (filename, path) = fs_path(dest)?;
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
            Sink::File { .. } => WRITE_CHUNK,
            Sink::Handle { .. } => 64 << 10,
        };
        Ok(Writer {
            sink,
            chunk: Vec::with_capacity(2 * limit),
            limit,
            count: 0,
        })
    }

    /// Add a block, writing the chunk once it is full.
    fn push(&mut self, py: Python<'_>, block: &PySDBlock) -> PyResult<()> {
        let count = self.count;
        let noted = |error| block.noted(py, error, || format!("Writing block {count} (from 0),"));
        let alone = self.sink.accept(py, block).map_err(noted)?;
        if alone {
            self.flush(py)?;
            // Text blocks that are not ASCII, which include any that are not
            // UTF-8, are written alone, so a decode error is this block's.
            self.sink
                .write(py, block.block.serialized())
                .map_err(noted)?;
        } else {
            self.chunk.extend_from_slice(block.block.serialized());
            if self.chunk.len() >= self.limit {
                self.flush(py)?;
            }
        }
        self.count += 1;
        Ok(())
    }

    /// Write the blocks still in the chunk.
    fn flush(&mut self, py: Python<'_>) -> PyResult<()> {
        if self.chunk.is_empty() {
            return Ok(());
        }
        let written = self.sink.write(py, &self.chunk);
        self.chunk.clear();
        written
    }
}

/// Write SDBlocks to a path or a file object, and return how many.
///
/// Blocks are written as `str(block)` or `bytes(block)`, but 64-256 KiB at
/// a time: to a path without holding the GIL, or to a file object with one
/// `write` call per chunk, as str for text blocks or bytes for byte-backed
/// ones. Text blocks must be valid UTF-8.
///
/// An error stops writing: the blocks before the failing one are written,
/// a path keeps them, and the error propagates. A failed write may leave
/// part of its chunk written.
///
/// A file object's `write` is called once per chunk and its return value
/// is not checked. Buffered and text file objects, such as `sys.stdout` or
/// those `open` returns by default, write everything or raise; an
/// unbuffered (`buffering=0`) or nonblocking one may write part of a chunk.
///
/// Args:
///     dest: a str, bytes, or os.PathLike path to create or truncate, or a
///         text or binary file object
///     blocks: an iterable of SDBlocks, all str or all bytes for a file
///         object
///
/// Returns:
///     The number of blocks written
///
/// Raises:
///     TypeError: if `blocks` is not iterable, produces something other
///         than an SDBlock, or mixes str and bytes blocks for a file object
///     UnicodeDecodeError: for a str block that is not UTF-8
///     OSError: if `dest` cannot be created or written
#[pyfunction]
#[pyo3(signature = (
    dest: "str | bytes | os.PathLike[str] | os.PathLike[bytes] | _SupportsWrite[Any]",
    blocks: "Iterable[SDBlock[Any]]",
))]
pub(crate) fn write_sdf(dest: &Bound<'_, PyAny>, blocks: &Bound<'_, PyAny>) -> PyResult<usize> {
    let py = dest.py();
    let blocks = blocks.try_iter()?;
    let mut writer = Writer::open(dest)?;
    let result = blocks.into_iter().try_for_each(|block| -> PyResult<()> {
        let block = block?;
        let block = block.cast::<PySDBlock>()?.try_borrow()?;
        writer.push(py, &block)
    });
    writer.flush(py)?;
    result.map(|()| writer.count)
}

/// Write SDBlocks one at a time, buffered as `write_sdf` writes them.
///
/// `close` writes the blocks still buffered and closes a path; a file
/// object is left open. Leaving a `with` block closes the writer, also on
/// an error. Writing to a closed writer raises ValueError.
#[pyclass(name = "_SDFWriter", module = "sendoff.native")]
pub(crate) struct PySDFWriter {
    writer: Option<Writer>,
}

#[pymethods]
impl PySDFWriter {
    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[pyo3(signature = (*_exc_info))]
    fn __exit__(&mut self, py: Python<'_>, _exc_info: &Bound<'_, PyTuple>) -> PyResult<()> {
        self.close(py)
    }

    /// Write a block, ending with $$$$, as `str(block)` or `bytes(block)`.
    ///
    /// Raises:
    ///     TypeError: for a str block after bytes ones to a file object, or
    ///         the reverse
    ///     UnicodeDecodeError: for a str block that is not UTF-8
    ///     ValueError: if the writer is closed
    ///     OSError: if a buffered chunk cannot be written
    #[pyo3(signature = (block: "SDBlock[AnyStr]", /))]
    fn write(&mut self, block: &Bound<'_, PySDBlock>) -> PyResult<()> {
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| PyValueError::new_err("I/O operation on a closed sdf_writer"))?;
        let py = block.py();
        let block = block.try_borrow()?;
        writer.push(py, &block)
    }

    /// Write the buffered blocks and close a path. Closing twice does nothing.
    ///
    /// Raises:
    ///     OSError: if the buffered blocks cannot be written; the writer is
    ///         closed anyway
    fn close(&mut self, py: Python<'_>) -> PyResult<()> {
        self.writer
            .take()
            .map_or(Ok(()), |mut writer| writer.flush(py))
    }

    #[expect(
        clippy::needless_pass_by_value,
        reason = "PyO3 passes the visitor by value"
    )]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        match &self.writer {
            Some(Writer {
                sink: Sink::Handle { write, .. },
                ..
            }) => visit.call(write),
            _ => Ok(()),
        }
    }
}

impl Drop for PySDFWriter {
    /// Write what is buffered, as an unclosed Python file does, reporting
    /// an error to `sys.unraisablehook`.
    fn drop(&mut self) {
        if let Some(mut writer) = self.writer.take() {
            Python::attach(|py| {
                if let Err(error) = writer.flush(py) {
                    error.write_unraisable(py, None);
                }
            });
        }
    }
}

/// Write SDBlocks one at a time to a path or a file object, buffered.
///
/// Use it as a context manager: `with sdf_writer(dest) as out:` then
/// `out.write(block)` per block, so one pass can write several files.
/// Blocks are written as `write_sdf` writes them; an error leaves the
/// blocks before it written.
///
/// Args:
///     dest: a str, bytes, or os.PathLike path to create or truncate, or a
///         text or binary file object, which stays open
///
/// Returns:
///     A writer whose `close` writes the blocks it still buffers
///
/// Raises:
///     OSError: if `dest` cannot be created
#[pyfunction]
#[pyo3(signature = (
    dest: "str | bytes | os.PathLike[str] | os.PathLike[bytes] | _SupportsWrite[Any]",
))]
pub(crate) fn sdf_writer(dest: &Bound<'_, PyAny>) -> PyResult<PySDFWriter> {
    Ok(PySDFWriter {
        writer: Some(Writer::open(dest)?),
    })
}
