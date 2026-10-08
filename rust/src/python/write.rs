//! Writing SDBlocks to paths and file objects a chunk at a time.

use std::fs::File;
use std::io::Write;

use pyo3::exceptions::PyTypeError;
use pyo3::intern;
use pyo3::prelude::*;

use crate::python::kind::{Kind, decode, same_kind};
use crate::python::os::{fs_path, os_error};
use crate::python::sdblock::PySDBlock;

/// Bytes buffered before each write to a path. 64 KiB to 1 MiB timed the
/// same, and 4 MiB was slower.
const WRITE_CHUNK: usize = 256 << 10;

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
                if !same_kind(expected, kind) {
                    return Err(PyTypeError::new_err(
                        "write_sdf does not accept mixed str and bytes blocks",
                    ));
                }
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
/// Blocks are written as `SDBlock.write` would, but 64-256 KiB at a time:
/// to a path without holding the GIL, or to a file object with one
/// `write` call per chunk, as str for text blocks or bytes for byte-backed
/// ones. Text blocks must be valid UTF-8.
///
/// An error stops writing as it would a loop of `SDBlock.write` calls:
/// the blocks before the failing one are written, a path keeps them, and
/// the error propagates. A failed write may leave part of its chunk
/// written.
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
    let mut sink = match dest.getattr(intern!(py, "write")) {
        Ok(write) => Sink::Handle { write, kind: None },
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
    let mut chunk = Vec::with_capacity(2 * limit);
    let mut count = 0;
    // `chunk` holds only blocks not yet passed to a write, so after an
    // error it holds the blocks before the failing one.
    let result = (|| -> PyResult<()> {
        for block in blocks {
            let block = block?;
            let block = block.cast::<PySDBlock>()?.try_borrow()?;
            let start = chunk.len();
            block.block.serialize_into(block.kind.strip(), &mut chunk);
            let noted =
                |error| block.noted(py, error, || format!("Writing block {count} (from 0),"));
            let alone = sink
                .accept(py, block.kind, &chunk[start..])
                .map_err(|error| {
                    chunk.truncate(start);
                    noted(error)
                })?;
            if alone {
                let mut written = Ok(());
                if start > 0 {
                    written = sink.write(py, &chunk[..start]);
                }
                // Text blocks that are not ASCII, which include any that
                // are not UTF-8, are written alone, so a decode error is
                // this block's.
                let written = written.and_then(|()| sink.write(py, &chunk[start..]).map_err(noted));
                chunk.clear();
                written?;
            } else if chunk.len() >= limit {
                let written = sink.write(py, &chunk);
                chunk.clear();
                written?;
            }
            count += 1;
        }
        Ok(())
    })();
    if !chunk.is_empty() {
        sink.write(py, &chunk)?;
    }
    result.map(|()| count)
}
