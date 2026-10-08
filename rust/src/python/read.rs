//! Reading SDBlocks lazily from paths, file objects, lines, and buffers.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{ErrorKind, Read};
use std::path::PathBuf;

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyIterator, PyString};
use pyo3::{PyTraverseError, PyVisit};

use crate::core::block::Block;
use crate::core::mapped::MappedBlocks;
use crate::core::newline::{TRANSLATE_SLACK, missing_ending, title_newline, translate_newlines};
use crate::core::split::{Scanner, Tail};
use crate::python::kind::{Kind, bytes_of, same_kind};
use crate::python::os::{fs_path, os_error};
use crate::python::sdblock::PySDBlock;

/// Characters or bytes requested per `read` call, or split per batch from a
/// mapped file. Blocks larger than this grow the scanner's buffer. Against
/// 1 MiB, 128 KiB was as fast for paths and buffers, 25-30% faster for text
/// file objects, and used less memory.
const CHUNK: usize = 128 << 10;

/// Settle `kind` on the kind of `from_lines` input found, which must not
/// change.
fn set_kind(kind: &mut Option<Kind>, found: Kind) -> PyResult<()> {
    if same_kind(kind, found) {
        Ok(())
    } else {
        Err(PyTypeError::new_err(
            "SDBlock.from_lines does not accept mixed str and bytes lines",
        ))
    }
}

/// Where `_SDBlockIterator` reads SDF from.
pub(super) enum Source {
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
    /// A file mapped from `Path`, split a batch of blocks at a time without
    /// holding the GIL.
    Mapped {
        blocks: MappedBlocks,
        batch: VecDeque<Block>,
    },
    /// A file opened from `Path` that does not map, read without holding
    /// the GIL.
    File {
        file: File,
        filename: Py<PyAny>,
        text: bool,
        pending_cr: bool,
    },
    Closed,
}

impl Source {
    /// Open a `Path`, mapping the file when it maps.
    fn open(&mut self, py: Python<'_>, kind: &mut Option<Kind>) -> PyResult<()> {
        let Source::Path {
            path,
            filename,
            text,
        } = self
        else {
            return Ok(());
        };
        let text_mode = *text;
        // Opening a FIFO waits for a writer, which may be a Python thread.
        let (file, mapped) = py
            .detach(|| {
                let file = File::open(&*path)?;
                let mapped = MappedBlocks::open(&file, text_mode);
                Ok((file, mapped))
            })
            .map_err(|error| os_error(py, error, filename))?;
        *kind = Some(Kind::of(*text));
        *self = match mapped {
            Some(blocks) => Source::Mapped {
                blocks,
                batch: VecDeque::new(),
            },
            None => Source::File {
                file,
                filename: filename.clone_ref(py),
                text: *text,
                pending_cr: false,
            },
        };
        Ok(())
    }

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
                // Lines given without an ending get the block's title line's.
                let mut newline = None;
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
                    let newline = *newline.get_or_insert_with(|| title_newline(bytes));
                    if !bytes.ends_with(b"\n") {
                        tail.extend(missing_ending(bytes, newline));
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
            Source::Path { .. } | Source::Mapped { .. } => {
                unreachable!("PySDBlockIterator opens paths and splits mapped files")
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
    pub(super) fn new(source: Source) -> Self {
        Self {
            source,
            kind: None,
            scanner: Scanner::default(),
        }
    }

    pub(super) fn of_kind(source: Source, kind: Kind) -> Self {
        Self {
            kind: Some(kind),
            ..Self::new(source)
        }
    }

    fn next_block(&mut self, py: Python<'_>) -> PyResult<Option<Block>> {
        let Self {
            source,
            kind,
            scanner,
        } = self;
        source.open(py, kind)?;
        if let Source::Mapped { blocks, batch } = source {
            if batch.is_empty() {
                py.detach(|| blocks.split(CHUNK, batch));
            }
            return Ok(batch.pop_front());
        }
        scanner.next_block(|tail| source.fill(py, kind, tail))
    }
}

#[pymethods]
impl PySDBlockIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PySDBlock>> {
        match self.next_block(py) {
            Ok(Some(block)) => Ok(Some(PySDBlock {
                block,
                kind: self.kind.unwrap_or(Kind::Text),
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

    /// Let the garbage collector find cycles through the source, such as a
    /// reader that keeps its parser.
    ///
    /// There is no `__clear__`: a cycle back to the parser runs through a
    /// container changed after the parser was made, and clearing that
    /// container breaks it.
    #[expect(
        clippy::needless_pass_by_value,
        reason = "PyO3 passes the visitor by value"
    )]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        // The other sources hold only str and bytes, which are never in cycles.
        match &self.source {
            Source::Lines(lines) => visit.call(lines),
            Source::Reader(read) => visit.call(read),
            _ => Ok(()),
        }
    }
}

/// Lazily parse SDBlocks from an SDF path.
///
/// The file is opened on the first `next()`. Regular files on local Linux
/// filesystems are memory-mapped and split a batch of blocks at a time
/// without holding the GIL; pipes, network filesystems, and other files are
/// read in chunks. As with any mapped file, truncating the file while
/// iterating can kill the process with SIGBUS, so write results to a new
/// path. str paths give str values decoded as UTF-8 on access, with
/// universal newlines translated to "\n", as `open(path)` reads; bytes
/// paths give bytes split on b"\n" only. Use a bytes path for feeds that
/// are not valid UTF-8: from a str path, their text raises
/// UnicodeDecodeError with a note naming the block and record.
///
/// Args:
///     path: a str, bytes, or os.PathLike path
///
/// Returns:
///     An iterator of SDBlocks; `close()` releases the file early
#[pyfunction]
#[pyo3(signature = (path: "str | bytes | os.PathLike[str] | os.PathLike[bytes]") -> "_SDBlockIterator[Any]")]
pub(crate) fn parse_sdf(path: &Bound<'_, PyAny>) -> PyResult<PySDBlockIterator> {
    let (filename, path) = fs_path(path)?;
    let text = filename.is_instance_of::<PyString>();
    Ok(PySDBlockIterator::new(Source::Path {
        path,
        filename: filename.unbind(),
        text,
    }))
}
