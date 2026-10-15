//! Record and record-name iterators over an `SDBlock`.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use crate::python::names::Namer;
use crate::python::sdblock::{PySDBlock, reading_record};

/// Shared cursor state for iterating over a block's records or record names,
/// raising `RuntimeError` if the block's records are edited during iteration.
struct RecordCursor {
    block: Py<PySDBlock>,
    edits: u64,
    pos: usize,
    finished: bool,
}

impl RecordCursor {
    fn new(block: Bound<'_, PySDBlock>) -> Self {
        let edits = block.borrow().edits;
        Self {
            block: block.unbind(),
            edits,
            pos: 0,
            finished: false,
        }
    }

    fn step<'py, T>(
        &mut self,
        py: Python<'py>,
        f: impl FnOnce(&PySDBlock, &mut usize) -> PyResult<Option<T>>,
    ) -> PyResult<Option<T>> {
        if self.finished {
            return Ok(None);
        }
        let block = self.block.bind(py).try_borrow()?;
        if block.edits != self.edits {
            self.finished = true;
            return Err(PyRuntimeError::new_err(
                "SDBlock records changed during iteration",
            ));
        }
        let result = f(&block, &mut self.pos);
        if !matches!(result, Ok(Some(_))) {
            self.finished = true;
        }
        result
    }
}

/// Lazily iterate over one block's SD data records.
#[pyclass(name = "_SDRecordIterator", module = "sendoff.native")]
pub(crate) struct PySDRecordIterator {
    cursor: RecordCursor,
}

impl PySDRecordIterator {
    pub(super) fn new(block: Bound<'_, PySDBlock>) -> Self {
        Self {
            cursor: RecordCursor::new(block),
        }
    }
}

#[pymethods]
impl PySDRecordIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyTuple>>> {
        self.cursor.step(py, |block, pos| {
            let mut records = block.block.records_at(*pos);
            let item = records.next();
            *pos = records.pos();
            match item {
                Some(Ok(record)) => Namer::with(block.kind, |names| names.lookup(py, record.name))
                    .and_then(|name| PyTuple::new(py, [name, block.kind.value(py, &record.value)?]))
                    .map_err(|error| block.noted(py, error, || reading_record(record.name)))
                    .map(Some),
                None => Ok(None),
                Some(Err(error)) => Err(error.into()),
            }
        })
    }
}

/// Lazily iterate over one block's record names.
#[pyclass(name = "_SDNameIterator", module = "sendoff.native")]
pub(crate) struct PySDNameIterator {
    cursor: RecordCursor,
}

impl PySDNameIterator {
    pub(super) fn new(block: Bound<'_, PySDBlock>) -> Self {
        Self {
            cursor: RecordCursor::new(block),
        }
    }
}

#[pymethods]
impl PySDNameIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&mut self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.cursor.step(py, |block, pos| {
            let mut headers = block.block.records_at(*pos).headers();
            match headers.next() {
                Some(Ok(header)) => {
                    *pos = header.end;
                    Namer::with(block.kind, |names| names.lookup(py, header.name))
                        .map_err(|error| block.noted(py, error, || reading_record(header.name)))
                        .map(Some)
                }
                None => Ok(None),
                Some(Err(error)) => Err(error.into()),
            }
        })
    }
}
