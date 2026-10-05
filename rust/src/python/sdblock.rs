use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyIterator, PyList, PyString, PyTuple};

use crate::core::sdblock;
use crate::python::ctable::PyCTable;

/// Handle a single molecule block from an SD file.
#[pyclass(name = "SDBlock", module = "sendoff.native")]
pub(crate) struct PySDBlock {
    inner: SDBlockInner,
}

enum SDBlockInner {
    Text(sdblock::SDBlock),
    Bytes(sdblock::BytesSDBlock),
}

enum SDRecord {
    Text(String, String),
    Bytes(Vec<u8>, Vec<u8>),
}

#[derive(Clone, Copy)]
enum InputKind {
    Text,
    Bytes,
}

enum BufferedBlock {
    Empty,
    Text(Vec<String>),
    Bytes(Vec<Vec<u8>>),
}

#[pyclass(name = "_SDBlockIterator", module = "sendoff.native")]
pub(crate) struct PySDBlockIterator {
    lines: Py<PyIterator>,
    kind: Option<InputKind>,
    block_lines: BufferedBlock,
    finished: bool,
}

#[pyclass(name = "_SDRecordIterator", module = "sendoff.native")]
pub(crate) struct PySDRecordIterator {
    block: Py<PySDBlock>,
    index: usize,
    finished: bool,
}

#[pymethods]
impl PySDBlockIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PySDBlock>>> {
        if self.finished {
            return Ok(None);
        }

        loop {
            let mut lines = self.lines.bind(py).clone();
            let line = match lines.next() {
                Some(Ok(line)) => line,
                Some(Err(error)) => {
                    self.finish();
                    return Err(error);
                }
                None => return self.final_block(py),
            };

            let block = if line.is_instance_of::<PyString>() {
                let text = match line.extract::<String>() {
                    Ok(text) => text,
                    Err(error) => {
                        self.finish();
                        return Err(error);
                    }
                };
                self.push_text_line(text)?
            } else if line.is_instance_of::<PyBytes>() {
                let bytes = match line.extract::<Vec<u8>>() {
                    Ok(bytes) => bytes,
                    Err(error) => {
                        self.finish();
                        return Err(error);
                    }
                };
                self.push_bytes_line(bytes)?
            } else {
                self.finish();
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "SDBlock.from_lines lines must be str or bytes",
                ));
            };

            if let Some(block) = block {
                return Py::new(py, block).map(Some);
            }
        }
    }
}

#[pymethods]
impl PySDRecordIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        if self.finished {
            return Ok(None);
        }

        let mut index = self.index;
        let record = {
            let block = self.block.borrow(py);
            match &block.inner {
                SDBlockInner::Text(block) => {
                    crate::core::sddata::next_record(&block.metadata, &mut index)
                        .map(|record| record.map(|(name, value)| SDRecord::Text(name, value)))
                }
                SDBlockInner::Bytes(block) => {
                    crate::core::sddata::next_record_bytes(&block.metadata, &mut index)
                        .map(|record| record.map(|(name, value)| SDRecord::Bytes(name, value)))
                }
            }
        };
        self.index = index;

        match record {
            Ok(Some(record)) => {
                let tuple = match record {
                    SDRecord::Text(name, value) => PyTuple::new(
                        py,
                        [
                            PyString::new(py, &name).into_any(),
                            PyString::new(py, &value).into_any(),
                        ],
                    )?,
                    SDRecord::Bytes(name, value) => PyTuple::new(
                        py,
                        [
                            PyBytes::new(py, &name).into_any(),
                            PyBytes::new(py, &value).into_any(),
                        ],
                    )?,
                };
                Ok(Some(tuple.into_any().unbind()))
            }
            Ok(None) => {
                self.finished = true;
                Ok(None)
            }
            Err(error) => {
                self.finished = true;
                Err(pyo3::exceptions::PyValueError::new_err(error.to_string()))
            }
        }
    }
}

impl PySDBlockIterator {
    fn push_text_line(&mut self, line: String) -> PyResult<Option<PySDBlock>> {
        if matches!(self.kind, Some(InputKind::Bytes)) {
            self.finish();
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "SDBlock.from_lines does not accept mixed str and bytes lines",
            ));
        }
        self.kind = Some(InputKind::Text);

        let is_end = line.starts_with("$$$$");
        self.block_lines = match std::mem::replace(&mut self.block_lines, BufferedBlock::Empty) {
            BufferedBlock::Empty => BufferedBlock::Text(vec![line]),
            BufferedBlock::Text(mut lines) => {
                lines.push(line);
                BufferedBlock::Text(lines)
            }
            BufferedBlock::Bytes(_) => {
                self.finish();
                return Err(pyo3::exceptions::PyRuntimeError::new_err(
                    "SDBlock line buffer does not match its input type",
                ));
            }
        };
        if !is_end {
            return Ok(None);
        }

        match std::mem::replace(&mut self.block_lines, BufferedBlock::Empty) {
            BufferedBlock::Text(lines) => Ok(Some(PySDBlock {
                inner: SDBlockInner::Text(sdblock::SDBlock::from_block_lines(lines)),
            })),
            BufferedBlock::Empty | BufferedBlock::Bytes(_) => {
                self.finish();
                Err(pyo3::exceptions::PyRuntimeError::new_err(
                    "SDBlock line buffer does not match its input type",
                ))
            }
        }
    }

    fn push_bytes_line(&mut self, line: Vec<u8>) -> PyResult<Option<PySDBlock>> {
        if matches!(self.kind, Some(InputKind::Text)) {
            self.finish();
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "SDBlock.from_lines does not accept mixed str and bytes lines",
            ));
        }
        self.kind = Some(InputKind::Bytes);

        let is_end = line.starts_with(b"$$$$");
        self.block_lines = match std::mem::replace(&mut self.block_lines, BufferedBlock::Empty) {
            BufferedBlock::Empty => BufferedBlock::Bytes(vec![line]),
            BufferedBlock::Bytes(mut lines) => {
                lines.push(line);
                BufferedBlock::Bytes(lines)
            }
            BufferedBlock::Text(_) => {
                self.finish();
                return Err(pyo3::exceptions::PyRuntimeError::new_err(
                    "SDBlock line buffer does not match its input type",
                ));
            }
        };
        if !is_end {
            return Ok(None);
        }

        match std::mem::replace(&mut self.block_lines, BufferedBlock::Empty) {
            BufferedBlock::Bytes(lines) => Ok(Some(PySDBlock {
                inner: SDBlockInner::Bytes(sdblock::BytesSDBlock::from_block_lines(lines)),
            })),
            BufferedBlock::Empty | BufferedBlock::Text(_) => {
                self.finish();
                Err(pyo3::exceptions::PyRuntimeError::new_err(
                    "SDBlock line buffer does not match its input type",
                ))
            }
        }
    }

    fn final_block(&mut self, py: Python<'_>) -> PyResult<Option<Py<PySDBlock>>> {
        self.finished = true;
        let block = match std::mem::replace(&mut self.block_lines, BufferedBlock::Empty) {
            BufferedBlock::Empty => return Ok(None),
            BufferedBlock::Text(lines) if !lines.is_empty() => PySDBlock {
                inner: SDBlockInner::Text(sdblock::SDBlock::from_block_lines(lines)),
            },
            BufferedBlock::Bytes(lines) if !lines.is_empty() => PySDBlock {
                inner: SDBlockInner::Bytes(sdblock::BytesSDBlock::from_block_lines(lines)),
            },
            BufferedBlock::Text(_) | BufferedBlock::Bytes(_) => return Ok(None),
        };
        Py::new(py, block).map(Some)
    }

    fn finish(&mut self) {
        self.finished = true;
        self.block_lines = BufferedBlock::Empty;
    }
}

#[pymethods]
impl PySDBlock {
    /// Parse SDBlocks from homogeneous str or bytes lines.
    ///
    /// Args:
    ///     lines: an Iterable of str or bytes lines separated by the $$$$ delimiter
    ///
    /// Returns:
    ///     An iterator that produces owned SDBlocks as lines are consumed
    #[staticmethod]
    fn from_lines(py: Python<'_>, lines: &Bound<'_, PyAny>) -> PyResult<Py<PySDBlockIterator>> {
        Py::new(
            py,
            PySDBlockIterator {
                lines: lines.try_iter()?.unbind(),
                kind: None,
                block_lines: BufferedBlock::Empty,
                finished: false,
            },
        )
    }

    /// Return the parsed molecule title without changing its input type.
    #[getter]
    fn title<'py>(&self, py: Python<'py>) -> Bound<'py, PyAny> {
        match &self.inner {
            SDBlockInner::Text(block) => PyString::new(py, &block.title).into_any(),
            SDBlockInner::Bytes(block) => PyBytes::new(py, &block.title).into_any(),
        }
    }

    /// Generate SD metadata records.
    ///
    /// Returns:
    ///     Tuples of record name and value in the input type
    fn records<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let records = match &self.inner {
            SDBlockInner::Text(block) => block
                .records()
                .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?
                .into_iter()
                .map(|(name, value)| {
                    let tuple = PyTuple::new(
                        py,
                        [
                            PyString::new(py, &name).into_any(),
                            PyString::new(py, &value).into_any(),
                        ],
                    )?;
                    Ok(tuple.into_any().unbind())
                })
                .collect::<PyResult<Vec<_>>>()?,
            SDBlockInner::Bytes(block) => block
                .records()
                .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?
                .into_iter()
                .map(|(name, value)| {
                    let tuple = PyTuple::new(
                        py,
                        [
                            PyBytes::new(py, &name).into_any(),
                            PyBytes::new(py, &value).into_any(),
                        ],
                    )?;
                    Ok(tuple.into_any().unbind())
                })
                .collect::<PyResult<Vec<_>>>()?,
        };
        PyList::new(py, records)
    }

    /// Lazily iterate over the SD metadata records.
    ///
    /// Returned names and values are owned and have the same type as the block.
    fn iter_records(slf: &Bound<'_, Self>, py: Python<'_>) -> PyResult<Py<PySDRecordIterator>> {
        Py::new(
            py,
            PySDRecordIterator {
                block: slf.clone().unbind(),
                index: 0,
                finished: false,
            },
        )
    }

    /// Write an SDBlock to a file-like handle.
    ///
    /// Args:
    ///     outh: a text or binary file-like handle matching the input type
    ///     with_newlines: each line is written with a trailing newline when absent
    #[pyo3(signature = (outh, with_newlines = true))]
    fn write(&self, outh: &Bound<'_, PyAny>, with_newlines: bool) -> PyResult<()> {
        match &self.inner {
            SDBlockInner::Text(block) => {
                outh.call_method1("write", (block.serialize(with_newlines),))?;
            }
            SDBlockInner::Bytes(block) => {
                let output = block.serialize(with_newlines);
                outh.call_method1("write", (PyBytes::new(outh.py(), &output),))?;
            }
        }
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
        match &mut self.inner {
            SDBlockInner::Text(block) => {
                if !record_name.is_instance_of::<PyString>() || !value.is_instance_of::<PyString>()
                {
                    return Err(pyo3::exceptions::PyTypeError::new_err(
                        "text SDBlock values must be appended as str",
                    ));
                }
                block.append_record(
                    &record_name.extract::<String>()?,
                    &value.extract::<String>()?,
                );
            }
            SDBlockInner::Bytes(block) => {
                if !record_name.is_instance_of::<PyBytes>() || !value.is_instance_of::<PyBytes>() {
                    return Err(pyo3::exceptions::PyTypeError::new_err(
                        "byte-backed SDBlock values must be appended as bytes",
                    ));
                }
                let record_name = record_name.extract::<Vec<u8>>()?;
                let value = value.extract::<Vec<u8>>()?;
                block.append_record(&record_name, &value);
            }
        }
        Ok(())
    }

    /// Parse out the underlying CTable object from a text SDBlock.
    ///
    /// Byte-backed SDBlocks cannot be converted without decoding their title.
    fn ctable(&self) -> PyResult<PyCTable> {
        match &self.inner {
            SDBlockInner::Text(block) => PyCTable::from_lines(block.ctable_lines()),
            SDBlockInner::Bytes(_) => Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "ctable operations are not supported for byte-backed SDBlocks",
            )),
        }
    }

    /// Get the number of atoms as indicated in the MDL block.
    ///
    /// Byte-backed blocks raise `NotImplementedError`.
    fn num_atoms(&self) -> PyResult<i64> {
        Ok(self.ctable()?.num_atoms_value())
    }

    /// Get the number of bonds as indicated in the MDL block.
    ///
    /// Byte-backed blocks raise `NotImplementedError`.
    fn num_bonds(&self) -> PyResult<i64> {
        Ok(self.ctable()?.num_bonds_value())
    }

    /// Renumber atom and bond indices to be 1-indexed and in order.
    ///
    /// Byte-backed blocks raise `NotImplementedError`.
    fn renumber_indices(&mut self, py: Python<'_>) -> PyResult<()> {
        if matches!(&self.inner, SDBlockInner::Bytes(_)) {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "ctable operations are not supported for byte-backed SDBlocks",
            ));
        }
        let mut ctable = self.ctable()?;
        ctable.renumber_indices_for_block(py)?;
        if let SDBlockInner::Text(block) = &mut self.inner {
            block.mdl = ctable.lines_after_title();
        }
        Ok(())
    }
}
