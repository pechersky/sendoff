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

#[pymethods]
impl PySDBlock {
    /// Parse SDBlocks from homogeneous str or bytes lines.
    ///
    /// Args:
    ///     lines: an Iterable of str or bytes lines separated by the $$$$ delimiter
    ///
    /// Returns:
    ///     An iterator of SDBlocks whose values have the same type as the input
    #[staticmethod]
    fn from_lines<'py>(py: Python<'py>, lines: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let mut text_lines = Vec::new();
        let mut byte_lines = Vec::new();
        let mut is_bytes = None;

        for line in lines.try_iter()? {
            let line = line?;
            if line.is_instance_of::<PyString>() {
                if is_bytes == Some(true) {
                    return Err(pyo3::exceptions::PyTypeError::new_err(
                        "SDBlock.from_lines does not accept mixed str and bytes lines",
                    ));
                }
                is_bytes = Some(false);
                text_lines.push(line.extract::<String>()?);
            } else if line.is_instance_of::<PyBytes>() {
                if is_bytes == Some(false) {
                    return Err(pyo3::exceptions::PyTypeError::new_err(
                        "SDBlock.from_lines does not accept mixed str and bytes lines",
                    ));
                }
                is_bytes = Some(true);
                byte_lines.push(line.extract::<Vec<u8>>()?);
            } else {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "SDBlock.from_lines lines must be str or bytes",
                ));
            }
        }

        let blocks: Vec<Self> = if is_bytes == Some(true) {
            let mut lines = byte_lines.into_iter().map(Ok::<_, PyErr>);
            sdblock::BytesSDBlock::from_lines(&mut lines)?
                .into_iter()
                .map(|inner| Self {
                    inner: SDBlockInner::Bytes(inner),
                })
                .collect()
        } else {
            let mut lines = text_lines.into_iter().map(Ok::<_, PyErr>);
            sdblock::SDBlock::from_lines(&mut lines)?
                .into_iter()
                .map(|inner| Self {
                    inner: SDBlockInner::Text(inner),
                })
                .collect()
        };
        let objects = blocks
            .into_iter()
            .map(|block| Py::new(py, block))
            .collect::<PyResult<Vec<_>>>()?;
        let list = PyList::new(py, objects)?;
        let iterator = PyIterator::from_object(&list.into_any())?;
        Ok(iterator.into_any())
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
