use pyo3::prelude::*;
use pyo3::types::{PyIterator, PyList};

use crate::core::sdblock;
use crate::python::ctable::PyCTable;

/// Handle a single molecule block from an SD file.
#[pyclass(name = "SDBlock", module = "sendoff.native")]
pub(crate) struct PySDBlock {
    inner: sdblock::SDBlock,
}

#[pymethods]
impl PySDBlock {
    /// Parse in SDBlocks from a sequence of lines.
    ///
    /// Args:
    ///     lines: an Iterable of str that optionally have whitespace that is stripped,
    ///         likely separated by the $$$$ SD delimiter
    ///
    /// Returns:
    ///     An iterator of SDBlocks parsed in from the lines
    #[staticmethod]
    fn from_lines<'py>(py: Python<'py>, lines: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let mut lines = lines
            .try_iter()?
            .map(|line| line.and_then(|line| line.extract::<String>()));
        let blocks = sdblock::SDBlock::from_lines(&mut lines)?;
        let objects = blocks
            .into_iter()
            .map(|inner| Py::new(py, Self { inner }))
            .collect::<PyResult<Vec<_>>>()?;
        let list = PyList::new(py, objects)?;
        let iterator = PyIterator::from_object(&list.into_any())?;
        Ok(iterator.into_any())
    }

    /// Return the parsed molecule title.
    #[getter]
    fn title(&self) -> &str {
        &self.inner.title
    }

    /// Generate SD metadata records one by one.
    ///
    /// Returns:
    ///     Tuples of record name and value. Multiline values include newlines.
    fn records(&self) -> Vec<(String, String)> {
        self.inner.records()
    }

    /// Write an SDBlock to a file-like handle.
    ///
    /// Args:
    ///     outh: file-like handle, such as what is returned by open(..., "w")
    ///     with_newlines: each line is written with a trailing "\\n".
    ///         The newline character is appended only if it was not in the line.
    #[pyo3(signature = (outh, with_newlines = true))]
    fn write(&self, outh: &Bound<'_, PyAny>, with_newlines: bool) -> PyResult<()> {
        outh.call_method1("write", (self.inner.serialize(with_newlines),))?;
        Ok(())
    }

    /// Append a field and value to the SD data record.
    ///
    /// Args:
    ///     record_name: field name
    ///     value: field value, with no terminating newline
    fn append_record(&mut self, record_name: &str, value: &str) {
        self.inner.append_record(record_name, value);
    }

    /// Parse out the underlying CTable object from the SDBlock.
    fn ctable(&self) -> PyResult<PyCTable> {
        PyCTable::from_lines(self.inner.ctable_lines())
    }

    /// Get the number of atoms indicated in the MDL block.
    fn num_atoms(&self) -> PyResult<i64> {
        Ok(self.ctable()?.num_atoms_value())
    }

    /// Get the number of bonds indicated in the MDL block.
    fn num_bonds(&self) -> PyResult<i64> {
        Ok(self.ctable()?.num_bonds_value())
    }

    /// Renumber atom and bond indices to be 1-indexed and in order.
    ///
    /// Utilizes the underlying CTable's renumbering and replaces this object's MDL lines.
    fn renumber_indices(&mut self, py: Python<'_>) -> PyResult<()> {
        let mut ctable = self.ctable()?;
        ctable.renumber_indices_for_block(py)?;
        self.inner.mdl = ctable.lines_after_title();
        Ok(())
    }
}

#[pyfunction]
pub(crate) fn parse_mdl(lines: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
    let mut lines = lines
        .try_iter()?
        .map(|line| line.and_then(|line| line.extract::<String>()));
    sdblock::parse_mdl(&mut lines)
}

#[pyfunction]
pub(crate) fn parse_metadata(lines: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
    let mut lines = lines
        .try_iter()?
        .map(|line| line.and_then(|line| line.extract::<String>()));
    sdblock::parse_metadata(&mut lines)
}

#[pyfunction]
pub(crate) fn serialize_sdfblock(
    title: &str,
    mdl_lines: Vec<String>,
    metadata_lines: Vec<String>,
    with_newlines: bool,
) -> String {
    sdblock::serialize(title, &mdl_lines, &metadata_lines, with_newlines)
}
