//! A5-owned SDData exports. No metadata/write algorithm is implemented in A3.

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;

#[pyfunction]
fn _records_iter(block: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = block;
    Err(PyNotImplementedError::new_err(
        "_records_iter: unreleased A5 scaffold",
    ))
}

#[pyfunction]
fn _write(
    block: &Bound<'_, PyAny>,
    outh: &Bound<'_, PyAny>,
    with_newlines: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let _ = (block, outh, with_newlines);
    Err(PyNotImplementedError::new_err(
        "_write: unreleased A5 scaffold",
    ))
}

#[pyfunction]
fn _append_record(
    block: &Bound<'_, PyAny>,
    record_name: &Bound<'_, PyAny>,
    value: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let _ = (block, record_name, value);
    Err(PyNotImplementedError::new_err(
        "_append_record: unreleased A5 scaffold",
    ))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_records_iter, m)?)?;
    m.add_function(wrap_pyfunction!(_write, m)?)?;
    m.add_function(wrap_pyfunction!(_append_record, m)?)?;
    Ok(())
}
