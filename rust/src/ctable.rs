//! A6-owned CTable exports. No header/count/traversal algorithm is implemented in A3.

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;

#[pyfunction]
fn _ctable_init(
    table: &Bound<'_, PyAny>,
    lines: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let _ = (table, lines, v3000);
    Err(PyNotImplementedError::new_err(
        "_ctable_init: unreleased A6 scaffold",
    ))
}

#[pyfunction]
fn _parse_format(line: &Bound<'_, PyAny>, formats: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = (line, formats);
    Err(PyNotImplementedError::new_err(
        "_parse_format: unreleased A6 scaffold",
    ))
}

#[pyfunction]
fn _parse_v2000_counts(line: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = line;
    Err(PyNotImplementedError::new_err(
        "_parse_v2000_counts: unreleased A6 scaffold",
    ))
}

#[pyfunction]
fn _parse_v3000_counts(line: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = line;
    Err(PyNotImplementedError::new_err(
        "_parse_v3000_counts: unreleased A6 scaffold",
    ))
}

#[pyfunction]
fn _atomlines(table: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = table;
    Err(PyNotImplementedError::new_err(
        "_atomlines: unreleased A6 scaffold",
    ))
}

#[pyfunction]
fn _bondlines(table: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = table;
    Err(PyNotImplementedError::new_err(
        "_bondlines: unreleased A6 scaffold",
    ))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_ctable_init, m)?)?;
    m.add_function(wrap_pyfunction!(_parse_format, m)?)?;
    m.add_function(wrap_pyfunction!(_parse_v2000_counts, m)?)?;
    m.add_function(wrap_pyfunction!(_parse_v3000_counts, m)?)?;
    m.add_function(wrap_pyfunction!(_atomlines, m)?)?;
    m.add_function(wrap_pyfunction!(_bondlines, m)?)?;
    Ok(())
}
