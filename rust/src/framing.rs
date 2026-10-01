//! A4-owned framing exports. No framing algorithm is implemented in A3.

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;

#[pyfunction]
fn _mdl_iter(lines: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = lines;
    Err(PyNotImplementedError::new_err(
        "_mdl_iter: unreleased A4 scaffold",
    ))
}

#[pyfunction]
fn _metadata_iter(lines: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = lines;
    Err(PyNotImplementedError::new_err(
        "_metadata_iter: unreleased A4 scaffold",
    ))
}

#[pyfunction]
fn _from_block_lines(
    cls: &Bound<'_, PyAny>,
    block_type: &Bound<'_, PyAny>,
    lines: &Bound<'_, PyAny>,
) -> PyResult<Py<PyAny>> {
    let _ = (cls, block_type, lines);
    Err(PyNotImplementedError::new_err(
        "_from_block_lines: unreleased A4 scaffold",
    ))
}

#[pyfunction]
fn _blocks_iter(cls: &Bound<'_, PyAny>, lines: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = (cls, lines);
    Err(PyNotImplementedError::new_err(
        "_blocks_iter: unreleased A4 scaffold",
    ))
}

#[pyfunction]
fn _read_sdf_lines(sdfpath: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let _ = sdfpath;
    Err(PyNotImplementedError::new_err(
        "_read_sdf_lines: unreleased A4 scaffold",
    ))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_mdl_iter, m)?)?;
    m.add_function(wrap_pyfunction!(_metadata_iter, m)?)?;
    m.add_function(wrap_pyfunction!(_from_block_lines, m)?)?;
    m.add_function(wrap_pyfunction!(_blocks_iter, m)?)?;
    m.add_function(wrap_pyfunction!(_read_sdf_lines, m)?)?;
    Ok(())
}
