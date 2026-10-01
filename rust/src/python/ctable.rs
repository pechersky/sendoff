use pyo3::prelude::*;

use crate::core::ctable;

#[pyfunction]
fn parse_v2000_counts(line: &str) -> PyResult<(i64, i64)> {
    ctable::parse_v2000_counts(line)
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))
}

#[pyfunction]
fn parse_v3000_counts(line: &str) -> PyResult<(i64, i64)> {
    ctable::parse_v3000_counts(line)
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid V3000 counts"))
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(parse_v2000_counts, module)?)?;
    module.add_function(wrap_pyfunction!(parse_v3000_counts, module)?)?;
    Ok(())
}
