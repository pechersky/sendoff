use pyo3::prelude::*;

use crate::core::ctable;

#[pyfunction]
fn _parse_v2000_counts(line: &str) -> PyResult<(i64, i64)> {
    ctable::parse_v2000_counts(line)
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))
}

#[pyfunction]
fn _parse_v3000_counts(line: &str) -> Option<(i64, i64)> {
    ctable::parse_v3000_counts(line)
}

#[pyfunction]
fn _parse_format(line: &str) -> Option<String> {
    ctable::parse_format(line).map(str::to_owned)
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(_parse_v2000_counts, module)?)?;
    module.add_function(wrap_pyfunction!(_parse_v3000_counts, module)?)?;
    module.add_function(wrap_pyfunction!(_parse_format, module)?)?;
    Ok(())
}
