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

#[pyfunction]
fn parse_format(line: &str) -> String {
    ctable::parse_format(line)
        .expect("expected a final CTable format token")
        .to_owned()
}

#[pyfunction]
fn atom_line_range(lines: Vec<String>) -> (usize, usize) {
    ctable::atom_line_range(&lines)
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(parse_v2000_counts, module)?)?;
    module.add_function(wrap_pyfunction!(parse_v3000_counts, module)?)?;
    module.add_function(wrap_pyfunction!(parse_format, module)?)?;
    module.add_function(wrap_pyfunction!(atom_line_range, module)?)?;
    Ok(())
}
