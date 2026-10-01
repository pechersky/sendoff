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

#[pyfunction]
fn bond_line_range(lines: Vec<String>) -> (usize, usize) {
    ctable::bond_line_range(&lines)
}

#[pyfunction]
fn validate_atom_indices(lines: Vec<String>, num_atoms: usize, strict: bool) -> Option<String> {
    ctable::validate_atom_indices(&lines, num_atoms, strict)
        .err()
        .map(str::to_owned)
}

#[pyfunction]
fn validate_bond_indices(lines: Vec<String>, num_bonds: usize, strict: bool) -> Option<String> {
    ctable::validate_bond_indices(&lines, num_bonds, strict)
        .err()
        .map(str::to_owned)
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(parse_v2000_counts, module)?)?;
    module.add_function(wrap_pyfunction!(parse_v3000_counts, module)?)?;
    module.add_function(wrap_pyfunction!(parse_format, module)?)?;
    module.add_function(wrap_pyfunction!(atom_line_range, module)?)?;
    module.add_function(wrap_pyfunction!(bond_line_range, module)?)?;
    module.add_function(wrap_pyfunction!(validate_atom_indices, module)?)?;
    module.add_function(wrap_pyfunction!(validate_bond_indices, module)?)?;
    Ok(())
}
