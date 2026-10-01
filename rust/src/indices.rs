//! A7-owned index exports. No validation/renumber algorithm is implemented in A3.

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;

#[pyfunction]
fn _valid_atom_indices(
    table: &Bound<'_, PyAny>,
    strict: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    errors: &Bound<'_, PyAny>,
) -> PyResult<bool> {
    let _ = (table, strict, v3000, errors);
    Err(PyNotImplementedError::new_err(
        "_valid_atom_indices: unreleased A7 scaffold",
    ))
}

#[pyfunction]
fn _valid_bond_indices(
    table: &Bound<'_, PyAny>,
    strict: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    errors: &Bound<'_, PyAny>,
) -> PyResult<bool> {
    let _ = (table, strict, v3000, errors);
    Err(PyNotImplementedError::new_err(
        "_valid_bond_indices: unreleased A7 scaffold",
    ))
}

#[pyfunction]
fn _renumber_ctable(
    table: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    duplicate_error: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let _ = (table, v3000, duplicate_error);
    Err(PyNotImplementedError::new_err(
        "_renumber_ctable: unreleased A7 scaffold",
    ))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_valid_atom_indices, m)?)?;
    m.add_function(wrap_pyfunction!(_valid_bond_indices, m)?)?;
    m.add_function(wrap_pyfunction!(_renumber_ctable, m)?)?;
    Ok(())
}
