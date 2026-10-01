use super::exact_text;
use crate::core::ctable;
use pyo3::prelude::*;
use pyo3::{exceptions::PyIndexError, types::PySlice};

fn integer<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyAny>> {
    if let Some(value) = ctable::integer(text) {
        return Ok(value.into_pyobject(py)?.into_any());
    }
    // Unicode digits, underscores, overflow and diagnostics retain Python's int rules.
    py.import("builtins")?.call_method1("int", (text,))
}

#[pyfunction]
fn parse_format<'py>(
    line: &Bound<'py, PyAny>,
    formats: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    if let Some(text) = exact_text(line)? {
        let format =
            ctable::format(text).ok_or_else(|| PyIndexError::new_err("list index out of range"))?;
        return formats.get_item(format);
    }
    formats.get_item(
        line.call_method0("strip")?
            .call_method0("split")?
            .get_item(-1)?,
    )
}

#[pyfunction]
fn parse_v2000_counts<'py>(
    py: Python<'py>,
    line: &Bound<'py, PyAny>,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
    if let Some(text) = exact_text(line)? {
        if let Some((atoms, bonds)) = ctable::parse_v2000_counts(text) {
            return Ok((
                atoms.into_pyobject(py)?.into_any(),
                bonds.into_pyobject(py)?.into_any(),
            ));
        }
        let (atoms, bonds) = ctable::v2000_fields(text);
        return Ok((integer(py, atoms)?, integer(py, bonds)?));
    }
    let builtins = py.import("builtins")?;
    let atoms = builtins
        .getattr("int")?
        .call1((line.get_item(py.get_type::<PySlice>().call1((py.None(), 3, py.None()))?)?,))?;
    let bonds = builtins
        .getattr("int")?
        .call1((line.get_item(py.get_type::<PySlice>().call1((3, 6, py.None()))?)?,))?;
    Ok((atoms, bonds))
}

#[pyfunction]
fn parse_v3000_counts<'py>(
    py: Python<'py>,
    line: &Bound<'py, PyAny>,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
    if let Some(text) = exact_text(line)? {
        if let Some((atoms, bonds)) = ctable::parse_v3000_counts(text) {
            return Ok((
                atoms.into_pyobject(py)?.into_any(),
                bonds.into_pyobject(py)?.into_any(),
            ));
        }
        let atoms = integer(
            py,
            ctable::token(text, 3)
                .ok_or_else(|| PyIndexError::new_err("list index out of range"))?,
        )?;
        let bonds = integer(
            py,
            ctable::token(text, 4)
                .ok_or_else(|| PyIndexError::new_err("list index out of range"))?,
        )?;
        return Ok((atoms, bonds));
    }
    let tokens = line.call_method0("split")?;
    let builtins = py.import("builtins")?;
    let atoms = builtins.getattr("int")?.call1((tokens.get_item(3)?,))?;
    let bonds = builtins.getattr("int")?.call1((tokens.get_item(4)?,))?;
    Ok((atoms, bonds))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_format, m)?)?;
    m.add_function(wrap_pyfunction!(parse_v2000_counts, m)?)?;
    m.add_function(wrap_pyfunction!(parse_v3000_counts, m)?)?;
    Ok(())
}
