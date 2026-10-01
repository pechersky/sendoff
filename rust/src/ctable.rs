use crate::{exact_text, unicode_text, whitespace};
use pyo3::prelude::*;
use pyo3::{
    exceptions::{PyIndexError, PyTypeError, PyValueError},
    ffi,
    types::{PyDict, PyList, PySlice, PyString, PyTuple},
};

fn strip<'py>(value: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    match exact_text(value)? {
        Some(text) => Ok(PyString::new(value.py(), text.trim_matches(whitespace)).into_any()),
        None => value.call_method0("strip"),
    }
}

fn integer<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyAny>> {
    if text.len() <= 40
        && let Ok(value) = text.trim().parse::<i128>()
    {
        return Ok(value.into_pyobject(py)?.into_any());
    }
    // Unicode digits, underscores, overflow and diagnostics retain Python's int rules.
    py.import("builtins")?.call_method1("int", (text,))
}

#[pyfunction]
fn ctable_init(
    py: Python<'_>,
    table: &Bound<'_, PyAny>,
    lines: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
) -> PyResult<()> {
    table.setattr(
        "lines",
        py.import("collections")?
            .getattr("deque")?
            .call1((lines,))?,
    )?;
    let builtins = py.import("builtins")?;
    let iterlines = builtins.call_method1("iter", (table.getattr("lines")?,))?;
    table.setattr(
        "title",
        strip(&builtins.call_method1("next", (&iterlines,))?)?,
    )?;
    for name in ["source", "comment", "counts"] {
        table.setattr(name, builtins.call_method1("next", (&iterlines,))?)?;
    }
    table.setattr(
        "format",
        table
            .getattr("parse_format")?
            .call1((table.getattr("counts")?,))?,
    )?;
    let parser = if table.getattr("format")?.is(v3000) {
        builtins.call_method1("next", (&iterlines,))?;
        table.setattr(
            "counts",
            strip(&builtins.call_method1("next", (&iterlines,))?)?,
        )?;
        "parse_v3000_counts"
    } else {
        "parse_v2000_counts"
    };
    let counts = table.getattr(parser)?.call1((table.getattr("counts")?,))?;
    let (atoms, bonds) = unpack_counts(py, &counts)?;
    drop(counts);
    table.setattr("num_atoms", &atoms)?;
    drop(atoms);
    table.setattr("num_bonds", &bonds)?;
    Ok(())
}

fn unpack_counts<'py>(
    py: Python<'py>,
    counts: &Bound<'py, PyAny>,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
    if counts.is_exact_instance_of::<PyTuple>() && counts.len()? == 2 {
        return Ok((counts.get_item(0)?, counts.get_item(1)?));
    }
    let mut items = counts.try_iter().map_err(|error| {
        // CPython unpacking only replaces TypeError when no iteration protocol exists.
        let non_iterable = unsafe {
            ffi::PyType_GetSlot(counts.get_type().as_type_ptr(), ffi::Py_tp_iter).is_null()
                && ffi::PySequence_Check(counts.as_ptr()) == 0
        };
        if error.is_instance_of::<PyTypeError>(py) && non_iterable {
            match error.value(py).call_method0("__str__").and_then(|message| {
                let name = message
                    .call_method1("removeprefix", ("'",))?
                    .call_method1("removesuffix", ("' object is not iterable",))?;
                PyString::new(py, "cannot unpack non-iterable {} object")
                    .call_method1("format", (name,))
            }) {
                Ok(message) => PyTypeError::new_err(message.unbind()),
                Err(error) => error,
            }
        } else {
            error
        }
    })?;
    let atoms = match items.next() {
        Some(value) => value?,
        None => {
            return Err(PyValueError::new_err(
                "not enough values to unpack (expected 2, got 0)",
            ));
        }
    };
    let bonds = match items.next() {
        Some(value) => value?,
        None => {
            return Err(PyValueError::new_err(
                "not enough values to unpack (expected 2, got 1)",
            ));
        }
    };
    if let Some(extra) = items.next() {
        extra?;
        if py.version_info().minor >= 14
            && (counts.is_exact_instance_of::<PyList>()
                || counts.is_exact_instance_of::<PyTuple>()
                || counts.is_exact_instance_of::<PyDict>())
        {
            return Err(PyValueError::new_err(format!(
                "too many values to unpack (expected 2, got {})",
                counts.len()?
            )));
        }
        return Err(PyValueError::new_err(
            "too many values to unpack (expected 2)",
        ));
    }
    drop(items);
    Ok((atoms, bonds))
}

#[pyfunction]
fn parse_format<'py>(
    line: &Bound<'py, PyAny>,
    formats: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    if let Some(text) = exact_text(line)? {
        let format = text
            .split(whitespace)
            .rfind(|token| !token.is_empty())
            .ok_or_else(|| PyIndexError::new_err("list index out of range"))?;
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
        let middle = text
            .char_indices()
            .nth(3)
            .map_or(text.len(), |(offset, _)| offset);
        let end = text
            .char_indices()
            .nth(6)
            .map_or(text.len(), |(offset, _)| offset);
        return Ok((
            integer(py, &text[..middle])?,
            integer(py, &text[middle..end])?,
        ));
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
        let mut tokens = text.split(whitespace).filter(|token| !token.is_empty());
        let atoms = integer(
            py,
            tokens
                .nth(3)
                .ok_or_else(|| PyIndexError::new_err("list index out of range"))?,
        )?;
        let bonds = integer(
            py,
            tokens
                .next()
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

#[pyfunction]
fn ctable_not_end_atom(py: Python<'_>, line: &Bound<'_, PyAny>) -> PyResult<bool> {
    not_prefix(py, line, "M  V30 END ATOM")
}

#[pyfunction]
fn ctable_not_begin_bond(py: Python<'_>, line: &Bound<'_, PyAny>) -> PyResult<bool> {
    not_prefix(py, line, "M  V30 BEGIN BOND")
}

#[pyfunction]
fn ctable_not_end_bond(py: Python<'_>, line: &Bound<'_, PyAny>) -> PyResult<bool> {
    not_prefix(py, line, "M  V30 END BOND")
}

fn not_prefix(py: Python<'_>, line: &Bound<'_, PyAny>, prefix: &str) -> PyResult<bool> {
    if line.is_instance_of::<PyString>()
        && let Some(text) = unicode_text(line.cast::<PyString>()?)?
    {
        return Ok(!text.starts_with(prefix));
    }
    Ok(!py
        .get_type::<PyString>()
        .getattr("startswith")?
        .call1((line, prefix))?
        .is_truthy()?)
}

#[pyfunction]
fn atomlines<'py>(py: Python<'py>, table: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    let itt = py.import("itertools")?;
    itt.getattr("takewhile")?.call1((
        wrap_pyfunction!(ctable_not_end_atom, py)?,
        itt.getattr("islice")?
            .call1((table.getattr("lines")?, 7, py.None()))?,
    ))
}

#[pyfunction]
fn bondlines<'py>(py: Python<'py>, table: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    let itt = py.import("itertools")?;
    itt.getattr("takewhile")?.call1((
        wrap_pyfunction!(ctable_not_end_bond, py)?,
        itt.getattr("islice")?.call1((
            itt.getattr("dropwhile")?.call1((
                wrap_pyfunction!(ctable_not_begin_bond, py)?,
                table.getattr("lines")?,
            ))?,
            1,
            py.None(),
        ))?,
    ))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ctable_init, m)?)?;
    m.add_function(wrap_pyfunction!(parse_format, m)?)?;
    m.add_function(wrap_pyfunction!(parse_v2000_counts, m)?)?;
    m.add_function(wrap_pyfunction!(parse_v3000_counts, m)?)?;
    m.add_function(wrap_pyfunction!(ctable_not_end_atom, m)?)?;
    m.add_function(wrap_pyfunction!(ctable_not_begin_bond, m)?)?;
    m.add_function(wrap_pyfunction!(ctable_not_end_bond, m)?)?;
    m.add_function(wrap_pyfunction!(atomlines, m)?)?;
    m.add_function(wrap_pyfunction!(bondlines, m)?)?;
    Ok(())
}
