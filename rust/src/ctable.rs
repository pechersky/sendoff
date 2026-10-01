use pyo3::prelude::*;
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    ffi,
    types::{PyDict, PyList, PySlice, PyString, PyTuple},
};

#[pyfunction]
fn _ctable_init(
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
        builtins
            .call_method1("next", (&iterlines,))?
            .call_method0("strip")?,
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
            builtins
                .call_method1("next", (&iterlines,))?
                .call_method0("strip")?,
        )?;
        "parse_v3000_counts"
    } else {
        "parse_v2000_counts"
    };
    let counts = table.getattr(parser)?.call1((table.getattr("counts")?,))?;
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
    drop(counts);
    table.setattr("num_atoms", &atoms)?;
    drop(atoms);
    table.setattr("num_bonds", &bonds)?;
    Ok(())
}

#[pyfunction]
fn _parse_format<'py>(
    line: &Bound<'py, PyAny>,
    formats: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    formats.get_item(
        line.call_method0("strip")?
            .call_method0("split")?
            .get_item(-1)?,
    )
}

#[pyfunction]
fn _parse_v2000_counts<'py>(
    py: Python<'py>,
    line: &Bound<'py, PyAny>,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
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
fn _parse_v3000_counts<'py>(
    py: Python<'py>,
    line: &Bound<'py, PyAny>,
) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyAny>)> {
    let tokens = line.call_method0("split")?;
    let builtins = py.import("builtins")?;
    let atoms = builtins.getattr("int")?.call1((tokens.get_item(3)?,))?;
    let bonds = builtins.getattr("int")?.call1((tokens.get_item(4)?,))?;
    Ok((atoms, bonds))
}

#[pyfunction]
fn _ctable_not_end_atom(py: Python<'_>, line: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(!py
        .get_type::<PyString>()
        .getattr("startswith")?
        .call1((line, "M  V30 END ATOM"))?
        .is_truthy()?)
}

#[pyfunction]
fn _ctable_not_begin_bond(py: Python<'_>, line: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(!py
        .get_type::<PyString>()
        .getattr("startswith")?
        .call1((line, "M  V30 BEGIN BOND"))?
        .is_truthy()?)
}

#[pyfunction]
fn _ctable_not_end_bond(py: Python<'_>, line: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(!py
        .get_type::<PyString>()
        .getattr("startswith")?
        .call1((line, "M  V30 END BOND"))?
        .is_truthy()?)
}

#[pyfunction]
fn _atomlines<'py>(py: Python<'py>, table: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    let itt = py.import("itertools")?;
    itt.getattr("takewhile")?.call1((
        wrap_pyfunction!(_ctable_not_end_atom, py)?,
        itt.getattr("islice")?
            .call1((table.getattr("lines")?, 7, py.None()))?,
    ))
}

#[pyfunction]
fn _bondlines<'py>(py: Python<'py>, table: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    let itt = py.import("itertools")?;
    itt.getattr("takewhile")?.call1((
        wrap_pyfunction!(_ctable_not_end_bond, py)?,
        itt.getattr("islice")?.call1((
            itt.getattr("dropwhile")?.call1((
                wrap_pyfunction!(_ctable_not_begin_bond, py)?,
                table.getattr("lines")?,
            ))?,
            1,
            py.None(),
        ))?,
    ))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_ctable_init, m)?)?;
    m.add_function(wrap_pyfunction!(_parse_format, m)?)?;
    m.add_function(wrap_pyfunction!(_parse_v2000_counts, m)?)?;
    m.add_function(wrap_pyfunction!(_parse_v3000_counts, m)?)?;
    m.add_function(wrap_pyfunction!(_ctable_not_end_atom, m)?)?;
    m.add_function(wrap_pyfunction!(_ctable_not_begin_bond, m)?)?;
    m.add_function(wrap_pyfunction!(_ctable_not_end_bond, m)?)?;
    m.add_function(wrap_pyfunction!(_atomlines, m)?)?;
    m.add_function(wrap_pyfunction!(_bondlines, m)?)?;
    Ok(())
}
