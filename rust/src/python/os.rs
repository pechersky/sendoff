//! Paths and OS errors, as Python's `os` module and `open` handle them.

use std::path::PathBuf;

use pyo3::exceptions::PyOSError;
use pyo3::prelude::*;

/// The `os.fspath` of a path argument, which errors name, and the path.
pub(super) fn fs_path<'py>(path: &Bound<'py, PyAny>) -> PyResult<(Bound<'py, PyAny>, PathBuf)> {
    let os = PyModule::import(path.py(), "os")?;
    let filename = os.call_method1("fspath", (path,))?;
    // fsdecode round-trips bytes paths losslessly through surrogateescape.
    let path = os.call_method1("fsdecode", (&filename,))?.extract()?;
    Ok((filename, path))
}

/// Build the `OSError` subclass Python's `open` would raise.
pub(super) fn os_error(py: Python<'_>, error: std::io::Error, filename: &Py<PyAny>) -> PyErr {
    let Some(code) = error.raw_os_error() else {
        return error.into();
    };
    os_error_from_code(py, code, &error, filename)
}

#[cfg(not(windows))]
fn os_error_from_code(
    py: Python<'_>,
    errno: i32,
    _: &std::io::Error,
    filename: &Py<PyAny>,
) -> PyErr {
    match PyModule::import(py, "os").and_then(|os| os.call_method1("strerror", (errno,))) {
        Ok(strerror) => PyOSError::new_err((errno, strerror.unbind(), filename.clone_ref(py))),
        Err(error) => error,
    }
}

/// `code` is a Windows error, which `OSError` takes as its fourth argument
/// and maps to an errno and subclass, as `PyErr_SetFromWindowsErr` does.
#[cfg(windows)]
fn os_error_from_code(
    py: Python<'_>,
    code: i32,
    error: &std::io::Error,
    filename: &Py<PyAny>,
) -> PyErr {
    let message = error.to_string();
    let message = message
        .strip_suffix(&format!(" (os error {code})"))
        .unwrap_or(&message)
        .to_owned();
    PyOSError::new_err((0, message, filename.clone_ref(py), code))
}
