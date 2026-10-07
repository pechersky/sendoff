//! Python exceptions for core errors.

use pyo3::exceptions::{PyNotImplementedError, PyOSError, PyValueError};
use pyo3::prelude::*;

use crate::core::block::InvalidMolblock;
use crate::core::ctable::error::CTableError;
use crate::core::record::MissingFieldName;

pyo3::import_exception!(sendoff.ctable, IndicesDuplicateError);
pyo3::import_exception!(sendoff.ctable, IndicesMismatchError);
pyo3::import_exception!(sendoff.ctable, IndicesOutOfOrderError);
pyo3::import_exception!(sendoff.ctable, MalformedCTableError);

impl From<CTableError> for PyErr {
    fn from(error: CTableError) -> Self {
        let message = error.to_string();
        match error {
            CTableError::Malformed { .. } => MalformedCTableError::new_err(message),
            CTableError::NotImplemented(_) | CTableError::Unsupported { .. } => {
                PyNotImplementedError::new_err(message)
            }
            CTableError::OutOfOrder { .. } => IndicesOutOfOrderError::new_err(message),
            CTableError::Duplicate { .. } | CTableError::AmbiguousAtom { .. } => {
                IndicesDuplicateError::new_err(message)
            }
            CTableError::Fewer { .. } | CTableError::More { .. } => {
                IndicesMismatchError::new_err(message)
            }
        }
    }
}

impl From<MissingFieldName> for PyErr {
    fn from(error: MissingFieldName) -> Self {
        PyValueError::new_err(error.to_string())
    }
}

impl From<InvalidMolblock> for PyErr {
    fn from(error: InvalidMolblock) -> Self {
        PyValueError::new_err(error.to_string())
    }
}

/// Build the `OSError` subclass Python's `open` would raise.
pub(super) fn os_error(py: Python<'_>, error: std::io::Error, filename: &Py<PyAny>) -> PyErr {
    let Some(errno) = error.raw_os_error() else {
        return error.into();
    };
    match PyModule::import(py, "os").and_then(|os| os.call_method1("strerror", (errno,))) {
        Ok(strerror) => PyOSError::new_err((errno, strerror.unbind(), filename.clone_ref(py))),
        Err(error) => error,
    }
}
