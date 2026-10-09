//! Python exceptions for core errors.

use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;

use crate::core::block::{InvalidMolblock, InvalidTitle};
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

impl From<InvalidTitle> for PyErr {
    fn from(error: InvalidTitle) -> Self {
        PyValueError::new_err(error.to_string())
    }
}
