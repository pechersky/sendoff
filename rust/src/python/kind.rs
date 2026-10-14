//! Whether a block holds str or bytes, and converting between the two.

use pyo3::exceptions::PyUnicodeDecodeError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

/// Whether a block was read from str or bytes input, which its names,
/// values, and lines have.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Kind {
    Text,
    Bytes,
}

impl Kind {
    /// The kind of str input if `text`, else of bytes input.
    pub(super) fn of(text: bool) -> Self {
        if text { Kind::Text } else { Kind::Bytes }
    }

    /// A str or bytes object of this kind holding `bytes`.
    // TODO: explore caching short values ("0", "N/A", vendor names) as the
    // names cache does; it pays only if values repeat across blocks, so
    // measure repeat rates over many vendors' SDFs first.
    pub(super) fn value<'py>(self, py: Python<'py>, bytes: &[u8]) -> PyResult<Bound<'py, PyAny>> {
        match self {
            Kind::Text => PyString::from_bytes(py, bytes).map(Bound::into_any),
            Kind::Bytes => Ok(PyBytes::new(py, bytes).into_any()),
        }
    }
}

/// Return the kind and UTF-8 or raw bytes of a str or bytes object.
pub(super) fn bytes_of<'a>(value: &'a Bound<'_, PyAny>) -> PyResult<Option<(Kind, &'a [u8])>> {
    if let Ok(text) = value.cast::<PyString>() {
        Ok(Some((Kind::Text, text.to_str()?.as_bytes())))
    } else if let Ok(bytes) = value.cast::<PyBytes>() {
        Ok(Some((Kind::Bytes, bytes.as_bytes())))
    } else {
        Ok(None)
    }
}

/// Settle `expected` on the first kind found, returning whether `found`
/// matches it.
pub(super) fn same_kind(expected: &mut Option<Kind>, found: Kind) -> bool {
    *expected.get_or_insert(found) == found
}

/// `bytes` as UTF-8, or the `UnicodeDecodeError` Python would raise.
pub(super) fn decode<'a>(py: Python<'_>, bytes: &'a [u8]) -> PyResult<&'a str> {
    std::str::from_utf8(bytes)
        .map_err(|error| PyUnicodeDecodeError::new_err_from_utf8(py, bytes, error))
}
