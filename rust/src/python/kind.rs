//! Whether a block holds str or bytes, and converting between the two.

use std::borrow::Cow;
use std::ffi::CStr;

use pyo3::exceptions::{PyUnicodeDecodeError, PyValueError};
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

use crate::core::strip::Strip;

/// How a text block decodes bytes that are not UTF-8, as `bytes.decode`'s
/// `errors` handlers do. There is no "ignore": dropping bytes can make lines
/// blank or delimiters, which the block was not split on.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Errors {
    Strict,
    Replace,
    BackslashReplace,
    SurrogateEscape,
}

impl Errors {
    const ALL: [Errors; 4] = [
        Errors::Strict,
        Errors::Replace,
        Errors::BackslashReplace,
        Errors::SurrogateEscape,
    ];

    /// The handler `bytes.decode` calls `name`.
    pub(super) fn named(name: &str) -> PyResult<Self> {
        Self::ALL
            .into_iter()
            .find(|errors| errors.name() == name)
            .ok_or_else(|| {
                PyValueError::new_err(format!(
                    "errors must be 'strict', 'replace', 'backslashreplace', or \
                     'surrogateescape', not {name:?}"
                ))
            })
    }

    pub(super) fn name(self) -> &'static str {
        self.c_name().to_str().expect("handler names are ASCII")
    }

    fn c_name(self) -> &'static CStr {
        match self {
            Errors::Strict => c"strict",
            Errors::Replace => c"replace",
            Errors::BackslashReplace => c"backslashreplace",
            Errors::SurrogateEscape => c"surrogateescape",
        }
    }
}

/// Whether a block was read from str or bytes input, which its names,
/// values, and lines have, and for str how it decodes bytes.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Kind {
    Text(Errors),
    Bytes,
}

impl Kind {
    /// Text that must be UTF-8.
    pub(super) const TEXT: Kind = Kind::Text(Errors::Strict);

    pub(super) fn is_text(self) -> bool {
        matches!(self, Kind::Text(_))
    }

    /// Whether this kind gives a str for any bytes, without raising.
    pub(super) fn is_lenient(self) -> bool {
        matches!(self, Kind::Text(errors) if errors != Errors::Strict)
    }

    /// A str or bytes object of this kind holding `bytes`.
    pub(super) fn value<'py>(self, py: Python<'py>, bytes: &[u8]) -> PyResult<Bound<'py, PyAny>> {
        match self {
            Kind::Text(Errors::Strict) => PyString::from_bytes(py, bytes).map(Bound::into_any),
            Kind::Text(errors) => match std::str::from_utf8(bytes) {
                Ok(text) => Ok(PyString::new(py, text).into_any()),
                Err(_) => decode_with(py, bytes, errors).map(Bound::into_any),
            },
            Kind::Bytes => Ok(PyBytes::new(py, bytes).into_any()),
        }
    }

    /// The Python `strip()` that values of this kind follow.
    pub(super) fn strip(self) -> Strip {
        match self {
            Kind::Text(_) => Strip::Str,
            Kind::Bytes => Strip::Bytes,
        }
    }

    /// The bytes of a str or bytes argument of this kind, or `None` for
    /// another type. A str encodes as UTF-8, with surrogateescape's
    /// surrogates back to their bytes for blocks that decode with it.
    pub(super) fn arg_bytes<'a>(
        self,
        value: &'a Bound<'_, PyAny>,
    ) -> PyResult<Option<Cow<'a, [u8]>>> {
        match self {
            Kind::Text(errors) => {
                let Ok(text) = value.cast::<PyString>() else {
                    return Ok(None);
                };
                match text.to_str() {
                    Ok(text) => Ok(Some(Cow::Borrowed(text.as_bytes()))),
                    Err(_) if errors == Errors::SurrogateEscape => {
                        let py = value.py();
                        let encoded = text.call_method1(
                            intern!(py, "encode"),
                            (intern!(py, "utf-8"), intern!(py, "surrogateescape")),
                        )?;
                        Ok(Some(Cow::Owned(
                            encoded.cast::<PyBytes>()?.as_bytes().to_vec(),
                        )))
                    }
                    Err(error) => Err(error),
                }
            }
            Kind::Bytes => Ok(value
                .cast::<PyBytes>()
                .ok()
                .map(|bytes| Cow::Borrowed(bytes.as_bytes()))),
        }
    }
}

/// `bytes` decoded as UTF-8 with the `errors` handler.
pub(super) fn decode_with<'py>(
    py: Python<'py>,
    bytes: &[u8],
    errors: Errors,
) -> PyResult<Bound<'py, PyString>> {
    PyString::from_encoded_object(
        &PyBytes::new(py, bytes),
        Some(c"utf-8"),
        Some(errors.c_name()),
    )
}

/// Return the kind and UTF-8 or raw bytes of a str or bytes object.
pub(super) fn bytes_of<'a>(value: &'a Bound<'_, PyAny>) -> PyResult<Option<(Kind, &'a [u8])>> {
    if let Ok(text) = value.cast::<PyString>() {
        Ok(Some((Kind::TEXT, text.to_str()?.as_bytes())))
    } else if let Ok(bytes) = value.cast::<PyBytes>() {
        Ok(Some((Kind::Bytes, bytes.as_bytes())))
    } else {
        Ok(None)
    }
}

/// Settle `expected` on the first kind found, returning whether `found` is
/// also str or also bytes.
pub(super) fn same_kind(expected: &mut Option<Kind>, found: Kind) -> bool {
    expected.get_or_insert(found).is_text() == found.is_text()
}

/// `bytes` as UTF-8, or the `UnicodeDecodeError` Python would raise.
pub(super) fn decode<'a>(py: Python<'_>, bytes: &'a [u8]) -> PyResult<&'a str> {
    std::str::from_utf8(bytes)
        .map_err(|error| PyUnicodeDecodeError::new_err_from_utf8(py, bytes, error))
}
