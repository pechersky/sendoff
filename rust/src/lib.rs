use pyo3::{exceptions::PyUnicodeEncodeError, prelude::*, types::PyString};

mod ctable;
mod framing;
mod indices;
mod sddata;

fn exact_text<'a>(value: &'a Bound<'_, PyAny>) -> PyResult<Option<&'a str>> {
    if !value.is_exact_instance_of::<PyString>() {
        return Ok(None);
    }
    unicode_text(value.cast::<PyString>()?)
}

fn unicode_text<'a>(value: &'a Bound<'_, PyString>) -> PyResult<Option<&'a str>> {
    match value.to_str() {
        Ok(text) => Ok(Some(text)),
        // Lone surrogates retain their Python Unicode semantics at the boundary.
        Err(error) if error.is_instance_of::<PyUnicodeEncodeError>(value.py()) => Ok(None),
        Err(error) => Err(error),
    }
}

fn whitespace(c: char) -> bool {
    c.is_whitespace() || matches!(c, '\u{1c}'..='\u{1f}')
}

#[pymodule]
fn native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    // ponytail: leaves add their own exports; no algorithm placeholders.
    framing::register(m)?;
    sddata::register(m)?;
    ctable::register(m)?;
    indices::register(m)?;
    Ok(())
}
