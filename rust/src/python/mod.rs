use pyo3::{exceptions::PyUnicodeEncodeError, prelude::*, types::PyString};

mod ctable;
mod framing;
mod indices;
mod sddata;

#[cfg(feature = "coverage")]
#[pyfunction]
fn flush_coverage() -> PyResult<()> {
    unsafe extern "C" {
        #[link_name = "__llvm_profile_dump"]
        fn dump_profile() -> std::ffi::c_int;
    }
    // LLVM's dump also prevents another profile write after test cleanup.
    if unsafe { dump_profile() } != 0 {
        return Err(pyo3::exceptions::PyOSError::new_err(
            "failed to write the Rust coverage profile",
        ));
    }
    Ok(())
}

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

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    #[cfg(feature = "coverage")]
    m.add_function(wrap_pyfunction!(flush_coverage, m)?)?;
    framing::register(m)?;
    sddata::register(m)?;
    ctable::register(m)?;
    indices::register(m)
}
