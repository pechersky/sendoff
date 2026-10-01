use pyo3::prelude::*;

mod ctable;
mod sddata;

#[cfg(feature = "coverage")]
#[pyfunction]
fn flush_coverage() -> PyResult<()> {
    unsafe extern "C" {
        #[link_name = "__llvm_profile_dump"]
        fn dump_profile() -> std::ffi::c_int;
    }
    if unsafe { dump_profile() } != 0 {
        return Err(pyo3::exceptions::PyOSError::new_err(
            "failed to write the Rust coverage profile",
        ));
    }
    Ok(())
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    ctable::register(module)?;
    sddata::register(module)?;
    #[cfg(feature = "coverage")]
    module.add_function(wrap_pyfunction!(flush_coverage, module)?)?;
    Ok(())
}
