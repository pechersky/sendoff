use pyo3::prelude::*;

pub mod core;
mod python;

#[pymodule]
mod native {
    #[cfg(feature = "coverage")]
    use pyo3::prelude::*;

    #[pymodule_export]
    use crate::python::ctable::PyCTable as CTable;

    #[pymodule_export]
    use crate::python::sdblock::PySDBlock as SDBlock;

    #[pymodule_export]
    #[allow(non_upper_case_globals)]
    const __version__: &str = env!("CARGO_PKG_VERSION");

    #[cfg(feature = "coverage")]
    #[pyfunction]
    fn flush_coverage() -> pyo3::PyResult<()> {
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
}
