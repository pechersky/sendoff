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
    use crate::python::read::PySDBlockIterator as _SDBlockIterator;

    #[pymodule_export]
    use crate::python::sdblock::PySDRecordIterator as _SDRecordIterator;

    #[pymodule_export]
    use crate::python::read::parse_sdf;

    #[pymodule_export]
    use crate::python::write::write_sdf;

    #[pymodule_export]
    #[allow(non_upper_case_globals)]
    const __version__: &str = env!("CARGO_PKG_VERSION");

    /// True when compiled without optimizations; never benchmark such a build.
    #[pymodule_export]
    const DEBUG_BUILD: bool = cfg!(debug_assertions);

    #[cfg(feature = "coverage")]
    #[pyfunction]
    #[allow(unsafe_code)]
    fn flush_coverage() -> pyo3::PyResult<()> {
        unsafe extern "C" {
            #[link_name = "__llvm_profile_dump"]
            fn dump_profile() -> std::ffi::c_int;
        }
        // SAFETY: the coverage runtime defines this, taking no arguments.
        if unsafe { dump_profile() } != 0 {
            return Err(pyo3::exceptions::PyOSError::new_err(
                "failed to write the Rust coverage profile",
            ));
        }
        Ok(())
    }
}
