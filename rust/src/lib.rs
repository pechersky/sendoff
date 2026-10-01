use pyo3::prelude::*;

pub mod core;
mod python;

#[pymodule]
fn native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    python::register(module)
}
