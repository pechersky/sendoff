use pyo3::prelude::*;

mod core;
mod python;

#[pymodule]
fn native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    python::register(m)
}
