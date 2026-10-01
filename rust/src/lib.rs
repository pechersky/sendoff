use pyo3::prelude::*;

mod ctable;
mod framing;
mod indices;
mod sddata;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add(
        "__doc__",
        "Mandatory private bridge; unreleased algorithm scaffolding, not a completed port.",
    )?;
    framing::register(m)?;
    sddata::register(m)?;
    ctable::register(m)?;
    indices::register(m)?;
    Ok(())
}
