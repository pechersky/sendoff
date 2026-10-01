use pyo3::prelude::*;

mod ctable;
mod framing;
mod indices;
mod sddata;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    // ponytail: leaves add their own exports; no algorithm placeholders.
    framing::register(m)?;
    sddata::register(m)?;
    ctable::register(m)?;
    indices::register(m)?;
    Ok(())
}
