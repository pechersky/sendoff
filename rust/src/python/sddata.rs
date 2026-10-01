use pyo3::prelude::*;

use crate::core::sddata;

#[pyfunction]
fn append_record(record_name: &str, value: &str) -> (String, String, String) {
    sddata::append_record(record_name, value)
}

pub fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(append_record, module)?)?;
    Ok(())
}
