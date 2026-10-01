use pyo3::prelude::*;

use crate::core::sddata;

#[pyfunction]
pub(crate) fn append_record(record_name: &str, value: &str) -> (String, String, String) {
    sddata::append_record(record_name, value)
}

#[pyfunction]
pub(crate) fn records(metadata: Vec<String>) -> Vec<(String, String)> {
    sddata::records(&metadata)
}
