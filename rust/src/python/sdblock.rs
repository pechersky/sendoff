use pyo3::prelude::*;

use crate::core::sdblock;

#[pyfunction]
pub(crate) fn parse_mdl(lines: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
    let mut lines = lines
        .try_iter()?
        .map(|line| line.and_then(|line| line.extract::<String>()));
    sdblock::parse_mdl(&mut lines)
}

#[pyfunction]
pub(crate) fn parse_metadata(lines: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
    let mut lines = lines
        .try_iter()?
        .map(|line| line.and_then(|line| line.extract::<String>()));
    sdblock::parse_metadata(&mut lines)
}

#[pyfunction]
pub(crate) fn serialize_sdfblock(
    title: &str,
    mdl_lines: Vec<String>,
    metadata_lines: Vec<String>,
    with_newlines: bool,
) -> String {
    sdblock::serialize(title, &mdl_lines, &metadata_lines, with_newlines)
}
