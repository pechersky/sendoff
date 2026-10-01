use pyo3::prelude::*;

use crate::core::sdblock;

#[pyfunction]
pub(crate) fn parse_mdl(lines: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
    let mut lines = lines
        .try_iter()?
        .map(|line| line.and_then(|line| line.extract::<String>()));
    sdblock::parse_mdl(&mut lines)
}
