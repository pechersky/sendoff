use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::core::ctable;

/// Handle a Connection table (CTAB) from a single molecule block.
///
/// Parsed in based on https://en.wikipedia.org/wiki/Chemical_table_file.
/// As part of parsing, record the title, the source line, the comment line,
/// and the counts line. The rest of the lines are kept but not parsed.
/// The counts line is parsed to infer the format and the number of atoms.
/// Note, the number of atoms in the counts line might not actually
/// be the number of lines of atoms.
#[pyclass(name = "CTable", module = "sendoff.native")]
pub(crate) struct PyCTable {
    inner: ctable::CTable,
}

#[pymethods]
impl PyCTable {
    /// Parse in lines representing the CTAB, including the title line.
    ///
    /// Args:
    ///     lines: an Iterable of str that comprises the CTAB
    ///         The title line is parsed in by the SDBlock, and should
    ///         be supplied as the first in lines.
    #[new]
    #[pyo3(signature = (lines: "Iterable[str]"))]
    fn new(lines: &Bound<'_, PyAny>) -> PyResult<Self> {
        let lines = lines
            .try_iter()?
            .map(|line| line?.extract::<String>())
            .collect::<PyResult<Vec<_>>>()?;
        Self::from_lines(lines)
    }

    /// Return a copy of the original connection-table lines.
    #[getter]
    fn lines(&self) -> &[String] {
        &self.inner.lines
    }

    /// Return the atom count from the counts line.
    #[getter]
    fn num_atoms(&self) -> i64 {
        self.inner.num_atoms
    }

    /// Return the bond count from the counts line.
    #[getter]
    fn num_bonds(&self) -> i64 {
        self.inner.num_bonds
    }

    /// Validate that the atom lines match the counts line.
    ///
    /// If strict, make sure they are 1-indexed and in order.
    /// This can break if the V3000 block has "-" terminated lines,
    ///     which means that the next line is a continuation of the previous.
    ///
    /// Args:
    ///     strict: the indices start with 1, and increment by one.
    ///
    /// Raises:
    ///     IndicesDuplicateError: if an atom line met has an index seen before
    ///     IndicesMismatchError: if number of atom lines does not match count line
    ///     IndicesOutOfOrderError: if strict, and indices are not in 1-indexed order
    ///     NotImplementedError: if trying to validate a V2000 format table
    ///
    /// Returns:
    ///     If all the checks pass, return True.
    #[pyo3(signature = (strict: "bool" = false))]
    fn valid_atom_indices(&self, py: Python<'_>, strict: bool) -> PyResult<bool> {
        self.inner
            .valid_atom_indices(strict)
            .map(|()| true)
            .map_err(|error| index_error(py, error, "atom"))
    }

    /// Validate that the bond lines match the counts line.
    ///
    /// If strict, make sure they are 1-indexed and in order.
    /// This can break if the V3000 block has "-" terminated lines,
    ///     which means that the next line is a continuation of the previous.
    ///
    /// Args:
    ///     strict: the indices start with 1, and increment by one.
    ///
    /// Raises:
    ///     IndicesDuplicateError: if an bond line met has an index seen before
    ///     IndicesMismatchError: if number of bond lines does not match count line
    ///     IndicesOutOfOrderError: if strict, and indices are not in 1-indexed order
    ///     NotImplementedError: if trying to validate a V2000 format table
    ///
    /// Returns:
    ///     If all the checks pass, return True.
    #[pyo3(signature = (strict: "bool" = false))]
    fn valid_bond_indices(&self, py: Python<'_>, strict: bool) -> PyResult<bool> {
        self.inner
            .valid_bond_indices(strict)
            .map(|()| true)
            .map_err(|error| index_error(py, error, "bond"))
    }

    /// Renumber the atom and bond indices in the block.
    ///
    /// Iterating through the atom lines, replace the indexes into a
    /// 1-indexed order, keeping track of what maps to what.
    /// Then, renumber the bond lines, making sure to map the new atom indices
    /// properly, and keeping the rest of each line, such as CFG, as written.
    /// The counts line, num_atoms, and num_bonds stay as they are; call
    /// fix_counts after adding or removing atom or bond lines.
    /// This will run regardless of whether the indices are already valid.
    /// It cannot fix duplicate atom indices properly if the duplicate one
    /// occurs in the bond lines, and will raise an error.
    /// The lines are then in-place replaced within the CTable object.
    /// Only implemented for V3000 tables.
    ///
    /// Indices that other data refers to (atom ATTCHORD, bond ENDPTS,
    /// SGROUP, COLLECTION, OBJ3D, and LINKNODE lines) are not remapped, so
    /// renumbering refuses to change them. It also refuses "-" continued
    /// atom or bond lines and CTABs nested in RGROUP or TEMPLATE blocks.
    ///
    /// Raises:
    ///     IndicesDuplicateError: if there was a duplicate atom index, and it was used
    ///         somewhere in a bond line. Raised, because it is not clear which of the
    ///         original indices to use in the remapping.
    ///     NotImplementedError: if trying to renumber in a V2000 format table,
    ///         or a table that renumbering would corrupt, as above
    fn renumber_indices(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner.renumber_indices().map_err(|error| match error {
            ctable::CTableError::Duplicate => {
                custom_error(py, "IndicesDuplicateError", "atom index mapping in bond")
            }
            other => index_error(py, other, "atom"),
        })
    }

    /// Set the counts line to the number of atom and bond lines.
    ///
    /// Updates num_atoms and num_bonds to match, and keeps the rest of the
    /// counts line, such as the Sgroup count, as written. Use it after
    /// adding or removing atom or bond lines.
    ///
    /// Raises:
    ///     NotImplementedError: for a V2000 format table, "-" continued atom
    ///         or bond lines, or CTABs nested in RGROUP or TEMPLATE blocks
    ///     ValueError: if an ATOM or BOND section has no END line
    fn fix_counts(&mut self) -> PyResult<()> {
        self.inner.fix_counts().map_err(error_to_py)
    }
}

impl PyCTable {
    pub(crate) fn from_lines(lines: Vec<String>) -> PyResult<Self> {
        let inner = ctable::CTable::parse(lines).map_err(error_to_py)?;
        Ok(Self { inner })
    }

    pub(crate) fn lines_after_title(&self) -> &[String] {
        &self.inner.lines[1..]
    }

    pub(crate) fn renumber_indices_for_block(&mut self, py: Python<'_>) -> PyResult<()> {
        self.renumber_indices(py)
    }

    pub(crate) fn fix_counts_for_block(&mut self) -> PyResult<()> {
        self.fix_counts()
    }
}

/// Map a structural CTable error to the exception main's Python raised.
///
/// Truncated tables raise ValueError rather than main's StopIteration, which
/// PEP 479 would turn into RuntimeError inside generators.
pub(crate) fn error_to_py(error: ctable::CTableError) -> PyErr {
    use pyo3::exceptions::{PyIndexError, PyKeyError, PyNotImplementedError, PyValueError};
    match error {
        ctable::CTableError::Truncated { needed, found } => PyValueError::new_err(format!(
            "truncated CTable: needs at least {needed} lines, found {found}"
        )),
        ctable::CTableError::MissingField(name) => {
            PyIndexError::new_err(format!("CTable line has no {name} field"))
        }
        ctable::CTableError::InvalidInteger(message) => PyValueError::new_err(message),
        ctable::CTableError::UnknownFormat(token) => PyKeyError::new_err(token),
        ctable::CTableError::NotImplemented => PyNotImplementedError::new_err(()),
        ctable::CTableError::MissingLine(marker) => {
            PyValueError::new_err(format!("CTable has no {marker:?} line"))
        }
        ctable::CTableError::Unsupported {
            operation,
            construct,
        } => {
            PyNotImplementedError::new_err(format!("{operation} does not support {construct} yet"))
        }
        other => PyValueError::new_err(format!("{other:?}")),
    }
}

fn index_error(py: Python<'_>, error: ctable::CTableError, noun: &str) -> PyErr {
    match error {
        ctable::CTableError::OutOfOrder => {
            custom_error(py, "IndicesOutOfOrderError", &format!("{noun}s"))
        }
        ctable::CTableError::Duplicate => {
            custom_error(py, "IndicesDuplicateError", &format!("{noun}s"))
        }
        ctable::CTableError::Fewer => custom_error(
            py,
            "IndicesMismatchError",
            &format!("fewer {noun} lines than count line"),
        ),
        ctable::CTableError::More => custom_error(
            py,
            "IndicesMismatchError",
            &format!("more {noun} lines than count line"),
        ),
        other => error_to_py(other),
    }
}

fn custom_error(py: Python<'_>, name: &str, message: &str) -> PyErr {
    PyModule::import(py, "sendoff.ctable")
        .and_then(|module| module.getattr(name))
        .and_then(|class| class.call1((message,)))
        .map_or_else(|error| error, PyErr::from_value)
}
