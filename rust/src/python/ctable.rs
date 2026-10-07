use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;

use crate::core::ctable::{self, CTableError, Section};

pyo3::import_exception!(sendoff.ctable, IndicesDuplicateError);
pyo3::import_exception!(sendoff.ctable, IndicesMismatchError);
pyo3::import_exception!(sendoff.ctable, IndicesOutOfOrderError);
pyo3::import_exception!(sendoff.ctable, MalformedCTableError);

/// Handle a Connection table (CTAB) from a single molecule block.
///
/// Parsed in based on https://en.wikipedia.org/wiki/Chemical_table_file.
/// As part of parsing, record the title, the source line, the comment line,
/// and the counts line. The rest of the lines are kept but not parsed.
/// The counts line is parsed to infer the format and the number of atoms.
/// Note, the number of atoms in the counts line might not actually
/// be the number of lines of atoms.
///
/// Errors about the table's contents are sendoff.ctable.CTableError, a
/// ValueError; MalformedCTableError names the CTable line (the title is
/// line 1). Operations a table's format or constructs do not support raise
/// NotImplementedError.
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
    ///
    /// Raises:
    ///     MalformedCTableError: if the table is shorter than its header,
    ///         or its format or counts cannot be read
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
    /// A line ending in "-" continues onto the next, and the two are read
    /// as one atom line, numbered by its first line in errors.
    ///
    /// Args:
    ///     strict: the indices start with 1, and increment by one.
    ///
    /// Raises:
    ///     MalformedCTableError: if an atom line has no integer index, the
    ///         ATOM section has no END line or its last line ends in "-", or
    ///         the root CTAB has more than one ATOM section
    ///     IndicesDuplicateError: if an atom line has an index seen before
    ///     IndicesMismatchError: if number of atom lines does not match count line
    ///     IndicesOutOfOrderError: if strict, and indices are not in 1-indexed order
    ///     NotImplementedError: if trying to validate a V2000 format table
    ///
    /// Returns:
    ///     If all the checks pass, return True.
    #[pyo3(signature = (strict: "bool" = false))]
    fn valid_atom_indices(&self, strict: bool) -> PyResult<bool> {
        self.inner
            .valid_indices(Section::Atom, strict)
            .map(|()| true)
            .map_err(PyErr::from)
    }

    /// Validate that the bond lines match the counts line.
    ///
    /// If strict, make sure they are 1-indexed and in order.
    /// A line ending in "-" continues onto the next, and the two are read
    /// as one bond line, numbered by its first line in errors.
    ///
    /// Args:
    ///     strict: the indices start with 1, and increment by one.
    ///
    /// Raises:
    ///     MalformedCTableError: if a bond line has no integer index, the
    ///         BOND section has no END line or its last line ends in "-", or
    ///         the root CTAB has more than one BOND section
    ///     IndicesDuplicateError: if a bond line has an index seen before
    ///     IndicesMismatchError: if number of bond lines does not match count line
    ///     IndicesOutOfOrderError: if strict, and indices are not in 1-indexed order
    ///     NotImplementedError: if trying to validate a V2000 format table
    ///
    /// Returns:
    ///     If all the checks pass, return True.
    #[pyo3(signature = (strict: "bool" = false))]
    fn valid_bond_indices(&self, strict: bool) -> PyResult<bool> {
        self.inner
            .valid_indices(Section::Bond, strict)
            .map(|()| true)
            .map_err(PyErr::from)
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
    ///     MalformedCTableError: if an atom or bond line is missing a field
    ///         or has a non-integer index, a bond refers to an atom index that
    ///         no atom line has, an ATOM or BOND section has no END line or
    ///         its last line ends in "-", or the root CTAB has more than one
    ///         of either
    ///     IndicesDuplicateError: if there was a duplicate atom index, and it was used
    ///         somewhere in a bond line. Raised, because it is not clear which of the
    ///         original indices to use in the remapping.
    ///     NotImplementedError: if trying to renumber in a V2000 format table,
    ///         or a table that renumbering would corrupt, as above
    fn renumber_indices(&mut self) -> PyResult<()> {
        Ok(self.inner.renumber_indices()?)
    }

    /// Set the counts line to the number of atom and bond lines.
    ///
    /// Updates num_atoms and num_bonds to match, and keeps the rest of the
    /// counts line, such as the Sgroup count, as written. Use it after
    /// adding or removing atom or bond lines. A line ending in "-" continues
    /// onto the next, and the two count as one.
    ///
    /// Raises:
    ///     MalformedCTableError: if the counts line has no atom or bond count
    ///         field, an ATOM or BOND section has no END line or its last
    ///         line ends in "-", or the root CTAB has more than one of either
    ///     NotImplementedError: for a V2000 format table, or CTABs nested in
    ///         RGROUP or TEMPLATE blocks
    fn fix_counts(&mut self) -> PyResult<()> {
        Ok(self.inner.fix_counts()?)
    }
}

impl PyCTable {
    pub(crate) fn from_lines(lines: Vec<String>) -> PyResult<Self> {
        let inner = ctable::CTable::parse(lines)?;
        Ok(Self { inner })
    }

    pub(crate) fn lines_after_title(&self) -> &[String] {
        &self.inner.lines[1..]
    }

    pub(crate) fn renumber_indices_for_block(&mut self) -> PyResult<()> {
        self.renumber_indices()
    }

    pub(crate) fn fix_counts_for_block(&mut self) -> PyResult<()> {
        self.fix_counts()
    }
}

impl From<CTableError> for PyErr {
    fn from(error: CTableError) -> Self {
        let message = error.to_string();
        match error {
            CTableError::Malformed { .. } => MalformedCTableError::new_err(message),
            CTableError::NotImplemented(_) | CTableError::Unsupported { .. } => {
                PyNotImplementedError::new_err(message)
            }
            CTableError::OutOfOrder { .. } => IndicesOutOfOrderError::new_err(message),
            CTableError::Duplicate { .. } | CTableError::AmbiguousAtom { .. } => {
                IndicesDuplicateError::new_err(message)
            }
            CTableError::Fewer { .. } | CTableError::More { .. } => {
                IndicesMismatchError::new_err(message)
            }
        }
    }
}
