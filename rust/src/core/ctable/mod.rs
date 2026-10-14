//! Connection tables: the MDL header's format and counts, and the V3000
//! atom and bond sections that validation, renumbering, and count fixing
//! read and edit.

pub mod error;
pub mod header;
mod v3000;

use error::CTableError;
use header::{Counts, HEADER_LINES, V3000_COUNTS_LINE, parse_header, with_v3000_counts};
use v3000::{logical_lines, renumbered, section_range, sections, validate_indices};

/// The supported connection-table formats.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CTableFormat {
    /// The fixed-width V2000 format.
    V2000,
    /// The token-based V3000 format.
    V3000,
}

/// A V3000 section of the root CTAB that holds one line per atom or bond.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Section {
    Atom,
    Bond,
}

impl Section {
    fn noun(self) -> &'static str {
        match self {
            Self::Atom => "atom",
            Self::Bond => "bond",
        }
    }

    fn index_field(self) -> &'static str {
        match self {
            Self::Atom => "atom index",
            Self::Bond => "bond index",
        }
    }

    fn markers(self) -> (&'static str, &'static str) {
        match self {
            Self::Atom => ("M  V30 BEGIN ATOM", "M  V30 END ATOM"),
            Self::Bond => ("M  V30 BEGIN BOND", "M  V30 END BOND"),
        }
    }
}

/// A parsed connection table with its original source lines and parsed header.
///
/// The PyO3 binding exposes this state and its operations directly to Python.
pub struct CTable {
    /// The original lines, including the title.
    pub lines: Vec<String>,
    /// The parsed format.
    pub format: CTableFormat,
    /// The atom count from the counts line.
    pub num_atoms: i64,
    /// The bond count from the counts line.
    pub num_bonds: i64,
}

impl CTable {
    /// Parse the format and atom/bond counts from source lines.
    ///
    /// V3000 input skips the CTAB BEGIN line and reads its counts from the
    /// following line, matching the Python API's existing behavior. The
    /// stored atom and bond counts come from the counts line, not from counting
    /// the raw atom and bond lines later in the table.
    pub fn parse(lines: Vec<String>) -> Result<Self, CTableError> {
        let header: [&str; HEADER_LINES] =
            std::array::from_fn(|index| lines.get(index).map_or("", String::as_str));
        let Counts {
            format,
            atoms: num_atoms,
            bonds: num_bonds,
        } = parse_header(&header[..lines.len().min(HEADER_LINES)])?;
        Ok(Self {
            lines,
            format,
            num_atoms,
            num_bonds,
        })
    }

    /// Validate the root CTAB's atom or bond lines against the counts line.
    ///
    /// Each logical line's index must be a unique integer and, with
    /// `strict`, its one-based position; the section must have as many
    /// logical lines as the counts line declares.
    pub fn valid_indices(&self, section: Section, strict: bool) -> Result<(), CTableError> {
        let operation = match section {
            Section::Atom => "validating atom indices",
            Section::Bond => "validating bond indices",
        };
        self.require_v3000(operation)?;
        let declared = match section {
            Section::Atom => self.num_atoms,
            Section::Bond => self.num_bonds,
        };
        let range = section_range(&self.lines, section)?;
        validate_indices(&self.lines, range, section, declared, strict)
    }

    fn require_v3000(&self, operation: &'static str) -> Result<(), CTableError> {
        match self.format {
            CTableFormat::V3000 => Ok(()),
            CTableFormat::V2000 => Err(CTableError::NotImplemented(operation)),
        }
    }

    /// Renumber atom and bond indices, replacing the stored lines on success.
    ///
    /// Only V3000 tables are supported. Atom and bond lines are numbered in
    /// their current order and bond atom references are remapped; the counts
    /// line and parsed counts are left as they are (see `fix_counts`). Lines
    /// are only replaced on success.
    pub fn renumber_indices(&mut self) -> Result<(), CTableError> {
        self.require_v3000("renumbering")?;
        let (atoms, bonds) = sections("renumbering", &self.lines)?;
        let (mut new_atoms, mut new_bonds) = renumbered(&self.lines, atoms.clone(), bonds.clone())?;
        self.lines[atoms].swap_with_slice(&mut new_atoms);
        self.lines[bonds].swap_with_slice(&mut new_bonds);
        Ok(())
    }

    /// Set the counts line's atom and bond counts to the number of logical
    /// atom and bond lines, and the parsed counts with them.
    ///
    /// Only V3000 tables are supported. The rest of the counts line, spacing
    /// and ending included, stays as it was.
    pub fn fix_counts(&mut self) -> Result<(), CTableError> {
        self.require_v3000("fixing counts")?;
        let (atoms, bonds) = sections("fixing counts", &self.lines)?;
        let atoms = logical_lines(&self.lines, atoms).count();
        let bonds = logical_lines(&self.lines, bonds).count();
        self.lines[V3000_COUNTS_LINE] =
            with_v3000_counts(&self.lines[V3000_COUNTS_LINE], atoms, bonds)?;
        self.num_atoms = atoms as i64;
        self.num_bonds = bonds as i64;
        Ok(())
    }
}

/// The offset in `line` of `field`, a subslice of it.
fn start_of(line: &str, field: &str) -> usize {
    field.as_ptr() as usize - line.as_ptr() as usize
}

/// The offset in `line` just past `field`, a subslice of it.
fn end_of(line: &str, field: &str) -> usize {
    start_of(line, field) + field.len()
}
