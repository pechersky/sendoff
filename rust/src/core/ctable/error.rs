//! How a connection table fails to parse, validate, or change.

use std::fmt;

use crate::core::ctable::Section;

/// How a CTable line breaks its format.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Malformed {
    /// The table ends before the line, after this many lines.
    Missing(usize),
    /// The line has no field with this name.
    NoField(&'static str),
    /// The named field, with this text, is not an integer.
    NotInteger(&'static str, String),
    /// The counts line's last field is this, not V2000 or V3000.
    UnknownFormat(String),
    /// This begin line of a section has no end line.
    Unterminated {
        begin: &'static str,
        end: &'static str,
    },
    /// The root CTAB has an earlier line with this begin marker.
    Repeated(&'static str),
    /// A bond names this atom index, which no atom line has.
    NoAtom(i64),
    /// The section's last line is continued (`-`) onto this end line.
    Dangling(&'static str),
}

/// Errors from parsing, validating, or renumbering a connection table.
///
/// Lines are numbered from 1, the title being line 1, so a block's CTable
/// numbers lines as the block does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CTableError {
    /// A line breaks the format.
    Malformed { line: usize, problem: Malformed },
    /// The named operation only supports V3000 tables.
    NotImplemented(&'static str),
    /// The operation would corrupt or miscount the named construct, so it
    /// refuses.
    Unsupported {
        operation: &'static str,
        construct: &'static str,
    },
    /// Strict validation found an index out of one-based order.
    OutOfOrder {
        section: Section,
        line: usize,
        index: i64,
    },
    /// An index occurs on an earlier line of the section.
    Duplicate {
        section: Section,
        line: usize,
        index: i64,
    },
    /// The section has fewer lines than the counts line declares.
    Fewer {
        section: Section,
        found: usize,
        declared: i64,
    },
    /// The section has more lines than the counts line declares.
    More {
        section: Section,
        found: usize,
        declared: i64,
    },
    /// A bond names an atom index that several atoms share, so renumbering
    /// cannot tell which atom it means.
    AmbiguousAtom { line: usize, index: i64 },
}

impl fmt::Display for Malformed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing(found) => write!(f, "missing; the table has only {found} lines"),
            Self::NoField(field) => write!(f, "no {field} field"),
            Self::NotInteger(field, text) => write!(f, "{field} {text:?} is not an integer"),
            Self::UnknownFormat(format) => {
                write!(f, "format {format:?} is not V2000 or V3000")
            }
            Self::Unterminated { begin, end } => write!(f, "{begin:?} has no {end:?} line"),
            Self::Repeated(begin) => {
                write!(f, "more than one {begin:?} line in the root CTAB")
            }
            Self::NoAtom(index) => write!(f, "bond atom {index} is not an atom index"),
            Self::Dangling(end) => write!(f, "continues (\"-\") onto the {end:?} line"),
        }
    }
}

impl fmt::Display for CTableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Malformed { line, problem } => write!(f, "CTable line {line}: {problem}"),
            Self::NotImplemented(operation) => {
                write!(f, "{operation} is only implemented for V3000 tables")
            }
            Self::Unsupported {
                operation,
                construct,
            } => write!(f, "{operation} does not support {construct} yet"),
            Self::OutOfOrder {
                section,
                line,
                index,
            } => write!(
                f,
                "CTable line {line}: {} {index} is out of order",
                section.index_field()
            ),
            Self::Duplicate {
                section,
                line,
                index,
            } => write!(
                f,
                "CTable line {line}: {} {index} is a duplicate",
                section.index_field()
            ),
            Self::Fewer {
                section,
                found,
                declared,
            } => write!(
                f,
                "fewer {} lines ({found}) than the counts line declares ({declared})",
                section.noun()
            ),
            Self::More {
                section,
                found,
                declared,
            } => write!(
                f,
                "more {} lines ({found}) than the counts line declares ({declared})",
                section.noun()
            ),
            Self::AmbiguousAtom { line, index } => write!(
                f,
                "CTable line {line}: bond atom {index} is the index of several atoms"
            ),
        }
    }
}

/// A `Malformed` error for CTable line `line`.
pub(super) fn malformed(line: usize, problem: Malformed) -> CTableError {
    CTableError::Malformed { line, problem }
}

/// Parse `text`, the named field of CTable line `line`, as an integer.
pub(super) fn integer(text: &str, line: usize, field: &'static str) -> Result<i64, CTableError> {
    text.trim()
        .parse()
        .map_err(|_| malformed(line, Malformed::NotInteger(field, text.to_owned())))
}
