//! How a connection table fails to parse, validate, or change.

use crate::core::ctable::Section;

/// How a CTable line breaks its format.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Malformed {
    /// The table ends before the line, after this many lines.
    #[error("missing; the table has only {0} lines")]
    Missing(usize),
    /// The line has no field with this name.
    #[error("no {0} field")]
    NoField(&'static str),
    /// The named field, with this text, is not an integer.
    #[error("{0} {1:?} is not an integer")]
    NotInteger(&'static str, String),
    /// The counts line's last field is this, not V2000 or V3000.
    #[error("format {0:?} is not V2000 or V3000")]
    UnknownFormat(String),
    /// This begin line of a section has no end line.
    #[error("{begin:?} has no {end:?} line")]
    Unterminated {
        begin: &'static str,
        end: &'static str,
    },
    /// The root CTAB has an earlier line with this begin marker.
    #[error("more than one {0:?} line in the root CTAB")]
    Repeated(&'static str),
    /// A bond names this atom index, which no atom line has.
    #[error("bond atom {0} is not an atom index")]
    NoAtom(i64),
    /// The section's last line is continued (`-`) onto this end line.
    #[error("continues (\"-\") onto the {0:?} line")]
    Dangling(&'static str),
}

/// Errors from parsing, validating, or renumbering a connection table.
///
/// Lines are numbered from 1, the title being line 1, so a block's CTable
/// numbers lines as the block does.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CTableError {
    /// A line breaks the format.
    #[error("CTable line {line}: {problem}")]
    Malformed { line: usize, problem: Malformed },
    /// The named operation only supports V3000 tables.
    #[error("{0} is only implemented for V3000 tables")]
    NotImplemented(&'static str),
    /// The operation would corrupt or miscount the named construct, so it
    /// refuses.
    #[error("{operation} does not support {construct} yet")]
    Unsupported {
        operation: &'static str,
        construct: &'static str,
    },
    /// Strict validation found an index out of one-based order.
    #[error("CTable line {line}: {} {index} is out of order", .section.index_field())]
    OutOfOrder {
        section: Section,
        line: usize,
        index: i64,
    },
    /// An index occurs on an earlier line of the section.
    #[error("CTable line {line}: {} {index} is a duplicate", .section.index_field())]
    Duplicate {
        section: Section,
        line: usize,
        index: i64,
    },
    /// The section has fewer lines than the counts line declares.
    #[error(
        "fewer {} lines ({found}) than the counts line declares ({declared})",
        .section.noun()
    )]
    Fewer {
        section: Section,
        found: usize,
        declared: i64,
    },
    /// The section has more lines than the counts line declares.
    #[error(
        "more {} lines ({found}) than the counts line declares ({declared})",
        .section.noun()
    )]
    More {
        section: Section,
        found: usize,
        declared: i64,
    },
    /// A bond names an atom index that several atoms share, so renumbering
    /// cannot tell which atom it means.
    #[error("CTable line {line}: bond atom {index} is the index of several atoms")]
    AmbiguousAtom { line: usize, index: i64 },
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
