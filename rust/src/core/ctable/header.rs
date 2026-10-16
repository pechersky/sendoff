//! The MDL header: the format and the atom and bond counts.

use crate::core::ctable::error::{CTableError, Malformed, integer, malformed};
use crate::core::ctable::{CTableFormat, end_of, start_of};

/// The most lines, title included, that `parse_header` reads.
pub const HEADER_LINES: usize = 6;

/// The index of the V3000 counts line, title included.
pub(super) const V3000_COUNTS_LINE: usize = 5;

/// The index of the counts line, which ends with the format.
const COUNTS_LINE: usize = 3;

/// The line at `index`, or an error naming the first line the table lacks.
fn line_at<'a>(lines: &[&'a str], index: usize) -> Result<&'a str, CTableError> {
    lines
        .get(index)
        .copied()
        .ok_or_else(|| malformed(lines.len() + 1, Malformed::Missing(lines.len())))
}

/// A CTable's format and the atom and bond counts its header declares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Counts {
    pub format: CTableFormat,
    pub atoms: i64,
    pub bonds: i64,
}

/// Parse the format and the declared atom and bond counts.
///
/// `lines` starts with the title and needs at most six lines. V3000 input
/// skips the CTAB BEGIN line and reads its counts from the following line.
pub fn parse_header(lines: &[&str]) -> Result<Counts, CTableError> {
    let first_counts = line_at(lines, COUNTS_LINE)?;
    let format = match fields(first_counts).next_back() {
        Some("V2000") => CTableFormat::V2000,
        Some("V3000") => CTableFormat::V3000,
        Some(other) => {
            return Err(malformed(
                COUNTS_LINE + 1,
                Malformed::UnknownFormat(other.to_owned()),
            ));
        }
        None => return Err(malformed(COUNTS_LINE + 1, Malformed::NoField("format"))),
    };
    let (atoms, bonds) = match format {
        CTableFormat::V2000 => parse_v2000_counts(first_counts)?,
        CTableFormat::V3000 => {
            let (atoms, bonds) = v3000_count_fields(line_at(lines, V3000_COUNTS_LINE)?)?;
            let line = V3000_COUNTS_LINE + 1;
            (
                integer(atoms, line, "atom count")?,
                integer(bonds, line, "bond count")?,
            )
        }
    };
    Ok(Counts {
        format,
        atoms,
        bonds,
    })
}

/// Parse the atom and bond counts from a V2000 counts line, its first two
/// fixed-width three-character fields.
fn parse_v2000_counts(line: &str) -> Result<(i64, i64), CTableError> {
    // Like the Python slices line[:3] and line[3:6], count characters.
    let offset = |chars| {
        line.char_indices()
            .nth(chars)
            .map_or(line.len(), |(at, _)| at)
    };
    let (atoms, bonds) = (&line[..offset(3)], &line[offset(3)..offset(6)]);
    let at = COUNTS_LINE + 1;
    Ok((
        integer(atoms, at, "atom count")?,
        integer(bonds, at, "bond count")?,
    ))
}

/// Split a counts line into fields as Python's argument-less `str.split()`
/// does for ASCII text.
fn fields(line: &str) -> impl DoubleEndedIterator<Item = &str> {
    line.split(|character: char| {
        character.is_ascii_whitespace() || ('\u{1c}'..='\u{1f}').contains(&character)
    })
    .filter(|field| !field.is_empty())
}

/// The format `counts_line` ends with, if it ends with one.
pub fn counts_line_format(counts_line: &str) -> Option<CTableFormat> {
    match fields(counts_line).next_back()? {
        "V2000" => Some(CTableFormat::V2000),
        "V3000" => Some(CTableFormat::V3000),
        _ => None,
    }
}

/// The atom and bond count fields of the V3000 `COUNTS` line, its fourth and
/// fifth fields.
fn v3000_count_fields(line: &str) -> Result<(&str, &str), CTableError> {
    let at = V3000_COUNTS_LINE + 1;
    let mut fields = fields(line).skip(3);
    let mut field = |name| {
        fields
            .next()
            .ok_or_else(|| malformed(at, Malformed::NoField(name)))
    };
    Ok((field("atom count")?, field("bond count")?))
}

/// The V3000 `COUNTS` line `line` with its atom and bond counts replaced;
/// the rest of the line, spacing and ending included, stays as it was.
pub(super) fn with_v3000_counts(
    line: &str,
    atoms: usize,
    bonds: usize,
) -> Result<String, CTableError> {
    let (atom_field, bond_field) = v3000_count_fields(line)?;
    Ok(format!(
        "{}{atoms}{}{bonds}{}",
        &line[..start_of(line, atom_field)],
        &line[end_of(line, atom_field)..start_of(line, bond_field)],
        &line[end_of(line, bond_field)..],
    ))
}
