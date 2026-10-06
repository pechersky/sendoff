use std::collections::{HashMap, HashSet};
use std::num::ParseIntError;

/// The supported connection-table formats.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CTableFormat {
    /// The fixed-width V2000 format.
    V2000,
    /// The token-based V3000 format.
    V3000,
}

/// Errors from parsing, validating, or renumbering a connection table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CTableError {
    /// The table has fewer lines than its header needs.
    Truncated { needed: usize, found: usize },
    /// A line lacks a required whitespace-delimited field.
    MissingField(&'static str),
    /// A field that must be an integer is not one.
    InvalidInteger(String),
    /// The final counts-line token does not name a supported format.
    UnknownFormat(String),
    /// The operation only supports V3000 tables.
    NotImplemented,
    /// Strict validation found an index out of one-based order.
    OutOfOrder,
    /// An index occurs twice, or a bond references a duplicated atom index.
    Duplicate,
    /// There are fewer unique indices than the counts line declares.
    Fewer,
    /// There are more unique indices than the counts line declares.
    More,
}

impl From<ParseIntError> for CTableError {
    fn from(error: ParseIntError) -> Self {
        Self::InvalidInteger(error.to_string())
    }
}

fn line_at<'a>(lines: &[&'a str], index: usize) -> Result<&'a str, CTableError> {
    lines.get(index).copied().ok_or(CTableError::Truncated {
        needed: index + 1,
        found: lines.len(),
    })
}

/// Parse the format, counts line, and declared atom and bond counts.
///
/// `lines` starts with the title and needs at most six lines. V3000 input
/// skips the CTAB BEGIN line and reads its stripped counts from the
/// following line.
pub fn parse_header<'a>(
    lines: &[&'a str],
) -> Result<(CTableFormat, &'a str, i64, i64), CTableError> {
    for index in 0..3 {
        line_at(lines, index)?;
    }
    let first_counts = line_at(lines, 3)?;
    let format = match parse_format(first_counts) {
        Some("V2000") => CTableFormat::V2000,
        Some("V3000") => CTableFormat::V3000,
        Some(other) => return Err(CTableError::UnknownFormat(other.to_owned())),
        None => return Err(CTableError::MissingField("format")),
    };
    let (counts, (atoms, bonds)) = match format {
        CTableFormat::V2000 => (first_counts, parse_v2000_counts(first_counts)?),
        CTableFormat::V3000 => {
            let counts = line_at(lines, 5)?.trim();
            (counts, parse_v3000_counts(counts)?)
        }
    };
    Ok((format, counts, atoms, bonds))
}

fn field<'a>(fields: &[&'a str], index: usize, name: &'static str) -> Result<&'a str, CTableError> {
    fields
        .get(index)
        .copied()
        .ok_or(CTableError::MissingField(name))
}

fn integer(text: &str) -> Result<i64, CTableError> {
    Ok(text.parse()?)
}

/// A parsed connection table with its original source lines and parsed header.
///
/// The PyO3 binding exposes this state and its operations directly to Python.
pub struct CTable {
    /// The original lines, including the title.
    pub lines: Vec<String>,
    /// The stripped molecule title.
    pub title: String,
    /// The unmodified source line.
    pub source: String,
    /// The unmodified comment line.
    pub comment: String,
    /// The original or stripped V3000 counts line.
    pub counts: String,
    /// The parsed format.
    pub format: CTableFormat,
    /// The atom count from the counts line.
    pub num_atoms: i64,
    /// The bond count from the counts line.
    pub num_bonds: i64,
}

impl CTable {
    /// Parse the title, header, format, and atom/bond counts from source lines.
    ///
    /// V3000 input skips the CTAB BEGIN line and reads its counts from the
    /// following line, matching the Python API's existing behavior. The
    /// stored atom and bond counts come from the counts line, not from counting
    /// the raw atom and bond lines later in the table.
    pub fn parse(lines: Vec<String>) -> Result<Self, CTableError> {
        let header: Vec<&str> = lines.iter().take(6).map(String::as_str).collect();
        let (format, counts, num_atoms, num_bonds) = parse_header(&header)?;
        let counts = counts.to_owned();
        let title = lines[0].trim().to_owned();
        let source = lines[1].clone();
        let comment = lines[2].clone();
        Ok(Self {
            lines,
            title,
            source,
            comment,
            counts,
            format,
            num_atoms,
            num_bonds,
        })
    }

    /// Validate raw atom lines against the parsed atom count.
    ///
    /// With `strict` enabled, indices must be one-based and in increasing
    /// order. The result distinguishes unsupported V2000 validation,
    /// out-of-order and duplicate indices, and fewer or more lines than the
    /// count line. Lines are selected using the raw V3000 atom-section
    /// boundaries.
    pub fn valid_atom_indices(&self, strict: bool) -> Result<(), CTableError> {
        if self.format != CTableFormat::V3000 {
            return Err(CTableError::NotImplemented);
        }
        validate_atom_indices(&self.lines, self.num_atoms as usize, strict)
    }

    /// Validate raw bond lines against the parsed bond count.
    ///
    /// With `strict` enabled, indices must be one-based and in increasing
    /// order. The result distinguishes unsupported V2000 validation,
    /// out-of-order and duplicate indices, and fewer or more lines than the
    /// count line. Lines are selected using the raw V3000 bond-section
    /// boundaries.
    pub fn valid_bond_indices(&self, strict: bool) -> Result<(), CTableError> {
        if self.format != CTableFormat::V3000 {
            return Err(CTableError::NotImplemented);
        }
        validate_bond_indices(&self.lines, self.num_bonds as usize, strict)
    }

    /// Renumber atom and bond indices, replacing the stored lines on success.
    ///
    /// Only V3000 tables are supported. Atom and bond lines are numbered in
    /// their current order, bond atom references are remapped, and the
    /// rewritten counts line reflects the number of table lines. Lines are
    /// only replaced on success.
    pub fn renumber_indices(&mut self) -> Result<(), CTableError> {
        if self.format != CTableFormat::V3000 {
            return Err(CTableError::NotImplemented);
        }
        self.lines = renumber_indices(&self.lines, &self.counts)?;
        Ok(())
    }
}

/// Parse the fixed-width atom and bond counts from a V2000 counts line.
///
/// The atom count uses the first three characters and the bond count uses the
/// next three. The returned counts describe the header, not the number of raw
/// table lines.
pub fn parse_v2000_counts(line: &str) -> Result<(i64, i64), ParseIntError> {
    let atom_field = line.get(..3).unwrap_or(line);
    let bond_field = line
        .get(3..6)
        .unwrap_or_else(|| line.get(3..).unwrap_or(""));
    Ok((atom_field.trim().parse()?, bond_field.trim().parse()?))
}

/// Split a counts line into fields as Python's argument-less `str.split()`
/// does for ASCII text.
fn fields(line: &str) -> impl DoubleEndedIterator<Item = &str> {
    line.split(|character: char| {
        character.is_ascii_whitespace() || ('\u{1c}'..='\u{1f}').contains(&character)
    })
    .filter(|field| !field.is_empty())
}

/// Parse atom and bond counts from the V3000 `COUNTS` line.
///
/// The atom and bond values are the fourth and fifth whitespace-delimited
/// fields.
pub fn parse_v3000_counts(line: &str) -> Result<(i64, i64), CTableError> {
    let mut fields = fields(line);
    let atoms = fields
        .nth(3)
        .ok_or(CTableError::MissingField("atom count"))?;
    let bonds = fields
        .next()
        .ok_or(CTableError::MissingField("bond count"))?;
    Ok((integer(atoms)?, integer(bonds)?))
}

/// Return the final whitespace-delimited token from a counts line.
pub fn parse_format(line: &str) -> Option<&str> {
    fields(line).next_back()
}

/// Return the half-open range of raw atom lines in a CTable.
///
/// The atom range starts at line seven and ends before `END ATOM` when that
/// marker exists. This preserves the legacy V2000 behavior of taking the
/// remaining lines after the header.
pub fn atom_line_range(lines: &[String]) -> (usize, usize) {
    let start = 7.min(lines.len());
    let end = lines
        .iter()
        .enumerate()
        .skip(start)
        .find_map(|(index, line)| line.starts_with("M  V30 END ATOM").then_some(index))
        .unwrap_or(lines.len());
    (start, end)
}

/// Return the half-open range of raw bond lines in a V3000 CTable.
///
/// If no `BEGIN BOND` marker exists, return an empty range at the end.
pub fn bond_line_range(lines: &[String]) -> (usize, usize) {
    let Some(start) = lines
        .iter()
        .position(|line| line.starts_with("M  V30 BEGIN BOND"))
        .map(|index| index + 1)
    else {
        return (lines.len(), lines.len());
    };
    let end = lines
        .iter()
        .enumerate()
        .skip(start)
        .find_map(|(index, line)| line.starts_with("M  V30 END BOND").then_some(index))
        .unwrap_or(lines.len());
    (start, end)
}

/// Validate atom-line count and unique indices and, when strict, their order.
///
/// Each raw line's third whitespace-delimited token is parsed as its index.
/// The number of unique indices is compared with `num_atoms`.
pub fn validate_atom_indices(
    lines: &[String],
    num_atoms: usize,
    strict: bool,
) -> Result<(), CTableError> {
    let (start, end) = atom_line_range(lines);
    validate_indices(&lines[start..end], num_atoms, strict)
}

/// Validate bond-line count and unique indices and, when strict, their order.
///
/// Each raw line's third whitespace-delimited token is parsed as its index.
/// The number of unique indices is compared with `num_bonds`.
pub fn validate_bond_indices(
    lines: &[String],
    num_bonds: usize,
    strict: bool,
) -> Result<(), CTableError> {
    let (start, end) = bond_line_range(lines);
    validate_indices(&lines[start..end], num_bonds, strict)
}

fn validate_indices(
    lines: &[String],
    expected_count: usize,
    strict: bool,
) -> Result<(), CTableError> {
    let mut seen_indices = HashSet::new();
    for (line_index, line) in lines.iter().enumerate() {
        let index = integer(
            line.split_whitespace()
                .nth(2)
                .ok_or(CTableError::MissingField("index"))?,
        )?;
        if strict && (line_index + 1) as i64 != index {
            return Err(CTableError::OutOfOrder);
        }
        if !seen_indices.insert(index) {
            return Err(CTableError::Duplicate);
        }
    }
    if seen_indices.len() < expected_count {
        return Err(CTableError::Fewer);
    }
    if seen_indices.len() > expected_count {
        return Err(CTableError::More);
    }
    Ok(())
}

/// Return lines with atom and bond indices renumbered in table order.
///
/// Atom and bond lines are numbered from one in their existing order, and
/// bond atom references are remapped through the atom-index map. The counts
/// line is updated to reflect the number of raw atom and bond lines. Fail with
/// `Duplicate` when duplicate atom indices make a bond's atom mapping
/// ambiguous, and with `MissingField` when a bond references an absent atom.
pub fn renumber_indices(lines: &[String], counts: &str) -> Result<Vec<String>, CTableError> {
    let prefix = "M  V30 ";
    let (atom_start, atom_end) = atom_line_range(lines);
    let atom_lines = &lines[atom_start..atom_end];
    let mut atom_index_mapping: HashMap<i64, Vec<usize>> = HashMap::new();
    let mut new_atom_lines = Vec::new();
    for (line_index, line) in atom_lines.iter().enumerate() {
        let fields: Vec<_> = line.split_whitespace().collect();
        let old_field = field(&fields, 2, "atom index")?;
        let old_index = integer(old_field)?;
        let new_index = line_index + 1;
        atom_index_mapping
            .entry(old_index)
            .or_default()
            .push(new_index);
        // Like the Python slice line[n:], skip characters, not bytes.
        let skip_chars = prefix.chars().count() + old_field.chars().count();
        let rest = line
            .char_indices()
            .nth(skip_chars)
            .map_or("", |(offset, _)| &line[offset..]);
        new_atom_lines.push(format!("{prefix}{new_index}{rest}"));
    }

    let (bond_start, bond_end) = bond_line_range(lines);
    let bond_lines = &lines[bond_start..bond_end];
    let mut new_bond_lines = Vec::new();
    for (line_index, line) in bond_lines.iter().enumerate() {
        let fields: Vec<_> = line.split_whitespace().collect();
        let from_index = integer(field(&fields, 4, "bond atom index")?)?;
        let to_index = integer(field(&fields, 5, "bond atom index")?)?;
        let from_mapping = atom_index_mapping
            .get(&from_index)
            .map_or(&[][..], Vec::as_slice);
        let to_mapping = atom_index_mapping
            .get(&to_index)
            .map_or(&[][..], Vec::as_slice);
        if from_mapping.len() > 1 || to_mapping.len() > 1 {
            return Err(CTableError::Duplicate);
        }
        let (Some(new_from), Some(new_to)) = (from_mapping.first(), to_mapping.first()) else {
            return Err(CTableError::MissingField("bonded atom"));
        };
        let trailing = if line.ends_with('\n') { "\n" } else { "" };
        new_bond_lines.push(format!(
            "{prefix}{} {} {new_from} {new_to}{trailing}",
            line_index + 1,
            field(&fields, 3, "bond type")?,
        ));
    }

    let suffix = format!(
        " {}",
        counts
            .split_whitespace()
            .skip(5)
            .collect::<Vec<_>>()
            .join(" ")
    );
    let new_counts = format!(
        "{prefix}COUNTS {} {}{suffix}",
        new_atom_lines.len(),
        new_bond_lines.len()
    );
    let mut new_lines = Vec::new();
    let mut appending = true;
    for line in lines {
        if appending {
            new_lines.push(line.clone());
        }
        if line.starts_with("M  V30 COUNTS") {
            new_lines.pop();
            new_lines.push(new_counts.clone());
        }
        if line.starts_with("M  V30 BEGIN ATOM") {
            appending = false;
        }
        if line.starts_with("M  V30 END ATOM") {
            appending = true;
            new_lines.extend(new_atom_lines.iter().cloned());
            new_lines.push(line.clone());
        }
        if line.starts_with("M  V30 BEGIN BOND") {
            appending = false;
        }
        if line.starts_with("M  V30 END BOND") {
            appending = true;
            new_lines.extend(new_bond_lines.iter().cloned());
            new_lines.push(line.clone());
        }
    }
    Ok(new_lines)
}
