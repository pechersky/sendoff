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

/// Errors raised while parsing a connection table.
#[derive(Debug)]
pub enum CTableParseError {
    /// The final token does not name a supported CTable format.
    UnknownFormat(String),
    /// A V2000 atom or bond count could not be parsed.
    InvalidV2000Counts(ParseIntError),
    /// A V3000 atom or bond count could not be parsed.
    InvalidV3000Counts,
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
    pub fn parse(lines: Vec<String>) -> Result<Self, CTableParseError> {
        let title = lines[0].trim().to_owned();
        let source = lines[1].clone();
        let comment = lines[2].clone();
        let first_counts = lines[3].clone();
        let format_token =
            parse_format(&first_counts).expect("expected a final CTable format token");
        let format = match format_token {
            "V2000" => CTableFormat::V2000,
            "V3000" => CTableFormat::V3000,
            other => return Err(CTableParseError::UnknownFormat(other.to_owned())),
        };
        let counts = match format {
            CTableFormat::V2000 => first_counts,
            CTableFormat::V3000 => lines[5].trim().to_owned(),
        };
        let (num_atoms, num_bonds) = match format {
            CTableFormat::V2000 => {
                parse_v2000_counts(&counts).map_err(CTableParseError::InvalidV2000Counts)?
            }
            CTableFormat::V3000 => {
                parse_v3000_counts(&counts).ok_or(CTableParseError::InvalidV3000Counts)?
            }
        };

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
    pub fn valid_atom_indices(&self, strict: bool) -> Result<(), &'static str> {
        if self.format != CTableFormat::V3000 {
            return Err("not_implemented");
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
    pub fn valid_bond_indices(&self, strict: bool) -> Result<(), &'static str> {
        if self.format != CTableFormat::V3000 {
            return Err("not_implemented");
        }
        validate_bond_indices(&self.lines, self.num_bonds as usize, strict)
    }

    /// Renumber atom and bond indices, replacing the stored lines on success.
    ///
    /// Only V3000 tables are supported. Atom and bond lines are numbered in
    /// their current order, bond atom references are remapped, and the
    /// rewritten counts line reflects the number of table lines. Return
    /// `"duplicate"` when a duplicate atom index makes a bond reference
    /// ambiguous, or `"not_implemented"` for V2000 input.
    pub fn renumber_indices(&mut self) -> Result<(), &'static str> {
        if self.format != CTableFormat::V3000 {
            return Err("not_implemented");
        }
        let Some(lines) = renumber_indices(&self.lines, &self.counts) else {
            return Err("duplicate");
        };
        self.lines = lines;
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

/// Parse atom and bond counts from the V3000 `COUNTS` line.
///
/// The atom and bond values are the fourth and fifth whitespace-delimited
/// fields. Return `None` if either field is absent or is not an integer.
pub fn parse_v3000_counts(line: &str) -> Option<(i64, i64)> {
    let mut fields = line
        .split(|character: char| {
            character.is_ascii_whitespace() || ('\u{1c}'..='\u{1f}').contains(&character)
        })
        .filter(|field| !field.is_empty());
    let atoms = fields.nth(3)?.parse().ok()?;
    let bonds = fields.next()?.parse().ok()?;
    Some((atoms, bonds))
}

/// Return the final whitespace-delimited token from a counts line.
pub fn parse_format(line: &str) -> Option<&str> {
    let mut fields = line
        .split(|character: char| {
            character.is_ascii_whitespace() || ('\u{1c}'..='\u{1f}').contains(&character)
        })
        .filter(|field| !field.is_empty());
    fields.next_back()
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
) -> Result<(), &'static str> {
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
) -> Result<(), &'static str> {
    let (start, end) = bond_line_range(lines);
    validate_indices(&lines[start..end], num_bonds, strict)
}

fn validate_indices(
    lines: &[String],
    expected_count: usize,
    strict: bool,
) -> Result<(), &'static str> {
    let mut seen_indices = HashSet::new();
    for (line_index, line) in lines.iter().enumerate() {
        let index = line
            .split_whitespace()
            .nth(2)
            .expect("index line has no index")
            .parse::<i64>()
            .expect("index is not an integer");
        if strict && (line_index + 1) as i64 != index {
            return Err("out_of_order");
        }
        if !seen_indices.insert(index) {
            return Err("duplicate");
        }
    }
    if seen_indices.len() < expected_count {
        return Err("fewer");
    }
    if seen_indices.len() > expected_count {
        return Err("more");
    }
    Ok(())
}

/// Return lines with atom and bond indices renumbered in table order.
///
/// Atom and bond lines are numbered from one in their existing order, and
/// bond atom references are remapped through the atom-index map. The counts
/// line is updated to reflect the number of raw atom and bond lines. Return
/// `None` when duplicate atom indices make a bond's atom mapping ambiguous.
pub fn renumber_indices(lines: &[String], counts: &str) -> Option<Vec<String>> {
    let prefix = "M  V30 ";
    let (atom_start, atom_end) = atom_line_range(lines);
    let atom_lines = &lines[atom_start..atom_end];
    let mut atom_index_mapping: HashMap<i64, Vec<usize>> = HashMap::new();
    let mut new_atom_lines = Vec::new();
    for (line_index, line) in atom_lines.iter().enumerate() {
        let fields: Vec<_> = line.split_whitespace().collect();
        let old_index = fields[2]
            .parse::<i64>()
            .expect("atom index is not an integer");
        let new_index = line_index + 1;
        atom_index_mapping
            .entry(old_index)
            .or_default()
            .push(new_index);
        let skip_chars = prefix.len() + fields[2].len();
        let mut new_line = format!("{prefix}{new_index}");
        new_line.push_str(
            line.get(skip_chars..)
                .expect("atom index offset is not a character boundary"),
        );
        new_atom_lines.push(new_line);
    }

    let (bond_start, bond_end) = bond_line_range(lines);
    let bond_lines = &lines[bond_start..bond_end];
    let mut new_bond_lines = Vec::new();
    for (line_index, line) in bond_lines.iter().enumerate() {
        let fields: Vec<_> = line.split_whitespace().collect();
        let from_index = fields[4]
            .parse::<i64>()
            .expect("bond atom index is not an integer");
        let to_index = fields[5]
            .parse::<i64>()
            .expect("bond atom index is not an integer");
        let from_mapping = atom_index_mapping
            .get(&from_index)
            .expect("bond references a missing atom index");
        let to_mapping = atom_index_mapping
            .get(&to_index)
            .expect("bond references a missing atom index");
        if from_mapping.len() > 1 || to_mapping.len() > 1 {
            return None;
        }
        let trailing = if line.ends_with('\n') { "\n" } else { "" };
        new_bond_lines.push(format!(
            "{prefix}{} {} {} {}{trailing}",
            line_index + 1,
            fields[3],
            from_mapping[0],
            to_mapping[0]
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
    Some(new_lines)
}
