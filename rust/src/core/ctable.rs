use std::collections::{HashMap, HashSet};
use std::num::ParseIntError;
use std::ops::Range;

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
    /// The operation would corrupt or miscount the named construct, so it
    /// refuses.
    Unsupported {
        operation: &'static str,
        construct: &'static str,
    },
    /// A section has no line with this end marker.
    MissingLine(&'static str),
    /// The root CTAB has more than one line with this begin marker.
    RepeatedLine(&'static str),
    /// Strict validation found an atom or bond index, as named, out of
    /// one-based order.
    OutOfOrder(&'static str),
    /// An atom or bond index, as named, occurs twice.
    Duplicate(&'static str),
    /// There are fewer unique atom or bond indices, as named, than the counts
    /// line declares.
    Fewer(&'static str),
    /// There are more unique atom or bond indices, as named, than the counts
    /// line declares.
    More(&'static str),
    /// A bond references an atom index that several atoms share, so
    /// renumbering cannot tell which atom it means.
    AmbiguousAtom,
}

impl From<ParseIntError> for CTableError {
    fn from(error: ParseIntError) -> Self {
        Self::InvalidInteger(error.to_string())
    }
}

/// The most lines, title included, that `parse_header` reads.
pub const HEADER_LINES: usize = 6;

/// The index of the V3000 counts line, title included.
const V3000_COUNTS_LINE: usize = 5;

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
            let counts = line_at(lines, V3000_COUNTS_LINE)?.trim();
            (counts, parse_v3000_counts(counts)?)
        }
    };
    Ok((format, counts, atoms, bonds))
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
        let header: [&str; HEADER_LINES] =
            std::array::from_fn(|index| lines.get(index).map_or("", String::as_str));
        let (format, _, num_atoms, num_bonds) =
            parse_header(&header[..lines.len().min(HEADER_LINES)])?;
        let title = lines[0].trim().to_owned();
        let source = lines[1].clone();
        let comment = lines[2].clone();
        Ok(Self {
            lines,
            title,
            source,
            comment,
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
    /// count line. The atom lines are those of the root CTAB's ATOM section.
    pub fn valid_atom_indices(&self, strict: bool) -> Result<(), CTableError> {
        if self.format != CTableFormat::V3000 {
            return Err(CTableError::NotImplemented);
        }
        let atoms = atom_section(&self.lines)?;
        validate_indices(&self.lines[atoms], self.num_atoms as usize, strict, "atom")
    }

    /// Validate raw bond lines against the parsed bond count.
    ///
    /// With `strict` enabled, indices must be one-based and in increasing
    /// order. The result distinguishes unsupported V2000 validation,
    /// out-of-order and duplicate indices, and fewer or more lines than the
    /// count line. The bond lines are those of the root CTAB's BOND section.
    pub fn valid_bond_indices(&self, strict: bool) -> Result<(), CTableError> {
        if self.format != CTableFormat::V3000 {
            return Err(CTableError::NotImplemented);
        }
        let bonds = bond_section(&self.lines)?;
        validate_indices(&self.lines[bonds], self.num_bonds as usize, strict, "bond")
    }

    /// Renumber atom and bond indices, replacing the stored lines on success.
    ///
    /// Only V3000 tables are supported. Atom and bond lines are numbered in
    /// their current order and bond atom references are remapped; the counts
    /// line and parsed counts are left as they are (see `fix_counts`). Lines
    /// are only replaced on success.
    pub fn renumber_indices(&mut self) -> Result<(), CTableError> {
        if self.format != CTableFormat::V3000 {
            return Err(CTableError::NotImplemented);
        }
        let (atoms, bonds) = sections("renumbering", &self.lines)?;
        let (mut new_atoms, mut new_bonds) = renumbered(&self.lines, atoms.clone(), bonds.clone())?;
        self.lines[atoms].swap_with_slice(&mut new_atoms);
        self.lines[bonds].swap_with_slice(&mut new_bonds);
        Ok(())
    }

    /// Set the counts line's atom and bond counts to the number of atom and
    /// bond lines, and the parsed counts with them.
    ///
    /// Only V3000 tables are supported. The rest of the counts line, spacing
    /// and ending included, stays as it was.
    pub fn fix_counts(&mut self) -> Result<(), CTableError> {
        if self.format != CTableFormat::V3000 {
            return Err(CTableError::NotImplemented);
        }
        let (atoms, bonds) = sections("fixing counts", &self.lines)?;
        let (atoms, bonds) = (atoms.len(), bonds.len());
        let line = &self.lines[V3000_COUNTS_LINE];
        let (atom_field, bond_field) = {
            let mut counts = fields(line).skip(3);
            (
                counts
                    .next()
                    .ok_or(CTableError::MissingField("atom count"))?,
                counts
                    .next()
                    .ok_or(CTableError::MissingField("bond count"))?,
            )
        };
        let fixed = format!(
            "{}{atoms}{}{bonds}{}",
            &line[..start_of(line, atom_field)],
            &line[end_of(line, atom_field)..start_of(line, bond_field)],
            &line[end_of(line, bond_field)..],
        );
        self.lines[V3000_COUNTS_LINE] = fixed;
        self.num_atoms = atoms as i64;
        self.num_bonds = bonds as i64;
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

/// Whether `line` starts with `prefix`, ignoring ASCII case as V3000 does.
fn starts_with_keyword(line: &str, prefix: &str) -> bool {
    line.as_bytes()
        .get(..prefix.len())
        .is_some_and(|start| start.eq_ignore_ascii_case(prefix.as_bytes()))
}

/// Find the root CTAB's atom and bond sections for `operation`, which edits
/// a table line by line. Refuse tables it would get wrong: CTABs nested in
/// RGROUP or TEMPLATE blocks, whose own atom and bond lines it would leave as
/// they are, and `-` continued atom or bond lines, which it would take for
/// atoms or bonds.
fn sections(
    operation: &'static str,
    lines: &[String],
) -> Result<(Range<usize>, Range<usize>), CTableError> {
    let unsupported = |construct| {
        Err(CTableError::Unsupported {
            operation,
            construct,
        })
    };
    let ctabs = lines
        .iter()
        .filter(|line| starts_with_keyword(line, "M  V30 BEGIN CTAB"))
        .count();
    if ctabs > 1 {
        return unsupported("CTABs nested in RGROUP or TEMPLATE blocks");
    }
    let (atoms, bonds) = (atom_section(lines)?, bond_section(lines)?);
    if lines[atoms.clone()]
        .iter()
        .chain(&lines[bonds.clone()])
        .any(|line| line.trim_end_matches(['\r', '\n']).ends_with('-'))
    {
        return unsupported("continued (-) atom or bond lines");
    }
    Ok((atoms, bonds))
}

/// Refuse to change indices that other data refers to, which renumbering
/// does not remap: atom ATTCHORD and bond ENDPTS properties, and Sgroup,
/// collection, 3D object, and link node lines.
fn check_references(
    lines: &[String],
    atom_lines: &[String],
    bond_lines: &[String],
) -> Result<(), CTableError> {
    let has_property = |lines: &[String], keyword: &str| {
        lines.iter().any(|line| {
            line.split_whitespace()
                .any(|field| starts_with_keyword(field, keyword))
        })
    };
    let refuse = |construct| {
        Err(CTableError::Unsupported {
            operation: "renumbering",
            construct,
        })
    };
    if has_property(atom_lines, "ATTCHORD=") {
        return refuse("ATTCHORD atom references");
    }
    if has_property(bond_lines, "ENDPTS=") {
        return refuse("ENDPTS atom references");
    }
    for (marker, construct) in [
        ("M  V30 BEGIN SGROUP", "SGROUP atom and bond references"),
        (
            "M  V30 BEGIN COLLECTION",
            "COLLECTION atom and bond references",
        ),
        ("M  V30 BEGIN OBJ3D", "OBJ3D atom references"),
        ("M  V30 LINKNODE", "LINKNODE atom references"),
    ] {
        if lines.iter().any(|line| starts_with_keyword(line, marker)) {
            return refuse(construct);
        }
    }
    Ok(())
}

/// The offset in `line` of `field`, a subslice of it.
fn start_of(line: &str, field: &str) -> usize {
    field.as_ptr() as usize - line.as_ptr() as usize
}

/// The offset in `line` just past `field`, a subslice of it.
fn end_of(line: &str, field: &str) -> usize {
    start_of(line, field) + field.len()
}

/// Return the range of root CTAB lines between the `begin` and `end` marker
/// lines, empty when there is no `begin` line.
///
/// The root CTAB ends at the first `END CTAB` line; CTABs nested in RGROUP or
/// TEMPLATE blocks follow it.
fn section(
    lines: &[String],
    begin: &'static str,
    end: &'static str,
) -> Result<Range<usize>, CTableError> {
    let root = lines
        .iter()
        .position(|line| starts_with_keyword(line, "M  V30 END CTAB"))
        .map_or(lines, |end| &lines[..end]);
    let mut begins = root
        .iter()
        .enumerate()
        .filter(|(_, line)| starts_with_keyword(line, begin));
    let Some((index, _)) = begins.next() else {
        return Ok(0..0);
    };
    if begins.next().is_some() {
        return Err(CTableError::RepeatedLine(begin));
    }
    let start = index + 1;
    let len = root[start..]
        .iter()
        .position(|line| starts_with_keyword(line, end))
        .ok_or(CTableError::MissingLine(end))?;
    Ok(start..start + len)
}

/// The range of the root CTAB's atom lines.
fn atom_section(lines: &[String]) -> Result<Range<usize>, CTableError> {
    section(lines, "M  V30 BEGIN ATOM", "M  V30 END ATOM")
}

/// The range of the root CTAB's bond lines.
fn bond_section(lines: &[String]) -> Result<Range<usize>, CTableError> {
    section(lines, "M  V30 BEGIN BOND", "M  V30 END BOND")
}

/// Validate a section's count and unique indices and, when strict, their
/// order. Each line's third whitespace-delimited token is its index.
fn validate_indices(
    lines: &[String],
    expected_count: usize,
    strict: bool,
    noun: &'static str,
) -> Result<(), CTableError> {
    let mut seen_indices = HashSet::new();
    for (line_index, line) in lines.iter().enumerate() {
        let index = integer(
            line.split_whitespace()
                .nth(2)
                .ok_or(CTableError::MissingField("index"))?,
        )?;
        if strict && (line_index + 1) as i64 != index {
            return Err(CTableError::OutOfOrder(noun));
        }
        if !seen_indices.insert(index) {
            return Err(CTableError::Duplicate(noun));
        }
    }
    if seen_indices.len() < expected_count {
        return Err(CTableError::Fewer(noun));
    }
    if seen_indices.len() > expected_count {
        return Err(CTableError::More(noun));
    }
    Ok(())
}

/// The new index of an old atom index, unless several atoms share it.
#[derive(Clone, Copy)]
enum NewIndex {
    Unique(usize),
    Duplicate,
}

/// Renumber the atom lines in `atoms` and the bond lines in `bonds`, returning
/// the new lines for each range.
///
/// Atom and bond lines are numbered from one in their existing order, and
/// bond atom references are remapped through the atom-index map; the rest of
/// each line, such as bond `CFG` and the ending, stays as it was. Fail with
/// `AmbiguousAtom` when duplicate atom indices make a bond's atom mapping
/// ambiguous, and with `MissingField` when a bond references an absent atom.
fn renumbered(
    lines: &[String],
    atoms: Range<usize>,
    bonds: Range<usize>,
) -> Result<(Vec<String>, Vec<String>), CTableError> {
    let prefix = "M  V30 ";
    let (atom_lines, bond_lines) = (&lines[atoms], &lines[bonds]);
    let mut changed = false;
    let mut atom_index_mapping = HashMap::with_capacity(atom_lines.len());
    let mut new_atom_lines = Vec::with_capacity(atom_lines.len());
    for (line_index, line) in atom_lines.iter().enumerate() {
        let old_field = line
            .split_whitespace()
            .nth(2)
            .ok_or(CTableError::MissingField("atom index"))?;
        let old_index = integer(old_field)?;
        let new_index = line_index + 1;
        changed |= usize::try_from(old_index) != Ok(new_index);
        atom_index_mapping
            .entry(old_index)
            .and_modify(|index| *index = NewIndex::Duplicate)
            .or_insert(NewIndex::Unique(new_index));
        // Like the Python slice line[n:], skip characters, not bytes.
        let skip_chars = prefix.chars().count() + old_field.chars().count();
        let rest = line
            .char_indices()
            .nth(skip_chars)
            .map_or("", |(offset, _)| &line[offset..]);
        new_atom_lines.push(format!("{prefix}{new_index}{rest}"));
    }

    let mut new_bond_lines = Vec::with_capacity(bond_lines.len());
    for (line_index, line) in bond_lines.iter().enumerate() {
        let missing = || CTableError::MissingField("bond atom index");
        let mut fields = line.split_whitespace().skip(2);
        changed |= fields.next().and_then(|field| field.parse().ok()) != Some(line_index + 1);
        let (Some(bond_type), Some(from_field)) = (fields.next(), fields.next()) else {
            return Err(missing());
        };
        let from_index = integer(from_field)?;
        let to_field = fields.next().ok_or_else(missing)?;
        let to_index = integer(to_field)?;
        let tail = &line[end_of(line, to_field)..];
        let (new_from, new_to) = match (
            atom_index_mapping.get(&from_index),
            atom_index_mapping.get(&to_index),
        ) {
            (Some(NewIndex::Duplicate), _) | (_, Some(NewIndex::Duplicate)) => {
                return Err(CTableError::AmbiguousAtom);
            }
            (Some(&NewIndex::Unique(from)), Some(&NewIndex::Unique(to))) => (from, to),
            _ => return Err(CTableError::MissingField("bonded atom")),
        };
        new_bond_lines.push(format!(
            "{prefix}{} {bond_type} {new_from} {new_to}{tail}",
            line_index + 1,
        ));
    }

    if changed {
        check_references(lines, atom_lines, bond_lines)?;
    }
    Ok((new_atom_lines, new_bond_lines))
}
