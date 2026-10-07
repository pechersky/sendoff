use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::Range;

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

fn malformed(line: usize, problem: Malformed) -> CTableError {
    CTableError::Malformed { line, problem }
}

/// Parse `text`, the named field of CTable line `line`, as an integer.
fn integer(text: &str, line: usize, field: &'static str) -> Result<i64, CTableError> {
    text.trim()
        .parse()
        .map_err(|_| malformed(line, Malformed::NotInteger(field, text.to_owned())))
}

/// The most lines, title included, that `parse_header` reads.
pub const HEADER_LINES: usize = 6;

/// The index of the V3000 counts line, title included.
const V3000_COUNTS_LINE: usize = 5;

/// The index of the counts line, which ends with the format.
const COUNTS_LINE: usize = 3;

/// The line at `index`, or an error naming the first line the table lacks.
fn line_at<'a>(lines: &[&'a str], index: usize) -> Result<&'a str, CTableError> {
    lines
        .get(index)
        .copied()
        .ok_or_else(|| malformed(lines.len() + 1, Malformed::Missing(lines.len())))
}

/// Parse the format and the declared atom and bond counts.
///
/// `lines` starts with the title and needs at most six lines. V3000 input
/// skips the CTAB BEGIN line and reads its counts from the following line.
pub fn parse_header(lines: &[&str]) -> Result<(CTableFormat, i64, i64), CTableError> {
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
    Ok((format, atoms, bonds))
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
        let (format, num_atoms, num_bonds) =
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

    /// Validate the root CTAB's atom or bond lines against the counts line.
    ///
    /// Each line's index must be a unique integer and, with `strict`, its
    /// one-based position; the section must have as many lines as the
    /// counts line declares.
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

    /// Set the counts line's atom and bond counts to the number of atom and
    /// bond lines, and the parsed counts with them.
    ///
    /// Only V3000 tables are supported. The rest of the counts line, spacing
    /// and ending included, stays as it was.
    pub fn fix_counts(&mut self) -> Result<(), CTableError> {
        self.require_v3000("fixing counts")?;
        let (atoms, bonds) = sections("fixing counts", &self.lines)?;
        let (atoms, bonds) = (atoms.len(), bonds.len());
        let line = &self.lines[V3000_COUNTS_LINE];
        let (atom_field, bond_field) = v3000_count_fields(line)?;
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
    let (atoms, bonds) = (
        section_range(lines, Section::Atom)?,
        section_range(lines, Section::Bond)?,
    );
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

/// Return the range of the root CTAB's lines between `section`'s begin and
/// end lines, empty when there is no begin line.
///
/// The root CTAB ends at the first `END CTAB` line; CTABs nested in RGROUP or
/// TEMPLATE blocks follow it.
fn section_range(lines: &[String], section: Section) -> Result<Range<usize>, CTableError> {
    let (begin, end) = section.markers();
    let root = lines
        .iter()
        .position(|line| starts_with_keyword(line, "M  V30 END CTAB"))
        .map_or(lines, |end| &lines[..end]);
    let mut begins = root
        .iter()
        .enumerate()
        .filter(|(_, line)| starts_with_keyword(line, begin))
        .map(|(index, _)| index);
    let Some(index) = begins.next() else {
        return Ok(0..0);
    };
    if let Some(repeat) = begins.next() {
        return Err(malformed(repeat + 1, Malformed::Repeated(begin)));
    }
    let start = index + 1;
    let len = root[start..]
        .iter()
        .position(|line| starts_with_keyword(line, end))
        .ok_or_else(|| malformed(start, Malformed::Unterminated { begin, end }))?;
    Ok(start..start + len)
}

/// The index field of `line`, CTable line `at`: its third field.
fn index_field(line: &str, at: usize, section: Section) -> Result<&str, CTableError> {
    line.split_whitespace()
        .nth(2)
        .ok_or_else(|| malformed(at, Malformed::NoField(section.index_field())))
}

/// Validate the count and unique indices of the `section` lines in `range`
/// and, when strict, their order.
fn validate_indices(
    lines: &[String],
    range: Range<usize>,
    section: Section,
    declared: i64,
    strict: bool,
) -> Result<(), CTableError> {
    let mut seen_indices = HashSet::new();
    for (position, line_index) in range.clone().enumerate() {
        let at = line_index + 1;
        let field = index_field(&lines[line_index], at, section)?;
        let index = integer(field, at, section.index_field())?;
        if strict && usize::try_from(index) != Ok(position + 1) {
            return Err(CTableError::OutOfOrder {
                section,
                line: at,
                index,
            });
        }
        if !seen_indices.insert(index) {
            return Err(CTableError::Duplicate {
                section,
                line: at,
                index,
            });
        }
    }
    let found = range.len();
    match (found as i64).cmp(&declared) {
        std::cmp::Ordering::Less => Err(CTableError::Fewer {
            section,
            found,
            declared,
        }),
        std::cmp::Ordering::Greater => Err(CTableError::More {
            section,
            found,
            declared,
        }),
        std::cmp::Ordering::Equal => Ok(()),
    }
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
/// ambiguous, and with `NoAtom` when a bond references an absent atom.
fn renumbered(
    lines: &[String],
    atoms: Range<usize>,
    bonds: Range<usize>,
) -> Result<(Vec<String>, Vec<String>), CTableError> {
    let prefix = "M  V30 ";
    let (atom_lines, bond_lines) = (&lines[atoms.clone()], &lines[bonds.clone()]);
    let mut changed = false;
    let mut atom_index_mapping = HashMap::with_capacity(atom_lines.len());
    let mut new_atom_lines = Vec::with_capacity(atom_lines.len());
    for (line_index, line) in atom_lines.iter().enumerate() {
        let at = atoms.start + line_index + 1;
        let old_field = index_field(line, at, Section::Atom)?;
        let old_index = integer(old_field, at, "atom index")?;
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
        let at = bonds.start + line_index + 1;
        let mut fields = line.split_whitespace().skip(2);
        let mut field = |name| {
            let text = fields
                .next()
                .ok_or_else(|| malformed(at, Malformed::NoField(name)))?;
            Ok::<_, CTableError>(text)
        };
        let index = field("bond index")?;
        let bond_type = field("bond type")?;
        let (from_field, to_field) = (field("first bond atom")?, field("second bond atom")?);
        changed |= usize::try_from(integer(index, at, "bond index")?) != Ok(line_index + 1);
        let from_index = integer(from_field, at, "first bond atom")?;
        let to_index = integer(to_field, at, "second bond atom")?;
        let tail = &line[end_of(line, to_field)..];
        let new_index = |index| match atom_index_mapping.get(&index) {
            Some(&NewIndex::Unique(new)) => Ok(Some(new)),
            Some(NewIndex::Duplicate) => Err(CTableError::AmbiguousAtom { line: at, index }),
            None => Ok(None),
        };
        // A duplicated atom index outranks a missing one in the same bond.
        let (new_from, new_to) = match (new_index(from_index)?, new_index(to_index)?) {
            (Some(from), Some(to)) => (from, to),
            (None, _) => return Err(malformed(at, Malformed::NoAtom(from_index))),
            (_, None) => return Err(malformed(at, Malformed::NoAtom(to_index))),
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
