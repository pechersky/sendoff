//! The root CTAB's V3000 atom and bond sections: finding them, reading
//! their logical lines, and validating and renumbering their indices.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::ops::Range;

use crate::core::ctable::error::{CTableError, Malformed, integer, malformed};
use crate::core::ctable::{Section, end_of};

/// Whether `line` starts with `prefix`, ignoring ASCII case as V3000 does.
fn starts_with_keyword(line: &str, prefix: &str) -> bool {
    line.as_bytes()
        .get(..prefix.len())
        .is_some_and(|start| start.eq_ignore_ascii_case(prefix.as_bytes()))
}

/// Find the root CTAB's atom and bond sections for `operation`, which edits
/// a table line by line. Refuse CTABs nested in RGROUP or TEMPLATE blocks,
/// whose own atom and bond lines it would leave as they are.
pub(super) fn sections(
    operation: &'static str,
    lines: &[String],
) -> Result<(Range<usize>, Range<usize>), CTableError> {
    let ctabs = lines
        .iter()
        .filter(|line| starts_with_keyword(line, "M  V30 BEGIN CTAB"))
        .count();
    if ctabs > 1 {
        return Err(CTableError::Unsupported {
            operation,
            construct: "CTABs nested in RGROUP or TEMPLATE blocks",
        });
    }
    Ok((
        section_range(lines, Section::Atom)?,
        section_range(lines, Section::Bond)?,
    ))
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

/// Return the range of the root CTAB's lines between `section`'s begin and
/// end lines, empty when there is no begin line.
///
/// The root CTAB ends at the first `END CTAB` line; CTABs nested in RGROUP or
/// TEMPLATE blocks follow it.
pub(super) fn section_range(
    lines: &[String],
    section: Section,
) -> Result<Range<usize>, CTableError> {
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
    let range = start..start + len;
    if let Some(last) = range
        .clone()
        .last()
        .filter(|&last| continued(&lines[last]).is_some())
    {
        return Err(malformed(last + 1, Malformed::Dangling(end)));
    }
    Ok(range)
}

/// The text of `line` before a terminal `-` that continues it onto the next
/// line, if it has one.
pub(super) fn continued(line: &str) -> Option<&str> {
    line.trim_end_matches(['\r', '\n']).strip_suffix('-')
}

/// The logical lines of the lines in `range`, each with the CTable line
/// number of its first line.
///
/// A line ending in `-` continues onto the next: the `-` and the next line's
/// `M  V30 ` marker are dropped and the two joined. `section_range` ensures a
/// section's last line does not continue.
pub(super) fn logical_lines(
    lines: &[String],
    range: Range<usize>,
) -> impl Iterator<Item = (usize, Cow<'_, str>)> {
    let mut index = range.start;
    std::iter::from_fn(move || {
        if index >= range.end {
            return None;
        }
        let first = index;
        let mut line = Cow::Borrowed(lines[first].as_str());
        index += 1;
        while index < range.end {
            let Some(head) = continued(&line).map(str::len) else {
                break;
            };
            let next = &lines[index];
            let payload = if starts_with_keyword(next, "M  V30 ") {
                &next["M  V30 ".len()..]
            } else {
                next
            };
            let joined = line.to_mut();
            joined.truncate(head);
            joined.push_str(payload);
            index += 1;
        }
        Some((first + 1, line))
    })
}

/// The index field of `line`, CTable line `at`: its third field.
fn index_field(line: &str, at: usize, section: Section) -> Result<&str, CTableError> {
    line.split_whitespace()
        .nth(2)
        .ok_or_else(|| malformed(at, Malformed::NoField(section.index_field())))
}

/// Validate the count and unique indices of the `section` lines in `range`,
/// read as logical lines, and, when strict, their order.
pub(super) fn validate_indices(
    lines: &[String],
    range: Range<usize>,
    section: Section,
    declared: i64,
    strict: bool,
) -> Result<(), CTableError> {
    let mut seen_indices = HashSet::new();
    let mut found = 0;
    for (position, (at, line)) in logical_lines(lines, range).enumerate() {
        found = position + 1;
        let field = index_field(&line, at, section)?;
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
pub(super) fn renumbered(
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
        let rest = &line[end_of(line, old_field)..];
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
