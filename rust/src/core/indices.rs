use super::ctable;
use super::whitespace;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Eq, Hash, PartialEq)]
pub enum Index {
    Small(i128),
    Large(Vec<u8>),
}

pub fn parse_index(text: &str) -> Option<Index> {
    ctable::integer(text).map(Index::Small)
}

pub enum ValidationFailure {
    OutOfOrder,
    Duplicate,
    Fewer,
    More,
}

pub struct Validator {
    seen: HashSet<Index>,
    position: usize,
}

impl Validator {
    pub fn new() -> Self {
        Self {
            seen: HashSet::new(),
            position: 0,
        }
    }

    pub fn accept(&mut self, index: Index, strict: bool) -> Result<(), ValidationFailure> {
        if strict && index != Index::Small((self.position + 1) as i128) {
            return Err(ValidationFailure::OutOfOrder);
        }
        if !self.seen.insert(index) {
            return Err(ValidationFailure::Duplicate);
        }
        self.position += 1;
        Ok(())
    }

    pub fn finish(self, expected: i128) -> Result<(), ValidationFailure> {
        let actual = self.seen.len() as i128;
        if actual < expected {
            Err(ValidationFailure::Fewer)
        } else if actual > expected {
            Err(ValidationFailure::More)
        } else {
            Ok(())
        }
    }

    pub fn len(&self) -> usize {
        self.seen.len()
    }
}

pub enum RenumberFailure {
    InvalidAtomIndex,
    InvalidBondIndex,
    DuplicateMapping,
    MissingEndpoint,
    EmptyCounts,
}

pub struct Renumbered {
    pub lines: Vec<String>,
}

pub fn renumber(
    all_lines: &[String],
    atom_lines: &[String],
    bond_lines: &[String],
    counts: &str,
) -> Result<Renumbered, RenumberFailure> {
    let mut atoms = Vec::with_capacity(atom_lines.len());
    let mut mapping: HashMap<Index, Vec<usize>> = HashMap::new();
    for (position, line) in atom_lines.iter().enumerate() {
        let token = ctable::token(line, 2).ok_or(RenumberFailure::InvalidAtomIndex)?;
        let index = parse_index(token).ok_or(RenumberFailure::InvalidAtomIndex)?;
        let new_index = position + 1;
        mapping.entry(index).or_default().push(new_index);
        let suffix = ctable::after_characters(line, 7 + token.chars().count());
        atoms.push(format!("M  V30 {new_index}{suffix}"));
    }

    let mut bonds = Vec::with_capacity(bond_lines.len());
    for (position, line) in bond_lines.iter().enumerate() {
        let order = ctable::token(line, 3).ok_or(RenumberFailure::InvalidBondIndex)?;
        let from = ctable::token(line, 4)
            .and_then(parse_index)
            .ok_or(RenumberFailure::InvalidBondIndex)?;
        let to = ctable::token(line, 5)
            .and_then(parse_index)
            .ok_or(RenumberFailure::InvalidBondIndex)?;
        if mapping.get(&from).is_some_and(|values| values.len() > 1)
            || mapping.get(&to).is_some_and(|values| values.len() > 1)
        {
            return Err(RenumberFailure::DuplicateMapping);
        }
        let new_from = mapping
            .get(&from)
            .and_then(|values| values.first())
            .ok_or(RenumberFailure::MissingEndpoint)?;
        let new_to = mapping
            .get(&to)
            .and_then(|values| values.first())
            .ok_or(RenumberFailure::MissingEndpoint)?;
        let newline = if line.ends_with('\n') { "\n" } else { "" };
        bonds.push(format!(
            "M  V30 {} {order} {new_from} {new_to}{newline}",
            position + 1
        ));
    }

    if counts.is_empty() {
        return Err(RenumberFailure::EmptyCounts);
    }
    let suffix = ctable::token(counts, 5)
        .map(|_| {
            counts
                .split(whitespace)
                .filter(|token| !token.is_empty())
                .skip(5)
                .collect::<Vec<_>>()
                .join(" ")
        })
        .unwrap_or_default();
    let new_counts = format!("M  V30 COUNTS {} {} {suffix}", atoms.len(), bonds.len());
    let mut lines = Vec::with_capacity(
        all_lines.len() + atoms.len().saturating_sub(1) + bonds.len().saturating_sub(1),
    );
    let mut appending = true;
    for line in all_lines {
        if appending {
            lines.push(line.clone());
        }
        if line.starts_with("M  V30 COUNTS") {
            lines.pop();
            lines.push(new_counts.clone());
        }
        if line.starts_with("M  V30 BEGIN ATOM") {
            appending = false;
        }
        if line.starts_with("M  V30 END ATOM") {
            appending = true;
            lines.extend(atoms.iter().cloned());
            lines.push(line.clone());
        }
        if line.starts_with("M  V30 BEGIN BOND") {
            appending = false;
        }
        if line.starts_with("M  V30 END BOND") {
            appending = true;
            lines.extend(bonds.iter().cloned());
            lines.push(line.clone());
        }
    }
    Ok(Renumbered { lines })
}
