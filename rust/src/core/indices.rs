use super::ctable;
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

pub enum MappingFailure {
    Duplicate,
    Missing,
}

pub struct Mapping {
    values: HashMap<Index, Vec<usize>>,
    next_index: usize,
}

impl Mapping {
    pub fn new() -> Self {
        Self {
            values: HashMap::new(),
            next_index: 1,
        }
    }

    pub fn add_atom(&mut self, index: Index) -> usize {
        let new_index = self.next_index;
        self.next_index += 1;
        self.values.entry(index).or_default().push(new_index);
        new_index
    }

    pub fn remap_bond(&self, from: Index, to: Index) -> Result<(usize, usize), MappingFailure> {
        if self
            .values
            .get(&from)
            .is_some_and(|values| values.len() > 1)
            || self.values.get(&to).is_some_and(|values| values.len() > 1)
        {
            return Err(MappingFailure::Duplicate);
        }
        let new_from = self
            .values
            .get(&from)
            .and_then(|values| values.first())
            .ok_or(MappingFailure::Missing)?;
        let new_to = self
            .values
            .get(&to)
            .and_then(|values| values.first())
            .ok_or(MappingFailure::Missing)?;
        Ok((*new_from, *new_to))
    }
}

#[derive(Clone, Copy)]
pub struct LineFlags {
    pub counts: bool,
    pub begin_atom: bool,
    pub end_atom: bool,
    pub begin_bond: bool,
    pub end_bond: bool,
}

pub struct Assembler<T> {
    lines: Vec<T>,
    atoms: Vec<T>,
    bonds: Vec<T>,
    counts: T,
    appending: bool,
}

impl<T: Clone> Assembler<T> {
    pub fn new(atoms: Vec<T>, bonds: Vec<T>, counts: T) -> Self {
        Self {
            lines: Vec::new(),
            atoms,
            bonds,
            counts,
            appending: true,
        }
    }

    pub fn push(&mut self, line: T, flags: LineFlags) {
        if self.appending {
            self.lines.push(line.clone());
        }
        if flags.counts {
            self.lines.pop();
            self.lines.push(self.counts.clone());
        }
        if flags.begin_atom {
            self.appending = false;
        }
        if flags.end_atom {
            self.appending = true;
            self.lines.extend(self.atoms.iter().cloned());
            self.lines.push(line.clone());
        }
        if flags.begin_bond {
            self.appending = false;
        }
        if flags.end_bond {
            self.appending = true;
            self.lines.extend(self.bonds.iter().cloned());
            self.lines.push(line);
        }
    }

    pub fn finish(self) -> Vec<T> {
        self.lines
    }
}

pub fn atom_line(text: &str, new_index: usize) -> String {
    let token = ctable::token(text, 2).unwrap_or_default();
    let suffix = ctable::after_characters(text, 7 + token.chars().count());
    format!("M  V30 {new_index}{suffix}")
}

pub fn bond_line(new_index: usize, order: &str, from: usize, to: usize, newline: bool) -> String {
    let newline = if newline { "\n" } else { "" };
    format!("M  V30 {new_index} {order} {from} {to}{newline}")
}

pub fn counts_line(atom_count: usize, bond_count: usize, suffix: &str) -> String {
    format!("M  V30 COUNTS {atom_count} {bond_count} {suffix}")
}
