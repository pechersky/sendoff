use std::collections::HashSet;
use std::num::ParseIntError;

pub fn parse_v2000_counts(line: &str) -> Result<(i64, i64), ParseIntError> {
    let atom_field = line.get(..3).unwrap_or(line);
    let bond_field = line
        .get(3..6)
        .unwrap_or_else(|| line.get(3..).unwrap_or(""));
    Ok((atom_field.trim().parse()?, bond_field.trim().parse()?))
}

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

pub fn parse_format(line: &str) -> Option<&str> {
    let mut fields = line
        .split(|character: char| {
            character.is_ascii_whitespace() || ('\u{1c}'..='\u{1f}').contains(&character)
        })
        .filter(|field| !field.is_empty());
    fields.next_back()
}

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

pub fn validate_atom_indices(
    lines: &[String],
    num_atoms: usize,
    strict: bool,
) -> Result<(), &'static str> {
    let (start, end) = atom_line_range(lines);
    validate_indices(&lines[start..end], num_atoms, strict)
}

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
