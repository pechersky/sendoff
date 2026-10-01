use std::collections::{HashMap, HashSet};
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
