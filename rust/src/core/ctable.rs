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
