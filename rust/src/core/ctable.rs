use super::whitespace;

pub fn integer(text: &str) -> Option<i128> {
    (text.len() <= 40)
        .then(|| text.trim().parse().ok())
        .flatten()
}

pub fn format(text: &str) -> Option<&str> {
    text.split(whitespace).rfind(|token| !token.is_empty())
}

pub fn v2000_fields(text: &str) -> (&str, &str) {
    let middle = text
        .char_indices()
        .nth(3)
        .map_or(text.len(), |(offset, _)| offset);
    let end = text
        .char_indices()
        .nth(6)
        .map_or(text.len(), |(offset, _)| offset);
    (&text[..middle], &text[middle..end])
}

pub fn parse_v2000_counts(text: &str) -> Option<(i128, i128)> {
    let (atoms, bonds) = v2000_fields(text);
    Some((integer(atoms)?, integer(bonds)?))
}

pub fn parse_v3000_counts(text: &str) -> Option<(i128, i128)> {
    let mut tokens = text.split(whitespace).filter(|token| !token.is_empty());
    Some((integer(tokens.nth(3)?)?, integer(tokens.next()?)?))
}

pub fn token(text: &str, position: usize) -> Option<&str> {
    text.split(whitespace)
        .filter(|token| !token.is_empty())
        .nth(position)
}

pub fn tokens(text: &str) -> Vec<&str> {
    text.split(whitespace)
        .filter(|token| !token.is_empty())
        .collect()
}

pub fn token_suffix(text: &str, position: usize) -> String {
    tokens(text)
        .into_iter()
        .skip(position)
        .collect::<Vec<_>>()
        .join(" ")
}

pub fn after_characters(text: &str, count: usize) -> &str {
    text.char_indices()
        .nth(count)
        .map_or("", |(offset, _)| &text[offset..])
}
