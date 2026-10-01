use super::whitespace;

pub fn record_name(text: &str) -> Option<&str> {
    let (_, remainder) = text.split_once("> ")?;
    let remainder = remainder.trim_matches(whitespace);
    let first = remainder
        .rsplit_once('>')
        .map_or(remainder, |(first, _)| first);
    first.split_once('<').map(|(_, name)| name)
}

pub fn header(name: &str) -> String {
    format!("> <{name}>\n")
}

pub fn value_line(value: &str) -> String {
    format!("{value}\n")
}
