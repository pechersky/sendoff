use std::fmt;

/// An invalid SD data record encountered while parsing.
#[derive(Debug)]
pub enum RecordParseError {
    /// A field header did not contain a field name.
    MissingFieldName,
}

impl fmt::Display for RecordParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingFieldName => formatter.write_str("SD data field header is missing a name"),
        }
    }
}

pub fn append_record(record_name: &str, value: &str) -> (String, String, String) {
    (
        format!("> <{record_name}>\n"),
        format!("{value}\n"),
        "\n".to_owned(),
    )
}

pub fn records(metadata: &[String]) -> Result<Vec<(String, String)>, RecordParseError> {
    let mut records = Vec::new();
    let mut index = 0;
    while let Some(record) = next_record(metadata, &mut index)? {
        records.push(record);
    }
    Ok(records)
}

pub fn records_bytes(metadata: &[Vec<u8>]) -> Result<Vec<(Vec<u8>, Vec<u8>)>, RecordParseError> {
    let mut records = Vec::new();
    let mut index = 0;
    while let Some(record) = next_record_bytes(metadata, &mut index)? {
        records.push(record);
    }
    Ok(records)
}

pub fn next_record(
    metadata: &[String],
    index: &mut usize,
) -> Result<Option<(String, String)>, RecordParseError> {
    while *index < metadata.len() {
        let start = *index;
        while *index < metadata.len() && !metadata[*index].trim().is_empty() {
            *index += 1;
        }
        let end = *index;
        if *index < metadata.len() {
            *index += 1;
        }

        let Some(header) = metadata[start..end].first() else {
            continue;
        };
        let header = header.trim();
        if !header.starts_with("> ") {
            continue;
        }
        let Some((_, field)) = header.split_once("> ") else {
            continue;
        };
        let field = field.trim();
        let field = field.rsplit_once('>').map_or(field, |(name, _)| name);
        let Some((_, field)) = field.split_once('<') else {
            return Err(RecordParseError::MissingFieldName);
        };
        let field = field.to_owned();

        let value_lines = &metadata[start + 1..end];
        let value_capacity = value_lines
            .iter()
            .fold(0usize, |size, line| size.saturating_add(line.trim().len()))
            .saturating_add(value_lines.len().saturating_sub(1));
        let mut value = String::with_capacity(value_capacity);
        for (line_index, line) in value_lines.iter().enumerate() {
            if line_index > 0 {
                value.push('\n');
            }
            value.push_str(line.trim());
        }
        return Ok(Some((field, value)));
    }
    Ok(None)
}

pub fn next_record_bytes(
    metadata: &[Vec<u8>],
    index: &mut usize,
) -> Result<Option<(Vec<u8>, Vec<u8>)>, RecordParseError> {
    while *index < metadata.len() {
        let start = *index;
        while *index < metadata.len() && !metadata[*index].trim_ascii().is_empty() {
            *index += 1;
        }
        let end = *index;
        if *index < metadata.len() {
            *index += 1;
        }

        let Some(header) = metadata[start..end].first() else {
            continue;
        };
        let Some(header) = header.trim_ascii().strip_prefix(b"> ") else {
            continue;
        };
        let field = header.trim_ascii();
        let field = field
            .iter()
            .rposition(|byte| *byte == b'>')
            .map_or(field, |end| &field[..end]);
        let Some(name_start) = field.iter().position(|byte| *byte == b'<') else {
            return Err(RecordParseError::MissingFieldName);
        };
        let field = field[name_start + 1..].to_vec();

        let value_lines = &metadata[start + 1..end];
        let value_capacity = value_lines
            .iter()
            .fold(0usize, |size, line| {
                size.saturating_add(line.trim_ascii().len())
            })
            .saturating_add(value_lines.len().saturating_sub(1));
        let mut value = Vec::with_capacity(value_capacity);
        for (line_index, line) in value_lines.iter().enumerate() {
            if line_index > 0 {
                value.push(b'\n');
            }
            value.extend_from_slice(line.trim_ascii());
        }
        return Ok(Some((field, value)));
    }
    Ok(None)
}
