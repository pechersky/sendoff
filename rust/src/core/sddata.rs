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

pub(crate) trait RecordValue: Sized {
    fn is_blank(&self) -> bool;
    fn record_name(&self) -> Result<Option<Self>, RecordParseError>;
    fn trimmed_len(&self) -> usize;
    fn append_trimmed(&self, output: &mut Self);
    fn with_capacity(capacity: usize) -> Self;
    fn push_line_separator(output: &mut Self);
}

impl RecordValue for String {
    fn is_blank(&self) -> bool {
        self.trim().is_empty()
    }

    fn record_name(&self) -> Result<Option<Self>, RecordParseError> {
        let Some(header) = self.trim().strip_prefix("> ") else {
            return Ok(None);
        };
        let field = header.trim();
        let field = field.rsplit_once('>').map_or(field, |(name, _)| name);
        let Some((_, name)) = field.split_once('<') else {
            return Err(RecordParseError::MissingFieldName);
        };
        Ok(Some(name.to_owned()))
    }

    fn trimmed_len(&self) -> usize {
        self.trim().len()
    }

    fn append_trimmed(&self, output: &mut Self) {
        output.push_str(self.trim());
    }

    fn with_capacity(capacity: usize) -> Self {
        String::with_capacity(capacity)
    }

    fn push_line_separator(output: &mut Self) {
        output.push('\n');
    }
}

impl RecordValue for Vec<u8> {
    fn is_blank(&self) -> bool {
        self.trim_ascii().is_empty()
    }

    fn record_name(&self) -> Result<Option<Self>, RecordParseError> {
        let Some(header) = self.trim_ascii().strip_prefix(b"> ") else {
            return Ok(None);
        };
        let field = header.trim_ascii();
        let field = field
            .iter()
            .rposition(|byte| *byte == b'>')
            .map_or(field, |end| &field[..end]);
        let Some(name_start) = field.iter().position(|byte| *byte == b'<') else {
            return Err(RecordParseError::MissingFieldName);
        };
        Ok(Some(field[name_start + 1..].to_vec()))
    }

    fn trimmed_len(&self) -> usize {
        self.trim_ascii().len()
    }

    fn append_trimmed(&self, output: &mut Self) {
        output.extend_from_slice(self.trim_ascii());
    }

    fn with_capacity(capacity: usize) -> Self {
        Vec::with_capacity(capacity)
    }

    fn push_line_separator(output: &mut Self) {
        output.push(b'\n');
    }
}

/// The name/value pairs parsed from one block's metadata lines.
pub type ParsedRecords<T> = Result<Vec<(T, T)>, RecordParseError>;

pub(crate) fn records<T: RecordValue>(metadata: &[T]) -> ParsedRecords<T> {
    let mut records = Vec::new();
    let mut index = 0;
    while let Some(record) = next_record(metadata, &mut index)? {
        records.push(record);
    }
    Ok(records)
}

pub(crate) fn next_record<T: RecordValue>(
    metadata: &[T],
    index: &mut usize,
) -> Result<Option<(T, T)>, RecordParseError> {
    while *index < metadata.len() {
        let start = *index;
        while *index < metadata.len() && !metadata[*index].is_blank() {
            *index += 1;
        }
        let end = *index;
        if *index < metadata.len() {
            *index += 1;
        }

        let Some(header) = metadata[start..end].first() else {
            continue;
        };
        let Some(field) = header.record_name()? else {
            continue;
        };

        let value_lines = &metadata[start + 1..end];
        let value_capacity = value_lines
            .iter()
            .fold(0usize, |size, line| size.saturating_add(line.trimmed_len()))
            .saturating_add(value_lines.len().saturating_sub(1));
        let mut value = T::with_capacity(value_capacity);
        for (line_index, line) in value_lines.iter().enumerate() {
            if line_index > 0 {
                T::push_line_separator(&mut value);
            }
            line.append_trimmed(&mut value);
        }
        return Ok(Some((field, value)));
    }
    Ok(None)
}
