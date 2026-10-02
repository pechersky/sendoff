pub fn append_record(record_name: &str, value: &str) -> (String, String, String) {
    (
        format!("> <{record_name}>\n"),
        format!("{value}\n"),
        "\n".to_owned(),
    )
}

pub fn records(metadata: &[String]) -> Vec<(String, String)> {
    let mut records = Vec::new();
    let mut chunk = Vec::new();

    for line in metadata.iter().map(|line| line.trim()) {
        if line.is_empty() {
            append_record_chunk(&mut records, &chunk);
            chunk.clear();
        } else {
            chunk.push(line);
        }
    }
    append_record_chunk(&mut records, &chunk);

    records
}

pub fn records_bytes(metadata: &[Vec<u8>]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut records = Vec::new();
    let mut chunk = Vec::new();

    for line in metadata.iter().map(|line| line.trim_ascii()) {
        if line.is_empty() {
            append_record_chunk_bytes(&mut records, &chunk);
            chunk.clear();
        } else {
            chunk.push(line);
        }
    }
    append_record_chunk_bytes(&mut records, &chunk);

    records
}

fn append_record_chunk(records: &mut Vec<(String, String)>, chunk: &[&str]) {
    let Some((record_line, value_lines)) = chunk.split_first() else {
        return;
    };
    if !record_line.starts_with("> ") {
        return;
    }

    let record_name = record_line.split_once("> ").unwrap().1.trim();
    let record_name = record_name
        .rsplit_once('>')
        .map_or(record_name, |(name, _)| name);
    let record_name = record_name.split_once('<').unwrap().1.to_owned();
    records.push((record_name, value_lines.join("\n")));
}

fn append_record_chunk_bytes(records: &mut Vec<(Vec<u8>, Vec<u8>)>, chunk: &[&[u8]]) {
    let Some((record_line, value_lines)) = chunk.split_first() else {
        return;
    };
    let Some(record_line) = record_line.strip_prefix(b"> ") else {
        return;
    };
    let record_line = record_line.trim_ascii();
    let record_line = record_line
        .iter()
        .rposition(|byte| *byte == b'>')
        .map_or(record_line, |end| &record_line[..end]);
    let Some(name_start) = record_line.iter().position(|byte| *byte == b'<') else {
        return;
    };
    let record_name = record_line[name_start + 1..].to_vec();

    let value_size = value_lines.iter().map(|line| line.len()).sum::<usize>()
        + value_lines.len().saturating_sub(1);
    let mut value = Vec::with_capacity(value_size);
    for (index, line) in value_lines.iter().enumerate() {
        if index > 0 {
            value.push(b'\n');
        }
        value.extend_from_slice(line);
    }
    records.push((record_name, value));
}
