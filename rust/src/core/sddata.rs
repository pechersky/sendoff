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
