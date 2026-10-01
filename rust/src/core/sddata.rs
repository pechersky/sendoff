pub fn append_record(record_name: &str, value: &str) -> (String, String, String) {
    (
        format!("> <{record_name}>\n"),
        format!("{value}\n"),
        "\n".to_owned(),
    )
}
