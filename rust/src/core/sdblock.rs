//! Pure Rust operations for SDF blocks.

/// Collect lines through the first `M  END` line, including that line.
pub fn parse_mdl<I, E>(lines: &mut I) -> Result<Vec<String>, E>
where
    I: Iterator<Item = Result<String, E>>,
{
    let mut mdl_lines = Vec::new();
    for line in lines {
        let line = line?;
        let is_end = line.starts_with("M  END");
        mdl_lines.push(line);
        if is_end {
            break;
        }
    }
    Ok(mdl_lines)
}

/// Collect metadata lines before the first `$$$$` delimiter.
pub fn parse_metadata<I, E>(lines: &mut I) -> Result<Vec<String>, E>
where
    I: Iterator<Item = Result<String, E>>,
{
    let mut metadata_lines = Vec::new();
    for line in lines {
        let line = line?;
        if line.starts_with("$$$$") {
            break;
        }
        metadata_lines.push(line);
    }
    Ok(metadata_lines)
}
