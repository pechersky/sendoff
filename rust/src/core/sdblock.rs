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
