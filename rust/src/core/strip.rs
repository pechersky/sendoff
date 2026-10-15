//! Python's `bytes.strip()` rules, for titles and SD data lines in either
//! mode, as RDKit and OEChem strip them.

/// Whether a byte is whitespace for Python's `bytes.strip()`; unlike
/// `u8::is_ascii_whitespace`, `\x0b` is.
#[inline]
fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | b'\x0b' | b'\x0c')
}

/// Strip whitespace from both ends like Python's `bytes.strip()`.
#[inline]
pub fn trim(bytes: &[u8]) -> &[u8] {
    let start = bytes
        .iter()
        .position(|&byte| !is_space(byte))
        .unwrap_or(bytes.len());
    trim_end(&bytes[start..])
}

/// Strip whitespace from the end like Python's `bytes.rstrip()`.
#[inline]
pub fn trim_end(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .iter()
        .rposition(|&byte| !is_space(byte))
        .map_or(0, |index| index + 1);
    &bytes[..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_python_bytes_whitespace_only() {
        assert_eq!(trim(b"\x0b\x0c \t\r\nx\x1c\n"), b"x\x1c");
        assert_eq!(trim_end(b" \tx\x1c\x0b\r\n"), b" \tx\x1c");
        assert_eq!(trim_end(b" \t\r\n"), b"");
        let wide = "\u{a0}\u{3000}x\u{2003}\x1f";
        assert_eq!(trim(wide.as_bytes()), wide.as_bytes());
    }
}
