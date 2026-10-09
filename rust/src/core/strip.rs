//! Python's `strip()` and `rstrip()` rules, for titles and SD data lines.

/// Which Python `strip()` a block's lines follow.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strip {
    /// `bytes.strip()`: ASCII whitespace only.
    Bytes,
    /// `str.strip()` over UTF-8: also `\x1c`-`\x1f` and Unicode spaces.
    Str,
}

impl Strip {
    #[inline]
    pub fn trim(self, bytes: &[u8]) -> &[u8] {
        match self {
            Strip::Bytes => trim(bytes),
            Strip::Str => trim_str(bytes),
        }
    }

    /// Strip only the end, like Python's `rstrip()`.
    #[inline]
    pub fn trim_end(self, bytes: &[u8]) -> &[u8] {
        match self {
            Strip::Bytes => trim_ascii_end_by(bytes, is_space),
            Strip::Str => trim_str_end(bytes),
        }
    }
}

/// Whether a byte is whitespace for Python's `bytes.strip()`.
#[inline]
fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | b'\x0b' | b'\x0c')
}

/// Strip whitespace from both ends like Python's `bytes.strip()`.
#[inline]
pub fn trim(bytes: &[u8]) -> &[u8] {
    trim_ascii_by(bytes, is_space)
}

#[inline]
fn trim_ascii_by(bytes: &[u8], is_space: impl Fn(u8) -> bool) -> &[u8] {
    let start = bytes
        .iter()
        .position(|&byte| !is_space(byte))
        .unwrap_or(bytes.len());
    trim_ascii_end_by(&bytes[start..], is_space)
}

#[inline]
fn trim_ascii_end_by(bytes: &[u8], is_space: impl Fn(u8) -> bool) -> &[u8] {
    let end = bytes
        .iter()
        .rposition(|&byte| !is_space(byte))
        .map_or(0, |index| index + 1);
    &bytes[..end]
}

/// Whether an ASCII byte is whitespace for Python's `str.strip()`.
#[inline]
fn is_str_space(byte: u8) -> bool {
    is_space(byte) || matches!(byte, b'\x1c'..=b'\x1f')
}

/// The UTF-8 length of a non-ASCII `str.isspace()` character starting
/// `bytes`, or 0.
#[inline]
fn wide_space_prefix(bytes: &[u8]) -> usize {
    match bytes {
        [0xc2, 0x85 | 0xa0, ..] => 2,
        [0xe1, 0x9a, 0x80, ..]
        | [0xe2, 0x80, 0x80..=0x8a | 0xa8 | 0xa9 | 0xaf, ..]
        | [0xe2, 0x81, 0x9f, ..]
        | [0xe3, 0x80, 0x80, ..] => 3,
        _ => 0,
    }
}

/// The UTF-8 length of a non-ASCII `str.isspace()` character ending
/// `bytes`, or 0.
#[inline]
fn wide_space_suffix(bytes: &[u8]) -> usize {
    match bytes {
        [.., 0xc2, 0x85 | 0xa0] => 2,
        [.., 0xe1, 0x9a, 0x80]
        | [.., 0xe2, 0x80, 0x80..=0x8a | 0xa8 | 0xa9 | 0xaf]
        | [.., 0xe2, 0x81, 0x9f]
        | [.., 0xe3, 0x80, 0x80] => 3,
        _ => 0,
    }
}

/// Strip whitespace from both ends of UTF-8 like Python's `str.strip()`.
///
/// Only complete whitespace sequences are removed, which are lead-byte
/// anchored, so stripping never changes whether the rest decodes.
#[inline]
pub fn trim_str(bytes: &[u8]) -> &[u8] {
    let trimmed = trim_ascii_by(bytes, is_str_space);
    match (trimmed.first(), trimmed.last()) {
        (Some(&first), Some(&last)) if first >= 0x80 || last >= 0x80 => trim_wide(trimmed),
        _ => trimmed,
    }
}

/// Strip whitespace from the end of UTF-8 like Python's `str.rstrip()`.
#[inline]
fn trim_str_end(bytes: &[u8]) -> &[u8] {
    let trimmed = trim_ascii_end_by(bytes, is_str_space);
    match trimmed.last() {
        Some(&last) if last >= 0x80 => trim_wide_end(trimmed),
        _ => trimmed,
    }
}

#[cold]
#[inline(never)]
fn trim_wide(bytes: &[u8]) -> &[u8] {
    let mut start = 0;
    loop {
        match wide_space_prefix(&bytes[start..]) {
            0 => break,
            len => start += len,
        }
        while start < bytes.len() && is_str_space(bytes[start]) {
            start += 1;
        }
    }
    trim_wide_end(&bytes[start..])
}

#[cold]
#[inline(never)]
fn trim_wide_end(bytes: &[u8]) -> &[u8] {
    let mut end = bytes.len();
    loop {
        match wide_space_suffix(&bytes[..end]) {
            0 => break,
            len => end -= len,
        }
        while end > 0 && is_str_space(bytes[end - 1]) {
            end -= 1;
        }
    }
    &bytes[..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_python_bytes_whitespace_only() {
        assert_eq!(trim(b"\x0b\x0c \t\r\nx\x1c\n"), b"x\x1c");
        assert_eq!(Strip::Bytes.trim_end(b" \tx\x1c\x0b\r\n"), b" \tx\x1c");
        assert_eq!(Strip::Bytes.trim_end(b" \t\r\n"), b"");
    }

    /// Every character for which Python's `str.isspace()` is true.
    const STR_SPACES: &str = "\t\n\x0b\x0c\r\x1c\x1d\x1e\x1f \u{85}\u{a0}\u{1680}\
        \u{2000}\u{2001}\u{2002}\u{2003}\u{2004}\u{2005}\u{2006}\u{2007}\u{2008}\
        \u{2009}\u{200a}\u{2028}\u{2029}\u{202f}\u{205f}\u{3000}";

    #[test]
    fn strips_python_str_whitespace() {
        for space in STR_SPACES.chars() {
            let padded = format!("{space}{space}x{space}y{space}{space}");
            let expected = format!("x{space}y");
            assert_eq!(
                trim_str(padded.as_bytes()),
                expected.as_bytes(),
                "{space:?}"
            );
            assert_eq!(trim_str(space.to_string().as_bytes()), b"");
            let kept = format!("{space}{space}x{space}y");
            assert_eq!(
                Strip::Str.trim_end(padded.as_bytes()),
                kept.as_bytes(),
                "{space:?}"
            );
        }
        assert_eq!(trim_str(STR_SPACES.as_bytes()), b"");
        assert_eq!(Strip::Str.trim_end(STR_SPACES.as_bytes()), b"");
        // Neighbours of the whitespace code points are kept.
        for kept in [
            "\u{84}", "\u{a1}", "\u{200b}", "\u{2027}", "\u{3001}", "\u{feff}",
        ] {
            assert_eq!(trim_str(kept.as_bytes()), kept.as_bytes(), "{kept:?}");
        }
        // Incomplete sequences are not whitespace.
        assert_eq!(trim_str(b"\xe3\x80x\xe3\x80"), b"\xe3\x80x\xe3\x80");
        assert_eq!(Strip::Str.trim_end(b"x\xe3\x80"), b"x\xe3\x80");
    }
}
