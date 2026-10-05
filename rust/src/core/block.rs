//! One molecule block stored as a single contiguous byte buffer.
//!
//! Every line in the buffer ends with `\n`, so serializing is a copy and
//! records, MDL lines, and counts are parsed from spans without per-line
//! allocation. Text blocks hold UTF-8 and are decoded only when values are
//! requested.

use std::borrow::Cow;
use std::fmt;
use std::sync::LazyLock;

use memchr::{memchr, memmem};

use crate::core::ctable::{self, CTableError, CTableFormat};

static MEND: LazyLock<memmem::Finder<'static>> = LazyLock::new(|| memmem::Finder::new(b"\nM  END"));

/// Whether a byte is whitespace for Python's `bytes.strip()`.
#[inline]
pub fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | b'\x0b' | b'\x0c')
}

/// Strip whitespace from both ends like Python's `bytes.strip()`.
#[inline]
pub fn trim(bytes: &[u8]) -> &[u8] {
    trim_ascii_by(bytes, is_space)
}

#[inline]
fn trim_ascii_by(bytes: &[u8], is_space: impl Fn(u8) -> bool) -> &[u8] {
    let Some(start) = bytes.iter().position(|&byte| !is_space(byte)) else {
        return &[];
    };
    let end = bytes
        .iter()
        .rposition(|&byte| !is_space(byte))
        .map_or(start, |index| index + 1);
    &bytes[start..end]
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
    let mut end = bytes.len();
    loop {
        match wide_space_suffix(&bytes[start..end]) {
            0 => break,
            len => end -= len,
        }
        while end > start && is_str_space(bytes[end - 1]) {
            end -= 1;
        }
    }
    &bytes[start..end]
}

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
}

/// Return the line starting at `*pos`, including its `\n`, and advance.
#[inline]
fn next_line<'a>(data: &'a [u8], pos: &mut usize) -> &'a [u8] {
    let rest = &data[*pos..];
    let len = memchr(b'\n', rest).map_or(rest.len(), |index| index + 1);
    *pos += len;
    &rest[..len]
}

/// Iterate over `\n`-terminated lines, yielding them without the `\n`.
pub fn lines(data: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut pos = 0;
    std::iter::from_fn(move || {
        (pos < data.len()).then(|| {
            let line = next_line(data, &mut pos);
            line.strip_suffix(b"\n").unwrap_or(line)
        })
    })
}

/// A metadata header that has no `<name>`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MissingFieldName;

impl fmt::Display for MissingFieldName {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("SD data field header is missing a name")
    }
}

/// One SD data record borrowed from a block where possible.
#[derive(Debug, PartialEq, Eq)]
pub struct Record<'a> {
    /// The field name between the header's first `<` and last `>`.
    pub name: &'a [u8],
    /// Stripped value lines joined with `\n`.
    pub value: Cow<'a, [u8]>,
}

/// One molecule block: a title line, MDL lines, and SD data lines.
#[derive(Clone, Debug)]
pub struct Block {
    data: Vec<u8>,
    title_end: usize,
    mdl_end: usize,
}

impl Block {
    /// Build a block from its raw lines.
    ///
    /// `delimiter` is the offset of a trailing `$$$$` line, if any. Like the
    /// Python reference, the MDL runs through the first `M  END` line and the
    /// SD data runs to the delimiter; without `M  END`, the MDL swallows the
    /// delimiter line and there is no SD data.
    pub fn new(mut data: Vec<u8>, delimiter: Option<usize>) -> Self {
        if data.last() != Some(&b'\n') {
            data.push(b'\n');
        }
        let title_end = memchr(b'\n', &data).map_or(data.len(), |index| index + 1);
        let body_end = delimiter.map_or(data.len(), |offset| offset.max(title_end));
        let body = &data[title_end..body_end];
        let mend = if body.starts_with(b"M  END") {
            Some(title_end)
        } else {
            MEND.find(body).map(|index| title_end + index + 1)
        };
        let Some(mut mdl_end) = mend else {
            let mdl_end = data.len();
            return Self {
                data,
                title_end,
                mdl_end,
            };
        };
        next_line(&data, &mut mdl_end);
        data.truncate(body_end);
        Self {
            data,
            title_end,
            mdl_end,
        }
    }

    /// The stripped title.
    pub fn title(&self, strip: Strip) -> &[u8] {
        strip.trim(&self.data[..self.title_end])
    }

    /// The MDL lines after the title, through `M  END`, each ending in `\n`.
    pub fn mdl(&self) -> &[u8] {
        &self.data[self.title_end..self.mdl_end]
    }

    /// The SD data lines, each ending in `\n`.
    pub fn metadata(&self) -> &[u8] {
        &self.data[self.mdl_end..]
    }

    /// Write the stripped title, the remaining lines, and `$$$$`.
    pub fn serialize_into(&self, output: &mut Vec<u8>, strip: Strip) {
        let title = self.title(strip);
        let body = &self.data[self.title_end..];
        output.reserve(title.len() + body.len() + 6);
        output.extend_from_slice(title);
        output.push(b'\n');
        output.extend_from_slice(body);
        output.extend_from_slice(b"$$$$\n");
    }

    /// The serialized size in bytes.
    pub fn serialized_len(&self, strip: Strip) -> usize {
        self.title(strip).len() + 1 + (self.data.len() - self.title_end) + 5
    }

    /// Append a `> <name>` record with a single value line.
    pub fn append_record(&mut self, name: &[u8], value: &[u8]) {
        self.data.reserve(name.len() + value.len() + 7);
        self.data.extend_from_slice(b"> <");
        self.data.extend_from_slice(name);
        self.data.extend_from_slice(b">\n");
        self.data.extend_from_slice(value);
        self.data.extend_from_slice(b"\n\n");
    }

    /// Replace the MDL lines, keeping the raw title line and SD data.
    pub fn replace_mdl<'a>(&mut self, mdl: impl IntoIterator<Item = &'a [u8]>) {
        let mut data = Vec::with_capacity(self.data.len() + 64);
        data.extend_from_slice(&self.data[..self.title_end]);
        for line in mdl {
            data.extend_from_slice(line);
            if !line.ends_with(b"\n") {
                data.push(b'\n');
            }
        }
        let mdl_end = data.len();
        data.extend_from_slice(self.metadata());
        self.data = data;
        self.mdl_end = mdl_end;
    }

    /// Parse the next record at or after `*pos` in the SD data.
    pub fn next_record(
        &self,
        pos: &mut usize,
        strip: Strip,
    ) -> Result<Option<Record<'_>>, MissingFieldName> {
        next_record(self.metadata(), pos, strip)
    }

    /// The format and atom and bond counts declared by the MDL header.
    ///
    /// Equivalent to parsing a `CTable` from the title and MDL lines, without
    /// copying them.
    pub fn counts(&self) -> Result<(CTableFormat, i64, i64), CTableError> {
        let header: Vec<Cow<'_, str>> = std::iter::once(&b""[..])
            .chain(lines(self.mdl()))
            .take(6)
            .map(String::from_utf8_lossy)
            .collect();
        let header: Vec<&str> = header.iter().map(AsRef::as_ref).collect();
        ctable::parse_header(&header).map(|(format, _, atoms, bonds)| (format, atoms, bonds))
    }
}

/// Parse the next record at or after `*pos`, as Python's groupby over
/// stripped lines does: a run of non-blank lines whose first line starts with
/// `> ` is a record; other runs are skipped.
pub fn next_record<'a>(
    metadata: &'a [u8],
    pos: &mut usize,
    strip: Strip,
) -> Result<Option<Record<'a>>, MissingFieldName> {
    match strip {
        Strip::Bytes => next_record_by(metadata, pos, trim),
        Strip::Str => next_record_by(metadata, pos, trim_str),
    }
}

#[inline]
fn next_record_by<'a>(
    metadata: &'a [u8],
    pos: &mut usize,
    trim: impl Fn(&[u8]) -> &[u8],
) -> Result<Option<Record<'a>>, MissingFieldName> {
    while *pos < metadata.len() {
        let header = trim(next_line(metadata, pos));
        if header.is_empty() {
            continue;
        }
        let Some(rest) = header.strip_prefix(b"> ") else {
            while *pos < metadata.len() && !trim(next_line(metadata, pos)).is_empty() {}
            continue;
        };
        let name = record_name(rest, &trim)?;
        let mut value: Cow<'a, [u8]> = Cow::Borrowed(&[]);
        let mut first = true;
        while *pos < metadata.len() {
            let line = trim(next_line(metadata, pos));
            if line.is_empty() {
                break;
            }
            if first {
                value = Cow::Borrowed(line);
                first = false;
            } else {
                let owned = value.to_mut();
                owned.push(b'\n');
                owned.extend_from_slice(line);
            }
        }
        return Ok(Some(Record { name, value }));
    }
    Ok(None)
}

/// Extract `name` from the text after `> ` in a `> <name>` header.
fn record_name<'a>(
    rest: &'a [u8],
    trim: &impl Fn(&[u8]) -> &[u8],
) -> Result<&'a [u8], MissingFieldName> {
    let field = trim(rest);
    let field = field
        .iter()
        .rposition(|&byte| byte == b'>')
        .map_or(field, |end| &field[..end]);
    let start = memchr(b'<', field).ok_or(MissingFieldName)?;
    Ok(&field[start + 1..])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn block(text: &str, delimited: bool) -> Block {
        let delimiter = delimited.then(|| {
            text.trim_end_matches('\n')
                .rfind('\n')
                .map_or(0, |index| index + 1)
        });
        Block::new(text.as_bytes().to_vec(), delimiter)
    }

    fn written(block: &Block) -> String {
        let mut output = Vec::new();
        block.serialize_into(&mut output, Strip::Bytes);
        assert_eq!(output.len(), block.serialized_len(Strip::Bytes));
        String::from_utf8(output).unwrap()
    }

    fn records(block: &Block) -> Result<Vec<(String, String)>, MissingFieldName> {
        let mut pos = 0;
        let mut found = Vec::new();
        while let Some(record) = block.next_record(&mut pos, Strip::Bytes)? {
            found.push((
                String::from_utf8(record.name.to_vec()).unwrap(),
                String::from_utf8(record.value.into_owned()).unwrap(),
            ));
        }
        Ok(found)
    }

    #[test]
    fn splits_title_mdl_and_metadata() {
        let block = block(" t \r\nsrc\nM  END\n> <A>\n a \nb\n\n> <B>\n\n$$$$\n", true);
        assert_eq!(block.title(Strip::Bytes), b"t");
        assert_eq!(block.mdl(), b"src\nM  END\n");
        assert_eq!(block.metadata(), b"> <A>\n a \nb\n\n> <B>\n\n");
        assert_eq!(
            records(&block).unwrap(),
            [("A".into(), "a\nb".into()), ("B".into(), String::new())]
        );
        assert_eq!(
            written(&block),
            "t\nsrc\nM  END\n> <A>\n a \nb\n\n> <B>\n\n$$$$\n"
        );
    }

    #[test]
    fn without_mend_the_mdl_swallows_the_delimiter() {
        let block = block("t\nx\n$$$$ tail\n", true);
        assert_eq!(block.mdl(), b"x\n$$$$ tail\n");
        assert_eq!(block.metadata(), b"");
        assert_eq!(written(&block), "t\nx\n$$$$ tail\n$$$$\n");
    }

    #[test]
    fn a_delimiter_title_is_the_whole_block() {
        let block = block("$$$$\n", true);
        assert_eq!(block.title(Strip::Bytes), b"$$$$");
        assert_eq!(written(&block), "$$$$\n$$$$\n");
    }

    #[test]
    fn unterminated_input_gets_a_newline() {
        let block = block("t\nM  END\n> <A>\nv", false);
        assert_eq!(records(&block).unwrap(), [("A".into(), "v".into())]);
        assert_eq!(written(&block), "t\nM  END\n> <A>\nv\n$$$$\n");
    }

    #[test]
    fn record_names_follow_python_slicing() {
        let names = |text: &str| {
            records(&block(&format!("t\nM  END\n{text}\nv\n"), false))
                .map(|found| found.into_iter().map(|(name, _)| name).collect::<Vec<_>>())
        };
        assert_eq!(names("> <a<b>c>").unwrap(), ["a<b>c"]);
        assert_eq!(names(">  < x > (1)").unwrap(), [" x "]);
        assert_eq!(names("> <a").unwrap(), ["a"]);
        assert_eq!(names("x\n> <skipped>").unwrap(), Vec::<String>::new());
        assert_eq!(names("> missing"), Err(MissingFieldName));
    }

    #[test]
    fn appended_records_follow_existing_metadata() {
        let mut block = block("t\nM  END\n> <A>\na\n\n$$$$\n", true);
        block.append_record(b"B", b"b");
        assert_eq!(
            records(&block).unwrap(),
            [("A".into(), "a".into()), ("B".into(), "b".into())]
        );
        assert_eq!(written(&block), "t\nM  END\n> <A>\na\n\n> <B>\nb\n\n$$$$\n");
    }

    #[test]
    fn strips_python_bytes_whitespace_only() {
        assert_eq!(trim(b"\x0b\x0c \t\r\nx\x1c\n"), b"x\x1c");
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
        }
        assert_eq!(trim_str(STR_SPACES.as_bytes()), b"");
        // Neighbours of the whitespace code points are kept.
        for kept in [
            "\u{84}", "\u{a1}", "\u{200b}", "\u{2027}", "\u{3001}", "\u{feff}",
        ] {
            assert_eq!(trim_str(kept.as_bytes()), kept.as_bytes(), "{kept:?}");
        }
        // Incomplete sequences are not whitespace.
        assert_eq!(trim_str(b"\xe3\x80x\xe3\x80"), b"\xe3\x80x\xe3\x80");
    }

    #[test]
    fn str_strip_ends_records_at_unicode_blank_lines() {
        let text = "\u{3000}t\u{a0}\nM  END\n> <A>\u{2003}\n\u{3000}\n\u{1c}tail\n\n";
        let block = block(text, false);
        assert_eq!(block.title(Strip::Str), b"t");
        let mut pos = 0;
        let record = block.next_record(&mut pos, Strip::Str).unwrap().unwrap();
        assert_eq!((record.name, &*record.value), (&b"A"[..], &b""[..]));
        assert_eq!(block.next_record(&mut pos, Strip::Str), Ok(None));
        assert_eq!(block.title(Strip::Bytes), "\u{3000}t\u{a0}".as_bytes());
        let record = block.next_record(&mut 0, Strip::Bytes).unwrap().unwrap();
        assert_eq!(record.name, b"A");
        assert_eq!(&*record.value, "\u{3000}\n\u{1c}tail".as_bytes());
    }

    #[test]
    fn counts_match_ctable_parsing() {
        let v2000 = block(
            "t\n\n\n  3  2  0  0  0  0            999 V2000\nM  END\n",
            false,
        );
        assert_eq!(v2000.counts(), Ok((CTableFormat::V2000, 3, 2)));
        let v3000 = block(
            "t\n\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 4 5 0 0 0\nM  END\n",
            false,
        );
        assert_eq!(v3000.counts(), Ok((CTableFormat::V3000, 4, 5)));
        assert_eq!(
            block("t\nM  END\n", false).counts(),
            Err(CTableError::Truncated {
                needed: 3,
                found: 2
            })
        );
    }

    #[test]
    fn replace_mdl_keeps_title_and_metadata() {
        let mut block = block(" t\nold\nM  END\n> <A>\na\n\n", false);
        block.replace_mdl([&b"new"[..], b"M  END\n"]);
        assert_eq!(block.mdl(), b"new\nM  END\n");
        assert_eq!(written(&block), "t\nnew\nM  END\n> <A>\na\n\n$$$$\n");
    }
}
