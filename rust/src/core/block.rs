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
    /// Offset of the header line in the SD data.
    pub start: usize,
}

/// The end of the first `M  END` line in `data[title_end..body_end]`.
fn mdl_end(data: &[u8], title_end: usize, body_end: usize) -> Option<usize> {
    let body = &data[title_end..body_end];
    let mut end = if body.starts_with(b"M  END") {
        title_end
    } else {
        title_end + MEND.find(body)? + 1
    };
    next_line(data, &mut end);
    Some(end)
}

/// Why a value cannot replace a block's title and MDL.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InvalidMolblock {
    MissingEnd,
    TrailingData,
    Delimiter,
}

impl fmt::Display for InvalidMolblock {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::MissingEnd => "molblock has no M  END line",
            Self::TrailingData => "molblock has data after its M  END line",
            Self::Delimiter => "molblock has a $$$$ line, which would end the block",
        })
    }
}

/// Append `> <name>`, the value line, and a blank line.
fn push_record(data: &mut Vec<u8>, name: &[u8], value: &[u8]) {
    data.reserve(name.len() + value.len() + 7);
    data.extend_from_slice(b"> <");
    data.extend_from_slice(name);
    data.extend_from_slice(b">\n");
    data.extend_from_slice(value);
    data.extend_from_slice(b"\n\n");
}

/// The delimiter line that ends every serialized block.
const TRAILER: &[u8] = b"$$$$\n";

/// One molecule block: a title line, MDL lines, and SD data lines.
///
/// `data` also ends with [`TRAILER`], so an unstripped block serializes
/// as a single copy of `data`.
#[derive(Clone, Debug)]
pub struct Block {
    data: Vec<u8>,
    title_end: usize,
    mdl_end: usize,
}

impl Block {
    /// Build a block from its raw lines.
    ///
    /// `delimiter` is the offset of a trailing `$$$$` line, if any. The MDL
    /// runs through the first `M  END` line and the SD data runs to the
    /// delimiter; without `M  END`, the MDL runs to the delimiter and there is
    /// no SD data. (sendoff < 1.0 kept the delimiter line in that MDL, so each
    /// write added a `$$$$` line that read back as a spurious block.)
    pub fn new(mut data: Vec<u8>, delimiter: Option<usize>) -> Self {
        if data.last() != Some(&b'\n') {
            data.push(b'\n');
        }
        let title_end = memchr(b'\n', &data).map_or(data.len(), |index| index + 1);
        let body_end = delimiter.map_or(data.len(), |offset| offset.max(title_end));
        let mdl_end = mdl_end(&data, title_end, body_end).unwrap_or(body_end);
        if data[body_end..] != *TRAILER {
            data.truncate(body_end);
            data.extend_from_slice(TRAILER);
        }
        Self {
            data,
            title_end,
            mdl_end,
        }
    }

    /// Where the SD data ends and the trailer starts.
    fn body_end(&self) -> usize {
        self.data.len() - TRAILER.len()
    }

    /// The stripped title.
    pub fn title(&self, strip: Strip) -> &[u8] {
        strip.trim(&self.data[..self.title_end])
    }

    /// The MDL lines after the title, through `M  END`, each ending in `\n`.
    pub fn mdl(&self) -> &[u8] {
        &self.data[self.title_end..self.mdl_end]
    }

    /// The stripped title line followed by the MDL lines: a mol block.
    ///
    /// Borrows from the block unless stripping changed the title.
    pub fn molblock(&self, strip: Strip) -> Cow<'_, [u8]> {
        let title = self.title(strip);
        // Trimming always drops the title line's `\n`; if that is all it
        // dropped, the title and MDL are already contiguous.
        if title.len() + 1 == self.title_end {
            return Cow::Borrowed(&self.data[..self.mdl_end]);
        }
        let mdl = self.mdl();
        let mut molblock = Vec::with_capacity(title.len() + 1 + mdl.len());
        molblock.extend_from_slice(title);
        molblock.push(b'\n');
        molblock.extend_from_slice(mdl);
        Cow::Owned(molblock)
    }

    /// The SD data lines, each ending in `\n`.
    pub fn metadata(&self) -> &[u8] {
        &self.data[self.mdl_end..self.body_end()]
    }

    /// The stripped title, the remaining lines, and `$$$$`.
    ///
    /// Borrows from the block unless stripping changed the title.
    pub fn serialized(&self, strip: Strip) -> Cow<'_, [u8]> {
        let title = self.title(strip);
        if title.len() + 1 == self.title_end {
            return Cow::Borrowed(&self.data);
        }
        let body = &self.data[self.title_end..];
        let mut output = Vec::with_capacity(title.len() + 1 + body.len());
        output.extend_from_slice(title);
        output.push(b'\n');
        output.extend_from_slice(body);
        Cow::Owned(output)
    }

    /// Append a `> <name>` record with a single value line.
    pub fn append_record(&mut self, name: &[u8], value: &[u8]) {
        self.data.truncate(self.body_end());
        push_record(&mut self.data, name, value);
        self.data.extend_from_slice(TRAILER);
    }

    /// Replace the title line with `title`, which must not contain `\n`.
    pub fn set_title(&mut self, title: &[u8]) {
        let old_end = self.title_end;
        let mut line = Vec::with_capacity(title.len() + 1);
        line.extend_from_slice(title);
        line.push(b'\n');
        self.data.splice(..old_end, line);
        self.title_end = title.len() + 1;
        self.mdl_end = self.mdl_end - old_end + self.title_end;
    }

    /// Replace the title line and MDL with `molblock`, keeping the SD data.
    ///
    /// `molblock` must run through an `M  END` line followed only by
    /// whitespace, with no `$$$$` line after the title.
    pub fn set_molblock(&mut self, molblock: &[u8]) -> Result<(), InvalidMolblock> {
        let title_end = memchr(b'\n', molblock).map_or(molblock.len(), |index| index + 1);
        let mdl_end =
            mdl_end(molblock, title_end, molblock.len()).ok_or(InvalidMolblock::MissingEnd)?;
        if !molblock[mdl_end..].iter().all(u8::is_ascii_whitespace) {
            return Err(InvalidMolblock::TrailingData);
        }
        if lines(&molblock[title_end..mdl_end]).any(|line| line.starts_with(b"$$$$")) {
            return Err(InvalidMolblock::Delimiter);
        }
        let metadata = &self.data[self.mdl_end..];
        let mut data = Vec::with_capacity(mdl_end + 1 + metadata.len());
        data.extend_from_slice(&molblock[..mdl_end]);
        if !data.ends_with(b"\n") {
            data.push(b'\n');
        }
        let new_mdl_end = data.len();
        data.extend_from_slice(metadata);
        self.data = data;
        self.title_end = title_end;
        self.mdl_end = new_mdl_end;
        Ok(())
    }

    /// Give the first record named `name` the single-line `value` and
    /// remove later ones, or append the record if there is none.
    pub fn set_record(
        &mut self,
        name: &[u8],
        value: &[u8],
        strip: Strip,
    ) -> Result<(), MissingFieldName> {
        if self.rewrite_records(name, Some(value), strip)? == 0 {
            self.append_record(name, value);
        }
        Ok(())
    }

    /// Remove every record named `name`, returning how many there were.
    pub fn remove_records(&mut self, name: &[u8], strip: Strip) -> Result<usize, MissingFieldName> {
        self.rewrite_records(name, None, strip)
    }

    /// Replace the whole SD data with `records`.
    pub fn replace_records<'a>(&mut self, records: impl IntoIterator<Item = (&'a [u8], &'a [u8])>) {
        self.data.truncate(self.mdl_end);
        for (name, value) in records {
            push_record(&mut self.data, name, value);
        }
        self.data.extend_from_slice(TRAILER);
    }

    /// Drop the records named `name`, putting `value` in place of the
    /// first, and return how many matched.
    fn rewrite_records(
        &mut self,
        name: &[u8],
        value: Option<&[u8]>,
        strip: Strip,
    ) -> Result<usize, MissingFieldName> {
        let metadata = self.metadata();
        let mut spans = Vec::new();
        let mut pos = 0;
        while let Some(record) = next_record(metadata, &mut pos, strip)? {
            if record.name == name {
                spans.push((record.start, pos));
            }
        }
        if spans.is_empty() {
            return Ok(0);
        }
        let mut data =
            Vec::with_capacity(self.data.len() + name.len() + value.map_or(0, <[u8]>::len) + 7);
        data.extend_from_slice(&self.data[..self.mdl_end]);
        let mut copied = 0;
        for (index, &(start, end)) in spans.iter().enumerate() {
            data.extend_from_slice(&metadata[copied..start]);
            if let (0, Some(value)) = (index, value) {
                push_record(&mut data, name, value);
            }
            copied = end;
        }
        data.extend_from_slice(&self.data[self.mdl_end + copied..]);
        self.data = data;
        Ok(spans.len())
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
        data.extend_from_slice(&self.data[self.mdl_end..]);
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
        let start = *pos;
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
        return Ok(Some(Record { name, value, start }));
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
        String::from_utf8(block.serialized(Strip::Bytes).into_owned()).unwrap()
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
    fn without_mend_the_mdl_runs_to_the_delimiter() {
        let block = block("t\nx\n$$$$ tail\n", true);
        assert_eq!(block.mdl(), b"x\n");
        assert_eq!(block.metadata(), b"");
        assert_eq!(written(&block), "t\nx\n$$$$\n");
        let undelimited = self::block("t\nx\n", false);
        assert_eq!(undelimited.mdl(), b"x\n");
        assert_eq!(written(&undelimited), "t\nx\n$$$$\n");
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
    fn set_molblock_replaces_title_and_mdl_only() {
        let mut block = block(" old \nsrc\nM  END\n> <A>\na\n\n$$$$\n", true);
        block.set_molblock(b"new\n  prog\n\nM  END").unwrap();
        assert_eq!(block.title(Strip::Bytes), b"new");
        assert_eq!(block.mdl(), b"  prog\n\nM  END\n");
        assert_eq!(written(&block), "new\n  prog\n\nM  END\n> <A>\na\n\n$$$$\n");
        block.set_molblock(b"M  END\nM  END\n \n\t\n").unwrap();
        assert_eq!(block.title(Strip::Bytes), b"M  END");
        assert_eq!(block.mdl(), b"M  END\n");
        assert_eq!(records(&block).unwrap(), [("A".into(), "a".into())]);
        for (bad, error) in [
            (&b"t\nx\n"[..], InvalidMolblock::MissingEnd),
            (b"M  END\n", InvalidMolblock::MissingEnd),
            (b"t\nM  END\n> <A>\n", InvalidMolblock::TrailingData),
            (b"t\n$$$$\nM  END\n", InvalidMolblock::Delimiter),
        ] {
            assert_eq!(block.set_molblock(bad), Err(error));
        }
        assert_eq!(block.title(Strip::Bytes), b"M  END");
    }

    #[test]
    fn set_title_keeps_the_rest() {
        let mut block = block(" old \nsrc\nM  END\n> <A>\na\n\n$$$$\n", true);
        block.set_title(b"a much longer title");
        assert_eq!(block.title(Strip::Bytes), b"a much longer title");
        assert_eq!(block.mdl(), b"src\nM  END\n");
        assert_eq!(records(&block).unwrap(), [("A".into(), "a".into())]);
        block.set_title(b"");
        assert_eq!(written(&block), "\nsrc\nM  END\n> <A>\na\n\n$$$$\n");
    }

    #[test]
    fn set_record_replaces_the_first_and_drops_the_rest() {
        let text = "t\nM  END\n> <A>  (1)\n1\nmore\n\njunk\n\n> <B>\nb\n\n> <A>\n2\n\n";
        let mut block = block(text, false);
        block.set_record(b"A", b"new", Strip::Bytes).unwrap();
        assert_eq!(
            written(&block),
            "t\nM  END\n> <A>\nnew\n\njunk\n\n> <B>\nb\n\n$$$$\n"
        );
        block.set_record(b"C", b"c", Strip::Bytes).unwrap();
        assert_eq!(
            records(&block).unwrap(),
            [
                ("A".into(), "new".into()),
                ("B".into(), "b".into()),
                ("C".into(), "c".into())
            ]
        );
        assert_eq!(block.remove_records(b"B", Strip::Bytes), Ok(1));
        assert_eq!(block.remove_records(b"B", Strip::Bytes), Ok(0));
        assert_eq!(
            written(&block),
            "t\nM  END\n> <A>\nnew\n\njunk\n\n> <C>\nc\n\n$$$$\n"
        );
        block.replace_records([(&b"Z"[..], &b"z"[..])]);
        assert_eq!(written(&block), "t\nM  END\n> <Z>\nz\n\n$$$$\n");
        let mut malformed = self::block("t\nM  END\n> bad\n\n", false);
        assert_eq!(
            malformed.remove_records(b"A", Strip::Bytes),
            Err(MissingFieldName)
        );
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
    fn molblock_is_the_stripped_title_and_mdl() {
        let plain = block("t\nsrc\nM  END\n> <A>\na\n\n", false);
        assert!(matches!(plain.molblock(Strip::Bytes), Cow::Borrowed(_)));
        assert_eq!(&*plain.molblock(Strip::Bytes), b"t\nsrc\nM  END\n");
        let padded = block(" t \nsrc\nM  END\n", false);
        assert_eq!(&*padded.molblock(Strip::Bytes), b"t\nsrc\nM  END\n");
        let empty = block("\nM  END\n", false);
        assert!(matches!(empty.molblock(Strip::Bytes), Cow::Borrowed(_)));
        assert_eq!(&*empty.molblock(Strip::Bytes), b"\nM  END\n");
    }

    #[test]
    fn replace_mdl_keeps_title_and_metadata() {
        let mut block = block(" t\nold\nM  END\n> <A>\na\n\n", false);
        block.replace_mdl([&b"new"[..], b"M  END\n"]);
        assert_eq!(block.mdl(), b"new\nM  END\n");
        assert_eq!(written(&block), "t\nnew\nM  END\n> <A>\na\n\n$$$$\n");
    }
}
