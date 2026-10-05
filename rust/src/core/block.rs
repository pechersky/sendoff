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
    let Some(start) = bytes.iter().position(|&byte| !is_space(byte)) else {
        return &[];
    };
    let end = bytes
        .iter()
        .rposition(|&byte| !is_space(byte))
        .map_or(start, |index| index + 1);
    &bytes[start..end]
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
    pub fn title(&self) -> &[u8] {
        trim(&self.data[..self.title_end])
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
    pub fn serialize_into(&self, output: &mut Vec<u8>) {
        let title = self.title();
        let body = &self.data[self.title_end..];
        output.reserve(title.len() + body.len() + 6);
        output.extend_from_slice(title);
        output.push(b'\n');
        output.extend_from_slice(body);
        output.extend_from_slice(b"$$$$\n");
    }

    /// The serialized size in bytes.
    pub fn serialized_len(&self) -> usize {
        self.title().len() + 1 + (self.data.len() - self.title_end) + 5
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
    pub fn next_record(&self, pos: &mut usize) -> Result<Option<Record<'_>>, MissingFieldName> {
        next_record(self.metadata(), pos)
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
        let name = record_name(rest)?;
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
fn record_name(rest: &[u8]) -> Result<&[u8], MissingFieldName> {
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
        block.serialize_into(&mut output);
        assert_eq!(output.len(), block.serialized_len());
        String::from_utf8(output).unwrap()
    }

    fn records(block: &Block) -> Result<Vec<(String, String)>, MissingFieldName> {
        let mut pos = 0;
        let mut found = Vec::new();
        while let Some(record) = block.next_record(&mut pos)? {
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
        assert_eq!(block.title(), b"t");
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
        assert_eq!(block.title(), b"$$$$");
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
