//! SD data records: `> <name>` header lines and their value lines.

use std::borrow::Cow;
use std::fmt;

use memchr::memchr;

use crate::core::lines::Lines;
use crate::core::strip::{Strip, trim, trim_str};

/// A metadata header that has no `<name>`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MissingFieldName;

impl fmt::Display for MissingFieldName {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("SD data field header is missing a name")
    }
}

/// Where an SD data record starts, and its name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Header<'a> {
    /// The field name between the header's first `<` and last `>`.
    pub name: &'a [u8],
    /// Offset of the header line in the SD data.
    pub start: usize,
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

/// Append `> <name>`, the value line, and a blank line, each ending in
/// `newline`.
pub fn push_record(data: &mut Vec<u8>, name: &[u8], value: &[u8], newline: &[u8]) {
    data.reserve(name.len() + value.len() + 4 + 3 * newline.len());
    data.extend_from_slice(b"> <");
    data.extend_from_slice(name);
    data.push(b'>');
    data.extend_from_slice(newline);
    data.extend_from_slice(value);
    data.extend_from_slice(newline);
    data.extend_from_slice(newline);
}

/// Parse the next record at or after `*pos` in `metadata`, as Python's
/// groupby over stripped lines does: a run of non-blank lines whose first
/// line starts with `> ` is a record; other runs are skipped.
pub fn next_record<'a>(
    metadata: &'a [u8],
    pos: &mut usize,
    strip: Strip,
) -> Result<Option<Record<'a>>, MissingFieldName> {
    match strip {
        Strip::Bytes => next_record_by::<true>(metadata, pos, trim),
        Strip::Str => next_record_by::<true>(metadata, pos, trim_str),
    }
}

/// Like `next_record`, but skip the value.
pub fn next_header<'a>(
    metadata: &'a [u8],
    pos: &mut usize,
    strip: Strip,
) -> Result<Option<Header<'a>>, MissingFieldName> {
    let found = match strip {
        Strip::Bytes => next_record_by::<false>(metadata, pos, trim),
        Strip::Str => next_record_by::<false>(metadata, pos, trim_str),
    };
    found.map(|record| {
        record.map(|record| Header {
            name: record.name,
            start: record.start,
        })
    })
}

/// Parse the next record; without `VALUE`, skip the value lines unbuilt.
#[inline]
fn next_record_by<'a, const VALUE: bool>(
    metadata: &'a [u8],
    pos: &mut usize,
    trim: impl Fn(&[u8]) -> &[u8],
) -> Result<Option<Record<'a>>, MissingFieldName> {
    let mut lines = Lines::new(metadata, *pos);
    let found = record_in::<VALUE>(&mut lines, &trim);
    *pos = lines.pos();
    found
}

#[inline(always)]
fn record_in<'a, const VALUE: bool>(
    lines: &mut Lines<'a>,
    trim: &impl Fn(&[u8]) -> &[u8],
) -> Result<Option<Record<'a>>, MissingFieldName> {
    loop {
        let start = lines.pos();
        let Some(line) = lines.next() else {
            return Ok(None);
        };
        let header = trim(line);
        if header.is_empty() {
            continue;
        }
        let Some(rest) = header.strip_prefix(b"> ") else {
            skip_run(lines, trim);
            continue;
        };
        let name = record_name(rest)?;
        let mut value: Cow<'a, [u8]> = Cow::Borrowed(&[]);
        if !VALUE {
            skip_run(lines, trim);
            return Ok(Some(Record { name, value, start }));
        }
        let mut first = true;
        for line in lines.by_ref() {
            let line = trim(line);
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
}

/// Consume lines through the next blank one.
#[inline(always)]
fn skip_run(lines: &mut Lines<'_>, trim: &impl Fn(&[u8]) -> &[u8]) {
    for line in lines {
        if is_blank(line, trim) {
            break;
        }
    }
}

/// Whether `line` strips to nothing, deciding most lines by their first
/// byte: printable ASCII is whitespace for neither `str` nor `bytes`.
#[inline(always)]
fn is_blank(line: &[u8], trim: &impl Fn(&[u8]) -> &[u8]) -> bool {
    match line.first() {
        Some(&byte) if byte > b' ' && byte < 0x80 => false,
        _ => trim(line).is_empty(),
    }
}

/// Extract `name` from the text after `> ` in a stripped `> <name>` header.
///
/// Python strips `rest` again first, but that can only remove leading
/// whitespace, which holds no `<` or `>`, so the name is the same.
fn record_name(field: &[u8]) -> Result<&[u8], MissingFieldName> {
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

    fn names(metadata: &str) -> Result<Vec<String>, MissingFieldName> {
        let mut pos = 0;
        let mut found = Vec::new();
        while let Some(record) = next_record(metadata.as_bytes(), &mut pos, Strip::Bytes)? {
            found.push(String::from_utf8(record.name.to_vec()).unwrap());
        }
        Ok(found)
    }

    #[test]
    fn record_names_follow_python_slicing() {
        assert_eq!(names("> <a<b>c>\nv\n").unwrap(), ["a<b>c"]);
        assert_eq!(names(">  < x > (1)\nv\n").unwrap(), [" x "]);
        assert_eq!(names("> <a\nv\n").unwrap(), ["a"]);
        assert_eq!(names("x\n> <skipped>\nv\n").unwrap(), Vec::<String>::new());
        assert_eq!(names("> missing\nv\n"), Err(MissingFieldName));
    }

    #[test]
    fn headers_match_records() {
        for metadata in [
            "> <A>\n1\n2\n\n> <B>\n\njunk\nmore\n\n> <A>\n3\n",
            "\n\n> <A> (1)\n\u{3000}\n> <B>\n",
            "",
        ] {
            for strip in [Strip::Bytes, Strip::Str] {
                let metadata = metadata.as_bytes();
                let (mut pos, mut headers) = (0, Vec::new());
                while let Some(header) = next_header(metadata, &mut pos, strip).unwrap() {
                    headers.push((header.name, header.start, pos));
                }
                let (mut pos, mut records) = (0, Vec::new());
                while let Some(record) = next_record(metadata, &mut pos, strip).unwrap() {
                    records.push((record.name, record.start, pos));
                }
                assert_eq!(headers, records);
            }
        }
        assert!(next_header(b"> A\n1\n", &mut 0, Strip::Bytes).is_err());
    }
}
