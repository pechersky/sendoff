//! SD data records: `> <name>` header lines and their value lines.

use std::borrow::Cow;
use std::ops::Range;

use memchr::memmem;

use crate::core::lines::{Lines, line_end};
use crate::core::strip::{Strip, trim, trim_str};

/// A metadata header that has no `<name>`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("SD data field header is missing a name")]
pub struct MissingFieldName;

/// Where an SD data record lies, and its name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Header<'a> {
    /// The field name between the header's first `<` and last `>`.
    pub name: &'a [u8],
    /// Offset of the header line in the SD data.
    pub start: usize,
    /// Offset just past the record's run of lines and the blank line
    /// ending it.
    pub end: usize,
}

impl Header<'_> {
    /// The ranges of `metadata`, which the header was read from, before and
    /// after its name on the header line, without the line ending; `None`
    /// when the name has no closing `>`.
    pub fn around_name(&self, metadata: &[u8]) -> Option<(Range<usize>, Range<usize>)> {
        let name_start = self.name.as_ptr().addr() - metadata.as_ptr().addr();
        let name_end = name_start + self.name.len();
        if metadata.get(name_end) != Some(&b'>') {
            return None;
        }
        let line = without_ending(&metadata[..line_end(metadata, self.start)]);
        Some((self.start..name_start, name_end..line.len()))
    }
}

/// One SD data record borrowed from a block where possible.
#[derive(Debug, PartialEq, Eq)]
pub struct Record<'a> {
    /// The field name between the header's first `<` and last `>`.
    pub name: &'a [u8],
    /// The value lines without their endings, joined with `\n`.
    pub value: Cow<'a, [u8]>,
    /// Offset of the header line in the SD data.
    pub start: usize,
}

/// How a new header line is written around its name.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum HeaderStyle {
    /// `> <name>`.
    #[default]
    Plain,
    /// `>  <name>`, as RDKit and Open Babel write.
    Padded,
    /// The text around the name on the first named header line, so
    /// `>  <B>  (1)` after `>  <A>  (1)`; `Plain` without one, or when that
    /// name has no closing `>`.
    Copy,
}

impl HeaderStyle {
    /// The text around a name, which `Copy` falls back to.
    pub fn fixed(self) -> (&'static [u8], &'static [u8]) {
        match self {
            Self::Padded => (b">  <", b">"),
            Self::Plain | Self::Copy => (b"> <", b">"),
        }
    }
}

/// Append a header line of `before`, `name`, and `after`, ending in
/// `newline`.
pub fn push_header(
    data: &mut Vec<u8>,
    (before, after): (&[u8], &[u8]),
    name: &[u8],
    newline: &[u8],
) {
    data.extend_from_slice(before);
    data.extend_from_slice(name);
    data.extend_from_slice(after);
    data.extend_from_slice(newline);
}

/// Append the value line and a blank line, each ending in `newline`.
pub fn push_value(data: &mut Vec<u8>, value: &[u8], newline: &[u8]) {
    data.extend_from_slice(value);
    data.extend_from_slice(newline);
    data.extend_from_slice(newline);
}

/// The named header lines of SD data, for new records to reuse by name.
pub struct HeaderPool<'a> {
    metadata: &'a [u8],
    headers: Headers<'a>,
    /// The next header, while records take them in order.
    next: Option<Header<'a>>,
    /// The names and offsets of the headers left, sorted, once a record
    /// comes out of order; a taken header's offset becomes `usize::MAX`.
    sorted: Option<Vec<(&'a [u8], usize)>>,
    /// What new header lines hold around their names.
    pub template: (&'a [u8], &'a [u8]),
}

impl<'a> HeaderPool<'a> {
    /// Pool the headers of `metadata`, skipping those without a name, for
    /// new headers in `style`.
    pub fn new(metadata: &'a [u8], strip: Strip, style: HeaderStyle) -> Self {
        let mut headers = Records::new(metadata, 0, strip).headers();
        let next = headers.find_map(Result::ok);
        let copied = match style {
            HeaderStyle::Copy => next.and_then(|header| header.around_name(metadata)),
            HeaderStyle::Plain | HeaderStyle::Padded => None,
        };
        let template = copied.map_or(style.fixed(), |(before, after)| {
            (&metadata[before], &metadata[after])
        });
        Self {
            metadata,
            headers,
            next,
            sorted: None,
            template,
        }
    }

    /// Take the line, with its ending, of the first header left named
    /// `name`.
    pub fn take(&mut self, name: &[u8]) -> Option<&'a [u8]> {
        let start = match &mut self.sorted {
            None if self.next.is_some_and(|header| header.name == name) => {
                let start = self.next?.start;
                self.next = self.headers.find_map(Result::ok);
                start
            }
            sorted => {
                let sorted = sorted.get_or_insert_with(|| {
                    let left = self.next.into_iter();
                    let left = left.chain(self.headers.by_ref().filter_map(Result::ok));
                    let mut sorted: Vec<_> =
                        left.map(|header| (header.name, header.start)).collect();
                    sorted.sort_unstable();
                    sorted
                });
                let first = sorted.partition_point(|&(other, _)| other < name);
                let (_, start) = sorted[first..]
                    .iter_mut()
                    .take_while(|(other, _)| *other == name)
                    .find(|(_, start)| *start != usize::MAX)?;
                std::mem::replace(start, usize::MAX)
            }
        };
        Some(&self.metadata[start..line_end(self.metadata, start)])
    }
}

/// The records of SD data from an offset on: a run of lines that do not
/// strip to nothing, whose first line strips to one starting with `> `, is
/// a record; other runs are skipped.
///
/// A header without a name yields an error; iterating further resumes
/// after its run.
pub struct Records<'a> {
    metadata: &'a [u8],
    pos: usize,
    strip: Strip,
}

impl<'a> Records<'a> {
    /// Iterate over the records of `metadata` at or after `pos`.
    #[inline]
    pub fn new(metadata: &'a [u8], pos: usize, strip: Strip) -> Self {
        Self {
            metadata,
            pos,
            strip,
        }
    }

    /// Where the search for the next record starts.
    #[inline]
    pub fn pos(&self) -> usize {
        self.pos
    }

    /// Iterate over the headers alone, skipping the value lines unbuilt.
    #[inline]
    pub fn headers(self) -> Headers<'a> {
        Headers(self)
    }

    #[inline(always)]
    fn parse<const VALUE: bool>(&mut self) -> Option<Result<Record<'a>, MissingFieldName>> {
        let mut lines = Lines::new(self.metadata, self.pos);
        let found = match self.strip {
            Strip::Bytes => record_in::<VALUE>(&mut lines, &trim),
            Strip::Str => record_in::<VALUE>(&mut lines, &trim_str),
        };
        self.pos = lines.pos();
        if found.is_err() {
            self.skip_nameless_run();
        }
        found.transpose()
    }

    /// Skip the rest of a nameless header's run, as groupby does. Out of
    /// line: inlined into `record_in`, this rare path grew the record loop
    /// by a fifth and slowed `get` and `in` up to 10%.
    #[cold]
    #[inline(never)]
    fn skip_nameless_run(&mut self) {
        let mut lines = Lines::new(self.metadata, self.pos);
        match self.strip {
            Strip::Bytes => skip_run(&mut lines, &trim),
            Strip::Str => skip_run(&mut lines, &trim_str),
        }
        self.pos = lines.pos();
    }
}

impl<'a> Iterator for Records<'a> {
    type Item = Result<Record<'a>, MissingFieldName>;

    fn next(&mut self) -> Option<Self::Item> {
        self.parse::<true>()
    }
}

/// The record headers of SD data, from `Records::headers`.
pub struct Headers<'a>(Records<'a>);

impl<'a> Iterator for Headers<'a> {
    type Item = Result<Header<'a>, MissingFieldName>;

    fn next(&mut self) -> Option<Self::Item> {
        let found = self.0.parse::<false>()?;
        Some(found.map(|record| Header {
            name: record.name,
            start: record.start,
            end: self.0.pos(),
        }))
    }
}

/// Parse the next record; without `VALUE`, skip the value lines unbuilt.
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
        if !VALUE {
            skip_run(lines, trim);
            let value = Cow::Borrowed(&[][..]);
            return Ok(Some(Record { name, value, start }));
        }
        let value_start = lines.pos();
        let (mut pos, mut value_end) = (value_start, value_start);
        // Whether a line before the last ends in `\r\n`, whose `\r` the
        // value drops.
        let (mut crlf, mut ended_crlf) = (false, false);
        for line in lines.by_ref() {
            if is_blank(line, trim) {
                break;
            }
            let ending = line.len() - without_ending(line).len();
            pos += line.len();
            value_end = pos - ending;
            crlf |= ended_crlf;
            ended_crlf = ending == 2;
        }
        let value = &lines.data()[value_start..value_end];
        let value = if crlf {
            Cow::Owned(crlf_to_lf(value))
        } else {
            Cow::Borrowed(value)
        };
        return Ok(Some(Record { name, value, start }));
    }
}

/// `line` without its `\n` and then one `\r`.
#[inline(always)]
fn without_ending(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

/// `value` with each `\r\n` as `\n`.
#[cold]
#[inline(never)]
fn crlf_to_lf(value: &[u8]) -> Vec<u8> {
    let mut lf = Vec::with_capacity(value.len());
    let mut start = 0;
    for cr in memmem::find_iter(value, b"\r\n") {
        lf.extend_from_slice(&value[start..cr]);
        start = cr + 1;
    }
    lf.extend_from_slice(&value[start..]);
    lf
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
    // Not memchr or memrchr: on headers a dozen bytes long their
    // dispatched calls cost more than these inlined loops, up to 7% of
    // `in`, `get`, and `len`.
    let field = field
        .iter()
        .rposition(|&byte| byte == b'>')
        .map_or(field, |end| &field[..end]);
    let start = field
        .iter()
        .position(|&byte| byte == b'<')
        .ok_or(MissingFieldName)?;
    Ok(&field[start + 1..])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(metadata: &str) -> Result<Vec<String>, MissingFieldName> {
        Records::new(metadata.as_bytes(), 0, Strip::Bytes)
            .map(|record| Ok(String::from_utf8(record?.name.to_vec()).unwrap()))
            .collect()
    }

    #[test]
    fn record_names_follow_python_slicing() {
        assert_eq!(names("> <a<b>c>\nv\n").unwrap(), ["a<b>c"]);
        assert_eq!(names(">  < x > (1)\nv\n").unwrap(), [" x "]);
        assert_eq!(names("> <a\nv\n").unwrap(), ["a"]);
        assert_eq!(names("x\n> <skipped>\nv\n").unwrap(), Vec::<String>::new());
        assert_eq!(names("> missing\nv\n"), Err(MissingFieldName));
    }

    fn values(metadata: &str, strip: Strip) -> Vec<(String, bool)> {
        Records::new(metadata.as_bytes(), 0, strip)
            .map(|record| {
                let value = record.unwrap().value;
                let borrowed = matches!(value, Cow::Borrowed(_));
                (String::from_utf8(value.into_owned()).unwrap(), borrowed)
            })
            .collect()
    }

    #[test]
    fn values_keep_their_lines_but_not_line_endings() {
        let metadata = "> <A>\n a \n\tb\x0c\n\n> <B>\n1\r\n 2\r\r\n\n> <C>\nc\r\n \n> <D>\n";
        let expected = [
            (" a \n\tb\x0c", true),
            ("1\n 2\r", false),
            ("c", true),
            ("", true),
        ];
        assert_eq!(
            values(metadata, Strip::Bytes),
            expected.map(|(v, b)| (v.into(), b))
        );
        // Lines that strip to nothing still end values, by the block's rules.
        let wide = "> <A>\nx\u{3000}\n\u{3000}\ny\n";
        assert_eq!(values(wide, Strip::Str), [("x\u{3000}".into(), true)]);
        assert_eq!(
            values(wide, Strip::Bytes),
            [("x\u{3000}\n\u{3000}\ny".into(), true)]
        );
        assert_eq!(values("> <A>\na\r", Strip::Bytes), [("a".into(), true)]);
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
                let headers: Vec<_> = Records::new(metadata, 0, strip)
                    .headers()
                    .map(|header| header.unwrap())
                    .collect();
                let mut records = Records::new(metadata, 0, strip);
                let mut expected = Vec::new();
                while let Some(record) = records.next() {
                    let record = record.unwrap();
                    expected.push(Header {
                        name: record.name,
                        start: record.start,
                        end: records.pos(),
                    });
                }
                assert_eq!(headers, expected);
            }
        }
        let mut headers = Records::new(b"> A\n> <1>\n\n> <B>\n", 0, Strip::Bytes).headers();
        assert_eq!(headers.next(), Some(Err(MissingFieldName)));
        assert_eq!(headers.next().unwrap().unwrap().name, b"B");
        assert_eq!(headers.next(), None);
    }

    #[test]
    fn header_pool_gives_each_name_its_lines_in_order() {
        let metadata = b"> bad\n> <A>\n\n>  <A> (1) \na\n\n> 7 <B>\r\nb\n\n> <A>\na\n\n";
        let mut pool = HeaderPool::new(metadata, Strip::Bytes, HeaderStyle::Copy);
        assert_eq!(pool.template, (&b">  <"[..], &b"> (1) "[..]));
        assert_eq!(pool.take(b"A"), Some(&b">  <A> (1) \n"[..]));
        assert_eq!(pool.take(b"A"), Some(&b"> <A>\n"[..]));
        assert_eq!(pool.take(b"B"), Some(&b"> 7 <B>\r\n"[..]));
        assert_eq!(pool.take(b"A"), None);
        assert_eq!(pool.take(b"C"), None);

        let mut pool = HeaderPool::new(metadata, Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(pool.template, (&b"> <"[..], &b">"[..]));
        for name in [&b"A"[..], b"B", b"A"] {
            assert!(pool.take(name).is_some());
        }
        assert_eq!(pool.take(b"A"), None);
        let pool = HeaderPool::new(metadata, Strip::Bytes, HeaderStyle::Padded);
        assert_eq!(pool.template, (&b">  <"[..], &b">"[..]));

        let mut empty = HeaderPool::new(b"x\n", Strip::Bytes, HeaderStyle::Copy);
        assert_eq!(empty.template, HeaderStyle::Plain.fixed());
        assert_eq!(empty.take(b"A"), None);
    }
}
