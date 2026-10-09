//! One molecule block stored as a single contiguous byte buffer.
//!
//! Every line in the buffer ends with `\n`, so serializing is a copy and
//! records, MDL lines, and counts are parsed from spans without per-line
//! allocation. Text blocks hold UTF-8 and are decoded only when values are
//! requested.

use std::borrow::Cow;
use std::ops::Range;
use std::sync::LazyLock;

use memchr::{memmem, memrchr};

use crate::core::ctable::error::CTableError;
use crate::core::ctable::header::{Counts, HEADER_LINES, parse_header};
use crate::core::lines::{find_newline, line_end, lines};
use crate::core::newline::{
    CRLF_TRAILER, TRAILER, missing_ending, newline, title_newline, trailer,
};
use crate::core::record::{
    HeaderPool, HeaderStyle, Headers, MissingFieldName, Records, push_header, push_value,
};
use crate::core::strip::Strip;

static MEND: LazyLock<memmem::Finder<'static>> = LazyLock::new(|| memmem::Finder::new(b"\nM  END"));

/// The end of the first `M  END` line in `data[title_end..body_end]`.
fn mdl_end(data: &[u8], title_end: usize, body_end: usize) -> Option<usize> {
    let body = &data[title_end..body_end];
    let start = if body.starts_with(b"M  END") {
        title_end
    } else {
        title_end + MEND.find(body)? + 1
    };
    Some(line_end(data, start))
}

/// Where the title, MDL, and SD data end in raw lines ending in `\n`, given
/// the offset of a `$$$$` line. A `$$$$` first line ends the record, which
/// then has no title line.
///
/// The MDL runs through the first `M  END` line and the SD data runs to the
/// delimiter; without `M  END`, the MDL runs to the delimiter and there is
/// no SD data.
fn layout(data: &[u8], delimiter: Option<usize>) -> Layout {
    let title_end = match delimiter {
        Some(0) => 0,
        _ => line_end(data, 0),
    };
    let body_end = delimiter.map_or(data.len(), |offset| offset.max(title_end));
    let mdl_end = mdl_end(data, title_end, body_end).unwrap_or(body_end);
    Layout {
        title_end,
        mdl_end,
        body_end,
    }
}

/// Where a block's title, MDL, and SD data end.
struct Layout {
    title_end: usize,
    mdl_end: usize,
    body_end: usize,
}

/// A block's data, ending with a `$$$$` line, and where its title and MDL
/// end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Parts<'a> {
    pub data: &'a [u8],
    pub title_end: usize,
    pub mdl_end: usize,
}

/// Why a value cannot replace a block's title and MDL.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum InvalidMolblock {
    #[error("molblock has no M  END line")]
    MissingEnd,
    #[error("molblock has data after its M  END line")]
    TrailingData,
    #[error("molblock has a $$$$ line, which would end the block")]
    Delimiter,
}

/// One molecule block: a title line, MDL lines, and SD data lines.
///
/// `data` also ends with a `$$$$` line, so an unstripped block serializes
/// as a single copy of `data`. Lines the block writes or ends itself (a
/// stripped title, new records, a missing or normalized `$$$$`, a last line
/// without `\n`) end like the title line, so CRLF input stays CRLF; a line
/// ending in a bare `\r` gets just `\n`.
///
/// A record that is only a `$$$$` line has no title line (`title_end` is 0),
/// since a line starting with `$$$$` always ends a record; it writes back as
/// that line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Block {
    data: Vec<u8>,
    title_end: usize,
    mdl_end: usize,
}

impl Block {
    /// Build a block from its raw lines.
    ///
    /// `delimiter` is the offset of a trailing `$$$$` line, if any, and
    /// `layout` finds the title, MDL, and SD data. (Pure-Python sendoff
    /// 0.1.x kept the delimiter line in an MDL without `M  END`, so each
    /// write added a `$$$$` line that read back as a spurious block.)
    pub fn new(mut data: Vec<u8>, delimiter: Option<usize>) -> Self {
        if data.last() != Some(&b'\n') {
            let ending = missing_ending(&data, title_newline(&data));
            data.extend_from_slice(ending);
        }
        let Layout {
            title_end,
            mdl_end,
            body_end,
        } = layout(&data, delimiter);
        if !matches!(&data[body_end..], TRAILER | CRLF_TRAILER) {
            let trailer = trailer(title_newline(&data));
            data.truncate(body_end);
            data.extend_from_slice(trailer);
        }
        Self {
            data,
            title_end,
            mdl_end,
        }
    }

    /// The data, ending with a `$$$$` line, and where the title and MDL end:
    /// what `from_parts` rebuilds the block from.
    ///
    /// Re-parsing the serialized block does not always restore it: a title
    /// set to start with `$$$$` reads back as a delimiter, and SD data can
    /// follow an MDL without `M  END`.
    pub fn parts(&self) -> Parts<'_> {
        Parts {
            data: &self.data,
            title_end: self.title_end,
            mdl_end: self.mdl_end,
        }
    }

    /// Rebuild a block from `parts`, or `None` if no block has them.
    pub fn from_parts(data: Vec<u8>, title_end: usize, mdl_end: usize) -> Option<Self> {
        let delimiter = trailer(newline(&data));
        let body_end = data.len().checked_sub(delimiter.len())?;
        let valid = data.ends_with(delimiter)
            && if title_end == 0 {
                mdl_end == 0 && body_end == 0
            } else {
                find_newline(&data) == Some(title_end - 1)
                    && title_end <= mdl_end
                    && mdl_end <= body_end
                    && data[mdl_end - 1] == b'\n'
            };
        valid.then_some(Self {
            data,
            title_end,
            mdl_end,
        })
    }

    /// The `$$$$` line that ends the data.
    fn trailer(&self) -> &'static [u8] {
        trailer(newline(&self.data))
    }

    /// Where the SD data ends and the trailer starts.
    fn body_end(&self) -> usize {
        self.data.len() - self.trailer().len()
    }

    /// The title line's ending, which lines the block writes also use; the
    /// `$$$$` line's for a block without a title line.
    fn newline(&self) -> &'static [u8] {
        match self.title_end {
            0 => newline(&self.data),
            end => newline(&self.data[..end]),
        }
    }

    /// Give a block that is only a `$$$$` line an empty title line, so that
    /// lines added to it do not read back as its title.
    fn ensure_title_line(&mut self) {
        if self.title_end == 0 {
            let newline = self.newline();
            self.data.splice(..0, newline.iter().copied());
            self.title_end = newline.len();
            self.mdl_end = newline.len();
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

    /// The stripped title line followed by the MDL lines: a mol block.
    ///
    /// Borrows from the block unless stripping changed the title.
    pub fn molblock(&self, strip: Strip) -> Cow<'_, [u8]> {
        self.retitled(strip, self.mdl_end)
    }

    /// The stripped title line and `data[title_end..end]`.
    ///
    /// Borrows from the block unless stripping changed more than the
    /// title line's ending.
    fn retitled(&self, strip: Strip, end: usize) -> Cow<'_, [u8]> {
        let title = self.title(strip);
        let newline = self.newline();
        if self.title_end == 0 || title.len() + newline.len() == self.title_end {
            return Cow::Borrowed(&self.data[..end]);
        }
        let rest = &self.data[self.title_end..end];
        let mut output = Vec::with_capacity(title.len() + newline.len() + rest.len());
        output.extend_from_slice(title);
        output.extend_from_slice(newline);
        output.extend_from_slice(rest);
        Cow::Owned(output)
    }

    /// The SD data lines, each ending in `\n`.
    pub fn metadata(&self) -> &[u8] {
        &self.data[self.mdl_end..self.body_end()]
    }

    /// The stripped title, the remaining lines, and `$$$$`.
    ///
    /// Borrows from the block unless stripping changed the title.
    pub fn serialized(&self, strip: Strip) -> Cow<'_, [u8]> {
        self.retitled(strip, self.data.len())
    }

    /// Append what `serialized` returns to `out`.
    pub fn serialize_into(&self, strip: Strip, out: &mut Vec<u8>) {
        if self.title_end > 0 {
            out.extend_from_slice(self.title(strip));
            out.extend_from_slice(self.newline());
        }
        out.extend_from_slice(&self.data[self.title_end..]);
    }

    /// Append a record with a single value line, after a blank line unless
    /// the SD data is empty or already ends with one, its header in `style`.
    pub fn append_record(&mut self, name: &[u8], value: &[u8], strip: Strip, style: HeaderStyle) {
        self.ensure_title_line();
        let template = match style {
            HeaderStyle::Copy => self.template(strip),
            HeaderStyle::Plain | HeaderStyle::Padded => None,
        };
        let (newline, trailer) = (self.newline(), self.trailer());
        let needs_separator = self.ends_in_text(strip);
        self.data.truncate(self.body_end());
        if needs_separator {
            self.data.extend_from_slice(newline);
        }
        match template {
            Some((before, after)) => {
                self.data.extend_from_within(before);
                self.data.extend_from_slice(name);
                self.data.extend_from_within(after);
                self.data.extend_from_slice(newline);
            }
            None => push_header(&mut self.data, style.fixed(), name, newline),
        }
        push_value(&mut self.data, value, newline);
        self.data.extend_from_slice(trailer);
    }

    /// The ranges of the data before and after the name on the first named
    /// header line, without its ending, for new headers to copy; `None`
    /// without one, or when that name has no closing `>`.
    fn template(&self, strip: Strip) -> Option<(Range<usize>, Range<usize>)> {
        let header = self.headers(strip).find_map(Result::ok)?;
        let (before, after) = header.around_name(self.metadata())?;
        let shift = |range: Range<usize>| range.start + self.mdl_end..range.end + self.mdl_end;
        Some((shift(before), shift(after)))
    }

    /// Whether the last SD data line is not blank, so that a header
    /// written after it would read as part of its run.
    fn ends_in_text(&self, strip: Strip) -> bool {
        let metadata = self.metadata();
        let body = metadata.strip_suffix(b"\n").unwrap_or(metadata);
        let start = memrchr(b'\n', body).map_or(0, |index| index + 1);
        !strip.trim(&body[start..]).is_empty()
    }

    /// Replace the title line with `title`, which must not contain `\n`,
    /// keeping the line's ending.
    pub fn set_title(&mut self, title: &[u8]) {
        let old_end = self.title_end;
        let newline = self.newline();
        let mut line = Vec::with_capacity(title.len() + newline.len());
        line.extend_from_slice(title);
        line.extend_from_slice(newline);
        self.data.splice(..old_end, line);
        self.title_end = title.len() + newline.len();
        self.mdl_end = self.mdl_end - old_end + self.title_end;
    }

    /// Replace the title line and MDL with `molblock`, keeping the SD data.
    ///
    /// `molblock` must run through an `M  END` line followed only by
    /// whitespace, with no `$$$$` line after the title.
    pub fn set_molblock(&mut self, molblock: &[u8]) -> Result<(), InvalidMolblock> {
        let title_end = line_end(molblock, 0);
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
            data.extend_from_slice(missing_ending(&data, title_newline(&data)));
        }
        let new_mdl_end = data.len();
        data.extend_from_slice(metadata);
        self.data = data;
        self.title_end = title_end;
        self.mdl_end = new_mdl_end;
        Ok(())
    }

    /// Give the first record named `name` the single-line `value`, keeping
    /// its header line, and remove later ones, or append the record, its
    /// header in `style`, if there is none.
    pub fn set_record(
        &mut self,
        name: &[u8],
        value: &[u8],
        strip: Strip,
        style: HeaderStyle,
    ) -> Result<(), MissingFieldName> {
        if self.rewrite_records(name, Some(value), strip)? == 0 {
            self.append_record(name, value, strip, style);
        }
        Ok(())
    }

    /// Remove every record named `name`, returning how many there were.
    pub fn remove_records(&mut self, name: &[u8], strip: Strip) -> Result<usize, MissingFieldName> {
        self.rewrite_records(name, None, strip)
    }

    /// Replace the whole SD data with `records`. The k-th record with a name
    /// keeps the k-th old header line with that name; other headers are in
    /// `style`, copying from the old SD data.
    pub fn replace_records<'a>(
        &mut self,
        records: impl IntoIterator<Item = (&'a [u8], &'a [u8])>,
        strip: Strip,
        style: HeaderStyle,
    ) {
        let mut records = records.into_iter().peekable();
        if records.peek().is_some() {
            self.ensure_title_line();
        }
        let (newline, trailer) = (self.newline(), self.trailer());
        let old = self.data.split_off(self.mdl_end);
        let mut headers = HeaderPool::new(&old[..old.len() - trailer.len()], strip, style);
        for (name, value) in records {
            match headers.take(name) {
                Some(line) => self.data.extend_from_slice(line),
                None => push_header(&mut self.data, headers.template, name, newline),
            }
            push_value(&mut self.data, value, newline);
        }
        self.data.extend_from_slice(trailer);
    }

    /// Drop the records named `name`, putting `value` under the first's
    /// header line, and return how many matched.
    fn rewrite_records(
        &mut self,
        name: &[u8],
        value: Option<&[u8]>,
        strip: Strip,
    ) -> Result<usize, MissingFieldName> {
        let metadata = self.metadata();
        let mut spans = Vec::new();
        for header in self.headers(strip) {
            let header = header?;
            if header.name == name {
                spans.push((header.start, header.end));
            }
        }
        if spans.is_empty() {
            return Ok(0);
        }
        let newline = self.newline();
        let mut data =
            Vec::with_capacity(self.data.len() + value.map_or(0, <[u8]>::len) + 2 * newline.len());
        data.extend_from_slice(&self.data[..self.mdl_end]);
        let mut copied = 0;
        for (index, &(start, end)) in spans.iter().enumerate() {
            data.extend_from_slice(&metadata[copied..start]);
            if let (0, Some(value)) = (index, value) {
                data.extend_from_slice(&metadata[start..line_end(metadata, start)]);
                push_value(&mut data, value, newline);
            }
            copied = end;
        }
        data.extend_from_slice(&self.data[self.mdl_end + copied..]);
        self.data = data;
        Ok(spans.len())
    }

    /// Replace the MDL lines, keeping the raw title line and SD data.
    pub fn replace_mdl<'a>(&mut self, mdl: impl IntoIterator<Item = &'a [u8]>) {
        let mut mdl = mdl.into_iter().peekable();
        if mdl.peek().is_some() {
            self.ensure_title_line();
        }
        let newline = self.newline();
        let mut data = Vec::with_capacity(self.data.len() + 64);
        data.extend_from_slice(&self.data[..self.title_end]);
        for line in mdl {
            data.extend_from_slice(line);
            if !line.ends_with(b"\n") {
                data.extend_from_slice(missing_ending(line, newline));
            }
        }
        let mdl_end = data.len();
        data.extend_from_slice(&self.data[self.mdl_end..]);
        self.data = data;
        self.mdl_end = mdl_end;
    }

    /// The records of the SD data at or after `pos`.
    pub fn records_at(&self, pos: usize, strip: Strip) -> Records<'_> {
        Records::new(self.metadata(), pos, strip)
    }

    /// The records of the SD data.
    pub fn records(&self, strip: Strip) -> Records<'_> {
        self.records_at(0, strip)
    }

    /// The record headers of the SD data, skipping the values.
    pub fn headers(&self, strip: Strip) -> Headers<'_> {
        self.records(strip).headers()
    }

    /// The format and atom and bond counts declared by the MDL header.
    ///
    /// Equivalent to parsing a `CTable` from the title and MDL lines, without
    /// copying them.
    pub fn counts(&self) -> Result<Counts, CTableError> {
        let mut source = std::iter::once(&b""[..])
            .chain(lines(self.mdl()))
            .map(String::from_utf8_lossy);
        let owned: [Option<Cow<'_, str>>; HEADER_LINES] = std::array::from_fn(|_| source.next());
        let found = owned.iter().take_while(|line| line.is_some()).count();
        let header = owned.each_ref().map(|line| line.as_deref().unwrap_or(""));
        parse_header(&header[..found])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::ctable::CTableFormat;
    use crate::core::ctable::error::Malformed;

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
        block
            .records(Strip::Bytes)
            .map(|record| {
                let record = record?;
                Ok((
                    String::from_utf8(record.name.to_vec()).unwrap(),
                    String::from_utf8(record.value.into_owned()).unwrap(),
                ))
            })
            .collect()
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
            "t\r\nsrc\nM  END\n> <A>\n a \nb\n\n> <B>\n\n$$$$\n"
        );
        let header_in_mdl = self::block("t\n> <A>\nM  END\n> <C>\nx\n$$$$\n", true);
        assert_eq!(records(&header_in_mdl).unwrap(), [("C".into(), "x".into())]);
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
    fn a_leading_delimiter_is_an_empty_record() {
        let block = block("$$$$ x\n", true);
        assert_eq!(block.title(Strip::Bytes), b"");
        assert_eq!(block.mdl(), b"");
        assert_eq!(block.molblock(Strip::Bytes).as_ref(), b"");
        assert_eq!(records(&block).unwrap(), []);
        assert_eq!(written(&block), "$$$$\n");
        let mut out = Vec::new();
        block.serialize_into(Strip::Bytes, &mut out);
        assert_eq!(out, b"$$$$\n");
        assert_eq!(written(&self::block("$$$$ x\r\n", true)), "$$$$\r\n");
        let Parts {
            data,
            title_end,
            mdl_end,
        } = block.parts();
        assert_eq!((title_end, mdl_end), (0, 0));
        assert_eq!(
            Block::from_parts(data.to_vec(), title_end, mdl_end),
            Some(block.clone())
        );
        assert_eq!(Block::from_parts(b"t\n$$$$\n".to_vec(), 0, 0), None);
        assert_eq!(Block::from_parts(b"$$$$\n".to_vec(), 0, 1), None);
    }

    #[test]
    fn lines_added_to_an_empty_record_follow_an_empty_title() {
        let empty = block("$$$$\r\n", true);
        let mut appended = empty.clone();
        appended.append_record(b"A", b"1", Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(written(&appended), "\r\n> <A>\r\n1\r\n\r\n$$$$\r\n");
        let mut set = empty.clone();
        set.set_record(b"A", b"1", Strip::Bytes, HeaderStyle::Plain)
            .unwrap();
        assert_eq!(set, appended);
        let mut replaced = empty.clone();
        replaced.replace_records([(&b"A"[..], &b"1"[..])], Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(replaced, appended);
        let mut cleared = empty.clone();
        cleared.replace_records([], Strip::Bytes, HeaderStyle::Plain);
        cleared.replace_mdl([]);
        assert_eq!(cleared, empty);
        let mut mdl = empty.clone();
        mdl.replace_mdl([&b"M  END"[..]]);
        assert_eq!(written(&mdl), "\r\nM  END\r\n$$$$\r\n");
        let mut titled = empty;
        titled.set_title(b"t");
        assert_eq!(written(&titled), "t\r\n$$$$\r\n");
    }

    #[test]
    fn crlf_blocks_keep_their_line_endings() {
        let tidy = block("t\r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$\r\n", true);
        assert!(matches!(tidy.serialized(Strip::Bytes), Cow::Borrowed(_)));
        assert!(matches!(
            tidy.molblock(Strip::Bytes),
            Cow::Borrowed(b"t\r\nM  END\r\n")
        ));
        assert_eq!(tidy.mdl(), b"M  END\r\n");

        // Lines the block writes itself end like the title line.
        let mut block = block(" t \r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$ x\r\n", true);
        assert_eq!(block.title(Strip::Bytes), b"t");
        assert_eq!(written(&block), "t\r\nM  END\r\n> <A>\r\n1\r\n\r\n$$$$\r\n");
        assert_eq!(
            written(&self::block("t\r\nx\r\n", false)),
            "t\r\nx\r\n$$$$\r\n"
        );
        // So does a last line without `\n`, unless a bare `\r` ends it.
        let unterminated = self::block("t\r\nM  END\r\n$$$$", true);
        assert_eq!(written(&unterminated), "t\r\nM  END\r\n$$$$\r\n");
        assert_eq!(written(&self::block("t\r\nx", false)), "t\r\nx\r\n$$$$\r\n");
        assert_eq!(
            written(&self::block("t\r\nx\r", false)),
            "t\r\nx\r\n$$$$\r\n"
        );
        assert_eq!(written(&self::block("t\r", false)), "t\r\n$$$$\r\n");
        assert_eq!(written(&self::block("t", false)), "t\n$$$$\n");
        block.append_record(b"B", b"2", Strip::Bytes, HeaderStyle::Plain);
        block.set_title(b"new");
        block
            .set_record(b"A", b"3", Strip::Bytes, HeaderStyle::Plain)
            .unwrap();
        let edited = "new\r\nM  END\r\n> <A>\r\n3\r\n\r\n> <B>\r\n2\r\n\r\n$$$$\r\n";
        assert_eq!(written(&block), edited);
        let Parts {
            data,
            title_end,
            mdl_end,
        } = block.parts();
        assert_eq!(
            Block::from_parts(data.to_vec(), title_end, mdl_end),
            Some(block.clone())
        );
        block.replace_records([(&b"C"[..], &b"c"[..])], Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(
            written(&block),
            "new\r\nM  END\r\n> <C>\r\nc\r\n\r\n$$$$\r\n"
        );

        // A delimiter line is kept as read, even if the title's ending differs.
        let mut mixed = self::block("t\nM  END\n$$$$\r\n", true);
        mixed.append_record(b"A", b"1", Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(written(&mixed), "t\nM  END\n> <A>\n1\n\n$$$$\r\n");
        assert_eq!(Block::from_parts(b"t\r\n$$$$\r".to_vec(), 3, 3), None);
        assert_eq!(Block::from_parts(b"\r\n".to_vec(), 2, 2), None);
    }

    #[test]
    fn unterminated_input_gets_a_newline() {
        let block = block("t\nM  END\n> <A>\nv", false);
        assert_eq!(records(&block).unwrap(), [("A".into(), "v".into())]);
        assert_eq!(written(&block), "t\nM  END\n> <A>\nv\n$$$$\n");
    }

    #[test]
    fn appended_records_follow_existing_metadata() {
        let mut block = block("t\nM  END\n> <A>\na\n\n$$$$\n", true);
        block.append_record(b"B", b"b", Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(
            records(&block).unwrap(),
            [("A".into(), "a".into()), ("B".into(), "b".into())]
        );
        assert_eq!(written(&block), "t\nM  END\n> <A>\na\n\n> <B>\nb\n\n$$$$\n");

        // SD data ending in a line that is not blank gets a blank line first.
        for (text, delimited) in [
            ("t\nM  END\n> <A>\na\n$$$$\n", true),
            ("t\nM  END\n> <A>\na", false),
        ] {
            let mut block = self::block(text, delimited);
            block.append_record(b"B", b"b", Strip::Bytes, HeaderStyle::Plain);
            assert_eq!(written(&block), "t\nM  END\n> <A>\na\n\n> <B>\nb\n\n$$$$\n");
        }
        let mut stray = self::block("t\r\nM  END\r\nstray\r\n$$$$\r\n", true);
        stray.append_record(b"B", b"b", Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(records(&stray).unwrap(), [("B".into(), "b".into())]);
        let mut spaced = self::block("t\nM  END\n> <A>\na\n \n$$$$\n", true);
        spaced.append_record(b"B", b"b", Strip::Bytes, HeaderStyle::Plain);
        assert_eq!(
            written(&spaced),
            "t\nM  END\n> <A>\na\n \n> <B>\nb\n\n$$$$\n"
        );
    }

    #[test]
    fn appended_headers_follow_the_style() {
        let text = "t\nM  END\n>  <A>  (1) \na\n\n";
        for (style, header) in [
            (HeaderStyle::Plain, "> <B>"),
            (HeaderStyle::Padded, ">  <B>"),
            (HeaderStyle::Copy, ">  <B>  (1) "),
        ] {
            let mut block = block(text, false);
            block.append_record(b"B", b"b", Strip::Bytes, style);
            assert_eq!(written(&block), format!("{text}{header}\nb\n\n$$$$\n"));
        }
    }

    #[test]
    fn copied_headers_follow_the_first_named_header() {
        for (text, header, newline) in [
            (
                "t\r\nM  END\r\n>  <A>  (1) \r\na\r\n\r\n",
                ">  <B>  (1) ",
                "\r\n",
            ),
            // The title's line ending, not the copied header's.
            ("t\nM  END\n > <A> (x)\r\na\n\n", " > <B> (x)", "\n"),
            // Not a line in a nameless header's run.
            (
                "t\nM  END\n> bad\n> <X> (v)\n\n> 7 <A> (M)\na\n\n",
                "> 7 <B> (M)",
                "\n",
            ),
            ("t\nM  END\n> <A\na\n\n> 7 <A> (M)\na\n\n", "> <B>", "\n"),
            // Opening with a blank line and a run that is not a record.
            ("t\nM  END\n\nnote\n\n> <A> x\na\n\n", "> <B> x", "\n"),
        ] {
            let mut block = block(text, false);
            block.append_record(b"B", b"b", Strip::Bytes, HeaderStyle::Copy);
            assert_eq!(
                written(&block),
                format!("{text}{header}{newline}b{newline}{newline}$$$${newline}")
            );
        }
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
        block
            .set_record(b"A", b"new", Strip::Bytes, HeaderStyle::Plain)
            .unwrap();
        assert_eq!(
            written(&block),
            "t\nM  END\n> <A>  (1)\nnew\n\njunk\n\n> <B>\nb\n\n$$$$\n"
        );
        block
            .set_record(b"C", b"c", Strip::Bytes, HeaderStyle::Plain)
            .unwrap();
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
            "t\nM  END\n> <A>  (1)\nnew\n\njunk\n\n> <C>\nc\n\n$$$$\n"
        );
        block.replace_records([(&b"Z"[..], &b"z"[..])], Strip::Bytes, HeaderStyle::Copy);
        assert_eq!(written(&block), "t\nM  END\n> <Z>  (1)\nz\n\n$$$$\n");
        let mut malformed = self::block("t\nM  END\n> bad\n\n", false);
        assert_eq!(
            malformed.remove_records(b"A", Strip::Bytes),
            Err(MissingFieldName)
        );
    }

    #[test]
    fn str_strip_ends_records_at_unicode_blank_lines() {
        let text = "\u{3000}t\u{a0}\nM  END\n> <A>\u{2003}\n\u{3000}\n\u{1c}tail\n\n";
        let block = block(text, false);
        assert_eq!(block.title(Strip::Str), b"t");
        let mut records = block.records(Strip::Str);
        let record = records.next().unwrap().unwrap();
        assert_eq!((record.name, &*record.value), (&b"A"[..], &b""[..]));
        assert_eq!(records.next(), None);
        assert_eq!(block.title(Strip::Bytes), "\u{3000}t\u{a0}".as_bytes());
        let record = block.records(Strip::Bytes).next().unwrap().unwrap();
        assert_eq!(record.name, b"A");
        assert_eq!(&*record.value, "\u{3000}\n\u{1c}tail".as_bytes());
    }

    #[test]
    fn counts_match_ctable_parsing() {
        let v2000 = block(
            "t\n\n\n  3  2  0  0  0  0            999 V2000\nM  END\n",
            false,
        );
        assert_eq!(
            v2000.counts(),
            Ok(Counts {
                format: CTableFormat::V2000,
                atoms: 3,
                bonds: 2
            })
        );
        let v3000 = block(
            "t\n\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 4 5 0 0 0\nM  END\n",
            false,
        );
        assert_eq!(
            v3000.counts(),
            Ok(Counts {
                format: CTableFormat::V3000,
                atoms: 4,
                bonds: 5
            })
        );
        assert_eq!(
            block("t\nM  END\n", false).counts(),
            Err(CTableError::Malformed {
                line: 3,
                problem: Malformed::Missing(2)
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
