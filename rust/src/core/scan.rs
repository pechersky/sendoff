//! Split a byte stream into molecule blocks at `$$$$` lines.
//!
//! The scanner owns a buffer that its caller refills on demand, so the same
//! search serves file reads, Python `read()` chunks, and Python line
//! iterables. Each block is copied out once into its own `Block`.

use std::sync::LazyLock;

use memchr::{memchr, memmem};

use crate::core::block::Block;

static DELIMITER: LazyLock<memmem::Finder<'static>> =
    LazyLock::new(|| memmem::Finder::new(b"\n$$$$"));

/// Incremental block splitter over a caller-filled buffer.
#[derive(Default)]
pub struct Scanner {
    buf: Vec<u8>,
    /// Start of the next block.
    pos: usize,
    /// Where to resume searching for `\n$$$$`.
    search_from: usize,
    /// Start of the current block's delimiter line, once found.
    delimiter: Option<usize>,
    eof: bool,
}

impl Scanner {
    /// Return the next block, calling `fill` to append more input as needed.
    ///
    /// `fill` appends any number of bytes and returns `true` at end of
    /// input. A block ends at the first line starting with `$$$$`, which may
    /// be its first line; trailing input without a delimiter is a final block.
    pub fn next_block<E>(
        &mut self,
        mut fill: impl FnMut(&mut Vec<u8>) -> Result<bool, E>,
    ) -> Result<Option<Block>, E> {
        loop {
            if let Some(block) = self.take_block() {
                return Ok(Some(block));
            }
            if self.eof {
                return Ok(None);
            }
            self.compact();
            self.eof = fill(&mut self.buf)?;
        }
    }

    fn take_block(&mut self) -> Option<Block> {
        let available = self.buf.len() - self.pos;
        if available == 0 || (available < 4 && !self.eof) {
            return None;
        }
        let start = match self.delimiter {
            Some(start) => start,
            None if self.buf[self.pos..].starts_with(b"$$$$") => self.pos,
            None => match DELIMITER.find(&self.buf[self.search_from..]) {
                Some(index) => self.search_from + index + 1,
                None if self.eof => return Some(self.emit(self.buf.len(), None)),
                None => {
                    self.search_from = self.buf.len().saturating_sub(4).max(self.pos);
                    return None;
                }
            },
        };
        self.delimiter = Some(start);
        let end = match memchr(b'\n', &self.buf[start..]) {
            Some(index) => start + index + 1,
            None if self.eof => self.buf.len(),
            None => return None,
        };
        Some(self.emit(end, Some(start)))
    }

    fn emit(&mut self, end: usize, delimiter: Option<usize>) -> Block {
        let data = self.buf[self.pos..end].to_vec();
        let delimiter = delimiter.map(|start| start - self.pos);
        self.pos = end;
        self.search_from = end;
        self.delimiter = None;
        Block::new(data, delimiter)
    }

    fn compact(&mut self) {
        if self.pos == 0 {
            return;
        }
        self.buf.drain(..self.pos);
        self.search_from -= self.pos;
        if let Some(start) = &mut self.delimiter {
            *start -= self.pos;
        }
        self.pos = 0;
    }
}

/// Translate `\r\n` and lone `\r` to `\n` in `buf[start..]`, in place.
///
/// `pending_cr` carries a chunk-final `\r` so that a `\n` starting the next
/// chunk is dropped, matching Python's universal newlines.
pub fn translate_newlines(buf: &mut Vec<u8>, start: usize, pending_cr: &mut bool) {
    let len = buf.len();
    if start == len {
        return;
    }
    let mut read = start;
    let mut write = start;
    if std::mem::take(pending_cr) && buf[read] == b'\n' {
        read += 1;
    }
    while let Some(offset) = memchr(b'\r', &buf[read..len]) {
        let cr = read + offset;
        if read != write {
            buf.copy_within(read..cr, write);
        }
        write += cr - read;
        buf[write] = b'\n';
        write += 1;
        read = cr + 1;
        if read == len {
            *pending_cr = true;
        } else if buf[read] == b'\n' {
            read += 1;
        }
    }
    if read != write {
        buf.copy_within(read..len, write);
    }
    buf.truncate(write + (len - read));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn split(input: &[u8], chunk: usize) -> Vec<(String, String)> {
        let mut scanner = Scanner::default();
        let mut offset = 0;
        let mut blocks = Vec::new();
        while let Some(block) = scanner
            .next_block(|buf| {
                let end = (offset + chunk).min(input.len());
                buf.extend_from_slice(&input[offset..end]);
                offset = end;
                Ok::<_, ()>(offset == input.len())
            })
            .unwrap()
        {
            let mut written = Vec::new();
            block.serialize_into(&mut written);
            blocks.push((
                String::from_utf8(block.title().to_vec()).unwrap(),
                String::from_utf8(written).unwrap(),
            ));
        }
        blocks
    }

    #[test]
    fn splits_at_delimiter_lines_for_every_chunk_size() {
        let input = b"a\nM  END\n> <X>\n1\n\n$$$$\n$$$$x\nb\nc\n$$$$\nd\nM  END";
        let expected = split(input, input.len());
        assert_eq!(
            expected,
            [
                ("a".into(), "a\nM  END\n> <X>\n1\n\n$$$$\n".into()),
                ("$$$$x".into(), "$$$$x\n$$$$\n".into()),
                ("b".into(), "b\nc\n$$$$\n$$$$\n".into()),
                ("d".into(), "d\nM  END\n$$$$\n".into()),
            ]
        );
        for chunk in 1..input.len() {
            assert_eq!(split(input, chunk), expected, "chunk {chunk}");
        }
    }

    #[test]
    fn empty_and_blank_input() {
        assert!(split(b"", 1).is_empty());
        assert_eq!(split(b"\n", 1), [(String::new(), "\n$$$$\n".into())]);
        assert_eq!(split(b"$$$", 1), [("$$$".into(), "$$$\n$$$$\n".into())]);
    }

    fn translated(input: &[u8], chunk: usize) -> Vec<u8> {
        let mut buf = Vec::new();
        let mut pending_cr = false;
        for piece in input.chunks(chunk) {
            let start = buf.len();
            buf.extend_from_slice(piece);
            translate_newlines(&mut buf, start, &mut pending_cr);
        }
        buf
    }

    #[test]
    fn translates_universal_newlines_across_chunks() {
        let input = b"a\r\nb\rc\n\r\r\nd\r";
        for chunk in 1..=input.len() {
            assert_eq!(
                translated(input, chunk),
                b"a\nb\nc\n\n\nd\n",
                "chunk {chunk}"
            );
        }
        assert_eq!(translated(b"plain\n", 3), b"plain\n");
    }
}
