//! Split a byte stream into molecule blocks at `$$$$` lines.
//!
//! The scanner owns a buffer that its caller refills on demand, so the same
//! search serves file reads, Python `read()` chunks, and Python line
//! iterables. Each block is copied out once into its own `Block`. The buffer
//! stays initialized up to its high-water mark, so refills write into it
//! without zero-filling first.

use std::sync::LazyLock;

use memchr::{memchr, memmem};

use crate::core::block::Block;
use crate::core::lines::{WINDOW, byte_mask};

static DELIMITER: LazyLock<memmem::Finder<'static>> =
    LazyLock::new(|| memmem::Finder::new(b"\n$$$$"));

/// Incremental block splitter over a caller-filled buffer.
#[derive(Default)]
pub struct Scanner {
    /// Input is `buf[..end]`; bytes after `end` are stale.
    buf: Vec<u8>,
    end: usize,
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
        mut fill: impl FnMut(&mut Tail<'_>) -> Result<bool, E>,
    ) -> Result<Option<Block>, E> {
        loop {
            if let Some(block) = self.take_block() {
                return Ok(Some(block));
            }
            if self.eof {
                return Ok(None);
            }
            self.compact();
            self.eof = fill(&mut Tail {
                buf: &mut self.buf,
                end: &mut self.end,
            })?;
        }
    }

    fn take_block(&mut self) -> Option<Block> {
        let buf = &self.buf[..self.end];
        let available = buf.len() - self.pos;
        if available == 0 || (available < 4 && !self.eof) {
            return None;
        }
        let start = match self.delimiter {
            Some(start) => start,
            None if buf[self.pos..].starts_with(b"$$$$") => self.pos,
            None => match DELIMITER.find(&buf[self.search_from..]) {
                Some(index) => self.search_from + index + 1,
                None if self.eof => return Some(self.emit(buf.len(), None)),
                None => {
                    self.search_from = buf.len().saturating_sub(4).max(self.pos);
                    return None;
                }
            },
        };
        self.delimiter = Some(start);
        let end = match memchr(b'\n', &buf[start..]) {
            Some(index) => start + index + 1,
            None if self.eof => buf.len(),
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
        self.buf.copy_within(self.pos..self.end, 0);
        self.end -= self.pos;
        self.search_from -= self.pos;
        if let Some(start) = &mut self.delimiter {
            *start -= self.pos;
        }
        self.pos = 0;
    }
}

/// Parse `data` as one block ending at its first `$$$$` line, if any.
///
/// Returns `None` when anything but ASCII whitespace follows that line.
pub fn single_block(mut data: Vec<u8>) -> Option<Block> {
    let delimiter = if data.starts_with(b"$$$$") {
        Some(0)
    } else {
        DELIMITER.find(&data).map(|index| index + 1)
    };
    if let Some(start) = delimiter {
        let end = memchr(b'\n', &data[start..]).map_or(data.len(), |index| start + index + 1);
        if !data[end..].iter().all(u8::is_ascii_whitespace) {
            return None;
        }
        data.truncate(end);
    }
    Some(Block::new(data, delimiter))
}

/// Writable space after a scanner's input.
pub struct Tail<'a> {
    buf: &'a mut Vec<u8>,
    end: &'a mut usize,
}

impl Tail<'_> {
    /// Append `bytes` to the input.
    pub fn extend(&mut self, bytes: &[u8]) {
        self.space(bytes.len())[..bytes.len()].copy_from_slice(bytes);
        self.commit(bytes.len());
    }

    /// Return at least `len` writable bytes with unspecified contents.
    pub fn space(&mut self, len: usize) -> &mut [u8] {
        let needed = *self.end + len;
        if self.buf.len() < needed {
            self.buf.resize(needed, 0);
        }
        &mut self.buf[*self.end..]
    }

    /// Append the first `len` bytes of `space` to the input.
    pub fn commit(&mut self, len: usize) {
        assert!(*self.end + len <= self.buf.len());
        *self.end += len;
    }
}

/// Spare bytes after the input that let `translate_newlines` move whole lanes.
pub const TRANSLATE_SLACK: usize = LANE;

const LANE: usize = 16;

/// Translate `\r\n` and lone `\r` in `space[..len]` to `\n` in place, and
/// return the translated length.
///
/// `pending_cr` carries a chunk-final `\r` so that a `\n` starting the next
/// chunk is dropped, matching Python's universal newlines. Translation is
/// fastest with `TRANSLATE_SLACK` bytes of `space` after the input.
pub fn translate_newlines(space: &mut [u8], len: usize, pending_cr: &mut bool) -> usize {
    if len == 0 {
        return 0;
    }
    let mut read = 0;
    let mut write = 0;
    if std::mem::take(pending_cr) && space[read] == b'\n' {
        read += 1;
    }
    // Scan each window once for `\r` instead of calling `memchr` per `\r`:
    // CRLF files have a `\r` every line, where per-call setup dominates.
    // Writes stay below `read`, so bytes at or after `read` are unmodified.
    let lane_end = space.len().saturating_sub(LANE);
    let mut base = 0;
    while base < len {
        let mut mask = byte_mask(&space[base..len.min(base + WINDOW)], b'\r');
        while mask != 0 {
            let cr = base + mask.trailing_zeros() as usize;
            mask &= mask - 1;
            if cr < read {
                continue;
            }
            if read - write >= LANE && cr <= lane_end {
                // Whole-lane copies are safe once writes trail reads by a
                // lane: they only overwrite consumed input or bytes that
                // later output replaces.
                let (mut from, mut to) = (read, write);
                while from < cr {
                    let lane: [u8; LANE] = space[from..from + LANE].try_into().unwrap();
                    space[to..to + LANE].copy_from_slice(&lane);
                    from += LANE;
                    to += LANE;
                }
            } else if read != write {
                space.copy_within(read..cr, write);
            }
            write += cr - read;
            space[write] = b'\n';
            write += 1;
            read = cr + 1;
            if read == len {
                *pending_cr = true;
            } else if space[read] == b'\n' {
                read += 1;
            }
        }
        base += WINDOW;
    }
    if read != write {
        space.copy_within(read..len, write);
    }
    write + (len - read)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::block::Strip;

    fn split(input: &[u8], chunk: usize) -> Vec<(String, String)> {
        let mut scanner = Scanner::default();
        let mut offset = 0;
        let mut blocks = Vec::new();
        while let Some(block) = scanner
            .next_block(|tail| {
                let end = (offset + chunk).min(input.len());
                tail.extend(&input[offset..end]);
                offset = end;
                Ok::<_, ()>(offset == input.len())
            })
            .unwrap()
        {
            blocks.push((
                String::from_utf8(block.title(Strip::Bytes).to_vec()).unwrap(),
                String::from_utf8(block.serialized(Strip::Bytes).into_owned()).unwrap(),
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
                ("b".into(), "b\nc\n$$$$\n".into()),
                ("d".into(), "d\nM  END\n$$$$\n".into()),
            ]
        );
        for chunk in 1..input.len() {
            assert_eq!(split(input, chunk), expected, "chunk {chunk}");
        }
    }

    #[test]
    fn single_block_allows_only_whitespace_after_its_delimiter() {
        let serialized = |input: &[u8]| {
            single_block(input.to_vec()).map(|block| block.serialized(Strip::Bytes).into_owned())
        };
        assert_eq!(
            serialized(b"t\nM  END"),
            Some(b"t\nM  END\n$$$$\n".to_vec())
        );
        assert_eq!(
            serialized(b"t\nM  END\n$$$$\n \n\t"),
            Some(b"t\nM  END\n$$$$\n".to_vec())
        );
        assert_eq!(serialized(b"$$$$"), Some(b"$$$$\n$$$$\n".to_vec()));
        assert_eq!(serialized(b""), Some(b"\n$$$$\n".to_vec()));
        assert_eq!(serialized(b"t\n$$$$\nu\n"), None);
        assert_eq!(serialized(b"$$$$\n$$$$\n"), None);
    }

    #[test]
    fn empty_and_blank_input() {
        assert!(split(b"", 1).is_empty());
        assert_eq!(split(b"\n", 1), [(String::new(), "\n$$$$\n".into())]);
        assert_eq!(split(b"$$$", 1), [("$$$".into(), "$$$\n$$$$\n".into())]);
    }

    fn translated(input: &[u8], chunk: usize, slack: usize) -> Vec<u8> {
        let mut out = Vec::new();
        let mut pending_cr = false;
        for piece in input.chunks(chunk) {
            let mut space = piece.to_vec();
            space.resize(piece.len() + slack, b'\r');
            let len = translate_newlines(&mut space, piece.len(), &mut pending_cr);
            out.extend_from_slice(&space[..len]);
        }
        out
    }

    #[test]
    fn translates_universal_newlines_across_chunks() {
        let input = b"a\r\nb\rc\n\r\r\nd\r";
        for chunk in 1..=input.len() {
            for slack in [0, TRANSLATE_SLACK] {
                assert_eq!(
                    translated(input, chunk, slack),
                    b"a\nb\nc\n\n\nd\n",
                    "chunk {chunk}"
                );
            }
        }
        assert_eq!(translated(b"plain\n", 3, 0), b"plain\n");
    }

    #[test]
    fn translation_matches_replace_for_random_input_and_chunks() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..1000 {
            // Dense and sparse line endings exercise short and multi-lane runs.
            let sparseness = [4, 40][(next() % 2) as usize];
            let len = (next() % 2000) as usize;
            let input: Vec<u8> = (0..len)
                .map(|_| match next() % sparseness {
                    0 => b'\r',
                    1 => b'\n',
                    _ => b'x',
                })
                .collect();
            let expected = String::from_utf8(input.clone())
                .unwrap()
                .replace("\r\n", "\n")
                .replace('\r', "\n");
            let chunk = 1 + (next() % 700) as usize;
            for slack in [0, TRANSLATE_SLACK] {
                let output = translated(&input, chunk, slack);
                assert_eq!(output, expected.as_bytes(), "{input:?}");
            }
        }
    }
}
