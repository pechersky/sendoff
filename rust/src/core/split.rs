//! Split a byte stream into molecule blocks at `$$$$` lines.
//!
//! The scanner owns a buffer that its caller refills on demand, so the same
//! search serves file reads, Python `read()` chunks, and Python line
//! iterables. Each block is copied out once into its own `Block`. The buffer
//! stays initialized up to its high-water mark, so refills write into it
//! without zero-filling first. `split_block` splits input that is already
//! complete, such as a mapped file, without a buffer, and `borrow_block`
//! splits a mapped file into blocks that borrow it.

use std::sync::{Arc, LazyLock};

use memchr::{memchr, memchr2, memmem, memrchr};
use memmap2::Mmap;

use crate::core::block::Block;
use crate::core::lines::{find_newline, line_end};
use crate::core::newline::translated_from;

static DELIMITER: LazyLock<memmem::Finder<'static>> =
    LazyLock::new(|| memmem::Finder::new(b"\n$$$$"));
static DOLLARS: LazyLock<memmem::Finder<'static>> = LazyLock::new(|| memmem::Finder::new(b"$$$$"));

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
        let end = match find_newline(&buf[start..]) {
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

/// Where the first block in complete input ends, and where its `$$$$` line
/// starts, if it has one.
fn bounds(data: &[u8]) -> (usize, Option<usize>) {
    let delimiter = if data.starts_with(b"$$$$") {
        Some(0)
    } else {
        DELIMITER.find(data).map(|index| index + 1)
    };
    let end = delimiter.map_or(data.len(), |start| line_end(data, start));
    (end, delimiter)
}

/// Where the first line starting with `$$$$` starts, with universal
/// newlines: a line starts after `\n`, or after a `\r` that `$` follows.
fn text_delimiter(data: &[u8]) -> Option<usize> {
    if data.starts_with(b"$$$$") {
        return Some(0);
    }
    let mut from = 1;
    while let Some(index) = DOLLARS.find(data.get(from..)?) {
        let start = from + index;
        if matches!(data[start - 1], b'\n' | b'\r') {
            return Some(start);
        }
        from = start + 1;
    }
    None
}

/// Where the first block in complete input ends, where its `$$$$` line
/// starts (if it has one), and how many leading bytes need no newline
/// translation; `text` reads universal newlines.
fn extent(data: &[u8], text: bool) -> (usize, Option<usize>, usize) {
    if !text {
        let (end, delimiter) = bounds(data);
        return (end, delimiter, end);
    }
    let delimiter = text_delimiter(data);
    let end = delimiter
        .and_then(|start| memchr2(b'\n', b'\r', &data[start..]).map(|index| start + index))
        .map_or(data.len(), |at| {
            at + 1 + usize::from(data[at..].starts_with(b"\r\n"))
        });
    let clean = memchr(b'\r', &data[..end]).unwrap_or(end);
    (end, delimiter, clean)
}

/// Split the first block from complete, nonempty input, returning it and
/// the input length it took.
///
/// `text` reads universal newlines, translating them as `translated` does,
/// so blocks match a `Scanner` over translated input.
pub fn split_block(data: &[u8], text: bool) -> (Block, usize) {
    let (end, delimiter, clean) = extent(data, text);
    (copied(&data[..end], delimiter, clean), end)
}

/// Split the first block from `map[pos..]`, which is nonempty, as
/// `split_block` does, returning it and the input length it took.
///
/// The block borrows the map when its bytes are already a block's: ending
/// in a `$$$$` line, and, for `text`, without `\r` to translate.
pub fn borrow_block(map: &Arc<Mmap>, pos: usize, text: bool) -> (Block, usize) {
    let data = &map[pos..];
    let (end, delimiter, clean) = extent(data, text);
    if let Some(start) = delimiter
        && clean == end
        && let Some(block) = Block::borrowed(map, pos, pos + end, start)
    {
        return (block, end);
    }
    (copied(&data[..end], delimiter, clean), end)
}

/// Copy one block's input, whose `$$$$` line, if any, starts at
/// `delimiter`, translating universal newlines from `clean`, as no `\r`
/// comes before it.
fn copied(data: &[u8], delimiter: Option<usize>, clean: usize) -> Block {
    if clean == data.len() {
        return Block::new(data.to_vec(), delimiter);
    }
    let block = translated_from(data, clean);
    // The delimiter line is the last line, found in the copy so that a file
    // changing under a map cannot set it past the block.
    let delimiter = delimiter.map(|_| {
        let lines = block.strip_suffix(b"\n").unwrap_or(&block);
        memrchr(b'\n', lines).map_or(0, |index| index + 1)
    });
    Block::new(block, delimiter)
}

/// Parse `data` as one block ending at its first `$$$$` line, if any.
///
/// Returns `None` when anything but ASCII whitespace follows that line.
pub fn single_block(mut data: Vec<u8>) -> Option<Block> {
    let (end, delimiter) = bounds(&data);
    if !data[end..].iter().all(u8::is_ascii_whitespace) {
        return None;
    }
    data.truncate(end);
    Some(Block::new(data, delimiter))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::newline::translated;
    use crate::core::testing::xorshift;

    fn scanned(input: &[u8], chunk: usize) -> Vec<Block> {
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
            blocks.push(block);
        }
        blocks
    }

    fn split(input: &[u8], chunk: usize) -> Vec<(String, String)> {
        scanned(input, chunk)
            .iter()
            .map(|block| {
                (
                    String::from_utf8(block.title().to_vec()).unwrap(),
                    String::from_utf8(block.serialized().to_vec()).unwrap(),
                )
            })
            .collect()
    }

    fn split_complete(input: &[u8], text: bool) -> Vec<Block> {
        let mut pos = 0;
        let mut blocks = Vec::new();
        while pos < input.len() {
            let (block, taken) = split_block(&input[pos..], text);
            assert!(taken > 0);
            blocks.push(block);
            pos += taken;
        }
        blocks
    }

    #[test]
    fn split_block_matches_the_scanner_on_random_input() {
        const TOKENS: [&[u8]; 10] = [
            b"\n",
            b"\r",
            b"\r\n",
            b"$$$$",
            b"$",
            b"M  END",
            b"> <X>",
            b"x",
            b" ",
            b"\n$$$$\r\n",
        ];
        let mut next = xorshift(0x2545_F491_4F6C_DD1D_u64);
        for _ in 0..3000 {
            let len = (next() % 60) as usize;
            let input: Vec<u8> = (0..len)
                .flat_map(|_| TOKENS[(next() % TOKENS.len() as u64) as usize])
                .copied()
                .collect();
            let chunk = 1 + (next() % 50) as usize;
            assert_eq!(
                split_complete(&input, false),
                scanned(&input, chunk),
                "{input:?}"
            );
            let translated = translated(&input);
            assert_eq!(
                split_complete(&input, true),
                scanned(&translated, chunk),
                "{input:?}"
            );
        }
    }

    #[test]
    fn split_block_reads_universal_newline_delimiters() {
        let blocks = split_complete(b"a\r$$$$\rb\r\n$$$$\r\nc\n$$$$x\r$$$$", true);
        let serialized: Vec<_> = blocks
            .iter()
            .map(|block| block.serialized().to_vec())
            .collect();
        assert_eq!(
            serialized,
            [&b"a\n$$$$\n"[..], b"b\n$$$$\n", b"c\n$$$$\n", b"$$$$\n"]
        );
        assert_eq!(split_complete(b"a\r$$$$\rb", false).len(), 1);
    }

    #[test]
    fn splits_at_delimiter_lines_for_every_chunk_size() {
        let input = b"a\nM  END\n> <X>\n1\n\n$$$$\n$$$$x\nb\nc\n$$$$\nd\nM  END";
        let expected = split(input, input.len());
        assert_eq!(
            expected,
            [
                ("a".into(), "a\nM  END\n> <X>\n1\n\n$$$$\n".into()),
                (String::new(), "$$$$\n".into()),
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
        let serialized =
            |input: &[u8]| single_block(input.to_vec()).map(|block| block.serialized().to_vec());
        assert_eq!(
            serialized(b"t\nM  END"),
            Some(b"t\nM  END\n$$$$\n".to_vec())
        );
        assert_eq!(
            serialized(b"t\nM  END\n$$$$\n \n\t"),
            Some(b"t\nM  END\n$$$$\n".to_vec())
        );
        assert_eq!(serialized(b"$$$$"), Some(b"$$$$\n".to_vec()));
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
}
