//! Line endings: which ending lines a block writes itself get, the `$$$$`
//! line that ends a serialized block, and translating universal newlines
//! as Python reads text.

use std::borrow::Cow;

use memchr::memchr;

use crate::core::lines::{WINDOW, byte_mask};

/// The delimiter lines that end serialized blocks, as read.
pub const TRAILER: &[u8] = b"$$$$\n";
pub const CRLF_TRAILER: &[u8] = b"$$$$\r\n";

/// `\r\n` if `line` ends with it, else `\n`.
pub fn newline(line: &[u8]) -> &'static [u8] {
    if line.ends_with(b"\r\n") {
        b"\r\n"
    } else {
        b"\n"
    }
}

/// The ending of `data`'s first line once it ends in `\n`, which lines a
/// block writes itself use.
pub fn title_newline(data: &[u8]) -> &'static [u8] {
    match memchr(b'\n', data) {
        Some(index) => newline(&data[..=index]),
        None if data.ends_with(b"\r") => b"\r\n",
        None => b"\n",
    }
}

/// What ends `line` if it lacks a `\n`: `\n` after a `\r`, else `newline`.
pub fn missing_ending(line: &[u8], newline: &'static [u8]) -> &'static [u8] {
    if line.ends_with(b"\r") {
        b"\n"
    } else {
        newline
    }
}

/// The `$$$$` line ending in `newline`.
pub fn trailer(newline: &[u8]) -> &'static [u8] {
    if newline == b"\r\n" {
        CRLF_TRAILER
    } else {
        TRAILER
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
    // Skip the `\r`-free prefix with vectorized `memchr` before entering the
    // windowed translator: LF-only inputs return immediately, and blocks with
    // LF molblocks and CRLF SD records jump straight to the first `\r`.
    let mut base = if read == 0 && space[0] != b'\r' {
        match memchr(b'\r', &space[..len]) {
            Some(first) => {
                read = first;
                write = first;
                first
            }
            None => return len,
        }
    } else {
        0
    };
    // Scan each window once for `\r` instead of calling `memchr` per `\r`:
    // CRLF files have a `\r` every line, where per-call setup dominates.
    // Writes stay below `read`, so bytes at or after `read` are unmodified,
    // and every `\r` a window's mask finds is at or after `read`.
    let lane_end = space.len().saturating_sub(LANE);
    while base < len {
        let mut mask = byte_mask(&space[base..len], b'\r');
        while mask != 0 {
            let cr = base + mask.trailing_zeros() as usize;
            mask &= mask - 1;
            debug_assert!(cr >= read);
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

/// `bytes` with `\r\n` and lone `\r` translated to `\n`, as Python reads
/// text.
pub fn translated(bytes: &[u8]) -> Vec<u8> {
    match memchr(b'\r', bytes) {
        Some(from) => translated_from(bytes, from),
        None => bytes.to_vec(),
    }
}

/// `translated`, for `bytes` known to have no `\r` before `from`, which it
/// copies without scanning.
pub fn translated_from(bytes: &[u8], from: usize) -> Vec<u8> {
    if from == bytes.len() {
        return bytes.to_vec();
    }
    let mut buffer = Vec::with_capacity(bytes.len() + TRANSLATE_SLACK);
    buffer.extend_from_slice(bytes);
    buffer.resize(bytes.len() + TRANSLATE_SLACK, 0);
    let len = from + translate_newlines(&mut buffer[from..], bytes.len() - from, &mut false);
    buffer.truncate(len);
    buffer
}

/// Like `translated`, but borrow `bytes` when it has no `\r`.
pub fn universal_newlines(bytes: &[u8]) -> Cow<'_, [u8]> {
    match memchr(b'\r', bytes) {
        Some(from) => Cow::Owned(translated_from(bytes, from)),
        None => Cow::Borrowed(bytes),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::testing::xorshift;

    fn chunked(input: &[u8], chunk: usize, slack: usize) -> Vec<u8> {
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
                    chunked(input, chunk, slack),
                    b"a\nb\nc\n\n\nd\n",
                    "chunk {chunk}"
                );
            }
        }
        assert_eq!(chunked(b"plain\n", 3, 0), b"plain\n");
    }

    #[test]
    fn translation_matches_replace_for_random_input_and_chunks() {
        let mut next = xorshift(0x9E37_79B9_7F4A_7C15_u64);
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
                let output = chunked(&input, chunk, slack);
                assert_eq!(output, expected.as_bytes(), "{input:?}");
            }
            let clean = memchr(b'\r', &input).unwrap_or(len);
            let from = (next() as usize) % (clean + 1);
            assert_eq!(translated_from(&input, from), expected.as_bytes());
        }
    }
}
