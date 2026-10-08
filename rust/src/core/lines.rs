//! Byte search a 64-byte window at a time, and line iteration on top of it.
//!
//! SD data and CRLF files have a match every few bytes, where the setup of
//! a search call per match costs more than the search: finding every match
//! in a window at once and stepping through a bitmask is cheaper.

/// Bytes covered by one `byte_mask`.
pub const WINDOW: usize = 64;

/// Bit `i` is set when `window[i]` is `byte`; `window` may be short.
#[inline]
pub fn byte_mask(window: &[u8], byte: u8) -> u64 {
    match <&[u8; WINDOW]>::try_from(window) {
        Ok(full) => full_mask(full, byte),
        Err(_) => {
            let mut padded = [0; WINDOW];
            padded[..window.len()].copy_from_slice(window);
            full_mask(&padded, byte) & low_bits(window.len())
        }
    }
}

/// The lowest `len` bits, for `len <= WINDOW`.
#[inline]
fn low_bits(len: usize) -> u64 {
    u64::MAX.checked_shr((WINDOW - len) as u32).unwrap_or(0)
}

#[inline]
fn full_mask(window: &[u8; WINDOW], byte: u8) -> u64 {
    #[cfg(target_arch = "x86_64")]
    {
        sse2_mask(window, byte)
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        swar_mask(window, byte)
    }
}

/// SSE2 is part of the x86_64 baseline, so no runtime detection is needed.
#[cfg(target_arch = "x86_64")]
#[inline]
fn sse2_mask(window: &[u8; WINDOW], byte: u8) -> u64 {
    use std::arch::x86_64::{_mm_cmpeq_epi8, _mm_loadu_si128, _mm_movemask_epi8, _mm_set1_epi8};
    window
        .chunks_exact(16)
        .enumerate()
        .fold(0, |mask, (index, lane)| {
            // SAFETY: SSE2 is enabled on every x86_64 target, and `lane` is
            // 16 readable bytes for the unaligned load.
            let bits = unsafe {
                let lane = _mm_loadu_si128(lane.as_ptr().cast());
                _mm_movemask_epi8(_mm_cmpeq_epi8(lane, _mm_set1_epi8(byte as i8)))
            } as u16;
            mask | (u64::from(bits) << (16 * index))
        })
}

/// Portable fallback: exact per-byte zero test on `word ^ byte`, 8 at a time.
#[cfg(any(test, not(target_arch = "x86_64")))]
fn swar_mask(window: &[u8; WINDOW], byte: u8) -> u64 {
    const LOW: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    let splat = u64::from_ne_bytes([byte; 8]);
    window
        .chunks_exact(8)
        .enumerate()
        .fold(0, |mask, (index, word)| {
            let word = u64::from_le_bytes(word.try_into().unwrap()) ^ splat;
            // Bit 7 of each byte is set exactly when that byte is zero.
            let zero = !(((word & LOW) + LOW) | word | LOW);
            // Gather those bits, byte `i`'s into bit `56 + i`.
            let bits = ((zero >> 7).wrapping_mul(0x0102_0408_1020_4080)) >> 56;
            mask | (bits << (8 * index))
        })
}

/// The lines of `data` from a position, each with its `\n` except perhaps
/// the last.
pub struct Lines<'a> {
    data: &'a [u8],
    pos: usize,
    /// Start of the window that `newlines` covers.
    base: usize,
    /// The window's newlines at or after `pos`, bit `i` for `data[base + i]`.
    newlines: u64,
}

impl<'a> Lines<'a> {
    /// Iterate over the lines of `data` starting at `pos`.
    #[inline]
    pub fn new(data: &'a [u8], pos: usize) -> Self {
        let newlines = data
            .get(pos..)
            .map_or(0, |rest| byte_mask(&rest[..rest.len().min(WINDOW)], b'\n'));
        Self {
            data,
            pos,
            base: pos,
            newlines,
        }
    }

    /// Where the next line starts.
    #[inline]
    pub fn pos(&self) -> usize {
        self.pos
    }
}

/// The lines of `data` without their `\n`.
pub fn lines(data: &[u8]) -> impl Iterator<Item = &[u8]> {
    Lines::new(data, 0).map(|line| line.strip_suffix(b"\n").unwrap_or(line))
}

/// Where the line starting at `start` ends, past its `\n`, or `data.len()`.
#[inline]
pub fn line_end(data: &[u8], start: usize) -> usize {
    find_newline(&data[start..]).map_or(data.len(), |index| start + index + 1)
}

/// The first `\n` in `data`. Not memchr: for titles, `M  END`, and `$$$$`
/// lines, a few bytes long, its out-of-line call costs more than this inlined
/// loop, which made block-level operations 1-2% faster on EC2.
#[inline]
pub fn find_newline(data: &[u8]) -> Option<usize> {
    data.iter().position(|&byte| byte == b'\n')
}

impl<'a> Iterator for Lines<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<&'a [u8]> {
        if self.pos >= self.data.len() {
            return None;
        }
        while self.newlines == 0 {
            self.base += WINDOW;
            let Some(rest) = self.data.get(self.base..).filter(|rest| !rest.is_empty()) else {
                let line = &self.data[self.pos..];
                self.pos = self.data.len();
                return Some(line);
            };
            self.newlines = byte_mask(&rest[..rest.len().min(WINDOW)], b'\n');
        }
        let end = self.base + self.newlines.trailing_zeros() as usize + 1;
        self.newlines &= self.newlines - 1;
        let line = &self.data[self.pos..end];
        self.pos = end;
        Some(line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::testing::xorshift;

    fn random_bytes(len: usize, seed: u64, alphabet: &[u8]) -> Vec<u8> {
        let mut next = xorshift(seed);
        (0..len)
            .map(|_| alphabet[(next() % alphabet.len() as u64) as usize])
            .collect()
    }

    #[test]
    fn masks_agree() {
        for seed in 1..50 {
            let data = random_bytes(WINDOW, seed, b"\n\rx\x8a\x8d\0 ");
            let window: &[u8; WINDOW] = data.as_slice().try_into().unwrap();
            for byte in [b'\n', b'\r', 0, b' ', 0x8a, b'y'] {
                let expected = data.iter().enumerate().fold(0, |mask, (index, &found)| {
                    mask | (u64::from(found == byte) << index)
                });
                assert_eq!(swar_mask(window, byte), expected);
                assert_eq!(byte_mask(window, byte), expected);
                for len in [0, 1, 10, 63] {
                    assert_eq!(byte_mask(&data[..len], byte), expected & low_bits(len));
                }
            }
        }
    }

    #[test]
    fn lines_match_split_inclusive_from_any_position() {
        for (seed, len) in (1..40).zip([0, 1, 2, 63, 64, 65, 127, 128, 129, 300].iter().cycle()) {
            let data = random_bytes(*len, seed, b"\n\nab \r");
            for pos in 0..=data.len() {
                let expected: Vec<&[u8]> =
                    data[pos..].split_inclusive(|&byte| byte == b'\n').collect();
                let mut lines = Lines::new(&data, pos);
                let mut found = Vec::new();
                let mut start = lines.pos();
                while let Some(line) = lines.next() {
                    assert_eq!(&data[start..lines.pos()], line);
                    start = lines.pos();
                    found.push(line);
                }
                assert_eq!(found, expected, "len {len} pos {pos}");
                assert_eq!(lines.pos(), data.len());
                assert_eq!(lines.next(), None);
            }
        }
    }
}
