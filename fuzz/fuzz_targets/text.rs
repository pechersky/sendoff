//! Byte masks, line iteration, and newline translation against naive
//! versions, translating in arbitrary chunks over stale spare bytes.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sendoff_fuzz::core::lines::{Lines, WINDOW, byte_mask, find_newline};
use sendoff_fuzz::core::newline::{
    TRANSLATE_SLACK, translate_newlines, translated, universal_newlines,
};

fn naive_translated(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len());
    let mut bytes = data.iter().copied().peekable();
    while let Some(byte) = bytes.next() {
        if byte == b'\r' {
            bytes.next_if_eq(&b'\n');
            out.push(b'\n');
        } else {
            out.push(byte);
        }
    }
    out
}

fuzz_target!(|input: (u8, u16, Vec<u16>, &[u8])| {
    let (byte, pos, cuts, data) = input;
    let pos = usize::from(pos) % (data.len() + 1);

    for window in data[pos..].chunks(WINDOW) {
        let expected = window.iter().enumerate().fold(0, |mask, (index, &found)| {
            mask | (u64::from(found == byte) << index)
        });
        assert_eq!(byte_mask(window, byte), expected);
    }

    let rest = &data[pos..];
    assert_eq!(find_newline(rest), rest.iter().position(|&b| b == b'\n'));
    let mut lines = Lines::new(data, pos);
    let mut start = pos;
    for expected in rest.split_inclusive(|&b| b == b'\n') {
        assert_eq!(lines.next(), Some(expected));
        assert_eq!(&data[start..lines.pos()], expected);
        start = lines.pos();
    }
    assert_eq!(lines.next(), None);
    assert_eq!(lines.pos(), data.len());

    let expected = naive_translated(data);
    assert_eq!(translated(data), expected);
    assert_eq!(*universal_newlines(data), *expected);

    // Each cut picks a chunk length and how many stale `\r` bytes follow it.
    let (mut out, mut pending, mut rest) = (Vec::new(), false, data);
    let mut cuts = cuts.into_iter();
    while !rest.is_empty() {
        let cut = cuts.next().unwrap_or(u16::MAX);
        let len = usize::from(cut & 0xff) % rest.len() + 1;
        let slack = usize::from(cut >> 8) % (2 * TRANSLATE_SLACK + 1);
        let (chunk, tail) = rest.split_at(len);
        let mut space = chunk.to_vec();
        space.resize(len + slack, b'\r');
        let translated = translate_newlines(&mut space, len, &mut pending);
        out.extend_from_slice(&space[..translated]);
        rest = tail;
    }
    assert_eq!(out, expected);
});
