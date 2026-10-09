//! Splitting complete input matches the scanner over arbitrary chunks, as
//! Python reads them, and a serialized block reads back as itself.

#![no_main]

use std::convert::Infallible;

use libfuzzer_sys::fuzz_target;
use sendoff_fuzz::check_block;
use sendoff_fuzz::core::block::Block;
use sendoff_fuzz::core::newline::{TRANSLATE_SLACK, translate_newlines};
use sendoff_fuzz::core::split::{Scanner, Tail, split_block};

/// Blocks from a scanner filled in chunks of the cut lengths, translating
/// newlines per chunk in text mode as the Python reader does.
fn scanned(data: &[u8], text: bool, cuts: &[u16]) -> Vec<Block> {
    let (mut scanner, mut blocks) = (Scanner::default(), Vec::new());
    let (mut rest, mut cuts, mut pending_cr) = (data, cuts.iter(), false);
    let mut fill = |tail: &mut Tail<'_>| {
        let len = cuts
            .next()
            .map_or(rest.len(), |&cut| usize::from(cut) % (rest.len() + 1));
        let (chunk, after) = rest.split_at(len);
        if text {
            let space = tail.space(len + TRANSLATE_SLACK);
            space[..len].copy_from_slice(chunk);
            let translated = translate_newlines(space, len, &mut pending_cr);
            tail.commit(translated);
        } else {
            tail.extend(chunk);
        }
        rest = after;
        Ok::<_, Infallible>(rest.is_empty())
    };
    while let Some(block) = scanner.next_block(&mut fill).unwrap() {
        blocks.push(block);
    }
    blocks
}

fuzz_target!(|input: (bool, Vec<u16>, &[u8])| {
    let (text, cuts, data) = input;
    let (mut blocks, mut rest) = (Vec::new(), data);
    while !rest.is_empty() {
        let (block, taken) = split_block(rest, text);
        assert!(0 < taken && taken <= rest.len());
        blocks.push(block);
        rest = &rest[taken..];
    }
    assert_eq!(scanned(data, text, &cuts), blocks);

    for block in &blocks {
        check_block(block);
        let once = block.serialized();
        let (again, taken) = split_block(once, false);
        assert_eq!(taken, once.len());
        assert_eq!(again.serialized(), once);
    }
});
