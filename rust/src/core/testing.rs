//! Helpers for unit tests.

use crate::core::block::Block;
use crate::core::record::HeaderStyle;

/// A xorshift64 generator of deterministic pseudo-random test input.
pub fn xorshift(seed: u64) -> impl FnMut() -> u64 {
    let mut state = seed | 1;
    move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    }
}

/// Block edits, each returning a summary of its result to compare.
pub const EDITS: [fn(&mut Block) -> String; 8] = [
    |block| format!("{:?}", block.set_title(b"M  END")),
    |block| format!("{:?}", block.set_title(b"")),
    |block| {
        block.append_record(b"A", b"M  END", HeaderStyle::Copy);
        String::new()
    },
    |block| {
        format!(
            "{:?}",
            block.set_record(b"B", b"M  END", HeaderStyle::Plain)
        )
    },
    |block| format!("{:?}", block.remove_records(b"A")),
    |block| {
        block.replace_records([(&b"A"[..], &b"M  END"[..])], HeaderStyle::Copy);
        String::new()
    },
    |block| {
        block.replace_mdl([&b"M  END"[..], b"> <B>"]);
        String::new()
    },
    |block| format!("{:?}", block.set_molblock(b"t\nM  END\n")),
];
