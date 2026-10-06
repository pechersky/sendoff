//! Record-name hashing and counting across blocks.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use crate::core::block::{BlockView, MissingFieldName, Strip};

/// A fast hash for short record names, eight bytes at a time.
#[derive(Default)]
pub struct NameHasher(u64);

impl Hasher for NameHasher {
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut word = [0; 8];
            word[..chunk.len()].copy_from_slice(chunk);
            self.0 = (self.0.rotate_left(5) ^ u64::from_le_bytes(word))
                .wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }

    fn finish(&self) -> u64 {
        // Fold the well-mixed high bits into the low bits buckets use.
        self.0 ^ (self.0 >> 32)
    }
}

/// A map keyed by record names.
pub type NameMap<V> = HashMap<Box<[u8]>, V, BuildHasherDefault<NameHasher>>;

/// Marks no name: the start of a block, or an unknown successor.
const NO_NAME: u32 = u32::MAX;

struct Counted {
    name: Box<[u8]>,
    blocks: usize,
    /// The last block counted, so repeats within a block count once.
    last_block: usize,
    /// The name that followed this one when last seen, or `NO_NAME`.
    next: u32,
}

/// How many blocks have each record name, in first-seen order.
///
/// Each name links to the name that last followed it, so blocks repeating
/// an earlier record order find their names by comparing bytes instead of
/// hashing them, even when some records are missing or new.
pub struct NameCounts {
    index: NameMap<u32>,
    names: Vec<Counted>,
    /// The first name of the last block.
    first: u32,
    blocks: usize,
}

impl Default for NameCounts {
    fn default() -> Self {
        Self {
            index: NameMap::default(),
            names: Vec::new(),
            first: NO_NAME,
            blocks: 0,
        }
    }
}

impl NameCounts {
    /// Count the record names of one more block.
    pub fn add(&mut self, block: &BlockView<'_>, strip: Strip) -> Result<(), MissingFieldName> {
        self.blocks += 1;
        let (mut pos, mut prev) = (0, NO_NAME);
        while let Some(name) = block.next_record_name(&mut pos, strip)? {
            let predicted = match prev {
                NO_NAME => self.first,
                prev => self.names[prev as usize].next,
            };
            let found = match self.names.get(predicted as usize) {
                Some(counted) if *counted.name == *name => predicted,
                _ => {
                    let found = self.find(name);
                    match prev {
                        NO_NAME => self.first = found,
                        prev => self.names[prev as usize].next = found,
                    }
                    found
                }
            };
            let counted = &mut self.names[found as usize];
            if counted.last_block != self.blocks {
                counted.last_block = self.blocks;
                counted.blocks += 1;
            }
            prev = found;
        }
        Ok(())
    }

    /// Return the index of `name`, adding it uncounted if new.
    fn find(&mut self, name: &[u8]) -> u32 {
        if let Some(&found) = self.index.get(name) {
            return found;
        }
        let found = u32::try_from(self.names.len()).expect("fewer than 2**32 record names");
        self.index.insert(name.into(), found);
        self.names.push(Counted {
            name: name.into(),
            blocks: 0,
            last_block: 0,
            next: NO_NAME,
        });
        found
    }

    /// The names and how many blocks have each, in first-seen order.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], usize)> {
        self.names
            .iter()
            .map(|counted| (&*counted.name, counted.blocks))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(blocks: &[&[u8]]) -> Vec<(String, usize)> {
        let mut counts = NameCounts::default();
        for data in blocks {
            counts
                .add(&BlockView::new(data, None), Strip::Bytes)
                .unwrap();
        }
        counts
            .iter()
            .map(|(name, blocks)| (String::from_utf8(name.to_vec()).unwrap(), blocks))
            .collect()
    }

    #[test]
    fn counts_blocks_per_name_in_first_seen_order() {
        let a_b_a: &[u8] = b"t\nM  END\n> <A>\n1\n\n> <B>\n\n> <A>\n2\n\n";
        let b_c: &[u8] = b"t\nM  END\n> <B>\n1\n\n> <C>\n\n";
        let none: &[u8] = b"t\nM  END\n";
        assert_eq!(
            counts(&[a_b_a, b_c, none, a_b_a, b_c]),
            [("A".into(), 2), ("B".into(), 4), ("C".into(), 2)]
        );
        assert_eq!(counts(&[none]), []);
    }

    #[test]
    fn counts_reports_malformed_headers() {
        let mut counts = NameCounts::default();
        let malformed = BlockView::new(b"t\nM  END\n> no name\n1\n\n", None);
        assert!(counts.add(&malformed, Strip::Bytes).is_err());
    }
}
