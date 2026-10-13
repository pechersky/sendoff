//! An index of SD data records, for repeated lookups by name.

use std::borrow::Cow;

use crate::core::record::{MissingFieldName, Records, record_value};
use crate::core::strip::Strip;

/// Where each record of some SD data lies, so lookups by name need not
/// parse it again.
///
/// Offsets are `u32`, half the size of `usize`, so SD data over 4 GiB goes
/// unindexed.
pub struct Index(Box<[Entry]>);

struct Entry {
    name: u32,
    name_end: u32,
    value: u32,
    value_end: u32,
    crlf: bool,
}

impl Entry {
    fn name<'a>(&self, metadata: &'a [u8]) -> &'a [u8] {
        &metadata[self.name as usize..self.name_end as usize]
    }
}

impl Index {
    /// Index the records of `metadata`, unless it is too large.
    ///
    /// # Errors
    ///
    /// `MissingFieldName` if a record header has no `<name>`.
    pub fn build(metadata: &[u8], strip: Strip) -> Result<Option<Self>, MissingFieldName> {
        if u32::try_from(metadata.len()).is_err() {
            return Ok(None);
        }
        // The length fits in u32, so do the offsets.
        let offset = |at: usize| at as u32;
        let entries = Records::new(metadata, 0, strip).spans().map(|span| {
            let span = span?;
            let name = span.name.as_ptr().addr() - metadata.as_ptr().addr();
            Ok(Entry {
                name: offset(name),
                name_end: offset(name + span.name.len()),
                value: offset(span.value.start),
                value_end: offset(span.value.end),
                crlf: span.crlf,
            })
        });
        Ok(Some(Self(entries.collect::<Result<_, _>>()?)))
    }

    /// The number of records, duplicates included.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Whether there are no records.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Whether any record of `metadata`, which the index was built from, is
    /// named `name`.
    pub fn contains(&self, metadata: &[u8], name: &[u8]) -> bool {
        self.0.iter().any(|entry| entry.name(metadata) == name)
    }

    /// The value of the last record of `metadata`, which the index was
    /// built from, named `name`.
    pub fn last<'a>(&self, metadata: &'a [u8], name: &[u8]) -> Option<Cow<'a, [u8]>> {
        let entry = self.0.iter().rfind(|entry| entry.name(metadata) == name)?;
        let value = entry.value as usize..entry.value_end as usize;
        Some(record_value(metadata, value, entry.crlf))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;

    #[test]
    fn lookups_match_records() {
        for metadata in [
            "> <A>\n1\n2\n\n> <B>\n\njunk\nmore\n\n> <A>\n3\n",
            "\n\n> <A> (1)\n\u{3000}\n> <B>\n",
            "> <A>\n1\r\n 2\r\r\n\n> <B>\nc\r\n \n> <C>\n",
            "> <a<b>c>\nv\n>  < x > (1)\nw\n\n> <>\n\n",
            "",
        ] {
            for strip in [Strip::Bytes, Strip::Str] {
                let data = metadata.as_bytes();
                let index = Index::build(data, strip).unwrap().unwrap();
                let records: Vec<_> = Records::new(data, 0, strip).map(Result::unwrap).collect();
                assert_eq!(index.len(), records.len());
                assert_eq!(index.is_empty(), records.is_empty());
                let mut last = HashMap::new();
                for record in &records {
                    last.insert(record.name, record.value.clone());
                }
                for (name, value) in &last {
                    assert!(index.contains(data, name));
                    assert_eq!(index.last(data, name).as_ref(), Some(value));
                }
                assert!(!index.contains(data, b"missing"));
                assert_eq!(index.last(data, b"missing"), None);
            }
        }
    }

    #[test]
    fn a_header_without_a_name_fails_the_index() {
        let metadata = b"> <A>\n1\n\n> missing\nv\n\n> <B>\n2\n";
        assert_eq!(
            Index::build(metadata, Strip::Bytes).err(),
            Some(MissingFieldName)
        );
    }
}
