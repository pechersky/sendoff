//! The native extension's pure-Rust core, compiled without PyO3 for fuzzing,
//! and invariants that every block keeps.

#[path = "../../rust/src/core/mod.rs"]
pub mod core;

use crate::core::block::Block;
use crate::core::ctable::CTable;
use crate::core::ctable::header::counts_line_format;
use crate::core::lines::lines;
use crate::core::record::Record;

/// Check what holds for any block, however it was read or edited.
pub fn check_block(block: &Block) {
    let parts = block.parts();
    let rebuilt = Block::from_parts(parts.data.to_vec(), parts.title_end, parts.mdl_end);
    assert_eq!(
        rebuilt.as_ref(),
        Some(block),
        "parts do not rebuild the block"
    );

    assert!(block.serialized().starts_with(block.molblock()));
    assert!(block.molblock().starts_with(block.title()));
    let records: Vec<_> = block.records().collect();
    let headers: Vec<_> = block.headers().collect();
    assert_eq!(records.len(), headers.len());
    for (record, header) in records.iter().zip(&headers) {
        match (record, header) {
            (Ok(record), Ok(header)) => {
                assert_eq!((record.name, record.start), (header.name, header.start));
                assert!(header.start < header.end && header.end <= block.metadata().len());
            }
            (Err(_), Err(_)) => {}
            _ => panic!("records and headers disagree: {record:?} {header:?}"),
        }
    }

    let mut source = vec![String::new()];
    source.extend(lines(block.mdl()).map(|line| String::from_utf8_lossy(line).into_owned()));
    let parsed =
        CTable::parse(source).map(|table| (table.format, table.num_atoms, table.num_bonds));
    let counts = block
        .counts()
        .map(|counts| (counts.format, counts.atoms, counts.bonds));
    assert_eq!(parsed, counts);
    assert!(block.mdl().ends_with(block.ctab()));
    let format = block.format();
    assert_eq!(
        format.and_then(|format| counts_line_format(format.name())),
        format
    );
    if let Ok((counts_format, ..)) = counts {
        assert_eq!(format, Some(counts_format));
    }
}

/// The named records, as `(name, value)` pairs.
pub fn named(block: &Block) -> Vec<(Vec<u8>, Vec<u8>)> {
    block
        .records()
        .filter_map(Result::ok)
        .map(|Record { name, value, .. }| (name.to_vec(), value.into_owned()))
        .collect()
}
