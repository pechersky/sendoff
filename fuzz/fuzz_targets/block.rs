//! Edit sequences keep every block invariant, and with names and values
//! that read back unchanged, leave the records they promise.

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use sendoff_fuzz::core::block::Block;
use sendoff_fuzz::core::record::HeaderStyle;
use sendoff_fuzz::core::split::split_block;
use sendoff_fuzz::core::strip::{trim, trim_end};
use sendoff_fuzz::{check_block, named};

#[derive(Arbitrary, Debug, Clone, Copy)]
enum Style {
    Plain,
    Padded,
    Copy,
}

impl From<Style> for HeaderStyle {
    fn from(style: Style) -> Self {
        match style {
            Style::Plain => HeaderStyle::Plain,
            Style::Padded => HeaderStyle::Padded,
            Style::Copy => HeaderStyle::Copy,
        }
    }
}

#[derive(Arbitrary, Debug)]
enum Edit<'a> {
    Append(&'a [u8], &'a [u8], Style),
    Set(&'a [u8], &'a [u8], Style),
    Remove(&'a [u8]),
    Replace(Vec<(&'a [u8], &'a [u8])>, Style),
    Title(&'a [u8]),
    Molblock(&'a [u8]),
    Mdl(Vec<&'a [u8]>),
}

/// A name that a header line reads back unchanged.
fn clean_name(name: &[u8]) -> bool {
    !name.is_empty()
        && name
            .iter()
            .all(|byte| byte.is_ascii_graphic() && !b"<>".contains(byte))
}

/// A single-line value that its record reads back unchanged.
fn clean_value(value: &[u8]) -> bool {
    !value.contains(&b'\n') && !value.ends_with(b"\r") && !trim(value).is_empty()
}

fn clean((name, value): (&[u8], &[u8])) -> bool {
    clean_name(name) && clean_value(value)
}

fn apply(block: &mut Block, edit: &Edit<'_>) {
    let before = block.clone();
    let owned = |name: &[u8], value: &[u8]| (name.to_vec(), value.to_vec());
    match *edit {
        Edit::Append(name, value, style) => {
            block.append_record(name, value, style.into());
            if clean((name, value)) {
                assert_eq!(named(block).last(), Some(&owned(name, value)));
            }
        }
        Edit::Set(name, value, style) => match block.set_record(name, value, style.into()) {
            Ok(()) if clean((name, value)) => {
                let found: Vec<_> = named(block)
                    .into_iter()
                    .filter(|(n, _)| n == name)
                    .collect();
                assert_eq!(found, [owned(name, value)]);
            }
            Ok(()) => {}
            Err(_) => assert_eq!(*block, before),
        },
        Edit::Remove(name) => match block.remove_records(name) {
            Ok(removed) => {
                let count = |block: &Block| named(block).iter().filter(|(n, _)| n == name).count();
                assert_eq!((removed, count(block)), (count(&before), 0));
            }
            Err(_) => assert_eq!(*block, before),
        },
        Edit::Replace(ref new, style) => {
            block.replace_records(new.iter().copied(), style.into());
            if new.iter().all(|&record| clean(record)) {
                let expected: Vec<_> = new
                    .iter()
                    .map(|&(name, value)| owned(name, value))
                    .collect();
                assert_eq!(named(block), expected);
            }
        }
        Edit::Title(title) => match block.set_title(title) {
            Ok(()) => assert_eq!(block.title(), trim_end(title)),
            Err(_) => assert_eq!(*block, before),
        },
        Edit::Molblock(molblock) => {
            if block.set_molblock(molblock).is_err() {
                assert_eq!(*block, before);
            }
        }
        Edit::Mdl(ref mdl) => {
            if mdl.iter().all(|line| !line.contains(&b'\n')) {
                block.replace_mdl(mdl.iter().copied());
            }
        }
    }
}

fuzz_target!(|input: (bool, Vec<Edit<'_>>, &[u8])| {
    let (text, edits, data) = input;
    if data.is_empty() {
        return;
    }
    let mut block = split_block(data, text).0;
    check_block(&block);
    for edit in edits.iter().take(8) {
        apply(&mut block, edit);
        check_block(&block);
    }
});
