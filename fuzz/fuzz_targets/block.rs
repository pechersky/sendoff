//! Edit sequences keep every block invariant, and with names and values
//! that read back unchanged, leave the records they promise.

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use sendoff_fuzz::core::block::Block;
use sendoff_fuzz::core::record::HeaderStyle;
use sendoff_fuzz::core::split::split_block;
use sendoff_fuzz::core::strip::Strip;
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
fn clean_value(value: &[u8], strip: Strip) -> bool {
    !value.contains(&b'\n') && !value.contains(&b'\r') && {
        let trimmed = strip.trim(value);
        !trimmed.is_empty() && trimmed.len() == value.len()
    }
}

fn clean((name, value): (&[u8], &[u8]), strip: Strip) -> bool {
    clean_name(name) && clean_value(value, strip)
}

fn apply(block: &mut Block, edit: &Edit<'_>, strip: Strip) {
    let before = block.clone();
    let records = |block: &Block| named(block, strip);
    let owned = |name: &[u8], value: &[u8]| (name.to_vec(), value.to_vec());
    match *edit {
        Edit::Append(name, value, style) => {
            block.append_record(name, value, strip, style.into());
            if clean((name, value), strip) {
                assert_eq!(records(block).last(), Some(&owned(name, value)));
            }
        }
        Edit::Set(name, value, style) => match block.set_record(name, value, strip, style.into()) {
            Ok(()) if clean((name, value), strip) => {
                let found: Vec<_> = records(block)
                    .into_iter()
                    .filter(|(n, _)| n == name)
                    .collect();
                assert_eq!(found, [owned(name, value)]);
            }
            Ok(()) => {}
            Err(_) => assert_eq!(*block, before),
        },
        Edit::Remove(name) => match block.remove_records(name, strip) {
            Ok(removed) => {
                let count =
                    |block: &Block| records(block).iter().filter(|(n, _)| n == name).count();
                assert_eq!((removed, count(block)), (count(&before), 0));
            }
            Err(_) => assert_eq!(*block, before),
        },
        Edit::Replace(ref new, style) => {
            block.replace_records(new.iter().copied(), strip, style.into());
            if new.iter().all(|&record| clean(record, strip)) {
                let expected: Vec<_> = new
                    .iter()
                    .map(|&(name, value)| owned(name, value))
                    .collect();
                assert_eq!(records(block), expected);
            }
        }
        Edit::Title(title) => {
            if !title.contains(&b'\n') && !title.contains(&b'\r') {
                block.set_title(title);
                assert_eq!(block.title(strip), strip.trim(title));
            }
        }
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

fuzz_target!(|input: (bool, bool, Vec<Edit<'_>>, &[u8])| {
    let (text, str_strip, edits, data) = input;
    if data.is_empty() {
        return;
    }
    let strip = if str_strip { Strip::Str } else { Strip::Bytes };
    let mut block = split_block(data, text).0;
    check_block(&block);
    for edit in edits.iter().take(8) {
        apply(&mut block, edit, strip);
        check_block(&block);
    }
});
