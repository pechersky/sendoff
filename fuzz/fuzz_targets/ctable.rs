//! CTable parsing, index validation, renumbering, and count fixing on raw
//! lines and on generated V3000 tables with arbitrary indices and noise.

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use sendoff_fuzz::core::ctable::{CTable, Section};

#[derive(Arbitrary, Debug)]
struct Atom<'a> {
    index: i8,
    element: &'a str,
    continued: bool,
}

#[derive(Arbitrary, Debug)]
struct Bond {
    index: i8,
    atoms: (i8, i8),
    continued: bool,
}

#[derive(Arbitrary, Debug)]
enum Input<'a> {
    Raw(Vec<String>),
    V3000 {
        counts: (i8, i8),
        atoms: Vec<Atom<'a>>,
        bonds: Vec<Bond>,
        noise: Vec<(u8, &'a str)>,
    },
}

/// A V3000 line, split over a `-` continuation when `continued`.
fn push(lines: &mut Vec<String>, line: String, continued: bool) {
    match line.split_once(' ').filter(|_| continued) {
        Some((first, rest)) => {
            lines.push(format!("M  V30 {first} -"));
            lines.push(format!("M  V30 {rest}"));
        }
        None => lines.push(format!("M  V30 {line}")),
    }
}

fn lines(input: Input<'_>) -> Vec<String> {
    let (counts, atoms, bonds, noise) = match input {
        Input::Raw(lines) => return lines,
        Input::V3000 {
            counts,
            atoms,
            bonds,
            noise,
        } => (counts, atoms, bonds, noise),
    };
    let mut lines: Vec<String> = [
        "title",
        "  program",
        "",
        "  0  0  0     0  0            999 V3000",
    ]
    .map(String::from)
    .into();
    lines.push("M  V30 BEGIN CTAB".into());
    lines.push(format!("M  V30 COUNTS {} {} 0 0 0", counts.0, counts.1));
    lines.push("M  V30 BEGIN ATOM".into());
    for Atom {
        index,
        element,
        continued,
    } in atoms
    {
        push(&mut lines, format!("{index} {element} 0 0 0 0"), continued);
    }
    lines.push("M  V30 END ATOM".into());
    lines.push("M  V30 BEGIN BOND".into());
    for Bond {
        index,
        atoms: (a, b),
        continued,
    } in bonds
    {
        push(&mut lines, format!("{index} 1 {a} {b}"), continued);
    }
    lines.push("M  V30 END BOND".into());
    lines.push("M  V30 END CTAB".into());
    lines.push("M  END".into());
    for (at, line) in noise {
        lines.insert(usize::from(at) % (lines.len() + 1), line.into());
    }
    lines
}

fn validate(table: &CTable) {
    for section in [Section::Atom, Section::Bond] {
        for strict in [false, true] {
            let _ = table.valid_indices(section, strict);
        }
    }
}

fuzz_target!(|input: Input<'_>| {
    let lines = lines(input);
    let Ok(mut table) = CTable::parse(lines.clone()) else {
        return;
    };
    validate(&table);
    if table.renumber_indices().is_err() {
        assert_eq!(table.lines, lines);
    }
    validate(&table);
    let before = table.lines.clone();
    if table.fix_counts().is_ok() {
        let reparsed = CTable::parse(table.lines.clone()).unwrap();
        assert_eq!(
            (reparsed.num_atoms, reparsed.num_bonds),
            (table.num_atoms, table.num_bonds)
        );
    } else {
        assert_eq!(table.lines, before);
    }
    validate(&table);
});
