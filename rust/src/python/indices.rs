use super::exact_text;
use crate::core::{
    ctable,
    indices::{self, Index},
    whitespace,
};
use pyo3::prelude::*;
use pyo3::{
    exceptions::{PyIndexError, PyNotImplementedError},
    types::{PyDict, PyInt, PyList, PySlice, PyString},
};
use std::collections::HashMap;

enum Tokens<'py> {
    Text(Vec<String>),
    Protocol(Bound<'py, PyAny>),
}

impl<'py> Tokens<'py> {
    fn split(line: &Bound<'py, PyAny>) -> PyResult<Self> {
        if let Some(text) = exact_text(line)? {
            Ok(Self::Text(
                text.split(whitespace)
                    .filter(|token| !token.is_empty())
                    .map(str::to_owned)
                    .collect(),
            ))
        } else {
            Ok(Self::Protocol(line.call_method0("split")?))
        }
    }

    fn get(&self, py: Python<'py>, position: usize) -> PyResult<Bound<'py, PyAny>> {
        match self {
            Self::Text(tokens) => tokens
                .get(position)
                .map(|token| PyString::new(py, token).into_any())
                .ok_or_else(|| PyIndexError::new_err("list index out of range")),
            Self::Protocol(tokens) => tokens.get_item(position),
        }
    }

    fn index(&self, py: Python<'py>, position: usize) -> PyResult<Index> {
        if let Self::Text(tokens) = self
            && let Some(token) = tokens.get(position)
            && let Some(index) = indices::parse_index(token)
        {
            return Ok(index);
        }
        // Python handles decimal Unicode, underscores, digit limits and int protocols.
        let value = py
            .import("builtins")?
            .getattr("int")?
            .call1((self.get(py, position)?,))?;
        if let Ok(value) = value.extract::<i128>() {
            return Ok(Index::Small(value));
        }
        let size = value.call_method0("bit_length")?.extract::<usize>()? / 8 + 1;
        let keywords = PyDict::new(py);
        keywords.set_item("signed", true)?;
        Ok(Index::Large(
            value
                .call_method("to_bytes", (size, "little"), Some(&keywords))?
                .extract()?,
        ))
    }
}

fn domain_error(errors: &Bound<'_, PyAny>, name: &str, message: &str) -> PyResult<PyErr> {
    Ok(PyErr::from_value(errors.getattr(name)?.call1((message,))?))
}

fn validation_error(
    errors: &Bound<'_, PyAny>,
    failure: indices::ValidationFailure,
    label: &str,
    singular: &str,
) -> PyResult<PyErr> {
    let (name, message) = match failure {
        indices::ValidationFailure::OutOfOrder => ("IndicesOutOfOrderError", label.to_owned()),
        indices::ValidationFailure::Duplicate => ("IndicesDuplicateError", label.to_owned()),
        indices::ValidationFailure::Fewer => (
            "IndicesMismatchError",
            format!("fewer {singular} lines than count line"),
        ),
        indices::ValidationFailure::More => (
            "IndicesMismatchError",
            format!("more {singular} lines than count line"),
        ),
    };
    domain_error(errors, name, &message)
}

fn plain_index(line: &Bound<'_, PyAny>, position: usize) -> PyResult<Option<Index>> {
    if let Some(text) = exact_text(line)?
        && let Some(token) = ctable::token(text, position)
        && let Some(index) = indices::parse_index(token)
    {
        return Ok(Some(index));
    }
    Ok(None)
}

fn exact_small_int(value: &Bound<'_, PyAny>) -> Option<i128> {
    value
        .is_exact_instance_of::<PyInt>()
        .then(|| value.extract().ok())
        .flatten()
}

fn validate(
    py: Python<'_>,
    table: &Bound<'_, PyAny>,
    strict: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    errors: &Bound<'_, PyAny>,
    atoms: bool,
) -> PyResult<bool> {
    if !table.getattr("format")?.is(v3000) {
        return Err(PyNotImplementedError::new_err(()));
    }
    let mut validator = indices::Validator::new();
    let (method, count, label, singular) = if atoms {
        ("atomlines", "num_atoms", "atoms", "atom")
    } else {
        ("bondlines", "num_bonds", "bonds", "bond")
    };
    let mut current_line = None;
    let mut current_tokens = None;
    for line in table.call_method0(method)?.try_iter()? {
        let line = current_line.insert(line?);
        let index = match plain_index(line, 2)? {
            Some(index) => index,
            None => current_tokens.insert(Tokens::split(line)?).index(py, 2)?,
        };
        let strict = strict.is_truthy()?;
        if let Err(failure) = validator.accept(index, strict) {
            return Err(validation_error(errors, failure, label, singular)?);
        }
    }
    let expected = table.getattr(count)?;
    if let Some(expected) = exact_small_int(&expected) {
        if let Err(failure) = validator.finish(expected) {
            return Err(validation_error(errors, failure, label, singular)?);
        }
    } else {
        let actual = validator.len().into_pyobject(py)?;
        if actual.lt(&expected)? {
            return Err(domain_error(
                errors,
                "IndicesMismatchError",
                &format!("fewer {singular} lines than count line"),
            )?);
        }
        if actual.gt(&expected)? {
            return Err(domain_error(
                errors,
                "IndicesMismatchError",
                &format!("more {singular} lines than count line"),
            )?);
        }
    }
    drop(current_line);
    drop(current_tokens);
    Ok(true)
}

#[pyfunction]
fn valid_atom_indices(
    py: Python<'_>,
    table: &Bound<'_, PyAny>,
    strict: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    errors: &Bound<'_, PyAny>,
) -> PyResult<bool> {
    validate(py, table, strict, v3000, errors, true)
}

#[pyfunction]
fn valid_bond_indices(
    py: Python<'_>,
    table: &Bound<'_, PyAny>,
    strict: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    errors: &Bound<'_, PyAny>,
) -> PyResult<bool> {
    validate(py, table, strict, v3000, errors, false)
}

fn startswith(line: &Bound<'_, PyAny>, prefix: &str) -> PyResult<bool> {
    match exact_text(line)? {
        Some(text) => Ok(text.starts_with(prefix)),
        None => line.call_method1("startswith", (prefix,))?.is_truthy(),
    }
}

fn collect_plain_lines<'py>(iterator: Bound<'py, PyAny>) -> PyResult<Option<Vec<String>>> {
    let mut lines = Vec::new();
    for line in iterator.try_iter()? {
        let line = line?;
        let Some(text) = exact_text(&line)? else {
            return Ok(None);
        };
        lines.push(text.to_owned());
    }
    Ok(Some(lines))
}

fn plain_table_lines(table: &Bound<'_, PyAny>) -> PyResult<Option<Vec<String>>> {
    collect_plain_lines(table.getattr("lines")?)
}

fn section_lines(lines: &[String]) -> (Vec<String>, Vec<String>) {
    let atoms = lines
        .iter()
        .skip(7)
        .take_while(|line| !line.starts_with("M  V30 END ATOM"))
        .cloned()
        .collect();
    let bonds = lines
        .iter()
        .skip_while(|line| !line.starts_with("M  V30 BEGIN BOND"))
        .skip(1)
        .take_while(|line| !line.starts_with("M  V30 END BOND"))
        .cloned()
        .collect();
    (atoms, bonds)
}

fn plain_indices(lines: &[String], positions: &[usize]) -> bool {
    positions.iter().all(|&position| {
        lines.iter().all(|line| {
            ctable::token(line, position).is_some_and(|token| indices::parse_index(token).is_some())
        })
    })
}

fn fast_renumber(
    py: Python<'_>,
    table: &Bound<'_, PyAny>,
    duplicate_error: &Bound<'_, PyAny>,
) -> PyResult<Option<Vec<String>>> {
    let Some(all_lines) = plain_table_lines(table)? else {
        return Ok(None);
    };
    let (plain_atoms, plain_bonds) = section_lines(&all_lines);
    if !plain_indices(&plain_atoms, &[2]) || !plain_indices(&plain_bonds, &[3, 4, 5]) {
        return Ok(None);
    }
    let atom_lines = table.call_method0("atomlines")?;
    let bond_lines = table.call_method0("bondlines")?;
    let Some(atom_lines) = collect_plain_lines(atom_lines)? else {
        return Ok(None);
    };
    let Some(bond_lines) = collect_plain_lines(bond_lines)? else {
        return Ok(None);
    };
    let counts = table.getattr("counts")?;
    let Some(counts) = exact_text(&counts)? else {
        return Ok(None);
    };
    match indices::renumber(&all_lines, &atom_lines, &bond_lines, counts) {
        Ok(result) => Ok(Some(result.lines)),
        Err(indices::RenumberFailure::DuplicateMapping) => Err(PyErr::from_value(
            duplicate_error.call1(("atom index mapping in bond",))?,
        )),
        Err(indices::RenumberFailure::MissingEndpoint) => {
            Err(PyIndexError::new_err("list index out of range"))
        }
        Err(indices::RenumberFailure::EmptyCounts) => {
            Err(PyIndexError::new_err("string index out of range"))
        }
        Err(
            indices::RenumberFailure::InvalidAtomIndex | indices::RenumberFailure::InvalidBondIndex,
        ) => {
            let _ = py;
            Ok(None)
        }
    }
}

fn renumber_compat(
    py: Python<'_>,
    table: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    duplicate_error: &Bound<'_, PyAny>,
) -> PyResult<()> {
    if !table.getattr("format")?.is(v3000) {
        return Err(PyNotImplementedError::new_err(()));
    }
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut mapping: HashMap<Index, Vec<usize>> = HashMap::new();
    let mut current_atomline = None;
    let mut current_tokens = None;
    let mut current_bondline = None;
    for (position, line) in table.call_method0("atomlines")?.try_iter()?.enumerate() {
        let line = current_atomline.insert(line?);
        let tokens = current_tokens.insert(Tokens::split(line)?);
        let index = tokens.index(py, 2)?;
        let new_index = position + 1;
        mapping.entry(index).or_default().push(new_index);
        let skip = 7 + match tokens {
            Tokens::Text(tokens) => tokens[2].chars().count(),
            Tokens::Protocol(_) => tokens.get(py, 2)?.len()?,
        };
        let prefix = format!("M  V30 {new_index}");
        let replacement = if let Some(text) = exact_text(line)? {
            let offset = text
                .char_indices()
                .nth(skip)
                .map_or(text.len(), |(offset, _)| offset);
            PyString::new(py, &(prefix + &text[offset..])).into_any()
        } else {
            PyString::new(py, &prefix).add(
                line.get_item(
                    py.get_type::<PySlice>()
                        .call1((skip, py.None(), py.None()))?,
                )?,
            )?
        };
        atoms.push(replacement);
    }
    for (position, line) in table.call_method0("bondlines")?.try_iter()?.enumerate() {
        let line = current_bondline.insert(line?);
        let tokens = current_tokens.insert(Tokens::split(line)?);
        let from = tokens.index(py, 4)?;
        let to = tokens.index(py, 5)?;
        if mapping.get(&from).is_some_and(|values| values.len() > 1)
            || mapping.get(&to).is_some_and(|values| values.len() > 1)
        {
            return Err(PyErr::from_value(
                duplicate_error.call1(("atom index mapping in bond",))?,
            ));
        }
        // A missing endpoint is the legacy empty-list subscript error.
        let new_from = *mapping
            .get(&from)
            .and_then(|values| values.first())
            .ok_or_else(|| PyIndexError::new_err("list index out of range"))?;
        let new_to = *mapping
            .get(&to)
            .and_then(|values| values.first())
            .ok_or_else(|| PyIndexError::new_err("list index out of range"))?;
        let trailing = if let Some(text) = exact_text(line)? {
            text.ends_with('\n')
        } else {
            line.get_item(-1)?.eq(PyString::new(py, "\n"))?
        };
        let trailing = if trailing { "\n" } else { "" };
        let replacement = match tokens {
            Tokens::Text(tokens) => PyString::new(
                py,
                &format!(
                    "M  V30 {} {} {new_from} {new_to}{trailing}",
                    position + 1,
                    tokens[3]
                ),
            )
            .into_any(),
            Tokens::Protocol(_) => PyString::new(py, "M  V30 {} {} {} {}{}").call_method1(
                "format",
                (position + 1, tokens.get(py, 3)?, new_from, new_to, trailing),
            )?,
        };
        bonds.push(replacement);
    }
    // The unused newline check is still observable on protocol operands.
    {
        let counts = table.getattr("counts")?;
        if let Some(text) = exact_text(&counts)? {
            if text.is_empty() {
                return Err(PyIndexError::new_err("string index out of range"));
            }
        } else {
            counts.get_item(-1)?.eq(PyString::new(py, "\n"))?;
        }
    }
    let count_tokens = Tokens::split(&table.getattr("counts")?)?;
    let prefix = format!("M  V30 COUNTS {} {} ", atoms.len(), bonds.len());
    let counts = match &count_tokens {
        Tokens::Text(tokens) => PyString::new(
            py,
            &(prefix + &tokens.get(5..).unwrap_or_default().join(" ")),
        )
        .into_any(),
        Tokens::Protocol(tokens) => {
            PyString::new(py, &prefix).add(py.get_type::<PyString>().call_method1(
                "join",
                (
                    " ",
                    tokens.get_item(py.get_type::<PySlice>().call1((
                        5,
                        py.None(),
                        py.None(),
                    ))?)?,
                ),
            )?)?
        }
    };
    let mut lines = Vec::new();
    let mut appending = true;
    let mut current_line = None;
    for line in table.getattr("lines")?.try_iter()? {
        let line = current_line.insert(line?);
        if appending {
            lines.push(line.clone());
        }
        if startswith(line, "M  V30 COUNTS")? {
            lines
                .pop()
                .ok_or_else(|| PyIndexError::new_err("pop from an empty deque"))?;
            lines.push(counts.clone());
        }
        if startswith(line, "M  V30 BEGIN ATOM")? {
            appending = false;
        }
        if startswith(line, "M  V30 END ATOM")? {
            appending = true;
            lines.extend(atoms.iter().cloned());
            lines.push(line.clone());
        }
        if startswith(line, "M  V30 BEGIN BOND")? {
            appending = false;
        }
        if startswith(line, "M  V30 END BOND")? {
            appending = true;
            lines.extend(bonds.iter().cloned());
            lines.push(line.clone());
        }
    }
    table.setattr(
        "lines",
        py.import("collections")?
            .getattr("deque")?
            .call1((PyList::new(py, lines)?,))?,
    )?;
    // Preserve observable protocol-object release order after successful replacement.
    drop(current_atomline);
    drop(current_tokens);
    drop(current_bondline);
    drop(count_tokens);
    drop(current_line);
    Ok(())
}

#[pyfunction]
fn renumber_ctable(
    py: Python<'_>,
    table: &Bound<'_, PyAny>,
    v3000: &Bound<'_, PyAny>,
    duplicate_error: &Bound<'_, PyAny>,
) -> PyResult<()> {
    if !table.getattr("format")?.is(v3000) {
        return Err(PyNotImplementedError::new_err(()));
    }
    if let Some(lines) = fast_renumber(py, table, duplicate_error)? {
        table.setattr(
            "lines",
            py.import("collections")?
                .getattr("deque")?
                .call1((PyList::new(py, lines)?,))?,
        )?;
        return Ok(());
    }
    renumber_compat(py, table, v3000, duplicate_error)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(valid_atom_indices, m)?)?;
    m.add_function(wrap_pyfunction!(valid_bond_indices, m)?)?;
    m.add_function(wrap_pyfunction!(renumber_ctable, m)?)?;
    Ok(())
}
