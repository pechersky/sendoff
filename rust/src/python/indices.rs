use super::exact_text;
use crate::core::{
    ctable,
    indices::{self, Index},
};
use pyo3::prelude::*;
use pyo3::{
    exceptions::{PyIndexError, PyNotImplementedError},
    types::{PyDict, PyInt, PyList, PySlice, PyString},
};

enum Tokens<'py> {
    Text(Vec<String>),
    Protocol(Bound<'py, PyAny>),
}

impl<'py> Tokens<'py> {
    fn split(line: &Bound<'py, PyAny>) -> PyResult<Self> {
        if let Some(text) = exact_text(line)? {
            Ok(Self::Text(
                ctable::tokens(text)
                    .into_iter()
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

    fn width(&self, position: usize) -> PyResult<usize> {
        match self {
            Self::Text(tokens) => tokens
                .get(position)
                .map(|token| token.chars().count())
                .ok_or_else(|| PyIndexError::new_err("list index out of range")),
            Self::Protocol(tokens) => tokens.get_item(position)?.len(),
        }
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

fn atom_replacement<'py>(
    py: Python<'py>,
    line: &Bound<'py, PyAny>,
    tokens: &Tokens<'py>,
    new_index: usize,
) -> PyResult<Bound<'py, PyAny>> {
    let skip = 7 + tokens.width(2)?;
    let replacement = if let Some(text) = exact_text(line)? {
        PyString::new(py, &indices::atom_line(text, new_index)).into_any()
    } else {
        PyString::new(py, &format!("M  V30 {new_index}")).add(
            line.get_item(
                py.get_type::<PySlice>()
                    .call1((skip, py.None(), py.None()))?,
            )?,
        )?
    };
    Ok(replacement)
}

fn bond_replacement<'py>(
    py: Python<'py>,
    line: &Bound<'py, PyAny>,
    tokens: &Tokens<'py>,
    position: usize,
    mapping: &indices::Mapping,
    duplicate_error: &Bound<'_, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    let from = tokens.index(py, 4)?;
    let to = tokens.index(py, 5)?;
    let (new_from, new_to) = match mapping.remap_bond(from, to) {
        Ok(values) => values,
        Err(indices::MappingFailure::Duplicate) => {
            return Err(PyErr::from_value(
                duplicate_error.call1(("atom index mapping in bond",))?,
            ));
        }
        Err(indices::MappingFailure::Missing) => {
            return Err(PyIndexError::new_err("list index out of range"));
        }
    };
    let replacement = match tokens {
        Tokens::Text(tokens) => PyString::new(
            py,
            &indices::bond_line(
                position + 1,
                &tokens[3],
                new_from,
                new_to,
                exact_text(line)?.is_some_and(|text| text.ends_with('\n')),
            ),
        )
        .into_any(),
        Tokens::Protocol(_) => PyString::new(py, "M  V30 {} {} {} {}{}").call_method1(
            "format",
            (
                position + 1,
                tokens.get(py, 3)?,
                new_from,
                new_to,
                if line.get_item(-1)?.eq(PyString::new(py, "\n"))? {
                    "\n"
                } else {
                    ""
                },
            ),
        )?,
    };
    Ok(replacement)
}

fn counts_replacement<'py>(
    py: Python<'py>,
    counts: &Bound<'py, PyAny>,
    tokens: &Tokens<'py>,
    atom_count: usize,
    bond_count: usize,
) -> PyResult<Bound<'py, PyAny>> {
    let prefix = format!("M  V30 COUNTS {atom_count} {bond_count} ");
    let replacement = match tokens {
        Tokens::Text(tokens) => {
            let suffix = exact_text(counts)?
                .map(|text| ctable::token_suffix(text, 5))
                .unwrap_or_else(|| tokens.get(5..).unwrap_or_default().join(" "));
            PyString::new(py, &indices::counts_line(atom_count, bond_count, &suffix)).into_any()
        }
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
    Ok(replacement)
}

fn line_flags(line: &Bound<'_, PyAny>) -> PyResult<indices::LineFlags> {
    Ok(indices::LineFlags {
        counts: startswith(line, "M  V30 COUNTS")?,
        begin_atom: startswith(line, "M  V30 BEGIN ATOM")?,
        end_atom: startswith(line, "M  V30 END ATOM")?,
        begin_bond: startswith(line, "M  V30 BEGIN BOND")?,
        end_bond: startswith(line, "M  V30 END BOND")?,
    })
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
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut mapping = indices::Mapping::new();
    let mut current_atomline = None;
    let mut current_tokens = None;
    let mut current_bondline = None;
    for line in table.call_method0("atomlines")?.try_iter()? {
        let line = current_atomline.insert(line?);
        let tokens = current_tokens.insert(Tokens::split(line)?);
        let index = tokens.index(py, 2)?;
        let new_index = mapping.add_atom(index);
        atoms.push(atom_replacement(py, line, tokens, new_index)?);
    }
    for (position, line) in table.call_method0("bondlines")?.try_iter()?.enumerate() {
        let line = current_bondline.insert(line?);
        let tokens = current_tokens.insert(Tokens::split(line)?);
        bonds.push(bond_replacement(
            py,
            line,
            tokens,
            position,
            &mapping,
            duplicate_error,
        )?);
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
    let counts = table.getattr("counts")?;
    let count_tokens = Tokens::split(&counts)?;
    let counts = counts_replacement(py, &counts, &count_tokens, atoms.len(), bonds.len())?;
    let mut assembler = indices::Assembler::new(atoms, bonds, counts);
    let mut current_line = None;
    for line in table.getattr("lines")?.try_iter()? {
        let line = current_line.insert(line?);
        assembler.push(line.clone(), line_flags(line)?);
    }
    let lines = assembler.finish();
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

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(valid_atom_indices, m)?)?;
    m.add_function(wrap_pyfunction!(valid_bond_indices, m)?)?;
    m.add_function(wrap_pyfunction!(renumber_ctable, m)?)?;
    Ok(())
}
