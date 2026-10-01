use pyo3::{
    class::{PyTraverseError, PyVisit},
    exceptions::{PyIndexError, PyRuntimeError, PyStopIteration, PyValueError},
    ffi,
    prelude::*,
    types::{PyDict, PyIterator, PyList, PyString, PyTuple},
};

use crate::framing::generator_error;

enum RecordInput {
    Item(Py<PyAny>, bool),
    End,
    KeyStop,
}

#[pyclass(name = "RecordsIter", module = "sendoff.native")]
struct RecordsIter {
    block: Option<Py<PyAny>>,
    iterator: Option<Py<PyAny>>,
    raw_line: Option<Py<PyAny>>,
    line: Option<Py<PyAny>>,
    pending: Option<(Py<PyAny>, bool)>,
    drain_key: Option<bool>,
    exhausted: bool,
    done: bool,
    running: bool,
}

impl RecordsIter {
    fn finish(&mut self) {
        self.block = None;
        self.iterator = None;
        self.raw_line = None;
        self.line = None;
        self.pending = None;
        self.drain_key = None;
        self.exhausted = true;
        self.done = true;
    }
}

#[pymethods]
impl RecordsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(slf: &Bound<'_, Self>) -> PyResult<Option<Py<PyAny>>> {
        {
            let mut state = slf
                .try_borrow_mut()
                .map_err(|_| PyValueError::new_err("generator already executing"))?;
            if state.running {
                return Err(PyValueError::new_err("generator already executing"));
            }
            if state.done {
                return Ok(None);
            }
            state.running = true;
        }

        let result = records_step(slf);
        let mut state = slf.borrow_mut();
        state.running = false;
        if result.is_err() || matches!(&result, Ok(None)) {
            state.finish();
        }
        result
    }

    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.block)?;
        visit.call(&self.iterator)?;
        visit.call(&self.raw_line)?;
        visit.call(&self.line)?;
        if let Some((line, _)) = &self.pending {
            visit.call(line)?;
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        self.finish();
    }
}

fn source_iterator(slf: &Bound<'_, RecordsIter>) -> PyResult<Py<PyAny>> {
    let py = slf.py();
    if let Some(iterator) = slf.borrow().iterator.as_ref() {
        return Ok(iterator.clone_ref(py));
    }

    let block = slf
        .borrow()
        .block
        .as_ref()
        .ok_or_else(|| PyRuntimeError::new_err("native iterator state was cleared"))?
        .clone_ref(py);
    let metadata = block
        .bind(py)
        .getattr("metadata")
        .map_err(|error| generator_error(py, error))?;
    let iterator = metadata
        .try_iter()
        .map_err(|error| generator_error(py, error))?
        .into_any()
        .unbind();
    slf.borrow_mut().iterator = Some(iterator.clone_ref(py));
    Ok(iterator)
}

fn next_item(py: Python<'_>, iterator: &Py<PyAny>) -> PyResult<Option<Py<PyAny>>> {
    let mut iterator = iterator.bind(py).cast::<PyIterator>()?.clone();
    match iterator.next() {
        None => Ok(None),
        Some(Ok(item)) => Ok(Some(item.unbind())),
        Some(Err(error)) => Err(error),
    }
}

fn next_record_input(slf: &Bound<'_, RecordsIter>) -> PyResult<RecordInput> {
    let py = slf.py();
    let iterator = source_iterator(slf)?;
    let Some(raw_line) = next_item(py, &iterator)? else {
        return Ok(RecordInput::End);
    };
    slf.borrow_mut().raw_line = Some(raw_line.clone_ref(py));
    let raw_line_bound = raw_line.bind(py);
    let line = if let Some(text) = crate::exact_text(raw_line_bound)? {
        let stripped = text.trim_matches(crate::whitespace);
        if stripped == text {
            raw_line.clone_ref(py)
        } else {
            PyString::new(py, stripped).into_any().unbind()
        }
    } else {
        raw_line_bound
            .call_method0("strip")
            .map_err(|error| generator_error(py, error))?
            .unbind()
    };
    slf.borrow_mut().line = Some(line.clone_ref(py));
    let key = if let Some(text) = crate::exact_text(line.bind(py))? {
        Ok(!text.is_empty())
    } else {
        line.bind(py).is_truthy()
    };
    match key {
        Ok(key) => Ok(RecordInput::Item(line, key)),
        Err(error) if error.is_instance_of::<PyStopIteration>(py) => Ok(RecordInput::KeyStop),
        Err(error) => Err(error),
    }
}

fn record_name(py: Python<'_>, line: &Py<PyAny>) -> PyResult<Py<PyAny>> {
    let line = line.bind(py);
    if let Some(text) = crate::exact_text(line)? {
        let remainder = text
            .split_once("> ")
            .map(|(_, remainder)| remainder)
            .ok_or_else(|| PyIndexError::new_err("list index out of range"))?;
        let remainder = remainder.trim_matches(crate::whitespace);
        let first = remainder
            .rsplit_once('>')
            .map_or(remainder, |(first, _)| first);
        let name = first
            .split_once('<')
            .map(|(_, name)| name)
            .ok_or_else(|| PyIndexError::new_err("list index out of range"))?;
        return Ok(PyString::new(py, name).into_any().unbind());
    }
    let split = line
        .call_method1("split", ("> ", 1))
        .map_err(|error| generator_error(py, error))?;
    let part = split
        .get_item(1)
        .map_err(|error| generator_error(py, error))?;
    let stripped = part
        .call_method0("strip")
        .map_err(|error| generator_error(py, error))?;
    let rsplit = stripped
        .call_method1("rsplit", (">", 1))
        .map_err(|error| generator_error(py, error))?;
    let first = rsplit
        .get_item(0)
        .map_err(|error| generator_error(py, error))?;
    let split = first
        .call_method1("split", ("<", 1))
        .map_err(|error| generator_error(py, error))?;
    split
        .get_item(1)
        .map(Bound::unbind)
        .map_err(|error| generator_error(py, error))
}

fn is_header(py: Python<'_>, line: &Py<PyAny>) -> PyResult<bool> {
    let line = line.bind(py);
    if let Some(text) = crate::exact_text(line)? {
        return Ok(text.starts_with("> "));
    }
    let result = line
        .call_method1("startswith", ("> ",))
        .map_err(|error| generator_error(py, error))?;
    result
        .is_truthy()
        .map_err(|error| generator_error(py, error))
}

fn drain_group(slf: &Bound<'_, RecordsIter>, key: bool, key_stop_ends: bool) -> PyResult<bool> {
    loop {
        match next_record_input(slf)? {
            RecordInput::Item(_, next_key) if next_key == key => {}
            RecordInput::Item(line, next_key) => {
                slf.borrow_mut().pending = Some((line, next_key));
                return Ok(true);
            }
            RecordInput::End => return Ok(false),
            RecordInput::KeyStop if key_stop_ends => return Ok(false),
            RecordInput::KeyStop => {}
        }
    }
}

fn join_values(py: Python<'_>, values: &[Py<PyAny>]) -> PyResult<Py<PyAny>> {
    if values.len() == 1 && crate::exact_text(values[0].bind(py))?.is_some() {
        return Ok(values[0].clone_ref(py));
    }
    let mut joined = String::new();
    for (index, value) in values.iter().enumerate() {
        let Some(text) = crate::exact_text(value.bind(py))? else {
            let python_values = PyList::empty(py);
            for value in values {
                python_values.append(value.bind(py))?;
            }
            return PyString::new(py, "\n")
                .call_method1("join", (python_values,))
                .map(Bound::unbind);
        };
        if index > 0 {
            joined.push('\n');
        }
        joined.push_str(text);
    }
    Ok(PyString::new(py, &joined).into_any().unbind())
}

fn records_step(slf: &Bound<'_, RecordsIter>) -> PyResult<Option<Py<PyAny>>> {
    let py = slf.py();
    if slf.borrow().exhausted {
        return Ok(None);
    }
    let drain_key = slf.borrow_mut().drain_key.take();
    if let Some(key) = drain_key
        && !drain_group(slf, key, false)?
    {
        return Ok(None);
    }
    loop {
        let next = slf.borrow_mut().pending.take();
        let (line, key) = match next {
            Some(item) => item,
            None => match next_record_input(slf)? {
                RecordInput::Item(line, key) => (line, key),
                RecordInput::End | RecordInput::KeyStop => return Ok(None),
            },
        };

        if !is_header(py, &line)? {
            if !drain_group(slf, key, true)? {
                return Ok(None);
            }
            continue;
        }

        let name = record_name(py, &line)?;
        let mut values = Vec::new();
        loop {
            match next_record_input(slf)? {
                RecordInput::Item(value, next_key) if next_key == key => values.push(value),
                RecordInput::Item(line, next_key) => {
                    slf.borrow_mut().pending = Some((line, next_key));
                    break;
                }
                RecordInput::End => {
                    slf.borrow_mut().exhausted = true;
                    break;
                }
                RecordInput::KeyStop => {
                    slf.borrow_mut().drain_key = Some(key);
                    break;
                }
            }
        }
        let value = join_values(py, &values)?;
        return PyTuple::new(py, [name, value])
            .map(Bound::into_any)
            .map(Bound::unbind)
            .map(Some);
    }
}

#[pyfunction]
fn records_iter(block: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    Py::new(
        block.py(),
        RecordsIter {
            block: Some(block.clone().unbind()),
            iterator: None,
            raw_line: None,
            line: None,
            pending: None,
            drain_key: None,
            exhausted: false,
            done: false,
            running: false,
        },
    )
    .map(Py::into_any)
}

fn write_lines(
    py: Python<'_>,
    lines: &Bound<'_, PyAny>,
    outh: &Bound<'_, PyAny>,
    with_newlines: &Bound<'_, PyAny>,
) -> PyResult<bool> {
    let iterator = match lines.try_iter() {
        Ok(iterator) => iterator,
        Err(error) if error.is_instance_of::<PyStopIteration>(py) => return Ok(false),
        Err(error) => return Err(error),
    };

    for item in iterator {
        let line = item?;
        outh.call_method1("write", (&line,))?;
        if with_newlines.is_truthy()? {
            let ends_with_newline = match crate::exact_text(&line)? {
                Some(text) => text.ends_with('\n'),
                None => line.call_method1("endswith", ("\n",))?.is_truthy()?,
            };
            if !ends_with_newline {
                outh.call_method1("write", ("\n",))?;
            }
        }
    }
    Ok(true)
}

#[pyfunction]
fn write(
    block: &Bound<'_, PyAny>,
    outh: &Bound<'_, PyAny>,
    with_newlines: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let py = block.py();
    let print = py.import("builtins")?.getattr("print")?;
    let title = block.getattr("title")?;
    let kwargs = PyDict::new(py);
    kwargs.set_item("file", outh)?;
    print.call((title,), Some(&kwargs))?;

    let mdl = block.getattr("mdl")?;
    let metadata = block.getattr("metadata")?;
    if write_lines(py, &mdl, outh, with_newlines)? {
        write_lines(py, &metadata, outh, with_newlines)?;
    }
    outh.call_method1("write", ("$$$$\n",))?;
    Ok(())
}

fn formatted_header<'py>(
    py: Python<'py>,
    record_name: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    if let Some(record_name) = crate::exact_text(record_name)? {
        let mut header = String::with_capacity(record_name.len() + 6);
        header.push_str("> <");
        header.push_str(record_name);
        header.push_str(">\n");
        return Ok(PyString::new(py, &header).into_any());
    }
    let spec = PyString::new(py, "");
    // Both operands are live Python objects under the GIL.
    let formatted = unsafe {
        Bound::from_owned_ptr_or_err(
            py,
            ffi::PyObject_Format(record_name.as_ptr(), spec.as_ptr()),
        )?
    };
    let parts = PyList::empty(py);
    parts.append(PyString::new(py, "> <"))?;
    parts.append(formatted)?;
    parts.append(PyString::new(py, ">\n"))?;
    PyString::new(py, "").call_method1("join", (parts,))
}

fn python_add<'py>(
    py: Python<'py>,
    left: &Bound<'py, PyAny>,
    right: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyAny>> {
    // Both operands are live Python objects under the GIL.
    let ptr = unsafe { ffi::PyNumber_Add(left.as_ptr(), right.as_ptr()) };
    unsafe { Bound::from_owned_ptr_or_err(py, ptr) }
}

#[pyfunction]
fn append_record(
    block: &Bound<'_, PyAny>,
    record_name: &Bound<'_, PyAny>,
    value: &Bound<'_, PyAny>,
) -> PyResult<()> {
    let py = block.py();

    let metadata = block.getattr("metadata")?;
    let append = metadata.getattr("append")?;
    let header = formatted_header(py, record_name)?;
    append.call1((header,))?;

    let metadata = block.getattr("metadata")?;
    let append = metadata.getattr("append")?;
    let value = if let Some(value) = crate::exact_text(value)? {
        let mut line = String::with_capacity(value.len() + 1);
        line.push_str(value);
        line.push('\n');
        PyString::new(py, &line).into_any()
    } else {
        python_add(py, value, PyString::new(py, "\n").as_any())?
    };
    append.call1((value,))?;

    let metadata = block.getattr("metadata")?;
    metadata
        .getattr("append")?
        .call1((PyString::new(py, "\n"),))?;
    Ok(())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<RecordsIter>()?;
    m.add_function(wrap_pyfunction!(records_iter, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    m.add_function(wrap_pyfunction!(append_record, m)?)?;
    Ok(())
}
