use pyo3::{
    class::{PyTraverseError, PyVisit},
    exceptions::{PyRuntimeError, PyStopIteration, PyValueError},
    prelude::*,
    types::PyIterator,
};

use crate::core::framing::{BlockState, FramingMode, FramingState};

#[pyclass(name = "FramingIter", module = "sendoff.native")]
struct FramingIter {
    framing: FramingState,
    lines: Option<Py<PyAny>>,
    iterator: Option<Py<PyAny>>,
    current: Option<Py<PyAny>>,
    done: bool,
    running: bool,
}

impl FramingIter {
    fn finish(&mut self) {
        self.lines = None;
        self.iterator = None;
        self.current = None;
        self.done = true;
    }
}

#[pymethods]
impl FramingIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(slf: &Bound<'_, Self>) -> PyResult<Option<Py<PyAny>>> {
        {
            let mut state = slf.borrow_mut();
            if state.running {
                return Err(PyValueError::new_err("generator already executing"));
            }
            if state.done {
                return Ok(None);
            }
            state.running = true;
        }

        let result = framing_step(slf);
        let mut state = slf.borrow_mut();
        state.running = false;
        if result.is_err() || matches!(&result, Ok(None)) {
            state.finish();
        }
        result
    }

    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.lines)?;
        visit.call(&self.iterator)?;
        visit.call(&self.current)
    }

    fn __clear__(&mut self) {
        self.finish();
    }
}

pub(crate) fn generator_error(py: Python<'_>, error: PyErr) -> PyErr {
    if !error.is_instance_of::<PyStopIteration>(py) {
        return error;
    }

    let runtime = pyo3::exceptions::PyRuntimeError::new_err("generator raised StopIteration");
    runtime.set_cause(py, Some(error.clone_ref(py)));
    runtime.set_context(py, Some(error));
    if let Err(error) = runtime.value(py).setattr("__suppress_context__", true) {
        return error;
    }
    runtime
}

fn has_prefix(py: Python<'_>, line: &Py<PyAny>, prefix: &str) -> PyResult<bool> {
    let line = line.bind(py);
    if let Some(text) = super::exact_text(line)? {
        return Ok(text.starts_with(prefix));
    }
    let matched = line
        .call_method1("startswith", (prefix,))
        .map_err(|error| generator_error(py, error))?;
    matched
        .is_truthy()
        .map_err(|error| generator_error(py, error))
}

fn source_iterator(slf: &Bound<'_, FramingIter>) -> PyResult<Py<PyAny>> {
    let py = slf.py();
    if let Some(iterator) = slf.borrow().iterator.as_ref() {
        return Ok(iterator.clone_ref(py));
    }

    let lines = slf
        .borrow()
        .lines
        .as_ref()
        .ok_or_else(|| PyRuntimeError::new_err("native iterator state was cleared"))?
        .clone_ref(py);
    let iterator = lines
        .bind(py)
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

fn framing_step(slf: &Bound<'_, FramingIter>) -> PyResult<Option<Py<PyAny>>> {
    let py = slf.py();
    let current = slf.borrow().current.as_ref().map(|line| line.clone_ref(py));
    if let Some(line) = current {
        let (checks_current, prefix) = {
            let state = slf.borrow();
            (state.framing.checks_current(), state.framing.prefix())
        };
        if checks_current {
            let matches = has_prefix(py, &line, prefix)?;
            let stop = slf.borrow().framing.stop_before_next(matches);
            if stop {
                slf.borrow_mut().finish();
                return Ok(None);
            }
        }
    }

    let iterator = source_iterator(slf)?;
    let Some(line) = next_item(py, &iterator)? else {
        slf.borrow_mut().finish();
        return Ok(None);
    };

    let (is_mdl, prefix) = {
        let state = slf.borrow();
        (state.framing.checks_current(), state.framing.prefix())
    };
    let should_yield = if is_mdl {
        true
    } else {
        let matches = has_prefix(py, &line, prefix)?;
        slf.borrow().framing.yield_line(matches)
    };
    slf.borrow_mut().current = Some(line.clone_ref(py));
    if !should_yield {
        slf.borrow_mut().finish();
        return Ok(None);
    }
    Ok(Some(line))
}

#[pyfunction]
fn mdl_iter(lines: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    Py::new(
        lines.py(),
        FramingIter {
            framing: FramingState::new(FramingMode::Mdl),
            lines: Some(lines.clone().unbind()),
            iterator: None,
            current: None,
            done: false,
            running: false,
        },
    )
    .map(Py::into_any)
}

#[pyfunction]
fn metadata_iter(lines: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    Py::new(
        lines.py(),
        FramingIter {
            framing: FramingState::new(FramingMode::Metadata),
            lines: Some(lines.clone().unbind()),
            iterator: None,
            current: None,
            done: false,
            running: false,
        },
    )
    .map(Py::into_any)
}

#[pyclass(name = "BlocksIter", module = "sendoff.native")]
struct BlocksIter {
    cls: Option<Py<PyAny>>,
    lines: Option<Py<PyAny>>,
    iterator: Option<Py<PyAny>>,
    blocks: BlockState<Py<PyAny>>,
    block: Option<Py<PyAny>>,
    current: Option<Py<PyAny>>,
    done: bool,
    running: bool,
}

impl BlocksIter {
    fn finish(&mut self) {
        self.cls = None;
        self.lines = None;
        self.iterator = None;
        self.blocks = BlockState::new();
        self.block = None;
        self.current = None;
        self.done = true;
    }
}

#[pymethods]
impl BlocksIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(slf: &Bound<'_, Self>) -> PyResult<Option<Py<PyAny>>> {
        {
            let mut state = slf.borrow_mut();
            if state.running {
                return Err(PyValueError::new_err("generator already executing"));
            }
            if state.done {
                return Ok(None);
            }
            state.running = true;
        }

        let result = blocks_step(slf);
        let mut state = slf.borrow_mut();
        state.running = false;
        if result.is_err() || matches!(&result, Ok(None)) {
            state.finish();
        }
        result
    }

    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.cls)?;
        visit.call(&self.lines)?;
        visit.call(&self.iterator)?;
        visit.call(&self.block)?;
        visit.call(&self.current)?;
        for line in self.blocks.lines() {
            visit.call(line)?;
        }
        Ok(())
    }

    fn __clear__(&mut self) {
        self.finish();
    }
}

fn blocks_iterator(slf: &Bound<'_, BlocksIter>) -> PyResult<Py<PyAny>> {
    let py = slf.py();
    if let Some(iterator) = slf.borrow().iterator.as_ref() {
        return Ok(iterator.clone_ref(py));
    }

    let lines = slf
        .borrow()
        .lines
        .as_ref()
        .ok_or_else(|| PyRuntimeError::new_err("native iterator state was cleared"))?
        .clone_ref(py);
    let iterator = lines
        .bind(py)
        .try_iter()
        .map_err(|error| generator_error(py, error))?
        .into_any()
        .unbind();
    slf.borrow_mut().iterator = Some(iterator.clone_ref(py));
    Ok(iterator)
}

fn blocks_step(slf: &Bound<'_, BlocksIter>) -> PyResult<Option<Py<PyAny>>> {
    let py = slf.py();
    let previous = slf.borrow_mut().block.take();
    drop(previous);
    let iterator = blocks_iterator(slf)?;
    loop {
        let Some(line) = next_item(py, &iterator)? else {
            slf.borrow_mut().finish();
            return Ok(None);
        };
        slf.borrow_mut().current = Some(line.clone_ref(py));
        let delimiter = has_prefix(py, &line, "$$$$")?;
        let Some(lines) = slf.borrow_mut().blocks.push(line, delimiter) else {
            continue;
        };

        let deque = py.import("collections")?.getattr("deque")?;
        let block_lines = deque.call0()?;
        for line in lines {
            block_lines.call_method1("append", (line,))?;
        }
        slf.borrow_mut().block = Some(block_lines.clone().unbind());
        let cls = slf
            .borrow()
            .cls
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("native iterator state was cleared"))?
            .clone_ref(py);
        return cls
            .bind(py)
            .call_method1("from_block_lines", (block_lines,))
            .map_err(|error| generator_error(py, error))
            .map(Bound::unbind)
            .map(Some);
    }
}

#[pyfunction]
fn blocks_iter(cls: &Bound<'_, PyAny>, lines: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    Py::new(
        lines.py(),
        BlocksIter {
            cls: Some(cls.clone().unbind()),
            lines: Some(lines.clone().unbind()),
            iterator: None,
            blocks: BlockState::new(),
            block: None,
            current: None,
            done: false,
            running: false,
        },
    )
    .map(Py::into_any)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FramingIter>()?;
    m.add_class::<BlocksIter>()?;
    m.add_function(wrap_pyfunction!(mdl_iter, m)?)?;
    m.add_function(wrap_pyfunction!(metadata_iter, m)?)?;
    m.add_function(wrap_pyfunction!(blocks_iter, m)?)?;
    Ok(())
}
