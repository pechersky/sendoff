use pyo3::{
    class::{PyTraverseError, PyVisit},
    exceptions::{PyRuntimeError, PyStopIteration, PyValueError},
    prelude::*,
    types::{PyIterator, PyString},
};

enum FramingMode {
    Mdl,
    Metadata,
}

#[pyclass(name = "FramingIter", module = "sendoff.native")]
struct FramingIter {
    mode: FramingMode,
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

    let runtime = PyRuntimeError::new_err("generator raised StopIteration");
    runtime.set_cause(py, Some(error.clone_ref(py)));
    runtime.set_context(py, Some(error));
    if let Err(error) = runtime.value(py).setattr("__suppress_context__", true) {
        return error;
    }
    runtime
}

fn has_prefix(py: Python<'_>, line: &Py<PyAny>, prefix: &str) -> PyResult<bool> {
    let line = line.bind(py);
    if let Some(text) = crate::exact_text(line)? {
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
    let mode_is_mdl = matches!(slf.borrow().mode, FramingMode::Mdl);

    if mode_is_mdl {
        let current = slf.borrow().current.as_ref().map(|line| line.clone_ref(py));
        if let Some(line) = current
            && has_prefix(py, &line, "M  END")?
        {
            slf.borrow_mut().finish();
            return Ok(None);
        }
    }

    let iterator = source_iterator(slf)?;
    let Some(line) = next_item(py, &iterator)? else {
        slf.borrow_mut().finish();
        return Ok(None);
    };

    slf.borrow_mut().current = Some(line.clone_ref(py));
    if !mode_is_mdl && has_prefix(py, &line, "$$$$")? {
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
            mode: FramingMode::Mdl,
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
            mode: FramingMode::Metadata,
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
    block: Option<Py<PyAny>>,
    current: Option<Py<PyAny>>,
    reset_after_yield: bool,
    done: bool,
    running: bool,
}

impl BlocksIter {
    fn finish(&mut self) {
        self.cls = None;
        self.lines = None;
        self.iterator = None;
        self.block = None;
        self.current = None;
        self.reset_after_yield = false;
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
        visit.call(&self.current)
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

fn empty_deque(py: Python<'_>) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("collections")?
        .getattr("deque")?
        .call0()?
        .unbind())
}

fn blocks_step(slf: &Bound<'_, BlocksIter>) -> PyResult<Option<Py<PyAny>>> {
    let py = slf.py();
    if slf.borrow().reset_after_yield {
        let block = empty_deque(py)?;
        let mut state = slf.borrow_mut();
        state.block = Some(block);
        state.reset_after_yield = false;
    }

    if slf.borrow().block.is_none() {
        let block = empty_deque(py)?;
        slf.borrow_mut().block = Some(block);
    }

    let iterator = blocks_iterator(slf)?;
    loop {
        let Some(line) = next_item(py, &iterator)? else {
            slf.borrow_mut().finish();
            return Ok(None);
        };
        slf.borrow_mut().current = Some(line.clone_ref(py));

        let block = slf
            .borrow()
            .block
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("native iterator state was cleared"))?
            .clone_ref(py);
        block.bind(py).call_method1("append", (line.bind(py),))?;
        if !has_prefix(py, &line, "$$$$")? {
            continue;
        }

        let cls = slf
            .borrow()
            .cls
            .as_ref()
            .ok_or_else(|| PyRuntimeError::new_err("native iterator state was cleared"))?
            .clone_ref(py);
        let parsed = cls
            .bind(py)
            .call_method1("from_block_lines", (block.bind(py),))
            .map_err(|error| generator_error(py, error))?
            .unbind();
        slf.borrow_mut().reset_after_yield = true;
        return Ok(Some(parsed));
    }
}

#[pyfunction]
fn from_block_lines(
    cls: &Bound<'_, PyAny>,
    block_type: &Bound<'_, PyAny>,
    lines: &Bound<'_, PyAny>,
) -> PyResult<Py<PyAny>> {
    let py = lines.py();
    let iterator = lines.try_iter()?;
    let builtins = py.import("builtins")?;
    let title = builtins
        .getattr("next")?
        .call1((iterator.as_any(),))?
        .into_any();
    let title = if let Some(text) = crate::exact_text(&title)? {
        let stripped = text.trim_matches(crate::whitespace);
        if stripped == text {
            title.clone()
        } else {
            PyString::new(py, stripped).into_any()
        }
    } else {
        title.call_method0("strip")?
    };
    let deque = py.import("collections")?.getattr("deque")?;
    let mdl = deque.call1((cls.call_method1("parse_mdl", (iterator.as_any(),))?,))?;
    let metadata = deque.call1((cls.call_method1("parse_metadata", (iterator.as_any(),))?,))?;
    block_type.call1((title, mdl, metadata)).map(Bound::unbind)
}

#[pyfunction]
fn blocks_iter(cls: &Bound<'_, PyAny>, lines: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    Py::new(
        lines.py(),
        BlocksIter {
            cls: Some(cls.clone().unbind()),
            lines: Some(lines.clone().unbind()),
            iterator: None,
            block: None,
            current: None,
            reset_after_yield: false,
            done: false,
            running: false,
        },
    )
    .map(Py::into_any)
}

#[pyfunction]
fn read_sdf_lines(sdfpath: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
    let builtins = sdfpath.py().import("builtins")?;
    let file = builtins.getattr("open")?.call1((sdfpath,))?;
    file.call_method0("readlines").map(Bound::unbind)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FramingIter>()?;
    m.add_class::<BlocksIter>()?;
    m.add_function(wrap_pyfunction!(mdl_iter, m)?)?;
    m.add_function(wrap_pyfunction!(metadata_iter, m)?)?;
    m.add_function(wrap_pyfunction!(from_block_lines, m)?)?;
    m.add_function(wrap_pyfunction!(blocks_iter, m)?)?;
    m.add_function(wrap_pyfunction!(read_sdf_lines, m)?)?;
    Ok(())
}
