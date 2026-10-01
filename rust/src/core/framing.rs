#[derive(Clone, Copy)]
pub enum FramingMode {
    Mdl,
    Metadata,
}

impl FramingMode {
    pub const fn prefix(self) -> &'static str {
        match self {
            Self::Mdl => "M  END",
            Self::Metadata => "$$$$",
        }
    }
}

pub struct FramingState {
    mode: FramingMode,
}

impl FramingState {
    pub const fn new(mode: FramingMode) -> Self {
        Self { mode }
    }

    pub const fn prefix(&self) -> &'static str {
        self.mode.prefix()
    }

    pub const fn checks_current(&self) -> bool {
        matches!(self.mode, FramingMode::Mdl)
    }

    pub const fn stop_before_next(&self, current_matches: bool) -> bool {
        matches!(self.mode, FramingMode::Mdl) && current_matches
    }

    pub const fn yield_line(&self, line_matches: bool) -> bool {
        !matches!(self.mode, FramingMode::Metadata) || !line_matches
    }
}

pub struct BlockState<T> {
    lines: Vec<T>,
}

impl<T> BlockState<T> {
    pub fn new() -> Self {
        Self { lines: Vec::new() }
    }

    pub fn lines(&self) -> &[T] {
        &self.lines
    }

    pub fn push(&mut self, line: T, delimiter: bool) -> Option<Vec<T>> {
        self.lines.push(line);
        delimiter.then(|| std::mem::take(&mut self.lines))
    }
}

impl<T> Default for BlockState<T> {
    fn default() -> Self {
        Self::new()
    }
}
