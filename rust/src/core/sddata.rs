use super::whitespace;

pub fn record_name(text: &str) -> Option<&str> {
    let (_, remainder) = text.split_once("> ")?;
    let remainder = remainder.trim_matches(whitespace);
    let first = remainder
        .rsplit_once('>')
        .map_or(remainder, |(first, _)| first);
    first.split_once('<').map(|(_, name)| name)
}

pub fn header(name: &str) -> String {
    format!("> <{name}>\n")
}

pub fn value_line(value: &str) -> String {
    format!("{value}\n")
}

pub fn trimmed(text: &str) -> &str {
    text.trim_matches(whitespace)
}

pub fn is_header(text: &str) -> bool {
    text.starts_with("> ")
}

pub fn needs_newline(text: &str) -> bool {
    !text.ends_with('\n')
}

pub struct Line {
    pub key: bool,
    pub header: bool,
    pub name: Option<Name>,
    pub value: Value,
}

pub struct Value {
    pub source: usize,
    pub text: Option<String>,
}

pub enum Name {
    Source(usize),
    Text(String),
}

pub enum Input {
    Item(Line),
    End,
    KeyStop,
}

pub enum JoinedValue {
    Source(usize),
    Text(String),
    Fallback(Vec<usize>),
}

pub struct Record {
    pub name: Name,
    pub value: JoinedValue,
}

pub struct RecordsState {
    pending: Option<Line>,
    drain_key: Option<bool>,
    exhausted: bool,
}

impl RecordsState {
    pub const fn new() -> Self {
        Self {
            pending: None,
            drain_key: None,
            exhausted: false,
        }
    }

    pub fn pending_sources(&self) -> (Option<usize>, Option<usize>) {
        self.pending.as_ref().map_or((None, None), |line| {
            let name = match &line.name {
                Some(Name::Source(source)) => Some(*source),
                _ => None,
            };
            (Some(line.value.source), name)
        })
    }

    pub fn next_record<E, F>(&mut self, mut next: F) -> Result<Option<Record>, E>
    where
        F: FnMut() -> Result<Input, E>,
    {
        if self.exhausted {
            return Ok(None);
        }

        if let Some(key) = self.drain_key.take()
            && !self.drain_group(key, false, &mut next)?
        {
            self.exhausted = true;
            return Ok(None);
        }

        loop {
            let input = match self.pending.take() {
                Some(line) => Input::Item(line),
                None => next()?,
            };
            let line = match input {
                Input::Item(line) => line,
                Input::End | Input::KeyStop => {
                    self.exhausted = true;
                    return Ok(None);
                }
            };

            if !line.header {
                if !self.drain_group(line.key, true, &mut next)? {
                    self.exhausted = true;
                    return Ok(None);
                }
                continue;
            }

            let name = line
                .name
                .expect("header names are prepared at the boundary");
            let key = line.key;
            let mut values = Vec::new();
            loop {
                match next()? {
                    Input::Item(next_line) if next_line.key == key => {
                        values.push(next_line.value);
                    }
                    Input::Item(next_line) => {
                        self.pending = Some(next_line);
                        break;
                    }
                    Input::End => {
                        self.exhausted = true;
                        break;
                    }
                    Input::KeyStop => {
                        self.drain_key = Some(key);
                        break;
                    }
                }
            }
            return Ok(Some(Record {
                name,
                value: join_values(&values),
            }));
        }
    }

    fn drain_group<E, F>(&mut self, key: bool, key_stop_ends: bool, next: &mut F) -> Result<bool, E>
    where
        F: FnMut() -> Result<Input, E>,
    {
        loop {
            match next()? {
                Input::Item(line) if line.key == key => {}
                Input::Item(line) => {
                    self.pending = Some(line);
                    return Ok(true);
                }
                Input::End => return Ok(false),
                Input::KeyStop if key_stop_ends => return Ok(false),
                Input::KeyStop => {}
            }
        }
    }
}

impl Default for RecordsState {
    fn default() -> Self {
        Self::new()
    }
}

fn join_values(values: &[Value]) -> JoinedValue {
    if values.len() == 1 && values[0].text.is_some() {
        return JoinedValue::Source(values[0].source);
    }

    let mut joined = String::new();
    for (index, value) in values.iter().enumerate() {
        let Some(text) = value.text.as_deref() else {
            return JoinedValue::Fallback(values.iter().map(|value| value.source).collect());
        };
        if index > 0 {
            joined.push('\n');
        }
        joined.push_str(text);
    }
    JoinedValue::Text(joined)
}
