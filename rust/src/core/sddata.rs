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

pub struct Value<T> {
    pub payload: T,
    pub text: Option<String>,
}

pub struct Line<T> {
    pub key: bool,
    pub value: Value<T>,
}

pub enum Name<T> {
    Source(T),
    Text(String),
}

pub enum Input<T> {
    Item(Line<T>),
    End,
    KeyStop,
}

pub enum Head<T> {
    No,
    Yes(Name<T>),
}

pub enum JoinedValue<T> {
    Source(T),
    Text(String),
    Fallback(Vec<T>),
}

pub struct Record<T> {
    pub name: Name<T>,
    pub value: JoinedValue<T>,
}

pub struct RecordsState<T> {
    pending: Option<Line<T>>,
    drain_key: Option<bool>,
    exhausted: bool,
}

impl<T> RecordsState<T> {
    pub const fn new() -> Self {
        Self {
            pending: None,
            drain_key: None,
            exhausted: false,
        }
    }

    pub fn pending(&self) -> Option<&Line<T>> {
        self.pending.as_ref()
    }

    pub fn next_record<E, F, H>(&mut self, mut next: F, mut head: H) -> Result<Option<Record<T>>, E>
    where
        F: FnMut() -> Result<Input<T>, E>,
        H: FnMut(&Line<T>) -> Result<Head<T>, E>,
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

            let name = match head(&line)? {
                Head::No => {
                    if !self.drain_group(line.key, true, &mut next)? {
                        self.exhausted = true;
                        return Ok(None);
                    }
                    continue;
                }
                Head::Yes(name) => name,
            };

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
                value: join_values(values),
            }));
        }
    }

    fn drain_group<E, F>(&mut self, key: bool, key_stop_ends: bool, next: &mut F) -> Result<bool, E>
    where
        F: FnMut() -> Result<Input<T>, E>,
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

impl<T> Default for RecordsState<T> {
    fn default() -> Self {
        Self::new()
    }
}

fn join_values<T>(values: Vec<Value<T>>) -> JoinedValue<T> {
    if values.len() == 1 && values[0].text.is_some() {
        let mut values = values.into_iter();
        let value = match values.next() {
            Some(value) => value,
            None => unreachable!("one value was checked above"),
        };
        return JoinedValue::Source(value.payload);
    }

    let mut joined = String::new();
    let mut first = true;
    let mut values = values.into_iter();
    while let Some(value) = values.next() {
        let Some(text) = value.text else {
            let mut fallback = vec![value.payload];
            fallback.extend(values.map(|value| value.payload));
            return JoinedValue::Fallback(fallback);
        };
        if !first {
            joined.push('\n');
        }
        first = false;
        joined.push_str(&text);
    }
    JoinedValue::Text(joined)
}
