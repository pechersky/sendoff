//! Pure Rust operations for SDF blocks.

use std::convert::Infallible;

/// The parsed data stored in one molecule block.
pub struct SDBlock {
    /// The stripped molecule title.
    pub title: String,
    /// The MDL lines, including the `M  END` marker.
    pub mdl: Vec<String>,
    /// The SD metadata lines, excluding the `$$$$` delimiter.
    pub metadata: Vec<String>,
}

impl SDBlock {
    /// Parse SDBlocks from a sequence of lines.
    ///
    /// Blocks are separated by lines starting with `$$$$`. A nonempty final
    /// block is parsed when the input ends without a delimiter.
    pub fn from_lines<I, E>(lines: &mut I) -> Result<Vec<Self>, E>
    where
        I: Iterator<Item = Result<String, E>>,
    {
        split_blocks(lines).map(|blocks| blocks.into_iter().map(Self::from_block_lines).collect())
    }

    pub(crate) fn from_block_lines(block_lines: Vec<String>) -> Self {
        let (title, mdl, metadata) = parse_block_lines(block_lines);
        Self {
            title,
            mdl,
            metadata,
        }
    }

    /// Generate SD metadata records.
    ///
    /// Each value includes embedded newlines when it spans multiple lines.
    pub fn records(&self) -> Result<Vec<(String, String)>, crate::core::sddata::RecordParseError> {
        crate::core::sddata::records(&self.metadata)
    }

    /// Append a field and value to the SD data record.
    pub fn append_record(&mut self, record_name: &str, value: &str) {
        let (name, value, separator) = crate::core::sddata::append_record(record_name, value);
        self.metadata.extend([name, value, separator]);
    }

    /// Return the title and MDL lines expected by `CTable`.
    ///
    /// The returned lines are parsed as a CTable without additional validation.
    pub fn ctable_lines(&self) -> Vec<String> {
        let mut lines = Vec::with_capacity(self.mdl.len() + 1);
        lines.push(self.title.clone());
        lines.extend(self.mdl.iter().cloned());
        lines
    }

    /// Serialize this block, preserving existing line endings when requested.
    pub fn serialize(&self, with_newlines: bool) -> String {
        serialize(&self.title, &self.mdl, &self.metadata, with_newlines)
    }
}

/// The parsed data stored in one molecule block, without text decoding.
pub struct BytesSDBlock {
    /// The stripped molecule title.
    pub title: Vec<u8>,
    /// The MDL lines, including the `M  END` marker.
    pub mdl: Vec<Vec<u8>>,
    /// The SD metadata lines, excluding the `$$$$` delimiter.
    pub metadata: Vec<Vec<u8>>,
}

impl BytesSDBlock {
    /// Parse SDBlocks from a sequence of byte lines.
    ///
    /// Blocks are separated by lines starting with `$$$$`. A nonempty final
    /// block is parsed when the input ends without a delimiter.
    pub fn from_lines<I, E>(lines: &mut I) -> Result<Vec<Self>, E>
    where
        I: Iterator<Item = Result<Vec<u8>, E>>,
    {
        split_blocks(lines).map(|blocks| blocks.into_iter().map(Self::from_block_lines).collect())
    }

    pub(crate) fn from_block_lines(block_lines: Vec<Vec<u8>>) -> Self {
        let (title, mdl, metadata) = parse_block_lines(block_lines);
        Self {
            title,
            mdl,
            metadata,
        }
    }

    /// Generate SD metadata records without decoding their bytes.
    pub fn records(
        &self,
    ) -> Result<Vec<(Vec<u8>, Vec<u8>)>, crate::core::sddata::RecordParseError> {
        crate::core::sddata::records(&self.metadata)
    }

    /// Append a byte-valued SD data record.
    pub fn append_record(&mut self, record_name: &[u8], value: &[u8]) {
        let mut name = Vec::with_capacity(record_name.len() + 5);
        name.extend_from_slice(b"> <");
        name.extend_from_slice(record_name);
        name.extend_from_slice(b">\n");

        let mut value_line = Vec::with_capacity(value.len() + 1);
        value_line.extend_from_slice(value);
        value_line.push(b'\n');
        self.metadata.extend([name, value_line, b"\n".to_vec()]);
    }

    /// Return the title and MDL lines expected by `CTable`.
    pub fn ctable_lines(&self) -> Vec<Vec<u8>> {
        let mut lines = Vec::with_capacity(self.mdl.len() + 1);
        lines.push(self.title.clone());
        lines.extend(self.mdl.iter().cloned());
        lines
    }

    /// Serialize this block, preserving existing line endings when requested.
    pub fn serialize(&self, with_newlines: bool) -> Vec<u8> {
        let mut output = Vec::new();
        output.extend_from_slice(&self.title);
        output.push(b'\n');
        for line in self.mdl.iter().chain(&self.metadata) {
            output.extend_from_slice(line);
            if with_newlines && !line.ends_with(b"\n") {
                output.push(b'\n');
            }
        }
        output.extend_from_slice(b"$$$$\n");
        output
    }
}

trait SdfLine: Sized {
    fn starts_with_ascii(&self, prefix: &[u8]) -> bool;
    fn trim_title(self) -> Self;
}

impl SdfLine for String {
    fn starts_with_ascii(&self, prefix: &[u8]) -> bool {
        self.as_bytes().starts_with(prefix)
    }

    fn trim_title(self) -> Self {
        self.trim().to_owned()
    }
}

impl SdfLine for Vec<u8> {
    fn starts_with_ascii(&self, prefix: &[u8]) -> bool {
        self.starts_with(prefix)
    }

    fn trim_title(self) -> Self {
        self.trim_ascii().to_vec()
    }
}

fn split_blocks<T, E, I>(lines: &mut I) -> Result<Vec<Vec<T>>, E>
where
    T: SdfLine,
    I: Iterator<Item = Result<T, E>>,
{
    let mut blocks = Vec::new();
    let mut block_lines = Vec::new();
    for line in lines {
        let line = line?;
        let is_end = line.starts_with_ascii(b"$$$$");
        block_lines.push(line);
        if is_end {
            blocks.push(std::mem::take(&mut block_lines));
        }
    }
    if !block_lines.is_empty() {
        blocks.push(block_lines);
    }
    Ok(blocks)
}

fn parse_block_lines<T: SdfLine>(mut block_lines: Vec<T>) -> (T, Vec<T>, Vec<T>) {
    let title = block_lines.remove(0).trim_title();
    let mut lines = block_lines.into_iter().map(Ok::<T, Infallible>);
    let mdl = collect_until(&mut lines, b"M  END", true).expect("infallible lines");
    let metadata = collect_until(&mut lines, b"$$$$", false).expect("infallible lines");
    (title, mdl, metadata)
}

fn collect_until<I, T, E>(lines: &mut I, marker: &[u8], include_marker: bool) -> Result<Vec<T>, E>
where
    T: SdfLine,
    I: Iterator<Item = Result<T, E>>,
{
    let mut collected = Vec::new();
    for line in lines {
        let line = line?;
        let found_marker = line.starts_with_ascii(marker);
        if found_marker {
            if include_marker {
                collected.push(line);
            }
            break;
        }
        collected.push(line);
    }
    Ok(collected)
}

/// Collect lines through the first `M  END` line, including that line.
pub fn parse_mdl<I, E>(lines: &mut I) -> Result<Vec<String>, E>
where
    I: Iterator<Item = Result<String, E>>,
{
    collect_until(lines, b"M  END", true)
}

/// Collect metadata lines before the first `$$$$` delimiter.
pub fn parse_metadata<I, E>(lines: &mut I) -> Result<Vec<String>, E>
where
    I: Iterator<Item = Result<String, E>>,
{
    collect_until(lines, b"$$$$", false)
}

/// Serialize an SDF block, preserving existing line endings when requested.
pub fn serialize(
    title: &str,
    mdl_lines: &[String],
    metadata_lines: &[String],
    with_newlines: bool,
) -> String {
    let mut output = String::new();
    output.push_str(title);
    output.push('\n');
    for line in mdl_lines.iter().chain(metadata_lines) {
        output.push_str(line);
        if with_newlines && !line.ends_with('\n') {
            output.push('\n');
        }
    }
    output.push_str("$$$$\n");
    output
}
