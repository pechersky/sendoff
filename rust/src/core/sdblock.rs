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
        let mut blocks = Vec::new();
        let mut block_lines = Vec::new();
        for line in lines {
            let line = line?;
            let is_end = line.starts_with("$$$$");
            block_lines.push(line);
            if is_end {
                blocks.push(Self::from_block_lines(std::mem::take(&mut block_lines)));
            }
        }
        if !block_lines.is_empty() {
            blocks.push(Self::from_block_lines(block_lines));
        }
        Ok(blocks)
    }

    pub(crate) fn from_block_lines(mut block_lines: Vec<String>) -> Self {
        let title = block_lines.remove(0).trim().to_owned();
        let mut lines = block_lines.into_iter().map(Ok::<String, Infallible>);
        let mdl = parse_mdl(&mut lines).expect("infallible lines");
        let metadata = parse_metadata(&mut lines).expect("infallible lines");
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
        let mut blocks = Vec::new();
        let mut block_lines = Vec::new();
        for line in lines {
            let line = line?;
            let is_end = line.starts_with(b"$$$$");
            block_lines.push(line);
            if is_end {
                blocks.push(Self::from_block_lines(std::mem::take(&mut block_lines)));
            }
        }
        if !block_lines.is_empty() {
            blocks.push(Self::from_block_lines(block_lines));
        }
        Ok(blocks)
    }

    pub(crate) fn from_block_lines(mut block_lines: Vec<Vec<u8>>) -> Self {
        let title = block_lines.remove(0).trim_ascii().to_vec();
        let mut lines = block_lines.into_iter();
        let mdl = parse_mdl_bytes(&mut lines);
        let metadata = parse_metadata_bytes(&mut lines);
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
        crate::core::sddata::records_bytes(&self.metadata)
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

/// Collect lines through the first `M  END` line, including that line.
pub fn parse_mdl<I, E>(lines: &mut I) -> Result<Vec<String>, E>
where
    I: Iterator<Item = Result<String, E>>,
{
    let mut mdl_lines = Vec::new();
    for line in lines {
        let line = line?;
        let is_end = line.starts_with("M  END");
        mdl_lines.push(line);
        if is_end {
            break;
        }
    }
    Ok(mdl_lines)
}

/// Collect metadata lines before the first `$$$$` delimiter.
pub fn parse_metadata<I, E>(lines: &mut I) -> Result<Vec<String>, E>
where
    I: Iterator<Item = Result<String, E>>,
{
    let mut metadata_lines = Vec::new();
    for line in lines {
        let line = line?;
        if line.starts_with("$$$$") {
            break;
        }
        metadata_lines.push(line);
    }
    Ok(metadata_lines)
}

fn parse_mdl_bytes<I>(lines: &mut I) -> Vec<Vec<u8>>
where
    I: Iterator<Item = Vec<u8>>,
{
    let mut mdl_lines = Vec::new();
    for line in lines {
        let is_end = line.starts_with(b"M  END");
        mdl_lines.push(line);
        if is_end {
            break;
        }
    }
    mdl_lines
}

fn parse_metadata_bytes<I>(lines: &mut I) -> Vec<Vec<u8>>
where
    I: Iterator<Item = Vec<u8>>,
{
    let mut metadata_lines = Vec::new();
    for line in lines {
        if line.starts_with(b"$$$$") {
            break;
        }
        metadata_lines.push(line);
    }
    metadata_lines
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
