pub mod ctable;
pub mod framing;
pub mod indices;
pub mod sddata;

pub fn whitespace(c: char) -> bool {
    c.is_whitespace() || matches!(c, '\u{1c}'..='\u{1f}')
}
