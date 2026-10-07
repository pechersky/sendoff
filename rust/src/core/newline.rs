//! Line endings: which ending lines a block writes itself get, and the
//! `$$$$` line that ends a serialized block.

/// The delimiter lines that end serialized blocks, as read.
pub const TRAILER: &[u8] = b"$$$$\n";
pub const CRLF_TRAILER: &[u8] = b"$$$$\r\n";

/// `\r\n` if `line` ends with it, else `\n`.
pub fn newline(line: &[u8]) -> &'static [u8] {
    if line.ends_with(b"\r\n") {
        b"\r\n"
    } else {
        b"\n"
    }
}

/// The ending of `data`'s first line once it ends in `\n`, which lines a
/// block writes itself use.
pub fn title_newline(data: &[u8]) -> &'static [u8] {
    match memchr::memchr(b'\n', data) {
        Some(index) => newline(&data[..=index]),
        None if data.ends_with(b"\r") => b"\r\n",
        None => b"\n",
    }
}

/// What ends `line` if it lacks a `\n`: `\n` after a `\r`, else `newline`.
pub fn missing_ending(line: &[u8], newline: &'static [u8]) -> &'static [u8] {
    if line.ends_with(b"\r") {
        b"\n"
    } else {
        newline
    }
}

/// The `$$$$` line ending in `newline`.
pub fn trailer(newline: &[u8]) -> &'static [u8] {
    if newline == b"\r\n" {
        CRLF_TRAILER
    } else {
        TRAILER
    }
}
