//! Compressed input, recognized by magic numbers no SD title line starts
//! with, so callers decompress it first rather than get garbage blocks.

/// A compression format: its name and how to decompress it.
pub struct Compression {
    pub name: &'static str,
    pub decompress: &'static str,
}

/// The compression format whose magic number `data` starts with, if any.
pub fn compression(data: &[u8]) -> Option<Compression> {
    let (name, decompress) = match data {
        [0x1f, 0x8b, ..] => ("gzip", "gzip.open() or zcat"),
        // "BZh", the block size, then the first block's magic.
        [b'B', b'Z', b'h', b'1'..=b'9', block @ ..] if block.starts_with(b"1AY&SY") => {
            ("bzip2", "bz2.open() or bzcat")
        }
        [0xfd, b'7', b'z', b'X', b'Z', 0, ..] => ("xz", "lzma.open() or xzcat"),
        [0x28, 0xb5, 0x2f, 0xfd, ..] => ("zstd", "compression.zstd.open() or zstdcat"),
        _ => return None,
    };
    Some(Compression { name, decompress })
}

#[cfg(test)]
mod tests {
    use super::compression;

    #[test]
    fn names_compressed_input() {
        let name = |data: &[u8]| compression(data).map(|found| found.name);
        assert_eq!(name(b"\x1f\x8b\x08\0rest"), Some("gzip"));
        assert_eq!(name(b"BZh91AY&SYrest"), Some("bzip2"));
        assert_eq!(name(b"\xfd7zXZ\0rest"), Some("xz"));
        assert_eq!(name(b"\x28\xb5\x2f\xfdrest"), Some("zstd"));
    }

    #[test]
    fn passes_sd_titles() {
        for title in [
            &b""[..],
            b"\x1f",
            b"\n",
            b"BZh9 is a title\n",
            b"BZh01AY&SY\n",
            b"\xfd7zXZ\n",
            b"aspirin\n",
        ] {
            assert!(compression(title).is_none(), "{title:?}");
        }
    }
}
