//! Split a memory-mapped file into blocks.
//!
//! Mapping skips the read calls and the buffer copies that chunked reads
//! need. Pages behind the cursor are dropped as it advances, so resident
//! memory stays bounded however large the file.

use std::collections::VecDeque;
use std::fs::File;

use memmap2::Mmap;

use crate::core::block::Block;
use crate::core::split::split_block;

/// Input kept mapped in behind the cursor before older pages are dropped.
///
/// An 8 MiB window parsed no faster but added its size to peak RSS.
const RESIDENT: usize = 1 << 20;

/// Blocks split from a mapped file.
pub struct MappedBlocks {
    map: Mmap,
    /// Start of the next block.
    pos: usize,
    /// Pages before this offset have been dropped.
    dropped: usize,
    text: bool,
}

impl MappedBlocks {
    /// Map `file` if it is on a filesystem where mapping beats reading.
    pub fn open(file: &File, text: bool) -> Option<Self> {
        if is_local(file) {
            Self::new(file, text)
        } else {
            None
        }
    }

    /// Map `file`, or return `None` unless it is a nonempty regular file
    /// that maps: pipes, devices, and procfs files are better read.
    ///
    /// `text` reads universal newlines, as `split_block` does.
    pub fn new(file: &File, text: bool) -> Option<Self> {
        let metadata = file.metadata().ok()?;
        if !metadata.is_file() || metadata.len() == 0 {
            return None;
        }
        #[allow(unsafe_code)]
        // SAFETY: the map is read-only and private to this struct, which
        // only copies blocks out of it. As with any mapped file, another
        // process truncating the file raises SIGBUS on the next read of a
        // page past its new end.
        let map = unsafe { Mmap::map(file) }.ok()?;
        #[cfg(unix)]
        let _ = map.advise(memmap2::Advice::Sequential);
        Some(Self {
            map,
            pos: 0,
            dropped: 0,
            text,
        })
    }

    /// Whether every block has been split.
    pub fn is_done(&self) -> bool {
        self.pos == self.map.len()
    }

    /// Split the next block into `out`, and more until `len` bytes are split
    /// or input ends.
    pub fn split(&mut self, len: usize, out: &mut VecDeque<Block>) {
        let stop = self.pos.saturating_add(len);
        while !self.is_done() {
            let (block, taken) = split_block(&self.map[self.pos..], self.text);
            out.push_back(block);
            self.pos += taken;
            if self.pos >= stop {
                break;
            }
        }
        self.drop_behind();
    }

    #[cfg(unix)]
    fn drop_behind(&mut self) {
        // Aligned for any page size up to 64 KiB.
        let upto = self.pos.saturating_sub(RESIDENT) & !0xFFFF;
        if upto < self.dropped + RESIDENT {
            return;
        }
        #[allow(unsafe_code)]
        // SAFETY: dropped pages of a read-only file mapping read back from
        // the file, and blocks behind the cursor were already copied out.
        let dropped = unsafe {
            self.map.unchecked_advise_range(
                memmap2::UncheckedAdvice::DontNeed,
                self.dropped,
                upto - self.dropped,
            )
        };
        if dropped.is_ok() {
            self.dropped = upto;
        }
    }

    #[cfg(not(unix))]
    fn drop_behind(&mut self) {}
}

/// Whether `file` is on a local disk or memory filesystem.
///
/// Mapping parsed faster than reading on ext4 and tmpfs, but 5-6x slower
/// over the 9p mount WSL gives Windows drives, where page faults are round
/// trips; network and FUSE filesystems risk the same, so they are read.
#[cfg(target_os = "linux")]
fn is_local(file: &File) -> bool {
    use std::os::fd::AsRawFd;

    // The constants and f_type differ in type across targets, not value.
    const LOCAL: [u32; 6] = [
        libc::BTRFS_SUPER_MAGIC as u32,
        libc::EXT4_SUPER_MAGIC as u32,
        libc::F2FS_SUPER_MAGIC as u32,
        libc::OVERLAYFS_SUPER_MAGIC as u32,
        libc::TMPFS_MAGIC as u32,
        libc::XFS_SUPER_MAGIC as u32,
    ];
    let mut stat = std::mem::MaybeUninit::<libc::statfs>::uninit();
    #[allow(unsafe_code)]
    // SAFETY: fstatfs fills `stat` when it returns 0, and only then is it read.
    let kind = unsafe {
        if libc::fstatfs(file.as_raw_fd(), stat.as_mut_ptr()) != 0 {
            return false;
        }
        stat.assume_init().f_type
    };
    LOCAL.contains(&(kind as u32))
}

/// Mapping is only measured on Linux, so other systems read.
#[cfg(not(target_os = "linux"))]
fn is_local(_: &File) -> bool {
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn blocks(input: &[u8], text: bool, batch: usize) -> Option<Vec<Block>> {
        let mut file = tempfile();
        file.write_all(input).unwrap();
        let mut mapped = MappedBlocks::new(&file, text)?;
        let mut out = VecDeque::new();
        while !mapped.is_done() {
            let before = out.len();
            mapped.split(batch, &mut out);
            assert!(out.len() > before);
        }
        Some(out.into())
    }

    fn tempfile() -> File {
        let path = std::env::temp_dir().join(format!(
            "sendoff-mapped-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let file = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .unwrap();
        std::fs::remove_file(&path).unwrap();
        file
    }

    #[test]
    fn splits_in_batches_of_whole_blocks() {
        let input = b"a\r\nM  END\r\n$$$$\r\nb\r$$$$\rc";
        let whole = blocks(input, true, usize::MAX).unwrap();
        assert_eq!(whole.len(), 3);
        for batch in [0, 1, 5, 100] {
            assert_eq!(blocks(input, true, batch).unwrap(), whole);
        }
        assert_eq!(blocks(input, false, 1).unwrap().len(), 2);
    }

    #[test]
    fn leaves_empty_files_to_reads() {
        assert!(blocks(b"", false, 1).is_none());
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn maps_only_local_filesystems() {
        assert!(!is_local(&File::open("/proc/self/status").unwrap()));
        if let Ok(file) = File::open("/dev/shm") {
            assert!(is_local(&file));
        }
    }

    #[cfg(unix)]
    #[test]
    fn drops_pages_behind_the_cursor() {
        let block = [&[b'x'; 4000][..], b"\n$$$$\n"].concat();
        let input = block.repeat(4 * RESIDENT / block.len());
        let mut file = tempfile();
        file.write_all(&input).unwrap();
        let mut mapped = MappedBlocks::new(&file, false).unwrap();
        let mut out = VecDeque::new();
        let mut count = 0;
        while !mapped.is_done() {
            mapped.split(1 << 17, &mut out);
            count += out.drain(..).count();
            assert!(mapped.pos - mapped.dropped < 2 * RESIDENT + (1 << 17) + (1 << 16));
        }
        assert_eq!(count, input.len() / block.len());
        assert!(mapped.dropped >= 2 * RESIDENT);
    }
}
