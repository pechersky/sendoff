//! Record-name objects cached across blocks, so blocks repeating a file's
//! record names reuse the same Python objects.

use std::cell::{RefCell, RefMut};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use pyo3::prelude::*;

use crate::python::kind::Kind;

/// Names cached per kind before the cache is emptied, which bounds what
/// one-off names can hold on to.
const MAX_NAMES: usize = 4096;
/// Longer names are rare and are not cached.
const MAX_CACHED_NAME: usize = 64;

/// A fast hash for short record names, eight bytes at a time. Not
/// rustc-hash's FxHasher, which ran record operations about 0.5% slower on
/// EC2 with the same instruction count.
#[derive(Default)]
struct NameHasher(u64);

impl Hasher for NameHasher {
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut word = [0; 8];
            word[..chunk.len()].copy_from_slice(chunk);
            self.0 = (self.0.rotate_left(5) ^ u64::from_le_bytes(word))
                .wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }

    fn finish(&self) -> u64 {
        // Fold the well-mixed high bits into the low bits buckets use.
        self.0 ^ (self.0 >> 32)
    }
}

/// Marks no cached name: the start of a block, or an unknown successor.
const NO_NAME: u32 = u32::MAX;

/// A cached record-name object.
struct Name {
    bytes: Box<[u8]>,
    object: Py<PyAny>,
    /// The name that followed this one when last seen, or `NO_NAME`.
    next: u32,
}

/// Cached record names of one kind. Each links to the name that last
/// followed it, so a block repeating an earlier record order finds its names
/// by comparing bytes instead of hashing them, even when some records are
/// missing or new.
struct Names {
    index: HashMap<Box<[u8]>, u32, BuildHasherDefault<NameHasher>>,
    names: Vec<Name>,
    /// The first name of the last block.
    first: u32,
}

impl Default for Names {
    fn default() -> Self {
        Names {
            index: HashMap::default(),
            names: Vec::new(),
            first: NO_NAME,
        }
    }
}

impl Names {
    /// Return the index of the name `bytes`, caching a new object for it.
    /// A full cache is emptied first, which resets `prev`.
    fn find(&mut self, py: Python<'_>, kind: Kind, bytes: &[u8], prev: &mut u32) -> PyResult<u32> {
        if let Some(&found) = self.index.get(bytes) {
            return Ok(found);
        }
        let object = kind.value(py, bytes)?.unbind();
        if self.names.len() >= MAX_NAMES {
            // Releasing exact str and bytes objects runs no Python code, so
            // it cannot re-enter the cache.
            self.index.clear();
            self.names.clear();
            self.first = NO_NAME;
            *prev = NO_NAME;
        }
        let found = self.names.len() as u32;
        self.index.insert(bytes.into(), found);
        self.names.push(Name {
            bytes: bytes.into(),
            object,
            next: NO_NAME,
        });
        Ok(found)
    }
}

thread_local! {
    /// Record-name objects by their bytes, for text and then bytes names.
    static NAMES: RefCell<[Names; 2]> = RefCell::default();
}

/// Record-name objects of one kind, from the cache when it is available.
pub(super) struct Namer<'a> {
    kind: Kind,
    cache: Option<RefMut<'a, [Names; 2]>>,
    /// The previous cached name in this block, or `NO_NAME` at its start.
    prev: u32,
}

impl Namer<'_> {
    /// Call `f` with a `Namer` for `kind` that reuses name objects across
    /// blocks, or creates new ones in a re-entrant call that finds the cache
    /// in use.
    pub(super) fn with<R>(kind: Kind, f: impl FnOnce(&mut Namer<'_>) -> R) -> R {
        NAMES.with(|cache| {
            f(&mut Namer {
                kind,
                cache: cache.try_borrow_mut().ok(),
                prev: NO_NAME,
            })
        })
    }

    /// Return the object for the next record name of a block, read in order.
    pub(super) fn name<'py>(
        &mut self,
        py: Python<'py>,
        bytes: &[u8],
    ) -> PyResult<Bound<'py, PyAny>> {
        let Some(cache) = self.cache.as_deref_mut() else {
            return self.kind.value(py, bytes);
        };
        if bytes.len() > MAX_CACHED_NAME {
            return self.kind.value(py, bytes);
        }
        let names = &mut cache[self.kind as usize];
        let predicted = match self.prev {
            NO_NAME => names.first,
            prev => names.names[prev as usize].next,
        };
        if let Some(name) = names.names.get(predicted as usize)
            && *name.bytes == *bytes
        {
            self.prev = predicted;
            return Ok(name.object.bind(py).clone());
        }
        let found = names.find(py, self.kind, bytes, &mut self.prev)?;
        match self.prev {
            NO_NAME => names.first = found,
            prev => names.names[prev as usize].next = found,
        }
        self.prev = found;
        Ok(names.names[found as usize].object.bind(py).clone())
    }

    /// Return the object for a record name read out of block order.
    pub(super) fn lookup<'py>(
        &mut self,
        py: Python<'py>,
        bytes: &[u8],
    ) -> PyResult<Bound<'py, PyAny>> {
        let Some(cache) = self.cache.as_deref_mut() else {
            return self.kind.value(py, bytes);
        };
        if bytes.len() > MAX_CACHED_NAME {
            return self.kind.value(py, bytes);
        }
        let names = &mut cache[self.kind as usize];
        let mut unlinked = NO_NAME;
        let found = names.find(py, self.kind, bytes, &mut unlinked)?;
        Ok(names.names[found as usize].object.bind(py).clone())
    }
}
