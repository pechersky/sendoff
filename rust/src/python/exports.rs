//! Blocks that refuse edits while Python buffers export their bytes, and
//! index their records for repeated lookups.

use std::borrow::Cow;
use std::ffi::c_void;
use std::ops::Deref;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, OnceLock};

use pyo3::PyResult;
use pyo3::exceptions::PyBufferError;

use crate::core::block::Block;
use crate::core::index::Index;
use crate::core::record::{MissingFieldName, Records};

/// Tracks active Python buffer exports via an `Arc<()>` reference count,
/// allocated on the first export so releasing a buffer needs no borrow of
/// the `SDBlock`.
#[derive(Default)]
struct Exports(OnceLock<Arc<()>>);

impl Exports {
    fn is_active(&self) -> bool {
        self.0
            .get()
            .is_some_and(|exports| Arc::strong_count(exports) > 1)
    }

    fn acquire(&self) -> *mut c_void {
        let exports = self.0.get_or_init(Arc::default);
        Arc::into_raw(Arc::clone(exports)).cast_mut().cast()
    }

    #[allow(unsafe_code)]
    unsafe fn release(token: *mut c_void) {
        // SAFETY: the caller passes a pointer from `Arc::into_raw`, once.
        drop(unsafe { Arc::from_raw(token.cast_const().cast::<()>()) });
    }
}

/// A block whose bytes buffers may export, so edits, which may move or free
/// them, raise until every export is released. Edits also drop the index
/// of its records.
pub(crate) struct Exported {
    block: Block,
    exports: Exports,
    /// The records' index, built on the second lookup since the last edit,
    /// so blocks looked up once allocate none.
    index: OnceLock<Option<Index>>,
    looked_up: AtomicBool,
}

impl Exported {
    pub(crate) fn new(block: Block) -> Self {
        Self {
            block,
            exports: Exports::default(),
            index: OnceLock::new(),
            looked_up: AtomicBool::new(false),
        }
    }

    /// The index of the records for a lookup: `None` for the first since
    /// the last edit, or SD data too large to index.
    fn index(&self) -> Result<Option<&Index>, MissingFieldName> {
        if let Some(index) = self.index.get() {
            return Ok(index.as_ref());
        }
        if !self.looked_up.swap(true, Ordering::Relaxed) {
            return Ok(None);
        }
        let index = Index::build(self.block.metadata())?;
        Ok(self.index.get_or_init(|| index).as_ref())
    }

    /// The value of the last record named `name`, if any.
    pub(crate) fn last<'a>(
        &'a self,
        name: &[u8],
    ) -> Result<Option<Cow<'a, [u8]>>, MissingFieldName> {
        let metadata = self.block.metadata();
        if let Some(index) = self.index()? {
            return Ok(index.last(metadata, name));
        }
        let mut found = None;
        for header in Records::new(metadata, 0).headers() {
            let header = header?;
            if header.name == name {
                found = Some(header.start);
            }
        }
        let Some(start) = found else {
            return Ok(None);
        };
        let record = Records::new(metadata, start).next();
        let record = record.and_then(Result::ok).expect("a record starts here");
        Ok(Some(record.value))
    }

    /// Whether any record is named `name`.
    pub(crate) fn contains(&self, name: &[u8]) -> Result<bool, MissingFieldName> {
        let metadata = self.block.metadata();
        if let Some(index) = self.index()? {
            return Ok(index.contains(metadata, name));
        }
        let mut seen = false;
        for header in Records::new(metadata, 0).headers() {
            seen |= header?.name == name;
        }
        Ok(seen)
    }

    /// The number of records, duplicates included.
    pub(crate) fn record_count(&self) -> Result<usize, MissingFieldName> {
        if let Some(index) = self.index()? {
            return Ok(index.len());
        }
        let mut count = 0;
        for header in self.block.headers() {
            header?;
            count += 1;
        }
        Ok(count)
    }

    /// The block, to edit, unless a buffer exports it; drops the index.
    pub(crate) fn get_mut(&mut self) -> PyResult<&mut Block> {
        if self.exports.is_active() {
            return Err(PyBufferError::new_err(
                "Existing exports of data: SDBlock cannot be edited",
            ));
        }
        self.index.take();
        *self.looked_up.get_mut() = false;
        Ok(&mut self.block)
    }

    /// Count an export, returning the token `release` takes.
    pub(crate) fn export(&self) -> *mut c_void {
        self.exports.acquire()
    }

    /// Count an export released.
    ///
    /// # Safety
    ///
    /// `token` comes from `export`, and is released once.
    #[allow(unsafe_code)]
    pub(crate) unsafe fn release(token: *mut c_void) {
        // SAFETY: forwarded from caller.
        unsafe { Exports::release(token) };
    }
}

impl Deref for Exported {
    type Target = Block;

    fn deref(&self) -> &Block {
        &self.block
    }
}

/// A copy is a new block, not exported.
impl Clone for Exported {
    fn clone(&self) -> Self {
        Self::new(self.block.clone())
    }
}

impl PartialEq for Exported {
    fn eq(&self, other: &Self) -> bool {
        self.block == other.block
    }
}
