//! Blocks that refuse edits while Python buffers export their bytes.

use std::ffi::c_void;
use std::ops::Deref;
use std::sync::{Arc, OnceLock};

use pyo3::PyResult;
use pyo3::exceptions::PyBufferError;

use crate::core::block::Block;

/// A block whose bytes buffers may export, so edits, which may move or free
/// them, raise until every export is released.
///
/// Each export holds a reference to `exports`, made on the first, so
/// releasing one needs no borrow of the block, which a mutable borrow
/// running Python code (a GC, say) could be holding.
pub(crate) struct Exported {
    block: Block,
    exports: OnceLock<Arc<()>>,
}

impl Exported {
    pub(crate) fn new(block: Block) -> Self {
        Self {
            block,
            exports: OnceLock::new(),
        }
    }

    /// The block, to edit, unless a buffer exports it.
    pub(crate) fn get_mut(&mut self) -> PyResult<&mut Block> {
        if self
            .exports
            .get()
            .is_some_and(|exports| Arc::strong_count(exports) > 1)
        {
            return Err(PyBufferError::new_err(
                "Existing exports of data: SDBlock cannot be edited",
            ));
        }
        Ok(&mut self.block)
    }

    /// Count an export, returning the token `release` takes.
    pub(crate) fn export(&self) -> *mut c_void {
        let exports = self.exports.get_or_init(Arc::default);
        Arc::into_raw(Arc::clone(exports)).cast_mut().cast()
    }

    /// Count an export released.
    ///
    /// # Safety
    ///
    /// `token` comes from `export`, and is released once.
    #[allow(unsafe_code)]
    pub(crate) unsafe fn release(token: *mut c_void) {
        // SAFETY: the caller passes a pointer from `Arc::into_raw`, once.
        drop(unsafe { Arc::from_raw(token.cast_const().cast::<()>()) });
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
