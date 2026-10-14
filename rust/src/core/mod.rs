//! Pure Rust algorithms, independent of PyO3.

pub mod block;
pub mod compressed;
pub mod ctable;
pub mod index;
pub mod lines;
pub mod mapped;
pub mod newline;
pub mod record;
pub mod split;
pub mod strip;
#[cfg(test)]
mod testing;
