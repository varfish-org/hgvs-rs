//! Code supporting mapping between coordinate systems.

pub mod alignment;
pub(crate) mod altseq;
pub mod assembly;
pub mod cigar;
mod error;
pub mod variant;

pub use error::Error;
