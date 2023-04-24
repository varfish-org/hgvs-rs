//! Error type definition.

use thiserror::Error;

/// Error type for data.
#[derive(Error, Debug)]
pub enum Error {
    // #[error("validation error")]
    // ValidationFailed(#[from] crate::validator::Error),
    // #[error("expected a GenomeVariant but received {0}")]
    // ExpectedGenomeVariant(String),
}
