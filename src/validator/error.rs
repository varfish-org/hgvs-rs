//! Error type definition.

use thiserror::Error;

/// Error type for validation of HGVS expressions.
#[derive(Error, Debug)]
pub enum Error {
    #[error("ref or alt must be non-empty in {0}")]
    RefOrAltMustBeNonEmpty(String),
    #[error("number of deleted bases must be positive in {0}")]
    NumDelBasesNotPositive(String),
    #[error("number of alternative bases must be positive in {0}")]
    NumAltBasesEmpty(String),
    #[error("number of inverted bases must be positive in {0}")]
    NumInvBasesNotPositive(String),

    #[error("Length implied by coordinates must equal count: {0}")]
    ImpliedLengthMismatch(String),
    #[error("start must be >=1 in {0}")]
    StartMustBePositive(String),
    #[error("end must be >=1 in {0}")]
    EndMustBePositive(String),
    #[error("sart <= end must hold in {0}")]
    StartMustBeLessThanEnd(String),
}
