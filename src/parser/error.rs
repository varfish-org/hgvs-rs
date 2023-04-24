//! Error type definition.

use thiserror::Error;

/// Error type for parsing of HGVS expressions.
#[derive(Error, Debug)]
pub enum Error {
    /// Invalid genome interval.
    #[error("{0} is not a valid genome interval")]
    InvalidGenomeInterval(String),
    /// Invalid transcript interval.
    #[error("{0} is not a valid tx interval")]
    InvalidTxInterval(String),
    /// Invalid CDS interval.
    #[error("{0} is not a valid CDS interval")]
    InvalidCdsInterval(String),
    /// Invalid HGVS expression.
    #[error("{0} is not a valid HGVS expression interval")]
    InvalidHgvsVariant(String),

    /// Ill-defined conversion.
    #[error("conversion of interval with different offsets (CDS start/end) is ill-defined: {0}")]
    IllDefinedConversion(String),
    /// Cannot None position into range.
    #[error("cannot convert interval with None position into range: {0}")]
    CannotNonePositionIntoRange(String),

    #[error("ref or alt must be non-empty in: {0}")]
    RefOrAltMustBeNonEmpty(String),
    #[error("number of deleted bases must be positive in: {0}")]
    NumDelBasesNotPositive(String),
    #[error("alternate bases must be non-empty in: {0}")]
    NumAltBasesEmpty(String),
    #[error("number of inverted bases must be positive in: {0}")]
    NumInvBasesNotPositive(String),
}
