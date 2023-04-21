//! Error type definition.

use thiserror::Error;

/// Error type for variant mapping.
#[derive(Error, Debug)]
pub enum Error {
    #[error("validation error")]
    ValidationFailed(#[from] crate::validator::Error),
    #[error("normalization error")]
    NormalizationFailed(#[from] crate::normalizer::Error),
    #[error("parsing failed")]
    ParsingFailed(#[from] crate::parser::Error),
    #[error("expected a GenomeVariant but received {0}")]
    ExpectedGenomeVariant(String),
    #[error("expected a TxVariant but received {0}")]
    ExpectedTxVariant(String),
    #[error("expected a CdsVariant but received {0}")]
    ExpectedCdsVariant(String),
    #[error("no NAEdit in HGVS.c variant: {0}")]
    NoNAEditInHgvsC(String),
    #[error("must have ProtVariant")]
    NotProtVariant,
    #[error("could not construct HGVS.p variant")]
    ProtVariantConstructionFailed,
    #[error("cannot get altered sequence for missing positions")]
    NoAlteredSequenceForMissingPositions,
    #[error("variant is missing nucleic acid edit")]
    NaEditMissing,
    #[error("can only update reference for c, g, m, n, r")]
    CannotUpdateReference,
    #[error("invalid CIGAR value: {0}")]
    InvalidCigarValue(char),
    #[error("invalid CIGAR value: {0}")]
    InvalidCigarCount(String),
    #[error("invalid CIGAR : {0}")]
    InvalidCigarOp(String),
}
