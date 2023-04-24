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
    #[error("sequence operation failed")]
    SequenceOperationFailed(#[from] crate::sequences::Error),
    #[error("problem accessing data")]
    DataError(#[from] crate::data::error::Error),
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
    #[error("invalid CIGAR op: {0}")]
    InvalidCigarOp(String),
    #[error("invalid CIGAR string: {0}")]
    InvalidCigarString(String),
    #[error(
        "position is beyond the bounds of transcript record (pos={0}, from_pos={1}, to_pos={2})"
    )]
    PositionBeyondTranscriptBounds(i32, String, String),
    #[error("algorithm error in CIGAR mapper")]
    CigarMapperError,
    #[error("not a GenomeVariant: {0}")]
    NotGenomeVariant(String),
    #[error("no alignments for {0} in {1} using {2}")]
    NoAlignments(String, String, String),
    #[error(
        "multiple chromosome alignments for {0} in {1} using {2} (non- \
        pseudoautosomal region) [{3}]"
    )]
    MultipleChromAlignsNonPar(String, String, String, String),
    #[error(
        "multiple chromosome alignments for {0} in {1} using {2} (likely \
        pseudoautosomal region)"
    )]
    MultipleChromAlignsLikelyPar(String, String, String),
    #[error(
        "multiple chromosome alignments for {0} in {1} using {2} \
        (in_par_assume={3} select {4} of them)"
    )]
    MultipleChromAlignsInParAssume(String, String, String, String, usize),
    #[error(
        "transcript {0} is not supported because its sequence length of
        {1} is not a multiple of 3"
    )]
    TranscriptLengthInvalid(String, usize),
    #[error("start pos ouf of range in reference sequence")]
    StartPosOutOfRange,
    #[error("got multiple AA variants which is not supported")]
    MultipleAAVariants,
    #[error("deletion sequence should not be empty")]
    DeletionSequenceEmpty,
    #[error("insertion sequence should not be empty")]
    InsertionSequenceEmpty,
    #[error("cannot build CIGAR string from empty exons")]
    EmptyExons,
    #[error("found no exons for tx_ac={0}, alt_ac={1}, alt_aln_method={2}")]
    NoExons(String, String, String),
    #[error("non-adjacent exons for tx_ac={0}, alt_ac={1}, alt_aln_method={2}: {3}")]
    NonAdjacentExons(String, String, String, String),
    #[error("CDS start and end must both be defined or undefined")]
    InconsistentCdsStartEnd,
    #[error("cannot project genome interval with missing start or end position: {0}")]
    MissingGenomeIntervalPosition(String),
    #[error("CDS is undefined for {0}; cannot map to c. coordinates (non-coding transcript?)")]
    CdsUndefined(String),
    #[error("coordinate is outside the bounds of the reference sequence")]
    CoordinateOutsideReference,
    #[error("c.{0} coordinate is out of bounds")]
    CoordinateOutOfBounds(String),
    #[error("cannot convert interval start: {0} to usize")]
    CannotConvertIntervalStart(i32),
    #[error("cannot convert interval end: {0} to usize")]
    CannotConvertIntervalEnd(i32),
}
