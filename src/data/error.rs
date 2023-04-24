//! Error type definition.

use thiserror::Error;

/// Error type for data.
#[derive(Error, Debug)]
pub enum Error {
    #[error("UTA Postgres access error")]
    UtaPostgresError(#[from] postgres::Error),
    #[error("sequence operation failed")]
    SequenceOperationFailed(#[from] crate::sequences::Error),
    #[error("problem with seqrepo access")]
    SeqRepoError(#[from] seqrepo::Error),
    #[error("no tx_exons for tx_ac={0}, alt_ac={1}, alt_aln_method={2}")]
    NoTxExons(String, String, String),
    #[error("could not get parent from {0}")]
    PathParent(String),
    #[error("could not get basename from {0}")]
    PathBasename(String),
    #[error("could not open cdot JSON file: {0}")]
    CdotJsonOpen(String),
    #[error("could not parse cdot JSON file: {0}")]
    CdotJsonParse(String),
    #[error("no gene found for {0}")]
    NoGeneFound(String),
    #[error("no transcript found for {0}")]
    NoTranscriptFound(String),
    #[error("no alignment found for {0} to {1}")]
    NoAlignmentFound(String, String),
    #[error("found no sequence record for accession {0}")]
    NoSequenceRecord(String),
}
