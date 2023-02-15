//! Data structures for representing HGVS variant descriptions.

/// Expression of "maybe uncertain".
#[derive(Clone, Debug, PartialEq)]
pub enum Mu<T> {
    /// Certain variant of `T`.
    Certain(T),
    /// Uncertain variant of `T`.
    Uncertain(T),
}

/// Representation of gene symbol, e.g., `TTN` or `Ttn`.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneSymbol {
    pub value: String,
}

/// Edit of nucleic acids.
#[derive(Clone, Debug, PartialEq)]
pub enum NaEdit {
    /// A substitution where both reference and alternative allele are nucleic acid strings
    /// (or empty).
    RefAlt {
        reference: String,
        alternative: String,
    },
    /// A substitution where the reference is a number and alternative is a count.
    NumAlt { count: i32, alternative: String },
    /// Insertion of one or more nucleic acid characters.
    Ins { alternative: String },
    /// Duplication of nucleic acid reference sequence.
    Dup { reference: String },
    /// Inversion of a (potentially empty) nucleic acid reference sequence.
    InvRef { reference: String },
    /// Inversion of a stretch given by its length.
    InvNum { count: i32 },
}

/// Uncertain change through extension.
#[derive(Clone, Debug, PartialEq)]
pub enum UncertainLengthChange {
    None,
    Unknown,
    Known(i32),
}

/// Representation of accession, e.g., `NM_01234.5`.
#[derive(Clone, Debug, PartialEq)]
pub struct Accession {
    pub value: String,
}

/// Protein edit with interval end edit.
#[derive(Clone, Debug, PartialEq)]
pub enum ProteinEdit {
    Fs {
        alternative: Option<String>,
        terminal: Option<String>,
        length: UncertainLengthChange,
    },
    Ext {
        /// Amino acid before "ext"
        aa_ext: Option<String>,
        /// Amino acid after "ext", terminal if shift is positive.
        ext_aa: Option<String>,
        /// Change in protein length.
        change: UncertainLengthChange,
    },
    Subst {
        alternative: String,
    },
    /// `delins`
    DelIns {
        alternative: String,
    },
    /// `ins`
    Ins {
        alternative: String,
    },
    /// `del`
    Del,
    /// `dup`
    Dup,
    /// `=`
    Ident,
}

/// A HGVS variant specification.
#[derive(Clone, Debug, PartialEq)]
pub enum HgvsVariant {
    /// Variant specification with `c.` position.
    CdsVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: CdsPosEdit,
    },
    /// Variant specification with `g.` position.
    GenomeVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: GenomePosEdit,
    },
    /// Variant specification with `m.` position.
    MtVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: MtPosEdit,
    },
    /// Variant specification with `n.` position.
    TxVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: TxPosEdit,
    },
    /// Variant specification with `p.` position.
    ProtVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: ProtPosEdit,
    },
    /// Variant specification with `r.` position.
    RnaVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: RnaPosEdit,
    },
}

/// Coding sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct CdsPosEdit {
    /// Interval on the CDS.
    pub pos: Mu<CdsInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// CDS position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct CdsInterval {
    /// Start position
    pub begin: CdsPos,
    /// End position
    pub end: CdsPos,
}

/// Specifies whether the CDS position is relative to the CDS start or
/// CDS end.
#[derive(Clone, Debug, PartialEq)]
pub enum CdsFrom {
    CdsStart,
    CdsEnd,
}

/// CDS position.
#[derive(Clone, Debug, PartialEq)]
pub struct CdsPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
    /// Whether starts at CDS start or end.
    pub cds_from: CdsFrom,
}

/// Genome sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomePosEdit {
    /// Interval on the genome.
    pub pos: Mu<GenomeInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// Genome position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeInterval {
    /// Start position
    pub begin: Option<i32>,
    /// End position
    pub end: Option<i32>,
}

/// Mitochondrial sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct MtPosEdit {
    /// Interval on the mitochondrium.
    pub pos: Mu<MtInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// Mitochondrial position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct MtInterval {
    /// Start position
    pub begin: Option<i32>,
    /// End position
    pub end: Option<i32>,
}

/// Transcript sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct TxPosEdit {
    /// Interval on a transcript.
    pub pos: Mu<TxInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// Transcript position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct TxInterval {
    /// Start position
    pub begin: TxPos,
    /// End position
    pub end: TxPos,
}

/// Transcript position.
#[derive(Clone, Debug, PartialEq)]
pub struct TxPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
}

/// RNA sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct RnaPosEdit {
    /// Interval on a transcript.
    pub pos: Mu<RnaInterval>,
    /// RNA change description.
    pub edit: Mu<NaEdit>,
}

/// RNA position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct RnaInterval {
    /// Start position
    pub begin: RnaPos,
    /// End position
    pub end: RnaPos,
}

/// RNA position.
#[derive(Clone, Debug, PartialEq)]
pub struct RnaPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
}

/// Protein sequence position with edit or special.
#[derive(Clone, Debug, PartialEq)]
pub enum ProtPosEdit {
    Ordinary {
        pos: Mu<ProtInterval>,
        edit: Mu<ProteinEdit>,
    },
    /// `=`
    NoChange,
    /// `(=)`
    NoChangeUncertain,
    /// `0`
    NoProtein,
    /// `0?`
    NoProteinUncertain,
}

/// Protein position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct ProtInterval {
    /// Start position
    pub begin: ProtPos,
    /// End position
    pub end: ProtPos,
}

/// Protein position.
#[derive(Clone, Debug, PartialEq)]
pub struct ProtPos {
    /// Amino acid value.
    pub aa: String,
    /// Number of `aa`.
    pub number: i32,
}
