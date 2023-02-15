//! Data structures for representing HGVS variant descriptions.

// impl HgvsVariant {
//     pub fn parse(input: &str) -> IResult<&str, Self> {

//     }
// }

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
    MitochondrialVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: MitochondriumPosEdit,
    },
    /// Variant specification with `n.` position.
    TranscriptVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: TranscriptPosEdit,
    },
    /// Variant specification with `p.` position.
    ProteinVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: ProteinPosEdit,
    },
    /// Variant specification with `r.` position.
    RnaVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        pos_edit: RnaPosEdit,
    },
}

/// Representation of gene symbol, e.g., `TTN` or `Ttn`.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneSymbol {
    pub value: String,
}

/// RNA sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct RnaPosEdit {
    /// Interval on a transcript.
    pub pos: Mu<Interval>,
    /// RNA change description.
    pub edit: Mu<NaEdit>,
}

/// The interval type
#[derive(Clone, Debug, PartialEq)]
pub enum PosType {
    Cds,
    Genome,
    Mitochorial,
    Transcript,
    Protein,
    Rna,
}

/// Interval.
#[derive(Clone, Debug, PartialEq)]
pub struct Interval {
    /// The type of the position
    pub pos_type: PosType,
    /// Start position
    pub pos: i32,
    /// End position
    pub end: i32,
}

/// Representation of accession, e.g., `NM_01234.5`.
#[derive(Clone, Debug, PartialEq)]
pub struct Accession {
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
    NumAlt { count: u32, alternative: String },
    /// Insertion of one or more nucleic acid characters.
    Ins { alternative: String },
    /// Duplication of nucleic acid reference sequence.
    Dup { reference: String },
    /// Inversion of a (potentially empty) nucleic acid reference sequence.
    InvRef { reference: String },
    /// Inversion of a stretch given by its length.
    InvNum { count: u32 },
}

/// Expression of "maybe uncertain".
#[derive(Clone, Debug, PartialEq)]
pub enum Mu<T> {
    /// Certain variant of `T`.
    Certain(T),
    /// Uncertain variant of `T`.
    Uncertain(T),
}

/// Coding sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct CdsPosEdit {
    /// Interval on the CDS.
    pub pos: Mu<Interval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// Genome sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomePosEdit {
    /// Interval on the genome.
    pub pos: Mu<Interval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// Mitochondrial sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct MitochondriumPosEdit {
    /// Interval on the mitochondrium.
    pub pos: Mu<Interval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// Transcript sequence position with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct TranscriptPosEdit {
    /// Interval on a transcript.
    pub pos: Mu<Interval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

/// Protein sequence position with edit or special.
#[derive(Clone, Debug, PartialEq)]
pub enum ProteinPosEdit {
    Ordinary {
        pos: Mu<Interval>,
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

/// Uncertain change through extension.
#[derive(Clone, Debug, PartialEq)]
pub enum UncertainChange {
    None,
    Unknown,
    Known(i32),
}

/// Protein edit with interval end edit.
#[derive(Clone, Debug, PartialEq)]
pub enum ProteinEdit {
    Fs {
        alternative: Option<String>,
        terminal: Option<String>,
        length: UncertainChange,
    },
    Ext {
        /// Amino acid before "ext"
        aa_ext: Option<String>,
        /// Amino acid after "ext", terminal if shift is positive.
        ext_aa: Option<String>,
        /// Change in protein length.
        change: UncertainChange,
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
