//! Data structures for representing HGVS variant descriptions.

use std::ops::{Deref, Range};

use crate::parser::error::Error;
use log::warn;
use crate::Sequence;

/// Expression of "maybe uncertain".
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Mu<T> {
    /// Certain variant of `T`.
    Certain(T),
    /// Uncertain variant of `T`.
    Uncertain(T),
}

impl<T> Mu<T> {
    pub fn from(value: T, is_certain: bool) -> Self {
        if is_certain {
            Mu::Certain(value)
        } else {
            Mu::Uncertain(value)
        }
    }

    pub fn is_certain(&self) -> bool {
        match &self {
            Mu::Certain(_) => true,
            Mu::Uncertain(_) => false,
        }
    }

    pub fn unwrap(self) -> T {
        match self {
            Mu::Certain(value) => value,
            Mu::Uncertain(value) => value,
        }
    }

    pub fn inner(&self) -> &T {
        match self {
            Mu::Certain(value) => value,
            Mu::Uncertain(value) => value,
        }
    }

    pub fn inner_mut(&mut self) -> &mut T {
        match self {
            Mu::Certain(value) => value,
            Mu::Uncertain(value) => value,
        }
    }
}

/// Representation of gene symbol, e.g., `TTN` or `Ttn`.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GeneSymbol {
    pub value: String,
}

impl GeneSymbol {
    pub fn new(value: &str) -> Self {
        Self {
            value: value.to_string(),
        }
    }

    pub fn from(value: String) -> Self {
        Self { value }
    }
}

impl Deref for GeneSymbol {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.value
    }
}

/// Edit of nucleic acids.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum NaEdit {
    /// A substitution where both reference and alternative allele are nucleic acid strings
    /// (or empty).
    RefAlt {
        reference: Sequence,
        alternative: Sequence,
    },
    /// A substitution where the reference is a number and alternative is a count.
    NumAlt { count: i32, alternative: Sequence },
    /// Deletion of one or more nucleic acid characters.
    DelRef { reference: Sequence },
    /// Deletion of a number of characters.
    DelNum { count: i32 },
    /// Insertion of one or more nucleic acid characters.
    Ins { alternative: Sequence },
    /// Duplication of nucleic acid reference sequence.
    Dup { reference: Sequence },
    /// Inversion of a (potentially empty) nucleic acid reference sequence.
    InvRef { reference: Sequence },
    /// Inversion of a stretch given by its length.
    InvNum { count: i32 },
}

impl NaEdit {
    /// Returns whether the edit has a count (and not reference sequence).
    pub fn is_na_edit_num(&self) -> bool {
        match self {
            NaEdit::RefAlt { .. }
            | NaEdit::DelRef { .. }
            | NaEdit::Ins { .. }
            | NaEdit::Dup { .. }
            | NaEdit::InvRef { .. } => false,
            NaEdit::NumAlt { .. } | NaEdit::DelNum { .. } | NaEdit::InvNum { .. } => true,
        }
    }

    /// Returns whether the edit is an insertion.
    pub fn is_ins(&self) -> bool {
        matches!(self, NaEdit::Ins { .. })
    }

    /// Returns whether the edit is a duplication.
    pub fn is_dup(&self) -> bool {
        matches!(self, NaEdit::Dup { .. })
    }

    /// Ensures that the reference is a count and no reference bases.
    pub fn with_num(&self) -> Self {
        match self {
            NaEdit::DelRef { reference } => NaEdit::DelNum {
                count: reference.len() as i32,
            },
            NaEdit::InvRef { reference } => NaEdit::InvNum {
                count: reference.len() as i32,
            },
            NaEdit::RefAlt {
                reference,
                alternative,
            } => NaEdit::NumAlt {
                count: reference.len() as i32,
                alternative: alternative.clone(),
            },
            NaEdit::NumAlt { .. }
            | NaEdit::DelNum { .. }
            | NaEdit::InvNum { .. }
            | NaEdit::Ins { .. }
            | NaEdit::Dup { .. } => self.clone(),
        }
    }

    /// Return whether the reference equals the given value.
    pub fn reference_equals(&self, value: &[u8]) -> bool {
        match self {
            NaEdit::RefAlt { reference, .. }
            | NaEdit::DelRef { reference }
            | NaEdit::Dup { reference }
            | NaEdit::InvRef { reference } => reference == value,
            _ => false,
        }
    }

    /// Return an updated `NaEdit` that has the reference replaced with the given sequence.
    pub fn with_reference(self, reference: Sequence) -> Self {
        match self {
            NaEdit::RefAlt {
                alternative,
                reference: old_reference,
            } => {
                if old_reference.is_empty() && alternative.is_empty() {
                    NaEdit::RefAlt {
                        // sic!
                        alternative: reference.clone(), // sic!
                        reference,                      // sic!
                    }
                } else {
                    NaEdit::RefAlt {
                        reference,
                        alternative,
                    }
                }
            }
            NaEdit::NumAlt { alternative, .. } => NaEdit::RefAlt {
                reference,
                alternative,
            },
            NaEdit::DelRef { .. } => NaEdit::DelRef { reference },
            NaEdit::DelNum { .. } => NaEdit::DelRef { reference },
            NaEdit::Ins { alternative } => {
                warn!("Calling with_reference() on NaEdit::Ins");
                NaEdit::Ins { alternative }
            }
            NaEdit::Dup { .. } => NaEdit::Dup { reference },
            NaEdit::InvRef { .. } => NaEdit::InvRef { reference },
            NaEdit::InvNum { .. } => NaEdit::InvRef { reference },
        }
    }
}

/// Uncertain change through extension.
#[derive(Clone, Debug, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub enum UncertainLengthChange {
    #[default]
    None,
    Unknown,
    Known(i32),
}

/// Representation of accession, e.g., `NM_01234.5`.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Accession {
    pub value: String,
}

impl Deref for Accession {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.value
    }
}

impl Accession {
    pub fn new(value: &str) -> Self {
        Self {
            value: value.to_string(),
        }
    }

    pub fn from(value: String) -> Self {
        Self { value }
    }
}

/// Protein edit with interval end edit.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum ProteinEdit {
    Fs {
        alternative: Option<Sequence>,
        terminal: Option<Sequence>,
        length: UncertainLengthChange,
    },
    Ext {
        /// Amino acid before "ext"
        aa_ext: Option<Sequence>,
        /// Amino acid after "ext", terminal if shift is positive.
        ext_aa: Option<Sequence>,
        /// Change in protein length.
        change: UncertainLengthChange,
    },
    Subst {
        alternative: Sequence,
    },
    /// `delins`
    DelIns {
        alternative: Sequence,
    },
    /// `ins`
    Ins {
        alternative: Sequence,
    },
    /// `del`
    Del,
    /// `dup`
    Dup,
    /// `=`
    Ident,
}

/// A HGVS variant specification.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum HgvsVariant {
    /// Variant specification with `c.` location.
    CdsVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        loc_edit: CdsLocEdit,
    },
    /// Variant specification with `g.` location.
    GenomeVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        loc_edit: GenomeLocEdit,
    },
    /// Variant specification with `m.` location.
    MtVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        loc_edit: MtLocEdit,
    },
    /// Variant specification with `n.` location.
    TxVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        loc_edit: TxLocEdit,
    },
    /// Variant specification with `p.` location.
    ProtVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        loc_edit: ProtLocEdit,
    },
    /// Variant specification with `r.` location.
    RnaVariant {
        accession: Accession,
        gene_symbol: Option<GeneSymbol>,
        loc_edit: RnaLocEdit,
    },
}

impl HgvsVariant {
    /// Return whether has a nucleic acid change and that one uses counts.
    pub fn is_na_edit_num(&self) -> bool {
        match self {
            HgvsVariant::CdsVariant { loc_edit, .. } => loc_edit.edit.inner().is_na_edit_num(),
            HgvsVariant::GenomeVariant { loc_edit, .. } => loc_edit.edit.inner().is_na_edit_num(),
            HgvsVariant::MtVariant { loc_edit, .. } => loc_edit.edit.inner().is_na_edit_num(),
            HgvsVariant::TxVariant { loc_edit, .. } => loc_edit.edit.inner().is_na_edit_num(),
            HgvsVariant::RnaVariant { loc_edit, .. } => loc_edit.edit.inner().is_na_edit_num(),
            HgvsVariant::ProtVariant { .. } => false,
        }
    }

    /// Replace `NaEdit` having a reference with one having a count.
    pub fn with_na_ref_num(self) -> Self {
        match self {
            HgvsVariant::CdsVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::CdsVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_num(),
            },
            HgvsVariant::GenomeVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::GenomeVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_num(),
            },
            HgvsVariant::MtVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::MtVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_num(),
            },
            HgvsVariant::TxVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::TxVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_num(),
            },
            HgvsVariant::ProtVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::ProtVariant {
                accession,
                gene_symbol,
                loc_edit,
            },
            HgvsVariant::RnaVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::RnaVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_num(),
            },
        }
    }

    /// Replace reference sequence.
    pub fn with_reference(self, value: Sequence) -> Self {
        match self {
            HgvsVariant::CdsVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::CdsVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_reference(value),
            },
            HgvsVariant::GenomeVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::GenomeVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_reference(value),
            },
            HgvsVariant::MtVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::MtVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_reference(value),
            },
            HgvsVariant::TxVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::TxVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_reference(value),
            },
            HgvsVariant::ProtVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => {
                warn!("Calling with_reference on ProtVariant");
                HgvsVariant::ProtVariant {
                    accession,
                    gene_symbol,
                    loc_edit,
                }
            }
            HgvsVariant::RnaVariant {
                accession,
                gene_symbol,
                loc_edit,
            } => HgvsVariant::RnaVariant {
                accession,
                gene_symbol,
                loc_edit: loc_edit.with_reference(value),
            },
        }
    }

    /// Return the gene symbol.
    pub fn gene_symbol(&self) -> &Option<GeneSymbol> {
        match self {
            HgvsVariant::CdsVariant { gene_symbol, .. } => gene_symbol,
            HgvsVariant::GenomeVariant { gene_symbol, .. } => gene_symbol,
            HgvsVariant::MtVariant { gene_symbol, .. } => gene_symbol,
            HgvsVariant::TxVariant { gene_symbol, .. } => gene_symbol,
            HgvsVariant::ProtVariant { gene_symbol, .. } => gene_symbol,
            HgvsVariant::RnaVariant { gene_symbol, .. } => gene_symbol,
        }
    }

    /// Return the accession.
    pub fn accession(&self) -> &Accession {
        match self {
            HgvsVariant::CdsVariant { accession, .. } => accession,
            HgvsVariant::GenomeVariant { accession, .. } => accession,
            HgvsVariant::MtVariant { accession, .. } => accession,
            HgvsVariant::TxVariant { accession, .. } => accession,
            HgvsVariant::ProtVariant { accession, .. } => accession,
            HgvsVariant::RnaVariant { accession, .. } => accession,
        }
    }

    /// Return the 0-based range of the location, possibly wrapped into `Mu`
    pub fn mu_loc_range(&self) -> Option<Mu<Range<i32>>> {
        match self {
            HgvsVariant::CdsVariant { loc_edit, .. } => loc_edit
                .loc
                .inner()
                .clone()
                .try_into()
                .ok()
                .map(|l| Mu::from(l, loc_edit.loc.is_certain())),
            HgvsVariant::GenomeVariant { loc_edit, .. } => loc_edit
                .loc
                .inner()
                .clone()
                .try_into()
                .ok()
                .map(|l| Mu::from(l, loc_edit.loc.is_certain())),
            HgvsVariant::MtVariant { loc_edit, .. } => loc_edit
                .loc
                .inner()
                .clone()
                .try_into()
                .ok()
                .map(|l| Mu::from(l, loc_edit.loc.is_certain())),
            HgvsVariant::TxVariant { loc_edit, .. } => {
                Some(From::from(loc_edit.loc.inner().clone()))
                    .map(|l| Mu::from(l, loc_edit.loc.is_certain()))
            }
            HgvsVariant::RnaVariant { loc_edit, .. } => {
                Some(From::from(loc_edit.loc.inner().clone()))
                    .map(|l| Mu::from(l, loc_edit.loc.is_certain()))
            }
            HgvsVariant::ProtVariant {
                loc_edit: ProtLocEdit::Ordinary { loc, .. },
                ..
            } => Some(Mu::from(loc.inner().clone().into(), loc.is_certain())),
            _ => None,
        }
    }

    pub fn loc_range(&self) -> Option<Range<i32>> {
        self.mu_loc_range().map(|l| l.inner().clone())
    }

    /// Return the `NaEdit` wrapped in `Mu`, if any.
    pub fn mu_na_edit(&self) -> Option<&Mu<NaEdit>> {
        match self {
            HgvsVariant::CdsVariant { loc_edit, .. } => Some(&loc_edit.edit),
            HgvsVariant::GenomeVariant { loc_edit, .. } => Some(&loc_edit.edit),
            HgvsVariant::MtVariant { loc_edit, .. } => Some(&loc_edit.edit),
            HgvsVariant::TxVariant { loc_edit, .. } => Some(&loc_edit.edit),
            HgvsVariant::RnaVariant { loc_edit, .. } => Some(&loc_edit.edit),
            _ => None,
        }
    }

    /// Return the `NaEdit` if any.
    pub fn na_edit(&self) -> Option<&NaEdit> {
        self.mu_na_edit().map(|e| e.inner())
    }

    /// Return the `ProtLocEdit` if any.
    pub fn mu_prot_edit(&self) -> Option<&Mu<ProteinEdit>> {
        match self {
            HgvsVariant::ProtVariant {
                loc_edit: ProtLocEdit::Ordinary { edit, .. },
                ..
            } => Some(edit),
            _ => None,
        }
    }

    /// Return the `ProtLocEdit` wrapped in `Mu`, if any.
    pub fn prot_edit(&self) -> Option<&ProteinEdit> {
        self.mu_prot_edit().map(|e| e.inner())
    }

    /// Return whether start or end position is intronic (offset != 0).
    pub fn spans_intron(&self) -> bool {
        match self {
            HgvsVariant::CdsVariant { loc_edit, .. } => {
                loc_edit.loc.inner().start.offset.is_some()
                    || loc_edit.loc.inner().end.offset.is_some()
            }
            HgvsVariant::TxVariant { loc_edit, .. } => {
                loc_edit.loc.inner().start.offset.is_some()
                    || loc_edit.loc.inner().end.offset.is_some()
            }
            HgvsVariant::RnaVariant { loc_edit, .. } => {
                loc_edit.loc.inner().start.offset.is_some()
                    || loc_edit.loc.inner().end.offset.is_some()
            }
            _ => false,
        }
    }
}

/// Coding sequence location with edit.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CdsLocEdit {
    /// Location on the CDS.
    pub loc: Mu<CdsInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl CdsLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: Sequence) -> Self {
        CdsLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }

    /// Return the LocEdit and ensure that the `NaEdit` has a count and no reference bases.
    fn with_num(self) -> Self {
        CdsLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_num()),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_num()),
            },
        }
    }
}

/// CDS position interval.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CdsInterval {
    /// Start position
    pub start: CdsPos,
    /// End position
    pub end: CdsPos,
}

impl TryFrom<CdsInterval> for Range<i32> {
    type Error = Error;

    /// The CDS interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn try_from(value: CdsInterval) -> Result<Self, Self::Error> {
        if value.start.offset.is_some() || value.end.offset.is_some() {
            warn!("Converting interval {:?} with offset to range!", &value);
        }
        if value.start.cds_from != value.end.cds_from {
            Err(Error::IllDefinedConversion(format!("{:?}", value)))
        } else {
            Ok(if value.start.base > 0 {
                (value.start.base - 1)..value.end.base
            } else {
                value.start.base..(value.end.base + 1)
            })
        }
    }
}

/// Specifies whether the CDS position is relative to the CDS start or
/// CDS end.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum CdsFrom {
    Start,
    End,
}

/// CDS position.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CdsPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
    /// Whether starts at CDS start or end.
    pub cds_from: CdsFrom,
}

/// Genome sequence location with edit.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GenomeLocEdit {
    /// Location on the genome.
    pub loc: Mu<GenomeInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl GenomeLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: Sequence) -> Self {
        GenomeLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }

    /// Return the LocEdit and ensure that the `NaEdit` has a count and no reference bases.
    fn with_num(self) -> Self {
        GenomeLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_num()),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_num()),
            },
        }
    }
}

/// Genome position interval.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GenomeInterval {
    /// Start position
    pub start: Option<i32>,
    /// End position
    pub end: Option<i32>,
}

impl TryInto<Range<i32>> for GenomeInterval {
    type Error = Error;

    /// The genome interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn try_into(self) -> Result<Range<i32>, Self::Error> {
        if let (Some(start), Some(end)) = (self.start, self.end) {
            Ok((start - 1)..end)
        } else {
            Err(Error::CannotNonePositionIntoRange(format!("{:?}", self)))
        }
    }
}

/// Mitochondrial sequence location with edit.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MtLocEdit {
    /// Location on the mitochondrium.
    pub loc: Mu<MtInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl MtLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: Sequence) -> Self {
        MtLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }

    /// Return the LocEdit and ensure that the `NaEdit` has a count and no reference bases.
    fn with_num(self) -> Self {
        MtLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_num()),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_num()),
            },
        }
    }
}
/// Mitochondrial position interval.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MtInterval {
    /// Start position
    pub start: Option<i32>,
    /// End position
    pub end: Option<i32>,
}

impl TryInto<Range<i32>> for MtInterval {
    type Error = Error;

    /// The MT interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn try_into(self) -> Result<Range<i32>, Self::Error> {
        if let (Some(start), Some(end)) = (self.start, self.end) {
            Ok((start - 1)..end)
        } else {
            Err(Error::CannotNonePositionIntoRange(format!("{:?}", self)))
        }
    }
}

/// Transcript sequence location with edit.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TxLocEdit {
    /// Loction on a transcript.
    pub loc: Mu<TxInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl TxLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: Sequence) -> Self {
        TxLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }

    /// Return the LocEdit and ensure that the `NaEdit` has a count and no reference bases.
    fn with_num(self) -> Self {
        TxLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_num()),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_num()),
            },
        }
    }
}

/// Transcript position interval.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TxInterval {
    /// Start position
    pub start: TxPos,
    /// End position
    pub end: TxPos,
}

impl From<TxInterval> for Range<i32> {
    /// The transcript interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn from(val: TxInterval) -> Self {
        if val.start.offset.is_some() || val.end.offset.is_some() {
            warn!("Converting interval {:?} with offset to range!", &val);
        }
        if val.start.base > 0 {
            (val.start.base - 1)..val.end.base
        } else {
            val.start.base..(val.end.base + 1)
        }
    }
}

/// Transcript position.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TxPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
}

/// RNA sequence location with edit.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RnaLocEdit {
    /// Location on a transcript.
    pub loc: Mu<RnaInterval>,
    /// RNA change description.
    pub edit: Mu<NaEdit>,
}

impl RnaLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: Sequence) -> Self {
        RnaLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }

    /// Return the LocEdit and ensure that the `NaEdit` has a count and no reference bases.
    fn with_num(self) -> Self {
        RnaLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_num()),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_num()),
            },
        }
    }
}
/// RNA position interval.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RnaInterval {
    /// Start position
    pub start: RnaPos,
    /// End position
    pub end: RnaPos,
}

impl From<RnaInterval> for Range<i32> {
    /// The RNA interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn from(val: RnaInterval) -> Self {
        if val.start.offset.is_some() || val.end.offset.is_some() {
            warn!("Converting interval {:?} with offset to range!", &val);
        }
        if val.start.base > 0 {
            (val.start.base - 1)..val.end.base
        } else {
            val.start.base..(val.end.base + 1)
        }
    }
}

/// RNA position.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RnaPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
}

/// Protein sequence location with edit or special.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum ProtLocEdit {
    Ordinary {
        loc: Mu<ProtInterval>,
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
    /// `?`
    Unknown,
    /// `Met1?`
    InitiationUncertain,
}

/// Protein position interval.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ProtInterval {
    /// Start position
    pub start: ProtPos,
    /// End position
    pub end: ProtPos,
}

impl From<ProtInterval> for Range<i32> {
    fn from(val: ProtInterval) -> Self {
        if val.start.number > 0 {
            (val.start.number - 1)..val.end.number
        } else {
            val.start.number..(val.end.number + 1)
        }
    }
}

/// Protein position.
#[derive(Clone, Debug, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct ProtPos {
    /// Amino acid value.
    pub aa: u8,
    /// Number of `aa`.
    pub number: i32,
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use super::{TxInterval, TxPos};
    use crate::parser::Mu;

    #[test]
    fn mu_construct() {
        assert_eq!(format!("{:?}", Mu::Certain(1)), "Certain(1)");
        assert_eq!(format!("{:?}", Mu::Certain(Some(1))), "Certain(Some(1))");
        assert_eq!(format!("{:?}", Mu::Uncertain(1)), "Uncertain(1)");
        assert_eq!(
            format!("{:?}", Mu::Uncertain(Some(1))),
            "Uncertain(Some(1))"
        );
    }

    #[test]
    fn mu_unwrap() {
        assert_eq!(Mu::Certain(1).unwrap(), 1);
        assert_eq!(Mu::Uncertain(1).unwrap(), 1);
    }

    #[test]
    fn mu_inner() {
        assert_eq!(Mu::Certain(1).inner(), &1);
        assert_eq!(Mu::Uncertain(1).inner(), &1);
    }

    #[test]
    fn mu_from() {
        assert_eq!(Mu::from(1, true), Mu::Certain(1));
        assert_eq!(Mu::from(1, false), Mu::Uncertain(1));
        assert_eq!(Mu::from(Some(1), true), Mu::Certain(Some(1)));
        assert_eq!(Mu::from(Some(1), false), Mu::Uncertain(Some(1)));
    }

    #[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
    pub struct TestInterval {
        pub start: TestPos,
        pub end: TestPos,
    }

    #[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
    pub struct TestPos {
        pub base: i32,
        pub offset: Option<i32>,
    }

    #[test]
    fn mu_from_problematic() {
        let test_itv = TestInterval {
            start: TestPos {
                base: 34,
                offset: Some(1),
            },
            end: TestPos {
                base: 35,
                offset: Some(-1),
            },
        };

        assert_eq!(
            Mu::from(test_itv, true),
            Mu::Certain(TestInterval {
                start: TestPos {
                    base: 34,
                    offset: Some(1),
                },
                end: TestPos {
                    base: 35,
                    offset: Some(-1),
                }
            })
        );
    }

    #[test]
    fn mu_from_tx_interval() {
        let test_itv = TxInterval {
            start: TxPos {
                base: 34,
                offset: Some(1),
            },
            end: TxPos {
                base: 35,
                offset: Some(-1),
            },
        };

        assert_eq!(
            Mu::from(test_itv, true),
            Mu::Certain(TxInterval {
                start: TxPos {
                    base: 34,
                    offset: Some(1),
                },
                end: TxPos {
                    base: 35,
                    offset: Some(-1),
                }
            })
        );
    }
}

// <LICENSE>
// Copyright 2023 hgvs-rs Contributors
// Copyright 2014 Bioutils Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// </LICENSE>
