//! Data structures for representing HGVS variant descriptions.

use std::ops::{Deref, Range};

use log::warn;

/// Expression of "maybe uncertain".
#[derive(Clone, Debug, PartialEq)]
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
#[derive(Clone, Debug, PartialEq)]
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
    /// Deletion of one or more nucleic acid characters.
    DelRef { reference: String },
    /// Deletion of a number of characters.
    DelNum { count: i32 },
    /// Insertion of one or more nucleic acid characters.
    Ins { alternative: String },
    /// Duplication of nucleic acid reference sequence.
    Dup { reference: String },
    /// Inversion of a (potentially empty) nucleic acid reference sequence.
    InvRef { reference: String },
    /// Inversion of a stretch given by its length.
    InvNum { count: i32 },
}

impl NaEdit {
    /// Return whether the reference equals the given value.
    pub fn reference_equals(&self, value: &str) -> bool {
        match self {
            NaEdit::RefAlt { reference, .. }
            | NaEdit::DelRef { reference }
            | NaEdit::Dup { reference }
            | NaEdit::InvRef { reference } => reference == value,
            _ => false,
        }
    }

    /// Return an updated `NaEdit` that has the reference replaced with the given sequence.
    pub fn with_reference(self, reference: String) -> Self {
        match self {
            NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => {
                // let (_, reference, alternative) = trim_common_suffixes(&reference, &alternative);
                // let (_, reference, alternative) = trim_common_prefixes(&reference, &alternative);
                if reference == alternative {
                    NaEdit::RefAlt {
                        reference: "".to_string(),
                        alternative: "".to_string(),
                    }
                } else {
                    NaEdit::RefAlt {
                        reference,
                        alternative,
                    }
                }
            }
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

impl lazy_static::__Deref for Accession {
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
    // Replace reference sequence.
    pub fn with_reference(self, value: String) -> Self {
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
            HgvsVariant::TxVariant { loc_edit, .. } => loc_edit
                .loc
                .inner()
                .clone()
                .try_into()
                .ok()
                .map(|l| Mu::from(l, loc_edit.loc.is_certain())),
            HgvsVariant::RnaVariant { loc_edit, .. } => loc_edit
                .loc
                .inner()
                .clone()
                .try_into()
                .ok()
                .map(|l| Mu::from(l, loc_edit.loc.is_certain())),
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
                loc_edit.loc.inner().start.offset.unwrap_or_default() > 0
                    || loc_edit.loc.inner().end.offset.unwrap_or_default() > 0
            }
            HgvsVariant::TxVariant { loc_edit, .. } => {
                loc_edit.loc.inner().start.offset.unwrap_or_default() > 0
                    || loc_edit.loc.inner().end.offset.unwrap_or_default() > 0
            }
            HgvsVariant::RnaVariant { loc_edit, .. } => {
                loc_edit.loc.inner().start.offset.unwrap_or_default() > 0
                    || loc_edit.loc.inner().end.offset.unwrap_or_default() > 0
            }
            _ => false,
        }
    }
}

/// Coding sequence location with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct CdsLocEdit {
    /// Location on the CDS.
    pub loc: Mu<CdsInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl CdsLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: String) -> Self {
        CdsLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }
}

/// CDS position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct CdsInterval {
    /// Start position
    pub start: CdsPos,
    /// End position
    pub end: CdsPos,
}

impl TryFrom<CdsInterval> for Range<i32> {
    type Error = anyhow::Error;

    /// The CDS interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn try_from(value: CdsInterval) -> Result<Self, Self::Error> {
        if value.start.offset.is_some() || value.end.offset.is_some() {
            warn!("Converting interval {:?} with offset to range!", &value);
        }
        if value.start.cds_from != value.end.cds_from {
            Err(anyhow::anyhow!(
                "Conversion of interval with different offsets (CDS start/end) is ill-defined"
            ))
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
#[derive(Clone, Debug, PartialEq)]
pub enum CdsFrom {
    Start,
    End,
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

/// Genome sequence location with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeLocEdit {
    /// Location on the genome.
    pub loc: Mu<GenomeInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl GenomeLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: String) -> Self {
        GenomeLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }
}

/// Genome position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeInterval {
    /// Start position
    pub start: Option<i32>,
    /// End position
    pub end: Option<i32>,
}

impl TryInto<Range<i32>> for GenomeInterval {
    type Error = anyhow::Error;

    /// The genome interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn try_into(self) -> Result<Range<i32>, Self::Error> {
        if let (Some(start), Some(end)) = (self.start, self.end) {
            Ok((start - 1)..end)
        } else {
            Err(anyhow::anyhow!(
                "Cannot convert interval with None position into range: {:?}",
                self
            ))
        }
    }
}

/// Mitochondrial sequence location with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct MtLocEdit {
    /// Location on the mitochondrium.
    pub loc: Mu<MtInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl MtLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: String) -> Self {
        MtLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }
}
/// Mitochondrial position interval.
#[derive(Clone, Debug, PartialEq)]
pub struct MtInterval {
    /// Start position
    pub start: Option<i32>,
    /// End position
    pub end: Option<i32>,
}

impl TryInto<Range<i32>> for MtInterval {
    type Error = anyhow::Error;

    /// The MT interval will be converted from 1-based inclusive coordinates
    /// `[start, end]` to 0-based, half-open Rust range `[start - 1, end)`.
    fn try_into(self) -> Result<Range<i32>, Self::Error> {
        if let (Some(start), Some(end)) = (self.start, self.end) {
            Ok((start - 1)..end)
        } else {
            Err(anyhow::anyhow!(
                "Cannot convert interval with None position into range: {:?}",
                self
            ))
        }
    }
}

/// Transcript sequence location with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct TxLocEdit {
    /// Loction on a transcript.
    pub loc: Mu<TxInterval>,
    /// DNA change description.
    pub edit: Mu<NaEdit>,
}

impl TxLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: String) -> Self {
        TxLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }
}

/// Transcript position interval.
#[derive(Clone, Debug, PartialEq)]
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
#[derive(Clone, Debug, PartialEq)]
pub struct TxPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
}

/// RNA sequence location with edit.
#[derive(Clone, Debug, PartialEq)]
pub struct RnaLocEdit {
    /// Location on a transcript.
    pub loc: Mu<RnaInterval>,
    /// RNA change description.
    pub edit: Mu<NaEdit>,
}

impl RnaLocEdit {
    /// Return the LocEdit with the reference replaced by `reference`.
    fn with_reference(self, reference: String) -> Self {
        RnaLocEdit {
            loc: self.loc,
            edit: match self.edit {
                Mu::Certain(edit) => Mu::Certain(edit.with_reference(reference)),
                Mu::Uncertain(edit) => Mu::Uncertain(edit.with_reference(reference)),
            },
        }
    }
}
/// RNA position interval.
#[derive(Clone, Debug, PartialEq)]
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
#[derive(Clone, Debug, PartialEq)]
pub struct RnaPos {
    /// Base position.
    pub base: i32,
    /// Optional offset.
    pub offset: Option<i32>,
}

/// Protein sequence location with edit or special.
#[derive(Clone, Debug, PartialEq)]
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
}

/// Protein position interval.
#[derive(Clone, Debug, PartialEq)]
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
#[derive(Clone, Debug, PartialEq)]
pub struct ProtPos {
    /// Amino acid value.
    pub aa: String,
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

    #[derive(Clone, Debug, PartialEq)]
    pub struct TestInterval {
        pub start: TestPos,
        pub end: TestPos,
    }

    #[derive(Clone, Debug, PartialEq)]
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
