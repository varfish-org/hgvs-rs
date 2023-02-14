use nom::{
    branch::alt,
    character::complete::alphanumeric1,
    character::complete::char,
    combinator::{map, opt},
    sequence::{pair, tuple},
    IResult,
};

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

// impl HgvsVariant {
//     pub fn parse(input: &str) -> IResult<&str, Self> {

//     }
// }

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

mod protein {
    use nom::{
        bytes::complete::{take, take_while, take_while1},
        multi::{many0, many1},
    };

    pub static AA1: &str = "ACDEFGHIKLMNPQRSTVWYBZXU";

    pub fn aa1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        let (rest, c) = take(1usize)(input)?;
        if !AA1.contains(&c) {
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Fail,
            )))
        } else {
            Ok((rest, c))
        }
    }

    pub fn aa10(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many0(aa1)(input)
    }

    pub fn aa11(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many1(aa1)(input)
    }

    pub static AAT1: &str = "ACDEFGHIKLMNPQRSTVWYBZXU*";

    pub fn aat1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        let (rest, c) = take(1usize)(input)?;
        if !AAT1.contains(&c) {
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Fail,
            )))
        } else {
            Ok((rest, c))
        }
    }

    pub fn aat10(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many0(aat1)(input)
    }

    pub fn aat11(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many1(aat1)(input)
    }

    pub const AA3: &'static [&'static str] = &[
        "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro",
        "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "Asx", "Glx", "Xaa", "Sec",
    ];

    pub fn aa3(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        let (rest, triplet) = take(3usize)(input)?;
        if !AA3.contains(&triplet) {
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Fail,
            )))
        } else {
            Ok((rest, triplet))
        }
    }

    pub fn aa30(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many0(aa3)(input)
    }

    pub fn aa31(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many1(aa3)(input)
    }

    pub const AAT3: &'static [&'static str] = &[
        "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro",
        "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "Asx", "Glx", "Xaa", "Sec", "Ter",
    ];

    pub fn aat3(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        let (rest, triplet) = take(3usize)(input)?;
        if !AAT3.contains(&triplet) {
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Fail,
            )))
        } else {
            Ok((rest, triplet))
        }
    }

    pub fn aat30(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many0(aat3)(input)
    }

    pub fn aat31(input: &str) -> Result<(&str, Vec<&str>), nom::Err<nom::error::Error<&str>>> {
        many1(aat3)(input)
    }
}

mod protein_edit {
    use nom::branch::alt;
    use nom::character::complete::{digit0, digit1};
    use nom::combinator::opt;
    use nom::multi::many0;
    use nom::sequence::pair;
    use nom::{bytes::complete::tag, IResult};

    use super::protein::*;
    use super::ProteinEdit;

    pub fn subst_qm(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, _) = tag("?")(input)?;
        Ok((
            rest,
            ProteinEdit::Subst {
                alternative: "?".to_owned(),
            },
        ))
    }

    pub fn subst_aa(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, seq) = alt((aat3, aat1))(input)?;
        Ok((
            rest,
            ProteinEdit::Subst {
                alternative: seq.to_owned(),
            },
        ))
    }

    pub fn delins(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, (_, seq)) = pair(tag("delins"), alt((aat31, aat11)))(input)?;
        Ok((
            rest,
            ProteinEdit::DelIns {
                alternative: seq.join(""),
            },
        ))
    }

    pub fn del(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, _) = tag("del")(input)?;
        Ok((rest, ProteinEdit::Del))
    }

    pub fn ins(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, (_, seq)) = pair(tag("ins"), alt((aat31, aat11)))(input)?;
        Ok((
            rest,
            ProteinEdit::Ins {
                alternative: seq.join(""),
            },
        ))
    }

    pub fn dup(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, _) = tag("dup")(input)?;
        Ok((rest, ProteinEdit::Dup))
    }

    pub fn fs(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, alternative) = opt(alt((aat3, aat1)))(input)?;
        let (rest, (_, terminal)) =
            pair(tag("fs"), opt(alt((tag("Ter"), tag("X"), tag("*")))))(rest)?;

        if let Some(terminal) = terminal {
            let (rest, count) = digit0(rest)?;
            if count.is_empty() {
                let (rest, qm) = opt(tag("?"))(rest)?;
                Ok((
                    rest,
                    ProteinEdit::Fs {
                        alternative: alternative.map(str::to_owned),
                        terminal: Some(terminal.to_string()),
                        length: if qm.is_some() {
                            super::UncertainChange::Unknown
                        } else {
                            super::UncertainChange::None
                        },
                    },
                ))
            } else {
                Ok((
                    rest,
                    ProteinEdit::Fs {
                        alternative: alternative.map(str::to_owned),
                        terminal: Some(terminal.to_string()),
                        length: super::UncertainChange::Known(str::parse::<i32>(count).unwrap()),
                    },
                ))
            }
        } else {
            Ok((
                rest,
                ProteinEdit::Fs {
                    alternative: alternative.map(str::to_owned),
                    terminal: None,
                    length: super::UncertainChange::None,
                },
            ))
        }
    }

    pub fn ext_neg_shift(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, aa_ext) = opt(alt((aat3, aat1)))(input)?;
        let (rest, _) = tag("ext")(rest)?;
        let (rest, ext_aa) = opt(alt((aat3, aat1)))(rest)?;
        let (rest, (_, offset)) = pair(tag("-"), digit1)(rest)?;

        Ok((
            rest,
            ProteinEdit::Ext {
                aa_ext: aa_ext.map(str::to_owned),
                ext_aa: ext_aa.map(str::to_owned),
                change: super::UncertainChange::Known(-str::parse::<i32>(offset).unwrap()),
            },
        ))
    }

    pub fn ext_pos_shift(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, aa_ext) = opt(alt((aat3, aat1)))(input)?;
        let (rest, _) = tag("ext")(rest)?;
        let (rest, ext_aa) = alt((tag("Ter"), tag("X"), tag("*")))(rest)?;
        let (rest, offset) = opt(alt((tag("?"), digit1)))(rest)?;

        if let Some(offset) = offset {
            Ok((
                rest,
                ProteinEdit::Ext {
                    aa_ext: aa_ext.map(str::to_owned),
                    ext_aa: Some(ext_aa.to_string()),
                    change: if offset == "?" {
                        super::UncertainChange::Unknown
                    } else {
                        super::UncertainChange::Known(str::parse::<i32>(offset).unwrap())
                    },
                },
            ))
        } else {
            Ok((
                rest,
                ProteinEdit::Ext {
                    aa_ext: aa_ext.map(str::to_owned),
                    ext_aa: Some(ext_aa.to_string()),
                    change: super::UncertainChange::None,
                },
            ))
        }
    }

    pub fn ext_minimal(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, aa_ext) = opt(alt((aat3, aat1)))(input)?;
        let (rest, _) = tag("ext")(rest)?;

        Ok((
            rest,
            ProteinEdit::Ext {
                aa_ext: aa_ext.map(str::to_owned),
                ext_aa: None,
                change: super::UncertainChange::None,
            },
        ))
    }

    pub fn ident(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, _) = tag("=")(input)?;
        Ok((rest, ProteinEdit::Ident))
    }
}

impl ProteinEdit {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        Ok(alt((
            protein_edit::subst_qm,
            protein_edit::subst_aa,
            protein_edit::delins,
            protein_edit::del,
            protein_edit::ins,
            protein_edit::dup,
            protein_edit::fs,
            protein_edit::ext_neg_shift,
            protein_edit::ext_pos_shift,
            protein_edit::ext_minimal,
            protein_edit::ident,
        ))(input)?)
    }
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
    /// Whether the position is uncertain.
    pub uncertain: bool,
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

// cf. https://stackoverflow.com/a/73437782/84349
fn narrowed_alphanumeric1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
    alphanumeric1(input)
}

mod na {
    use nom::{
        bytes::complete::{take_while, take_while1},
        character::complete::one_of,
    };

    pub static NA_IUPAC: &str = "ACGTURYMKWSBDHVNacgturymkwsbdhvn";

    pub fn na(input: &str) -> Result<(&str, char), nom::Err<nom::error::Error<&str>>> {
        one_of(NA_IUPAC)(input)
    }

    pub fn na0(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        take_while(|c: char| NA_IUPAC.contains(c))(input)
    }

    pub fn na1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        take_while1(|c: char| NA_IUPAC.contains(c))(input)
    }
}

mod na_edit {
    use nom::bytes::complete::tag;
    use nom::character::complete::{char as nom_char, digit1};
    use nom::sequence::tuple;
    use nom::{multi::many0, sequence::pair, IResult};

    use super::{
        na::{na, na0, na1},
        NaEdit,
    };

    pub fn ident(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (dna_vec, _)) = pair(many0(na), nom_char('='))(input)?;
        let dna: String = dna_vec.into_iter().collect();
        Ok((
            rest,
            NaEdit::RefAlt {
                reference: dna.clone(),
                alternative: dna,
            },
        ))
    }

    pub fn subst(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (src, _, dst)) = tuple((na, nom_char('>'), na))(input)?;
        Ok((
            rest,
            NaEdit::RefAlt {
                reference: src.to_string(),
                alternative: dst.to_string(),
            },
        ))
    }

    pub fn delins_ref_alt(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (_, reference, _, alternative)) =
            tuple((tag("del"), na0, tag("ins"), na1))(input)?;
        Ok((
            rest,
            NaEdit::RefAlt {
                reference: reference.to_string(),
                alternative: alternative.to_string(),
            },
        ))
    }

    pub fn delins_num_alt(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (_, count, _, alternative)) =
            tuple((tag("del"), digit1, tag("ins"), na1))(input)?;
        Ok((
            rest,
            NaEdit::NumAlt {
                count: count.parse::<u32>().unwrap(),
                alternative: alternative.to_string(),
            },
        ))
    }

    pub fn ins(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (_, alternative)) = tuple((tag("ins"), na1))(input)?;
        Ok((
            rest,
            NaEdit::Ins {
                alternative: alternative.to_string(),
            },
        ))
    }

    pub fn dup(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (_, reference)) = tuple((tag("dup"), na0))(input)?;
        Ok((
            rest,
            NaEdit::Dup {
                reference: reference.to_string(),
            },
        ))
    }

    pub fn inv_num(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (_, count)) = tuple((tag("inv"), digit1))(input)?;
        Ok((
            rest,
            NaEdit::InvNum {
                count: count.parse::<u32>().unwrap(),
            },
        ))
    }
    pub fn inv_ref(input: &str) -> IResult<&str, NaEdit> {
        let (rest, (_, reference)) = tuple((tag("inv"), na0))(input)?;
        Ok((
            rest,
            NaEdit::InvRef {
                reference: reference.to_string(),
            },
        ))
    }
}

impl NaEdit {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        Ok(alt((
            na_edit::ident,
            na_edit::subst,
            na_edit::delins_ref_alt,
            na_edit::delins_num_alt,
            na_edit::ins,
            na_edit::dup,
            na_edit::inv_num,
            na_edit::inv_ref,
        ))(input)?)
    }
}

impl Accession {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        let parser_accession = tuple((
            narrowed_alphanumeric1,
            opt(pair(char('_'), alphanumeric1)),
            opt(pair(char('.'), alphanumeric1)),
        ));

        let mut parser = map(parser_accession, |(a, b, c)| {
            let mut value = a.to_string();
            if let Some(b) = b {
                value.push(b.0);
                value.push_str(b.1);
            }
            if let Some(c) = c {
                value.push(c.0);
                value.push_str(c.1);
            }

            Self { value }
        });

        parser(input)
    }
}

/// Representation of gene symbol, e.g., `TTN` or `Ttn`.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneSymbol {
    pub value: String,
}

impl GeneSymbol {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        let mut gene_symbol_parser = map(alphanumeric1, |symbol: &str| Self {
            value: symbol.to_owned(),
        });
        gene_symbol_parser(input)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn protein_aa1() {
        assert_eq!(protein::aa1("ACD"), Ok(("CD", "A")));
        assert_eq!(
            protein::aa1("*"),
            Err(nom::Err::Error(nom::error::Error::new(
                "*",
                nom::error::ErrorKind::Fail,
            )))
        );
    }

    #[test]
    fn protein_aa10() {
        assert_eq!(protein::aa10(""), Ok(("", vec![])));
        assert_eq!(protein::aa10("ACD"), Ok(("", vec!["A", "C", "D"])));
    }

    #[test]
    fn protein_aa11() {
        assert_eq!(
            protein::aa11(""),
            Err(nom::Err::Error(nom::error::Error::new(
                "",
                nom::error::ErrorKind::Eof,
            )))
        );
        assert_eq!(protein::aa11("ACD"), Ok(("", vec!["A", "C", "D"])));
    }

    #[test]
    fn protein_aat1() {
        assert_eq!(protein::aat1("*ACD"), Ok(("ACD", "*")));
        assert_eq!(
            protein::aat1("="),
            Err(nom::Err::Error(nom::error::Error::new(
                "=",
                nom::error::ErrorKind::Fail,
            )))
        );
    }

    #[test]
    fn protein_aat10() {
        assert_eq!(protein::aat10(""), Ok(("", vec![])));
        assert_eq!(protein::aat10("ACD*"), Ok(("", vec!["A", "C", "D", "*"])));
    }

    #[test]
    fn protein_aat11() {
        assert_eq!(
            protein::aat11(""),
            Err(nom::Err::Error(nom::error::Error::new(
                "",
                nom::error::ErrorKind::Eof,
            )))
        );
        assert_eq!(protein::aat11("ACD*"), Ok(("", vec!["A", "C", "D", "*"])));
    }

    #[test]
    fn protein_aa3() {
        assert_eq!(protein::aa3("LeuMet"), Ok(("Met", "Leu")));
        assert_eq!(
            protein::aa3("Ter"),
            Err(nom::Err::Error(nom::error::Error::new(
                "Ter",
                nom::error::ErrorKind::Fail,
            )))
        );
    }

    #[test]
    fn protein_aa30() {
        assert_eq!(protein::aa30(""), Ok(("", vec![])));
        assert_eq!(
            protein::aa30("MetValTrp"),
            Ok(("", vec!["Met", "Val", "Trp"]))
        );
    }

    #[test]
    fn protein_aa31() {
        assert_eq!(
            protein::aa31(""),
            Err(nom::Err::Error(nom::error::Error::new(
                "",
                nom::error::ErrorKind::Eof,
            )))
        );
        assert_eq!(
            protein::aa30("MetValTrp"),
            Ok(("", vec!["Met", "Val", "Trp"]))
        );
    }

    #[test]
    fn protein_aat3() {
        assert_eq!(protein::aat3("TerLeuMet"), Ok(("LeuMet", "Ter")));
        assert_eq!(
            protein::aat3("==="),
            Err(nom::Err::Error(nom::error::Error::new(
                "===",
                nom::error::ErrorKind::Fail,
            )))
        );
    }

    #[test]
    fn protein_aat30() {
        assert_eq!(protein::aat30(""), Ok(("", vec![])));
        assert_eq!(
            protein::aat30("TerLeuMet"),
            Ok(("", vec!["Ter", "Leu", "Met"]))
        );
    }

    #[test]
    fn protein_aat31() {
        assert_eq!(
            protein::aat31(""),
            Err(nom::Err::Error(nom::error::Error::new(
                "",
                nom::error::ErrorKind::Eof,
            )))
        );
        assert_eq!(
            protein::aat31("TerLeuMet"),
            Ok(("", vec!["Ter", "Leu", "Met"]))
        );
    }

    #[test]
    fn proteinedit_subst_qm() {
        assert_eq!(
            protein_edit::subst_qm("?"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "?".to_owned()
                }
            ))
        );
    }

    #[test]
    fn proteinedit_subst_aa_aat1() {
        assert_eq!(
            protein_edit::subst_aa("L"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "L".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::subst_aa("X"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "X".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::subst_aa("*"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "*".to_owned()
                }
            ))
        );
    }

    #[test]
    fn proteinedit_subst_aa_aat3() {
        assert_eq!(
            protein_edit::subst_aa("Leu"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "Leu".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::subst_aa("Ter"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "Ter".to_owned()
                }
            ))
        );
    }

    #[test]
    fn proteinedit_delins_aat3() {
        assert_eq!(
            protein_edit::delins("delinsLeu"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "Leu".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::delins("delinsLeuTer"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "LeuTer".to_owned()
                }
            ))
        );
    }

    #[test]
    fn proteinedit_delins_aat1() {
        assert_eq!(
            protein_edit::delins("delinsL"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "L".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::delins("delinsL*"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "L*".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::delins("delinsLX"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "LX".to_owned()
                }
            ))
        );
    }

    #[test]
    fn proteinedit_dup() {
        assert_eq!(protein_edit::dup("dup"), Ok(("", ProteinEdit::Dup)));
    }

    #[test]
    fn proteinedit_fs_aat1() {
        assert_eq!(
            protein_edit::fs("Lfs"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: Some("L".to_owned()),
                    terminal: None,
                    length: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::fs("fs"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: None,
                    length: UncertainChange::None
                }
            ))
        );
        assert_eq!(
            protein_edit::fs("fs*"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("*".to_owned()),
                    length: UncertainChange::None
                }
            ))
        );
        assert_eq!(
            protein_edit::fs("fsX"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("X".to_owned()),
                    length: UncertainChange::None
                }
            ))
        );
        assert_eq!(
            protein_edit::fs("fs*?"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("*".to_owned()),
                    length: UncertainChange::Unknown
                }
            ))
        );
        assert_eq!(
            protein_edit::fs("fsX?"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("X".to_owned()),
                    length: UncertainChange::Unknown
                }
            ))
        );
        assert_eq!(
            protein_edit::fs("fs*12"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("*".to_owned()),
                    length: UncertainChange::Known(12)
                }
            ))
        );
        assert_eq!(
            protein_edit::fs("fsX12"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("X".to_owned()),
                    length: UncertainChange::Known(12)
                }
            ))
        );
    }

    #[test]
    fn proteinedit_ext_neg_shift() {
        assert_eq!(
            protein_edit::ext_neg_shift("ext-1"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: None,
                    change: UncertainChange::Known(-1),
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_neg_shift("LextM-1"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("M".to_string()),
                    change: UncertainChange::Known(-1),
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_neg_shift("LeuextMet-1"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Met".to_string()),
                    change: UncertainChange::Known(-1),
                }
            ))
        );
    }

    #[test]
    fn proteinedit_ext_pos_shift() {
        assert_eq!(
            protein_edit::ext_pos_shift("LextX"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("X".to_string()),
                    change: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("Lext*"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("*".to_string()),
                    change: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("LeuextTer"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("extX"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("X".to_string()),
                    change: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("ext*"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("*".to_string()),
                    change: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("extTer"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainChange::None,
                }
            ))
        );

        assert_eq!(
            protein_edit::ext_pos_shift("LextX?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("X".to_string()),
                    change: UncertainChange::Unknown,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("Lext*?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("*".to_string()),
                    change: UncertainChange::Unknown,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("LeuextTer?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainChange::Unknown,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("extX?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("X".to_string()),
                    change: UncertainChange::Unknown,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("ext*?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("*".to_string()),
                    change: UncertainChange::Unknown,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("extTer?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainChange::Unknown,
                }
            ))
        );

        assert_eq!(
            protein_edit::ext_pos_shift("LextX10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("X".to_string()),
                    change: UncertainChange::Known(10),
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("Lext*10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("*".to_string()),
                    change: UncertainChange::Known(10),
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("LeuextTer10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainChange::Known(10),
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("extX10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("X".to_string()),
                    change: UncertainChange::Known(10),
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("ext*10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("*".to_string()),
                    change: UncertainChange::Known(10),
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_pos_shift("extTer10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainChange::Known(10),
                }
            ))
        );
    }

    #[test]
    fn proteinedit_ext_minimal() {
        assert_eq!(
            protein_edit::ext_minimal("Lext"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: None,
                    change: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_minimal("Leuext"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: None,
                    change: UncertainChange::None,
                }
            ))
        );
        assert_eq!(
            protein_edit::ext_minimal("ext"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: None,
                    change: UncertainChange::None,
                }
            ))
        );
    }

    #[test]
    fn proteinedit_ident() {
        assert_eq!(protein_edit::ident("="), Ok(("", ProteinEdit::Ident)));
    }
    #[test]
    fn proteinedit_del() {
        assert_eq!(protein_edit::del("del"), Ok(("", ProteinEdit::Del)));
    }

    #[test]
    fn proteinedit_ins_at1() {
        assert_eq!(
            protein_edit::ins("insL"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "L".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::ins("insL*"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "L*".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::ins("insLX"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "LX".to_owned()
                }
            ))
        );
    }

    #[test]
    fn proteinedit_ins_at3() {
        assert_eq!(
            protein_edit::ins("insLeu"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "Leu".to_owned()
                }
            ))
        );
        assert_eq!(
            protein_edit::ins("insLeuTer"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "LeuTer".to_owned()
                }
            ))
        );
    }

    #[test]
    fn naedit_ident() {
        assert_eq!(
            na_edit::ident("="),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "".to_owned(),
                    alternative: "".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::ident("C="),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "C".to_owned(),
                    alternative: "C".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::ident("CG="),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "CG".to_owned(),
                    alternative: "CG".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_subst() {
        assert_eq!(
            na_edit::subst("C>T"),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "C".to_owned(),
                    alternative: "T".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_delins_ref_alt() {
        assert_eq!(
            na_edit::delins_ref_alt("delinsT"),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "".to_owned(),
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::delins_ref_alt("delGinsT"),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "G".to_owned(),
                    alternative: "T".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_delins_num_alt() {
        assert_eq!(
            na_edit::delins_num_alt("del1insT"),
            Ok((
                "",
                NaEdit::NumAlt {
                    count: 1,
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::delins_num_alt("del1insTT"),
            Ok((
                "",
                NaEdit::NumAlt {
                    count: 1,
                    alternative: "TT".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_ins() {
        assert_eq!(
            na_edit::ins("insT"),
            Ok((
                "",
                NaEdit::Ins {
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::ins("insTT"),
            Ok((
                "",
                NaEdit::Ins {
                    alternative: "TT".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_dup() {
        assert_eq!(
            na_edit::dup("dup"),
            Ok((
                "",
                NaEdit::Dup {
                    reference: "".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::dup("dupT"),
            Ok((
                "",
                NaEdit::Dup {
                    reference: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::dup("dupTT"),
            Ok((
                "",
                NaEdit::Dup {
                    reference: "TT".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_inv_ref() {
        assert_eq!(
            na_edit::inv_ref("inv"),
            Ok((
                "",
                NaEdit::InvRef {
                    reference: "".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::inv_ref("invC"),
            Ok((
                "",
                NaEdit::InvRef {
                    reference: "C".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::inv_ref("invCT"),
            Ok((
                "",
                NaEdit::InvRef {
                    reference: "CT".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_inv_num() {
        assert_eq!(
            na_edit::inv_num("inv1"),
            Ok(("", NaEdit::InvNum { count: 1 }))
        );
        assert_eq!(
            na_edit::inv_num("inv12"),
            Ok(("", NaEdit::InvNum { count: 12 }))
        );
    }

    #[test]
    fn na_edit_parse() {
        assert_eq!(
            NaEdit::parse("="),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "".to_owned(),
                    alternative: "".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("C="),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "C".to_owned(),
                    alternative: "C".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("CG="),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "CG".to_owned(),
                    alternative: "CG".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("C>T"),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "C".to_owned(),
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("delinsT"),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "".to_owned(),
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("delGinsT"),
            Ok((
                "",
                NaEdit::RefAlt {
                    reference: "G".to_owned(),
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("del1insT"),
            Ok((
                "",
                NaEdit::NumAlt {
                    count: 1,
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("del1insTT"),
            Ok((
                "",
                NaEdit::NumAlt {
                    count: 1,
                    alternative: "TT".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("insT"),
            Ok((
                "",
                NaEdit::Ins {
                    alternative: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("insTT"),
            Ok((
                "",
                NaEdit::Ins {
                    alternative: "TT".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("dup"),
            Ok((
                "",
                NaEdit::Dup {
                    reference: "".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("dupT"),
            Ok((
                "",
                NaEdit::Dup {
                    reference: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("dupTT"),
            Ok((
                "",
                NaEdit::Dup {
                    reference: "TT".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("inv"),
            Ok((
                "",
                NaEdit::InvRef {
                    reference: "".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("invC"),
            Ok((
                "",
                NaEdit::InvRef {
                    reference: "C".to_owned(),
                }
            ))
        );
        assert_eq!(
            NaEdit::parse("invCT"),
            Ok((
                "",
                NaEdit::InvRef {
                    reference: "CT".to_owned(),
                }
            ))
        );
        assert_eq!(NaEdit::parse("inv1"), Ok(("", NaEdit::InvNum { count: 1 })));
        assert_eq!(
            NaEdit::parse("inv12"),
            Ok(("", NaEdit::InvNum { count: 12 }))
        );
    }

    #[test]
    fn accession_parse() {
        assert_eq!(
            Accession::parse("NM_01234"),
            Ok((
                "",
                Accession {
                    value: "NM_01234".to_owned()
                }
            )),
        );
        assert_eq!(
            Accession::parse("NM_01234.5"),
            Ok((
                "",
                Accession {
                    value: "NM_01234.5".to_owned()
                }
            )),
        );
        assert_eq!(
            Accession::parse("LRG_01234"),
            Ok((
                "",
                Accession {
                    value: "LRG_01234".to_owned()
                }
            )),
        );
        assert_eq!(
            Accession::parse("LRG_01234.1"),
            Ok((
                "",
                Accession {
                    value: "LRG_01234.1".to_owned()
                }
            )),
        );
    }

    #[test]
    fn gene_symbol_parse() {
        assert_eq!(
            GeneSymbol::parse("TTN"),
            Ok((
                "",
                GeneSymbol {
                    value: "TTN".to_owned()
                }
            ))
        );
        assert_eq!(
            GeneSymbol::parse("Ttn"),
            Ok((
                "",
                GeneSymbol {
                    value: "Ttn".to_owned()
                }
            ))
        );
    }
}
