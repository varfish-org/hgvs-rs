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

/// Protein edit with interval end edit.
#[derive(Clone, Debug, PartialEq)]
pub enum ProteinEdit {
    Fs {
        alt: Option<String>,
        length: Option<i32>,
    },
    Ext {
        alt: Option<String>,
        length: Option<i32>,
    },
    Subst {
        alt: Option<String>,
    },
    /// `delins`
    DelIns {
        alt: Option<String>,
    },
    /// `ins`
    Ins {
        alt: Option<String>,
    },
    /// `del`
    Del,
    /// `dup`
    Dup,
    /// `=`
    Ident,
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

mod dna {
    use nom::{
        bytes::complete::{take_while, take_while1},
        character::complete::one_of,
    };

    pub static DNA_IUPAC: &str = "ACGTRYMKWSBDHVNacgtrymkwsbdhvn";

    pub fn dna(input: &str) -> Result<(&str, char), nom::Err<nom::error::Error<&str>>> {
        one_of(DNA_IUPAC)(input)
    }

    pub fn dna0(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        take_while(|c: char| DNA_IUPAC.contains(c))(input)
    }

    pub fn dna1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        take_while1(|c: char| DNA_IUPAC.contains(c))(input)
    }
}

mod na_edit {
    use nom::bytes::complete::tag;
    use nom::character::complete::{char as nom_char, digit1};
    use nom::sequence::tuple;
    use nom::{multi::many0, sequence::pair, IResult};

    use super::{
        dna::{dna, dna0, dna1},
        NaEdit,
    };

    pub fn ident(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (dna_vec, _)) = pair(many0(dna), nom_char('='))(input)?;
        let dna: String = dna_vec.into_iter().collect();
        Ok((
            rest,
            NaEdit::RefAlt {
                reference: dna.clone(),
                alternative: dna,
            },
        ))
    }

    pub fn subst(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (src, _, dst)) = tuple((dna, nom_char('>'), dna))(input)?;
        Ok((
            rest,
            NaEdit::RefAlt {
                reference: src.to_string(),
                alternative: dst.to_string(),
            },
        ))
    }

    pub fn delins_ref_alt(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (_, reference, _, alternative)) =
            tuple((tag("del"), dna0, tag("ins"), dna1))(input)?;
        Ok((
            rest,
            NaEdit::RefAlt {
                reference: reference.to_string(),
                alternative: alternative.to_string(),
            },
        ))
    }

    pub fn delins_num_alt(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (_, count, _, alternative)) =
            tuple((tag("del"), digit1, tag("ins"), dna1))(input)?;
        Ok((
            rest,
            NaEdit::NumAlt {
                count: count.parse::<u32>().unwrap(),
                alternative: alternative.to_string(),
            },
        ))
    }

    pub fn ins(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (_, alternative)) = tuple((tag("ins"), dna1))(input)?;
        Ok((
            rest,
            NaEdit::Ins {
                alternative: alternative.to_string(),
            },
        ))
    }

    pub fn dup(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (_, reference)) = tuple((tag("dup"), dna0))(input)?;
        Ok((
            rest,
            NaEdit::Dup {
                reference: reference.to_string(),
            },
        ))
    }

    pub fn inv_num(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (_, count)) = tuple((tag("inv"), digit1))(input)?;
        Ok((
            rest,
            NaEdit::InvNum {
                count: count.parse::<u32>().unwrap(),
            },
        ))
    }
    pub fn inv_ref(input: &str) -> IResult<&str, super::NaEdit> {
        let (rest, (_, reference)) = tuple((tag("inv"), dna0))(input)?;
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

    // #[test]
    // fn parse_protein_edit() {}

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
