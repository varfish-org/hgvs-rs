//! Provide implementation of parsing to data structures.

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::alphanumeric1,
    character::complete::char,
    combinator::{map, opt, recognize},
    sequence::{pair, tuple},
    IResult,
};

use crate::parser::ds::*;
use crate::parser::parse_funcs::*;

/// Parsers for constructing `HgvsVariant` records.
mod hgvs_variant {}

// impl HgvsVariant {
//     fn parse_cds_variant(input: &str) -> IResult<&str, Self> {
//         let xs = tuple((Accession::parse, opt(GeneSymbol::parse), CdsPosEdit::parse))(input)?;
//     }

//     /// Parse a `HgvsVariant` from the given `str`.
//     pub fn parse(input: &str) -> IResult<&str, Self> {

//     }
// }

pub trait Parseable {
    fn parse(input: &str) -> IResult<&str, Self>
    where
        Self: Sized;
}

impl Parseable for ProteinEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        alt((
            protein_edit::fs,
            protein_edit::ext_neg_shift,
            protein_edit::ext_pos_shift,
            protein_edit::ext_minimal,
            protein_edit::ident,
            protein_edit::subst_qm,
            protein_edit::subst_aa,
            protein_edit::delins,
            protein_edit::del,
            protein_edit::ins,
            protein_edit::dup,
        ))(input)
    }
}

impl Parseable for NaEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        alt((
            na_edit::ident,
            na_edit::subst,
            na_edit::delins_ref_alt,
            na_edit::delins_num_alt,
            na_edit::ins,
            na_edit::dup,
            na_edit::inv_num,
            na_edit::inv_ref,
        ))(input)
    }
}

impl Parseable for Accession {
    fn parse(input: &str) -> IResult<&str, Self> {
        let parser_accession = recognize(tuple((
            alphanum::narrowed_alphanumeric1,
            opt(pair(char('_'), alphanumeric1)),
            opt(pair(char('.'), alphanumeric1)),
        )));

        let mut parser = map(parser_accession, |value| Self {
            value: value.to_string(),
        });

        parser(input)
    }
}

impl GeneSymbol {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        let mut gene_symbol_parser = map(alphanumeric1, |symbol: &str| Self {
            value: symbol.to_owned(),
        });
        gene_symbol_parser(input)
    }
}

impl<T> Mu<T>
where
    T: Parseable,
{
    fn parse_uncertain(input: &str) -> IResult<&str, Self>
    where
        Self: Sized,
    {
        map(tuple((tag("("), T::parse, tag(")"))), |(_, value, _)| {
            Mu::Uncertain(value)
        })(input)
    }

    fn parse_certain(input: &str) -> IResult<&str, Self>
    where
        Self: Sized,
    {
        map(T::parse, |value| Mu::Certain(value))(input)
    }
}

impl<T> Parseable for Mu<T>
where
    T: Parseable,
{
    fn parse(input: &str) -> IResult<&str, Self>
    where
        Self: Sized,
    {
        alt((Mu::<T>::parse_uncertain, Mu::<T>::parse_certain))(input)
    }
}

impl Parseable for CdsInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        cds_pos::int(input)
    }
}

impl Parseable for CdsPosEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<CdsInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| CdsPosEdit { pos, edit },
        )(input)
    }
}

impl Parseable for GenomeInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        genome_pos::int(input)
    }
}

impl Parseable for GenomePosEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<GenomeInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| GenomePosEdit { pos, edit },
        )(input)
    }
}

impl Parseable for MtInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        mt_pos::int(input)
    }
}

impl Parseable for MtPosEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<MtInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| MtPosEdit { pos, edit },
        )(input)
    }
}

impl Parseable for TxInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        tx_pos::int(input)
    }
}

impl Parseable for TxPosEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<TxInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| TxPosEdit { pos, edit },
        )(input)
    }
}

impl Parseable for RnaInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        rna_pos::int(input)
    }
}

impl Parseable for RnaPosEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<RnaInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| RnaPosEdit { pos, edit },
        )(input)
    }
}

impl Parseable for ProtInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        prot_pos::int(input)
    }
}

impl ProtPosEdit {
    fn parse_no_change(input: &str) -> IResult<&str, Self> {
        map(tag("="), |_| ProtPosEdit::NoChange)(input)
    }

    fn parse_no_change_uncertain(input: &str) -> IResult<&str, Self> {
        map(tag("(=)"), |_| ProtPosEdit::NoChangeUncertain)(input)
    }

    fn parse_no_protein(input: &str) -> IResult<&str, Self> {
        map(tag("0"), |_| ProtPosEdit::NoProtein)(input)
    }

    fn parse_no_protein_uncertain(input: &str) -> IResult<&str, Self> {
        map(tag("0?"), |_| ProtPosEdit::NoProteinUncertain)(input)
    }

    fn parse_ordinary(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<ProtInterval>::parse, Mu::<ProteinEdit>::parse),
            |(pos, edit)| ProtPosEdit::Ordinary { pos, edit },
        )(input)
    }
}

impl Parseable for ProtPosEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        alt((
            Self::parse_ordinary,
            Self::parse_no_protein_uncertain,
            Self::parse_no_protein,
            Self::parse_no_change_uncertain,
            Self::parse_no_change,
        ))(input)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn cdsposedit_parse() {
        assert_eq!(
            CdsPosEdit::parse("(123_123)C>T"),
            Ok((
                "",
                CdsPosEdit {
                    pos: Mu::Uncertain(CdsInterval {
                        begin: CdsPos {
                            base: 123,
                            offset: None,
                            cds_from: CdsFrom::Start
                        },
                        end: CdsPos {
                            base: 123,
                            offset: None,
                            cds_from: CdsFrom::Start
                        }
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: "C".to_string(),
                        alternative: "T".to_string()
                    })
                }
            ))
        )
    }

    #[test]
    fn genomeposedit_parse() {
        assert_eq!(
            GenomePosEdit::parse("(123_123)C>T"),
            Ok((
                "",
                GenomePosEdit {
                    pos: Mu::Uncertain(GenomeInterval {
                        begin: Some(123),
                        end: Some(123),
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: "C".to_string(),
                        alternative: "T".to_string()
                    })
                }
            ))
        )
    }

    #[test]
    fn mtposedit_parse() {
        assert_eq!(
            MtPosEdit::parse("(123_123)C>T"),
            Ok((
                "",
                MtPosEdit {
                    pos: Mu::Uncertain(MtInterval {
                        begin: Some(123),
                        end: Some(123),
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: "C".to_string(),
                        alternative: "T".to_string()
                    })
                }
            ))
        )
    }

    #[test]
    fn txposedit_parse() {
        assert_eq!(
            TxPosEdit::parse("(123_123)C>T"),
            Ok((
                "",
                TxPosEdit {
                    pos: Mu::Uncertain(TxInterval {
                        begin: TxPos {
                            base: 123,
                            offset: None,
                        },
                        end: TxPos {
                            base: 123,
                            offset: None,
                        }
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: "C".to_string(),
                        alternative: "T".to_string()
                    })
                }
            ))
        )
    }

    #[test]
    fn rnaposedit_parse() {
        assert_eq!(
            RnaPosEdit::parse("(123_123)C>T"),
            Ok((
                "",
                RnaPosEdit {
                    pos: Mu::Uncertain(RnaInterval {
                        begin: RnaPos {
                            base: 123,
                            offset: None,
                        },
                        end: RnaPos {
                            base: 123,
                            offset: None,
                        }
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: "C".to_string(),
                        alternative: "T".to_string()
                    })
                }
            ))
        )
    }

    #[test]
    fn protposedit_parse() {
        assert_eq!(
            ProtPosEdit::parse("(Leu123_Leu123)(Thr)"),
            Ok((
                "",
                ProtPosEdit::Ordinary {
                    pos: Mu::Uncertain(ProtInterval {
                        begin: ProtPos {
                            aa: "Leu".to_string(),
                            number: 123,
                        },
                        end: ProtPos {
                            aa: "Leu".to_string(),
                            number: 123,
                        }
                    }),
                    edit: Mu::Uncertain(ProteinEdit::Subst {
                        alternative: "Thr".to_string()
                    })
                }
            ))
        );
        assert_eq!(ProtPosEdit::parse("="), Ok(("", ProtPosEdit::NoChange)));
        assert_eq!(
            ProtPosEdit::parse("(=)"),
            Ok(("", ProtPosEdit::NoChangeUncertain))
        );
        assert_eq!(ProtPosEdit::parse("0"), Ok(("", ProtPosEdit::NoProtein)));
        assert_eq!(
            ProtPosEdit::parse("0?"),
            Ok(("", ProtPosEdit::NoProteinUncertain))
        );
    }

    #[test]
    fn cdsinterval_parse() {
        assert_eq!(
            CdsInterval::parse("123_123"),
            Ok((
                "",
                CdsInterval {
                    begin: CdsPos {
                        base: 123,
                        offset: None,
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: 123,
                        offset: None,
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
    }

    #[test]
    fn genomeinterval_parse() {
        assert_eq!(
            GenomeInterval::parse("123_123"),
            Ok((
                "",
                GenomeInterval {
                    begin: Some(123),
                    end: Some(123),
                }
            ))
        );
    }

    #[test]
    fn mtinterval_parse() {
        assert_eq!(
            MtInterval::parse("123_123"),
            Ok((
                "",
                MtInterval {
                    begin: Some(123),
                    end: Some(123),
                }
            ))
        );
    }

    #[test]
    fn txinterval_parse() {
        assert_eq!(
            TxInterval::parse("123_123"),
            Ok((
                "",
                TxInterval {
                    begin: TxPos {
                        base: 123,
                        offset: None
                    },
                    end: TxPos {
                        base: 123,
                        offset: None
                    },
                }
            ))
        );
    }

    #[test]
    fn rnainterval_parse() {
        assert_eq!(
            RnaInterval::parse("123_123"),
            Ok((
                "",
                RnaInterval {
                    begin: RnaPos {
                        base: 123,
                        offset: None
                    },
                    end: RnaPos {
                        base: 123,
                        offset: None
                    },
                }
            ))
        );
    }

    #[test]
    fn protinterval_parse() {
        assert_eq!(
            ProtInterval::parse("Leu123_Leu123"),
            Ok((
                "",
                ProtInterval {
                    begin: ProtPos {
                        aa: "Leu".to_string(),
                        number: 123
                    },
                    end: ProtPos {
                        aa: "Leu".to_string(),
                        number: 123
                    },
                }
            ))
        );
    }

    #[test]
    fn proteinedit_parse() {
        assert_eq!(
            ProteinEdit::parse("?"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "?".to_owned()
                }
            ))
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
    }
}
