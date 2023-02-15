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

impl HgvsVariant {
    fn parse_cds_variant(input: &str) -> IResult<&str, Self> {
        map(
            tuple((
                Accession::parse,
                opt(tuple((tag("("), GeneSymbol::parse, tag(")")))),
                tag(":c."),
                CdsLocEdit::parse,
            )),
            |(accession, opt_gs, _, pos_edit)| HgvsVariant::CdsVariant {
                accession,
                gene_symbol: opt_gs.map(|(_, gene_symbol, _)| gene_symbol),
                loc_edit: pos_edit,
            },
        )(input)
    }

    fn parse_genome_variant(input: &str) -> IResult<&str, Self> {
        map(
            tuple((
                Accession::parse,
                opt(tuple((tag("("), GeneSymbol::parse, tag(")")))),
                tag(":g."),
                GenomeLocEdit::parse,
            )),
            |(accession, opt_gs, _, pos_edit)| HgvsVariant::GenomeVariant {
                accession,
                gene_symbol: opt_gs.map(|(_, gene_symbol, _)| gene_symbol),
                loc_edit: pos_edit,
            },
        )(input)
    }

    fn parse_mt_variant(input: &str) -> IResult<&str, Self> {
        map(
            tuple((
                Accession::parse,
                opt(tuple((tag("("), GeneSymbol::parse, tag(")")))),
                tag(":m."),
                MtLocEdit::parse,
            )),
            |(accession, opt_gs, _, pos_edit)| HgvsVariant::MtVariant {
                accession,
                gene_symbol: opt_gs.map(|(_, gene_symbol, _)| gene_symbol),
                loc_edit: pos_edit,
            },
        )(input)
    }

    fn parse_tx_variant(input: &str) -> IResult<&str, Self> {
        map(
            tuple((
                Accession::parse,
                opt(tuple((tag("("), GeneSymbol::parse, tag(")")))),
                tag(":n."),
                TxLocEdit::parse,
            )),
            |(accession, opt_gs, _, pos_edit)| HgvsVariant::TxVariant {
                accession,
                gene_symbol: opt_gs.map(|(_, gene_symbol, _)| gene_symbol),
                loc_edit: pos_edit,
            },
        )(input)
    }

    fn parse_prot_variant(input: &str) -> IResult<&str, Self> {
        map(
            tuple((
                Accession::parse,
                opt(tuple((tag("("), GeneSymbol::parse, tag(")")))),
                tag(":p."),
                ProtLocEdit::parse,
            )),
            |(accession, opt_gs, _, pos_edit)| HgvsVariant::ProtVariant {
                accession,
                gene_symbol: opt_gs.map(|(_, gene_symbol, _)| gene_symbol),
                loc_edit: pos_edit,
            },
        )(input)
    }

    fn parse_rna_variant(input: &str) -> IResult<&str, Self> {
        map(
            tuple((
                Accession::parse,
                opt(tuple((tag("("), GeneSymbol::parse, tag(")")))),
                tag(":r."),
                RnaLocEdit::parse,
            )),
            |(accession, opt_gs, _, pos_edit)| HgvsVariant::RnaVariant {
                accession,
                gene_symbol: opt_gs.map(|(_, gene_symbol, _)| gene_symbol),
                loc_edit: pos_edit,
            },
        )(input)
    }
}

impl Parseable for HgvsVariant {
    /// Parse a `HgvsVariant` from the given `str`.
    fn parse(input: &str) -> IResult<&str, Self> {
        alt((
            Self::parse_cds_variant,
            Self::parse_genome_variant,
            Self::parse_mt_variant,
            Self::parse_tx_variant,
            Self::parse_prot_variant,
            Self::parse_rna_variant,
        ))(input)
    }
}

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
            na_edit::del_ref,
            na_edit::del_num,
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
        cds_pos::loc(input)
    }
}

impl Parseable for CdsLocEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<CdsInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| CdsLocEdit { loc: pos, edit },
        )(input)
    }
}

impl Parseable for GenomeInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        genome_pos::loc(input)
    }
}

impl Parseable for GenomeLocEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<GenomeInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| GenomeLocEdit { loc: pos, edit },
        )(input)
    }
}

impl Parseable for MtInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        mt_pos::loc(input)
    }
}

impl Parseable for MtLocEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<MtInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| MtLocEdit { loc: pos, edit },
        )(input)
    }
}

impl Parseable for TxInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        tx_pos::loc(input)
    }
}

impl Parseable for TxLocEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<TxInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| TxLocEdit { loc: pos, edit },
        )(input)
    }
}

impl Parseable for RnaInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        rna_pos::loc(input)
    }
}

impl Parseable for RnaLocEdit {
    fn parse(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<RnaInterval>::parse, Mu::<NaEdit>::parse),
            |(pos, edit)| RnaLocEdit { loc: pos, edit },
        )(input)
    }
}

impl Parseable for ProtInterval {
    fn parse(input: &str) -> IResult<&str, Self> {
        prot_pos::loc(input)
    }
}

impl ProtLocEdit {
    fn parse_no_change(input: &str) -> IResult<&str, Self> {
        map(tag("="), |_| ProtLocEdit::NoChange)(input)
    }

    fn parse_no_change_uncertain(input: &str) -> IResult<&str, Self> {
        map(tag("(=)"), |_| ProtLocEdit::NoChangeUncertain)(input)
    }

    fn parse_no_protein(input: &str) -> IResult<&str, Self> {
        map(tag("0"), |_| ProtLocEdit::NoProtein)(input)
    }

    fn parse_no_protein_uncertain(input: &str) -> IResult<&str, Self> {
        map(tag("0?"), |_| ProtLocEdit::NoProteinUncertain)(input)
    }

    fn parse_ordinary(input: &str) -> IResult<&str, Self> {
        map(
            pair(Mu::<ProtInterval>::parse, Mu::<ProteinEdit>::parse),
            |(pos, edit)| ProtLocEdit::Ordinary { loc: pos, edit },
        )(input)
    }
}

impl Parseable for ProtLocEdit {
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
    fn hgvsvariant_parse() {
        assert_eq!(
            HgvsVariant::parse("NR_01234.1(XYZ):c.123_123C>T"),
            Ok((
                "",
                HgvsVariant::CdsVariant {
                    accession: Accession {
                        value: "NR_01234.1".to_string()
                    },
                    gene_symbol: Some(GeneSymbol {
                        value: "XYZ".to_string()
                    }),
                    loc_edit: CdsLocEdit {
                        loc: Mu::Certain(CdsInterval {
                            begin: CdsPos {
                                base: 123,
                                offset: None,
                                cds_from: CdsFrom::Start
                            },
                            end: CdsPos {
                                base: 123,
                                offset: None,
                                cds_from: CdsFrom::Start
                            },
                        }),
                        edit: Mu::Certain(NaEdit::RefAlt {
                            reference: "C".to_string(),
                            alternative: "T".to_string()
                        })
                    }
                }
            ))
        );
        assert_eq!(
            HgvsVariant::parse("NR_01234.1(XYZ):g.123_123C>T"),
            Ok((
                "",
                HgvsVariant::GenomeVariant {
                    accession: Accession {
                        value: "NR_01234.1".to_string()
                    },
                    gene_symbol: Some(GeneSymbol {
                        value: "XYZ".to_string()
                    }),
                    loc_edit: GenomeLocEdit {
                        loc: Mu::Certain(GenomeInterval {
                            begin: Some(123),
                            end: Some(123),
                        }),
                        edit: Mu::Certain(NaEdit::RefAlt {
                            reference: "C".to_string(),
                            alternative: "T".to_string()
                        })
                    }
                }
            ))
        );
        assert_eq!(
            HgvsVariant::parse("NR_01234.1(XYZ):m.123_123C>T"),
            Ok((
                "",
                HgvsVariant::MtVariant {
                    accession: Accession {
                        value: "NR_01234.1".to_string()
                    },
                    gene_symbol: Some(GeneSymbol {
                        value: "XYZ".to_string()
                    }),
                    loc_edit: MtLocEdit {
                        loc: Mu::Certain(MtInterval {
                            begin: Some(123),
                            end: Some(123),
                        }),
                        edit: Mu::Certain(NaEdit::RefAlt {
                            reference: "C".to_string(),
                            alternative: "T".to_string()
                        })
                    }
                }
            ))
        );
        assert_eq!(
            HgvsVariant::parse("NR_01234.1(XYZ):n.123_123C>T"),
            Ok((
                "",
                HgvsVariant::TxVariant {
                    accession: Accession {
                        value: "NR_01234.1".to_string()
                    },
                    gene_symbol: Some(GeneSymbol {
                        value: "XYZ".to_string()
                    }),
                    loc_edit: TxLocEdit {
                        loc: Mu::Certain(TxInterval {
                            begin: TxPos {
                                base: 123,
                                offset: None
                            },
                            end: TxPos {
                                base: 123,
                                offset: None
                            },
                        }),
                        edit: Mu::Certain(NaEdit::RefAlt {
                            reference: "C".to_string(),
                            alternative: "T".to_string()
                        })
                    }
                }
            ))
        );
        assert_eq!(
            HgvsVariant::parse("NR_01234.1(XYZ):r.123_123C>T"),
            Ok((
                "",
                HgvsVariant::RnaVariant {
                    accession: Accession {
                        value: "NR_01234.1".to_string()
                    },
                    gene_symbol: Some(GeneSymbol {
                        value: "XYZ".to_string()
                    }),
                    loc_edit: RnaLocEdit {
                        loc: Mu::Certain(RnaInterval {
                            begin: RnaPos {
                                base: 123,
                                offset: None
                            },
                            end: RnaPos {
                                base: 123,
                                offset: None
                            },
                        }),
                        edit: Mu::Certain(NaEdit::RefAlt {
                            reference: "C".to_string(),
                            alternative: "T".to_string()
                        })
                    }
                }
            ))
        );
        assert_eq!(
            HgvsVariant::parse("NR_01234.1(XYZ):p.Leu3_Leu3Thr"),
            Ok((
                "",
                HgvsVariant::ProtVariant {
                    accession: Accession {
                        value: "NR_01234.1".to_string()
                    },
                    gene_symbol: Some(GeneSymbol {
                        value: "XYZ".to_string()
                    }),
                    loc_edit: ProtLocEdit::Ordinary {
                        loc: Mu::Certain(ProtInterval {
                            begin: ProtPos {
                                aa: "Leu".to_string(),
                                number: 3
                            },
                            end: ProtPos {
                                aa: "Leu".to_string(),
                                number: 3
                            },
                        }),
                        edit: Mu::Certain(ProteinEdit::Subst {
                            alternative: "Thr".to_string()
                        })
                    }
                }
            ))
        );
    }

    #[test]
    fn cdsposedit_parse() {
        assert_eq!(
            CdsLocEdit::parse("(123_123)C>T"),
            Ok((
                "",
                CdsLocEdit {
                    loc: Mu::Uncertain(CdsInterval {
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
            GenomeLocEdit::parse("(123_123)C>T"),
            Ok((
                "",
                GenomeLocEdit {
                    loc: Mu::Uncertain(GenomeInterval {
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
            MtLocEdit::parse("(123_123)C>T"),
            Ok((
                "",
                MtLocEdit {
                    loc: Mu::Uncertain(MtInterval {
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
            TxLocEdit::parse("(123_123)C>T"),
            Ok((
                "",
                TxLocEdit {
                    loc: Mu::Uncertain(TxInterval {
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
            RnaLocEdit::parse("(123_123)C>T"),
            Ok((
                "",
                RnaLocEdit {
                    loc: Mu::Uncertain(RnaInterval {
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
            ProtLocEdit::parse("(Leu123_Leu123)(Thr)"),
            Ok((
                "",
                ProtLocEdit::Ordinary {
                    loc: Mu::Uncertain(ProtInterval {
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
        assert_eq!(ProtLocEdit::parse("="), Ok(("", ProtLocEdit::NoChange)));
        assert_eq!(
            ProtLocEdit::parse("(=)"),
            Ok(("", ProtLocEdit::NoChangeUncertain))
        );
        assert_eq!(ProtLocEdit::parse("0"), Ok(("", ProtLocEdit::NoProtein)));
        assert_eq!(
            ProtLocEdit::parse("0?"),
            Ok(("", ProtLocEdit::NoProteinUncertain))
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
