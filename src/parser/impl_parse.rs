//! Provide implementation of parsing to data structures.

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::char,
    character::complete::{alphanumeric1, digit1, satisfy},
    combinator::{all_consuming, map, opt, recognize},
    sequence::{pair, tuple},
    AsChar, IResult,
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
        all_consuming(alt((
            Self::parse_cds_variant,
            Self::parse_genome_variant,
            Self::parse_mt_variant,
            Self::parse_tx_variant,
            Self::parse_prot_variant,
            Self::parse_rna_variant,
        )))(input)
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
            na_edit::delins_ref_alt,
            na_edit::delins_num_alt,
            na_edit::del_num,
            na_edit::del_ref,
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
            satisfy(|c| c.is_alpha()),
            alphanum::narrowed_alphanumeric1,
            opt(pair(char('_'), alphanumeric1)),
            opt(pair(char('.'), digit1)),
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
    fn parse_initiation_uncertain(input: &str) -> IResult<&str, Self> {
        map(tag("Met1?"), |_| ProtLocEdit::InitiationUncertain)(input)
    }

    fn parse_no_change(input: &str) -> IResult<&str, Self> {
        map(tag("="), |_| ProtLocEdit::NoChange)(input)
    }

    fn parse_no_change_uncertain(input: &str) -> IResult<&str, Self> {
        map(tag("(=)"), |_| ProtLocEdit::NoChangeUncertain)(input)
    }

    fn parse_unknown(input: &str) -> IResult<&str, Self> {
        map(tag("?"), |_| ProtLocEdit::Unknown)(input)
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
            Self::parse_initiation_uncertain,
            Self::parse_ordinary,
            Self::parse_no_protein_uncertain,
            Self::parse_no_protein,
            Self::parse_no_change_uncertain,
            Self::parse_no_change,
            Self::parse_unknown,
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
                            start: CdsPos {
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
                            start: Some(123),
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
                            start: Some(123),
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
                            start: TxPos {
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
                            start: RnaPos {
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
                            start: ProtPos {
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
                        start: CdsPos {
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
                        start: Some(123),
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
                        start: Some(123),
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
                        start: TxPos {
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
                        start: RnaPos {
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
                        start: ProtPos {
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
                    start: CdsPos {
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
                    start: Some(123),
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
                    start: Some(123),
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
                    start: TxPos {
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
                    start: RnaPos {
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
                    start: ProtPos {
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

    // The following is a port of the tests in `test_hgvs_grammar_full.py` of the Python
    // package

    /// Tests of the grammar
    ///
    /// Code takes a tab-delimited text file of the form:
    ///
    ///   Func    Test    Valid   InType  Expected
    ///   pm      -+      True    string
    ///   pm      *       False   one
    ///   num     1|+1    True    list    1|1
    ///
    /// Headers are defined as follows:
    /// Func: function name to call in the grammar
    /// Test: item(s) to test
    /// Valid: if the input is expected to be valid (True or False)
    /// InType: 3 type:
    /// - one: input is a single value
    /// - string: input is a string; test each character in the string separately
    /// - list: input is a list delimited by a pipe character ("|")
    /// Expected: expected result (if stringifying input does not return the same answer, e,g. "+1" -> "1")
    /// - if expected is left blank, then it is assumed that stringifying the parsed input returns the same answer.
    mod grammar_full {
        use nom::combinator::all_consuming;

        use crate::parser::{impl_parse::Parseable, Accession};

        #[test]
        fn parser_grammar() -> Result<(), anyhow::Error> {
            let path = "tests/data/parser/grammar_test.tsv";
            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .flexible(true)
                .from_path(path)?;

            for row in rdr.records() {
                let row = row?;
                if row.get(0) == Some("Func") {
                    continue; // skip header
                }

                // setup input
                let inputs = split_inputs(
                    row.get(1).expect("problem with test input file"),
                    row.get(3).expect("problem with test input file"),
                )?;
                let expected_results = if row.get(4).is_some() && row.get(4) != Some("") {
                    split_inputs(
                        row.get(4).expect("problem with test input file"),
                        row.get(3).expect("problem with test input file"),
                    )?
                } else {
                    inputs.clone()
                };
                let is_valid = row
                    .get(2)
                    .expect("problem with test input file")
                    .to_lowercase()
                    == "true";

                for (input, expected) in inputs.into_iter().zip(expected_results.into_iter()) {
                    let func = row.get(0).expect("problem with test input file");
                    match func {
                        "accn" => {
                            let x = input.clone();
                            let res = all_consuming(Accession::parse)(&x);
                            if is_valid {
                                let (r, acc) = res.expect("problem with test input file");
                                assert_eq!(r, "");
                                assert_eq!(acc.as_str(), expected);
                            } else {
                                assert!(res.is_err());
                            }
                        }
                        _ => println!("We need to implement the other cases as well"),
                    }
                }
            }

            Ok(())
        }

        fn split_inputs(in_string: &str, in_type: &str) -> Result<Vec<String>, anyhow::Error> {
            let inputs = if in_type == "list" {
                in_string
                    .split('|')
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>()
            } else if in_type == "string" {
                in_string.chars().map(|c| c.to_string()).collect::<Vec<_>>()
            } else if in_type == "one" {
                vec![in_string.to_string()]
            } else {
                panic!("should never reach here (in_type={:?})", &in_type)
            };

            Ok(inputs
                .into_iter()
                .filter(|s| s != "None")
                .collect::<Vec<_>>())
        }
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
