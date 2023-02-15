//! Provide implementation of parsing to data structures.

use nom::{
    branch::alt,
    character::complete::alphanumeric1,
    character::complete::char,
    combinator::{map, opt, recognize},
    sequence::{pair, tuple},
    IResult,
};

use crate::parser::ds::*;
use crate::parser::parse_funcs::*;

impl ProteinEdit {
    pub fn parse(input: &str) -> IResult<&str, Self> {
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

impl NaEdit {
    pub fn parse(input: &str) -> IResult<&str, Self> {
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

impl Accession {
    pub fn parse(input: &str) -> IResult<&str, Self> {
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

impl CdsInterval {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        cds_pos::int(input)
    }
}

impl GenomeInterval {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        genome_pos::int(input)
    }
}

impl MtInterval {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        mt_pos::int(input)
    }
}

impl TxInterval {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        tx_pos::int(input)
    }
}

impl RnaInterval {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        rna_pos::int(input)
    }
}

impl ProtInterval {
    pub fn parse(input: &str) -> IResult<&str, Self> {
        prot_pos::int(input)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use pretty_assertions::assert_eq;

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
