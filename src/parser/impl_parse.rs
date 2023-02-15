//! Provide implementation of parsing to data structures.

use nom::{
    branch::alt,
    character::complete::alphanumeric1,
    character::complete::char,
    combinator::{map, opt},
    sequence::{pair, tuple},
    IResult,
};

use crate::parser::parse_funcs::*;
use crate::parser::ds::*;

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
        let parser_accession = tuple((
            alphanum::narrowed_alphanumeric1,
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

        assert_eq!(
            ProteinEdit::parse("L"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "L".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("X"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "X".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("*"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "*".to_owned()
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("Leu"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "Leu".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("Ter"),
            Ok((
                "",
                ProteinEdit::Subst {
                    alternative: "Ter".to_owned()
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("delinsLeu"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "Leu".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("delinsLeuTer"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "LeuTer".to_owned()
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("delinsL"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "L".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("delinsL*"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "L*".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("delinsLX"),
            Ok((
                "",
                ProteinEdit::DelIns {
                    alternative: "LX".to_owned()
                }
            ))
        );

        assert_eq!(ProteinEdit::parse("dup"), Ok(("", ProteinEdit::Dup)));

        assert_eq!(
            ProteinEdit::parse("Lfs"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: Some("L".to_owned()),
                    terminal: None,
                    length: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("fs"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: None,
                    length: UncertainLengthChange::None
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("fs*"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("*".to_owned()),
                    length: UncertainLengthChange::None
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("fsX"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("X".to_owned()),
                    length: UncertainLengthChange::None
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("fs*?"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("*".to_owned()),
                    length: UncertainLengthChange::Unknown
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("fsX?"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("X".to_owned()),
                    length: UncertainLengthChange::Unknown
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("fs*12"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("*".to_owned()),
                    length: UncertainLengthChange::Known(12)
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("fsX12"),
            Ok((
                "",
                ProteinEdit::Fs {
                    alternative: None,
                    terminal: Some("X".to_owned()),
                    length: UncertainLengthChange::Known(12)
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("ext-1"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: None,
                    change: UncertainLengthChange::Known(-1),
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("LextM-1"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("M".to_string()),
                    change: UncertainLengthChange::Known(-1),
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("LeuextMet-1"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Met".to_string()),
                    change: UncertainLengthChange::Known(-1),
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("LextX"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("X".to_string()),
                    change: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("Lext*"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("*".to_string()),
                    change: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("LeuextTer"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("extX"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("X".to_string()),
                    change: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("ext*"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("*".to_string()),
                    change: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("extTer"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainLengthChange::None,
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("LextX?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("X".to_string()),
                    change: UncertainLengthChange::Unknown,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("Lext*?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("*".to_string()),
                    change: UncertainLengthChange::Unknown,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("LeuextTer?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainLengthChange::Unknown,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("extX?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("X".to_string()),
                    change: UncertainLengthChange::Unknown,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("ext*?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("*".to_string()),
                    change: UncertainLengthChange::Unknown,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("extTer?"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainLengthChange::Unknown,
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("LextX10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("X".to_string()),
                    change: UncertainLengthChange::Known(10),
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("Lext*10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: Some("*".to_string()),
                    change: UncertainLengthChange::Known(10),
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("LeuextTer10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainLengthChange::Known(10),
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("extX10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("X".to_string()),
                    change: UncertainLengthChange::Known(10),
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("ext*10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("*".to_string()),
                    change: UncertainLengthChange::Known(10),
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("extTer10"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: Some("Ter".to_string()),
                    change: UncertainLengthChange::Known(10),
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("Lext"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("L".to_owned()),
                    ext_aa: None,
                    change: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("Leuext"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: Some("Leu".to_owned()),
                    ext_aa: None,
                    change: UncertainLengthChange::None,
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("ext"),
            Ok((
                "",
                ProteinEdit::Ext {
                    aa_ext: None,
                    ext_aa: None,
                    change: UncertainLengthChange::None,
                }
            ))
        );

        assert_eq!(ProteinEdit::parse("="), Ok(("", ProteinEdit::Ident)));

        assert_eq!(ProteinEdit::parse("del"), Ok(("", ProteinEdit::Del)));

        assert_eq!(
            ProteinEdit::parse("insL"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "L".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("insL*"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "L*".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("insLX"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "LX".to_owned()
                }
            ))
        );

        assert_eq!(
            ProteinEdit::parse("insLeu"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "Leu".to_owned()
                }
            ))
        );
        assert_eq!(
            ProteinEdit::parse("insLeuTer"),
            Ok((
                "",
                ProteinEdit::Ins {
                    alternative: "LeuTer".to_owned()
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
