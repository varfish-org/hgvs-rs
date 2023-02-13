//! Code supporting parsing of HGVS variant specification based on the
//! `hgvs.pymeta` file from `biocommons/hgvs`.

use pest::error::Error as PestError;
use pest::Parser;

use pest_derive::Parser;

/// Parser for HGVS expressions.
///
/// The subset is limited to is limited to those rules that define sequence variants
/// precisely.  It does not current cover rules for translocations or conversions.
#[derive(Parser)]
#[grammar = "parser/hgvs.pest"]
pub struct HgvsParser;

/// The AST for HGVS expressions.
pub enum Variant<'a> {
    // CVariant { accession: Accession<'a> },
    GVariant { accession: Accession<'a> },
    // MVariant { accession: Accession<'a> },
    // NVariant { accession: Accession<'a> },
    // PVariant { accession: Accession<'a> },
    // RVariant { accession: Accession<'a> },
}

pub fn parse_variant(variant: &str) -> Result<Variant, PestError<Rule>> {
    let variant = HgvsParser::parse(Rule::hgvs_variant, variant)?
        .next()
        .unwrap();

    use pest::iterators::Pair;

    fn parse_value(pair: Pair<Rule>) -> Variant {
        match pair.as_rule() {
            Rule::g_variant => Variant::GVariant {
                accession: Accession {
                    value: pair.into_inner().next().unwrap().as_str(),
                },
            },
            Rule::hgvs_variant
            | Rule::accession
            | Rule::opt_gene_expr
            | Rule::paren_gene
            | Rule::gene_symbol => unreachable!(),
        }
    }

    Ok(parse_value(variant))
}

/// HGVS accession
pub struct Accession<'a> {
    /// The accession's string value.
    pub value: &'a str,
}

#[derive(Debug, PartialEq)]
pub enum Error {
    ParseError(String),
}

impl From<pest::error::Error<Rule>> for Error {
    fn from(value: pest::error::Error<Rule>) -> Self {
        Error::ParseError(format!("Parser error: {}", &value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pest::Parser;
    use pretty_assertions::assert_eq;

    #[test]
    fn parser_acccession() -> Result<(), Error> {
        assert_eq!(
            format!(
                "{:?}",
                HgvsParser::parse(Rule::accession, "NM_01234.5")?
            ),
            "[Pair { rule: accession, span: Span { str: \"NM_01234.5\", start: 0, end: 10 }, inner: [] }]"
        );

        assert_eq!(
            format!(
                "{:?}",
                HgvsParser::parse(Rule::accession, "LRG_01234.1")?
            ),
            "[Pair { rule: accession, span: Span { str: \"LRG_01234.1\", start: 0, end: 11 }, inner: [] }]"
        );

        Ok(())
    }

    #[test]
    fn parser_opt_gene_expr() -> Result<(), Error> {
        assert_eq!(
            format!(
                "{:?}",
                HgvsParser::parse(Rule::opt_gene_expr, "(TTN)")?
            ),
            "[Pair { rule: opt_gene_expr, span: Span { str: \"(TTN)\", start: 0, end: 5 }, inner: [] }]"
        );

        assert_eq!(
            format!("{:?}", HgvsParser::parse(Rule::opt_gene_expr, "")?),
            "[Pair { rule: opt_gene_expr, span: Span { str: \"\", start: 0, end: 0 }, inner: [] }]"
        );

        Ok(())
    }
}
