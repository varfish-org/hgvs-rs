//! Implementation of parser functions.

/// Code for parsing alpha/numeric strings.
pub mod alphanum {
    use nom::character::complete::alphanumeric1;

    // cf. https://stackoverflow.com/a/73437782/84349
    pub fn narrowed_alphanumeric1(
        input: &str,
    ) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        alphanumeric1(input)
    }
}

/// Code for parsing amino acid residues and protein sequences.
#[allow(dead_code)] // NB: code is tested welll
pub mod protein {
    use nom::{
        bytes::complete::take,
        multi::{many0, many1},
    };

    pub static AA1: &str = "ACDEFGHIKLMNPQRSTVWYBZXU";

    pub fn aa1(input: &str) -> Result<(&str, &str), nom::Err<nom::error::Error<&str>>> {
        let (rest, c) = take(1usize)(input)?;
        if !AA1.contains(c) {
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
        if !AAT1.contains(c) {
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

    pub const AA3: &[&str] = &[
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

    pub const AAT3: &[&str] = &[
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

/// Code for parsing protein edits.
pub mod protein_edit {
    use nom::branch::alt;
    use nom::character::complete::{digit0, digit1};
    use nom::combinator::opt;

    use nom::sequence::pair;
    use nom::{bytes::complete::tag, IResult};

    use super::protein::*;
    use crate::parser::ds::{ProteinEdit, UncertainLengthChange};

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
                            UncertainLengthChange::Unknown
                        } else {
                            UncertainLengthChange::None
                        },
                    },
                ))
            } else {
                Ok((
                    rest,
                    ProteinEdit::Fs {
                        alternative: alternative.map(str::to_owned),
                        terminal: Some(terminal.to_string()),
                        length: UncertainLengthChange::Known(str::parse::<i32>(count).unwrap()),
                    },
                ))
            }
        } else {
            Ok((
                rest,
                ProteinEdit::Fs {
                    alternative: alternative.map(str::to_owned),
                    terminal: None,
                    length: UncertainLengthChange::None,
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
                change: UncertainLengthChange::Known(-str::parse::<i32>(offset).unwrap()),
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
                        UncertainLengthChange::Unknown
                    } else {
                        UncertainLengthChange::Known(str::parse::<i32>(offset).unwrap())
                    },
                },
            ))
        } else {
            Ok((
                rest,
                ProteinEdit::Ext {
                    aa_ext: aa_ext.map(str::to_owned),
                    ext_aa: Some(ext_aa.to_string()),
                    change: UncertainLengthChange::None,
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
                change: UncertainLengthChange::None,
            },
        ))
    }

    pub fn ident(input: &str) -> IResult<&str, ProteinEdit> {
        let (rest, _) = tag("=")(input)?;
        Ok((rest, ProteinEdit::Ident))
    }
}

/// Functions for parsing nucleic acid residues and sequences.
pub mod na {
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

/// Functions for parsing nucleic acid edits.
pub mod na_edit {
    use nom::bytes::complete::tag;
    use nom::character::complete::{char as nom_char, digit1};
    use nom::combinator::map;
    use nom::sequence::tuple;
    use nom::{multi::many0, sequence::pair, IResult};

    use crate::parser::NaEdit;

    use super::na::{na, na0, na1};

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

    pub fn del_ref(input: &str) -> IResult<&str, NaEdit> {
        map(tuple((tag("del"), na0)), |(_, reference)| NaEdit::DelRef {
            reference: reference.to_string(),
        })(input)
    }

    pub fn del_num(input: &str) -> IResult<&str, NaEdit> {
        map(tuple((tag("del"), digit1)), |(_, count)| NaEdit::DelNum {
            count: str::parse::<i32>(count).unwrap(),
        })(input)
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
                count: count.parse::<i32>().unwrap(),
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
                count: count.parse::<i32>().unwrap(),
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

/// Parsing of CDS position and interval.
pub mod cds_pos {
    use nom::{
        branch::alt,
        bytes::complete::tag,
        character::complete::digit1,
        combinator::{map, opt, recognize},
        sequence::{pair, tuple},
        IResult,
    };

    use crate::parser::{CdsFrom, CdsInterval, CdsPos};

    fn pos_from_start(input: &str) -> IResult<&str, CdsPos> {
        let (rest, base) = recognize(pair(opt(alt((tag("+"), tag("-")))), digit1))(input)?;
        let (rest, offset) = opt(recognize(pair(alt((tag("+"), tag("-"))), digit1)))(rest)?;
        Ok((
            rest,
            CdsPos {
                base: str::parse::<i32>(base).unwrap(),
                offset: offset.map(|offset| str::parse::<i32>(offset).unwrap()),
                cds_from: CdsFrom::Start,
            },
        ))
    }

    fn pos_from_end(input: &str) -> IResult<&str, CdsPos> {
        let (rest, _) = tag("*")(input)?;
        let (rest, base) = digit1(rest)?;
        let (rest, offset) = opt(recognize(pair(alt((tag("+"), tag("-"))), digit1)))(rest)?;
        Ok((
            rest,
            CdsPos {
                base: str::parse::<i32>(base).unwrap(),
                offset: offset.map(|offset| str::parse::<i32>(offset).unwrap()),
                cds_from: CdsFrom::End,
            },
        ))
    }

    pub fn pos(input: &str) -> IResult<&str, CdsPos> {
        alt((pos_from_start, pos_from_end))(input)
    }

    pub fn int(input: &str) -> IResult<&str, CdsInterval> {
        let (rest, (begin, _, end)) = tuple((pos, tag("_"), pos))(input)?;
        Ok((rest, CdsInterval { start: begin, end }))
    }

    pub fn loc(input: &str) -> IResult<&str, CdsInterval> {
        alt((
            int,
            map(pos, |pos| CdsInterval {
                start: pos.clone(),
                end: pos,
            }),
        ))(input)
    }
}

/// Parsing of Genome position and interval.
pub mod genome_pos {
    use nom::{
        branch::alt, bytes::complete::tag, character::complete::digit1, combinator::map,
        sequence::tuple, IResult,
    };

    use crate::parser::GenomeInterval;

    pub fn pos(input: &str) -> IResult<&str, Option<i32>> {
        map(alt((digit1, tag("?"))), |value| {
            if value == "?" {
                None
            } else {
                Some(str::parse::<i32>(value).unwrap())
            }
        })(input)
    }

    pub fn int(input: &str) -> IResult<&str, GenomeInterval> {
        map(tuple((pos, tag("_"), pos)), |(begin, _, end)| {
            GenomeInterval { start: begin, end }
        })(input)
    }

    pub fn loc(input: &str) -> IResult<&str, GenomeInterval> {
        alt((
            int,
            map(pos, |pos| GenomeInterval {
                start: pos,
                end: pos,
            }),
        ))(input)
    }
}

/// Parsing of mt position and interval.
pub mod mt_pos {
    use nom::{
        branch::alt, bytes::complete::tag, character::complete::digit1, combinator::map,
        sequence::tuple, IResult,
    };

    use crate::parser::MtInterval;

    pub fn pos(input: &str) -> IResult<&str, Option<i32>> {
        map(alt((digit1, tag("?"))), |value| {
            if value == "?" {
                None
            } else {
                Some(str::parse::<i32>(value).unwrap())
            }
        })(input)
    }

    pub fn int(input: &str) -> IResult<&str, MtInterval> {
        map(tuple((pos, tag("_"), pos)), |(begin, _, end)| MtInterval {
            start: begin,
            end,
        })(input)
    }

    pub fn loc(input: &str) -> IResult<&str, MtInterval> {
        alt((
            int,
            map(pos, |pos| MtInterval {
                start: pos,
                end: pos,
            }),
        ))(input)
    }
}

/// Parsing of transcript position and interval.
pub mod tx_pos {
    use nom::{
        branch::alt,
        bytes::complete::tag,
        character::complete::digit1,
        combinator::{map, opt, recognize},
        sequence::{pair, tuple},
        IResult,
    };

    use crate::parser::{TxInterval, TxPos};

    pub fn pos(input: &str) -> IResult<&str, TxPos> {
        let (rest, base) = recognize(pair(opt(alt((tag("+"), tag("-")))), digit1))(input)?;
        let (rest, offset) = opt(recognize(pair(alt((tag("+"), tag("-"))), digit1)))(rest)?;
        Ok((
            rest,
            TxPos {
                base: str::parse::<i32>(base).unwrap(),
                offset: offset.map(|offset| str::parse::<i32>(offset).unwrap()),
            },
        ))
    }

    pub fn int(input: &str) -> IResult<&str, TxInterval> {
        let (rest, (begin, _, end)) = tuple((pos, tag("_"), pos))(input)?;
        Ok((rest, TxInterval { start: begin, end }))
    }

    pub fn loc(input: &str) -> IResult<&str, TxInterval> {
        alt((
            int,
            map(pos, |pos| TxInterval {
                start: pos.clone(),
                end: pos,
            }),
        ))(input)
    }
}

/// Parsing of RNA position and interval.
pub mod rna_pos {
    use nom::{
        branch::alt,
        bytes::complete::tag,
        character::complete::digit1,
        combinator::{map, opt, recognize},
        sequence::{pair, tuple},
        IResult,
    };

    use crate::parser::{RnaInterval, RnaPos};

    pub fn pos(input: &str) -> IResult<&str, RnaPos> {
        let (rest, base) = recognize(pair(opt(alt((tag("+"), tag("-")))), digit1))(input)?;
        let (rest, offset) = opt(recognize(pair(alt((tag("+"), tag("-"))), digit1)))(rest)?;
        Ok((
            rest,
            RnaPos {
                base: str::parse::<i32>(base).unwrap(),
                offset: offset.map(|offset| str::parse::<i32>(offset).unwrap()),
            },
        ))
    }

    pub fn int(input: &str) -> IResult<&str, RnaInterval> {
        let (rest, (begin, _, end)) = tuple((pos, tag("_"), pos))(input)?;
        Ok((rest, RnaInterval { start: begin, end }))
    }

    pub fn loc(input: &str) -> IResult<&str, RnaInterval> {
        alt((
            int,
            map(pos, |pos| RnaInterval {
                start: pos.clone(),
                end: pos,
            }),
        ))(input)
    }
}

/// Parsing of protein position and interval
pub mod prot_pos {
    use nom::{
        branch::alt,
        bytes::complete::tag,
        character::complete::digit1,
        combinator::map,
        sequence::{pair, tuple},
        IResult,
    };

    use crate::parser::{ProtInterval, ProtPos};

    use super::protein::{aat1, aat3};

    pub fn pos(input: &str) -> IResult<&str, ProtPos> {
        let (rest, (aa, number)) = pair(alt((aat3, aat1)), digit1)(input)?;
        Ok((
            rest,
            ProtPos {
                aa: aa.to_string(),
                number: str::parse::<i32>(number).unwrap(),
            },
        ))
    }

    pub fn int(input: &str) -> IResult<&str, ProtInterval> {
        let (rest, (begin, _, end)) = tuple((pos, tag("_"), pos))(input)?;
        Ok((rest, ProtInterval { start: begin, end }))
    }

    pub fn loc(input: &str) -> IResult<&str, ProtInterval> {
        alt((
            int,
            map(pos, |pos| ProtInterval {
                start: pos.clone(),
                end: pos,
            }),
        ))(input)
    }
}

#[cfg(test)]
mod test {
    use crate::parser::{
        ds::{ProteinEdit, UncertainLengthChange},
        CdsFrom, CdsInterval, CdsPos, GenomeInterval, MtInterval, NaEdit, ProtInterval, ProtPos,
        RnaInterval, RnaPos, TxInterval, TxPos,
    };

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
                    length: UncertainLengthChange::None,
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
                    length: UncertainLengthChange::None
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
                    length: UncertainLengthChange::None
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
                    length: UncertainLengthChange::None
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
                    length: UncertainLengthChange::Unknown
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
                    length: UncertainLengthChange::Unknown
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
                    length: UncertainLengthChange::Known(12)
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
                    length: UncertainLengthChange::Known(12)
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
                    change: UncertainLengthChange::Known(-1),
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
                    change: UncertainLengthChange::Known(-1),
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
                    change: UncertainLengthChange::Known(-1),
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::Unknown,
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
                    change: UncertainLengthChange::Unknown,
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
                    change: UncertainLengthChange::Unknown,
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
                    change: UncertainLengthChange::Unknown,
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
                    change: UncertainLengthChange::Unknown,
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
                    change: UncertainLengthChange::Unknown,
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
                    change: UncertainLengthChange::Known(10),
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
                    change: UncertainLengthChange::Known(10),
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
                    change: UncertainLengthChange::Known(10),
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
                    change: UncertainLengthChange::Known(10),
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
                    change: UncertainLengthChange::Known(10),
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
                    change: UncertainLengthChange::Known(10),
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::None,
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
                    change: UncertainLengthChange::None,
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
    fn naedit_del_ref() {
        assert_eq!(
            na_edit::del_ref("delT"),
            Ok((
                "",
                NaEdit::DelRef {
                    reference: "T".to_owned(),
                }
            ))
        );
        assert_eq!(
            na_edit::del_ref("del"),
            Ok((
                "",
                NaEdit::DelRef {
                    reference: "".to_owned(),
                }
            ))
        );
    }

    #[test]
    fn naedit_del_num() {
        assert_eq!(
            na_edit::del_num("del3"),
            Ok(("", NaEdit::DelNum { count: 3 }))
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
    fn mtpos_loc() {
        assert_eq!(
            mt_pos::loc("42"),
            Ok((
                "",
                MtInterval {
                    start: Some(42),
                    end: Some(42),
                }
            ))
        );
    }

    #[test]
    fn mtpos_int() {
        assert_eq!(
            mt_pos::int("42_100"),
            Ok((
                "",
                MtInterval {
                    start: Some(42),
                    end: Some(100),
                }
            ))
        );
        assert_eq!(
            mt_pos::int("?_100"),
            Ok((
                "",
                MtInterval {
                    start: None,
                    end: Some(100),
                }
            ))
        );
        assert_eq!(
            mt_pos::int("42_?"),
            Ok((
                "",
                MtInterval {
                    start: Some(42),
                    end: None,
                }
            ))
        );
        assert_eq!(
            mt_pos::int("?_?"),
            Ok((
                "",
                MtInterval {
                    start: None,
                    end: None,
                }
            ))
        );
    }

    #[test]
    fn cdspos_loc() {
        assert_eq!(
            cds_pos::loc("123"),
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
    fn cdspos_int() {
        assert_eq!(
            cds_pos::int("123_123"),
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
        assert_eq!(
            cds_pos::int("+123_+123"),
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
        assert_eq!(
            cds_pos::int("-123_-123"),
            Ok((
                "",
                CdsInterval {
                    start: CdsPos {
                        base: -123,
                        offset: None,
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: -123,
                        offset: None,
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
        assert_eq!(
            cds_pos::int("123+42_123+42"),
            Ok((
                "",
                CdsInterval {
                    start: CdsPos {
                        base: 123,
                        offset: Some(42),
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: 123,
                        offset: Some(42),
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
        assert_eq!(
            cds_pos::int("+123+42_+123+42"),
            Ok((
                "",
                CdsInterval {
                    start: CdsPos {
                        base: 123,
                        offset: Some(42),
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: 123,
                        offset: Some(42),
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
        assert_eq!(
            cds_pos::int("-123+42_-123+42"),
            Ok((
                "",
                CdsInterval {
                    start: CdsPos {
                        base: -123,
                        offset: Some(42),
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: -123,
                        offset: Some(42),
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
        assert_eq!(
            cds_pos::int("123-42_123-42"),
            Ok((
                "",
                CdsInterval {
                    start: CdsPos {
                        base: 123,
                        offset: Some(-42),
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: 123,
                        offset: Some(-42),
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
        assert_eq!(
            cds_pos::int("+123-42_+123-42"),
            Ok((
                "",
                CdsInterval {
                    start: CdsPos {
                        base: 123,
                        offset: Some(-42),
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: 123,
                        offset: Some(-42),
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
        assert_eq!(
            cds_pos::int("-123-42_-123-42"),
            Ok((
                "",
                CdsInterval {
                    start: CdsPos {
                        base: -123,
                        offset: Some(-42),
                        cds_from: CdsFrom::Start,
                    },
                    end: CdsPos {
                        base: -123,
                        offset: Some(-42),
                        cds_from: CdsFrom::Start,
                    },
                }
            ))
        );
    }

    #[test]
    fn cdspos_pos() {
        assert_eq!(
            cds_pos::pos("123"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: None,
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("+123"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: None,
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("-123"),
            Ok((
                "",
                CdsPos {
                    base: -123,
                    offset: None,
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("123+42"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: Some(42),
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("+123+42"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: Some(42),
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("-123+42"),
            Ok((
                "",
                CdsPos {
                    base: -123,
                    offset: Some(42),
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("123-42"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: Some(-42),
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("+123-42"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: Some(-42),
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("-123-42"),
            Ok((
                "",
                CdsPos {
                    base: -123,
                    offset: Some(-42),
                    cds_from: CdsFrom::Start,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("*123"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: None,
                    cds_from: CdsFrom::End,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("*123+42"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: Some(42),
                    cds_from: CdsFrom::End,
                }
            ))
        );
        assert_eq!(
            cds_pos::pos("*123-42"),
            Ok((
                "",
                CdsPos {
                    base: 123,
                    offset: Some(-42),
                    cds_from: CdsFrom::End,
                }
            ))
        );
    }

    #[test]
    fn genomepos_loc() {
        assert_eq!(
            genome_pos::loc("42"),
            Ok((
                "",
                GenomeInterval {
                    start: Some(42),
                    end: Some(42),
                }
            ))
        );
    }

    #[test]
    fn genomepos_int() {
        assert_eq!(
            genome_pos::int("42_100"),
            Ok((
                "",
                GenomeInterval {
                    start: Some(42),
                    end: Some(100),
                }
            ))
        );
        assert_eq!(
            genome_pos::int("?_100"),
            Ok((
                "",
                GenomeInterval {
                    start: None,
                    end: Some(100),
                }
            ))
        );
        assert_eq!(
            genome_pos::int("42_?"),
            Ok((
                "",
                GenomeInterval {
                    start: Some(42),
                    end: None,
                }
            ))
        );
        assert_eq!(
            genome_pos::int("?_?"),
            Ok((
                "",
                GenomeInterval {
                    start: None,
                    end: None,
                }
            ))
        );
    }

    #[test]
    fn rnapos_loc() {
        assert_eq!(
            tx_pos::loc("123"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: 123,
                        offset: None,
                    },
                    end: TxPos {
                        base: 123,
                        offset: None,
                    },
                }
            ))
        );
    }

    #[test]
    fn rnapos_int() {
        assert_eq!(
            tx_pos::int("123_123"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: 123,
                        offset: None,
                    },
                    end: TxPos {
                        base: 123,
                        offset: None,
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("+123_+123"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: 123,
                        offset: None,
                    },
                    end: TxPos {
                        base: 123,
                        offset: None,
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("-123_-123"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: -123,
                        offset: None,
                    },
                    end: TxPos {
                        base: -123,
                        offset: None,
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("123+42_123+42"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: 123,
                        offset: Some(42),
                    },
                    end: TxPos {
                        base: 123,
                        offset: Some(42),
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("+123+42_+123+42"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: 123,
                        offset: Some(42),
                    },
                    end: TxPos {
                        base: 123,
                        offset: Some(42),
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("-123+42_-123+42"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: -123,
                        offset: Some(42),
                    },
                    end: TxPos {
                        base: -123,
                        offset: Some(42),
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("123-42_123-42"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: 123,
                        offset: Some(-42),
                    },
                    end: TxPos {
                        base: 123,
                        offset: Some(-42),
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("+123-42_+123-42"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: 123,
                        offset: Some(-42),
                    },
                    end: TxPos {
                        base: 123,
                        offset: Some(-42),
                    },
                }
            ))
        );
        assert_eq!(
            tx_pos::int("-123-42_-123-42"),
            Ok((
                "",
                TxInterval {
                    start: TxPos {
                        base: -123,
                        offset: Some(-42),
                    },
                    end: TxPos {
                        base: -123,
                        offset: Some(-42),
                    },
                }
            ))
        );
    }

    #[test]
    fn txpos_pos() {
        assert_eq!(
            tx_pos::pos("123"),
            Ok((
                "",
                TxPos {
                    base: 123,
                    offset: None,
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("+123"),
            Ok((
                "",
                TxPos {
                    base: 123,
                    offset: None,
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("-123"),
            Ok((
                "",
                TxPos {
                    base: -123,
                    offset: None,
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("123+42"),
            Ok((
                "",
                TxPos {
                    base: 123,
                    offset: Some(42),
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("+123+42"),
            Ok((
                "",
                TxPos {
                    base: 123,
                    offset: Some(42),
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("-123+42"),
            Ok((
                "",
                TxPos {
                    base: -123,
                    offset: Some(42),
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("123-42"),
            Ok((
                "",
                TxPos {
                    base: 123,
                    offset: Some(-42),
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("+123-42"),
            Ok((
                "",
                TxPos {
                    base: 123,
                    offset: Some(-42),
                }
            ))
        );
        assert_eq!(
            tx_pos::pos("-123-42"),
            Ok((
                "",
                TxPos {
                    base: -123,
                    offset: Some(-42),
                }
            ))
        );
    }

    #[test]
    fn txpos_loc() {
        assert_eq!(
            rna_pos::loc("123"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: 123,
                        offset: None,
                    },
                    end: RnaPos {
                        base: 123,
                        offset: None,
                    },
                }
            ))
        );
    }

    #[test]
    fn txpos_int() {
        assert_eq!(
            rna_pos::int("123_123"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: 123,
                        offset: None,
                    },
                    end: RnaPos {
                        base: 123,
                        offset: None,
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("+123_+123"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: 123,
                        offset: None,
                    },
                    end: RnaPos {
                        base: 123,
                        offset: None,
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("-123_-123"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: -123,
                        offset: None,
                    },
                    end: RnaPos {
                        base: -123,
                        offset: None,
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("123+42_123+42"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: 123,
                        offset: Some(42),
                    },
                    end: RnaPos {
                        base: 123,
                        offset: Some(42),
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("+123+42_+123+42"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: 123,
                        offset: Some(42),
                    },
                    end: RnaPos {
                        base: 123,
                        offset: Some(42),
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("-123+42_-123+42"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: -123,
                        offset: Some(42),
                    },
                    end: RnaPos {
                        base: -123,
                        offset: Some(42),
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("123-42_123-42"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: 123,
                        offset: Some(-42),
                    },
                    end: RnaPos {
                        base: 123,
                        offset: Some(-42),
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("+123-42_+123-42"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: 123,
                        offset: Some(-42),
                    },
                    end: RnaPos {
                        base: 123,
                        offset: Some(-42),
                    },
                }
            ))
        );
        assert_eq!(
            rna_pos::int("-123-42_-123-42"),
            Ok((
                "",
                RnaInterval {
                    start: RnaPos {
                        base: -123,
                        offset: Some(-42),
                    },
                    end: RnaPos {
                        base: -123,
                        offset: Some(-42),
                    },
                }
            ))
        );
    }

    #[test]
    fn rnapos_pos() {
        assert_eq!(
            rna_pos::pos("123"),
            Ok((
                "",
                RnaPos {
                    base: 123,
                    offset: None,
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("+123"),
            Ok((
                "",
                RnaPos {
                    base: 123,
                    offset: None,
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("-123"),
            Ok((
                "",
                RnaPos {
                    base: -123,
                    offset: None,
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("123+42"),
            Ok((
                "",
                RnaPos {
                    base: 123,
                    offset: Some(42),
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("+123+42"),
            Ok((
                "",
                RnaPos {
                    base: 123,
                    offset: Some(42),
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("-123+42"),
            Ok((
                "",
                RnaPos {
                    base: -123,
                    offset: Some(42),
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("123-42"),
            Ok((
                "",
                RnaPos {
                    base: 123,
                    offset: Some(-42),
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("+123-42"),
            Ok((
                "",
                RnaPos {
                    base: 123,
                    offset: Some(-42),
                }
            ))
        );
        assert_eq!(
            rna_pos::pos("-123-42"),
            Ok((
                "",
                RnaPos {
                    base: -123,
                    offset: Some(-42),
                }
            ))
        );
    }

    #[test]
    fn protpos_pos() {
        assert_eq!(
            prot_pos::pos("Leu123"),
            Ok((
                "",
                ProtPos {
                    aa: "Leu".to_string(),
                    number: 123
                }
            ))
        );
        assert_eq!(
            prot_pos::pos("L123"),
            Ok((
                "",
                ProtPos {
                    aa: "L".to_string(),
                    number: 123
                }
            ))
        );
        assert_eq!(
            prot_pos::pos("Ter123"),
            Ok((
                "",
                ProtPos {
                    aa: "Ter".to_string(),
                    number: 123
                }
            ))
        );
        assert_eq!(
            prot_pos::pos("*123"),
            Ok((
                "",
                ProtPos {
                    aa: "*".to_string(),
                    number: 123
                }
            ))
        );
        assert_eq!(
            prot_pos::pos("X123"),
            Ok((
                "",
                ProtPos {
                    aa: "X".to_string(),
                    number: 123
                }
            ))
        );
    }

    #[test]
    fn protpos_loc() {
        assert_eq!(
            prot_pos::loc("Leu10"),
            Ok((
                "",
                ProtInterval {
                    start: ProtPos {
                        aa: "Leu".to_string(),
                        number: 10
                    },
                    end: ProtPos {
                        aa: "Leu".to_string(),
                        number: 10
                    }
                },
            ))
        );
    }

    #[test]
    fn protpos_int() {
        assert_eq!(
            prot_pos::int("Leu10_Leu123"),
            Ok((
                "",
                ProtInterval {
                    start: ProtPos {
                        aa: "Leu".to_string(),
                        number: 10
                    },
                    end: ProtPos {
                        aa: "Leu".to_string(),
                        number: 123
                    }
                },
            ))
        );
        assert_eq!(
            prot_pos::int("L10_L123"),
            Ok((
                "",
                ProtInterval {
                    start: ProtPos {
                        aa: "L".to_string(),
                        number: 10
                    },
                    end: ProtPos {
                        aa: "L".to_string(),
                        number: 123
                    }
                },
            ))
        );
        assert_eq!(
            prot_pos::int("Leu10_Ter123"),
            Ok((
                "",
                ProtInterval {
                    start: ProtPos {
                        aa: "Leu".to_string(),
                        number: 10
                    },
                    end: ProtPos {
                        aa: "Ter".to_string(),
                        number: 123
                    }
                },
            ))
        );
        assert_eq!(
            prot_pos::int("L10_*123"),
            Ok((
                "",
                ProtInterval {
                    start: ProtPos {
                        aa: "L".to_string(),
                        number: 10
                    },
                    end: ProtPos {
                        aa: "*".to_string(),
                        number: 123
                    }
                },
            ))
        );
        assert_eq!(
            prot_pos::int("L10_X123"),
            Ok((
                "",
                ProtInterval {
                    start: ProtPos {
                        aa: "L".to_string(),
                        number: 10
                    },
                    end: ProtPos {
                        aa: "X".to_string(),
                        number: 123
                    }
                },
            ))
        );
    }
}
