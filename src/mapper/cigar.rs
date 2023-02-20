//! Code supporting the `CigarMapper`

use std::fmt::Display;

use nom::{combinator::all_consuming, multi::many0};

/// CIGAR operation as parsed from UTA.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum CigarOp {
    /// =
    Eq,
    /// D
    Del,
    /// I
    Ins,
    /// M
    Match,
    /// S
    Skip,
    /// X
    Mismatch,
}

impl CigarOp {
    pub fn is_advance_ref(&self) -> bool {
        matches!(
            self,
            CigarOp::Eq | CigarOp::Match | CigarOp::Mismatch | CigarOp::Ins | CigarOp::Skip
        )
    }

    pub fn is_advance_tgt(&self) -> bool {
        matches!(
            self,
            CigarOp::Eq | CigarOp::Match | CigarOp::Mismatch | CigarOp::Del
        )
    }
}

impl TryFrom<char> for CigarOp {
    type Error = anyhow::Error;

    fn try_from(value: char) -> Result<Self, anyhow::Error> {
        Ok(match value {
            '=' => Self::Eq,
            'D' => Self::Del,
            'I' => Self::Ins,
            'M' => Self::Match,
            'N' => Self::Skip,
            'X' => Self::Mismatch,
            _ => return Err(anyhow::anyhow!("Invalid CIGAR character {}", value)),
        })
    }
}

impl From<CigarOp> for char {
    fn from(val: CigarOp) -> Self {
        match val {
            CigarOp::Eq => '=',
            CigarOp::Del => 'D',
            CigarOp::Ins => 'I',
            CigarOp::Match => 'M',
            CigarOp::Skip => 'N',
            CigarOp::Mismatch => 'X',
        }
    }
}

impl Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", std::convert::Into::<char>::into(*self))
    }
}

/// CIGAR element consisting of count and CIGAR operation.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct CigarElement {
    pub count: i32,
    pub op: CigarOp,
}

impl Display for CigarElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.count > 1 {
            write!(f, "{}", self.count)?;
        }
        write!(f, "{}", self.op)
    }
}

impl CigarElement {
    fn from_strs(count: &str, op: &str) -> CigarElement {
        CigarElement {
            count: if count.is_empty() {
                1
            } else {
                str::parse(count).unwrap()
            },
            op: op.chars().next().unwrap().try_into().unwrap(),
        }
    }
}

#[derive(Debug, PartialEq, PartialOrd, Default, Clone)]
pub struct CigarString {
    pub elems: Vec<CigarElement>,
}

impl CigarString {
    fn from(elems: Vec<CigarElement>) -> Self {
        Self { elems }
    }
}

impl std::ops::Deref for CigarString {
    type Target = Vec<CigarElement>;
    fn deref(&self) -> &Self::Target {
        &self.elems
    }
}
impl std::ops::DerefMut for CigarString {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.elems
    }
}

impl Display for CigarString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for item in &self.elems {
            write!(f, "{}", &item)?
        }
        Ok(())
    }
}

pub mod parse {
    use nom::{
        bytes::complete::take_while_m_n,
        character::complete::digit0,
        error::{context, VerboseError},
        sequence::pair,
        IResult,
    };

    type Res<T, U> = IResult<T, U, VerboseError<T>>;

    use super::CigarElement;

    pub fn is_cigar_op_char(c: char) -> bool {
        "=DIMNX".contains(c)
    }

    pub fn cigar_element(input: &str) -> Res<&str, CigarElement> {
        context(
            "cigar_element",
            pair(digit0, take_while_m_n(1, 1, is_cigar_op_char)),
        )(input)
        .map(|(rest, (count, op))| (rest, CigarElement::from_strs(count, op)))
    }
}

/// Parse a CIGAR `str` into a real one.
pub fn parse_cigar_string(input: &str) -> Result<CigarString, anyhow::Error> {
    Ok(CigarString::from(
        all_consuming(many0(parse::cigar_element))(input)
            .map_err(|e| anyhow::anyhow!("Problem with parsing: {:?}", e))?
            .1,
    ))
}

#[cfg(test)]
mod test {
    use crate::mapper::cigar::{parse_cigar_string, CigarElement, CigarOp};

    #[test]
    fn parse_cigar_string_simple() -> Result<(), anyhow::Error> {
        // assert_eq!(parse_cigar_string("")?, vec![]);
        assert_eq!(
            parse_cigar_string("M")?.elems,
            vec![CigarElement {
                count: 1,
                op: CigarOp::Match
            }]
        );
        assert_eq!(
            parse_cigar_string("MM")?.elems,
            vec![
                CigarElement {
                    count: 1,
                    op: CigarOp::Match
                },
                CigarElement {
                    count: 1,
                    op: CigarOp::Match
                }
            ]
        );
        assert_eq!(
            parse_cigar_string("1M")?.elems,
            vec![CigarElement {
                count: 1,
                op: CigarOp::Match,
            },]
        );
        assert_eq!(
            parse_cigar_string("1M2I3X")?.elems,
            vec![
                CigarElement {
                    count: 1,
                    op: CigarOp::Match,
                },
                CigarElement {
                    count: 2,
                    op: CigarOp::Ins,
                },
                CigarElement {
                    count: 3,
                    op: CigarOp::Mismatch,
                },
            ]
        );
        assert_eq!(
            parse_cigar_string("1MI3X")?.elems,
            vec![
                CigarElement {
                    count: 1,
                    op: CigarOp::Match,
                },
                CigarElement {
                    count: 1,
                    op: CigarOp::Ins,
                },
                CigarElement {
                    count: 3,
                    op: CigarOp::Mismatch,
                },
            ]
        );

        Ok(())
    }
}
