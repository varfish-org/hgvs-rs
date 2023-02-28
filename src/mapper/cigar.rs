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

/// Provide coordinate mapping between two sequences whose alignment is given by a CIGAR string.
///
/// CIGAR is about alignments between positions in two sequences.  It is base-centric.
///
/// Unfortunately, base-centric coordinate systems require additional complexity to refer to
/// zero-width positions.
///
/// This code uses interbase intervals.  Interbase positions are zero-width boundaries between
/// bases.  They often look similar to zero-based, right open coordinates. (But don't call them
/// that.  It upsets me deeply.)  The most important difference is that zero width intervals
/// neatly represent insertions between bases (or before or after the sequence).
#[derive(Default, Debug)]
pub struct CigarMapper {
    pub cigar_string: CigarString,
    pub ref_pos: Vec<i32>,
    pub tgt_pos: Vec<i32>,
    pub cigar_op: Vec<CigarOp>,
    pub ref_len: i32,
    pub tgt_len: i32,
}

#[derive(Debug, PartialEq)]
pub struct CigarMapperResult {
    pub pos: i32,
    pub offset: i32,
    pub cigar_op: CigarOp,
}

impl CigarMapper {
    pub fn new(cigar_string: &CigarString) -> Self {
        let (ref_pos, tgt_pos, cigar_op) = Self::init(cigar_string);

        Self {
            cigar_string: cigar_string.clone(),
            ref_len: *ref_pos.last().unwrap(),
            tgt_len: *tgt_pos.last().unwrap(),
            ref_pos,
            tgt_pos,
            cigar_op,
        }
    }

    /// For a given CIGAR string, return the start positions of each aligned segment in ref
    /// and tgt, and a list of CIGAR operators.
    fn init(cigar_string: &CigarString) -> (Vec<i32>, Vec<i32>, Vec<CigarOp>) {
        let cigar_len = cigar_string.len();

        let mut ref_pos = vec![-1; cigar_len];
        let mut tgt_pos = vec![-1; cigar_len];
        let mut cigar_op = vec![CigarOp::Mismatch; cigar_len];
        let mut ref_cur = 0;
        let mut tgt_cur = 0;
        for (i, CigarElement { count, op }) in cigar_string.iter().enumerate() {
            ref_pos[i] = ref_cur;
            tgt_pos[i] = tgt_cur;
            cigar_op[i] = *op;
            if op.is_advance_ref() {
                ref_cur += *count;
            }
            if op.is_advance_tgt() {
                tgt_cur += *count;
            }
        }
        ref_pos.push(ref_cur);
        tgt_pos.push(tgt_cur);

        (ref_pos, tgt_pos, cigar_op)
    }

    pub fn map_ref_to_tgt(
        &self,
        pos: i32,
        end: &str,
        strict_bounds: bool,
    ) -> Result<CigarMapperResult, anyhow::Error> {
        self.map(&self.ref_pos, &self.tgt_pos, pos, end, strict_bounds)
    }

    pub fn map_tgt_to_ref(
        &self,
        pos: i32,
        end: &str,
        strict_bounds: bool,
    ) -> Result<CigarMapperResult, anyhow::Error> {
        self.map(&self.tgt_pos, &self.ref_pos, pos, end, strict_bounds)
    }

    /// Map position between aligned segments.
    ///
    /// Positions in this function are 0-based, base-counting.
    fn map(
        &self,
        from_pos: &[i32],
        to_pos: &[i32],
        pos: i32,
        end: &str,
        strict_bounds: bool,
    ) -> Result<CigarMapperResult, anyhow::Error> {
        if strict_bounds && (pos < 0 || pos > *from_pos.last().unwrap()) {
            return Err(anyhow::anyhow!(
                "Position is beyond the bounds of transcript record (pos={}, from_pos={:?}, to_pos={:?}]",
                pos,
                from_pos,
                to_pos,
            ));
        }

        // Find aligned segment to use as basis for mapping.  It is okay for pos to be
        // before first element or after last.
        let pos_i = {
            let mut pos_i = 0;
            while pos_i < self.cigar_op.len() {
                if pos < from_pos[pos_i + 1] {
                    break;
                }
                pos_i += 1;
            }
            std::cmp::min(pos_i, self.cigar_op.len().saturating_sub(1))
        };

        let cigar_op = self.cigar_op[pos_i];

        if cigar_op == CigarOp::Eq || cigar_op == CigarOp::Match || cigar_op == CigarOp::Mismatch {
            Ok(CigarMapperResult {
                pos: to_pos[pos_i] + (pos - from_pos[pos_i]),
                offset: 0,
                cigar_op,
            })
        } else if cigar_op == CigarOp::Del || cigar_op == CigarOp::Ins {
            Ok(CigarMapperResult {
                pos: if end == "start" {
                    to_pos[pos_i] - 1
                } else {
                    to_pos[pos_i]
                },
                offset: 0,
                cigar_op,
            })
        } else if cigar_op == CigarOp::Skip {
            if pos - from_pos[pos_i] < from_pos[pos_i + 1] - pos {
                Ok(CigarMapperResult {
                    pos: to_pos[pos_i] - 1,
                    offset: pos - from_pos[pos_i] + 1,
                    cigar_op,
                })
            } else {
                Ok(CigarMapperResult {
                    pos: to_pos[pos_i],
                    offset: -(from_pos[pos_i + 1] - pos),
                    cigar_op,
                })
            }
        } else {
            Err(anyhow::anyhow!("Algorithm error in CIGAR mapper"))
        }
    }
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use super::{parse_cigar_string, CigarElement, CigarMapper, CigarMapperResult, CigarOp};

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

    #[test]
    fn cigar_mapper_simple() -> Result<(), anyhow::Error> {
        // 0   1   2           3   4   5               6       7   8   9  tgt
        // =   =   =   N   N   =   X   =   N   N   N   =   I   =   D   =
        // 0   1   2   3   4   5   6   7   8   9  10  11  12  13      14  ref
        let cigar = "3=2N=X=3N=I=D=".to_string();
        let cigar_str = parse_cigar_string(&cigar)?;
        let cigar_mapper = CigarMapper::new(&cigar_str);

        assert_eq!(cigar_mapper.ref_len, 15);
        assert_eq!(cigar_mapper.tgt_len, 10);
        assert_eq!(cigar_mapper.ref_pos.len(), cigar_mapper.tgt_pos.len());
        assert_eq!(
            cigar_mapper.ref_pos,
            vec![0, 3, 5, 6, 7, 8, 11, 12, 13, 14, 14, 15]
        );
        assert_eq!(
            cigar_mapper.tgt_pos,
            vec![0, 3, 3, 4, 5, 6, 6, 7, 7, 8, 9, 10]
        );

        // ref to tgt
        {
            let cases = vec![
                (0, "start", 0, 0, CigarOp::Eq),
                (0, "end", 0, 0, CigarOp::Eq),
                (1, "start", 1, 0, CigarOp::Eq),
                (1, "end", 1, 0, CigarOp::Eq),
                (2, "start", 2, 0, CigarOp::Eq),
                (2, "end", 2, 0, CigarOp::Eq),
                (3, "start", 2, 1, CigarOp::Skip),
                (3, "end", 2, 1, CigarOp::Skip),
                (4, "start", 3, -1, CigarOp::Skip),
                (4, "end", 3, -1, CigarOp::Skip),
                (5, "start", 3, 0, CigarOp::Eq),
                (5, "end", 3, 0, CigarOp::Eq),
                (6, "start", 4, 0, CigarOp::Mismatch),
                (6, "end", 4, 0, CigarOp::Mismatch),
                (7, "start", 5, 0, CigarOp::Eq),
                (7, "end", 5, 0, CigarOp::Eq),
                (8, "start", 5, 1, CigarOp::Skip),
                (8, "end", 5, 1, CigarOp::Skip),
                (9, "start", 5, 2, CigarOp::Skip),
                (9, "end", 5, 2, CigarOp::Skip),
                (10, "start", 6, -1, CigarOp::Skip),
                (10, "end", 6, -1, CigarOp::Skip),
                (11, "start", 6, 0, CigarOp::Eq),
                (11, "end", 6, 0, CigarOp::Eq),
                (12, "start", 6, 0, CigarOp::Ins),
                (12, "end", 7, 0, CigarOp::Ins),
                (13, "start", 7, 0, CigarOp::Eq),
                (13, "end", 7, 0, CigarOp::Eq),
                (14, "start", 9, 0, CigarOp::Eq),
                (14, "end", 9, 0, CigarOp::Eq),
            ];
            for (arg_pos, arg_end, pos, offset, cigar_op) in cases {
                assert_eq!(
                    cigar_mapper.map_ref_to_tgt(arg_pos, arg_end, true)?,
                    CigarMapperResult {
                        pos,
                        offset,
                        cigar_op
                    },
                    "case = {:?}",
                    (arg_pos, arg_end, pos, offset, cigar_op)
                );
            }
        }

        // tgt to ref
        {
            let cases = vec![
                (0, "start", 0, 0, CigarOp::Eq),
                (0, "end", 0, 0, CigarOp::Eq),
                (1, "start", 1, 0, CigarOp::Eq),
                (1, "end", 1, 0, CigarOp::Eq),
                (2, "start", 2, 0, CigarOp::Eq),
                (2, "end", 2, 0, CigarOp::Eq),
                (3, "start", 5, 0, CigarOp::Eq),
                (3, "end", 5, 0, CigarOp::Eq),
                (4, "start", 6, 0, CigarOp::Mismatch),
                (4, "end", 6, 0, CigarOp::Mismatch),
                (5, "start", 7, 0, CigarOp::Eq),
                (5, "end", 7, 0, CigarOp::Eq),
                (6, "start", 11, 0, CigarOp::Eq),
                (6, "end", 11, 0, CigarOp::Eq),
                (7, "start", 13, 0, CigarOp::Eq),
                (7, "end", 13, 0, CigarOp::Eq),
                (8, "start", 13, 0, CigarOp::Del),
                (8, "end", 14, 0, CigarOp::Del),
                (9, "start", 14, 0, CigarOp::Eq),
                (9, "end", 14, 0, CigarOp::Eq),
            ];
            for (arg_pos, arg_end, pos, offset, cigar_op) in cases {
                assert_eq!(
                    cigar_mapper.map_tgt_to_ref(arg_pos, arg_end, true)?,
                    CigarMapperResult {
                        pos,
                        offset,
                        cigar_op
                    },
                    "case = {:?}",
                    (arg_pos, arg_end, pos, offset, cigar_op)
                );
            }
        }

        Ok(())
    }

    #[test]
    fn cigar_mapper_strict_bounds() -> Result<(), anyhow::Error> {
        // 0   1   2           3   4   5               6       7   8   9  tgt
        // =   =   =   N   N   =   X   =   N   N   N   =   I   =   D   =
        // 0   1   2   3   4   5   6   7   8   9  10  11  12  13      14  ref
        let cigar = "3=2N=X=3N=I=D=".to_string();
        let cigar_str = parse_cigar_string(&cigar)?;
        let cigar_mapper = CigarMapper::new(&cigar_str);

        // error for out of bounds on left?
        assert!(cigar_mapper.map_ref_to_tgt(-1, "start", true).is_err());
        // ... and right?
        assert!(cigar_mapper
            .map_ref_to_tgt(cigar_mapper.ref_len + 1, "start", true)
            .is_err());

        // test whether 1 base outside bounds results in correct position
        assert_eq!(
            cigar_mapper.map_ref_to_tgt(0, "start", true)?,
            CigarMapperResult {
                pos: 0,
                offset: 0,
                cigar_op: CigarOp::Eq,
            }
        );
        assert_eq!(
            cigar_mapper.map_ref_to_tgt(-1, "start", false)?,
            CigarMapperResult {
                pos: -1,
                offset: 0,
                cigar_op: CigarOp::Eq,
            }
        );
        assert_eq!(
            cigar_mapper.map_ref_to_tgt(cigar_mapper.ref_len, "start", true)?,
            CigarMapperResult {
                pos: cigar_mapper.tgt_len,
                offset: 0,
                cigar_op: CigarOp::Eq,
            }
        );
        assert_eq!(
            cigar_mapper.map_ref_to_tgt(cigar_mapper.ref_len - 1, "start", false)?,
            CigarMapperResult {
                pos: cigar_mapper.tgt_len - 1,
                offset: 0,
                cigar_op: CigarOp::Eq,
            }
        );

        Ok(())
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
