//! Mapping positions between pairs of sequence alignments.
//!
//! `AlignmentMapper` is at the heart of mapping between aligned sequences.

// Implementation note re: "no-zero correction": HGVS doesn't have a
// 0. Counting is -3, -2, -1, 1, 2, 3 :-/ Coordinate calculations must
// take this discontinuity in c. positions into account.  The
// implementation of imaginary transcript positions creates a second
// discontinuity. (By analogy with c., n.0 is declared to not exist.)
// The strategy used in this code is to use internal c0 and n0
// coordinates, which include 0, for coordinate calculations and to
// translate these to c. and n. positions as needed.
//
//                imag.                                                 imag.
//              upstream     5' UTR           CDS          3' UTR      downstr
//                                     |>
//            - - - - - - ———————————— ||||||||||||||||| ——————————— - - - - - -
//                           a     b     C     D     E     f     g     h     i
//    c.        -4    -3    -2    -1  !  1     2     3  ! *1    *2    *3    *4
//    c0        -4    -3    -2    -1     0     1     2     3     4     5     6
//    n0        -2    -1     0     1     2     3     4     5     6     7     8
//    n.        -2    -1  !  1     2     3     4     5     6     7     8     9
//    g.   ... 123   124   125   126   127   128   129   130   131   132   133 ...

use std::{fmt::Display, rc::Rc};

use nom::{combinator::all_consuming, multi::many0};

use crate::{
    data::interface::{Provider, TxExonsRecord},
    parser::{GenomeInterval, Mu, TxInterval, TxPos},
};

/// Convert zero-based coordinate to hgvs (1 based, missing zero)
fn zbc_to_hgvs(i: i32) -> i32 {
    if i >= 0 {
        i + 1
    } else {
        i
    }
}

/// Convert hgvs (1 based, missing zero)
fn hgvs_to_zbc(i: i32) -> i32 {
    if i >= 1 {
        i - 1
    } else {
        i
    }
}

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
        match self {
            CigarOp::Eq | CigarOp::Match | CigarOp::Mismatch | CigarOp::Ins | CigarOp::Skip => true,
            _ => false,
        }
    }

    pub fn is_advance_tgt(&self) -> bool {
        match self {
            CigarOp::Eq | CigarOp::Match | CigarOp::Mismatch | CigarOp::Del => true,
            _ => false,
        }
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

impl Into<char> for CigarOp {
    fn into(self) -> char {
        match self {
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
        write!(f, "{}", std::convert::Into::<char>::into(self.clone()))
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
            count: if count.len() == 0 {
                1
            } else {
                str::parse(count).unwrap()
            },
            op: op.chars().nth(0).unwrap().try_into().unwrap(),
        }
    }
}

#[derive(Debug, PartialEq, PartialOrd, Default, Clone)]
pub struct CigarString{
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
        Ok(for item in &self.elems {
                write!(f, "{}", &item)?
            }
        )
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
    Ok(CigarString::from(all_consuming(many0(parse::cigar_element))(input)
        .map_err(|e| anyhow::anyhow!("Problem with parsing: {:?}", e))?
        .1))
}

/// Builds a single CIGAR string representing an alignment of the transcript sequence to a
/// reference sequence, including introns.
///
/// The input exons are expected to be in transcript order, and the resulting CIGAR is also
/// in transcript order.
pub fn build_tx_cigar(
    exons: &Vec<TxExonsRecord>,
    strand: i16,
) -> Result<CigarString, anyhow::Error> {
    if exons.is_empty() {
        return Err(anyhow::anyhow!(
            "Cannot build CIGAR string from empty exons"
        ));
    }

    // Parse CIGAR string and flip if on reverse strand.
    let exon_cigars: Result<Vec<CigarString>, anyhow::Error> = exons
        .iter()
        .map(|record| {
            let mut cigar = parse_cigar_string(&record.cigar)?;
            if strand == -1 {
                cigar.reverse();
            }
            Ok(cigar)
        })
        .collect();
    let exon_cigars = exon_cigars?;

    let mut result = exon_cigars[0].clone(); // exon 1
    for i in 1..exon_cigars.len() {
        result.push(CigarElement {
            count: exons[i].alt_start_i - exons[i - 1].alt_end_i,
            op: CigarOp::Skip,
        });
        result.append(&mut exon_cigars[i].clone());
    }
    Ok(result)
}

/// Helper function that wraps a value into `Option` but returns `None` for the default value.
pub fn none_if_default<T>(value: T) -> Option<T>
where
    T: Default + PartialEq,
{
    if value == T::default() {
        None
    } else {
        Some(value)
    }
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
#[derive(Default)]
struct CigarMapper {
    pub cigar_string: CigarString,
    pub ref_pos: Vec<i32>,
    pub tgt_pos: Vec<i32>,
    pub cigar_op: Vec<CigarOp>,
    pub ref_len: i32,
    pub tgt_len: i32,
}

#[derive(Debug, PartialEq)]
struct CigarMapperResult {
    pub pos: i32,
    pub offset: i32,
    pub cigar_op: CigarOp,
}

impl CigarMapper {
    fn new(cigar_string: &CigarString) -> Self {
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
                "Position is beyond the bounds of transcript record"
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
            if pos - from_pos[pos_i] + 1 <= from_pos[pos_i + 1] - pos {
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

/// Configuration for mapping.
struct Config {
    /// Require transcript variants to be within transcript sequence bounds.
    pub strict_bounds: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            strict_bounds: true,
        }
    }
}

/// Map HGVS location objects between genomic (g), non-coding (n) and cds (c)
/// coordinates according to a CIGAR string.
struct AlignmentMapper {
    /// Configuration for alignment mapping.
    pub config: Config,
    /// Data provider to use for the mapping.
    pub provider: Rc<dyn Provider>,

    /// The transcript accession.
    pub tx_ac: String,
    /// The reference sequence asccession.
    pub alt_ac: String,
    /// The alignment method.
    pub alt_aln_method: String,
    pub strand: i16,
    pub gc_offset: i32,
    pub cds_start_i: Option<i32>,
    pub cds_end_i: Option<i32>,
    pub tgt_len: i32,
    pub cigar_mapper: CigarMapper,
}

impl AlignmentMapper {
    pub fn new(
        provider: Rc<dyn Provider>,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<AlignmentMapper, anyhow::Error> {
        let (strand, gc_offset, cds_start_i, cds_end_i, tgt_len, cigar_mapper) =
            if alt_aln_method != "transcript" {
                let tx_info = provider.get_tx_info(tx_ac, alt_ac, alt_aln_method)?;
                let tx_exons = {
                    let mut tx_exons = provider.get_tx_exons(tx_ac, alt_ac, alt_aln_method)?;
                    if tx_exons.is_empty() {
                        return Err(anyhow::anyhow!(
                            "Found no exons for tx_ac={}, alt_ac={}, alt_aln_method={}",
                            tx_ac,
                            alt_ac,
                            alt_aln_method
                        ));
                    }

                    // Issue biocommons/hgvs#386: An assumption when building the CIGAR string is that
                    // exons are adjacent. Assert that here.
                    tx_exons.sort_by(|a, b| a.ord.partial_cmp(&b.ord).unwrap());
                    let offenders = tx_exons
                        .windows(2)
                        .filter(|pair| {
                            let lhs = &pair[0];
                            let rhs = &pair[1];
                            lhs.tx_end_i != rhs.tx_start_i
                        })
                        .collect::<Vec<_>>();
                    if !offenders.is_empty() {
                        return Err(anyhow::anyhow!(
                            "Non-adjacent exons for tx_acc={}, alt_acc={}, alt_aln_method={}: {:?}",
                            tx_ac,
                            alt_ac,
                            alt_aln_method,
                            &offenders
                        ));
                    }

                    tx_exons
                };

                let strand = tx_exons[0].alt_strand;
                let gc_offset = tx_exons[0].alt_start_i;
                let cds_start_i = tx_info.cds_start_i;
                let cds_end_i = tx_info.cds_end_i;

                let cigar_mapper = CigarMapper::new(&build_tx_cigar(&tx_exons, strand)?);
                let tgt_len = cigar_mapper.tgt_len;

                (
                    strand,
                    gc_offset,
                    cds_start_i,
                    cds_end_i,
                    tgt_len,
                    cigar_mapper,
                )
            } else {
                // this covers the identity cases n <-> c
                let tx_identity_info = provider.get_tx_identity_info(tx_ac)?;

                let cds_start_i = tx_identity_info.cds_start_i;
                let cds_end_i = tx_identity_info.cds_end_i;
                let tgt_len: i32 = tx_identity_info.lengths.iter().sum();

                (
                    1, // strand
                    0, // gc_offset
                    Some(cds_start_i),
                    Some(cds_end_i),
                    tgt_len,
                    Default::default(),
                )
            };

        if cds_start_i.is_none() != cds_end_i.is_none() {
            return Err(anyhow::anyhow!(
                "CDs start and end position but be both or neither defined."
            ));
        }

        Ok(AlignmentMapper {
            config: Default::default(),
            provider: provider,
            tx_ac: tx_ac.to_string(),
            alt_ac: alt_ac.to_string(),
            alt_aln_method: alt_aln_method.to_string(),
            cds_start_i,
            cds_end_i,
            tgt_len,
            cigar_mapper,
            strand: strand,
            gc_offset: gc_offset,
        })
    }

    /// Convert a genomic (g.) interval to a transcript cDNA (n.) interval.
    fn g_to_n(&self, g_interval: &GenomeInterval) -> Result<Mu<TxInterval>, anyhow::Error> {
        if let GenomeInterval {
            begin: Some(begin),
            end: Some(end),
        } = g_interval
        {
            let grs = begin - 1 - self.gc_offset;
            let gre = end - 1 - self.gc_offset;

            // Compute start/end on forward strand with respect to the genome.
            //
            // frs, fre = (f)orward (r)na (s)tart & (e)nd
            let frs = self
                .cigar_mapper
                .map_ref_to_tgt(grs, "start", self.config.strict_bounds)?;
            let fre = self
                .cigar_mapper
                .map_ref_to_tgt(gre, "end", self.config.strict_bounds)?;

            // Project to reverse strand if necessary.
            let (frs, fre) = if self.strand == -1 {
                (
                    CigarMapperResult {
                        pos: self.tgt_len - 1 - fre.pos,
                        offset: -fre.offset,
                        cigar_op: fre.cigar_op,
                    },
                    CigarMapperResult {
                        pos: self.tgt_len - 1 - frs.pos,
                        offset: -frs.offset,
                        cigar_op: frs.cigar_op,
                    },
                )
            } else {
                (frs, fre)
            };

            // The position is uncertain if the alignment ends in a gap.
            Ok(Mu::from(
                TxInterval {
                    begin: TxPos {
                        base: zbc_to_hgvs(frs.pos),
                        offset: none_if_default(frs.offset),
                    },
                    end: TxPos {
                        base: zbc_to_hgvs(fre.pos),
                        offset: none_if_default(fre.offset),
                    },
                },
                frs.cigar_op != CigarOp::Del
                    && frs.cigar_op != CigarOp::Ins
                    && fre.cigar_op != CigarOp::Del
                    && fre.cigar_op != CigarOp::Ins,
            ))
        } else {
            Err(anyhow::anyhow!(
                "Cannot project genome interval with missing start or end position: {}",
                g_interval
            ))
        }
    }
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use crate::{
        data::interface::TxExonsRecord,
        mapper::alignment::{parse_cigar_string, CigarElement, CigarMapperResult, CigarOp},
    };

    use super::{build_tx_cigar, none_if_default, CigarMapper};

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
    fn build_tx_cigar_empty() {
        assert!(!build_tx_cigar(&Vec::new(), 1).is_ok());
    }

    #[test]
    fn build_tx_cigar_forward() -> Result<(), anyhow::Error> {
        let exons = vec![
            TxExonsRecord {
                tx_start_i: 0,
                tx_end_i: 10,
                alt_start_i: 100,
                alt_end_i: 110,
                cigar: "5M1I4M".to_string(),
                ..Default::default()
            },
            TxExonsRecord {
                tx_start_i: 10,
                tx_end_i: 21,
                alt_start_i: 120,
                alt_end_i: 131,
                cigar: "7M1I2M".to_string(),
                ..Default::default()
            },
        ];

        assert_eq!(
            format!("{}", &build_tx_cigar(&exons, 1)?),
            "5MI4M10N7MI2M"
        );

        Ok(())
    }

    #[test]
    fn build_tx_cigar_reverse() -> Result<(), anyhow::Error> {
        let exons = vec![
            TxExonsRecord {
                tx_start_i: 0,
                tx_end_i: 10,
                alt_start_i: 100,
                alt_end_i: 110,
                cigar: "5M1I4M".to_string(),
                ..Default::default()
            },
            TxExonsRecord {
                tx_start_i: 10,
                tx_end_i: 21,
                alt_start_i: 120,
                alt_end_i: 131,
                cigar: "7M1I2M".to_string(),
                ..Default::default()
            },
        ];

        assert_eq!(
            format!("{}", &build_tx_cigar(&exons, -1)?),
            "4MI5M10N2MI7M"
        );

        Ok(())
    }

    #[test]
    fn run_none_if_default() {
        assert_eq!(none_if_default(0u32), None);
        assert_eq!(none_if_default(1u32), Some(1u32));
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
                    "case = {:?}", (arg_pos, arg_end, pos, offset, cigar_op)

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
                    "case = {:?}", (arg_pos, arg_end, pos, offset, cigar_op)
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
