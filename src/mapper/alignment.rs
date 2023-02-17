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

use std::rc::Rc;

use regex::Regex;

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

lazy_static::lazy_static! {
    static ref RE_CIGAR_ELEM: Regex = Regex::new(r"\d+[=DIMNX]").unwrap();
}

/// Builds a single CIGAR string representing an alignment of the transcript sequence to a
/// reference sequence, including introns.
///
/// The input exons are expected to be in transcript order, and the resulting CIGAR is also
/// in transcript order.
pub fn build_tx_cigar(exons: &Vec<TxExonsRecord>, strand: i16) -> Result<String, anyhow::Error> {
    fn reverse_cigar(cigar: &str) -> String {
        let mut v = RE_CIGAR_ELEM.find_iter(cigar).map(|m| m.as_str()).collect::<Vec<_>>();
        v.reverse();
        v.join("")
    }

    if exons.is_empty() {
        return Err(anyhow::anyhow!("Cannot build CIGAR string from empty exons"));
    }

    // Flip orientation of all CIGARs if on reverse strand.
    let exon_cigars = if strand == -1 {
        exons.iter().map(|record| reverse_cigar(&record.cigar)).collect::<Vec<_>>()
    } else {
        exons.iter().map(|record| record.cigar.clone()).collect::<Vec<_>>()
    };

    let mut tx_cigar = vec![exon_cigars[0].clone()];  // exon 1
    for i in 1..exon_cigars.len() {
        let intron = format!("{}N", exons[i].alt_start_i - exons[i - 1].alt_end_i);
        tx_cigar.push(intron);
        tx_cigar.push(exon_cigars[i].clone());
    }

    Ok(tx_cigar.join(""))
}

#[derive(Default)]
struct CigarMapper {
    cigar: String,
}

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

struct CigarMapperResult {
    pub pos: i32,
    pub offset: i32,
    pub cigar_op: char,
}

impl CigarMapper {
    fn new(cigar: String) -> Self {
        Self { cigar }
    }

    pub fn tgt_len(&self) -> i32 {
        todo!()
    }

    pub fn map_ref_to_tgt(
        &self,
        pos: i32,
        end: &str,
        strict_bounds: bool,
    ) -> Result<CigarMapperResult, anyhow::Error> {
        todo!()
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

                let cigar_mapper = CigarMapper::new(build_tx_cigar(&tx_exons, strand)?);
                let tgt_len = cigar_mapper.tgt_len();

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
                !"DI".contains(frs.cigar_op) && !"DI".contains(fre.cigar_op),
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

    use super::{none_if_default, build_tx_cigar};

    #[test]
    fn build_tx_cigar_empty() {
        assert!(!build_tx_cigar(&Vec::new(), 1).is_ok());
    }

    #[test]
    fn run_none_if_default() {
        assert_eq!(none_if_default(0u32), None);
        assert_eq!(none_if_default(1u32), Some(1u32));
    }
}
