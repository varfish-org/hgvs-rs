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

use crate::{
    data::interface::{Provider, TxExonsRecord},
    parser::{CdsInterval, CdsPos, GenomeInterval, Mu, TxInterval, TxPos},
};

use super::cigar::{
    parse_cigar_string, CigarElement, CigarMapper, CigarMapperResult, CigarOp, CigarString,
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

/// Builds a single CIGAR string representing an alignment of the transcript sequence to a
/// reference sequence, including introns.
///
/// The input exons are expected to be in transcript order, and the resulting CIGAR is also
/// in transcript order.
pub fn build_tx_cigar(exons: &Vec<TxExonsRecord>, strand: i16) -> Result<CigarString, Error> {
    if exons.is_empty() {
        return Err(anyhow::anyhow!(
            "Cannot build CIGAR string from empty exons"
        ));
    }

    // Parse CIGAR string and flip if on reverse strand.
    let exon_cigars: Result<Vec<CigarString>, Error> = exons
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

/// Configuration for mapping.
#[derive(Debug, Clone)]
pub struct Config {
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
pub struct Mapper {
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

impl Mapper {
    pub fn new(
        config: &Config,
        provider: Rc<dyn Provider>,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Mapper, Error> {
        let config = config.clone();
        let (strand, gc_offset, cds_start_i, cds_end_i, tgt_len, cigar_mapper) = if alt_aln_method
            != "transcript"
        {
            let tx_info = provider.get_tx_info(tx_ac, alt_ac, alt_aln_method)?;
            let tx_exons = {
                let tx_exons = provider.get_tx_exons(tx_ac, alt_ac, alt_aln_method)?;
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
                let mut sorted_exons = tx_exons.clone();
                sorted_exons
                    .sort_by(|a, b| a.ord.partial_cmp(&b.ord).expect("comparison failed / NaN?"));
                let offenders = sorted_exons
                    .windows(2)
                    .filter(|pair| {
                        let lhs = &pair[0];
                        let rhs = &pair[1];
                        lhs.tx_end_i != rhs.tx_start_i
                    })
                    .collect::<Vec<_>>();
                if !offenders.is_empty() {
                    return Err(anyhow::anyhow!(
                        "Non-adjacent exons for tx_acc={}, alt_acc={}, alt_aln_method={}: {:#?}",
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

            if cds_start_i.is_none() != cds_end_i.is_none() {
                return Err(anyhow::anyhow!(
                    "CDS start and end must both be defined or undefined"
                ));
            }

            let cigar = build_tx_cigar(&tx_exons, strand)?;
            let cigar_mapper = CigarMapper::new(&cigar);
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

        Ok(Mapper {
            config,
            provider,
            tx_ac: tx_ac.to_string(),
            alt_ac: alt_ac.to_string(),
            alt_aln_method: alt_aln_method.to_string(),
            cds_start_i,
            cds_end_i,
            tgt_len,
            cigar_mapper,
            strand,
            gc_offset,
        })
    }

    /// Convert a genomic (g.) interval to a transcript (n.) interval.
    pub fn g_to_n(&self, g_interval: &GenomeInterval) -> Result<Mu<TxInterval>, Error> {
        if let GenomeInterval {
            start: Some(begin),
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
            let n_interval = TxInterval {
                start: TxPos {
                    base: zbc_to_hgvs(frs.pos),
                    offset: none_if_default(frs.offset),
                },
                end: TxPos {
                    base: zbc_to_hgvs(fre.pos),
                    offset: none_if_default(fre.offset),
                },
            };
            let result = if frs.cigar_op != CigarOp::Del
                && frs.cigar_op != CigarOp::Ins
                && fre.cigar_op != CigarOp::Del
                && fre.cigar_op != CigarOp::Ins
            {
                Mu::Certain(n_interval)
            } else {
                Mu::Uncertain(n_interval)
            };
            Ok(result)
        } else {
            Err(anyhow::anyhow!(
                "Cannot project genome interval with missing start or end position: {}",
                g_interval
            ))
        }
    }

    /// Convert a transcript (n.) interval to a genomic (g.) interval.
    pub fn n_to_g(&self, n_interval: &TxInterval) -> Result<Mu<GenomeInterval>, Error> {
        let frs = hgvs_to_zbc(n_interval.start.base);
        let start_offset = n_interval.start.offset.unwrap_or(0);
        let fre = hgvs_to_zbc(n_interval.end.base);
        let end_offset = n_interval.end.offset.unwrap_or(0);

        let (fre, frs, start_offset, end_offset) = if self.strand == -1 {
            (
                self.tgt_len - 1 - frs,
                self.tgt_len - 1 - fre,
                -end_offset,
                -start_offset,
            )
        } else {
            (fre, frs, start_offset, end_offset)
        };

        // Obtain the genomic range start (grs) and end (gre).
        let grs = self
            .cigar_mapper
            .map_tgt_to_ref(frs, "start", self.config.strict_bounds)?;
        let gre = self
            .cigar_mapper
            .map_tgt_to_ref(fre, "end", self.config.strict_bounds)?;
        let (grs_pos, gre_pos) = (grs.pos + self.gc_offset + 1, gre.pos + self.gc_offset + 1);
        let (gs, ge) = (grs_pos + start_offset, gre_pos + end_offset);

        // The returned interval would be uncertain when locating at alignment gaps.
        Ok(Mu::from(
            GenomeInterval {
                start: Some(gs),
                end: Some(ge),
            },
            grs.cigar_op != CigarOp::Del
                && grs.cigar_op != CigarOp::Ins
                && gre.cigar_op != CigarOp::Del
                && gre.cigar_op != CigarOp::Ins,
        ))
    }

    fn pos_n_to_c(&self, pos: &TxPos) -> CdsPos {
        let cds_start_i = self
            .cds_start_i
            .expect("cannot convert n. to c. without CDS positions");
        let cds_end_i = self
            .cds_end_i
            .expect("cannot convert n. to c. without CDS positions");
        if pos.base <= cds_start_i {
            CdsPos {
                base: pos.base - cds_start_i - (if pos.base > 0 { 1 } else { 0 }),
                offset: pos.offset,
                cds_from: crate::parser::CdsFrom::Start,
            }
        } else if pos.base > cds_start_i && pos.base <= cds_end_i {
            CdsPos {
                base: pos.base - cds_start_i,
                offset: pos.offset,
                cds_from: crate::parser::CdsFrom::Start,
            }
        } else {
            CdsPos {
                base: pos.base - cds_end_i,
                offset: pos.offset,
                cds_from: crate::parser::CdsFrom::End,
            }
        }
    }

    /// Convert a transcript (n.) interval to a CDS (c.) interval.
    pub fn n_to_c(&self, n_interval: &TxInterval) -> Result<CdsInterval, Error> {
        if self.cds_start_i.is_none() {
            return Err(anyhow::anyhow!(
                "CDS is undefined for {}; cannot map to c. coordinate (non-coding transcript?)",
                self.tx_ac
            ));
        }

        if self.config.strict_bounds
            && (n_interval.start.base <= 0 || n_interval.end.base > self.tgt_len)
        {
            return Err(anyhow::anyhow!(
                "The given coordinate is outside the bounds of the reference sequence."
            ));
        }

        Ok(CdsInterval {
            start: self.pos_n_to_c(&n_interval.start),
            end: self.pos_n_to_c(&n_interval.end),
        })
    }

    pub fn pos_c_to_n(&self, pos: &CdsPos) -> Result<TxPos, Error> {
        let cds_start_i = self
            .cds_start_i
            .expect("cannot convert c. to n. without CDS positions");
        let cds_end_i = self
            .cds_end_i
            .expect("cannot convert n. to c. without CDS positions");

        let n = match pos.cds_from {
            crate::parser::CdsFrom::Start => {
                let n = pos.base + cds_start_i;
                if pos.base < 0 {
                    // correct for lack of c.0 coordinate
                    n + 1
                } else {
                    n
                }
            }
            crate::parser::CdsFrom::End => pos.base + cds_end_i,
        };

        let n = if n <= 0 {
            // correct for lack of n.0 coordinate
            n - 1
        } else {
            n
        };

        if (n <= 0 || n > self.tgt_len) && self.config.strict_bounds {
            Err(anyhow::anyhow!("c.{:?} coordinate is out of boounds", pos))
        } else {
            Ok(TxPos {
                base: n,
                offset: pos.offset,
            })
        }
    }

    /// Convert a a CDS (c.) interval to a transcript (n.) interval.
    pub fn c_to_n(&self, c_interval: &CdsInterval) -> Result<TxInterval, Error> {
        if self.cds_start_i.is_none() {
            return Err(anyhow::anyhow!(
                "CDS is undefined for {}; cannot map to c. coordinate (non-coding transcript?)",
                self.tx_ac
            ));
        }

        let n_start = self.pos_c_to_n(&c_interval.start);
        let n_end = self.pos_c_to_n(&c_interval.end);
        let result = TxInterval {
            start: n_start?,
            end: n_end?,
        };
        Ok(result)
    }

    /// Convert a genomic (g.) interval to a CDS (c.) interval.
    pub fn g_to_c(&self, g_interval: &GenomeInterval) -> Result<Mu<CdsInterval>, Error> {
        let n_interval = self.g_to_n(g_interval)?;
        Ok(match &n_interval {
            Mu::Certain(n_interval) => Mu::Certain(self.n_to_c(n_interval)?),
            Mu::Uncertain(n_interval) => Mu::Uncertain(self.n_to_c(n_interval)?),
        })
    }

    /// Convert a CDS (c.) interval to a genomic (g.) interval.
    pub fn c_to_g(&self, c_interval: &CdsInterval) -> Result<Mu<GenomeInterval>, Error> {
        let n_interval = self.c_to_n(c_interval)?;
        let g_interval = self.n_to_g(&n_interval)?;
        Ok(g_interval)
    }

    /// Return the transcript is coding.
    pub fn is_coding_transcript(&self) -> bool {
        self.cds_start_i.is_some()
    }

    /// Return whether the given genome interval is in bounds.
    pub fn is_g_interval_in_bounds(&self, g_interval: &GenomeInterval) -> bool {
        let grs = g_interval
            .start
            .expect("cannot check g. interval for non-concrete positions")
            - 1
            - self.gc_offset;
        let gre = g_interval
            .end
            .expect("cannot check g. interval for non-concrete positions")
            - 1
            - self.gc_offset;
        grs >= 0 && gre <= self.cigar_mapper.ref_len
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;

    use crate::{
        data::{interface::TxExonsRecord, uta_sr::test_helpers::build_provider},
        parser::{CdsFrom, CdsInterval, CdsPos, GenomeInterval, Mu, TxInterval, TxPos},
    };

    use super::{build_tx_cigar, none_if_default, Mapper};

    #[test]
    fn build_tx_cigar_empty() {
        assert!(build_tx_cigar(&Vec::new(), 1).is_err());
    }

    #[test]
    fn build_tx_cigar_forward() -> Result<(), Error> {
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

        assert_eq!(format!("{}", &build_tx_cigar(&exons, 1)?), "5MI4M10N7MI2M");

        Ok(())
    }

    #[test]
    fn build_tx_cigar_reverse() -> Result<(), Error> {
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

        assert_eq!(format!("{}", &build_tx_cigar(&exons, -1)?), "4MI5M10N2MI7M");

        Ok(())
    }

    #[test]
    fn run_none_if_default() {
        assert_eq!(none_if_default(0u32), None);
        assert_eq!(none_if_default(1u32), Some(1u32));
        assert_eq!(none_if_default(1i32), Some(1i32));
        assert_eq!(none_if_default(-1i32), Some(-1i32));
    }

    #[test]
    fn construction() -> Result<(), Error> {
        let provider = build_provider()?;

        assert_eq!(provider.data_version(), "uta_20210129");
        assert_eq!(provider.schema_version(), "1.1");

        Ok(())
    }

    #[test]
    fn failures() -> Result<(), Error> {
        let provider = build_provider()?;

        // unknown sequences

        assert!(Mapper::new(
            &Default::default(),
            provider.clone(),
            "bogus",
            "NM_033089.6",
            "splign"
        )
        .is_err());
        assert!(Mapper::new(
            &Default::default(),
            provider.clone(),
            "bogus",
            "NM_033089.6",
            "transcript"
        )
        .is_err());
        assert!(Mapper::new(
            &Default::default(),
            provider.clone(),
            "NM_033089.6",
            "bogus",
            "splign"
        )
        .is_err());
        assert!(Mapper::new(
            &Default::default(),
            provider.clone(),
            "NM_000051.3",
            "NC_000011.9",
            "bogus"
        )
        .is_err());

        // invalid intervals

        {
            let am = Mapper::new(
                &Default::default(),
                provider.clone(),
                "NM_000348.3",
                "NC_000002.11",
                "splign",
            )?;
            assert!(am
                .n_to_c(&TxInterval {
                    start: TxPos {
                        base: -1,
                        offset: None
                    },
                    end: TxPos {
                        base: -1,
                        offset: None
                    },
                },)
                .is_err());
        }

        {
            let am = Mapper::new(
                &Default::default(),
                provider,
                "NM_000348.3",
                "NC_000002.11",
                "splign",
            )?;
            assert!(am
                .c_to_n(&CdsInterval {
                    start: CdsPos {
                        base: 99999,
                        offset: None,
                        cds_from: CdsFrom::Start
                    },
                    end: CdsPos {
                        base: 99999,
                        offset: None,
                        cds_from: CdsFrom::Start
                    },
                },)
                .is_err());
        }

        Ok(())
    }

    /// Helper for running multiple projection cases.
    fn run_test_cases(
        tx_ac: &str,
        alt_ac: &str,
        cases: &Vec<(GenomeInterval, TxInterval, CdsInterval)>,
    ) -> Result<(), Error> {
        let provider = build_provider()?;
        let mapper = Mapper::new(&Default::default(), provider, tx_ac, alt_ac, "splign")?;

        for (g_interval, n_interval, c_interval) in cases {
            assert_eq!(
                c_interval,
                mapper.g_to_c(g_interval)?.inner(),
                "{}~{} {} mapper.g_to_c",
                tx_ac,
                alt_ac,
                &g_interval
            );
            assert_eq!(
                c_interval,
                &mapper.n_to_c(n_interval)?,
                "{}~{} {} mapper.n_to_c",
                tx_ac,
                alt_ac,
                &n_interval
            );

            assert_eq!(
                g_interval,
                mapper.c_to_g(c_interval)?.inner(),
                "{}~{} {} mapper.c_to_g",
                tx_ac,
                alt_ac,
                &c_interval
            );
            assert_eq!(
                Mu::Certain(g_interval.clone()),
                mapper.n_to_g(n_interval)?,
                "{}~{} {} mapper.n_to_g",
                tx_ac,
                alt_ac,
                &g_interval
            );

            assert_eq!(
                n_interval,
                &mapper.c_to_n(c_interval)?,
                "{}~{} {} mapper.c_to_n",
                tx_ac,
                alt_ac,
                &c_interval
            );
            assert_eq!(
                Mu::Certain(n_interval.clone()),
                mapper.g_to_n(g_interval)?,
                "{}~{} {} mapper.g_to_n",
                tx_ac,
                alt_ac,
                &n_interval
            );
        }

        Ok(())
    }

    /// Use NM_178434.2 tests to test mapping with uncertain positions
    // #[test]  // not yet supported (same as Python package)
    fn _test_lce3c_uncertain() -> Result<(), Error> {
        let tx_ac = "NM_178434.2";
        let alt_ac = "NC_000001.10";
        let test_cases = vec![
            (
                GenomeInterval::from_str("?_152573139")?,
                TxInterval::from_str("?_2")?,
                CdsInterval::from_str("?_-69")?,
            ),
            (
                GenomeInterval::from_str("152573138_?")?,
                TxInterval::from_str("1_?")?,
                CdsInterval::from_str("-70_?")?,
            ),
        ];

        run_test_cases(tx_ac, alt_ac, &test_cases)?;

        Ok(())
    }

    /// NM_178434.2: LCE3C single exon, strand = +1, all coordinate input/output are in HGVS
    #[test]
    fn test_lce3c() -> Result<(), Error> {
        let tx_ac = "NM_178434.2";
        let alt_ac = "NC_000001.10";
        let test_cases = vec![
            // 5'
            (
                GenomeInterval::from_str("152573138")?,
                TxInterval::from_str("1")?,
                CdsInterval::from_str("-70")?,
            ),
            (
                GenomeInterval::from_str("152573140")?,
                TxInterval::from_str("3")?,
                CdsInterval::from_str("-68")?,
            ),
            // CDS
            (
                GenomeInterval::from_str("152573207")?,
                TxInterval::from_str("70")?,
                CdsInterval::from_str("-1")?,
            ),
            (
                GenomeInterval::from_str("152573208")?,
                TxInterval::from_str("71")?,
                CdsInterval::from_str("1")?,
            ),
            // 3'
            (
                GenomeInterval::from_str("152573492")?,
                TxInterval::from_str("355")?,
                CdsInterval::from_str("285")?,
            ),
            (
                GenomeInterval::from_str("152573493")?,
                TxInterval::from_str("356")?,
                CdsInterval::from_str("*1")?,
            ),
            (
                GenomeInterval::from_str("152573560")?,
                TxInterval::from_str("423")?,
                CdsInterval::from_str("*68")?,
            ),
            (
                GenomeInterval::from_str("152573562")?,
                TxInterval::from_str("425")?,
                CdsInterval::from_str("*70")?,
            ),
        ];

        run_test_cases(tx_ac, alt_ac, &test_cases)?;

        Ok(())
    }

    /// NM_033445.2: H2AW single exon, strand = -1, all coordinate input/output are in HGVS
    #[test]
    fn test_h2a2() -> Result<(), Error> {
        let tx_ac = "NM_033445.2";
        let alt_ac = "NC_000001.10";
        let test_cases = vec![
            // 3'
            (
                GenomeInterval::from_str("228645560")?,
                TxInterval::from_str("1")?,
                CdsInterval::from_str("-42")?,
            ),
            (
                GenomeInterval::from_str("228645558")?,
                TxInterval::from_str("3")?,
                CdsInterval::from_str("-40")?,
            ),
            // CDS
            (
                GenomeInterval::from_str("228645519")?,
                TxInterval::from_str("42")?,
                CdsInterval::from_str("-1")?,
            ),
            (
                GenomeInterval::from_str("228645518")?,
                TxInterval::from_str("43")?,
                CdsInterval::from_str("1")?,
            ),
            // 5'
            (
                GenomeInterval::from_str("228645126")?,
                TxInterval::from_str("435")?,
                CdsInterval::from_str("393")?,
            ),
            (
                GenomeInterval::from_str("228645125")?,
                TxInterval::from_str("436")?,
                CdsInterval::from_str("*1")?,
            ),
            (
                GenomeInterval::from_str("228645124")?,
                TxInterval::from_str("437")?,
                CdsInterval::from_str("*2")?,
            ),
            (
                GenomeInterval::from_str("228645065")?,
                TxInterval::from_str("496")?,
                CdsInterval::from_str("*61")?,
            ),
        ];

        run_test_cases(tx_ac, alt_ac, &test_cases)?;

        Ok(())
    }

    /// NM_014357.4: LCE2B, two exons, strand = +1, all coordinate input/output are in HGVS
    #[test]
    fn test_lce2b() -> Result<(), Error> {
        let tx_ac = "NM_014357.4";
        let alt_ac = "NC_000001.10";
        let test_cases = vec![
            // 5'
            (
                GenomeInterval::from_str("152658599")?,
                TxInterval::from_str("1")?,
                CdsInterval::from_str("-54")?,
            ),
            (
                GenomeInterval::from_str("152658601")?,
                TxInterval::from_str("3")?,
                CdsInterval::from_str("-52")?,
            ),
            // CDS
            (
                GenomeInterval::from_str("152659319")?,
                TxInterval::from_str("54")?,
                CdsInterval::from_str("-1")?,
            ),
            (
                GenomeInterval::from_str("152659320")?,
                TxInterval::from_str("55")?,
                CdsInterval::from_str("1")?,
            ),
            // around end of exon 1
            (
                GenomeInterval::from_str("152658632")?,
                TxInterval::from_str("34")?,
                CdsInterval::from_str("-21")?,
            ),
            (
                GenomeInterval::from_str("152658633")?,
                TxInterval::from_str("34+1")?,
                CdsInterval::from_str("-21+1")?,
            ),
            // span
            (
                GenomeInterval::from_str("152658633_152659299")?,
                TxInterval::from_str("34+1_35-1")?,
                CdsInterval::from_str("-21+1_-20-1")?,
            ),
            // around start of exon 2
            (
                GenomeInterval::from_str("152659300")?,
                TxInterval::from_str("35")?,
                CdsInterval::from_str("-20")?,
            ),
            (
                GenomeInterval::from_str("152659299")?,
                TxInterval::from_str("35-1")?,
                CdsInterval::from_str("-20-1")?,
            ),
            // around end of exon 2
            (
                GenomeInterval::from_str("152659652")?,
                TxInterval::from_str("387")?,
                CdsInterval::from_str("333")?,
            ),
            (
                GenomeInterval::from_str("152659653")?,
                TxInterval::from_str("388")?,
                CdsInterval::from_str("*1")?,
            ),
            // span
            (
                GenomeInterval::from_str("152659651_152659654")?,
                TxInterval::from_str("386_389")?,
                CdsInterval::from_str("332_*2")?,
            ),
            // 3'
            (
                GenomeInterval::from_str("152659877")?,
                TxInterval::from_str("612")?,
                CdsInterval::from_str("*225")?,
            ),
        ];

        run_test_cases(tx_ac, alt_ac, &test_cases)?;

        Ok(())
    }

    /// NM_178449.3: PTH2, two exons, strand = -1, all coordinate input/output are in HGVS
    #[test]
    fn test_pth2() -> Result<(), Error> {
        let tx_ac = "NM_178449.3";
        let alt_ac = "NC_000019.9";
        let test_cases = vec![
            // 3'
            (
                GenomeInterval::from_str("49926698")?,
                TxInterval::from_str("1")?,
                CdsInterval::from_str("-102")?,
            ),
            // CDS
            (
                GenomeInterval::from_str("49926597")?,
                TxInterval::from_str("102")?,
                CdsInterval::from_str("-1")?,
            ),
            (
                GenomeInterval::from_str("49926596")?,
                TxInterval::from_str("103")?,
                CdsInterval::from_str("1")?,
            ),
            // around end of exon 1
            (
                GenomeInterval::from_str("49926469")?,
                TxInterval::from_str("230")?,
                CdsInterval::from_str("128")?,
            ),
            (
                GenomeInterval::from_str("49926468")?,
                TxInterval::from_str("230+1")?,
                CdsInterval::from_str("128+1")?,
            ),
            // span
            (
                GenomeInterval::from_str("49925901_49926467")?,
                TxInterval::from_str("230+2_231-2")?,
                CdsInterval::from_str("128+2_129-2")?,
            ),
            // around start of exons 2
            (
                GenomeInterval::from_str("49925900")?,
                TxInterval::from_str("231-1")?,
                CdsInterval::from_str("129-1")?,
            ),
            (
                GenomeInterval::from_str("49925899")?,
                TxInterval::from_str("231")?,
                CdsInterval::from_str("129")?,
            ),
            // around end of exon2
            (
                GenomeInterval::from_str("49925725")?,
                TxInterval::from_str("405")?,
                CdsInterval::from_str("303")?,
            ),
            (
                GenomeInterval::from_str("49925724")?,
                TxInterval::from_str("406")?,
                CdsInterval::from_str("*1")?,
            ),
            (
                GenomeInterval::from_str("49925671")?,
                TxInterval::from_str("459")?,
                CdsInterval::from_str("*54")?,
            ),
        ];

        run_test_cases(tx_ac, alt_ac, &test_cases)?;

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
