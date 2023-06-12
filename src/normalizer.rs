//! Variant normalization.

use std::{cmp::Ordering, ops::Range, rc::Rc, sync::Arc};

pub use crate::normalizer::error::Error;
use crate::{
    data::interface::Provider,
    mapper::variant,
    parser::{
        GenomeInterval, GenomeLocEdit, HgvsVariant, MtInterval, MtLocEdit, Mu, NaEdit, RnaInterval,
        RnaLocEdit, RnaPos, TxInterval, TxLocEdit, TxPos,
    },
    sequences::{revcomp, trim_common_prefixes, trim_common_suffixes},
    validator::Validator,
};

mod error {
    /// Error type for normalization of HGVS expressins.
    #[derive(thiserror::Error, Debug)]
    pub enum Error {
        #[error("integer conversion failed")]
        IntegerConversion(#[from] std::num::TryFromIntError),
        #[error("validation error")]
        ValidationFailed(#[from] crate::validator::Error),
        #[error("problem accessing data")]
        DataError(#[from] crate::data::error::Error),
        #[error("replacing reference failed: {0}")]
        ReplaceReferenceFailed(String),
        #[error("c_to_n mapping failed for {0}")]
        CToNMappingFailed(String),
        #[error("n_to_c mapping failed for {0}")]
        NToCMappingFailed(String),
        #[error("validation failed in normalization: {0}")]
        Validation(String),
        #[error("cannot normalize protein-level variant: {0}")]
        ProteinVariant(String),
        #[error("cannot normalize intronic variant: {0}")]
        IntronicVariant(String),
        #[error("coordinates are out of bound in: {0}")]
        CoordinatesOutOfBounds(String),
        #[error("cannot find exon for normalization of start of: {0}")]
        ExonNotFoundForStart(String),
        #[error("cannot find exon for normalization of end of: {0}")]
        ExonNotFoundForEnd(String),
        #[error("normalization unsupported when spanning exon-intron boundary: {0}")]
        ExonIntronBoundary(String),
        #[error("normalization unsupported when spanning UTR-exon boundary: {0}")]
        UtrExonBoundary(String),
        #[error("variant span is outside of sequence bounds: {0}")]
        VariantSpanOutsideSequenceBounds(String),
    }
}

/// A direction with respect to a sequence.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Direction {
    ThreeToFive,
    FiveToThree,
}

/// Configuration for the normalizer.
#[derive(Debug, Clone)]
pub struct Config {
    pub alt_aln_method: String,
    pub cross_boundaries: bool,
    pub shuffle_direction: Direction,
    pub replace_reference: bool,
    /// Whether to validate sequence length.  Can be disabled for mapping from genome
    // TODO: inconsistent with passing in the validator...
    #[allow(dead_code)]
    pub validate: bool,
    pub window_size: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            alt_aln_method: "splign".to_string(),
            cross_boundaries: false,
            shuffle_direction: Direction::FiveToThree,
            replace_reference: true,
            validate: true,
            window_size: 20,
        }
    }
}

/// Normalizes variants (5' and 3' shifting).
pub struct Normalizer<'a> {
    pub provider: Arc<dyn Provider + Send + Sync>,
    pub validator: Arc<dyn Validator + Send + Sync>,
    pub config: Config,
    pub mapper: &'a variant::Mapper,
}

/// Helper type used in `Normalizer::check_and_guard()`.
struct CheckAndGuardResult {
    var: HgvsVariant,
    as_is: bool,
    cds_to_tx: bool,
}

impl<'a> Normalizer<'a> {
    pub fn new(
        mapper: &'a variant::Mapper,
        provider: Arc<dyn Provider + Send + Sync>,
        validator: Arc<dyn Validator + Send + Sync>,
        config: Config,
    ) -> Self {
        Self {
            mapper,
            provider,
            validator,
            config,
        }
    }

    pub fn normalize(&self, var: &HgvsVariant) -> Result<HgvsVariant, Error> {
        let is_genome = matches!(&var, HgvsVariant::GenomeVariant { .. });

        // Run the pre-normalization checks (a) whether trying to normalize the variant is an
        // error, (b) whether the variant should be returned as is, and (c) whether the variant
        // was translated from CDS to transcript and needs translation into the other direction.
        let CheckAndGuardResult {
            var,
            as_is,
            cds_to_tx,
        } = self.check_and_guard(var, is_genome)?;
        if as_is {
            return Ok(var);
        }

        // Compute boundary for the shuffling depending on the configuration and normalize
        // the alleles.
        let boundary = self.get_boundary(&var)?;
        let (start, end, reference, alternative) =
            self.normalize_alleles(&var, boundary.clone())?;

        // Build a new variant from the given alleles and possibly project back to the CDS.
        self.build_result(var, start, end, reference, alternative, boundary, cds_to_tx)
    }

    // # Args
    //
    // * `is_genome` -- allows for disabling length validation for genome (where contigs are likely
    //   missing from the provider to save memory)
    fn check_and_guard(
        &self,
        orig_var: &HgvsVariant,
        is_genome: bool,
    ) -> Result<CheckAndGuardResult, Error> {
        let var = orig_var.clone();
        self.validator
            .validate(&var)
            .map_err(|e| Error::Validation(e.to_string()))?;

        // Bail out if we cannot reliably normalize.
        if var.mu_na_edit().map(|e| !e.is_certain()).unwrap_or(true)
            || var.mu_loc_range().map(|e| !e.is_certain()).unwrap_or(false)
        {
            return Ok(CheckAndGuardResult {
                var,
                as_is: true,
                cds_to_tx: false,
            });
        }

        // Guard against unsupported normalizations.
        if let HgvsVariant::ProtVariant { .. } = var {
            return Err(Error::ProteinVariant(format!("{}", var)));
        }

        // NB: once we support gene conversions, guard against this here as well.

        let var = if self.config.replace_reference {
            self.mapper
                .replace_reference(var.clone())
                .map_err(|_e| Error::ReplaceReferenceFailed(format!("{}", var)))?
        } else {
            var
        };

        // Guard against identity variants.
        if let Some(na_edit) = var.na_edit() {
            if format!("{na_edit}") == "=" {
                return Ok(CheckAndGuardResult {
                    var,
                    as_is: true,
                    cds_to_tx: false,
                });
            }
        }

        // For CDS variants, first convert to transcript variatn and perform normalization on this
        // level.  Then, convert normalized variant back.
        let (var, cds_to_tx) = if let HgvsVariant::CdsVariant { .. } = var {
            (
                self.mapper
                    .c_to_n(&var)
                    .map_err(|_e| Error::CToNMappingFailed(format!("{}", var)))?,
                true,
            )
        } else {
            (var, false)
        };

        // Guard against intronic variants.
        match &var {
            HgvsVariant::TxVariant { .. } | HgvsVariant::RnaVariant { .. } => {
                if var.spans_intron() {
                    Err(Error::IntronicVariant(format!("{}", &var)))
                } else {
                    Ok(())
                }
            }
            _ => Ok(()),
        }?;

        fn valid_seq_len(provider: &dyn Provider, ac: &str, len: usize) -> Result<bool, Error> {
            let res = provider.get_seq_part(ac, Some(len - 1), Some(len))?;
            Ok(!res.is_empty())
        }

        // Bail out when coordinates are out of bounds.
        if let Some(var_loc_range) = var.loc_range() {
            if var_loc_range.start < 0
                || !is_genome
                    && !valid_seq_len(
                        self.provider.as_ref(),
                        var.accession(),
                        var_loc_range.end as usize,
                    )?
            {
                return Err(Error::CoordinatesOutOfBounds(format!("{}", &var)));
            }
        } else {
            panic!("Cannot happen; guarded against above")
        }

        // Assertion against any unsupported variant type.
        assert!(
            matches!(
                &var,
                HgvsVariant::GenomeVariant { .. }
                    | HgvsVariant::MtVariant { .. }
                    | HgvsVariant::TxVariant { .. }
                    | HgvsVariant::RnaVariant { .. }
            ),
            "Variant must have one of the types gmnr"
        );

        Ok(CheckAndGuardResult {
            var,
            as_is: false,
            cds_to_tx,
        })
    }

    /// Obtain position of exon-intron boundary for the current variant.
    fn get_boundary(&self, var: &HgvsVariant) -> Result<Range<i32>, Error> {
        if !self.config.cross_boundaries
            && matches!(
                &var,
                HgvsVariant::RnaVariant { .. } | HgvsVariant::TxVariant { .. }
            )
        {
            // Obtain genomic accession.
            let map_info = self
                .provider
                .as_ref()
                .get_tx_mapping_options(var.accession())?;
            let map_info = map_info
                .into_iter()
                .filter(|r| r.alt_aln_method == self.config.alt_aln_method)
                .collect::<Vec<_>>();
            let alt_ac = &map_info[0].alt_ac;

            // Obtain tx info.
            let tx_info = self.provider.as_ref().get_tx_info(
                var.accession(),
                alt_ac,
                &self.config.alt_aln_method,
            )?;
            let cds_start = tx_info.cds_start_i;
            let cds_end = tx_info.cds_end_i;

            // Obtain exon info.
            let exon_info = self.provider.as_ref().get_tx_exons(
                var.accession(),
                alt_ac,
                &self.config.alt_aln_method,
            )?;
            let mut exon_starts = exon_info.iter().map(|r| r.tx_start_i).collect::<Vec<_>>();
            exon_starts.sort();
            let mut exon_ends = exon_info.iter().map(|r| r.tx_end_i).collect::<Vec<_>>();
            exon_ends.sort();
            exon_starts.push(
                *exon_ends
                    .last()
                    .expect("should not happen; must have at least one exon"),
            );
            exon_ends.push(std::i32::MAX);

            // Find the end pos of the exon where the var locates.
            let _left = 0;
            let _right = std::i32::MAX;

            // NB: the following content is from the original Python code.
            // TODO: #242: implement methods to find tx regions
            let loc_range = var
                .loc_range()
                .expect("location must have a concrete position");
            let i = exon_starts
                .iter()
                .enumerate()
                .find(|(idx, _x)| {
                    loc_range.start >= exon_starts[*idx] && loc_range.start < exon_ends[*idx]
                })
                .ok_or(Error::ExonNotFoundForStart(format!("{}", &var)))?
                .0;
            let j = exon_starts
                .iter()
                .enumerate()
                .find(|(idx, _x)| {
                    loc_range.end > exon_starts[*idx] && loc_range.end - 1 < exon_ends[*idx]
                })
                .ok_or(Error::ExonNotFoundForEnd(format!("{}", &var)))?
                .0;

            if i != j {
                return Err(Error::ExonIntronBoundary(format!("{}", &var)));
            }

            let mut left = exon_starts[i];
            let mut right = exon_ends[i];

            if let Some(cds_start) = cds_start {
                if loc_range.end - 1 < cds_start {
                    right = std::cmp::min(right, cds_start);
                } else if loc_range.start >= cds_start {
                    left = std::cmp::max(left, cds_start);
                } else {
                    return Err(Error::UtrExonBoundary(format!("{}", &var)));
                }
            }

            if let Some(cds_end) = cds_end {
                if loc_range.start >= cds_end {
                    left = std::cmp::max(left, cds_end);
                } else if loc_range.end - 1 < cds_end {
                    right = std::cmp::min(right, cds_end);
                } else {
                    return Err(Error::UtrExonBoundary(format!("{}", &var)));
                }
            }

            Ok(left..right)
        } else {
            // For variant types g., m., etc.
            Ok(0..std::i32::MAX)
        }
    }

    /// NB: The returned start/end are 1-based!
    fn normalize_alleles(
        &self,
        var: &HgvsVariant,
        boundary: Range<i32>,
    ) -> Result<(i32, i32, String, String), Error> {
        let (reference, alternative) = self.get_ref_alt(var, &boundary)?;
        let win_size = self.config.window_size;

        if self.config.shuffle_direction == Direction::FiveToThree {
            self.normalize_alleles_5_to_3(
                reference,
                alternative,
                win_size.try_into()?,
                var,
                boundary,
            )
        } else {
            self.normalize_alleles_3_to_5(
                reference,
                alternative,
                win_size.try_into()?,
                var,
                boundary,
            )
        }
    }

    fn normalize_alleles_5_to_3(
        &self,
        reference: String,
        alternative: String,
        win_size: i32,
        var: &HgvsVariant,
        boundary: Range<i32>,
    ) -> Result<(i32, i32, String, String), Error> {
        let mut reference = reference;
        let mut alternative = alternative;
        let loc_range = var
            .loc_range()
            .expect("must have a concrete base pair location");
        let (mut base, mut start, mut stop) = match var
            .na_edit()
            .expect("should not happen; must have a NAEdit")
        {
            NaEdit::Ins { .. } => (loc_range.start + 1, 1, 1),
            NaEdit::Dup { .. } => (loc_range.end, 1, 1),
            _ => (loc_range.start + 1, 0, loc_range.end - loc_range.start),
        };

        loop {
            let ref_seq = self.fetch_bounded_seq(
                var,
                base - 1,
                base + stop - 1 + win_size,
                win_size,
                &boundary,
            )?;
            if ref_seq.is_empty() {
                break;
            }
            let (orig_start, _orig_stop) = (start, stop);
            (start, stop, reference, alternative) = normalize_alleles(
                &ref_seq,
                start,
                stop,
                reference,
                alternative,
                ref_seq.len(),
                win_size,
                false,
            )?;
            if stop < ref_seq.len().try_into()? || start == orig_start {
                break;
            }
            // If stop at the end of the window, try to extend the shuffling to the right.
            base += start - orig_start;
            stop -= start - orig_start;
            start = orig_start;
        }

        Ok((base + start, base + stop, reference, alternative))
    }

    fn normalize_alleles_3_to_5(
        &self,
        reference: String,
        alternative: String,
        win_size: i32,
        var: &HgvsVariant,
        boundary: Range<i32>,
    ) -> Result<(i32, i32, String, String), Error> {
        let mut reference = reference;
        let mut alternative = alternative;
        let loc_range = var
            .loc_range()
            .expect("must have a concrete base pair location");
        let mut base = std::cmp::max(loc_range.start + 1 - win_size, 1);
        let (mut start, mut stop) = match var
            .na_edit()
            .expect("should not happen; must have a NAEdit")
        {
            NaEdit::Ins { .. } => (loc_range.end - base, loc_range.end - base),
            NaEdit::Dup { .. } => (loc_range.end - base + 1, loc_range.end - base + 1),
            _ => (loc_range.start + 1 - base, loc_range.end - base + 1),
        };

        loop {
            if base < boundary.start + 1 {
                start -= boundary.start + 1 - base;
                stop -= boundary.start + 1 - base;
                base = boundary.start + 1;
            }
            let ref_seq =
                self.fetch_bounded_seq(var, base - 1, base + stop - 1, start, &boundary)?;
            if ref_seq.is_empty() {
                break;
            }
            let (_orig_start, orig_stop) = (start, stop);
            (start, stop, reference, alternative) = normalize_alleles(
                &ref_seq,
                start,
                stop,
                reference,
                alternative,
                0,
                win_size,
                true,
            )?;
            if start > 0 || stop == orig_stop {
                break;
            }
            // If stop at the end of the window, try to extend the shuffling to the left.
            base -= orig_stop - stop;
            start += orig_stop - stop;
            stop = orig_stop;
        }
        Ok((base + start, base + stop, reference, alternative))
    }

    /// NB: The parameter start/end are 1-based!
    #[allow(clippy::too_many_arguments)]
    fn build_result(
        &self,
        var: HgvsVariant,
        start: i32,
        end: i32,
        reference: String,
        alternative: String,
        boundary: Range<i32>,
        cds_to_tx: bool,
    ) -> Result<HgvsVariant, Error> {
        let ref_len = reference.len() as i32;
        let alt_len = alternative.len() as i32;

        let (ref_start, ref_end, edit) = match alt_len.cmp(&ref_len) {
            Ordering::Equal => {
                self.build_result_len_eq(start, end, ref_len, &reference, &alternative, alt_len)
            }
            Ordering::Less => {
                self.build_result_len_less(start, end, alt_len, &reference, &alternative)
            }
            Ordering::Greater => self.build_result_len_gt(
                ref_len,
                &var,
                start,
                alt_len,
                end,
                &boundary,
                &alternative,
                &reference,
            )?,
        };

        // Ensure the start is not 0.
        let (ref_start, ref_end, edit, _reference, alternative) = if ref_start == 0 {
            let reference = self.fetch_bounded_seq(&var, 0, 1, 0, &boundary)?;
            let alternative = format!("{}{}", alternative, &reference);

            (
                1,
                1,
                NaEdit::RefAlt {
                    reference: reference.clone(),
                    alternative: alternative.clone(),
                },
                reference,
                alternative,
            )
        } else {
            (ref_start, ref_end, edit, reference, alternative)
        };

        // Ensure the end is not outside of reference sequence.
        let tgt_len = self.get_tgt_len(&var)?;
        let (ref_start, ref_end, edit) = if ref_end == tgt_len.saturating_add(1) {
            let reference = self.fetch_bounded_seq(&var, tgt_len - 1, tgt_len, 0, &boundary)?;
            let alternative = format!("{}{}", &reference, alternative);
            (
                tgt_len,
                tgt_len,
                NaEdit::RefAlt {
                    reference,
                    alternative,
                },
            )
        } else {
            (ref_start, ref_end, edit)
        };

        let res = self.build_result_construct(var, ref_start, ref_end, edit, cds_to_tx)?;
        Ok(res)
    }

    fn build_result_construct(
        &self,
        var: HgvsVariant,
        ref_start: i32,
        ref_end: i32,
        edit: NaEdit,
        cds_to_tx: bool,
    ) -> Result<HgvsVariant, Error> {
        Ok(match var {
            HgvsVariant::GenomeVariant {
                accession,
                gene_symbol,
                ..
            } => HgvsVariant::GenomeVariant {
                accession,
                gene_symbol,
                loc_edit: GenomeLocEdit {
                    loc: Mu::Certain(GenomeInterval {
                        start: Some(ref_start),
                        end: Some(ref_end),
                    }),
                    edit: Mu::Certain(edit),
                },
            },
            HgvsVariant::MtVariant {
                accession,
                gene_symbol,
                ..
            } => HgvsVariant::MtVariant {
                accession,
                gene_symbol,
                loc_edit: MtLocEdit {
                    loc: Mu::Certain(MtInterval {
                        start: Some(ref_start),
                        end: Some(ref_end),
                    }),
                    edit: Mu::Certain(edit),
                },
            },
            HgvsVariant::TxVariant {
                accession,
                gene_symbol,
                ..
            } => {
                let var_t = HgvsVariant::TxVariant {
                    accession,
                    gene_symbol,
                    loc_edit: TxLocEdit {
                        loc: Mu::Certain(TxInterval {
                            start: TxPos {
                                base: ref_start,
                                offset: None,
                            },
                            end: TxPos {
                                base: ref_end,
                                offset: None,
                            },
                        }),
                        edit: Mu::Certain(edit),
                    },
                };

                if cds_to_tx {
                    self.mapper
                        .n_to_c(&var_t)
                        .map_err(|_e| Error::NToCMappingFailed(format!("{}", var_t)))?
                } else {
                    var_t
                }
            }
            HgvsVariant::RnaVariant {
                accession,
                gene_symbol,
                ..
            } => HgvsVariant::RnaVariant {
                accession,
                gene_symbol,
                loc_edit: RnaLocEdit {
                    loc: Mu::Certain(RnaInterval {
                        start: RnaPos {
                            base: ref_start,
                            offset: None,
                        },
                        end: RnaPos {
                            base: ref_end,
                            offset: None,
                        },
                    }),
                    edit: Mu::Certain(edit),
                },
            },
            _ => panic!("Cannot happen; variant types guarded above"),
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn build_result_len_gt(
        &self,
        ref_len: i32,
        var: &HgvsVariant,
        start: i32,
        alt_len: i32,
        end: i32,
        boundary: &Range<i32>,
        alternative: &String,
        reference: &str,
    ) -> Result<(i32, i32, NaEdit), Error> {
        Ok(if ref_len == 0 {
            let adj_seq = if self.config.shuffle_direction == Direction::FiveToThree {
                self.fetch_bounded_seq(var, start - alt_len - 1, end - 1, 0, boundary)?
            } else {
                self.fetch_bounded_seq(var, start - 1, start + alt_len - 1, 0, boundary)?
            };

            if *alternative != adj_seq {
                // ins
                (
                    start - 1,
                    end,
                    NaEdit::Ins {
                        alternative: alternative.clone(),
                    },
                )
            } else {
                // dup
                if self.config.shuffle_direction == Direction::FiveToThree {
                    (
                        start - alt_len,
                        end - 1,
                        NaEdit::Dup {
                            reference: alternative.clone(),
                        },
                    )
                } else {
                    (
                        start,
                        start + alt_len - 1,
                        NaEdit::Dup {
                            reference: alternative.clone(),
                        },
                    )
                }
            }
        } else {
            // delins
            (
                start,
                end - 1,
                NaEdit::RefAlt {
                    reference: reference.to_owned(),
                    alternative: alternative.clone(),
                },
            )
        })
    }

    /// Build the result for the del/delins case.
    fn build_result_len_less(
        &self,
        start: i32,
        end: i32,
        alt_len: i32,
        reference: &str,
        alternative: &str,
    ) -> (i32, i32, NaEdit) {
        (
            start,
            end - 1,
            if alt_len == 0 {
                NaEdit::DelRef {
                    reference: reference.to_owned(),
                }
            } else {
                NaEdit::RefAlt {
                    reference: reference.to_owned(),
                    alternative: alternative.to_owned(),
                }
            },
        )
    }

    /// Build the result for the case of ref.len() == alt.len()
    fn build_result_len_eq(
        &self,
        start: i32,
        end: i32,
        ref_len: i32,
        reference: &str,
        alternative: &str,
        alt_len: i32,
    ) -> (i32, i32, NaEdit) {
        let ref_start = start;
        let ref_end = end - 1;

        if ref_len > 1 && *reference == revcomp(alternative) {
            // inversion
            (
                ref_start,
                ref_end,
                NaEdit::InvRef {
                    reference: reference.to_owned(),
                },
            )
        } else if ref_len == 0 && alt_len == 0 {
            // identity
            (
                ref_end,
                ref_end,
                NaEdit::RefAlt {
                    reference: reference.to_owned(),
                    alternative: alternative.to_owned(),
                },
            )
        } else {
            // substitution or delins
            (
                ref_start,
                ref_end,
                if alt_len == 0 {
                    NaEdit::DelRef {
                        reference: reference.to_owned(),
                    }
                } else {
                    NaEdit::RefAlt {
                        reference: reference.to_owned(),
                        alternative: alternative.to_owned(),
                    }
                },
            )
        }
    }

    /// Fetch reference sequence from HGVS data provider.
    ///
    /// The start position is 0 and the interval is half-open.
    fn fetch_bounded_seq(
        &self,
        var: &HgvsVariant,
        start: i32,
        end: i32,
        window_size: i32,
        boundary: &Range<i32>,
    ) -> Result<String, Error> {
        let var_len = end - start - window_size;

        let start = std::cmp::max(start, boundary.start);
        let end = std::cmp::min(end, boundary.end);
        if start >= end {
            return Ok("".to_string());
        }

        let seq = self.provider.get_seq_part(
            var.accession(),
            Some(start.try_into()?),
            Some(end.try_into()?),
        )?;
        let seq_len: i32 = seq.len().try_into()?;

        if seq_len < end - start && seq_len < var_len {
            Err(Error::VariantSpanOutsideSequenceBounds(format!("{}", &var)))
        } else {
            Ok(seq)
        }
    }

    fn get_tgt_len(&self, var: &HgvsVariant) -> Result<i32, Error> {
        Ok(
            if matches!(
                var,
                HgvsVariant::GenomeVariant { .. } | HgvsVariant::MtVariant { .. }
            ) {
                std::i32::MAX
            } else {
                let id_info = self.provider.get_tx_identity_info(var.accession())?;
                id_info.lengths.into_iter().sum()
            },
        )
    }

    fn get_ref_alt(
        &self,
        var: &HgvsVariant,
        boundary: &Range<i32>,
    ) -> Result<(String, String), Error> {
        // Get reference allele.
        let reference = match var
            .na_edit()
            .expect("should not happen; must have a NAEdit")
        {
            NaEdit::Ins { .. } | NaEdit::Dup { .. } => "".to_string(),
            NaEdit::RefAlt { .. } | NaEdit::InvRef { .. } | NaEdit::DelRef { .. } => {
                let loc_range = var
                    .loc_range()
                    .expect("must have a concrete base pair location");
                self.fetch_bounded_seq(var, loc_range.start, loc_range.end, 0, boundary)?
            }
            _ => panic!("Cannot work with NumAlt,DelNum/InvNum"),
        };

        // Get alternative allele.
        let alternative = match var
            .na_edit()
            .expect("should not happen; must have a NAEdit")
        {
            NaEdit::DelRef { .. } => "".to_string(),
            NaEdit::RefAlt { alternative, .. } | NaEdit::Ins { alternative } => {
                alternative.to_string()
            }
            NaEdit::Dup { .. } => {
                let loc_range = var
                    .loc_range()
                    .expect("must have a concrete base pair location");
                self.fetch_bounded_seq(var, loc_range.start, loc_range.end, 0, boundary)?
            }
            NaEdit::InvRef { .. } => revcomp(&reference),
            _ => panic!("Cannot work with NumAlt,DelNum/InvNum"),
        };

        Ok((reference, alternative))
    }
}

#[allow(clippy::too_many_arguments)]
fn normalize_alleles(
    ref_seq: &str,
    start: i32,
    stop: i32,
    reference: String,
    alternative: String,
    bound: usize,
    ref_step: i32,
    left: bool,
) -> Result<(i32, i32, String, String), Error> {
    if left {
        normalize_alleles_left(
            ref_seq,
            start.try_into()?,
            stop.try_into()?,
            reference,
            alternative,
            bound,
            ref_step.try_into()?,
        )
    } else {
        normalize_alleles_right(
            ref_seq,
            start.try_into()?,
            stop.try_into()?,
            reference,
            alternative,
            bound,
            ref_step.try_into()?,
        )
    }
}

fn normalize_alleles_left(
    ref_seq: &str,
    start: usize,
    stop: usize,
    reference: String,
    alternative: String,
    bound: usize,
    ref_step: usize,
) -> Result<(i32, i32, String, String), Error> {
    // Step 1: Trim common suffix./
    let (trimmed, reference, alternative) = trim_common_suffixes(&reference, &alternative);
    let mut stop = stop - trimmed;

    // Step 2: Trim common prefix.
    let (trimmed, mut reference, mut alternative) = trim_common_prefixes(&reference, &alternative);
    let mut start = start + trimmed;

    // Step 3: While a null allele exists, right shuffle by appending alleles with reference
    //         and trimming common prefixes.

    let shuffle = true;
    while shuffle && (reference.is_empty() || alternative.is_empty()) && start > bound {
        let step = std::cmp::min(ref_step, start - bound);

        let r = ref_seq[(start - step)..(start - bound)].to_uppercase();
        let new_reference = format!("{r}{reference}");
        let new_alternative = format!("{r}{alternative}");

        let (trimmed, new_reference, new_alternative) =
            trim_common_suffixes(&new_reference, &new_alternative);

        if trimmed == 0 {
            break;
        }

        start -= trimmed;
        stop -= trimmed;

        if trimmed == step {
            reference = new_reference;
            alternative = new_alternative;
        } else {
            let left = step - trimmed;
            let r = left..;
            reference = new_reference[r.clone()].to_string();
            alternative = new_alternative[r].to_string();
            break;
        }
    }

    Ok((start.try_into()?, stop.try_into()?, reference, alternative))
}

fn normalize_alleles_right(
    ref_seq: &str,
    start: usize,
    stop: usize,
    reference: String,
    alternative: String,
    bound: usize,
    ref_step: usize,
) -> Result<(i32, i32, String, String), Error> {
    // Step 1: Trim common prefix.
    let (trimmed, reference, alternative) = trim_common_prefixes(&reference, &alternative);
    let mut start = start + trimmed;

    // Step 2: Trim common suffix.
    let (trimmed, mut reference, mut alternative) = trim_common_suffixes(&reference, &alternative);
    let mut stop = stop - trimmed;

    // Step 3: While a null allele exists, right shuffle by appending alleles with reference
    //         and trimming common prefixes.

    let shuffle = true;
    while shuffle && (reference.is_empty() || alternative.is_empty()) && stop < bound {
        let step = std::cmp::min(ref_step, bound - stop);

        let r = ref_seq[stop..(stop + step)].to_uppercase();
        let new_reference = format!("{reference}{r}");
        let new_alternative = format!("{alternative}{r}");

        let (trimmed, new_reference, new_alternative) =
            trim_common_prefixes(&new_reference, &new_alternative);

        if trimmed == 0 {
            break;
        }

        start += trimmed;
        stop += trimmed;

        if trimmed == step {
            reference = new_reference;
            alternative = new_alternative;
        } else {
            let left = step - trimmed;
            let r = ..(new_reference.len() - left);
            reference = new_reference[r].to_string();
            let r = ..(new_alternative.len() - left);
            alternative = new_alternative[r].to_string();
            break;
        }
    }

    Ok((start.try_into()?, stop.try_into()?, reference, alternative))
}

#[cfg(test)]
mod test {
    use test_log::test;

    use anyhow::Error;
    use std::{rc::Rc, str::FromStr};

    use pretty_assertions::assert_eq;

    use super::{Config, Direction, Normalizer};
    use crate::{
        data::uta_sr::test_helpers::build_provider,
        mapper::variant::Mapper,
        parser::{HgvsVariant, NoRef},
        validator::IntrinsicValidator,
    };

    fn normalizers(
        mapper: &Mapper,
    ) -> Result<(Normalizer, Normalizer, Normalizer, Normalizer), Error> {
        let provider = mapper.provider();
        let validator = Arc::new(IntrinsicValidator::new(true));

        Ok((
            Normalizer::new(
                mapper,
                provider.clone(),
                validator.clone(),
                Config {
                    shuffle_direction: Direction::FiveToThree,
                    cross_boundaries: true,
                    ..Default::default()
                },
            ),
            Normalizer::new(
                mapper,
                provider.clone(),
                validator.clone(),
                Config {
                    shuffle_direction: Direction::ThreeToFive,
                    cross_boundaries: true,
                    ..Default::default()
                },
            ),
            Normalizer::new(
                mapper,
                provider.clone(),
                validator.clone(),
                Config {
                    shuffle_direction: Direction::FiveToThree,
                    cross_boundaries: false,
                    ..Default::default()
                },
            ),
            Normalizer::new(
                mapper,
                provider.clone(),
                validator,
                Config {
                    shuffle_direction: Direction::ThreeToFive,
                    cross_boundaries: false,
                    ..Default::default()
                },
            ),
        ))
    }

    #[test]
    fn normalize_cds_3_prime_shuffling() -> Result<(), Error> {
        let mapper = Mapper::new(&Default::default(), build_provider()?);
        let (norm, _norm5, _normc, _norm5c) = normalizers(&mapper)?;

        // The following tests correspond to the Python ones "3' shuffling" and "5' shuffling".

        let cases3 = vec![
            // gene COL1A
            ("NM_000088.3:c.589_600inv", "NM_000088.3:c.590_599inv"),
            // gene DEFB133
            ("NM_001166478.1:c.31del", "NM_001166478.1:c.35del"),
            ("NM_001166478.1:c.35_36insT", "NM_001166478.1:c.35dup"),
            ("NM_001166478.1:c.36_37insTC", "NM_001166478.1:c.36_37dup"),
            ("NM_001166478.1:c.35_36dup", "NM_001166478.1:c.36_37dup"),
            (
                "NM_001166478.1:c.2_7delinsTTTAGA",
                "NM_001166478.1:c.3_4delinsTT",
            ),
            ("NM_001166478.1:c.30_31insT", "NM_001166478.1:c.35dup"),
            ("NM_001166478.1:c.59delG", "NM_001166478.1:c.61del"),
            (
                "NM_001166478.1:c.36_37insTCTCTC",
                "NM_001166478.1:c.37_38insCTCTCT",
            ),
            // gene ATM
            ("NM_000051.3:c.14_15insT", "NM_000051.3:c.15dup"),
        ];

        for (input, exp_3) in cases3 {
            let raw = HgvsVariant::from_str(input)?;

            // 5' -> 3' shuffling
            let res_3 = norm.normalize(&raw)?;
            assert_eq!(
                format!("{}", &NoRef(&res_3)),
                exp_3,
                "{:?} ~ {:?} ~ {:?}",
                &res_3,
                &raw,
                &exp_3
            );
        }

        Ok(())
    }

    #[test]
    fn normalize_cds_5_prime_shuffling() -> Result<(), Error> {
        let mapper = Mapper::new(&Default::default(), build_provider()?);
        let (_norm, norm5, _normc, _norm5c) = normalizers(&mapper)?;

        let cases5 = vec![
            // gene COL1A
            ("NM_000088.3:c.589_600inv", "NM_000088.3:c.590_599inv"),
            // gene DEFB133
            ("NM_001166478.1:c.34del", "NM_001166478.1:c.31del"),
            ("NM_001166478.1:c.35_36insT", "NM_001166478.1:c.31dup"),
            ("NM_001166478.1:c.36_37insTC", "NM_001166478.1:c.35_36dup"),
            ("NM_001166478.1:c.35_36dup", "NM_001166478.1:c.35_36dup"),
            (
                "NM_001166478.1:c.2_7delinsTTTAGA",
                "NM_001166478.1:c.3_4delinsTT",
            ),
            ("NM_001166478.1:c.30_31insT", "NM_001166478.1:c.31dup"),
            ("NM_001166478.1:c.61delG", "NM_001166478.1:c.59del"),
            (
                "NM_001166478.1:c.36_37insTCTCTC",
                "NM_001166478.1:c.34_35insTCTCTC",
            ),
            // gene ATM
            ("NM_000051.3:c.14_15insT", "NM_000051.3:c.14dup"),
        ];

        for (input, exp_5) in cases5 {
            let raw = HgvsVariant::from_str(input)?;

            // 3' -> 5' shuffling
            let res_5 = norm5.normalize(&raw)?;
            assert_eq!(
                format!("{}", &NoRef(&res_5)),
                exp_5,
                "{:?} ~ {:?} ~ {:?}",
                &res_5,
                &raw,
                &exp_5
            );
        }

        Ok(())
    }

    #[test]
    fn normalize_cds_around_exon_intron_boundary() -> Result<(), Error> {
        let mapper = Mapper::new(&Default::default(), build_provider()?);
        let (_norm, _norm5, normc, norm5c) = normalizers(&mapper)?;

        // The following tests correspond to the Python ones "Around exon-intron boundary".
        {
            let raw = HgvsVariant::from_str("NM_001166478.1:c.59delG")?;
            let exp = "NM_001166478.1:c.60del";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001166478.1:c.61delG")?;
            let exp = "NM_001166478.1:c.61del";
            let res = norm5c.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            // gene MECP2
            let raw = HgvsVariant::from_str("NM_001110792.1:c.1030_1035del")?;
            let exp = "NM_001110792.1:c.1029_1034del";
            let res = norm5c.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001166478.1:c.59_61del")?;
            // with self.assertRaises(HGVSUnsupportedOperationError):
            assert!(normc.normalize(&raw).is_err());
        }
        Ok(())
    }

    #[test]
    fn normalize_cds_utr_variant() -> Result<(), Error> {
        let mapper = Mapper::new(&Default::default(), build_provider()?);
        let (norm, norm5, normc, _norm5c) = normalizers(&mapper)?;

        // The following tests correspond to the Python ones "UTR variants".
        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.-5_-4insA")?;
            let exp = "NM_000051.3:c.-3dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.-4_-3insAC")?;
            let exp = "NM_000051.3:c.-3_-2dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.-2_-1insCA")?;
            let exp = "NM_000051.3:c.-1_1dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.-2_-1insCA")?;
            let exp = "NM_000051.3:c.-1_1insAC";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.-4_-3insA")?;
            let exp = "NM_000051.3:c.-4dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.1_2insCA")?;
            let exp = "NM_000051.3:c.-1_1dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.*2_*3insT")?;
            let exp = "NM_000051.3:c.*4dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.9170_9171insAT")?;
            let exp = "NM_000051.3:c.9171_*1dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.*4_*5insT")?;
            let exp = "NM_000051.3:c.*3dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_000051.3:c.9171_*1insA")?;
            let exp = "NM_000051.3:c.9171dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            // gene BRCA2
            //
            // NB: this used to fail in Python hgvs but works here now as we do not
            // perform comprehensive validation yet (1 bp interval/position, but 3bp
            // deleted).
            let raw = HgvsVariant::from_str("NM_000059.3:c.7790delAAG")?;
            assert_eq!(
                "NM_000059.3:c.7791delA",
                format!("{}", norm.normalize(&raw)?)
            );
        }

        Ok(())
    }

    #[test]
    fn normalize_genome_3_prime_shuffling() -> Result<(), Error> {
        let mapper = Mapper::new(&Default::default(), build_provider()?);
        let (norm, _norm5, _normc, _norm5c) = normalizers(&mapper)?;

        // 3' shuffling
        {
            // GRCh38:chr6
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917122_49917123insA")?;
            let exp = "NC_000006.11:g.49917127dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917121_49917122insGA")?;
            let exp = "NC_000006.11:g.49917122_49917123dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917122_49917123dup")?;
            let exp = "NC_000006.11:g.49917122_49917123dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917122_49917123dupGA")?;
            let exp = "NC_000006.11:g.49917122_49917123dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917098delC")?;
            let exp = "NC_000006.11:g.49917099del";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917151_49917156delinsTCTAAA")?;
            let exp = "NC_000006.11:g.49917154_49917155delinsAA";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            // GRCh37.p13:chr1
            let raw = HgvsVariant::from_str("NC_000001.10:g.1647893delinsCTTTCTT")?;
            let exp = "NC_000001.10:g.1647895_1647900dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            // GRCh38.p12:chr3
            let raw = HgvsVariant::from_str(
                "NC_000003.12:g.46709584_46709610del27insAAGAAGAAGAAGAAGAAGAAGAAGAAG",
            )?;
            let exp = "NC_000003.12:g.46709584_46709610=";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        Ok(())
    }

    #[test]
    fn normalize_genome_5_prime_shuffling() -> Result<(), Error> {
        let mapper = Mapper::new(&Default::default(), build_provider()?);
        let (norm, norm5, _normc, _norm5c) = normalizers(&mapper)?;

        // 5' shuffling
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917122_49917123insA")?;
            let exp = "NC_000006.11:g.49917123dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917121_49917122insGA")?;
            let exp = "NC_000006.11:g.49917121_49917122dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917122_49917123dup")?;
            let exp = "NC_000006.11:g.49917121_49917122dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917122_49917123dupGA")?;
            let exp = "NC_000006.11:g.49917121_49917122dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917099delC")?;
            let exp = "NC_000006.11:g.49917098del";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NC_000006.11:g.49917151_49917156delinsTCTAAA")?;
            let exp = "NC_000006.11:g.49917154_49917155delinsAA";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str(
                "NC_000003.12:g.46709584_46709610del27insAAGAAGAAGAAGAAGAAGAAGAAGAAG",
            )?;
            let exp = "NC_000003.12:g.46709584_46709610=";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NC_000009.11:g.36233991_36233992delCAinsTG")?;
            let exp = "NC_000009.11:g.36233991_36233992inv";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NG_032871.1:g.32476_53457delinsAATTAAGGTATA")?;
            // with self.assertRaises(HGVSInvalidVariantError):
            assert!(norm.normalize(&raw).is_err());
        }
        {
            let raw = HgvsVariant::from_str("NG_032871.1:g.32476_53457delinsAATTAAGGTATA")?;
            // with self.assertRaises(HGVSInvalidVariantError):
            assert!(norm.normalize(&raw).is_err());
        }

        Ok(())
    }

    #[test]
    fn normalize_near_tx_start_end() -> Result<(), Error> {
        let mapper = Mapper::new(&Default::default(), build_provider()?);
        let (norm, norm5, normc, norm5c) = normalizers(&mapper)?;

        {
            // gene OR9A4
            let raw = HgvsVariant::from_str("NM_001001656.1:c.935T>C")?;
            let exp = "NM_001001656.1:c.935T>C";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.945G>C")?;
            let exp = "NM_001001656.1:c.945G>C";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.945dup")?;
            let exp = "NM_001001656.1:c.945dup";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.935_945del")?;
            let exp = "NM_001001656.1:c.935_945del";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.946G>C")?;
            // with self.assertRaises(HGVSError):
            assert!(norm.normalize(&raw).is_err());
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.946dup")?;
            // with self.assertRaises(HGVSError):
            assert!(norm.normalize(&raw).is_err());
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.935_946del")?;
            // with self.assertRaises(HGVSError):
            assert!(norm.normalize(&raw).is_err());
        }

        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.935T>C")?;
            let exp = "NM_001001656.1:c.935T>C";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.945G>C")?;
            let exp = "NM_001001656.1:c.945G>C";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.945dup")?;
            let exp = "NM_001001656.1:c.945dup";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.935_945del")?;
            let exp = "NM_001001656.1:c.935_945del";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.946G>C")?;
            // with self.assertRaises(HGVSError):
            assert!(normc.normalize(&raw).is_err());
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.946dup")?;
            // with self.assertRaises(HGVSError):
            assert!(normc.normalize(&raw).is_err());
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.935_946del")?;
            // with self.assertRaises(HGVSError):
            assert!(normc.normalize(&raw).is_err());
        }

        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.1A>T")?;
            let exp = "NM_001001656.1:c.1A>T";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.1del")?;
            let exp = "NM_001001656.1:c.1del";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.1dup")?;
            let exp = "NM_001001656.1:c.1dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.1A>T")?;
            let exp = "NM_001001656.1:c.1A>T";
            let res = norm5c.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.1del")?;
            let exp = "NM_001001656.1:c.1del";
            let res = norm5c.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_001001656.1:c.1dup")?;
            let exp = "NM_001001656.1:c.1dup";
            let res = norm5c.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1_2insCA")?;
            let exp = "NM_212556.2:c.1_2insCA";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.2_3insCAT")?;
            let exp = "NM_212556.2:c.2_3insCAT";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1delinsCA")?;
            let exp = "NM_212556.2:c.1delinsCA";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1_2insCA")?;
            let exp = "NM_212556.2:c.1delinsACA";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.2_3insCAT")?;
            let exp = "NM_212556.2:c.1delinsATCA";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1delinsCA")?;
            let exp = "NM_212556.2:c.1delinsCA";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1delinsAA")?;
            let exp = "NM_212556.2:c.1dup";
            let res = norm5.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1400_1401insAC")?;
            let exp = "NM_212556.2:c.1401delinsACA";
            let res = norm.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1400_1401insAC")?;
            let exp = "NM_212556.2:c.1401delinsACA";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }
        {
            let raw = HgvsVariant::from_str("NM_212556.2:c.1401delinsAA")?;
            let exp = "NM_212556.2:c.1401dup";
            let res = normc.normalize(&raw)?;

            assert_eq!(
                format!("{}", &NoRef(&res)),
                exp,
                "{:?} ~ {:?} ~ {:?}",
                &res,
                &raw,
                &exp
            );
        }

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
