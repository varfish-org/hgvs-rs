//! Code for mapping variants between sequences.

use std::{ops::Range, rc::Rc};

use log::{debug, info};

use super::alignment::{Config as AlignmentConfig, Mapper as AlignmentMapper};
use crate::{
    data::interface::Provider,
    mapper::Error,
    normalizer::{self, Config as NormalizerConfig, Normalizer},
    parser::{
        Accession, CdsInterval, CdsLocEdit, CdsPos, GeneSymbol, GenomeInterval, GenomeLocEdit,
        HgvsVariant, Mu, NaEdit, TxInterval, TxLocEdit, TxPos,
    },
    sequences::revcomp,
    validator::{ValidationLevel, Validator},
};
use std::ops::Deref;

/// Configuration for Mapper.
///
/// Defaults are taken from `hgvs` Python library.
#[derive(Debug, PartialEq, Clone)]
pub struct Config {
    pub replace_reference: bool,
    pub strict_validation: bool,
    pub prevalidation_level: ValidationLevel,
    pub add_gene_symbol: bool,
    pub strict_bounds: bool,
    /// Re-normalize out of bounds genome variants on minus strand.  This can be
    /// switched off so genome sequence does not have to be available in provider.
    pub renormalize_g: bool,
    /// Use the genome sequence in case of uncertain g-to-n projections.  This
    /// can be switched off so genome sequence does not have to be available.
    pub genome_seq_available: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            replace_reference: true,
            strict_validation: false,
            prevalidation_level: ValidationLevel::Full,
            add_gene_symbol: false,
            strict_bounds: true,
            renormalize_g: true,
            genome_seq_available: true,
        }
    }
}

/// Projects variants between sequences using `alignment::Mapper`.
pub struct Mapper {
    config: Config,
    provider: Rc<dyn Provider>,
    validator: Rc<dyn Validator>,
}

/// Maps SequenceVariant objects between g., n., r., c., and p. representations.
///
/// g⟷{c,n,r} projections are similar in that c, n, and r variants
/// may use intronic coordinates. There are two essential differences
/// that distinguish the three types:
///
/// * Sequence start: In n and r variants, position 1 is the sequence
///   start; in c variants, 1 is the transcription start site.
/// * Alphabet: In n and c variants, sequences are DNA; in
///   r. variants, sequences are RNA.
///
/// This differences are summarized in this diagram::
///
/// ```text
/// g ----acgtatgcac--gtctagacgt----      ----acgtatgcac--gtctagacgt----      ----acgtatgcac--gtctagacgt----
///         \         \/         /              \         \/         /              \         \/         /
/// c      acgtATGCACGTCTAGacgt         n      acgtatgcacgtctagacgt         r      acguaugcacgucuagacgu
///            1                               1                                   1
/// p          MetHisValTer
///
/// The g excerpt and exon structures are identical. The g⟷n
/// transformation, which is the most basic, accounts for the offset
/// of the aligned sequences (shown with "1") and the exon structure.
/// The g⟷c transformation is akin to g⟷n transformation, but
/// requires an addition offset to account for the translation start
/// site (c.1).  The CDS in uppercase. The g⟷c transformation is
/// akin to g⟷n transformation with a change of alphabet.
///
/// Therefore, this this code uses g⟷n as the core transformation
/// between genomic and c, n, and r variants: All c⟷g and r⟷g
/// transformations use n⟷g after accounting for the above
/// differences. For example, c_to_g accounts for the transcription
/// start site offset, then calls n_to_g.
impl Mapper {
    pub fn new(config: &Config, provider: Rc<dyn Provider>) -> Mapper {
        let validator = config
            .prevalidation_level
            .validator(config.strict_validation, provider.clone());

        let _n_config = normalizer::Config::default();

        Mapper {
            config: config.clone(),
            provider: provider.clone(),
            validator: validator.clone(),
        }
    }

    pub fn config(&self) -> &Config {
        &self.config
    }

    /// Return a copy of the internal provider.
    pub fn provider(&self) -> Rc<dyn Provider> {
        self.provider.clone()
    }

    /// Obtain new `alignment::Mapper` for the given arguments, possibly caching results.
    fn build_alignment_mapper(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<AlignmentMapper, Error> {
        // TODO: implement caching
        AlignmentMapper::new(
            &AlignmentConfig {
                strict_bounds: self.config.strict_bounds,
            },
            self.provider.clone(),
            tx_ac,
            alt_ac,
            alt_aln_method,
        )
    }

    /// Construct a new normalizer for the variant mapper.
    pub fn normalizer(&self) -> Result<Normalizer, Error> {
        Ok(Normalizer::new(
            self,
            self.provider.clone(),
            self.validator.clone(),
            NormalizerConfig {
                replace_reference: self.config.replace_reference,
                ..Default::default()
            },
        ))
    }

    /// Convert from genome (g.) variant to transcript variant (g. or n.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    /// * `alt_al_method` -- alignment method, e.g., `splign`
    pub fn g_to_t(
        &self,
        var_g: &HgvsVariant,
        tx_ac: &str,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, Error> {
        self.validator.validate(var_g)?;
        let mapper = self.build_alignment_mapper(tx_ac, var_g.accession(), alt_aln_method)?;
        if mapper.is_coding_transcript() {
            self.g_to_c(var_g, tx_ac, alt_aln_method)
        } else {
            self.g_to_n(var_g, tx_ac, alt_aln_method)
        }
    }

    /// Convert from genome (g.) variant to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn g_to_n(
        &self,
        var_g: &HgvsVariant,
        tx_ac: &str,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, Error> {
        self.validator.validate(var_g)?;
        let var_g = if self.config.replace_reference {
            self.replace_reference(var_g.clone())?
        } else {
            var_g.clone()
        };
        if let HgvsVariant::GenomeVariant {
            accession,
            loc_edit,
            gene_symbol,
        } = &var_g
        {
            let mapper = self.build_alignment_mapper(tx_ac, &accession.value, alt_aln_method)?;

            let var_g = if mapper.strand == -1
                && !self.config.strict_bounds
                && !mapper.is_g_interval_in_bounds(loc_edit.loc.inner())
                && self.config.renormalize_g
            {
                info!("Renormalizing out-of-bounds minus strand variant on genomic sequenc");
                self.normalizer()?.normalize(&var_g)?
            } else {
                var_g.clone()
            };

            let pos_n = mapper.g_to_n(loc_edit.loc.inner())?;
            let pos_n = Mu::from(
                pos_n.inner(),
                loc_edit.loc.is_certain() && pos_n.is_certain(),
            );
            // The original Python code falls back to the genome for uncertain positions.  This
            // cannot be done if we do not have the original genome sequence.
            let pos_n_certain = pos_n.is_certain();
            let pos_n = pos_n.inner();
            let (pos_n, edit_n) = if pos_n_certain || !self.config.genome_seq_available {
                let edit_n = self.convert_edit_check_strand(mapper.strand, &loc_edit.edit)?;
                if let NaEdit::Ins { alternative } = edit_n.inner() {
                    if pos_n.start.offset.is_none()
                        && pos_n.end.offset.is_none()
                        && pos_n.end.base - pos_n.start.base > 1
                    {
                        (
                            Mu::Certain(TxInterval {
                                start: TxPos {
                                    base: pos_n.start.base + 1,
                                    ..pos_n.start
                                },
                                end: TxPos {
                                    base: pos_n.end.base - 1,
                                    ..pos_n.end
                                },
                            }),
                            Mu::Certain(NaEdit::RefAlt {
                                reference: "".to_string(),
                                alternative: alternative.clone(),
                            }),
                        )
                    } else {
                        (Mu::Certain((*pos_n).clone()), edit_n)
                    }
                } else {
                    (Mu::Certain((*pos_n).clone()), edit_n)
                }
            } else {
                // This is the how the original code handles uncertain positions.  We will reach
                // here if the position is uncertain and we have the genome sequence.
                let pos_g = mapper.n_to_g(pos_n)?;
                let edit_n = NaEdit::RefAlt {
                    reference: "".to_string(),
                    alternative: self.get_altered_sequence(
                        mapper.strand,
                        pos_g.inner().clone().try_into()?,
                        &var_g,
                    )?,
                };
                (Mu::Certain((*pos_n).clone()), Mu::Certain(edit_n))
            };

            // the following is not needed?
            // pos_n.uncertain = var_g.posedit.pos.uncertain

            let var_n = HgvsVariant::TxVariant {
                accession: Accession::new(tx_ac),
                gene_symbol: self.fetch_gene_symbol(tx_ac, gene_symbol)?,
                loc_edit: TxLocEdit {
                    loc: pos_n.clone(),
                    edit: edit_n,
                },
            };

            let var_n = if self.config.replace_reference
                && pos_n.inner().start.base >= 0
                && pos_n.inner().end.base < mapper.tgt_len
            {
                self.replace_reference(var_n)?
            } else {
                var_n
            };

            Ok(var_n)
        } else {
            Err(Error::ExpectedGenomeVariant(format!("{}", &var_g)))
        }
    }

    /// Convert from transcript variant (g.) variant to genome variant (n.).
    ///
    /// # Args
    ///
    /// * `var_n` -- `HgvsVariant::TxVariant` to project
    /// * `alt_ac` -- alternative contig accession
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn n_to_g(
        &self,
        var_n: &HgvsVariant,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, Error> {
        self.validator.validate(var_n)?;
        let var_n = self.replace_reference(var_n.clone())?;
        if let HgvsVariant::TxVariant {
            accession,
            gene_symbol: _,
            loc_edit,
        } = &var_n
        {
            let mapper = self.build_alignment_mapper(&accession.value, alt_ac, alt_aln_method)?;
            let pos_g = mapper.n_to_g(loc_edit.loc.inner())?;

            let (pos_g, edit_g) = if let Mu::Certain(pos_g) = pos_g {
                let edit_g = self.convert_edit_check_strand(mapper.strand, &loc_edit.edit)?;
                if let (NaEdit::Ins { alternative }, Some(end), Some(start)) =
                    (edit_g.inner(), pos_g.end, pos_g.start)
                {
                    if end - start > 1 {
                        (
                            Mu::Certain(GenomeInterval {
                                start: Some(start + 1),
                                end: Some(end - 1),
                            }),
                            Mu::from(
                                NaEdit::RefAlt {
                                    reference: "".to_string(),
                                    alternative: alternative.to_owned(),
                                },
                                edit_g.is_certain(),
                            ),
                        )
                    } else {
                        (Mu::Certain(pos_g), edit_g)
                    }
                } else {
                    (Mu::Certain(pos_g), edit_g)
                }
            } else {
                // variant at alignment gap
                let pos_n = mapper.g_to_n(pos_g.inner())?;
                let edit_g = NaEdit::RefAlt {
                    reference: "".to_string(),
                    alternative: self.get_altered_sequence(
                        mapper.strand,
                        pos_n.inner().clone().into(),
                        &var_n,
                    )?,
                };
                (pos_g, Mu::Certain(edit_g))
            };

            // the following is not needed?
            // pos_g.uncertain = var_n.posedit.pos.uncertain

            let var_g = HgvsVariant::GenomeVariant {
                accession: Accession::new(alt_ac),
                gene_symbol: None,
                loc_edit: GenomeLocEdit {
                    loc: pos_g,
                    edit: edit_g,
                },
            };

            let var_g = if self.config.replace_reference {
                self.replace_reference(var_g)?
            } else {
                var_g
            };

            // No gene symbol for g. variants (actually, *should* for NG, but no way to distinguish)

            Ok(var_g)
        } else {
            Err(Error::ExpectedTxVariant(format!("{}", &var_n)))
        }
    }

    /// Convert from genome (g.) variant to CDS variant (c.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn g_to_c(
        &self,
        var_g: &HgvsVariant,
        tx_ac: &str,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, Error> {
        self.validator.validate(var_g)?;
        let var_g = if self.config.replace_reference {
            self.replace_reference(var_g.clone())?
        } else {
            var_g.clone()
        };
        if let HgvsVariant::GenomeVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = &var_g
        {
            let mapper = self.build_alignment_mapper(tx_ac, &accession.value, alt_aln_method)?;
            let pos_c = mapper.g_to_c(loc_edit.loc.inner())?;

            let (pos_c, edit_c) = if let Mu::Certain(pos_c) = pos_c {
                let edit_c = self.convert_edit_check_strand(mapper.strand, &loc_edit.edit)?;
                if let NaEdit::Ins { alternative } = edit_c.inner() {
                    if pos_c.start.offset.is_none()
                        && pos_c.end.offset.is_none()
                        && pos_c.end.base - pos_c.start.base > 1
                    {
                        (
                            Mu::Certain(CdsInterval {
                                start: CdsPos {
                                    base: pos_c.start.base + 1,
                                    ..pos_c.start
                                },
                                end: CdsPos {
                                    base: pos_c.end.base - 1,
                                    ..pos_c.end
                                },
                            }),
                            Mu::Certain(NaEdit::RefAlt {
                                reference: "".to_string(),
                                alternative: alternative.clone(),
                            }),
                        )
                    } else {
                        (Mu::Certain(pos_c), edit_c)
                    }
                } else {
                    (Mu::Certain(pos_c), edit_c)
                }
            } else {
                let pos_g = mapper.c_to_g(pos_c.inner())?;
                let edit_c = NaEdit::RefAlt {
                    reference: "".to_string(),
                    alternative: self.get_altered_sequence(
                        mapper.strand,
                        pos_g.inner().clone().try_into()?,
                        &var_g,
                    )?,
                };
                (Mu::Certain((*pos_c.inner()).clone()), Mu::Certain(edit_c))
            };

            let var_c = HgvsVariant::CdsVariant {
                accession: Accession::new(tx_ac),
                gene_symbol: self.fetch_gene_symbol(tx_ac, gene_symbol)?,
                loc_edit: CdsLocEdit {
                    loc: pos_c,
                    edit: edit_c,
                },
            };

            let var_c = if self.config.replace_reference {
                self.replace_reference(var_c)?
            } else {
                var_c
            };

            Ok(var_c)
        } else {
            Err(Error::ExpectedGenomeVariant(format!("{}", &var_g)))
        }
    }

    /// Convert from CDS variant (c.) to genome variant (g.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::CdsVariant` to project
    /// * `alt_ac` -- alternative contig accession
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn c_to_g(
        &self,
        var_c: &HgvsVariant,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, Error> {
        self.validator.validate(var_c)?;
        let var_c = if self.config.replace_reference {
            self.replace_reference(var_c.clone())?
        } else {
            var_c.clone()
        };
        if let HgvsVariant::CdsVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = &var_c
        {
            let mapper = self.build_alignment_mapper(&accession.value, alt_ac, alt_aln_method)?;
            let pos_g = mapper.c_to_g(loc_edit.loc.inner())?;

            let (pos_g, edit_g) = if let Mu::Certain(pos_g) = pos_g {
                let edit_g = self.convert_edit_check_strand(mapper.strand, &loc_edit.edit)?;
                if let (NaEdit::Ins { alternative }, Some(end), Some(start)) =
                    (edit_g.inner(), pos_g.end, pos_g.start)
                {
                    if end - start > 1 {
                        (
                            Mu::Certain(GenomeInterval {
                                start: Some(start + 1),
                                end: Some(end - 1),
                            }),
                            Mu::from(
                                NaEdit::RefAlt {
                                    reference: "".to_string(),
                                    alternative: alternative.clone(),
                                },
                                edit_g.is_certain(),
                            ),
                        )
                    } else {
                        (Mu::Certain(pos_g), edit_g)
                    }
                } else {
                    (Mu::Certain(pos_g), edit_g)
                }
            } else {
                // variant at alignment gap
                let pos_n = mapper.g_to_n(pos_g.inner())?;
                let var_n = HgvsVariant::TxVariant {
                    accession: var_c.accession().clone(),
                    gene_symbol: var_c.gene_symbol().clone(),
                    loc_edit: TxLocEdit {
                        loc: pos_n.clone(),
                        edit: Mu::Certain(
                            var_c
                                .na_edit()
                                .ok_or(Error::NoNAEditInHgvsC(format!("{}", &var_c)))?
                                .clone(),
                        ),
                    },
                };
                let edit_n = NaEdit::RefAlt {
                    reference: "".to_string(),
                    alternative: self.get_altered_sequence(
                        mapper.strand,
                        pos_n.inner().clone().into(),
                        &var_n,
                    )?,
                };
                (pos_g, Mu::Certain(edit_n))
            };

            let var_g = HgvsVariant::GenomeVariant {
                accession: Accession::from(alt_ac.to_string()),
                gene_symbol: self.fetch_gene_symbol(accession.deref().as_str(), gene_symbol)?,
                loc_edit: GenomeLocEdit {
                    loc: pos_g,
                    edit: edit_g,
                },
            };

            let var_g = if self.config.replace_reference {
                self.replace_reference(var_g)?
            } else {
                var_g
            };

            Ok(var_g)
        } else {
            Err(Error::ExpectedCdsVariant(format!("{}", &var_c)))
        }
    }

    /// Convert from transcript (c. or n.) to genome (g.) variant.
    ///
    /// # Args
    ///
    /// * `var_t` -- `HgvsVariant::TxVariant` or `HgvsVariant::CdsVariant` to project
    /// * `alt_ac` -- accession of alternativ esequence
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn t_to_g(
        &self,
        var_t: &HgvsVariant,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, Error> {
        self.validator.validate(var_t)?;
        let var_t = if self.config.replace_reference {
            self.replace_reference(var_t.clone())?
        } else {
            var_t.clone()
        };
        match var_t {
            HgvsVariant::TxVariant { .. } => self.n_to_g(&var_t, alt_ac, alt_aln_method),
            HgvsVariant::CdsVariant { .. } => self.c_to_g(&var_t, alt_ac, alt_aln_method),
            _ => Err(Error::ExpectedCdsVariant(format!("{}", &var_t))),
        }
    }

    /// Convert from CDS variant (c.) to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::CdsVariant` to project
    pub fn c_to_n(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, Error> {
        log::debug!("c_to_n({})", var_c);
        self.validator.validate(var_c)?;
        let var_c = if self.config.replace_reference {
            self.replace_reference(var_c.clone())?
        } else {
            var_c.clone()
        };
        if let HgvsVariant::CdsVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = &var_c
        {
            let mapper =
                self.build_alignment_mapper(&accession.value, &accession.value, "transcript")?;

            let pos_n = mapper.c_to_n(loc_edit.loc.inner())?;
            let var_n = HgvsVariant::TxVariant {
                accession: accession.clone(),
                gene_symbol: self.fetch_gene_symbol(accession.deref().as_str(), gene_symbol)?,
                loc_edit: TxLocEdit {
                    loc: Mu::from(pos_n, loc_edit.loc.is_certain()),
                    edit: loc_edit.edit.clone(),
                },
            };

            let var_n = if self.config.replace_reference {
                self.replace_reference(var_n)?
            } else {
                var_n
            };

            log::debug!("c_to_n({}) = {}", var_c, &var_n);
            Ok(var_n)
        } else {
            Err(Error::ExpectedCdsVariant(format!("{}", &var_c)))
        }
    }

    /// Convert from transcript variant (n.) to CDS variant (c.).
    ///
    /// # Args
    ///
    /// * `var_n` -- `HgvsVariant::TxVariant` to project
    pub fn n_to_c(&self, var_n: &HgvsVariant) -> Result<HgvsVariant, Error> {
        self.validator.validate(var_n)?;
        let var_n = if self.config.replace_reference {
            self.replace_reference(var_n.clone())?
        } else {
            var_n.clone()
        };
        if let HgvsVariant::TxVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = &var_n
        {
            let mapper =
                self.build_alignment_mapper(&accession.value, &accession.value, "transcript")?;
            let pos_c = mapper.n_to_c(loc_edit.loc.inner())?;

            let var_c = HgvsVariant::CdsVariant {
                accession: accession.clone(),
                gene_symbol: self.fetch_gene_symbol(accession.deref().as_str(), gene_symbol)?,
                loc_edit: CdsLocEdit {
                    loc: Mu::from(pos_c, loc_edit.loc.is_certain()),
                    edit: loc_edit.edit.clone(),
                },
            };

            let var_c = if self.config.replace_reference {
                self.replace_reference(var_c)?
            } else {
                var_c
            };

            Ok(var_c)
        } else {
            Err(Error::ExpectedTxVariant(format!("{}", &var_n)))
        }
    }

    /// Convert from CDS variant (c.) to protein variant (p.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::TxVariant` to project
    /// * `pro_ac` -- Protein accession
    pub fn c_to_p(&self, var_c: &HgvsVariant, prot_ac: Option<&str>) -> Result<HgvsVariant, Error> {
        use super::altseq::*;

        if let HgvsVariant::CdsVariant {
            accession,
            gene_symbol: _,
            loc_edit: _,
        } = &var_c
        {
            self.validator.validate(var_c)?;

            let var_c = if self.config.replace_reference {
                self.replace_reference(var_c.clone())?
            } else {
                var_c.clone()
            };

            let reference_data = RefTranscriptData::new(
                self.provider.clone(),
                accession.deref(),
                prot_ac.map(|s| s.to_string()).as_deref(),
            )?;
            let builder = AltSeqBuilder::new(var_c, reference_data.clone());

            // NB: the following comment is from the original code.
            // TODO: handle case where you get 2+ alt sequences back;  currently get list of 1 element
            // loop structure implemented to handle this, but doesn't really do anything currently.

            let var_ps: Result<Vec<_>, Error> = builder
                .build_altseq()?
                .into_iter()
                .map(|alt_data| {
                    let builder = AltSeqToHgvsp::new(reference_data.clone(), alt_data);
                    builder.build_hgvsp()
                })
                .collect();
            let var_p = var_ps?
                .into_iter()
                .next()
                .ok_or(Error::ProtVariantConstructionFailed)?;

            let var_p = if let HgvsVariant::ProtVariant {
                accession,
                gene_symbol,
                loc_edit,
            } = var_p
            {
                HgvsVariant::ProtVariant {
                    gene_symbol: self
                        .fetch_gene_symbol(accession.deref().as_str(), &gene_symbol)?,
                    accession,
                    loc_edit,
                }
            } else {
                return Err(Error::NotProtVariant);
            };

            Ok(var_p)
        } else {
            Err(Error::ExpectedCdsVariant(format!("{}", &var_c)))
        }
    }

    fn get_altered_sequence(
        &self,
        strand: i16,
        interval: Range<i32>,
        var: &HgvsVariant,
    ) -> Result<String, Error> {
        let mut seq = self.provider.as_ref().get_seq_part(
            var.accession(),
            Some(
                interval
                    .start
                    .try_into()
                    .map_err(|_e| Error::CannotConvertIntervalStart(interval.start))?,
            ),
            Some(
                interval
                    .end
                    .try_into()
                    .map_err(|_e| Error::CannotConvertIntervalEnd(interval.end))?,
            ),
        )?;

        let r = var
            .loc_range()
            .ok_or(Error::NoAlteredSequenceForMissingPositions)?;
        let r = ((r.start - interval.start) as usize)..((r.end - interval.start) as usize);

        let na_edit = var.na_edit().ok_or(Error::NaEditMissing)?;

        match na_edit {
            NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => {
                seq.replace_range(r, alternative)
            }
            NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => seq.replace_range(r, ""),
            NaEdit::Ins { alternative } => {
                seq.replace_range((r.start + 1)..(r.start + 1), alternative)
            }
            NaEdit::Dup { .. } => {
                let seg = seq[r.clone()].to_string();
                seq.replace_range(r.end..r.end, &seg);
            }
            NaEdit::InvRef { .. } | NaEdit::InvNum { .. } => {
                let rc = revcomp(&seq[r.clone()]);
                seq.replace_range(r, &rc);
            }
        }

        Ok(if strand == -1 { revcomp(&seq) } else { seq })
    }

    /// Convert an edit from one type to another, based on the strand and type.
    fn convert_edit_check_strand(
        &self,
        strand: i16,
        edit: &Mu<NaEdit>,
    ) -> Result<Mu<NaEdit>, Error> {
        let result = if strand == 1 {
            edit.inner().clone()
        } else {
            match edit.inner() {
                NaEdit::RefAlt {
                    reference,
                    alternative,
                } => NaEdit::RefAlt {
                    reference: revcomp(reference),
                    alternative: revcomp(alternative),
                },
                NaEdit::NumAlt { count, alternative } => NaEdit::NumAlt {
                    count: *count,
                    alternative: revcomp(alternative),
                },
                NaEdit::DelRef { reference } => NaEdit::DelRef {
                    reference: revcomp(reference),
                },
                NaEdit::Ins { alternative } => NaEdit::Ins {
                    alternative: revcomp(alternative),
                },
                NaEdit::Dup { reference } => NaEdit::Dup {
                    reference: revcomp(reference),
                },
                NaEdit::InvRef { reference } => NaEdit::InvRef {
                    reference: revcomp(reference),
                },
                NaEdit::DelNum { count } => NaEdit::DelNum { count: *count },
                NaEdit::InvNum { count } => NaEdit::InvNum { count: *count },
            }
        };
        Ok(Mu::from(result, edit.is_certain()))
    }

    /// Fetch reference sequence for variant and return updated `HgvsVariant` if necessary.
    pub fn replace_reference(&self, var: HgvsVariant) -> Result<HgvsVariant, Error> {
        match &var {
            HgvsVariant::ProtVariant { .. } => Err(Error::CannotUpdateReference),
            _ => Ok(()),
        }?;

        if let Some(NaEdit::Ins { .. }) = var.na_edit() {
            // Insertions have no reference sequence (zero-width); return as-is.
            return Ok(var);
        }

        if var.spans_intron() {
            debug!(
                "Can't update reference sequence for intronic variant {}",
                var
            );
            return Ok(var);
        }

        // For c. variants, we need coordinates on underlying sequence.
        let (r, ac): (Range<_>, _) = match &var {
            HgvsVariant::CdsVariant {
                accession,
                loc_edit,
                ..
            } => {
                let mapper = self.build_alignment_mapper(accession, accession, "transcript")?;
                (mapper.c_to_n(loc_edit.loc.inner())?.into(), accession)
            }
            HgvsVariant::GenomeVariant {
                accession,
                loc_edit,
                ..
            } => (loc_edit.loc.inner().clone().try_into()?, accession),
            HgvsVariant::MtVariant {
                accession,
                loc_edit,
                ..
            } => (loc_edit.loc.inner().clone().try_into()?, accession),
            HgvsVariant::TxVariant {
                accession,
                loc_edit,
                ..
            } => (loc_edit.loc.inner().clone().into(), accession),
            HgvsVariant::RnaVariant {
                accession,
                loc_edit,
                ..
            } => (loc_edit.loc.inner().clone().into(), accession),
            _ => panic!("Cases excluded above; cannot happen"),
        };

        // NB: The following comment is from the original code (no strict_bounds there either).
        // When strict_bounds is false and an error occurs, return variant as-is.
        if r.start < 0 {
            // This is an out-of-bounds variant.
            return Ok(var);
        }
        log::debug!("get_seq_part({}, {}, {})", ac, r.start, r.end);
        let seq = self.provider.as_ref().get_seq_part(
            ac,
            Some(r.start as usize),
            Some(r.end as usize),
        )?;
        if seq.len() != r.len() {
            // Tried to read beyond seq end; this is an out-of-bounds variant.
            log::debug!("Bailing out on out-of-bounds variant: {}", &var);
            return Ok(var);
        }

        let na_edit = var
            .na_edit()
            .expect("Variant must be of nucleic acid type here");
        if !na_edit.reference_equals(&seq) {
            Ok(var.with_reference(seq))
        } else {
            Ok(var)
        }
    }

    fn fetch_gene_symbol(
        &self,
        tx_ac: &str,
        gene_symbol: &Option<GeneSymbol>,
    ) -> Result<Option<GeneSymbol>, Error> {
        if !self.config.add_gene_symbol {
            Ok(gene_symbol.clone())
        } else if let Some(gene_symbol) = gene_symbol {
            Ok(Some(gene_symbol.clone()))
        } else {
            let hgnc = self.provider.as_ref().get_tx_identity_info(tx_ac)?.hgnc;
            if hgnc.is_empty() {
                Ok(None)
            } else {
                Ok(Some(GeneSymbol::from(hgnc)))
            }
        }
    }
}

#[cfg(test)]
mod test {
    use anyhow::Error;
    use pretty_assertions::assert_eq;
    use regex::Regex;
    use std::{
        path::{Path, PathBuf},
        str::FromStr,
    };
    use test_log::test;

    use crate::{
        data::uta_sr::test_helpers::build_provider,
        parser::{HgvsVariant, NoRef},
    };

    use super::{Config, Mapper};

    fn build_mapper() -> Result<Mapper, Error> {
        let provider = build_provider()?;
        let config = Config::default();
        Ok(Mapper::new(&config, provider))
    }

    #[test]
    fn fail_for_invalid_variant_types() -> Result<(), Error> {
        let mapper = build_mapper()?;

        let hgvs_g = "NC_000007.13:g.36561662C>T";
        let hgvs_c = "NM_001637.3:c.1582G>A"; // gene AOAH

        let var_g = HgvsVariant::from_str(hgvs_g)?;
        let var_c = HgvsVariant::from_str(hgvs_c)?;

        assert!(mapper.g_to_c(&var_c, "NM_001637.3", "splign").is_err());
        assert!(mapper.g_to_t(&var_c, "NM_001637.3", "splign").is_err());
        assert!(mapper.n_to_g(&var_c, "NM_001637.3", "splign").is_err());
        assert!(mapper.c_to_g(&var_g, "NM_001637.3", "splign").is_err());
        assert!(mapper.t_to_g(&var_g, "NM_001637.3", "splign").is_err());
        assert!(mapper.c_to_n(&var_g).is_err());
        assert!(mapper.n_to_c(&var_g).is_err());
        assert!(mapper.c_to_p(&var_g, None).is_err());

        Ok(())
    }

    #[test]
    fn fail_c_to_p_on_invalid_nm_accession() -> Result<(), Error> {
        let mapper = build_mapper()?;

        let hgvs_g = "NC_000007.13:g.36561662C>T";
        let var_g = HgvsVariant::from_str(hgvs_g)?;

        assert!(mapper.c_to_p(&var_g, Some("NM_999999.1")).is_err());

        Ok(())
    }

    #[test]
    fn fail_on_undefined_cds() -> Result<(), Error> {
        let mapper = build_mapper()?;

        let hgvs_n = "NR_111984.1:n.44G>A"; // legit
        let hgvs_c = "NR_111984.1:c.44G>A"; // bogus: c. with non-coding tx accession

        let var_n = HgvsVariant::from_str(hgvs_n)?;
        let var_c = HgvsVariant::from_str(hgvs_c)?;

        // n_to_c: transcript is non-coding
        assert!(mapper.n_to_c(&var_n).is_err());

        // c_to_n: var_c is bogus
        assert!(mapper.c_to_n(&var_c).is_err());

        Ok(())
    }

    #[test]
    fn map_var_of_unsupported_validation() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_003777.3:c.13552_*36del57"; // gene DNAH11
        let var_c = HgvsVariant::from_str(hgvs_c)?;

        let var_g = mapper.c_to_g(&var_c, "NC_000007.13", "splign")?;
        assert_eq!(
            format!("{}", &NoRef(&var_g)),
            "NC_000007.13:g.21940852_21940908del"
        );

        Ok(())
    }

    #[test]
    fn map_to_unknown_p_effect() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_020975.4:c.625+9C>T"; // gene RET
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_066124.1:p.?");

        Ok(())
    }

    // TODO(#17): Need to implement validation.
    // #[test]
    // fn map_of_c_out_of_cds_bound() -> Result<(), Error> {
    //     let mapper = build_mapper()?;
    //     let hgvs_c = "NM_145901.2:c.343T>C"; // gene HMGA1
    //     let var_c = HgvsVariant::from_str(hgvs_c)?;
    //     assert!(mapper.c_to_p(&var_c, None).is_err());

    //     Ok(())
    // }

    #[test]
    fn map_of_dup_at_cds_end() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_001051.2:c.1257dupG"; // gene SSTR3
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_001042.1:p.=");

        Ok(())
    }

    // TODO(#17): Need to implement validation.
    // #[test]
    // fn map_of_c_out_of_reference_bound() -> Result<(), Error> {
    //     let mapper = build_mapper()?;
    //     let hgvs_c = "NM_000249.3:c.-73960_*46597del"; // gene MLH1
    //     let var_c = HgvsVariant::from_str(hgvs_c)?;
    //     assert!(mapper.c_to_p(&var_c, None).is_err());

    //     Ok(())
    // }

    /// The following tests corresponds to the `test_hgvs_variantmapper_cp_sanity.py`
    /// test suite of the Python package.  It uses a mock data provider, defined
    /// in the `sanity_mock` module.

    mod sanity_mock {
        use anyhow::Error;
        use std::{
            path::{Path, PathBuf},
            rc::Rc,
        };

        use crate::data::interface::Provider as ProviderInterface;
        use crate::{
            data::interface::TxIdentityInfo,
            mapper::variant::{Config, Mapper},
        };

        #[derive(Debug, serde::Deserialize)]
        struct ProviderRecord {
            pub accession: String,
            pub transcript_sequence: String,
            pub cds_start_i: i32,
            pub cds_end_i: i32,
        }

        pub struct Provider {
            records: Vec<ProviderRecord>,
        }

        impl Provider {
            pub fn new(path: &Path) -> Result<Self, Error> {
                let mut records = Vec::new();

                let mut rdr = csv::ReaderBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(true)
                    .from_path(path)?;
                for record in rdr.deserialize() {
                    records.push(record?);
                }

                Ok(Self { records })
            }
        }

        impl ProviderInterface for Provider {
            fn data_version(&self) -> &str {
                panic!("for test use only");
            }

            fn schema_version(&self) -> &str {
                panic!("for test use only");
            }

            fn get_assembly_map(
                &self,
                _assembly: crate::static_data::Assembly,
            ) -> linked_hash_map::LinkedHashMap<String, String> {
                panic!("for test use only");
            }

            fn get_gene_info(
                &self,
                _hgnc: &str,
            ) -> Result<crate::data::interface::GeneInfoRecord, crate::data::error::Error>
            {
                panic!("for test use only");
            }

            fn get_pro_ac_for_tx_ac(
                &self,
                _tx_ac: &str,
            ) -> Result<Option<String>, crate::data::error::Error> {
                panic!("for test use only");
            }

            fn get_seq_part(
                &self,
                tx_ac: &str,
                begin: Option<usize>,
                end: Option<usize>,
            ) -> Result<String, crate::data::error::Error> {
                for record in &self.records {
                    if record.accession == tx_ac {
                        let seq = &record.transcript_sequence;
                        return match (begin, end) {
                            (None, None) => Ok(seq.to_string()),
                            (None, Some(end)) => Ok(seq[..end].to_string()),
                            (Some(begin), None) => Ok(seq[begin..].to_string()),
                            (Some(begin), Some(end)) => Ok(seq[begin..end].to_string()),
                        };
                    }
                }
                Err(crate::data::error::Error::NoSequenceRecord(
                    tx_ac.to_string(),
                ))
            }

            fn get_acs_for_protein_seq(
                &self,
                _seq: &str,
            ) -> Result<Vec<String>, crate::data::error::Error> {
                panic!("for test use only");
            }

            fn get_similar_transcripts(
                &self,
                _tx_ac: &str,
            ) -> Result<Vec<crate::data::interface::TxSimilarityRecord>, crate::data::error::Error>
            {
                panic!("for test use only");
            }

            fn get_tx_exons(
                &self,
                _tx_ac: &str,
                _alt_ac: &str,
                _alt_aln_method: &str,
            ) -> Result<Vec<crate::data::interface::TxExonsRecord>, crate::data::error::Error>
            {
                todo!()
            }

            fn get_tx_for_gene(
                &self,
                _gene: &str,
            ) -> Result<Vec<crate::data::interface::TxInfoRecord>, crate::data::error::Error>
            {
                panic!("for test use only");
            }

            fn get_tx_for_region(
                &self,
                _alt_ac: &str,
                _alt_aln_method: &str,
                _start_i: i32,
                _end_i: i32,
            ) -> Result<Vec<crate::data::interface::TxForRegionRecord>, crate::data::error::Error>
            {
                panic!("for test use only");
            }

            fn get_tx_identity_info(
                &self,
                tx_ac: &str,
            ) -> Result<TxIdentityInfo, crate::data::error::Error> {
                for record in &self.records {
                    if record.accession == tx_ac {
                        return Ok(TxIdentityInfo {
                            tx_ac: record.accession.clone(),
                            alt_ac: record.accession.clone(),
                            alt_aln_method: "splign".to_string(),
                            cds_start_i: record.cds_start_i,
                            cds_end_i: record.cds_end_i,
                            lengths: Vec::new(),
                            hgnc: "MOCK".to_string(),
                        });
                    }
                }
                Err(crate::data::error::Error::NoSequenceRecord(
                    tx_ac.to_string(),
                ))
            }

            fn get_tx_info(
                &self,
                _tx_ac: &str,
                _alt_ac: &str,
                _alt_aln_method: &str,
            ) -> Result<crate::data::interface::TxInfoRecord, crate::data::error::Error>
            {
                panic!("for test use only");
            }

            fn get_tx_mapping_options(
                &self,
                _tx_ac: &str,
            ) -> Result<
                Vec<crate::data::interface::TxMappingOptionsRecord>,
                crate::data::error::Error,
            > {
                panic!("for test use only");
            }
        }

        pub fn build_mapper(strict_bounds: bool) -> Result<Mapper, Error> {
            let path = PathBuf::from("tests/data/mapper/sanity_cp.tsv");
            let provider = Rc::new(Provider::new(&path)?);
            let config = Config {
                strict_bounds,
                ..Default::default()
            };
            Ok(Mapper::new(&config, provider))
        }
    }

    fn test_hgvs_c_to_p_conversion(hgvsc: &str, hgvsp_expected: &str) -> Result<(), Error> {
        let mapper = sanity_mock::build_mapper(false)?;

        let var_c = HgvsVariant::from_str(hgvsc)?;
        let ac_p = "MOCK";

        let var_p = mapper.c_to_p(&var_c, Some(ac_p))?;
        let hgvsp_actual = format!("{}", &var_p);

        assert_eq!(hgvsp_actual, hgvsp_expected);

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_silent() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.6A>G";
        let hgvsp_expected = "MOCK:p.Lys2=";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_substitution() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.6A>T";
        let hgvsp_expected = "MOCK:p.Lys2Asn";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_substitution_introduces_stop_codon() -> Result<(), Error> {
        let hgvsc = "NM_999996.1:c.8C>A";
        let hgvsp_expected = "MOCK:p.Ser3Ter";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_substitution_removes_stop_codon() -> Result<(), Error> {
        let hgvsc = "NM_999998.1:c.30G>T";
        let hgvsp_expected = "MOCK:p.Ter10TyrextTer3";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    //xx
    #[test]
    fn hgvs_c_to_p_insertion_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.6_7insGGG";
        let hgvsp_expected = "MOCK:p.Lys2_Ala3insGly";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_insertion_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.22_23insT";
        let hgvsp_expected = "MOCK:p.Ala8ValfsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_adds_stop() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.8_9insTT";
        let hgvsp_expected = "MOCK:p.Lys4Ter";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.10_12del";
        let hgvsp_expected = "MOCK:p.Lys4del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion2_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.4_15del";
        let hgvsp_expected = "MOCK:p.Lys2_Ala5del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion3_no_frameshift_c_term() -> Result<(), Error> {
        let hgvsc = "NM_999995.1:c.4_6del";
        let hgvsp_expected = "MOCK:p.Lys3del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion4_no_frameshift_c_term() -> Result<(), Error> {
        let hgvsc = "NM_999994.1:c.4_9del";
        let hgvsp_expected = "MOCK:p.Lys3_Lys4del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion5_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999994.1:c.20_25del";
        let hgvsp_expected = "MOCK:p.Ala7_Arg9delinsGly";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion6_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.5_7del";
        let hgvsp_expected = "MOCK:p.Lys2_Ala3delinsThr";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion7_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999993.1:c.13_24del";
        let hgvsp_expected = "MOCK:p.Arg5_Ala8del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_frameshift_nostop() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.11_12del";
        let hgvsp_expected = "MOCK:p.Lys4SerfsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_frameshift_adds_stop() -> Result<(), Error> {
        let hgvsc = "NM_999997.1:c.7del";
        let hgvsp_expected = "MOCK:p.Ala3ArgfsTer6";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_no_frameshift_removes_stop_plus_previous() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.25_30del";
        let hgvsp_expected = "MOCK:p.Lys9_Ter10delinsGly";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_indel_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.11_12delinsTCCCA";
        let hgvsp_expected = "MOCK:p.Lys4delinsIlePro";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_indel2_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.11_18delinsTCCCA";
        let hgvsp_expected = "MOCK:p.Lys4_Phe6delinsIlePro";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_indel_frameshift_nostop() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.8delinsGG";
        let hgvsp_expected = "MOCK:p.Ala3GlyfsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_1aa_no_frameshift_2() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.10_12dup";
        let hgvsp_expected = "MOCK:p.Lys4dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_1aa_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.16_18dup";
        let hgvsp_expected = "MOCK:p.Phe6dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_2aa_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.16_21dup";
        let hgvsp_expected = "MOCK:p.Phe6_Arg7dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_2aa2_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999995.1:c.4_6dup";
        let hgvsp_expected = "MOCK:p.Lys3dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_3aa_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.16_24dup";
        let hgvsp_expected = "MOCK:p.Phe6_Ala8dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.12_13dup";
        let hgvsp_expected = "MOCK:p.Ala5GlufsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_intron() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.12+1G>A";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_five_prime_utr() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.-2A>G";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_three_prime_utr() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.*3G>A";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_into_three_prime_utr_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.27_*3del";
        let hgvsp_expected = "MOCK:p.Lys9XaafsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_into_three_prime_utr_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999995.1:c.28_*3del";
        let hgvsp_expected = "MOCK:p.Lys10_Ter11delinsArgGlnPheArg";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_delins_into_three_prime_utr_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999995.1:c.28_*3delinsGGG";
        let hgvsp_expected = "MOCK:p.Lys10_Ter11delinsGlyArgGlnPheArg";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    /// See recommendations re p.? (p.Met1?) at:
    /// http://varnomen.hgvs.org/recommendations/protein/variant/substitution/
    #[test]
    fn hgvs_c_to_p_substitution_removes_start_codon() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.1A>G";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_from_five_prime_utr_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.-3_1del";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_from_five_prime_utr_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.-3_3del";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_delins_from_five_prime_utr_no_frameshift() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.-3_3delinsAAA";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_delete_entire_gene() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.-3_*1del";
        let hgvsp_expected = "MOCK:p.0?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    /// Check the case with multiple stop codons.  We introduced a change in hgvs-rs
    /// that does not handle multiple stop codons in the transcript sequence as
    /// conservatively as the Python version.
    #[test]
    fn hgvs_c_to_p_multiple_stop_codons() -> Result<(), Error> {
        let hgvsc = "NM_999992.1:c.4G>A";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    // The following tests correspond to the tests in `test_hgvs_variantmapper_cp_real.py`.
    //
    // For adding tests, you will have to
    //
    // - add a record to `real_cp.tsv`
    // - update `bootstrap.sh` with the HGNC symbol if necessary
    // - re-run `bootstrap.sh` so the records are pulled into the subset
    // - re-create the local database and import the subset
    // - re-run the test with `TEST_SEQREPO_CACHE_MODE=write` so the relevant queries to
    //   the seqrepo are cached

    #[test]
    fn hgvs_c_to_p_format() -> Result<(), Error> {
        let mapper = build_mapper()?;
        // gene SIL1
        let hgvs_c = "NM_022464.4:c.3G>A";
        // let hgvsp_expected_alternative = "NP_071909.1:p.?";

        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_071909.1:p.Met1?");

        // TODO(#25): implement formatting of display and uncomment
        // alt_format_p = var_p.format(conf={"p_init_met": False})
        // self.assertEqual(hgvsp_expected_alternative, alt_format_p)

        Ok(())
    }

    mod gcp_tests {
        use anyhow::Error;
        use std::path::Path;

        #[derive(Debug, serde::Deserialize)]
        pub struct Record {
            pub id: String,
            #[serde(alias = "HGVSg")]
            pub hgvs_g: String,
            #[serde(alias = "HGVSc")]
            pub hgvs_c: String,
            #[serde(alias = "HGVSp")]
            pub hgvs_p: Option<String>,
            pub description: Option<String>,
            pub alternatives: Option<String>,
        }

        pub fn load_records(path: &Path) -> Result<Vec<Record>, Error> {
            let mut records = Vec::new();

            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .flexible(true)
                .comment(Some(b'#'))
                .from_path(path)?;
            for record in rdr.deserialize() {
                let mut record: Record = record?;
                // p.(*) => p.
                record.hgvs_p = record.hgvs_p.map(|s| s.replace(['(', ')'], ""));
                records.push(record);
            }

            Ok(records)
        }
    }

    #[test]
    fn cp_real() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let path = PathBuf::from("tests/data/mapper/real_cp.tsv");
        let records = gcp_tests::load_records(&path)?;

        for record in records {
            let var_c = HgvsVariant::from_str(&record.hgvs_c)?;
            let prot_ac = record
                .hgvs_p
                .as_ref()
                .expect("problem with result in test")
                .split(':')
                .next()
                .map(|s| s.to_string());
            let var_p = mapper.c_to_p(&var_c, prot_ac.as_deref())?;
            let result = format!("{}", &var_p);
            let expected = &record.hgvs_p.expect("problem with result in test");

            let expected = if &result != expected {
                expected.replace('*', "Ter")
            } else {
                expected.clone()
            };
            assert_eq!(result, expected);
        }

        Ok(())
    }

    // The following tests correspond to those in `test_hgvs_variantmapper_gcp.py`.

    fn run_gxp_test(path: &str, noref: bool) -> Result<(), Error> {
        fn rm_del_seq(var: &HgvsVariant, noref: bool) -> String {
            let tmp = if noref {
                format!("{}", &NoRef(var))
            } else {
                format!("{var}")
            };
            let re = Regex::new(r"del\w+ins").expect("problem with regex in test");
            re.replace(&tmp, "delins").to_string()
        }

        let mapper = build_mapper()?;
        let records = gcp_tests::load_records(Path::new(path))?;

        for record in &records {
            let var_g = HgvsVariant::from_str(&record.hgvs_g)?;
            let var_x = HgvsVariant::from_str(&record.hgvs_c)?;
            let var_p = record
                .hgvs_p
                .as_ref()
                .map(|s| HgvsVariant::from_str(s))
                .transpose()?;

            // g -> x
            let var_x_test = match &var_x {
                HgvsVariant::CdsVariant { accession, .. } => {
                    mapper.g_to_c(&var_g, accession, "splign")?
                }
                HgvsVariant::TxVariant { accession, .. } => {
                    mapper.g_to_n(&var_g, accession, "splign")?
                }
                _ => panic!("cannot happen"),
            };

            // Use `del<COUNT>` syntax in output when we saw this in the input.  The original
            // Python library implements this by always storing the count in the nucleic acid
            // edit.
            let var_x_test = if var_x.is_na_edit_num() {
                var_x_test.with_na_ref_num()
            } else {
                var_x_test
            };

            assert_eq!(
                rm_del_seq(&var_x, noref),
                rm_del_seq(&var_x_test, noref),
                "{} != {} (g>t; {}; HGVSg={})",
                var_x,
                var_x_test,
                &record.id,
                &record.hgvs_g
            );

            // c, n -> g
            let var_g_test = match &var_x {
                HgvsVariant::CdsVariant { .. } => {
                    mapper.c_to_g(&var_x, var_g.accession(), "splign")?
                }
                HgvsVariant::TxVariant { .. } => {
                    mapper.n_to_g(&var_x, var_g.accession(), "splign")?
                }
                _ => panic!("cannot happen"),
            };

            // Use `del<COUNT>` syntax in output when we saw this in the input.  The original
            // Python library implements this by always storing the count in the nucleic acid
            // edit.
            let var_g_test = if var_g.is_na_edit_num() {
                var_g_test.with_na_ref_num()
            } else {
                var_g_test
            };

            assert_eq!(
                rm_del_seq(&var_g, noref),
                rm_del_seq(&var_g_test, noref),
                "{} != {} (t>g; {}; HGVSc={})",
                var_g,
                var_g_test,
                &record.id,
                &record.hgvs_c
            );

            if let Some(var_p) = &var_p {
                // c -> p
                let hgvs_p_exp = format!("{var_p}");
                let var_p_test = mapper.c_to_p(&var_x, Some(var_p.accession()))?;

                // TODO: if expected value isn't uncertain, strip uncertain from test
                // if var_p.posedit and not var_p.posedit.uncertain:
                //     # if expected value isn't uncertain, strip uncertain from test
                //     var_p_test.posedit.uncertain = False

                let mut hgvs_p_test = format!("{}", &var_p_test);

                if hgvs_p_exp.ends_with("Ter") {
                    let re = Regex::new(r"Ter\d+$").expect("problem with regex in test");
                    hgvs_p_test = re.replace(&hgvs_p_test, "Ter").to_string();
                }

                assert_eq!(
                    hgvs_p_exp, hgvs_p_test,
                    "{} != {} ({})",
                    &hgvs_p_exp, &hgvs_p_test, &record.id,
                );
            }
        }

        Ok(())
    }

    #[test]
    fn zcchc3_dbsnp() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/ZCCHC3-dbSNP.tsv", false)
    }

    #[test]
    fn orai1_dbsnp() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/ORAI1-dbSNP.tsv", false)
    }

    #[test]
    fn folr3_dbsnp() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/FOLR3-dbSNP.tsv", false)
    }

    #[test]
    fn adra2b_dbsnp() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/ADRA2B-dbSNP.tsv", false)
    }

    #[test]
    fn jrk_dbsnp() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/JRK-dbSNP.tsv", false)
    }

    #[test]
    fn nefl_dbsnp() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/NEFL-dbSNP.tsv", false)
    }

    #[test]
    fn dnah11_hgmd() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/DNAH11-HGMD.tsv", true)
    }

    #[test]
    fn dnah11_dbsnp_nm_003777() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/DNAH11-dbSNP-NM_003777.tsv", false)
    }

    #[test]
    fn dnah11_db_snp_nm_001277115() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/DNAH11-dbSNP-NM_001277115.tsv", false)
    }

    #[test]
    fn regression() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/regression.tsv", false)
    }

    #[ignore]
    #[test]
    fn dnah11_db_snp_full() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/DNAH11-dbSNP.tsv", false)
    }

    #[test]
    fn real() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/real.tsv", false)
    }

    /// Check for issues with variants affecting `Met1` leading to `p.Met1?`.
    #[test]
    fn real_met1() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/real-met1.tsv", false)
    }

    #[test]
    fn noncoding() -> Result<(), Error> {
        run_gxp_test("tests/data/mapper/gcp/noncoding.tsv", false)
    }

    // #[test]
    // fn case() -> Result<(), Error> {
    //     let mapper = build_mapper()?;

    //     let s_c = "NM_000425.3:c.3772dupT";
    //     let s_p = "NP_000416.1:p.Ter1258Leuext*96";

    //     let var_c = HgvsVariant::from_str(s_c)?;
    //     let var_p = mapper.c_to_p(&var_c, None)?;

    //     let hgvsp_actual = format!("{}", &var_p);
    //     assert_eq!(hgvsp_actual, s_p);

    //     Ok(())
    // }
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
