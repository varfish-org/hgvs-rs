//! Code for mapping variants between sequences.

use log::{debug, info};
use quick_cache::{unsync::Cache, Equivalent};
use std::cell::RefCell;
use std::ops::Deref;
use std::{ops::Range, sync::Arc};

use crate::mapper::alignment;
use crate::mapper::altseq::{ref_transcript_data_shared, AltSeqBuilder, VariantLocation};
use crate::{
    data::interface::Provider,
    mapper::Error,
    normalizer::{self, Direction, Normalizer},
    parser::{
        Accession, CdsInterval, CdsLocEdit, CdsPos, GeneSymbol, GenomeInterval, GenomeLocEdit,
        HgvsVariant, Mu, NaEdit, TxInterval, TxLocEdit, TxPos,
    },
    sequences::revcomp,
    validator::{ValidationLevel, Validator},
};

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
    pub shuffle_direction: Direction,
    pub window_size: usize,
    /// Treat insertions and duplications that touch an exon boundary as intronic in `c_to_p`,
    /// giving `p.?`.  With `false`, e.g., `c.26_26+1insC` and `c.26+1dup` change the protein.
    /// Before c.1 and after the stop codon, they change no base of the CDS and give `p.?` like
    /// a variant in the UTR.
    pub ins_at_boundary_is_intronic: bool,
    /// Let `c_to_p` move an insertion or duplication in an exon or intron to the other side
    /// of an exon boundary, as preferred by `shift_over_boundary_preference`.  The shifts need
    /// the genome alignment.  `assembly::Mapper::c_to_p()` provides it if `genome_seq_available`
    /// is set.  Without it, `c_to_p` fails with `Error::ShiftOverBoundaryNeedsAlignment` for
    /// such variants, unless the preference is `Default`.
    pub shift_over_boundary: bool,
    /// See `shift_over_boundary`.
    pub shift_over_boundary_preference: ShiftOverBoundaryPreference,
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
            shuffle_direction: Default::default(),
            window_size: 20,
            ins_at_boundary_is_intronic: true,
            shift_over_boundary: false,
            shift_over_boundary_preference: Default::default(),
        }
    }
}

/// Where `c_to_p` with `Config::shift_over_boundary` puts an insertion or duplication that fits
/// on both sides of an exon boundary.
///
/// `c_to_p` shifts the variant fully 3' and fully 5' on the genome.  If one of these shifts
/// lies on the other side of an exon boundary, `c_to_p` uses it for the protein if it is in
/// the preferred region.  The c. variant itself does not change.
#[derive(Debug, PartialEq, Eq, Clone, Copy, Default)]
#[non_exhaustive]
pub enum ShiftOverBoundaryPreference {
    /// Keep the variant where the HGVS rules put it.
    #[default]
    Default,
    /// Prefer the intron, which gives `p.?`.
    Intron,
    /// Prefer the exon, which gives a protein change.
    Exon,
}

/// Projects variants between sequences using `alignment::Mapper`.
pub struct Mapper {
    pub config: Config,
    provider: Arc<dyn Provider + Send + Sync>,
    validator: Arc<dyn Validator + Send + Sync>,
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
    pub fn new(config: &Config, provider: Arc<dyn Provider + Send + Sync>) -> Mapper {
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
    pub fn provider(&self) -> Arc<dyn Provider + Send + Sync> {
        self.provider.clone()
    }

    /// Obtain new `alignment::Mapper` for the given arguments, possibly caching results.
    fn build_alignment_mapper(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Arc<alignment::Mapper>, Error> {
        build_alignment_mapper_cached(
            &self.provider,
            self.config.strict_bounds,
            tx_ac,
            alt_ac,
            alt_aln_method,
        )
    }

    /// Construct a new normalizer for the variant mapper.
    pub fn normalizer(&self) -> Result<Normalizer<'_>, Error> {
        Ok(Normalizer::new(
            self,
            self.provider.clone(),
            self.validator.clone(),
            normalizer::Config {
                replace_reference: self.config.replace_reference,
                shuffle_direction: self.config.shuffle_direction,
                window_size: self.config.window_size,
                ..Default::default()
            },
        ))
    }

    /// Construct a new normalizer for the variant mapper.
    pub fn left_normalizer(&self) -> Result<Normalizer<'_>, Error> {
        Ok(Normalizer::new(
            self,
            self.provider.clone(),
            self.validator.clone(),
            normalizer::Config {
                replace_reference: self.config.replace_reference,
                shuffle_direction: Direction::ThreeToFive,
                ..Default::default()
            },
        ))
    }

    /// Construct a new normalizer for the variant mapper.
    pub fn right_normalizer(&self) -> Result<Normalizer<'_>, Error> {
        Ok(Normalizer::new(
            self,
            self.provider.clone(),
            self.validator.clone(),
            normalizer::Config {
                replace_reference: self.config.replace_reference,
                shuffle_direction: Direction::FiveToThree,
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

        let (mapper, var_g) = match &var_g {
            HgvsVariant::GenomeVariant {
                accession,
                loc_edit,
                gene_symbol: _,
            } => {
                let mapper =
                    self.build_alignment_mapper(tx_ac, &accession.value, alt_aln_method)?;
                if mapper.strand == -1
                    && !self.config.strict_bounds
                    && !mapper.is_g_interval_in_bounds(loc_edit.loc.inner())
                    && self.config.renormalize_g
                {
                    info!("Renormalizing out-of-bounds minus strand variant on genomic sequence");
                    (mapper, self.left_normalizer()?.normalize(&var_g)?)
                } else {
                    (mapper, var_g)
                }
            }
            _ => unreachable!(),
        };

        if let HgvsVariant::GenomeVariant {
            accession: _,
            loc_edit,
            gene_symbol,
        } = &var_g
        {
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
                        && distance_without_zero(pos_n.start.base, pos_n.end.base) > 1
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
                // This is how the original code handles uncertain positions.  We will reach
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
            Err(Error::ExpectedGenomeVariant(format!("{}", var_g)))
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
            Err(Error::ExpectedTxVariant(format!("{}", var_n)))
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
                        && pos_c.start.cds_from == pos_c.end.cds_from
                        && distance_without_zero(pos_c.start.base, pos_c.end.base) > 1
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
            Err(Error::ExpectedGenomeVariant(format!("{}", var_g)))
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
                                .ok_or(Error::NoNAEditInHgvsC(format!("{}", var_c)))?
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
            Err(Error::ExpectedCdsVariant(format!("{}", var_c)))
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
            _ => Err(Error::ExpectedCdsVariant(format!("{}", var_t))),
        }
    }

    /// Convert from CDS variant (c.) to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::CdsVariant` to project
    pub fn c_to_n(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, Error> {
        tracing::debug!("c_to_n({})", var_c);
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

            tracing::debug!("c_to_n({}) = {}", var_c, &var_n);
            Ok(var_n)
        } else {
            Err(Error::ExpectedCdsVariant(format!("{}", var_c)))
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
            Err(Error::ExpectedTxVariant(format!("{}", var_n)))
        }
    }

    /// Convert from CDS variant (c.) to protein variant (p.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::TxVariant` to project
    /// * `pro_ac` -- Protein accession
    pub fn c_to_p(&self, var_c: &HgvsVariant, prot_ac: Option<&str>) -> Result<HgvsVariant, Error> {
        self.c_to_p_with_alt(var_c, prot_ac, None)
    }

    /// Like `c_to_p()`, but takes the reference of intron bases from the genome alignment
    /// `alt`, given as `(alt_ac, alt_aln_method)`.
    pub(crate) fn c_to_p_with_alt(
        &self,
        var_c: &HgvsVariant,
        prot_ac: Option<&str>,
        alt: Option<(&str, &str)>,
    ) -> Result<HgvsVariant, Error> {
        use super::altseq::*;

        if let HgvsVariant::CdsVariant {
            accession,
            gene_symbol: _,
            loc_edit: _,
        } = &var_c
        {
            self.validator.validate(var_c)?;

            let var_c = if self.config.replace_reference {
                self.replace_reference_with_alt(var_c.clone(), alt)?
            } else {
                var_c.clone()
            };

            let reference_data = ref_transcript_data_shared(
                &self.provider,
                accession.deref(),
                prot_ac.map(|s| s.to_string()).as_deref(),
            )?;
            let builder = AltSeqBuilder {
                ins_at_boundary_is_intronic: self.config.ins_at_boundary_is_intronic,
                ..AltSeqBuilder::new(var_c, &reference_data)
            };
            let builder = self.shift_over_boundary(builder, alt)?;

            // NB: the following comment is from the original code.
            // TODO: handle case where you get 2+ alt sequences back;  currently get list of 1 element
            // loop structure implemented to handle this, but doesn't really do anything currently.

            let var_ps: Result<Vec<_>, Error> = builder
                .build_altseq()?
                .into_iter()
                .map(|alt_data| {
                    let builder = AltSeqToHgvsp::new(&reference_data, alt_data);
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
            Err(Error::ExpectedCdsVariant(format!("{}", var_c)))
        }
    }

    /// Apply `Config::shift_over_boundary` to the variant of `builder`, as `c_to_p` of
    /// biocommons/hgvs does.
    ///
    /// For an insertion or duplication in an exon or intron, take the first shift, 3' then
    /// 5', that lies on the other side of an exon boundary.  Return its builder if that side
    /// is the preferred one, else `builder`.  `alt` is the genome alignment as
    /// `(alt_ac, alt_aln_method)`.
    fn shift_over_boundary<'a>(
        &self,
        builder: AltSeqBuilder<'a>,
        alt: Option<(&str, &str)>,
    ) -> Result<AltSeqBuilder<'a>, Error> {
        let Some(preferred) = self.shift_preference(&builder.var_c)? else {
            return Ok(builder);
        };
        let region = builder.variant_location()?;
        if !matches!(region, VariantLocation::Intron | VariantLocation::Cds) {
            return Ok(builder);
        }
        let (alt_ac, alt_aln_method) =
            alt.ok_or_else(|| Error::ShiftOverBoundaryNeedsAlignment(builder.var_c.to_string()))?;

        let tx_ac = builder.var_c.accession().to_string();
        let var_g = self.c_to_g(&builder.var_c, alt_ac, alt_aln_method)?;
        let mut seen = vec![var_g.to_string()];
        for direction in [Direction::FiveToThree, Direction::ThreeToFive] {
            // Skip a shift that fails, unless data is missing, see `DIFFERENCES.md`.
            let shifted_var_g = match self.shift_var_g(&var_g, direction, alt_aln_method) {
                Ok(shifted_var_g) => shifted_var_g,
                Err(e) if is_data_error(&e) => return Err(e),
                Err(e) => {
                    debug!("Cannot shift {} over the boundary: {}", var_g, e);
                    continue;
                }
            };
            let key = shifted_var_g.to_string();
            if seen.contains(&key) {
                continue;
            }
            seen.push(key);
            let shifted_var_c = match self.g_to_c(&shifted_var_g, &tx_ac, alt_aln_method) {
                Ok(shifted_var_c) => shifted_var_c,
                Err(e) if is_data_error(&e) => return Err(e),
                Err(e) => {
                    debug!("Cannot map {} to {}: {}", shifted_var_g, tx_ac, e);
                    continue;
                }
            };

            let shifted = AltSeqBuilder {
                ins_at_boundary_is_intronic: self.config.ins_at_boundary_is_intronic,
                ..AltSeqBuilder::new(shifted_var_c, builder.reference_data)
            };
            let shifted_region = shifted.variant_location()?;
            if matches!(
                shifted_region,
                VariantLocation::Intron | VariantLocation::Cds
            ) && shifted_region != region
            {
                return Ok(if shifted_region == preferred {
                    shifted
                } else {
                    builder
                });
            }
        }
        Ok(builder)
    }

    /// The region that `c_to_p` prefers for a shift of `var_c` over an exon boundary, or `None`
    /// if it does not shift.  Only insertions and duplications shift.
    fn shift_preference(&self, var_c: &HgvsVariant) -> Result<Option<VariantLocation>, Error> {
        if !self.config.shift_over_boundary
            || !var_c
                .na_edit()
                .is_some_and(|edit| edit.is_ins() || edit.is_dup())
        {
            return Ok(None);
        }
        // `Default` keeps the variant, so it needs no shift and no genome, see `DIFFERENCES.md`.
        Ok(match self.config.shift_over_boundary_preference {
            ShiftOverBoundaryPreference::Default => None,
            ShiftOverBoundaryPreference::Intron => Some(VariantLocation::Intron),
            ShiftOverBoundaryPreference::Exon => Some(VariantLocation::Cds),
        })
    }

    /// Whether `c_to_p` shifts `var_c` over an exon boundary, which needs the genome alignment.
    /// This is the case for an insertion or duplication in an exon or intron if the preference
    /// is not `Default`.
    pub(crate) fn shift_needs_alignment(&self, var_c: &HgvsVariant) -> Result<bool, Error> {
        if self.shift_preference(var_c)?.is_none() {
            return Ok(false);
        }
        let reference_data =
            ref_transcript_data_shared(&self.provider, &var_c.accession().value, None)?;
        let region = AltSeqBuilder {
            ins_at_boundary_is_intronic: self.config.ins_at_boundary_is_intronic,
            ..AltSeqBuilder::new(var_c.clone(), &reference_data)
        }
        .variant_location()?;
        Ok(matches!(
            region,
            VariantLocation::Intron | VariantLocation::Cds
        ))
    }

    /// Shift `var_g` fully in `direction` and rewrite a duplication as an insertion, as
    /// `_var_g_shift_with_rewrite` of biocommons/hgvs does.
    ///
    /// The insertion goes to the end of the duplicated bases that the shift reached: after
    /// them for a 3' shift and before them for a 5' shift on the genome.  The dup is fully
    /// shifted, so the insertion cannot move further in `direction` either.
    ///
    /// biocommons/hgvs picks the end by the strand of the transcript instead, see
    /// `DIFFERENCES.md`.
    fn shift_var_g(
        &self,
        var_g: &HgvsVariant,
        direction: Direction,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, Error> {
        // As in biocommons/hgvs, the normalizer does not validate.
        let normalizer = Normalizer::new(
            self,
            self.provider.clone(),
            ValidationLevel::Null.validator(false, self.provider.clone()),
            normalizer::Config {
                alt_aln_method: alt_aln_method.to_string(),
                shuffle_direction: direction,
                window_size: self.config.window_size,
                ..Default::default()
            },
        );
        let var_g = normalizer.normalize(var_g)?;
        if !var_g.na_edit().is_some_and(NaEdit::is_dup) {
            return Ok(var_g);
        }

        let var_g = self.replace_reference(var_g)?;
        let HgvsVariant::GenomeVariant {
            accession,
            gene_symbol,
            loc_edit,
        } = &var_g
        else {
            return Err(Error::ExpectedGenomeVariant(var_g.to_string()));
        };
        let (Some(start), Some(end), NaEdit::Dup { reference }) = (
            loc_edit.loc.inner().start,
            loc_edit.loc.inner().end,
            loc_edit.edit.inner(),
        ) else {
            return Err(Error::MissingGenomeIntervalPosition(var_g.to_string()));
        };
        let (start, end) = match direction {
            Direction::FiveToThree => (end, end + 1),
            Direction::ThreeToFive => (start - 1, start),
        };
        Ok(HgvsVariant::GenomeVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: GenomeLocEdit {
                loc: Mu::Certain(GenomeInterval {
                    start: Some(start),
                    end: Some(end),
                }),
                edit: Mu::Certain(NaEdit::Ins {
                    alternative: reference.clone(),
                }),
            },
        })
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

        let (start, end) = if interval.start > r.start || interval.start > r.end {
            tracing::warn!(
                "Altered sequence range start {} is greater than variant range start {} or end {}, clamping. Variant description is {}",
                interval.start,
                r.end,
                r.start,
                &var
            );
            // FIXME: workaround for cases that would end up with negative start positions
            //  this needs more investigation
            let n = r.len();
            let mut start = r.start - interval.start;
            let mut end = r.end - interval.start;
            if start < 0 {
                start += n as i32;
            }
            if end < 0 {
                end += n as i32;
            }
            (start, end)
        } else {
            (r.start - interval.start, r.end - interval.start)
        };

        let start = usize::try_from(start).map_err(|_| Error::CannotConvertIntervalStart(start))?;
        let end = usize::try_from(end).map_err(|_| Error::CannotConvertIntervalEnd(end))?;
        let r = start..end;
        let r = if r.end > seq.len() {
            tracing::warn!(
                    "Altered sequence range {:?} is incompatible with sequence length {:?}, clamping. Variant description is {}",
                    r,
                    seq.len(),
                    &var
                );
            r.start..seq.len()
        } else {
            r
        };

        let na_edit = var.na_edit().ok_or(Error::NaEditMissing)?;

        match na_edit {
            NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => {
                seq.replace_range(r, alternative)
            }
            NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => {
                // FIXME the original code in python simply does `del seq[pos_start:pos_end]`,
                //  which does not error if `pos_end > len(seq)`. Check if this is intended or not.
                seq.replace_range(r, "")
            }
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
        self.replace_reference_with_alt(var, None)
    }

    /// Like `replace_reference()`, but takes the reference of an intronic c. variant from the
    /// genome alignment `alt`, given as `(alt_ac, alt_aln_method)`.
    fn replace_reference_with_alt(
        &self,
        var: HgvsVariant,
        alt: Option<(&str, &str)>,
    ) -> Result<HgvsVariant, Error> {
        match &var {
            HgvsVariant::ProtVariant { .. } => Err(Error::CannotUpdateReference),
            _ => Ok(()),
        }?;

        if let Some(NaEdit::Ins { .. }) = var.na_edit() {
            // Insertions have no reference sequence (zero-width); return as-is.
            return Ok(var);
        }

        let alt = alt.filter(|_| matches!(var, HgvsVariant::CdsVariant { .. }));
        if var.spans_intron() && alt.is_none() {
            debug!(
                "Can't update reference sequence for intronic variant {}",
                var
            );
            return Ok(var);
        }

        // For c. variants, we need coordinates on underlying sequence.
        let mut strand = 1;
        let (r, ac): (Range<_>, &str) = match &var {
            HgvsVariant::CdsVariant {
                accession,
                loc_edit,
                ..
            } => match alt {
                Some((alt_ac, alt_aln_method)) if var.spans_intron() => {
                    let mapper = self.build_alignment_mapper(accession, alt_ac, alt_aln_method)?;
                    strand = mapper.strand;
                    let pos_g = mapper.c_to_g(loc_edit.loc.inner())?;
                    (pos_g.inner().clone().try_into()?, alt_ac)
                }
                _ => {
                    let mapper = self.build_alignment_mapper(accession, accession, "transcript")?;
                    (mapper.c_to_n(loc_edit.loc.inner())?.into(), accession)
                }
            },
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
        tracing::trace!("get_seq_part({}, {}, {})", ac, r.start, r.end);
        let seq = self.provider.as_ref().get_seq_part(
            ac,
            Some(r.start as usize),
            Some(r.end as usize),
        )?;
        if seq.len() != r.len() {
            // Tried to read beyond seq end; this is an out-of-bounds variant.
            tracing::debug!("Bailing out on out-of-bounds variant: {}", &var);
            return Ok(var);
        }
        let seq = if strand == -1 { revcomp(&seq) } else { seq };

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

/// The distance from the position `start` to `end`, as `BaseOffsetPosition.__sub__` of
/// biocommons/hgvs computes it for positions without offset.  There is no position 0, e.g.,
/// c.0, so the distance from c.-1 to c.1 is 1.
fn distance_without_zero(start: i32, end: i32) -> i32 {
    end - start - i32::from(end > 0 && start < 0)
}

/// Whether `e` reports missing data rather than a problem with the variant.  The normalizer
/// reports a failed lookup of the reference as `ReplaceReferenceFailed`.
fn is_data_error(e: &Error) -> bool {
    matches!(
        e,
        Error::DataError(_)
            | Error::NormalizationFailed(
                normalizer::Error::DataError(_) | normalizer::Error::ReplaceReferenceFailed(_)
            )
    )
}

/// Cache key of `build_alignment_mapper_cached`.
#[derive(Hash, PartialEq, Eq)]
struct AlignmentMapperKey {
    data_version: String,
    schema_version: String,
    strict_bounds: bool,
    tx_ac: String,
    alt_ac: String,
    alt_aln_method: String,
}

/// Borrowed form of `AlignmentMapperKey`, so that lookups need no allocation.
///
/// It must hash like `AlignmentMapperKey`, hence the same fields in the same order.
#[derive(Hash)]
struct AlignmentMapperKeyRef<'a> {
    data_version: &'a str,
    schema_version: &'a str,
    strict_bounds: bool,
    tx_ac: &'a str,
    alt_ac: &'a str,
    alt_aln_method: &'a str,
}

impl Equivalent<AlignmentMapperKey> for AlignmentMapperKeyRef<'_> {
    fn equivalent(&self, key: &AlignmentMapperKey) -> bool {
        self.data_version == key.data_version
            && self.schema_version == key.schema_version
            && self.strict_bounds == key.strict_bounds
            && self.tx_ac == key.tx_ac
            && self.alt_ac == key.alt_ac
            && self.alt_aln_method == key.alt_aln_method
    }
}

impl AlignmentMapperKeyRef<'_> {
    fn to_owned(&self) -> AlignmentMapperKey {
        AlignmentMapperKey {
            data_version: self.data_version.to_string(),
            schema_version: self.schema_version.to_string(),
            strict_bounds: self.strict_bounds,
            tx_ac: self.tx_ac.to_string(),
            alt_ac: self.alt_ac.to_string(),
            alt_aln_method: self.alt_aln_method.to_string(),
        }
    }
}

thread_local! {
    /// Each thread has its own cache.  Even hits update the eviction state, so a shared cache
    /// needs a lock for them, and threads would also contend on the reference counts of the
    /// shared values.
    static ALIGNMENT_MAPPER_CACHE: RefCell<Cache<AlignmentMapperKey, Result<Arc<alignment::Mapper>, Error>>> =
        RefCell::new(Cache::new(1000));
}

/// A cached version of `alignment::Mapper::new`.
/// All of this function's arguments contribute to the key;
/// that is why the supplied provider's `data_version` and `schema_version`
/// should return sensible values which allow distinguishing them from other providers.
///
/// Because the cached value must implement `Clone`
/// and the type is `Result<Arc<alignment::Mapper>, Error>`,
/// `Error`, too, must implement `Clone`.
/// Sadly, postgres Errors do not do that, so either:
/// 1. wrap non-clonable errors in `Arc`
/// 2. convert the result to an option instead
/// 3. return a generic error instead
fn build_alignment_mapper_cached(
    provider: &Arc<dyn Provider + Send + Sync>,
    strict_bounds: bool,
    tx_ac: &str,
    alt_ac: &str,
    alt_aln_method: &str,
) -> Result<Arc<alignment::Mapper>, Error> {
    let key = AlignmentMapperKeyRef {
        data_version: provider.data_version(),
        schema_version: provider.schema_version(),
        strict_bounds,
        tx_ac,
        alt_ac,
        alt_aln_method,
    };
    if let Some(hit) = ALIGNMENT_MAPPER_CACHE.with_borrow(|c| c.get(&key).cloned()) {
        return hit;
    }
    let key = key.to_owned();
    let result = alignment::Mapper::new(
        &alignment::Config { strict_bounds },
        provider.clone(),
        tx_ac,
        alt_ac,
        alt_aln_method,
    )
    .map(Arc::new);
    ALIGNMENT_MAPPER_CACHE.with_borrow_mut(|c| c.insert(key, result.clone()));
    result
}
#[cfg(test)]
mod test {
    use std::{
        path::{Path, PathBuf},
        str::FromStr,
    };

    use anyhow::Error;
    use pretty_assertions::assert_eq;
    use regex::Regex;
    use test_log::test;

    use crate::{
        data::uta_sr::test_helpers::build_provider,
        parser::{HgvsVariant, NoRef},
    };

    use super::{Config, Mapper};

    /// Cache lookups only hit if the borrowed key hashes like the owned one.
    #[test]
    fn alignment_mapper_key_ref_matches_key() {
        use quick_cache::Equivalent;
        use std::hash::{BuildHasher, RandomState};

        let key_ref = super::AlignmentMapperKeyRef {
            data_version: "uta_20210129",
            schema_version: "1.1",
            strict_bounds: true,
            tx_ac: "NM_001637.3",
            alt_ac: "NC_000007.13",
            alt_aln_method: "splign",
        };
        let key = key_ref.to_owned();
        let state = RandomState::new();

        assert_eq!(state.hash_one(&key_ref), state.hash_one(&key));
        assert!(key_ref.equivalent(&key));
    }

    #[test]
    fn issue_131() -> Result<(), Error> {
        let mapper = build_mapper()?;

        let var_c = HgvsVariant::from_str("NM_001253909.2:c.416_417insGTG")?;
        let var_p_test = mapper.c_to_p(&var_c, None)?;

        assert_eq!(format!("{}", &var_p_test), "NP_001240838.1:p.=");
        insta::assert_yaml_snapshot!(&var_p_test);

        Ok(())
    }

    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::Mapper>();
    }

    fn build_mapper() -> Result<Mapper, Error> {
        let provider = build_provider()?;
        let config = Config {
            renormalize_g: false,
            ..Default::default()
        };
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

    /// c.1257 is the last base of the stop codon, so the copy goes into the 3' UTR.
    /// biocommons/hgvs gives `p.=`.
    #[test]
    fn map_of_dup_at_cds_end() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_001051.2:c.1257dupG"; // gene SSTR3
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_001042.1:p.?");

        Ok(())
    }

    #[test]
    fn map_of_ins_three_prime_utr() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_004985.4:c.567_*1insCCC"; // gene KRAS
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_004976.2:p.?");

        Ok(())
    }

    #[test]
    fn map_of_dup_three_prime_utr() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_153223.3:c.2959_*1dup"; // gene CEP120
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_694955.2:p.?");

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
        use std::{
            path::{Path, PathBuf},
            sync::Arc,
        };

        use anyhow::Error;

        use crate::data::interface;
        use crate::{
            data,
            data::interface::TxIdentityInfo,
            mapper::variant::{Config, Mapper},
        };
        use indexmap::IndexMap;
        use std::sync::atomic::AtomicUsize;

        static PROVIDER_COUNT: AtomicUsize = AtomicUsize::new(0);

        #[derive(Debug, serde::Deserialize)]
        pub struct ProviderRecord {
            pub accession: String,
            pub transcript_sequence: String,
            pub cds_start_i: i32,
            pub cds_end_i: i32,
            #[serde(skip)]
            pub translation_exceptions: Vec<interface::TranslationException>,
            #[serde(skip)]
            pub lengths: Vec<i32>,
        }

        pub struct Provider {
            data_version: String,
            schema_version: String,
            records: Vec<ProviderRecord>,
            translation_table: crate::sequences::TranslationTable,
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
                Ok(Self::with_records(records))
            }

            pub fn with_records(records: Vec<ProviderRecord>) -> Self {
                let number = PROVIDER_COUNT.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                let dummy_version = format!("provider_{number}");
                Self {
                    records,
                    data_version: dummy_version.clone(),
                    schema_version: dummy_version,
                    translation_table: Default::default(),
                }
            }

            /// Use `translation_table` for every transcript.
            pub fn with_translation_table(
                self,
                translation_table: crate::sequences::TranslationTable,
            ) -> Self {
                Self {
                    translation_table,
                    ..self
                }
            }
        }

        impl interface::Provider for Provider {
            fn data_version(&self) -> &str {
                &self.data_version
            }

            fn schema_version(&self) -> &str {
                &self.schema_version
            }

            fn get_assembly_map(
                &self,
                _assembly: &str,
            ) -> Result<IndexMap<String, String>, data::error::Error> {
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
                        // As seqrepo, return the part of the sequence within its bounds.
                        let seq = &record.transcript_sequence;
                        let end = end.unwrap_or(seq.len()).min(seq.len());
                        let begin = begin.unwrap_or(0).min(end);
                        return Ok(seq[begin..end].to_string());
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
                            cds_start_i: Some(record.cds_start_i),
                            cds_end_i: Some(record.cds_end_i),
                            lengths: record.lengths.clone(),
                            hgnc: "MOCK".to_string(),
                            translation_table: self.translation_table,
                        });
                    }
                }
                Err(crate::data::error::Error::NoSequenceRecord(
                    tx_ac.to_string(),
                ))
            }

            fn get_tx_translation_exceptions(
                &self,
                tx_ac: &str,
            ) -> Result<Vec<interface::TranslationException>, crate::data::error::Error>
            {
                self.records
                    .iter()
                    .find(|record| record.accession == tx_ac)
                    .map(|record| record.translation_exceptions.clone())
                    .ok_or_else(|| crate::data::error::Error::NoSequenceRecord(tx_ac.to_string()))
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
            let provider = Arc::new(Provider::new(&path)?);
            let config = Config {
                strict_bounds,
                renormalize_g: false,
                ..Default::default()
            };
            Ok(Mapper::new(&config, provider))
        }
    }

    /// Transcripts whose proteins need translation exceptions:
    ///
    /// * `NM_999990.1`: `MAKLWEPURV*`, with selenocysteine at the UGA codon 8.
    /// * `NM_999991.1`: `MAKLWEPRV*`, with methionine at the ACG start codon.
    /// * `NM_999992.1`: `MAKXLWERV*`, with an unspecified amino acid at the UAG codon 4, i.e.,
    ///   a readthrough of this stop codon.
    fn translation_exceptions_provider() -> std::sync::Arc<sanity_mock::Provider> {
        use crate::data::interface::TranslationException;

        std::sync::Arc::new(sanity_mock::Provider::with_records(vec![
            sanity_mock::ProviderRecord {
                accession: "NM_999990.1".to_string(),
                transcript_sequence: "AAAATCAAAATGGCGAAACTGTGGGAACCGTGACGCGTGTAAGGGG".to_string(),
                cds_start_i: 9,
                cds_end_i: 42,
                translation_exceptions: vec![TranslationException {
                    position: 8,
                    amino_acid: 'U',
                }],
                lengths: Vec::new(),
            },
            sanity_mock::ProviderRecord {
                accession: "NM_999991.1".to_string(),
                transcript_sequence: "AAAATCAAAACGGCGAAACTGTGGGAACCGCGCGTGTAAGGGG".to_string(),
                cds_start_i: 9,
                cds_end_i: 39,
                translation_exceptions: vec![TranslationException {
                    position: 1,
                    amino_acid: 'M',
                }],
                lengths: Vec::new(),
            },
            sanity_mock::ProviderRecord {
                accession: "NM_999992.1".to_string(),
                transcript_sequence: "AAAATCAAAATGGCGAAATAGCTGTGGGAACGCGTGTAAGGGG".to_string(),
                cds_start_i: 9,
                cds_end_i: 39,
                translation_exceptions: vec![TranslationException {
                    position: 4,
                    amino_acid: 'X',
                }],
                lengths: Vec::new(),
            },
        ]))
    }

    #[test]
    fn translation_exceptions_in_reference_protein() -> Result<(), Error> {
        for (tx_ac, aa_sequence) in [
            ("NM_999990.1", "MAKLWEPURV*"),
            ("NM_999991.1", "MAKLWEPRV*"),
            ("NM_999992.1", "MAKXLWERV*"),
        ] {
            let ref_data = crate::mapper::altseq::RefTranscriptData::new(
                translation_exceptions_provider(),
                tx_ac,
                Some("MOCK"),
            )?;
            assert_eq!(ref_data.aa_sequence.as_ref(), aa_sequence, "{tx_ac}");
        }

        Ok(())
    }

    /// A translation exception applies only at the codons of the reference that the variant
    /// leaves unchanged and in frame.
    #[test]
    fn hgvs_c_to_p_translation_exceptions() -> Result<(), Error> {
        let config = Config {
            strict_bounds: false,
            renormalize_g: false,
            ..Default::default()
        };
        let mapper = Mapper::new(&config, translation_exceptions_provider());

        for (hgvsc, hgvsp_expected) in [
            // A new UGA is a stop codon.
            ("NM_999990.1:c.15G>A", "MOCK:p.Trp5Ter"),
            // The selenocysteine codon changes.
            ("NM_999990.1:c.24A>G", "MOCK:p.Sec8Trp"),
            // The selenocysteine codon stays before a variant ...
            ("NM_999990.1:c.26G>A", "MOCK:p.Arg9His"),
            // ... and after an in-frame deletion.
            ("NM_999990.1:c.16_18del", "MOCK:p.Glu6del"),
            // The ACG start codon stays methionine.
            ("NM_999991.1:c.16_18del", "MOCK:p.Glu6del"),
            // The readthrough stop codon does not end the protein, before a variant ...
            ("NM_999992.1:c.23G>A", "MOCK:p.Arg8His"),
            // ... and after an in-frame deletion.
            ("NM_999992.1:c.4_6del", "MOCK:p.Ala2del"),
        ] {
            let var_p = mapper.c_to_p(&HgvsVariant::from_str(hgvsc)?, Some("MOCK"))?;
            assert_eq!(format!("{}", &var_p), hgvsp_expected, "{hgvsc}");
        }

        Ok(())
    }

    /// Transcripts whose CDS has no stop codon, e.g., GENCODE `cds_end_NF`:
    ///
    /// * `NM_999993.1`: `ATG GCC AAG CTG TGG GAA` (`MAKLWE`).
    /// * `NM_999994.1`: `ATG GCC AAG CTG TGG GAA GC` (`MAKLWE` and an incomplete codon).
    /// * `NM_999996.1`: `NM_999993.1` followed by `GGCTAACTGATAGG`, with stop codons after
    ///   c.18 in two frames.
    ///
    /// The CDS of `NM_999995.1` has a stop codon, but its length is not a multiple of 3, as a
    /// +1 ribosomal frameshift skips the `T` after `TCC`: `ATG GCC TCC T GAT GGC TAA`
    /// (`MASDG*`).
    fn no_stop_codon_provider() -> std::sync::Arc<sanity_mock::Provider> {
        let record =
            |accession: &str, transcript_sequence: &str, cds_end_i| sanity_mock::ProviderRecord {
                accession: accession.to_string(),
                transcript_sequence: transcript_sequence.to_string(),
                cds_start_i: 9,
                cds_end_i,
                translation_exceptions: Vec::new(),
                lengths: Vec::new(),
            };
        std::sync::Arc::new(sanity_mock::Provider::with_records(vec![
            record("NM_999993.1", "AAAATCAAAATGGCCAAGCTGTGGGAA", 27),
            record("NM_999994.1", "AAAATCAAAATGGCCAAGCTGTGGGAAGC", 29),
            record("NM_999995.1", "AAAATCAAAATGGCCTCCTGATGGCTAAGGGG", 28),
            record(
                "NM_999996.1",
                "AAAATCAAAATGGCCAAGCTGTGGGAAGGCTAACTGATAGG",
                27,
            ),
        ]))
    }

    /// The reference protein holds the complete codons of the CDS, and `cds_end` tells how the
    /// CDS ends.  A CDS with an incomplete codon is rejected if a stop codon in frame shows
    /// that the CDS changes the frame.
    #[test]
    fn reference_protein_and_cds_end() -> Result<(), Error> {
        use crate::mapper::altseq::{CdsEnd, RefTranscriptData};

        for (provider, tx_ac, aa_sequence, cds_stop, cds_end) in [
            (
                translation_exceptions_provider(),
                "NM_999991.1",
                "MAKLWEPRV*",
                39,
                CdsEnd::StopCodon,
            ),
            (
                no_stop_codon_provider(),
                "NM_999993.1",
                "MAKLWE",
                27,
                CdsEnd::NoStopCodon,
            ),
            (
                no_stop_codon_provider(),
                "NM_999994.1",
                "MAKLWE",
                29,
                CdsEnd::IncompleteCodon("GC".to_string()),
            ),
        ] {
            let ref_data = RefTranscriptData::new(provider, tx_ac, Some("MOCK"))?;
            assert_eq!(
                (
                    ref_data.aa_sequence.as_ref(),
                    ref_data.cds_stop,
                    ref_data.cds_end
                ),
                (aa_sequence, cds_stop, cds_end),
                "{tx_ac}"
            );
        }

        assert!(matches!(
            crate::mapper::altseq::RefTranscriptData::new(
                no_stop_codon_provider(),
                "NM_999995.1",
                Some("MOCK"),
            ),
            Err(crate::mapper::error::Error::TranscriptLengthInvalid(tx_ac, 19))
                if tx_ac == "NM_999995.1"
        ));

        Ok(())
    }

    /// Without a stop codon, the last amino acid is an ordinary one, and the amino acids after
    /// it are unknown.
    #[test]
    fn hgvs_c_to_p_no_stop_codon() -> Result<(), Error> {
        let config = Config {
            strict_bounds: false,
            renormalize_g: false,
            ..Default::default()
        };
        let mapper = Mapper::new(&config, no_stop_codon_provider());

        for (hgvsc, hgvsp_expected) in [
            ("NM_999993.1:c.18A>G", "MOCK:p.Glu6="),
            ("NM_999993.1:c.18A>C", "MOCK:p.Glu6Asp"),
            ("NM_999993.1:c.16_18del", "MOCK:p.Glu6del"),
            ("NM_999993.1:c.13_15del", "MOCK:p.Trp5del"),
            ("NM_999993.1:c.15_16insGGG", "MOCK:p.Trp5_Glu6insGly"),
            ("NM_999993.1:c.15G>A", "MOCK:p.Trp5Ter"),
            ("NM_999993.1:c.13_14insCCTA", "MOCK:p.Trp5SerfsTer2"),
            ("NM_999993.1:c.17_18insC", "MOCK:p.Glu6AspfsTer?"),
            // A new stop codon gives the same result as with a stop codon at the end.
            ("NM_999993.1:c.15_16insGGGTAA", "MOCK:p.Trp5_Glu6insGlyTer"),
            ("NM_999993.1:c.4_6delinsGGGTAA", "MOCK:p.Ala2delinsGlyTer"),
            // The protein ends at c.18 also if a stop codon follows.
            ("NM_999996.1:c.18A>G", "MOCK:p.Glu6="),
            ("NM_999996.1:c.18A>C", "MOCK:p.Glu6Asp"),
            ("NM_999996.1:c.4_6dup", "MOCK:p.Ala2dup"),
            ("NM_999996.1:c.15_16insGGGTAA", "MOCK:p.Trp5_Glu6insGlyTer"),
            ("NM_999996.1:c.17_18insCC", "MOCK:p.Glu6AspfsTer?"),
            // The alt protein keeps the last amino acid after an in-frame dup.
            ("NM_999993.1:c.4_6dup", "MOCK:p.Ala2dup"),
            ("NM_999993.1:c.6_7insGCC", "MOCK:p.Ala2dup"),
            // The CDS goes on after c.18, so there is no 3' UTR, and the amino acids after
            // c.18 are unknown.  This also holds for a variant that starts in the CDS.
            ("NM_999996.1:c.18_*1insA", "MOCK:p.?"),
            ("NM_999996.1:c.*1G>A", "MOCK:p.?"),
            ("NM_999996.1:c.*2dup", "MOCK:p.?"),
            ("NM_999996.1:c.*1_*3inv", "MOCK:p.?"),
            ("NM_999996.1:c.17_*1del", "MOCK:p.?"),
            ("NM_999996.1:c.18_*1delinsTC", "MOCK:p.?"),
            // The CDS ends before the new frame changes a complete codon.
            ("NM_999993.1:c.18del", "MOCK:p.?"),
            ("NM_999993.1:c.16del", "MOCK:p.?"),
            // The copy follows the last amino acid.
            ("NM_999993.1:c.16_18dup", "MOCK:p.?"),
            ("NM_999994.1:c.18A>G", "MOCK:p.Glu6="),
            ("NM_999994.1:c.18A>C", "MOCK:p.Glu6Asp"),
            ("NM_999994.1:c.16_18del", "MOCK:p.Glu6del"),
            ("NM_999994.1:c.4_6dup", "MOCK:p.Ala2dup"),
            ("NM_999994.1:c.15G>A", "MOCK:p.Trp5Ter"),
            ("NM_999994.1:c.13_14insCCTA", "MOCK:p.Trp5SerfsTer2"),
            // The new frame completes the incomplete codon ...
            ("NM_999994.1:c.16del", "MOCK:p.Glu6LysfsTer?"),
            // ... but can keep the last amino acid.
            ("NM_999994.1:c.18del", "MOCK:p.?"),
            // Changes after the last complete codon.
            ("NM_999994.1:c.19G>A", "MOCK:p.?"),
            ("NM_999994.1:c.19_20insGTGGGAAC", "MOCK:p.?"),
            ("NM_999994.1:c.18_19insTAA", "MOCK:p.?"),
            ("NM_999994.1:c.20dup", "MOCK:p.?"),
            ("NM_999994.1:c.19_20dup", "MOCK:p.?"),
            // c.-1 is an A, so the deletion keeps the ATG.
            ("NM_999994.1:c.1del", "MOCK:p.?"),
        ] {
            let var_p = mapper.c_to_p(&HgvsVariant::from_str(hgvsc)?, Some("MOCK"))?;
            assert_eq!(format!("{}", &var_p), hgvsp_expected, "{hgvsc}");
        }

        Ok(())
    }

    /// Transcripts whose CDS ends with an incomplete stop codon, as many mitochondrial ones do:
    ///
    /// * `NM_999982.1`: `ATG GCC AAG CTG TGG T`, followed by `CCCC`.
    /// * `NM_999981.1`: `ATG GCC AAG CTG TGG TA`, followed by `CCCC`.
    /// * `NM_999980.1`: `ATG GCC AAG CTG TGG C`, followed by `CCCC`.
    fn poly_a_provider(
        translation_table: crate::sequences::TranslationTable,
    ) -> std::sync::Arc<sanity_mock::Provider> {
        let record =
            |accession: &str, transcript_sequence: &str, cds_end_i| sanity_mock::ProviderRecord {
                accession: accession.to_string(),
                transcript_sequence: transcript_sequence.to_string(),
                cds_start_i: 5,
                cds_end_i,
                translation_exceptions: Vec::new(),
                lengths: Vec::new(),
            };
        std::sync::Arc::new(
            sanity_mock::Provider::with_records(vec![
                record("NM_999982.1", "AAAACATGGCCAAGCTGTGGTCCCC", 21),
                record("NM_999981.1", "AAAACATGGCCAAGCTGTGGTACCCC", 22),
                record("NM_999980.1", "AAAACATGGCCAAGCTGTGGCCCCC", 21),
            ])
            .with_translation_table(translation_table),
        )
    }

    /// As in biocommons/hgvs, the poly(A) tail completes a stop codon that ends with `T` or
    /// `TA`, but only with the vertebrate mitochondrial translation table.  The completing
    /// `A`s go into the translated bases only, so the transcript sequence and `cds_stop` stay.
    #[test]
    fn poly_a_completes_mitochondrial_stop_codon() -> Result<(), Error> {
        use crate::mapper::altseq::{CdsEnd, RefTranscriptData};
        use crate::sequences::TranslationTable;

        for (translation_table, tx_ac, transcript_sequence, aa_sequence, cds_stop, cds_end) in [
            (
                TranslationTable::VertebrateMitochondrial,
                "NM_999982.1",
                "AAAACATGGCCAAGCTGTGGTCCCC",
                "MAKLW*",
                21,
                CdsEnd::StopCodon,
            ),
            (
                TranslationTable::VertebrateMitochondrial,
                "NM_999981.1",
                "AAAACATGGCCAAGCTGTGGTACCCC",
                "MAKLW*",
                22,
                CdsEnd::StopCodon,
            ),
            (
                TranslationTable::VertebrateMitochondrial,
                "NM_999980.1",
                "AAAACATGGCCAAGCTGTGGCCCCC",
                "MAKLW",
                21,
                CdsEnd::IncompleteCodon("C".to_string()),
            ),
            (
                TranslationTable::Standard,
                "NM_999982.1",
                "AAAACATGGCCAAGCTGTGGTCCCC",
                "MAKLW",
                21,
                CdsEnd::IncompleteCodon("T".to_string()),
            ),
        ] {
            let ref_data =
                RefTranscriptData::new(poly_a_provider(translation_table), tx_ac, Some("MOCK"))?;
            assert_eq!(
                (
                    ref_data.transcript_sequence.as_ref(),
                    ref_data.aa_sequence.as_ref(),
                    ref_data.cds_stop,
                    ref_data.cds_end,
                ),
                (transcript_sequence, aa_sequence, cds_stop, cds_end),
                "{tx_ac} with {translation_table:?}"
            );
        }

        Ok(())
    }

    /// `c_to_p` on the mitochondrial transcripts of `poly_a_provider()`.  The expectations are
    /// those of biocommons/hgvs, without the parentheses, except where marked.
    #[test]
    fn hgvs_c_to_p_mitochondrial_poly_a() -> Result<(), Error> {
        let config = Config {
            strict_bounds: false,
            renormalize_g: false,
            ..Default::default()
        };
        let mapper = Mapper::new(
            &config,
            poly_a_provider(crate::sequences::TranslationTable::VertebrateMitochondrial),
        );

        for (hgvsc, hgvsp_expected) in [
            ("NM_999982.1:c.4G>A", "MOCK:p.Ala2Thr"),
            // `TGA` codes Trp in the vertebrate mitochondrial table.
            ("NM_999982.1:c.15G>A", "MOCK:p.Trp5="),
            ("NM_999982.1:c.16T>C", "MOCK:p.Ter6GlnextTer?"),
            ("NM_999982.1:c.13_15del", "MOCK:p.Trp5Ter"),
            ("NM_999982.1:c.10_12dup", "MOCK:p.Leu4dup"),
            ("NM_999982.1:c.*1C>G", "MOCK:p.?"),
            // The c.* positions count the bases of the transcript, not the completing `A`s, so
            // the deletions end with the `C`s after c.16.  biocommons/hgvs deletes the `A`s
            // instead, see `DIFFERENCES.md`.
            ("NM_999982.1:c.16_*1del", "MOCK:p.Ter6ProextTer?"),
            ("NM_999982.1:c.16_*2del", "MOCK:p.Ter6XaaextTer?"),
            // Bases inserted after the CDS go in front of the completing `A`s.  The CDS has 16
            // bases, so `c.16_17insG` is not canonical.  It is used on purpose, as
            // biocommons/hgvs counts it as exonic.  The canonical `c.16_*1insG` is a 3' UTR
            // variant by its positions.
            ("NM_999982.1:c.16dup", "MOCK:p.Ter6LeuextTer?"),
            ("NM_999982.1:c.16_17insG", "MOCK:p.Ter6TrpextTer?"),
            ("NM_999982.1:c.16_*1insG", "MOCK:p.?"),
            ("NM_999981.1:c.4G>A", "MOCK:p.Ala2Thr"),
            ("NM_999981.1:c.17A>C", "MOCK:p.Ter6SerextTer?"),
            ("NM_999981.1:c.16_17del", "MOCK:p.Ter6ThrextTer?"),
        ] {
            let var_p = mapper.c_to_p(&HgvsVariant::from_str(hgvsc)?, Some("MOCK"))?;
            assert_eq!(format!("{}", &var_p), hgvsp_expected, "{hgvsc}");
        }

        // The alt transcript has the bases of the transcript only.
        let ref_data = crate::mapper::altseq::RefTranscriptData::new(
            poly_a_provider(crate::sequences::TranslationTable::VertebrateMitochondrial),
            "NM_999982.1",
            Some("MOCK"),
        )?;
        for (hgvsc, transcript_sequence) in [
            ("NM_999982.1:c.16_*1del", "AAAACATGGCCAAGCTGTGGCCC"),
            ("NM_999982.1:c.16T>C", "AAAACATGGCCAAGCTGTGGCCCCC"),
        ] {
            let alt_data =
                crate::mapper::altseq::AltSeqBuilder::new(HgvsVariant::from_str(hgvsc)?, &ref_data)
                    .build_altseq()?;
            assert_eq!(
                alt_data[0].transcript_sequence, transcript_sequence,
                "{hgvsc}"
            );
        }

        Ok(())
    }

    fn test_hgvs_c_to_p_conversion(hgvsc: &str, hgvsp_expected: &str) -> Result<(), Error> {
        let mapper = sanity_mock::build_mapper(false)?;

        let var_c = HgvsVariant::from_str(hgvsc)?;
        let ac_p = "MOCK";

        let var_p = mapper.c_to_p(&var_c, Some(ac_p))?;
        let hgvsp_actual = format!("{}", var_p);

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
    fn hgvs_c_to_p_sub_three_prime_ut() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.*3G>A";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_ins_three_prime_utr() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.39_*1insA";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_three_prime_utr() -> Result<(), Error> {
        let hgvsc = "NM_999999.1:c.12_*1dup";
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

    /// As in biocommons/hgvs, a reference protein with more than one stop codon gives `p.?`
    /// for every variant.  `NM_999992.1` codes `MG*ARQPAL*` after a 5' UTR of 9 bases, and
    /// `NM_999986.1` codes `MGKP*A*` without a 5' UTR.
    #[test]
    fn hgvs_c_to_p_multiple_stop_codons() -> Result<(), Error> {
        let sanity = sanity_mock::build_mapper(false)?;
        let config = Config {
            strict_bounds: false,
            renormalize_g: false,
            ..Default::default()
        };
        let no_utr = Mapper::new(
            &config,
            std::sync::Arc::new(sanity_mock::Provider::with_records(vec![
                sanity_mock::ProviderRecord {
                    accession: "NM_999986.1".to_string(),
                    transcript_sequence: "ATGGGGAAACCCTAGGCCTAGGG".to_string(),
                    cds_start_i: 0,
                    cds_end_i: 21,
                    translation_exceptions: Vec::new(),
                    lengths: Vec::new(),
                },
            ])),
        );

        for (mapper, hgvsc) in [
            (&sanity, "NM_999992.1:c.1A>G"),
            (&sanity, "NM_999992.1:c.4G>A"),
            (&sanity, "NM_999992.1:c.4_6dup"),
            (&sanity, "NM_999992.1:c.6_7insGGG"),
            (&sanity, "NM_999992.1:c.4_5inv"),
            (&sanity, "NM_999992.1:c.13C>T"),
            (&no_utr, "NM_999986.1:c.1A>G"),
            (&no_utr, "NM_999986.1:c.4G>A"),
            (&no_utr, "NM_999986.1:c.4_6dup"),
            (&no_utr, "NM_999986.1:c.6_7insGGG"),
            (&no_utr, "NM_999986.1:c.4_5inv"),
            (&no_utr, "NM_999986.1:c.10C>T"),
        ] {
            let var_p = mapper.c_to_p(&HgvsVariant::from_str(hgvsc)?, Some("MOCK"))?;
            assert_eq!(format!("{}", &var_p), "MOCK:p.?", "{hgvsc}");
        }

        Ok(())
    }

    /// As in biocommons/hgvs (biocommons/hgvs#830), a stop codon that an insertion or dup forms
    /// within its own bases ends the translation there, so the variant is in frame.  The CDS of
    /// `NM_999985.1` is `ATG AAA GCG AAA GCG GCT AAG GCT GAA TAG` (`MKAKAAKAE*`), and the one
    /// of `NM_999984.1` is `ATG AAA GCG AAA GCG TTT CGC GTA GAA TAG` (`MKAKAFRVE*`).  The
    /// expectations are those of biocommons/hgvs, without the parentheses.
    #[test]
    fn hgvs_c_to_p_stop_codon_in_insertion() -> Result<(), Error> {
        let record = |accession: &str, cds: &str| sanity_mock::ProviderRecord {
            accession: accession.to_string(),
            transcript_sequence: format!("AAAATCAAAC{cds}GGGGGCCCCC"),
            cds_start_i: 10,
            cds_end_i: 40,
            translation_exceptions: Vec::new(),
            lengths: Vec::new(),
        };
        let config = Config {
            strict_bounds: false,
            renormalize_g: false,
            ..Default::default()
        };
        let mapper = Mapper::new(
            &config,
            std::sync::Arc::new(sanity_mock::Provider::with_records(vec![
                record("NM_999985.1", "ATGAAAGCGAAAGCGGCTAAGGCTGAATAG"),
                record("NM_999984.1", "ATGAAAGCGAAAGCGTTTCGCGTAGAATAG"),
            ])),
        );

        for (hgvs_c, expected) in [
            // The copy holds the stop codons `TAA` and `TGA` in the new frame.
            ("NM_999985.1:c.17_29dup", "MOCK:p.Ter10TyrextTer1"),
            (
                "NM_999985.1:c.29_30insCTAAGGCTGAATA",
                "MOCK:p.Ter10TyrextTer1",
            ),
            ("NM_999985.1:c.6_7insTAAGG", "MOCK:p.Ala3Ter"),
            // The stop codon is not in the frame of the CDS.
            ("NM_999985.1:c.6_7insGGTAAGG", "MOCK:p.Ala3GlyfsTer7"),
            // The copy forms the stop codon `TAG` in front of `AA TAG`.
            ("NM_999984.1:c.21_25dup", "MOCK:p.Glu9Ala"),
            ("NM_999984.1:c.22_26dup", "MOCK:p.Glu9="),
        ] {
            let var_p = mapper.c_to_p(&HgvsVariant::from_str(hgvs_c)?, Some("MOCK"))?;
            assert_eq!(format!("{}", &var_p), expected, "{hgvs_c}");
        }

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

    // cf. https://github.com/biocommons/hgvs/pull/765
    #[test]
    fn test_765_vm_g_to_c() -> Result<(), Error> {
        let mapper = build_mapper()?;
        let var_g = HgvsVariant::from_str("NC_000007.13:g.106545832_106545833insAG")?;
        let var_c = mapper.g_to_c(&var_g, "NM_002649.2", "splign")?;
        assert_eq!(var_c.to_string(), "NM_002649.2:c.3309_*1insAG");
        Ok(())
    }

    /// There is no c.0, so an insertion between c.-1 and c.1 stays an insertion.
    #[test]
    fn g_to_c_ins_before_c1() -> Result<(), Error> {
        let mapper = build_mapper()?;
        for (hgvs_g, tx_ac, hgvs_c) in [
            // PTEN, plus strand.
            (
                "NC_000010.10:g.89624226_89624227insA",
                "NM_000314.4",
                "NM_000314.4:c.-1_1insA",
            ),
            // KRAS, minus strand.
            (
                "NC_000012.11:g.25398318_25398319insT",
                "NM_004985.4",
                "NM_004985.4:c.-1_1insA",
            ),
        ] {
            let var_g = HgvsVariant::from_str(hgvs_g)?;
            let var_c = mapper.g_to_c(&var_g, tx_ac, "splign")?;
            assert_eq!(var_c.to_string(), hgvs_c, "{hgvs_g}");
            assert_eq!(
                mapper.c_to_g(&var_c, var_g.accession(), "splign")?,
                var_g,
                "{hgvs_c}"
            );
        }
        Ok(())
    }

    // The following tests correspond to `tests/issues/test_655.py` of biocommons/hgvs.

    /// Check `c_to_p()` of `hgvs_c`, where `None` expects a dup without reference.
    fn check_c_to_p(
        mapper: &Mapper,
        hgvs_c: &str,
        prot_ac: Option<&str>,
        expected: Option<&str>,
    ) -> Result<(), Error> {
        let var_p = mapper.c_to_p(&HgvsVariant::from_str(hgvs_c)?, prot_ac);
        match expected {
            Some(expected) => assert_eq!(format!("{}", &var_p?), expected, "{hgvs_c}"),
            None => assert!(
                matches!(
                    var_p,
                    Err(crate::mapper::Error::DuplicationReferenceMissing(_))
                ),
                "{hgvs_c} gave {var_p:?}"
            ),
        }
        Ok(())
    }

    /// `NM_999999.1` of `sanity_cp.tsv` with exons of 10, 25, and 30 bases, so exon 2 spans
    /// c.2 to c.26.
    fn build_mapper_655_sanity(ins_at_boundary_is_intronic: bool) -> Mapper {
        let provider = std::sync::Arc::new(sanity_mock::Provider::with_records(vec![
            sanity_mock::ProviderRecord {
                accession: "NM_999999.1".to_string(),
                transcript_sequence: "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG".to_string(),
                cds_start_i: 9,
                cds_end_i: 39,
                translation_exceptions: Vec::new(),
                lengths: vec![10, 25, 30],
            },
        ]));
        let config = Config {
            prevalidation_level: crate::validator::ValidationLevel::Intrinsic,
            strict_bounds: false,
            renormalize_g: false,
            ins_at_boundary_is_intronic,
            ..Default::default()
        };
        Mapper::new(&config, provider)
    }

    /// With `false`, the dups of intron bases fail: there is no genome to take their reference
    /// from.
    #[test]
    fn issue_655_sanity_c_to_p() -> Result<(), Error> {
        let exonic = build_mapper_655_sanity(false);
        let intronic = build_mapper_655_sanity(true);
        for (hgvs_c, expected_exonic) in [
            ("NM_999999.1:c.2-1_2insG", Some("MOCK:p.Met1?")),
            ("NM_999999.1:c.26_26+1insC", Some("MOCK:p.Lys9AsnfsTer?")),
            ("NM_999999.1:c.2dup", Some("MOCK:p.Met1?")),
            ("NM_999999.1:c.2-1dup", None),
            ("NM_999999.1:c.26dup", Some("MOCK:p.Ter10IleextTer?")),
            ("NM_999999.1:c.26+1dup", None),
        ] {
            check_c_to_p(&exonic, hgvs_c, Some("MOCK"), expected_exonic)?;
            check_c_to_p(&intronic, hgvs_c, Some("MOCK"), Some("MOCK:p.?"))?;
        }
        Ok(())
    }

    /// KRAS, minus strand, exon 2 ends at c.111.  Without `replace_reference`, as the
    /// seqrepo cache lacks the transcript parts.
    #[test]
    fn issue_655_real_c_to_p() -> Result<(), Error> {
        let build = |ins_at_boundary_is_intronic| -> Result<Mapper, Error> {
            let config = Config {
                replace_reference: false,
                renormalize_g: false,
                ins_at_boundary_is_intronic,
                ..Default::default()
            };
            Ok(Mapper::new(&config, build_provider()?))
        };
        let exonic = build(false)?;
        let intronic = build(true)?;
        for (hgvs_c, expected_exonic) in [
            (
                "NM_004985.4:c.111dup",
                Some("NP_004976.2:p.Asp38GlyfsTer10"),
            ),
            (
                "NM_004985.4:c.110_111dup",
                Some("NP_004976.2:p.Asp38ArgfsTer8"),
            ),
            ("NM_004985.4:c.111+1dup", None),
            ("NM_004985.4:c.111+1_111+4dup", None),
        ] {
            check_c_to_p(&exonic, hgvs_c, None, expected_exonic)?;
            check_c_to_p(&intronic, hgvs_c, None, Some("NP_004976.2:p.?"))?;
        }
        Ok(())
    }

    /// `AltSeqBuilder::variant_location()` gives the region that `c_to_p` uses, if both have the
    /// same `ins_at_boundary_is_intronic`.  `NM_999989.1` has exons of 10, 15, and 25 bases:
    /// exon 1 is the 5' UTR, exon 2 spans c.1 to c.15, and exon 3 spans c.16 to c.*10.  The CDS
    /// is `ATG AAA GCG AAA GCG TTT CGC GCG AAA TAG` (`MKAKAFRAK*`), and c.-2 to c.3 and c.30 to
    /// c.*2 are `AC ATG` and `TAG GG`.  The proteins are those of biocommons/hgvs, without the
    /// parentheses, except where marked.
    #[test]
    fn variant_location_agrees_with_c_to_p() -> Result<(), Error> {
        use crate::mapper::altseq::{ref_transcript_data_cached, AltSeqBuilder, VariantLocation};
        use VariantLocation::*;

        let provider = std::sync::Arc::new(sanity_mock::Provider::with_records(vec![
            sanity_mock::ProviderRecord {
                accession: "NM_999989.1".to_string(),
                transcript_sequence: "AAAATCAAACATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGGGGCCCCC"
                    .to_string(),
                cds_start_i: 10,
                cds_end_i: 40,
                translation_exceptions: Vec::new(),
                lengths: vec![10, 15, 25],
            },
        ]));
        let reference_data =
            ref_transcript_data_cached(provider.clone(), "NM_999989.1", Some("MOCK"))?;

        for ins_at_boundary_is_intronic in [false, true] {
            let config = Config {
                strict_bounds: false,
                renormalize_g: false,
                ins_at_boundary_is_intronic,
                ..Default::default()
            };
            let mapper = Mapper::new(&config, provider.clone());
            // With `ins_at_boundary_is_intronic`, the insertions and dups at an exon boundary
            // within the CDS are intronic.
            let (boundary, boundary_p) = if ins_at_boundary_is_intronic {
                (Intron, "MOCK:p.?")
            } else {
                (Cds, "MOCK:p.Arg7SerfsTer?")
            };
            for (hgvs_c, location, hgvs_p) in [
                ("c.4A>G", Cds, "MOCK:p.Lys2Glu"),
                ("c.-2_2del", Cds, "MOCK:p.Met1?"),
                ("c.2dup", Cds, "MOCK:p.Met1?"),
                ("c.29_30insC", Cds, "MOCK:p.Ter10TyrextTer?"),
                // These differ from biocommons/hgvs on purpose: the same change fits in front
                // of c.1 or after the stop codon and keeps both intact, see `DIFFERENCES.md`.
                // biocommons/hgvs gives `p.Met1?` for the first four, except `p.?` for `c.1dup`
                // with `ins_at_boundary_is_intronic`, as c.1 starts an exon.  It gives `p.(=)` for
                // the dups at the end and `p.(Ter10Ter)` for `c.30del`.
                ("c.-1_1insGCC", FivePrimeUtr, "MOCK:p.?"),
                ("c.1_2insA", FivePrimeUtr, "MOCK:p.?"),
                ("c.-1_1del", FivePrimeUtr, "MOCK:p.?"),
                ("c.1dup", FivePrimeUtr, "MOCK:p.?"),
                ("c.30dup", ThreePrimeUtr, "MOCK:p.?"),
                ("c.28_30dup", ThreePrimeUtr, "MOCK:p.?"),
                ("c.30del", ThreePrimeUtr, "MOCK:p.?"),
                ("c.-5C>G", FivePrimeUtr, "MOCK:p.?"),
                ("c.-1+5del", FivePrimeUtr, "MOCK:p.?"),
                ("c.*3G>T", ThreePrimeUtr, "MOCK:p.?"),
                ("c.30_*1insA", ThreePrimeUtr, "MOCK:p.?"),
                ("c.-5_*2dup", ThreePrimeUtr, "MOCK:p.?"),
                ("c.15+5del", Intron, "MOCK:p.?"),
                ("c.15_15+3del", Intron, "MOCK:p.?"),
                ("c.16-1_16insT", boundary, boundary_p),
                ("c.16dup", boundary, boundary_p),
                ("c.-10_*10del", WholeGene, "MOCK:p.0?"),
                ("c.-5_*5inv", WholeGene, "MOCK:p.?"),
            ] {
                let var_c = HgvsVariant::from_str(&format!("NM_999989.1:{hgvs_c}"))?;
                let builder = AltSeqBuilder {
                    ins_at_boundary_is_intronic,
                    ..AltSeqBuilder::new(var_c.clone(), &reference_data)
                };
                let var_p = mapper.c_to_p(&var_c, Some("MOCK"))?;
                assert_eq!(
                    (builder.variant_location()?, var_p.to_string().as_str()),
                    (location, hgvs_p),
                    "{hgvs_c} with ins_at_boundary_is_intronic = {ins_at_boundary_is_intronic}"
                );
            }
        }

        Ok(())
    }

    /// `NM_999989.1` with exons of 10, 30, and 10 bases: exon 1 is the 5' UTR, exon 2 the CDS
    /// `ATG AAA GCG AAA GCG TTT CGC GCG AAA TAG` (`MKAKAFRAK*`), and exon 3 the 3' UTR.
    /// `NM_999988.1` has the same sequence in a single exon.  `NM_999987.1` has the CDS
    /// without the stop codon, `MKAKAFRAK`, in a single exon.
    fn build_mapper_cds_boundary(ins_at_boundary_is_intronic: bool) -> Mapper {
        let record = |accession: &str, transcript_sequence: &str, cds_end_i, lengths| {
            sanity_mock::ProviderRecord {
                accession: accession.to_string(),
                transcript_sequence: transcript_sequence.to_string(),
                cds_start_i: 10,
                cds_end_i,
                translation_exceptions: Vec::new(),
                lengths,
            }
        };
        let sequence = "AAAATCAAACATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGGGGCCCCC";
        let provider = std::sync::Arc::new(sanity_mock::Provider::with_records(vec![
            record("NM_999989.1", sequence, 40, vec![10, 30, 10]),
            record("NM_999988.1", sequence, 40, vec![50]),
            record(
                "NM_999987.1",
                "AAAATCAAACATGAAAGCGAAAGCGTTTCGCGCGAAAGGGGGCCCCC",
                37,
                vec![47],
            ),
        ]));
        let config = Config {
            prevalidation_level: crate::validator::ValidationLevel::Intrinsic,
            strict_bounds: false,
            renormalize_g: false,
            ins_at_boundary_is_intronic,
            ..Default::default()
        };
        Mapper::new(&config, provider)
    }

    /// An insertion that can go before c.1 or after the stop codon changes no base of the CDS,
    /// so it gives `p.?` like a variant in the UTR.  With `ins_at_boundary_is_intronic` set to
    /// `false`, this includes the insertions and dups at the exon boundaries next to c.1 and
    /// the stop codon.
    #[test]
    fn cds_boundary_ins_c_to_p() -> Result<(), Error> {
        let exonic = build_mapper_cds_boundary(false);
        let intronic = build_mapper_cds_boundary(true);
        for (hgvs_c, expected) in [
            // In the UTRs.
            ("NM_999989.1:c.-3_-2insA", "MOCK:p.?"),
            ("NM_999989.1:c.-1dup", "MOCK:p.?"),
            ("NM_999989.1:c.*1_*2insA", "MOCK:p.?"),
            ("NM_999989.1:c.*1dup", "MOCK:p.?"),
            // At the start of the CDS.  In `NM_999989.1`, c.1 is at an exon boundary, in
            // `NM_999988.1` it is not.
            ("NM_999989.1:c.-1_1insGCC", "MOCK:p.?"),
            ("NM_999989.1:c.1-1_1insGCC", "MOCK:p.?"),
            ("NM_999989.1:c.1-1dupG", "MOCK:p.?"),
            ("NM_999989.1:c.1-1dup", "MOCK:p.?"),
            ("NM_999989.1:c.1-2_1-1dupAG", "MOCK:p.?"),
            ("NM_999988.1:c.-1_1insGCC", "MOCK:p.?"),
            // The copy can go before c.1.
            ("NM_999989.1:c.1dup", "MOCK:p.?"),
            ("NM_999989.1:c.1_3dup", "MOCK:p.?"),
            ("NM_999988.1:c.1dup", "MOCK:p.?"),
            ("NM_999988.1:c.1_3dup", "MOCK:p.?"),
            ("NM_999988.1:c.-2_2dup", "MOCK:p.?"),
            // Also if the dup starts before the transcript.  `NM_999988.1` starts with c.-10.
            // With `GG` before it, `c.-12_2dup` is `c.2_3insGGAAAATCAAACAT`.
            ("NM_999988.1:c.-12_2dup", "MOCK:p.?"),
            ("NM_999988.1:c.2_3insGGAAAATCAAACAT", "MOCK:p.?"),
            // The inserted A can go before c.1, like `c.1dup`.
            ("NM_999988.1:c.1_2insA", "MOCK:p.?"),
            // At the end of the CDS, likewise.
            ("NM_999989.1:c.30_*1insA", "MOCK:p.?"),
            ("NM_999989.1:c.30_30+1insA", "MOCK:p.?"),
            ("NM_999989.1:c.30+1dupG", "MOCK:p.?"),
            ("NM_999989.1:c.30+1dup", "MOCK:p.?"),
            ("NM_999989.1:c.30+1_30+2dupGT", "MOCK:p.?"),
            ("NM_999988.1:c.30_*1insA", "MOCK:p.?"),
            // The copy can go after the stop codon.
            ("NM_999989.1:c.30dup", "MOCK:p.?"),
            ("NM_999989.1:c.28_30dup", "MOCK:p.?"),
            ("NM_999988.1:c.30dup", "MOCK:p.?"),
            ("NM_999988.1:c.28_30dup", "MOCK:p.?"),
            // Also if the dup ends after the transcript.  `NM_999988.1` ends with c.*10.  With
            // `TT` after it, `c.29_*12dup` is `c.28_29insAGGGGGGCCCCCTT`.
            ("NM_999988.1:c.29_*12dup", "MOCK:p.?"),
            ("NM_999988.1:c.28_29insAGGGGGGCCCCCTT", "MOCK:p.?"),
            // The inserted G can go after the stop codon, like `c.30dup`.
            ("NM_999988.1:c.29_30insG", "MOCK:p.?"),
            // Without a stop codon, the CDS goes on after c.27, so the copy changes the CDS.
            // The amino acids after the last complete codon are unknown.
            ("NM_999987.1:c.26dup", "MOCK:p.?"),
            ("NM_999987.1:c.27_*1insA", "MOCK:p.?"),
            // Likewise for the bases after c.27, e.g., for `c.*1dup`, which is `c.27_*1insG`.
            ("NM_999987.1:c.*1dup", "MOCK:p.?"),
            ("NM_999987.1:c.*1_*3dup", "MOCK:p.?"),
            // These cannot go outside the CDS.
            ("NM_999989.1:c.1_2insGCC", "MOCK:p.Met1?"),
            ("NM_999989.1:c.2dup", "MOCK:p.Met1?"),
            ("NM_999989.1:c.29_30insC", "MOCK:p.Ter10TyrextTer?"),
            ("NM_999989.1:c.29dup", "MOCK:p.Ter10="),
            ("NM_999988.1:c.1_2insGCC", "MOCK:p.Met1?"),
            ("NM_999988.1:c.2dup", "MOCK:p.Met1?"),
            ("NM_999988.1:c.29_30insC", "MOCK:p.Ter10TyrextTer?"),
            ("NM_999988.1:c.29dup", "MOCK:p.Ter10="),
            // The copy forms the stop codon TAA in front of TAG.  A dup moves `cds_stop` like an
            // insertion, so every placement gives the same alt protein.
            ("NM_999989.1:c.26_28dup", "MOCK:p.="),
            ("NM_999989.1:c.27_29dup", "MOCK:p.="),
            ("NM_999989.1:c.25_26insAAT", "MOCK:p.="),
            ("NM_999989.1:c.28_29insAAT", "MOCK:p.="),
            ("NM_999989.1:c.29_30insATA", "MOCK:p.="),
            ("NM_999989.1:c.29_30insATAG", "MOCK:p.="),
            ("NM_999988.1:c.26_28dup", "MOCK:p.="),
            ("NM_999988.1:c.28_29insAAT", "MOCK:p.="),
        ] {
            check_c_to_p(&exonic, hgvs_c, Some("MOCK"), Some(expected))?;
            check_c_to_p(&intronic, hgvs_c, Some("MOCK"), Some(expected))?;
        }
        Ok(())
    }

    /// Like an insertion, a deletion, delins, or inversion that can go before c.1 or after the
    /// stop codon changes no base of the CDS, so it gives `p.?`.
    #[test]
    fn cds_boundary_del_c_to_p() -> Result<(), Error> {
        let exonic = build_mapper_cds_boundary(false);
        let intronic = build_mapper_cds_boundary(true);
        for tx_ac in ["NM_999989.1", "NM_999988.1"] {
            for (edit, expected) in [
                // c.-2_-1 is `AC`, so deleting `CA` gives the same sequence.
                ("c.-1_1del", "MOCK:p.?"),
                // Changes the C at c.-1 to a G and keeps the A at c.1.
                ("c.-1_1delinsGA", "MOCK:p.?"),
                // c.*1 is a G, like c.30.
                ("c.30del", "MOCK:p.?"),
                ("c.30_*1del", "MOCK:p.?"),
                // Changes the G at c.*1 to a C and keeps the G at c.30.
                ("c.30_*1delinsGC", "MOCK:p.?"),
                // c.-6 is a T, so the inversion keeps the A at c.1, like `c.-5_-1inv`.
                ("c.-6_1inv", "MOCK:p.?"),
                ("c.-5_-1inv", "MOCK:p.?"),
                // c.*6 is a C, so the inversion keeps the G at c.30, like `c.*1_*5inv`.
                ("c.30_*6inv", "MOCK:p.?"),
                ("c.*1_*5inv", "MOCK:p.?"),
                // These cannot go outside the CDS.
                ("c.1del", "MOCK:p.Met1?"),
                ("c.-1_1delinsGG", "MOCK:p.Met1?"),
                ("c.-5_1inv", "MOCK:p.Met1?"),
                ("c.29del", "MOCK:p.Ter10TrpextTer?"),
                ("c.28_30del", "MOCK:p.Ter10GlyextTer?"),
                // An identity changes no base.  c.7_8 is `GC`, so its inversion is one, too.
                ("c.30=", "MOCK:p.Ter10="),
                ("c.7_8inv", "MOCK:p.Ala3="),
            ] {
                let hgvs_c = format!("{tx_ac}:{edit}");
                check_c_to_p(&exonic, &hgvs_c, Some("MOCK"), Some(expected))?;
                check_c_to_p(&intronic, &hgvs_c, Some("MOCK"), Some(expected))?;
            }
        }
        // Without a stop codon, the CDS goes on after c.27, and the amino acids there are
        // unknown.  `c.27_*1delinsAC` keeps the A at c.27.
        for edit in [
            "c.*1G>A",
            "c.*1G=",
            "c.*1_*3del",
            "c.*1_*3inv",
            "c.27_*1delinsAC",
        ] {
            let hgvs_c = format!("NM_999987.1:{edit}");
            check_c_to_p(&exonic, &hgvs_c, Some("MOCK"), Some("MOCK:p.?"))?;
            check_c_to_p(&intronic, &hgvs_c, Some("MOCK"), Some("MOCK:p.?"))?;
        }
        Ok(())
    }

    mod gcp_tests {
        use std::path::Path;

        use anyhow::Error;

        #[derive(Debug, serde::Deserialize)]
        pub struct Record {
            pub id: String,
            #[serde(alias = "HGVSg")]
            pub hgvs_g: String,
            #[serde(alias = "HGVSc")]
            pub hgvs_c: String,
            #[serde(alias = "HGVSp")]
            pub hgvs_p: Option<String>,
            #[allow(dead_code)]
            pub description: Option<String>,
            #[allow(dead_code)]
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
            let result = format!("{}", var_p);
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
                format!("{}", NoRef(var))
            } else {
                format!("{var}")
            };
            let re = Regex::new(r"del\w+ins").expect("problem with regex in test");
            re.replace(&tmp, "delins").to_string()
        }

        let mapper = build_mapper()?;
        let records = gcp_tests::load_records(Path::new(path))?;
        let re = Regex::new(r"Ter\d+$").expect("problem with regex in test");

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

                let mut hgvs_p_test = format!("{}", var_p_test);

                if hgvs_p_exp.ends_with("Ter") {
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
