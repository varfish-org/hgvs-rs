//! Code for mapping variants between sequences.

use std::{ops::Range, rc::Rc};

use log::{debug, info};

use super::alignment::Mapper as AlignmentMapper;
use crate::{
    data::interface::Provider,
    normalizer::{self, Normalizer},
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
}

impl Default for Config {
    fn default() -> Self {
        Self {
            replace_reference: true,
            strict_validation: false,
            prevalidation_level: ValidationLevel::Full,
            add_gene_symbol: false,
            strict_bounds: true,
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
    ) -> Result<AlignmentMapper, anyhow::Error> {
        // TODO: implement caching
        AlignmentMapper::new(self.provider.clone(), tx_ac, alt_ac, alt_aln_method)
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
        var_t: &HgvsVariant,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<HgvsVariant, anyhow::Error> {
        self.validator.validate(var_t)?;
        let var_t = self.replace_reference(var_t.clone())?;
        match var_t {
            HgvsVariant::TxVariant { .. } => self.n_to_g(&var_t, alt_ac, alt_aln_method),
            HgvsVariant::CdsVariant { .. } => self.c_to_g(&var_t, alt_ac, alt_aln_method),
            _ => {
                return Err(anyhow::anyhow!(
                    "Expected transcript or CDS variant but received {}",
                    &var_t
                ))
            }
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
    ) -> Result<HgvsVariant, anyhow::Error> {
        self.validator.validate(var_g)?;
        let var_g = self.replace_reference(var_g.clone())?;
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
            {
                info!("Renormalizing out-of-bounds minus strand variant on genomic sequenc");
                Normalizer::new(
                    self,
                    self.provider.clone(),
                    self.validator.clone(),
                    Default::default(),
                )
                .normalize(&var_g)?
            } else {
                var_g.clone()
            };

            let pos_n = mapper.g_to_n(loc_edit.loc.inner())?;
            let pos_n = Mu::from(
                pos_n.inner(),
                loc_edit.loc.is_certain() && pos_n.is_certain(),
            );
            let (pos_n, edit_n) = if let Mu::Certain(pos_n) = pos_n {
                let edit_n = self.convert_edit_check_strand(mapper.strand, &loc_edit.edit)?;
                if let NaEdit::Ins { .. } = edit_n.inner() {
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
                            Mu::Certain(NaEdit::Ins {
                                alternative: "".to_string(),
                            }),
                        )
                    } else {
                        (Mu::Certain((*pos_n).clone()), edit_n)
                    }
                } else {
                    (Mu::Certain((*pos_n).clone()), edit_n)
                }
            } else {
                let pos_g = mapper.n_to_g(pos_n.inner())?;
                let edit_n = NaEdit::RefAlt {
                    reference: "".to_string(),
                    alternative: self.get_altered_sequence(
                        mapper.strand,
                        pos_g.inner().clone().try_into()?,
                        &var_g,
                    )?,
                };
                (Mu::Uncertain((*pos_n.inner()).clone()), Mu::Certain(edit_n))
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
            Err(anyhow::anyhow!(
                "Expected a GenomeVariant but received {}",
                &var_g
            ))
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
    ) -> Result<HgvsVariant, anyhow::Error> {
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
                if let (NaEdit::Ins { .. }, Some(end), Some(start)) =
                    (edit_g.inner(), pos_g.end, pos_g.start)
                {
                    if end - start > 1 {
                        (
                            Mu::Certain(GenomeInterval {
                                start: Some(start + 1),
                                end: Some(end - 1),
                            }),
                            Mu::from(
                                NaEdit::Ins {
                                    alternative: "".to_string(),
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
            Err(anyhow::anyhow!(
                "Expected a TxVariant but received {}",
                &var_n
            ))
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
    ) -> Result<HgvsVariant, anyhow::Error> {
        self.validator.validate(var_g)?;
        let var_g = self.replace_reference(var_g.clone())?;
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
                if let NaEdit::Ins { .. } = edit_c.inner() {
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
                            Mu::Certain(NaEdit::Ins {
                                alternative: "".to_string(),
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
                (Mu::Uncertain((*pos_c.inner()).clone()), Mu::Certain(edit_c))
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
            Err(anyhow::anyhow!(
                "Expected a GenomeVariant but received {}",
                &var_g
            ))
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
    ) -> Result<HgvsVariant, anyhow::Error> {
        self.validator.validate(var_c)?;
        let var_c = self.replace_reference(var_c.clone())?;
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
                if let (NaEdit::Ins { .. }, Some(end), Some(start)) =
                    (edit_g.inner(), pos_g.end, pos_g.start)
                {
                    if end - start > 1 {
                        (
                            Mu::Certain(GenomeInterval {
                                start: Some(start + 1),
                                end: Some(end - 1),
                            }),
                            Mu::from(
                                NaEdit::Ins {
                                    alternative: "".to_string(),
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
                let edit_n = NaEdit::RefAlt {
                    reference: "".to_string(),
                    alternative: self.get_altered_sequence(
                        mapper.strand,
                        pos_n.inner().clone().into(),
                        &var_c,
                    )?,
                };
                (pos_g, Mu::Certain(edit_n))
            };

            let var_g = HgvsVariant::GenomeVariant {
                accession: accession.clone(),
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
            Err(anyhow::anyhow!(
                "Expected a CdsVariant but received {}",
                &var_c
            ))
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
    ) -> Result<HgvsVariant, anyhow::Error> {
        self.validator.validate(var_t)?;
        let var_t = self.replace_reference(var_t.clone())?;
        match var_t {
            HgvsVariant::TxVariant { .. } => self.n_to_g(&var_t, alt_ac, alt_aln_method),
            HgvsVariant::CdsVariant { .. } => self.c_to_g(&var_t, alt_ac, alt_aln_method),
            _ => {
                return Err(anyhow::anyhow!(
                    "Expected transcript or CDS variant but received {}",
                    &var_t
                ))
            }
        }
    }

    /// Convert from CDS variant (c.) to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::CdsVariant` to project
    pub fn c_to_n(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        log::debug!("c_to_n({})", var_c);
        self.validator.validate(var_c)?;
        let var_c = self.replace_reference(var_c.clone())?;
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
            Err(anyhow::anyhow!(
                "Expected a CdsVariant but received {}",
                &var_c
            ))
        }
    }

    /// Convert from transcript variant (n.) to CDS variant (c.).
    ///
    /// # Args
    ///
    /// * `var_n` -- `HgvsVariant::TxVariant` to project
    pub fn n_to_c(&self, var_n: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        self.validator.validate(var_n)?;
        let var_n = self.replace_reference(var_n.clone())?;
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
            Err(anyhow::anyhow!(
                "Expected a TxVariant but received {}",
                &var_n
            ))
        }
    }

    /// Convert from CDS variant (c.) to protein variant (p.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::TxVariant` to project
    /// * `pro_ac` -- Protein accession
    pub fn c_to_p(
        &self,
        var_c: &HgvsVariant,
        prot_ac: Option<&str>,
    ) -> Result<HgvsVariant, anyhow::Error> {
        use super::altseq::*;

        if let HgvsVariant::CdsVariant {
            accession,
            gene_symbol: _,
            loc_edit: _,
        } = &var_c
        {
            self.validator.validate(var_c)?;

            let reference_data = RefTranscriptData::new(
                self.provider.clone(),
                accession.deref(),
                prot_ac.map(|s| s.to_string()).as_deref(),
            )?;
            let builder = AltSeqBuilder::new(var_c.clone(), reference_data.clone());

            // NB: the following comment is from the original code.
            // TODO: handle case where you get 2+ alt sequences back;  currently get lis tof 1 element
            // loop structure implemented to handle this, but doesn't really do anything currently.

            let var_ps: Result<Vec<_>, anyhow::Error> = builder
                .build_altseq()?
                .into_iter()
                .map(|alt_data| {
                    let builder = AltSeqToHgvsp::new(reference_data.clone(), alt_data);
                    builder.build_hgvsp()
                })
                .collect();
            let var_p = var_ps?.into_iter().next().unwrap();

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
                return Err(anyhow::anyhow!(
                    "This cannot happen; must have ProtVariant here."
                ));
            };

            Ok(var_p)
        } else {
            Err(anyhow::anyhow!(
                "Expected a CdsVariant but received {}",
                &var_c
            ))
        }
    }

    fn get_altered_sequence(
        &self,
        strand: i16,
        interval: Range<i32>,
        var: &HgvsVariant,
    ) -> Result<String, anyhow::Error> {
        let mut seq = self.provider.as_ref().get_seq_part(
            var.accession(),
            Some(interval.start.try_into()?),
            Some(interval.end.try_into()?),
        )?;

        let r = var.loc_range().ok_or(anyhow::anyhow!(
            "Cannot get altered sequence for missing positions"
        ))?;
        let r = ((r.start - interval.start) as usize)..((r.end - interval.start) as usize);

        let na_edit = var
            .na_edit()
            .ok_or(anyhow::anyhow!("Variant is missing nucleic acid edit"))?;

        match na_edit {
            NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => {
                seq.replace_range(r, alternative)
            }
            NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => seq.replace_range(r, ""),
            NaEdit::Ins { alternative } => seq.replace_range(r, alternative),
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
    ) -> Result<Mu<NaEdit>, anyhow::Error> {
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
    pub(crate) fn replace_reference(&self, var: HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        log::debug!("replace_reference({})", &var);
        match &var {
            HgvsVariant::ProtVariant { .. } => Err(anyhow::anyhow!(
                "Can only update reference for c, g, m, n, r"
            )),
            _ => Ok(()),
        }?;

        log::debug!("???");

        if let Some(NaEdit::Ins { .. }) = var.na_edit() {
            // Insertions have no reference sequence (zero-width); return as-is.
            return Ok(var);
        }

        log::debug!("111");

        if var.spans_intron() {
            debug!(
                "Can't update reference sequence for intronic variant {}",
                var
            );
            return Ok(var);
        }

        log::debug!("???");

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

        log::debug!("foo");

        let na_edit = var
            .na_edit()
            .expect("Variant must be of nucleic acid type here");
        if !na_edit.reference_equals(&seq) {
            debug!("Replaced reference sequence in {} with {}", &var, &seq);
            Ok(var.with_reference(seq))
        } else {
            Ok(var)
        }
    }

    fn fetch_gene_symbol(
        &self,
        tx_ac: &str,
        gene_symbol: &Option<GeneSymbol>,
    ) -> Result<Option<GeneSymbol>, anyhow::Error> {
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
    use std::{rc::Rc, str::FromStr};
    use test_log::test;

    use crate::{
        data::uta::{Config as ProviderConfig, Provider},
        parser::HgvsVariant,
    };

    use super::{Config, Mapper};

    fn get_config() -> ProviderConfig {
        ProviderConfig {
            db_url: std::env::var("TEST_UTA_DATABASE_URL")
                .expect("Environment variable TEST_UTA_DATABASE_URL undefined!"),
            db_schema: std::env::var("TEST_UTA_DATABASE_SCHEMA")
                .expect("Environment variable TEST_UTA_DATABASE_SCHEMA undefined!"),
        }
    }

    fn build_mapper() -> Result<Mapper, anyhow::Error> {
        let provider = Rc::new(Provider::with_config(&get_config())?);
        let config = Config::default();
        Ok(Mapper::new(&config, provider))
    }

    #[test]
    fn fail_for_invalid_variant_types() -> Result<(), anyhow::Error> {
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
    fn fail_c_to_p_on_invalid_nm_accession() -> Result<(), anyhow::Error> {
        let mapper = build_mapper()?;

        let hgvs_g = "NC_000007.13:g.36561662C>T";
        let var_g = HgvsVariant::from_str(hgvs_g)?;

        assert!(mapper.c_to_p(&var_g, Some("NM_999999.1")).is_err());

        Ok(())
    }

    #[test]
    fn fail_on_undefined_cds() -> Result<(), anyhow::Error> {
        let mapper = build_mapper()?;

        let hgvs_n = "NR_111984.1:n.44G>A"; // legit
        let hgvs_c = "NR_111984.1:c.44G>A"; // bogus: c. with non-coding tx accession

        let var_n = HgvsVariant::from_str(hgvs_n)?;
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let _tx_ac = if let HgvsVariant::TxVariant { accession, .. } = &var_c {
            Some(accession.value.clone())
        } else {
            None
        }
        .unwrap();

        // n_to_c: transcript is non-coding
        assert!(mapper.n_to_c(&var_n).is_err());

        // c_to_n: var_c is bogus
        assert!(mapper.c_to_n(&var_c).is_err());

        Ok(())
    }

    #[test]
    fn map_var_of_unsupported_validation() -> Result<(), anyhow::Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_003777.3:c.13552_*36del57"; // gene DNAH11
        let var_c = HgvsVariant::from_str(hgvs_c)?;

        let var_g = mapper.c_to_g(&var_c, "NC_000007.13", "splign")?;
        assert_eq!(format!("{}", &var_g), "NC_000007.13:g.21940852_21940908del");

        Ok(())
    }

    #[test]
    fn map_to_unknown_p_effect() -> Result<(), anyhow::Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_020975.4:c.625+9C>T"; // gene RET
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_066124.1:p.?");

        Ok(())
    }

    #[test]
    fn map_of_c_out_of_cds_bound() -> Result<(), anyhow::Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_145901.2:c.343T>C"; // gene HMGA1
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        assert!(mapper.c_to_p(&var_c, None).is_err());

        Ok(())
    }

    #[test]
    fn map_of_dup_at_cds_end() -> Result<(), anyhow::Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_001051.2:c.1257dupG"; // gene SSTR3
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        let var_p = mapper.c_to_p(&var_c, None)?;
        assert_eq!(format!("{}", &var_p), "NP_001042.1:p.=)");

        Ok(())
    }

    #[test]
    fn map_of_c_out_of_reference_bound() -> Result<(), anyhow::Error> {
        let mapper = build_mapper()?;
        let hgvs_c = "NM_000249.3:c.-73960_*46597del"; // gene MLH1
        let var_c = HgvsVariant::from_str(hgvs_c)?;
        assert!(mapper.c_to_p(&var_c, None).is_err());

        Ok(())
    }

    /// The following tests corresponds to the `test_hgvs_variantmapper_cp_sanity.py`
    /// test suite of the Python package.  It uses a mock data provider, defined
    /// in the `sanity_mock` module.

    mod sanity_mock {
        use std::{
            path::{Path, PathBuf},
            rc::Rc,
        };



        use crate::{data::interface::Provider as ProviderInterface};
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
            pub fn new(path: &Path) -> Result<Self, anyhow::Error> {
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
            ) -> Result<crate::data::interface::GeneInfoRecord, anyhow::Error> {
                panic!("for test use only");
            }

            fn get_pro_ac_for_tx_ac(&self, _tx_ac: &str) -> Result<Option<String>, anyhow::Error> {
                panic!("for test use only");
            }

            fn get_seq_part(
                &self,
                tx_ac: &str,
                begin: Option<usize>,
                end: Option<usize>,
            ) -> Result<String, anyhow::Error> {
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
                Err(anyhow::anyhow!("Found no record for accession {}", &tx_ac))
            }

            fn get_acs_for_protein_seq(&self, _seq: &str) -> Result<Vec<String>, anyhow::Error> {
                panic!("for test use only");
            }

            fn get_similar_transcripts(
                &self,
                _tx_ac: &str,
            ) -> Result<Vec<crate::data::interface::TxSimilarityRecord>, anyhow::Error>
            {
                panic!("for test use only");
            }

            fn get_tx_exons(
                &self,
                _tx_ac: &str,
                _alt_ac: &str,
                _alt_aln_method: &str,
            ) -> Result<Vec<crate::data::interface::TxExonsRecord>, anyhow::Error> {
                todo!()
            }

            fn get_tx_for_gene(
                &self,
                _gene: &str,
            ) -> Result<Vec<crate::data::interface::TxInfoRecord>, anyhow::Error> {
                panic!("for test use only");
            }

            fn get_tx_for_region(
                &self,
                _alt_ac: &str,
                _alt_aln_method: &str,
                _start_i: i32,
                _end_i: i32,
            ) -> Result<Vec<crate::data::interface::TxForRegionRecord>, anyhow::Error> {
                panic!("for test use only");
            }

            fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, anyhow::Error> {
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
                Err(anyhow::anyhow!("Found no record for accession {}", &tx_ac))
            }

            fn get_tx_info(
                &self,
                _tx_ac: &str,
                _alt_ac: &str,
                _alt_aln_method: &str,
            ) -> Result<crate::data::interface::TxInfoRecord, anyhow::Error> {
                panic!("for test use only");
            }

            fn get_tx_mapping_options(
                &self,
                _tx_ac: &str,
            ) -> Result<Vec<crate::data::interface::TxMappingOptionsRecord>, anyhow::Error>
            {
                panic!("for test use only");
            }
        }

        pub fn build_mapper() -> Result<Mapper, anyhow::Error> {
            let path = PathBuf::from("tests/data/mapper/sanity_cp.tsv");
            let provider = Rc::new(Provider::new(&path)?);
            let config = Config::default();
            Ok(Mapper::new(&config, provider))
        }
    }

    fn test_hgvs_c_to_p_conversion(hgvsc: &str, hgvsp_expected: &str) -> Result<(), anyhow::Error> {
        let mapper = sanity_mock::build_mapper()?;

        let var_c = HgvsVariant::from_str(hgvsc)?;
        let ac_p = "MOCK";

        let var_p = mapper.c_to_p(&var_c, Some(ac_p))?;
        let hgvsp_actual = format!("{}", &var_p);

        assert_eq!(hgvsp_actual, hgvsp_expected);

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_silent() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.6A>G";
        let hgvsp_expected = "MOCK:p.Lys2=";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_substitution() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.6A>T";
        let hgvsp_expected = "MOCK:p.Lys2Asn";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_substitution_introduces_stop_codon() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999996.1:c.8C>A";
        let hgvsp_expected = "MOCK:p.Ser3Ter";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_substitution_removes_stop_codon() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999998.1:c.30G>T";
        let hgvsp_expected = "MOCK:p.Ter10TyrextTer3";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    //xx
    #[test]
    fn hgvs_c_to_p_insertion_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.6_7insGGG";
        let hgvsp_expected = "MOCK:p.Lys2_Ala3insGly";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_insertion_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.22_23insT";
        let hgvsp_expected = "MOCK:p.Ala8ValfsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_adds_stop() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.8_9insTT";
        let hgvsp_expected = "MOCK:p.Lys4Ter";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.10_12del";
        let hgvsp_expected = "MOCK:p.Lys4del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion2_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.4_15del";
        let hgvsp_expected = "MOCK:p.Lys2_Ala5del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion3_no_frameshift_c_term() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999995.1:c.4_6del";
        let hgvsp_expected = "MOCK:p.Lys3del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion4_no_frameshift_c_term() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999994.1:c.4_9del";
        let hgvsp_expected = "MOCK:p.Lys3_Lys4del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion5_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999994.1:c.20_25del";
        let hgvsp_expected = "MOCK:p.Ala7_Arg9delinsGly";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion6_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.5_7del";
        let hgvsp_expected = "MOCK:p.Lys2_Ala3delinsThr";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion7_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999993.1:c.13_24del";
        let hgvsp_expected = "MOCK:p.Arg5_Ala8del";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_frameshift_nostop() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.11_12del";
        let hgvsp_expected = "MOCK:p.Lys4SerfsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_frameshift_adds_stop() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999997.1:c.7del";
        let hgvsp_expected = "MOCK:p.Ala3ArgfsTer6";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_no_frameshift_removes_stop_plus_previous() -> Result<(), anyhow::Error>
    {
        let hgvsc = "NM_999999.1:c.25_30del";
        let hgvsp_expected = "MOCK:p.Lys9_Ter10delinsGly";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_indel_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.11_12delinsTCCCA";
        let hgvsp_expected = "MOCK:p.Lys4delinsIlePro";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_indel2_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.11_18delinsTCCCA";
        let hgvsp_expected = "MOCK:p.Lys4_Phe6delinsIlePro";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_indel_frameshift_nostop() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.8delinsGG";
        let hgvsp_expected = "MOCK:p.Ala3GlyfsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_1aa_no_frameshift_2() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.10_12dup";
        let hgvsp_expected = "MOCK:p.Lys4dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_1aa_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.16_18dup";
        let hgvsp_expected = "MOCK:p.Phe6dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_2aa_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.16_21dup";
        let hgvsp_expected = "MOCK:p.Phe6_Arg7dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_2aa2_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999995.1:c.4_6dup";
        let hgvsp_expected = "MOCK:p.Lys3dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_3aa_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.16_24dup";
        let hgvsp_expected = "MOCK:p.Phe6_Ala8dup";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_dup_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.12_13dup";
        let hgvsp_expected = "MOCK:p.Ala5GlufsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_intron() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.12+1G>A";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_five_prime_utr() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.-2A>G";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_three_prime_utr() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.*3G>A";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_into_three_prime_utr_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.27_*3del";
        let hgvsp_expected = "MOCK:p.Lys9XaafsTer?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_into_three_prime_utr_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999995.1:c.28_*3del";
        let hgvsp_expected = "MOCK:p.Lys10_Ter11delinsArgGlnPheArg";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_delins_into_three_prime_utr_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999995.1:c.28_*3delinsGGG";
        let hgvsp_expected = "MOCK:p.Lys10_Ter11delinsGlyArgGlnPheArg";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    /// See recommendations re p.? (p.Met1?) at:
    /// http://varnomen.hgvs.org/recommendations/protein/variant/substitution/
    #[test]
    fn hgvs_c_to_p_substitution_removes_start_codon() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.1A>G";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_from_five_prime_utr_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.-3_1del";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_deletion_from_five_prime_utr_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.-3_3del";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_delins_from_five_prime_utr_no_frameshift() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.-3_3delinsAAA";
        let hgvsp_expected = "MOCK:p.Met1?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_delete_entire_gene() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999999.1:c.-3_*1del";
        let hgvsp_expected = "MOCK:p.0?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }

    #[test]
    fn hgvs_c_to_p_multiple_stop_codons() -> Result<(), anyhow::Error> {
        let hgvsc = "NM_999992.1:c.4G>A";
        let hgvsp_expected = "MOCK:p.?";
        test_hgvs_c_to_p_conversion(hgvsc, hgvsp_expected)?;

        Ok(())
    }
}
