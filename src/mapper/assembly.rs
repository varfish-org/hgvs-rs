//! Implementation of `AssemblyMapper`.

use std::collections::{HashMap, HashSet};
use std::ops::Range;

use std::sync::Arc;

use crate::mapper::error::Error;
use crate::mapper::variant::{self, ShiftOverBoundaryPreference};
use crate::normalizer::Direction;
use crate::parser::HgvsVariant;
use crate::{data::interface::Provider, validator::ValidationLevel};

#[derive(Debug, PartialEq, Eq, Default, Clone, Copy)]
pub enum InParAssume {
    #[default]
    X,
    Y,
    None,
}

impl InParAssume {
    fn chrom_name(&self) -> &str {
        match &self {
            InParAssume::X => "X",
            InParAssume::Y => "Y",
            InParAssume::None => panic!("no chromosome name"),
        }
    }
}

/// Configuration for `Assemblymapper`.
#[derive(Debug)]
pub struct Config {
    /// Assembly name, e.g., "GRCh38".
    pub assembly: String,
    /// Optional custom assembly map from accession to chromosome name.
    /// If provided, this bypasses the provider's `get_assembly_map` method.
    pub asm_map: Option<HashMap<String, String>>,
    pub alt_aln_method: String,
    pub normalize: bool,
    pub in_par_assume: InParAssume,
    pub prevalidation_level: ValidationLevel,
    pub replace_reference: bool,
    pub strict_validation: bool,
    pub strict_bounds: bool,
    pub add_gene_symbol: bool,
    /// Always re-normalize genome variants during g-to-n projections.
    /// This can be switched off so genome sequence does not have to be available in provider.
    pub renormalize_g: bool,
    /// Use the genome sequence in case of uncertain g-to-n projections.  This
    /// can be switched off so genome sequence does not have to be available.
    pub genome_seq_available: bool,
    pub shuffle_direction: Direction,
    pub window_size: usize,
    /// See `variant::Config::ins_at_boundary_is_intronic`.
    pub ins_at_boundary_is_intronic: bool,
    /// See `variant::Config::shift_over_boundary`.  This needs `genome_seq_available`.
    pub shift_over_boundary: bool,
    /// See `variant::Config::shift_over_boundary_preference`.
    pub shift_over_boundary_preference: ShiftOverBoundaryPreference,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            assembly: "GRCh38".to_string(),
            asm_map: None,
            alt_aln_method: "splign".to_string(),
            normalize: true,
            in_par_assume: InParAssume::X,
            prevalidation_level: ValidationLevel::Intrinsic,
            replace_reference: true,
            strict_validation: true,
            strict_bounds: true,
            add_gene_symbol: false,
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

/// Provides simplified variant mapping for a single assembly and transcript-reference
/// alignment method.
///
/// `AssemblyMapper` uses `variant::Mapper`, to provide all projection functionality, and add:
///
/// * Automatic selection of genomic sequence accession
/// * Transcript selection from genomic coordinates
/// * Normalization after projection
/// * Special handling for PAR regions
///
/// `AssemblyMapper` is instantiated with an assembly name and `alt_aln_method`. These enable
/// the following conveniences over `variant::Mapper`:
///
/// * The assembly and alignment method are used to automatically select an appropriate
///   chromosomal reference sequence when mapping from a transcript to a genome (i.e.,
///   `c_to_g()` and `n_to_g()`).
///
/// * A new method, relevant_trancripts(g_variant), returns a list of transcript accessions
///   available for the specified variant. These accessions are candidates mapping from
///   genomic to trancript coordinates (i.e., g_to_c(...) and g_to_n(...)).
///
/// Note: AssemblyMapper supports only chromosomal references (e.g. NC_000006.11). It does
/// not support contigs or other genomic sequences (e.g., NT_167249.1).
pub struct Mapper {
    pub config: Config,
    provider: Arc<dyn Provider + Send + Sync>,
    pub inner: variant::Mapper,
    /// Accessions of contigs in assembly.
    asm_accessions: HashSet<String>,
    /// Map from accession to contig name.
    asm_map: HashMap<String, String>,
}

impl Mapper {
    /// Construct new assembly mapper from config and provider.
    pub fn new(config: Config, provider: Arc<dyn Provider + Send + Sync>) -> Self {
        let inner_config = variant::Config {
            replace_reference: config.replace_reference,
            strict_validation: config.strict_validation,
            prevalidation_level: config.prevalidation_level,
            add_gene_symbol: config.add_gene_symbol,
            strict_bounds: config.strict_bounds,
            renormalize_g: config.renormalize_g,
            genome_seq_available: config.genome_seq_available,
            shuffle_direction: config.shuffle_direction,
            window_size: config.window_size,
            ins_at_boundary_is_intronic: config.ins_at_boundary_is_intronic,
            shift_over_boundary: config.shift_over_boundary,
            shift_over_boundary_preference: config.shift_over_boundary_preference,
        };
        let inner = variant::Mapper::new(&inner_config, provider.clone());
        let asm_map = config.asm_map.clone().unwrap_or_else(|| {
            provider
                .as_ref()
                .get_assembly_map(&config.assembly)
                .unwrap_or_else(|_| {
                    panic!("assembly map for {} should be available", config.assembly)
                })
                .into_iter()
                .map(|(key, value)| (key.clone(), value.clone()))
                .collect()
        });

        let asm_accessions = asm_map.keys().cloned().collect::<HashSet<_>>();

        Self {
            config,
            provider,
            inner,
            asm_accessions,
            asm_map,
        }
    }

    pub fn variant_mapper(&self) -> &variant::Mapper {
        &self.inner
    }

    /// Convert from genome (g.) variant to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    pub fn g_to_n(&self, var_g: &HgvsVariant, tx_ac: &str) -> Result<HgvsVariant, Error> {
        let var = self
            .inner
            .g_to_n(var_g, tx_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from genome (g.) variant to CDS variant (c.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    pub fn g_to_c(&self, var_g: &HgvsVariant, tx_ac: &str) -> Result<HgvsVariant, Error> {
        let var = self
            .inner
            .g_to_c(var_g, tx_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from genome (g.) variant to transcript variant (g. or n.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    pub fn g_to_t(&self, var_g: &HgvsVariant, tx_ac: &str) -> Result<HgvsVariant, Error> {
        let var = self
            .inner
            .g_to_t(var_g, tx_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from transcript variant (g.) variant to genome variant (n.).
    ///
    /// # Args
    ///
    /// * `var_n` -- `HgvsVariant::TxVariant` to project
    pub fn n_to_g(&self, var_n: &HgvsVariant) -> Result<HgvsVariant, Error> {
        let alt_ac = self.alt_ac_for_tx_ac(var_n.accession())?;
        let var = self
            .inner
            .n_to_g(var_n, &alt_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from CDS variant (c.) to genome variant (g.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::CdsVariant` to project
    /// * `alt_ac` -- alternative contig accession
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn c_to_g(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, Error> {
        let alt_ac = self.alt_ac_for_tx_ac(var_c.accession())?;
        let var = self
            .inner
            .c_to_g(var_c, &alt_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from transcript (c. or n.) to genome (g.) variant.
    ///
    /// # Args
    ///
    /// * `var_t` -- `HgvsVariant::TxVariant` or `HgvsVariant::CdsVariant` to project
    /// * `alt_ac` -- accession of alternativ esequence
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn t_to_g(&self, var_t: &HgvsVariant) -> Result<HgvsVariant, Error> {
        let alt_ac = self.alt_ac_for_tx_ac(var_t.accession())?;
        let var = self
            .inner
            .t_to_g(var_t, &alt_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from CDS variant (c.) to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::CdsVariant` to project
    pub fn c_to_n(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, Error> {
        let var = self.inner.c_to_n(var_c)?;
        self.maybe_normalize(&var)
    }

    /// Convert from transcript variant (n.) to CDS variant (c.).
    ///
    /// # Args
    ///
    /// * `var_n` -- `HgvsVariant::TxVariant` to project
    pub fn n_to_c(&self, var_n: &HgvsVariant) -> Result<HgvsVariant, Error> {
        let var = self.inner.n_to_c(var_n)?;
        self.maybe_normalize(&var)
    }

    /// Convert from CDS variant (c.) to protein variant (p.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::TxVariant` to project
    pub fn c_to_p(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, Error> {
        // Only the genome has the reference of intron bases, e.g., of a dup next to an exon.
        // Shifting an insertion or duplication over an exon boundary needs the genome, too.
        let shifts = self.config.shift_over_boundary
            && self.config.shift_over_boundary_preference != ShiftOverBoundaryPreference::Default
            && var_c
                .na_edit()
                .is_some_and(|edit| edit.is_ins() || edit.is_dup());
        let alt_ac = if self.config.genome_seq_available && (var_c.spans_intron() || shifts) {
            Some(self.alt_ac_for_tx_ac(var_c.accession())?)
        } else {
            None
        };
        let alt = alt_ac
            .as_deref()
            .map(|alt_ac| (alt_ac, self.config.alt_aln_method.as_str()));
        let var = self.inner.c_to_p_with_alt(var_c, None, alt)?;
        self.maybe_normalize(&var)
    }

    /// Obtain relevant transcript accessions.
    ///
    /// # Args
    ///
    /// * `var_g` -- genome variant to obtain the transcript accessions for.
    ///
    /// # Returns
    ///
    /// `Vec` of relevant transcript accessions.
    pub fn relevant_transcripts(&self, var_g: &HgvsVariant) -> Result<Vec<String>, Error> {
        match var_g {
            HgvsVariant::GenomeVariant { loc_edit, .. } => {
                let r: Range<i32> = loc_edit.loc.inner().clone().try_into()?;
                Ok(self
                    .provider
                    .as_ref()
                    .get_tx_for_region(
                        var_g.accession(),
                        &self.config.alt_aln_method,
                        r.start,
                        r.end,
                    )?
                    .into_iter()
                    .map(|rec| rec.tx_ac)
                    .collect::<Vec<_>>())
            }
            _ => Err(Error::NotGenomeVariant(format!("{0}", var_g))),
        }
    }

    /// Normalize variant if requested and ignore errors.  This is better than checking whether
    /// the variant is intronic because future UTAs will support LRG, which will enable checking
    /// intronic variants.
    ///
    /// Includes a special round-trip normalization for intronic variants,
    /// similar to how VariantValidator does it.
    pub fn maybe_normalize(&self, var: &HgvsVariant) -> Result<HgvsVariant, Error> {
        if self.config.normalize {
            // spans_intron() checks offsets and returns false for g., m., p.
            // and returns true for c., n. (and r.) (depending on offsets).
            if var.spans_intron() {
                match self.normalize_intronic(var) {
                    Ok(norm_var) => return Ok(norm_var),
                    Err(e) => {
                        tracing::warn!("Intronic normalization of variant {} failed: {}", &var, e);
                    }
                }
            }

            // otherwise (or if intronic normalization failed): do default normalization
            let normalizer = self.inner.normalizer()?;
            normalizer
                .normalize(var)
                .map_err(|e| tracing::warn!("Normalization of variant {} failed: {}", &var, e))
                .or_else(|_| Ok(var.clone()))
        } else {
            Ok(var.clone())
        }
    }

    /// Normalize intronic variants via genomic round-trip (within bounds).
    fn normalize_intronic(&self, var: &HgvsVariant) -> Result<HgvsVariant, Error> {
        // get genomic variant
        let alt_ac = self.alt_ac_for_tx_ac(var.accession())?;
        let var_g = self
            .inner
            .t_to_g(var, &alt_ac, &self.config.alt_aln_method)?;

        // determine strand & intron bounds
        let exons =
            self.provider
                .get_tx_exons(var.accession(), &alt_ac, &self.config.alt_aln_method)?;
        let strand = exons
            .first()
            .ok_or_else(|| {
                Error::NoExons(
                    var.accession().to_string(),
                    alt_ac.clone(),
                    self.config.alt_aln_method.clone(),
                )
            })?
            .alt_strand;

        // get genomic start/end of variant
        let var_range = var_g
            .loc_range()
            .ok_or_else(|| Error::MissingGenomeIntervalPosition(format!("{}", var_g)))?;
        let var_start = var_range.start;
        let var_end = var_range.end;

        // find intron "containing" the variant to set boundaries
        let mut boundary = None;

        debug_assert!(
            exons.is_sorted_by_key(|e| e.alt_start_i),
            "Exons are not sorted by alt_start_i: {:?} (alt_ac={})",
            exons,
            alt_ac
        );

        for w in exons.windows(2) {
            let prev_exon_end = w[0].alt_end_i;
            let next_exon_start = w[1].alt_start_i;

            // check if variant is roughly within this gap
            if var_start >= prev_exon_end && var_end <= next_exon_start {
                boundary = Some(prev_exon_end..next_exon_start);
                break;
            }
        }

        let boundary = boundary.ok_or_else(|| Error::NoIntronBoundary(var.to_string()))?;

        // normalize on genomic sequence within bounds
        let norm_g = if strand == 1 {
            self.inner
                .right_normalizer()?
                .normalize_bounded(&var_g, boundary)?
        } else {
            self.inner
                .left_normalizer()?
                .normalize_bounded(&var_g, boundary)?
        };

        // map back to transcript
        match var {
            HgvsVariant::CdsVariant { .. } => {
                self.inner
                    .g_to_c(&norm_g, var.accession(), &self.config.alt_aln_method)
            }
            HgvsVariant::TxVariant { .. } => {
                self.inner
                    .g_to_n(&norm_g, var.accession(), &self.config.alt_aln_method)
            }
            _ => Err(Error::ExpectedTxVariant(format!("{}", var))),
        }
    }

    /// Return chromosomal accession for the given transcript accession and the assembly and
    /// aln_method from configuration.
    fn alt_ac_for_tx_ac(&self, tx_ac: &str) -> Result<String, Error> {
        let mut alt_acs = Vec::new();
        for record in self.provider.get_tx_mapping_options(tx_ac)? {
            if record.alt_aln_method == self.config.alt_aln_method
                && self.asm_accessions.contains(&record.alt_ac)
            {
                alt_acs.push(record.alt_ac);
            }
        }

        if alt_acs.is_empty() {
            return Err(Error::NoAlignments(
                tx_ac.to_string(),
                self.config.assembly.clone(),
                self.config.alt_aln_method.clone(),
            ));
        }

        if alt_acs.len() == 1 {
            Ok(alt_acs
                .first()
                .expect("should not happen; checked for one alt_ac above")
                .to_owned())
        } else {
            // alt_acs.len() > 1
            // Perform sanity check on result if more than one contig is returned.
            let mut alts = alt_acs
                .iter()
                .map(|ac| self.asm_map.get(ac))
                .filter(|x| x.is_some())
                .map(|x| x.expect("should not happen; empty alt_ac").clone())
                .collect::<Vec<_>>();
            alts.sort();
            if alts.join("") != "XY" {
                return Err(Error::MultipleChromAlignsNonPar(
                    tx_ac.to_string(),
                    self.config.assembly.clone(),
                    self.config.alt_aln_method.to_string(),
                    format!("{:?}", alts),
                ));
            }

            // Assume PAR
            if self.config.in_par_assume == InParAssume::None {
                return Err(Error::MultipleChromAlignsLikelyPar(
                    tx_ac.to_string(),
                    self.config.assembly.clone(),
                    self.config.alt_aln_method.to_string(),
                ));
            }

            let alt_acs: Vec<String> = if self.config.in_par_assume == InParAssume::None {
                alt_acs.into_iter().collect()
            } else {
                alt_acs
                    .into_iter()
                    .filter(|ac| {
                        // Resolve the raw accession (e.g., "NC_000023.11") to its chromosome name (e.g., "X")
                        if let Some(chrom_name) = self.asm_map.get(ac) {
                            chrom_name == self.config.in_par_assume.chrom_name()
                        } else {
                            // Consistently filter out accessions that are missing from the asm_map
                            false
                        }
                    })
                    .collect()
            };
            if alt_acs.len() != 1 {
                Err(Error::MultipleChromAlignsInParAssume(
                    tx_ac.to_string(),
                    self.config.assembly.clone(),
                    self.config.alt_aln_method.to_string(),
                    format!("{:?}", self.config.in_par_assume),
                    alt_acs.len(),
                ))
            } else {
                Ok(alt_acs
                    .first()
                    .expect("should not happen; checked for exactly one alt_ac above")
                    .to_owned())
            }
        }
    }

    pub fn replace_reference(&self, var: HgvsVariant) -> Result<HgvsVariant, Error> {
        self.inner.replace_reference(var)
    }
}

#[cfg(test)]
mod test {
    use crate::data::uta_sr::test_helpers::build_provider;
    use anyhow::Error;

    use super::{Config, Mapper};

    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::Mapper>();
    }

    fn build_mapper_38(normalize: bool) -> Result<Mapper, Error> {
        let provider = build_provider()?;
        let config = Config {
            assembly: "GRCh38".to_string(),
            normalize,
            renormalize_g: false,
            ..Config::default()
        };
        Ok(Mapper::new(config, provider))
    }

    fn build_mapper_37(normalize: bool) -> Result<Mapper, Error> {
        let provider = build_provider()?;
        let config = Config {
            assembly: "GRCh37".to_string(),
            normalize,
            renormalize_g: false,
            ..Config::default()
        };
        Ok(Mapper::new(config, provider))
    }

    #[test]
    fn smoke() -> Result<(), Error> {
        build_mapper_38(true)?;
        Ok(())
    }

    /// Corresponds to `tests/issues/test_655.py` of biocommons/hgvs.
    mod issue_655 {
        use anyhow::Error;
        use std::str::FromStr;

        use crate::data::uta_sr::test_helpers::build_provider;
        use crate::parser::HgvsVariant;

        use super::{Config, Mapper};

        fn build_mapper(ins_at_boundary_is_intronic: bool) -> Result<Mapper, Error> {
            let config = Config {
                assembly: "GRCh37".to_string(),
                renormalize_g: false,
                ins_at_boundary_is_intronic,
                ..Config::default()
            };
            Ok(Mapper::new(config, build_provider()?))
        }

        /// DNAH11, plus strand, exon 21 ends at c.4110.  The dup copies the genome base G at
        /// c.4110+1; the transcript base T at c.4110 would give `p.Glu1371Ter`.
        #[test]
        fn c_to_p_dup_of_intron_base() -> Result<(), Error> {
            let var_c = HgvsVariant::from_str("NM_003777.3:c.4110+1dup")?;

            let exonic = build_mapper(false)?;
            assert_eq!(
                format!("{}", exonic.c_to_p(&var_c)?),
                "NP_003768.2:p.Glu1371GlyfsTer23"
            );
            // Without the genome, the variant mapper lacks the reference of the dup.
            assert!(matches!(
                exonic.variant_mapper().c_to_p(&var_c, None),
                Err(crate::mapper::Error::DuplicationReferenceMissing(_))
            ));

            let intronic = build_mapper(true)?;
            assert_eq!(format!("{}", intronic.c_to_p(&var_c)?), "NP_003768.2:p.?");

            Ok(())
        }
    }

    /// Corresponds to `tests/issues/test_714.py` of biocommons/hgvs.  Its transcripts are not
    /// in the test data, so these cases use KRAS and PTEN instead.
    mod issue_714 {
        use anyhow::Error;
        use std::collections::HashMap;
        use std::str::FromStr;

        use crate::data::uta_sr::test_helpers::build_provider;
        use crate::mapper::variant::ShiftOverBoundaryPreference;
        use crate::parser::HgvsVariant;

        use super::{Config, Mapper};

        fn config(
            shift_over_boundary: bool,
            shift_over_boundary_preference: ShiftOverBoundaryPreference,
        ) -> Config {
            Config {
                assembly: "GRCh37".to_string(),
                renormalize_g: false,
                shift_over_boundary,
                shift_over_boundary_preference,
                ..Config::default()
            }
        }

        fn build_mapper(
            shift_over_boundary: bool,
            shift_over_boundary_preference: ShiftOverBoundaryPreference,
        ) -> Result<Mapper, Error> {
            Ok(Mapper::new(
                config(shift_over_boundary, shift_over_boundary_preference),
                build_provider()?,
            ))
        }

        /// `(var_c, var_p without shifts or with DEFAULT, var_p with EXON, var_p with INTRON,
        /// var_p with KEEP_SPLICE_SITE)`
        ///
        /// KRAS is on the minus strand.  Exon 2 ends with `c.106_111` `ATAGAG` (Ile36, Glu37)
        /// and intron 2 starts with `c.111+1_111+3` `GTA`.  Exon 3 starts with `c.112_123`
        /// `GATTCCTACAGG` (Asp38 to Arg41) after `c.112-2_112-1` `AG`.  PTEN is on the plus
        /// strand.  Exon 3 ends with `c.207_209` `TCT` and intron 3 starts with `c.209+1_209+2`
        /// `GT`.  Exon 8 starts with `c.802_810` `GACAAAATG` (Asp268 to Met270) after
        /// `c.802-17_802-1` `TTTTTTTTTTTTTTTAG`.
        ///
        /// KEEP_SPLICE_SITE gives the result of EXON unless a placement in the intron keeps the
        /// splice site.
        pub(super) const CASES: &[(&str, &str, &str, &str, &str)] = &[
            // An intron insertion, which fits into the exon as `c.110_111insGGC`.
            (
                "NM_004985.4:c.111+1_111+2insCGG",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Glu37_Asp38insAla",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Glu37_Asp38insAla",
            ),
            // The same insertion, given in the exon.
            (
                "NM_004985.4:c.110_111insGGC",
                "NP_004976.2:p.Glu37_Asp38insAla",
                "NP_004976.2:p.Glu37_Asp38insAla",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Glu37_Asp38insAla",
            ),
            // A dup over the exon end, which fits into the exon as `c.108_109insGAG`.
            (
                "NM_004985.4:c.110_111+1dup",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Glu37dup",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Glu37dup",
            ),
            // A dup over the acceptor site, which fits into the exon as `c.112_113insGG`.
            (
                "NM_004985.4:c.112-1_112dup",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer8",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer8",
            ),
            // A dup of the last exon base, which fits into the exon as `c.208_209insT`.
            (
                "NM_000314.4:c.209dup",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Cys71LeufsTer3",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Cys71LeufsTer3",
            ),
            // A dup over the acceptor site, which fits into the exon as `c.802_803insGG`.
            (
                "NM_000314.4:c.802-1_802dup",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Asp268GlyfsTer9",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Asp268GlyfsTer9",
            ),
            // A dup in the poly-T before the acceptor, which cannot reach the exon.
            (
                "NM_000314.4:c.802-3dup",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
            ),
            // A dup of FOLR3 exon 3, whose alignment has a gap (`150=2D37=`), and of the 20 intron
            // bases before it.
            (
                "NM_000804.2:c.169-20_357dup",
                "NP_000795.2:p.?",
                "NP_000795.2:p.?",
                "NP_000795.2:p.?",
                "NP_000795.2:p.?",
            ),
            // Pairs of the same situation on the minus strand (KRAS) and on the plus strand
            // (PTEN).  A shift rewrites a dup as an insertion at the end of the dup that the
            // shift reached, so both strands give the same result.  Each comment names the
            // KRAS variant first.
            //
            // An intron dup at the donor site.  It shifts to `c.111dup` and `c.209_209+1dup`,
            // and fits into the exon as `c.110_111insG` and `c.208_209insTG`.  The PTEN dup also
            // fits into the intron after the splice site as `c.209+2_209+3insGT`, so
            // KEEP_SPLICE_SITE gives `p.?`.
            (
                "NM_004985.4:c.111+1dup",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer10",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer10",
            ),
            (
                "NM_000314.4:c.209+1_209+2dup",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Cys71PhefsTer29",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
            ),
            // An intron dup at the acceptor site.  It shifts to `c.112dup` and `c.802dup`, and
            // fits into the exon as `c.112_113insG` and `c.802_803insG`.
            (
                "NM_004985.4:c.112-1dup",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer10",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer10",
            ),
            (
                "NM_000314.4:c.802-1dup",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Asp268GlyfsTer30",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Asp268GlyfsTer30",
            ),
            // An exon dup before the donor site.  It shifts to `c.110_111dup` and
            // `c.208_209dup`, and fits into the intron as `c.111_111+1insAG` and
            // `c.209_209+1insCT`.  These insertions at the exon edge count as intronic.
            (
                "NM_004985.4:c.109_110dup",
                "NP_004976.2:p.Asp38ArgfsTer8",
                "NP_004976.2:p.Asp38ArgfsTer8",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38ArgfsTer8",
            ),
            (
                "NM_000314.4:c.207_208dup",
                "NP_000305.3:p.Cys71PhefsTer29",
                "NP_000305.3:p.Cys71PhefsTer29",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Cys71PhefsTer29",
            ),
            // An exon dup after the acceptor site.  It shifts to `c.112_121dup` and
            // `c.802_809dup`, and fits into the intron as `c.112-1_112insGATTCCTACA` and
            // `c.802-1_802insGACAAAAT`.  These insertions at the exon edge count as intronic.
            (
                "NM_004985.4:c.113_122dup",
                "NP_004976.2:p.Lys42PhefsTer9",
                "NP_004976.2:p.Lys42PhefsTer9",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Lys42PhefsTer9",
            ),
            (
                "NM_000314.4:c.803_810dup",
                "NP_000305.3:p.Phe271ThrfsTer8",
                "NP_000305.3:p.Phe271ThrfsTer8",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Phe271ThrfsTer8",
            ),
        ];

        /// As `CASES`, with `ins_at_boundary_is_intronic` set to `false`, so that an insertion
        /// at an exon edge counts as exonic.
        pub(super) const CASES_INS_AT_BOUNDARY_IS_EXONIC: &[(&str, &str, &str, &str, &str)] = &[
            // A dup of the first intron base, which counts as exonic.  It fits into the intron as
            // `c.111+1_111+2insG`.
            (
                "NM_004985.4:c.111+1dup",
                "NP_004976.2:p.Asp38GlyfsTer10",
                "NP_004976.2:p.Asp38GlyfsTer10",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer10",
            ),
            // A dup of the last exon base.  Both of its rewrites, `c.208_209insT` and
            // `c.209_209+1insT`, count as exonic.
            (
                "NM_000314.4:c.209dup",
                "NP_000305.3:p.Cys71LeufsTer3",
                "NP_000305.3:p.Cys71LeufsTer3",
                "NP_000305.3:p.Cys71LeufsTer3",
                "NP_000305.3:p.Cys71LeufsTer3",
            ),
            // A dup of the first exon base on both strands.  It shifts to `c.112-1dup` and
            // `c.802-1dup`, and fits into the intron as `c.112-2_112-1insG` and
            // `c.802-2_802-1insG`.
            (
                "NM_004985.4:c.112dup",
                "NP_004976.2:p.Asp38GlyfsTer10",
                "NP_004976.2:p.Asp38GlyfsTer10",
                "NP_004976.2:p.?",
                "NP_004976.2:p.Asp38GlyfsTer10",
            ),
            (
                "NM_000314.4:c.802dup",
                "NP_000305.3:p.Asp268GlyfsTer30",
                "NP_000305.3:p.Asp268GlyfsTer30",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Asp268GlyfsTer30",
            ),
        ];

        #[test]
        fn c_to_p() -> Result<(), Error> {
            check_c_to_p(CASES, true)
        }

        #[test]
        fn c_to_p_ins_at_boundary_is_exonic() -> Result<(), Error> {
            check_c_to_p(CASES_INS_AT_BOUNDARY_IS_EXONIC, false)
        }

        /// Check `cases` with the option off and with each preference.
        pub(super) fn check_c_to_p(
            cases: &[(&str, &str, &str, &str, &str)],
            ins_at_boundary_is_intronic: bool,
        ) -> Result<(), Error> {
            let build = |shift_over_boundary, shift_over_boundary_preference| {
                let config = Config {
                    ins_at_boundary_is_intronic,
                    ..config(shift_over_boundary, shift_over_boundary_preference)
                };
                Ok::<_, Error>(Mapper::new(config, build_provider()?))
            };
            let disabled = build(false, ShiftOverBoundaryPreference::Exon)?;
            let default = build(true, ShiftOverBoundaryPreference::Default)?;
            let exon = build(true, ShiftOverBoundaryPreference::Exon)?;
            let intron = build(true, ShiftOverBoundaryPreference::Intron)?;
            let keep_splice_site = build(true, ShiftOverBoundaryPreference::KeepSpliceSite)?;
            for (hgvs_c, expected_default, expected_exon, expected_intron, expected_keep) in cases {
                let var_c = HgvsVariant::from_str(hgvs_c)?;
                for (label, mapper, expected) in [
                    ("disabled", &disabled, expected_default),
                    ("DEFAULT", &default, expected_default),
                    ("EXON", &exon, expected_exon),
                    ("INTRON", &intron, expected_intron),
                    ("KEEP_SPLICE_SITE", &keep_splice_site, expected_keep),
                ] {
                    assert_eq!(
                        format!("{}", mapper.c_to_p(&var_c)?),
                        *expected,
                        "{hgvs_c} with {label}"
                    );
                }
            }
            Ok(())
        }

        /// Only the shifts look up the genome alignment.  Without it, `Default` still works.
        #[test]
        fn c_to_p_without_alignment() -> Result<(), Error> {
            // An assembly without chromosome 12, so KRAS has no alignment.
            let asm_map = HashMap::from([("NC_000010.10".to_string(), "10".to_string())]);
            let build = |shift_over_boundary_preference| {
                let config = Config {
                    asm_map: Some(asm_map.clone()),
                    ..config(true, shift_over_boundary_preference)
                };
                Ok::<_, Error>(Mapper::new(config, build_provider()?))
            };
            let var_c = HgvsVariant::from_str("NM_004985.4:c.110_111insGGC")?;
            assert_eq!(
                format!(
                    "{}",
                    build(ShiftOverBoundaryPreference::Default)?.c_to_p(&var_c)?
                ),
                "NP_004976.2:p.Glu37_Asp38insAla"
            );
            assert!(matches!(
                build(ShiftOverBoundaryPreference::Exon)?.c_to_p(&var_c),
                Err(crate::mapper::Error::NoAlignments(..))
            ));
            Ok(())
        }

        /// The variant mapper lacks the genome alignment for the shifts.
        #[test]
        fn variant_mapper_c_to_p_needs_alignment() -> Result<(), Error> {
            let mapper = build_mapper(true, ShiftOverBoundaryPreference::Exon)?;
            let var_c = HgvsVariant::from_str("NM_004985.4:c.110_111insGGC")?;
            assert!(matches!(
                mapper.variant_mapper().c_to_p(&var_c, None),
                Err(crate::mapper::Error::ShiftOverBoundaryNeedsAlignment(_))
            ));
            Ok(())
        }
    }

    /// `ShiftOverBoundaryPreference::KeepSpliceSite`, which biocommons/hgvs does not have.
    mod keep_splice_site {
        use anyhow::Error;
        use std::str::FromStr;

        use crate::data::uta_sr::test_helpers::build_provider;
        use crate::mapper::variant::ShiftOverBoundaryPreference;
        use crate::parser::HgvsVariant;

        use super::issue_714::{check_c_to_p, CASES, CASES_INS_AT_BOUNDARY_IS_EXONIC};
        use super::{Config, Mapper};

        /// As `issue_714::CASES`, on the same transcripts.  PTEN intron 3 starts with
        /// `c.209+1_209+9` `GTAAGTATG`.  PTEN exon 2 ends with `c.160_164` `GTAAG` and intron 2
        /// starts with `c.164+1_164+5` `GTAAG`.
        const CASES_KEEP: &[(&str, &str, &str, &str, &str)] = &[
            // A dup of the first intron base.  It fits into the exon only at its edge, as
            // `c.209_209+1insG`.  This insertion keeps the splice site.
            (
                "NM_000314.4:c.209+1dup",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Cys71LeufsTer3",
            ),
            // Insertions between `+1` and `+2`.  They cannot move, so they change the splice
            // site.
            (
                "NM_004985.4:c.111+1_111+2insC",
                "NP_004976.2:p.?",
                "NP_004976.2:p.?",
                "NP_004976.2:p.?",
                "NP_004976.2:p.?",
            ),
            (
                "NM_000314.4:c.209+1_209+2insA",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
            ),
            // A dup of the intron start `GTAAG`, which repeats the exon end.  It fits into the
            // exon as a dup of `c.160_164`, and into the intron after the splice site as
            // `c.164+5_164+6insGTAAG`.
            (
                "NM_000314.4:c.164+1_164+5dup",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Phe56Ter",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
            ),
        ];

        /// As `CASES_KEEP`, with `ins_at_boundary_is_intronic` set to `false`.
        const CASES_KEEP_INS_AT_BOUNDARY_IS_EXONIC: &[(&str, &str, &str, &str, &str)] = &[
            // Dups of the essential splice-site bases.  They fit at the exon edge, which counts
            // as exonic here, and into the intron next to the essential bases.
            (
                "NM_004985.4:c.111+1_111+2dup",
                "NP_004976.2:p.Asp38ValfsTer8",
                "NP_004976.2:p.Asp38ValfsTer8",
                "NP_004976.2:p.?",
                "NP_004976.2:p.?",
            ),
            (
                "NM_004985.4:c.112-2_112-1dup",
                "NP_004976.2:p.Asp38ArgfsTer8",
                "NP_004976.2:p.Asp38ArgfsTer8",
                "NP_004976.2:p.?",
                "NP_004976.2:p.?",
            ),
            (
                "NM_000314.4:c.802-2_802-1dup",
                "NP_000305.3:p.Asp268ArgfsTer9",
                "NP_000305.3:p.Asp268ArgfsTer9",
                "NP_000305.3:p.?",
                "NP_000305.3:p.?",
            ),
            // As in `CASES_KEEP`.  The other preferences now count the dup as exonic.
            (
                "NM_000314.4:c.209+1dup",
                "NP_000305.3:p.Cys71LeufsTer3",
                "NP_000305.3:p.Cys71LeufsTer3",
                "NP_000305.3:p.?",
                "NP_000305.3:p.Cys71LeufsTer3",
            ),
        ];

        #[test]
        fn c_to_p() -> Result<(), Error> {
            check_c_to_p(CASES_KEEP, true)
        }

        #[test]
        fn c_to_p_ins_at_boundary_is_exonic() -> Result<(), Error> {
            check_c_to_p(CASES_KEEP_INS_AT_BOUNDARY_IS_EXONIC, false)
        }

        /// `KeepSpliceSite` decides by sequence, so `ins_at_boundary_is_intronic` does not
        /// change its result.
        #[test]
        fn c_to_p_ignores_ins_at_boundary_is_intronic() -> Result<(), Error> {
            let build = |ins_at_boundary_is_intronic| {
                let config = Config {
                    assembly: "GRCh37".to_string(),
                    renormalize_g: false,
                    ins_at_boundary_is_intronic,
                    shift_over_boundary: true,
                    shift_over_boundary_preference: ShiftOverBoundaryPreference::KeepSpliceSite,
                    ..Config::default()
                };
                Ok::<_, Error>(Mapper::new(config, build_provider()?))
            };
            let intronic = build(true)?;
            let exonic = build(false)?;
            for cases in [
                CASES,
                CASES_INS_AT_BOUNDARY_IS_EXONIC,
                CASES_KEEP,
                CASES_KEEP_INS_AT_BOUNDARY_IS_EXONIC,
            ] {
                for (hgvs_c, _, _, _, expected) in cases {
                    let var_c = HgvsVariant::from_str(hgvs_c)?;
                    for mapper in [&intronic, &exonic] {
                        assert_eq!(format!("{}", mapper.c_to_p(&var_c)?), *expected, "{hgvs_c}");
                    }
                }
            }
            Ok(())
        }
    }

    /// The following is a port of `Test_variant::Mapper` in
    /// `test_hgvs_variant::Mapper_near_discrepancies.py` (sic!)
    mod cases {
        use anyhow::Error;
        use std::str::FromStr;

        use rstest::rstest;

        use crate::parser::{HgvsVariant, NoRef};

        use super::{build_mapper_37, build_mapper_38};

        #[test]
        fn test_quick_aoah() -> Result<(), Error> {
            let mapper = build_mapper_38(true)?;
            let hgvs_g = "NC_000007.13:g.36561662C>T";
            let hgvs_c = "NM_001637.3:c.1582G>A";
            let hgvs_p = "NP_001628.1:p.Gly528Arg";

            let var_g = HgvsVariant::from_str(hgvs_g)?;
            let var_c = mapper.g_to_c(&var_g, "NM_001637.3")?;
            let var_p = mapper.c_to_p(&var_c)?;

            assert_eq!(format!("{var_c}"), hgvs_c);
            assert_eq!(format!("{var_p}"), hgvs_p);

            Ok(())
        }

        #[test]
        fn test_c_to_p_brca2() -> Result<(), Error> {
            let mapper = build_mapper_38(true)?;
            let hgvs_c = "NM_000059.3:c.7790delAAG";
            let var_c = HgvsVariant::from_str(hgvs_c)?;

            // NB: this used to fail in Python hgvs but works here now as we do not
            // perform comprehensive validation yet (1 bp interval/position, but 3bp
            // deleted).
            assert_eq!(
                "NP_000050.2:p.Glu2598LysfsTer50",
                format!("{}", mapper.c_to_p(&var_c)?)
            );

            Ok(())
        }

        #[rstest]
        #[case("NM_000059.3:c.7791A>G", "NP_000050.2:p.Lys2597=", "BRCA2", 38)]
        #[case("NM_000302.3:c.1594_1596del", "NP_000293.2:p.Glu532del", "PLOD1", 38)]
        #[case(
            "NM_000090.3:c.2490_2516del",
            "NP_000081.1:p.Glu832_Gly840del",
            "COL3A1",
            38
        )]
        #[case("NM_001292004.1:c.376=", "NC_000005.10:g.74715659=", "HEXB", 38)]
        #[case(
            "NM_000116.4:c.-8_-3inv6",
            "NC_000023.11:g.154411836_154411841inv",
            "TAZ",
            38
        )]
        #[case(
            "NM_000348.3:c.89_91inv3",
            "NC_000002.12:g.31580810_31580812inv",
            "SDR5A2",
            38
        )]
        #[case("NM_001637.3:c.1582_1583inv", "NP_001628.1:p.Gly528Pro", "AOAH", 38)]
        #[case("NM_025137.3:c.-20_*20inv", "NP_079413.3:p.?", "SPG11", 38)]
        fn project_c_to_x(
            #[case] hgvs_c: &str,
            #[case] hgvs_x: &str,
            #[case] gene: &str,
            #[case] build: i32,
        ) -> Result<(), Error> {
            let mapper = if build == 38 {
                build_mapper_38(true)?
            } else {
                build_mapper_37(true)?
            };

            let var_c = HgvsVariant::from_str(hgvs_c)?;
            let var_x = HgvsVariant::from_str(hgvs_x)?;
            let actual = match &var_x {
                HgvsVariant::GenomeVariant { .. } => mapper.c_to_g(&var_c)?,
                HgvsVariant::ProtVariant { .. } => mapper.c_to_p(&var_c)?,
                _ => panic!("not implemented"),
            };

            assert_eq!(format!("{}", NoRef(&actual)), hgvs_x, "gene={gene}");

            Ok(())
        }

        #[rstest]
        #[case(
            "NC_000019.10:g.50378563_50378564insTG",
            "NM_007121.5:n.796_798delinsTG",
            "NR1H2",
            38
        )]
        #[case(
            "NC_000007.14:g.149779575delC",
            "NM_198455.2:n.1115_1116insAG",
            "SSPO",
            38
        )]
        #[case(
            "NC_000012.11:g.122064775C>T",
            "NM_032790.3:c.127_128insTGCCAC",
            "ORAI1",
            38
        )]
        #[case(
            "NC_000002.11:g.73675227_73675228insCTC",
            "NM_015120.4:c.1574_1576=",
            "ALMS1",
            37
        )]
        #[case("NM_000116.4:c.-120_-119insT", "NC_000023.10:g.153640061=", "TAZ", 37)]
        #[case(
            "NM_000348.3:c.88del",
            "NC_000002.11:g.(31805882_31805883)=",
            "SRD5A2",
            37
        )]
        fn project_at_alignment_discrepancies(
            #[case] hgvs_lhs: &str,
            #[case] hgvs_rhs: &str,
            #[case] gene: &str,
            #[case] build: i32,
        ) -> Result<(), Error> {
            let mapper = if build == 38 {
                build_mapper_38(true)?
            } else {
                build_mapper_37(true)?
            };

            let var_lhs = HgvsVariant::from_str(hgvs_lhs)?;
            let var_rhs = HgvsVariant::from_str(hgvs_rhs)?;
            let actual = match (&var_lhs, &var_rhs) {
                (HgvsVariant::CdsVariant { .. }, HgvsVariant::GenomeVariant { .. }) => {
                    mapper.c_to_g(&var_lhs)?
                }
                (HgvsVariant::GenomeVariant { .. }, HgvsVariant::CdsVariant { .. }) => {
                    mapper.g_to_c(&var_lhs, var_rhs.accession())?
                }
                (HgvsVariant::GenomeVariant { .. }, HgvsVariant::TxVariant { .. }) => {
                    mapper.g_to_n(&var_lhs, var_rhs.accession())?
                }
                _ => panic!("not implemented"),
            };

            assert_eq!(
                format!("{}", crate::parser::NoRef(&actual)),
                hgvs_rhs,
                "gene={gene}"
            );

            Ok(())
        }

        #[rstest]
        #[case(
            "NM_080877.2:c.1733_1735delinsTTT",
            "NP_543153.1:p.Pro578_Lys579delinsLeuTer",
            "SLC34A3"
        )]
        #[case(
            "NM_001034853.1:c.2847_2848delAGinsCT",
            "NP_001030025.1:p.Glu949_Glu950delinsAspTer",
            "RPGR"
        )]
        #[case(
            "NM_001034853.1:c.2847_2848inv",
            "NP_001030025.1:p.Glu949_Glu950delinsAspTer",
            "RPGR"
        )]
        #[case("NM_080877.2:c.1735A>T", "NP_543153.1:p.Lys579Ter", "SLC34A3")]
        #[case("NM_080877.2:c.1795_*3delinsTAG", "NP_543153.1:p.Leu599Ter", "SLC34A3")]
        fn c_to_p_with_stop_gain(
            #[case] hgvs_c: &str,
            #[case] hgvs_p: &str,
            #[case] gene: &str,
        ) -> Result<(), Error> {
            let mapper = build_mapper_38(true)?;
            let var_c = HgvsVariant::from_str(hgvs_c)?;

            let actual = mapper.c_to_p(&var_c)?;

            assert_eq!(
                format!("{}", crate::parser::NoRef(&actual)),
                hgvs_p,
                "gene={gene}, hgvs_c={hgvs_c}, hgvs_p={hgvs_p}"
            );

            Ok(())
        }
    }

    /// The following is a port of the `Test_RefReplacement` in
    /// `test_hgvs_variant::Mapper_near_discrepancies.py` (sic!)
    mod ref_replacement {
        use anyhow::Error;
        use std::str::FromStr;

        use rstest::rstest;

        use crate::parser::HgvsVariant;

        use super::build_mapper_38;

        // These casese attempt to test reference update in four dimensions:
        //
        // - variant type: n, c, g
        // - major mapping paths: c<->n, c<->g, n<->g
        // - variant class: sub, del, ins, delins, dup
        // - strand: +/-
        //
        // ADRB2    │ NM_000024.5 │  239 │ 1481 │ NC_000005.9  │  1 │ 148206155,148208197 | 284=1X32=1X1724=
        //
        // cseq = hdp.fetch_seq("NM_000024.5")
        // gseq = hdp.fetch_seq("NC_000005.9",148206155,148208197)
        // cseq[280:290] = "CAATAGAAGC"
        // gseq[280:290] = "CAATGGAAGC"
        //                      ^ @ n.285
        #[rstest]
        // These variants are in and around the first sub:
        #[case(
            "NM_000024.5:c.42C>N",
            "NC_000005.9:g.148206436C>N",
            "NM_000024.5:n.281C>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.43A>N",
            "NC_000005.9:g.148206437A>N",
            "NM_000024.5:n.282A>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.44A>N",
            "NC_000005.9:g.148206438A>N",
            "NM_000024.5:n.283A>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.45T>N",
            "NC_000005.9:g.148206439T>N",
            "NM_000024.5:n.284T>N",
            "ADRB2"
        )]
        // ref repl
        #[case(
            "NM_000024.5:c.46A>N",
            "NC_000005.9:g.148206440G>N",
            "NM_000024.5:n.285A>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.47G>N",
            "NC_000005.9:g.148206441G>N",
            "NM_000024.5:n.286G>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.48A>N",
            "NC_000005.9:g.148206442A>N",
            "NM_000024.5:n.287A>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.49A>N",
            "NC_000005.9:g.148206443A>N",
            "NM_000024.5:n.288A>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.50G>N",
            "NC_000005.9:g.148206444G>N",
            "NM_000024.5:n.289G>N",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.51C>N",
            "NC_000005.9:g.148206445C>N",
            "NM_000024.5:n.281C>N",
            "ADRB2"
        )]
        // ins, del, delins, dup
        #[case(
            "NM_000024.5:c.46_47insNN",
            "NC_000005.9:g.148206440_148206441insNN",
            "NM_000024.5:n.285_286insNN",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.45_47delTAG",
            "NC_000005.9:g.148206439_148206441delTGG",
            "NM_000024.5:n.284_286delTAG",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.45_47delTAGinsNNNN",
            "NC_000005.9:g.148206439_148206441delTGGinsNNNN",
            "NM_000024.5:n.284_286delTAGinsNNNN",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.45_47delTAGinsNNNN",
            "NC_000005.9:g.148206439_148206441delTGGinsNNNN",
            "NM_000024.5:n.284_286delTAGinsNNNN",
            "ADRB2"
        )]
        #[case(
            "NM_000024.5:c.46dupA",
            "NC_000005.9:g.148206440dupG",
            "NM_000024.5:n.285dupA",
            "ADRB2"
        )]
        // IFNA16   │ NM_002173.2 │    6 │  576 │ NC_000009.11 │ -1 │  21216371, 21217310 | 691=2X246=
        // cseq = hdp.fetch_seq("NM_002173.2")
        // gseq = reverse_complement(hdp.fetch_seq("NC_000009.11",21216371,21217310))
        // cseq[685:695] = "AAATTTCAAA"
        // gseq[685:695] = "AAATTTTCAA"
        //                        ^^ @ n.692_693
        // These variants are in and around the 2X substitution
        #[case(
            "NM_002173.2:c.*110A>N",
            "NC_000009.11:g.21216625T>N",
            "NM_002173.2:n.686A>N",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*111A>N",
            "NC_000009.11:g.21216624T>N",
            "NM_002173.2:n.687A>N",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*112A>N",
            "NC_000009.11:g.21216623T>N",
            "NM_002173.2:n.688A>N",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*113T>N",
            "NC_000009.11:g.21216622A>N",
            "NM_002173.2:n.689T>N",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*114T>N",
            "NC_000009.11:g.21216621A>N",
            "NM_002173.2:n.690T>N",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*115T>N",
            "NC_000009.11:g.21216620A>N",
            "NM_002173.2:n.691T>N",
            "IFNA16"
        )]
        // ref repl
        #[case(
            "NM_002173.2:c.*116C>N",
            "NC_000009.11:g.21216619A>N",
            "NM_002173.2:n.692C>N",
            "IFNA16"
        )]
        // ref repl
        #[case(
            "NM_002173.2:c.*117A>N",
            "NC_000009.11:g.21216618G>N",
            "NM_002173.2:n.693A>N",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*118A>N",
            "NC_000009.11:g.21216617T>N",
            "NM_002173.2:n.694A>N",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*119A>N",
            "NC_000009.11:g.21216616T>N",
            "NM_002173.2:n.695A>N",
            "IFNA16"
        )]
        // ins, del, delins, dup:
        #[case(
            "NM_002173.2:c.*115_*117insNN",
            "NC_000009.11:g.21216618_21216620insNN",
            "NM_002173.2:n.691_693insNN",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*114_*117delTTCA",
            "NC_000009.11:g.21216618_21216621delGAAA",
            "NM_002173.2:n.690_693delTTCA",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*115_*117delTCAinsNN",
            "NC_000009.11:g.21216618_21216620delGAAinsNN",
            "NM_002173.2:n.691_693delTCAinsNN",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*115_*117delTCAinsNN",
            "NC_000009.11:g.21216618_21216620delGAAinsNN",
            "NM_002173.2:n.691_693delTCAinsNN",
            "IFNA16"
        )]
        #[case(
            "NM_002173.2:c.*115_*117dupTCA",
            "NC_000009.11:g.21216618_21216620dupGAA",
            "NM_002173.2:n.691_693dupTCA",
            "IFNA16"
        )]
        // genomic variation within an indel discrepancy
        // NM_032790.3  c        108 >
        // NM_032790.3  n        301 > CGGGGAGCCCCCGGGGGCC------CCGCCACCGCCGCCGT
        //                             |||||||||||||||||||------||||||||||||||||
        // NC_000012.11 g  122064755 > CGGGGAGCCCCCGGGGGCCCCGCCACCGCCACCGCCGCCGT
        //                                              **********
        #[case(
            "NM_032790.3:c.125_128delCCCC",
            "NC_000012.11:g.122064772_122064781delCCCCGCCACC",
            "NM_032790.3:n.318_321delCCCC",
            "ORAI1"
        )]
        fn cases(
            #[case] hgvs_c: &str,
            #[case] hgvs_g: &str,
            #[case] hgvs_n: &str,
            #[case] gene: &str,
        ) -> Result<(), Error> {
            let mapper = build_mapper_38(true)?;
            let var_c = HgvsVariant::from_str(hgvs_c)?;
            let var_g = HgvsVariant::from_str(hgvs_g)?;
            let var_n = HgvsVariant::from_str(hgvs_n)?;

            let var_c = var_c.with_reference("NNNNNN".to_string());
            let fixed_c = mapper.replace_reference(var_c)?;
            assert_eq!(format!("{}", fixed_c), hgvs_c, "gene={gene}");

            let var_g = var_g.with_reference("NNNNNN".to_string());
            let fixed_g = mapper.replace_reference(var_g)?;
            assert_eq!(format!("{}", fixed_g), hgvs_g, "gene={gene}");

            let var_n = var_n.with_reference("NNNNNN".to_string());
            let fixed_n = mapper.replace_reference(var_n)?;
            assert_eq!(format!("{}", fixed_n), hgvs_n, "gene={gene}");

            Ok(())
        }
    }

    /// The following is a port of `Test_AssemblyMapper` in
    /// `test_hgvs_variant::Mapper_near_discrepancies.py` (sic!)
    mod projections {
        use anyhow::Error;
        use std::str::FromStr;

        use crate::parser::HgvsVariant;

        use super::build_mapper_37;

        fn test_projections(
            hgvs_g: &str,
            hgvs_c: &str,
            hgvs_n: &str,
            hgvs_p: &str,
        ) -> Result<(), Error> {
            let mapper = build_mapper_37(false)?;

            let var_g = HgvsVariant::from_str(hgvs_g)?;
            let var_c = HgvsVariant::from_str(hgvs_c)?;
            let var_n = HgvsVariant::from_str(hgvs_n)?;

            let res_cg = mapper.c_to_g(&var_c)?;
            assert_eq!(format!("{res_cg}"), hgvs_g,);

            let res_gc = mapper.g_to_c(&var_g, var_c.accession())?;
            assert_eq!(format!("{res_gc}"), hgvs_c,);

            let res_ng = mapper.n_to_g(&var_n)?;
            assert_eq!(format!("{res_ng}"), hgvs_g,);

            let res_gn = mapper.g_to_n(&var_g, var_n.accession())?;
            assert_eq!(format!("{res_gn}"), hgvs_n,);

            let res_cn = mapper.c_to_n(&var_c)?;
            assert_eq!(format!("{res_cn}"), hgvs_n,);

            let res_nc = mapper.n_to_c(&var_n)?;
            assert_eq!(format!("{res_nc}"), hgvs_c,);

            let res_cp = mapper.c_to_p(&var_c)?;
            assert_eq!(format!("{res_cp}"), hgvs_p,);

            Ok(())
        }

        #[test]
        fn snv() -> Result<(), Error> {
            test_projections(
                "NC_000007.13:g.36561662C>T",
                "NM_001637.3:c.1582G>A",
                "NM_001637.3:n.1983G>A",
                "NP_001628.1:p.Gly528Arg",
            )
        }

        #[test]
        fn intronic() -> Result<(), Error> {
            test_projections(
                "NC_000010.10:g.89711873A>C",
                "NM_000314.4:c.493-2A>C",
                "NM_000314.4:n.1524-2A>C",
                "NP_000305.3:p.?",
            )
        }
    }

    /// The following is a port of the `test_hgvs_variant::Mapper_near_discrepancies.py` (sic!)
    mod near_discrepancies {
        use anyhow::Error;
        use std::str::FromStr;

        use crate::parser::{HgvsVariant, NoRef};

        use super::build_mapper_38;

        #[test]
        fn run() -> Result<(), Error> {
            let mapper = build_mapper_38(false)?;
            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_path("tests/data/mapper/proj-near-disc.tsv")?;

            for row in rdr.records() {
                let row = row?;
                let disc_type = row.get(0).expect("problem in test TSV file");
                let loc_type = row.get(1).expect("problem in test TSV file");
                let variant = row.get(2).expect("problem in test TSV file");
                let expected = row.get(3).expect("problem in test TSV file");

                let var_n = HgvsVariant::from_str(variant)?;
                let var_g = mapper.n_to_g(&var_n)?;

                let actual = format!("{}", NoRef(&var_g)).replace(['(', ')'], "");

                assert_eq!(
                    actual, expected,
                    "loc_type={loc_type} disc_type={disc_type}",
                )
            }

            Ok(())
        }
    }
}
