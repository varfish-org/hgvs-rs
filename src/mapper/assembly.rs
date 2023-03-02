//! Implementation of `AssemblyMapper`.

use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::rc::Rc;

use crate::mapper::variant::{Config as VariantMapperConfig, Mapper as VariantMapper};
use crate::parser::HgvsVariant;
use crate::{data::interface::Provider, static_data::Assembly, validator::ValidationLevel};

#[derive(Debug, PartialEq, Eq, Default)]
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
    pub assembly: Assembly,
    pub alt_aln_method: String,
    pub normalize: bool,
    pub in_par_assume: InParAssume,
    pub prevalidation_level: ValidationLevel,
    pub replace_reference: bool,
    pub strict_validation: bool,
    pub strict_bounds: bool,
    pub add_gene_symbol: bool,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            assembly: Assembly::Grch38,
            alt_aln_method: "splign".to_string(),
            normalize: true,
            in_par_assume: InParAssume::X,
            prevalidation_level: ValidationLevel::Intrinsic,
            replace_reference: true,
            strict_validation: true,
            strict_bounds: true,
            add_gene_symbol: false,
        }
    }
}

/// Provides simplified variant mapping for a single assembly and transcript-reference
/// alignment method.
///
/// `AssemblyMapper` uses `VariantMapper`, to provide all projection functionality, and add:
///
/// * Automatic selection of genomic sequence accession
/// * Transcript selection from genomic coordinates
/// * Normalization after projection
/// * Special handling for PAR regions
///
/// `AssemblyMapper` is instantiated with an assembly name and `alt_aln_method`. These enable
/// the following conveniences over `VariantMapper`:
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
    config: Config,
    provider: Rc<dyn Provider>,
    inner: VariantMapper,
    /// Accessions of contigs in assembly.
    asm_accessions: HashSet<String>,
    /// Map from accession to contig name.
    asm_map: HashMap<String, String>,
}

impl Mapper {
    /// Construct new assembly mapper from config and provider.
    pub fn new(config: Config, provider: Rc<dyn Provider>) -> Self {
        let inner_config = VariantMapperConfig {
            replace_reference: config.replace_reference,
            strict_validation: config.strict_validation,
            prevalidation_level: config.prevalidation_level,
            add_gene_symbol: config.add_gene_symbol,
            strict_bounds: config.strict_bounds,
        };
        let inner = VariantMapper::new(&inner_config, provider.clone());
        let asm_accessions = provider
            .as_ref()
            .get_assembly_map(config.assembly)
            .keys()
            .clone()
            .into_iter()
            .map(|s| s.to_string())
            .collect::<HashSet<_>>();
        let asm_map = HashMap::from_iter(
            provider
                .as_ref()
                .get_assembly_map(config.assembly)
                .iter()
                .map(|(key, value)| (key.clone(), value.clone())),
        );
        Self {
            config,
            provider,
            inner,
            asm_accessions,
            asm_map,
        }
    }

    /// Convert from genome (g.) variant to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    pub fn g_to_n(&self, var_g: &HgvsVariant, tx_ac: &str) -> Result<HgvsVariant, anyhow::Error> {
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
    pub fn g_to_c(&self, var_g: &HgvsVariant, tx_ac: &str) -> Result<HgvsVariant, anyhow::Error> {
        let var = self
            .inner
            .g_to_n(var_g, tx_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from genome (g.) variant to transcript variant (g. or n.).
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    pub fn g_to_t(&self, var_g: &HgvsVariant, tx_ac: &str) -> Result<HgvsVariant, anyhow::Error> {
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
    pub fn n_to_g(&self, var_n: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        let alt_ac = self.alt_ac_for_tx_ac(&var_n.accession())?;
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
    pub fn c_to_g(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        let alt_ac = self.alt_ac_for_tx_ac(&var_c.accession())?;
        let var = self
            .inner
            .n_to_g(var_c, &alt_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from transcript (c. or n.) to genome (g.) variant.
    ///
    /// # Args
    ///
    /// * `var_t` -- `HgvsVariant::TxVariant` or `HgvsVariant::CdsVariant` to project
    /// * `alt_ac` -- accession of alternativ esequence
    /// * `alt_al_method` -- alignment method, e.g., `"splign"`
    pub fn t_to_g(&self, var_t: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        let alt_ac = self.alt_ac_for_tx_ac(&var_t.accession())?;
        let var = self
            .inner
            .n_to_g(var_t, &alt_ac, &self.config.alt_aln_method)?;
        self.maybe_normalize(&var)
    }

    /// Convert from CDS variant (c.) to transcript variant (n.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::CdsVariant` to project
    pub fn c_to_n(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        let var = self.inner.c_to_n(var_c)?;
        self.maybe_normalize(&var)
    }

    /// Convert from transcript variant (n.) to CDS variant (c.).
    ///
    /// # Args
    ///
    /// * `var_n` -- `HgvsVariant::TxVariant` to project
    pub fn n_to_c(&self, var_n: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        let var = self.inner.n_to_c(var_n)?;
        self.maybe_normalize(&var)
    }

    /// Convert from CDS variant (c.) to protein variant (p.).
    ///
    /// # Args
    ///
    /// * `var_c` -- `HgvsVariant::TxVariant` to project
    pub fn c_to_p(&self, var_c: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        let var = self.inner.c_to_p(var_c, None)?;
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
    pub fn relevant_transcripts(&self, var_g: &HgvsVariant) -> Result<Vec<String>, anyhow::Error> {
        match var_g {
            HgvsVariant::GenomeVariant { loc_edit, .. } => {
                let r: Range<i32> = loc_edit.loc.inner().clone().try_into()?;
                Ok(self
                    .provider
                    .as_ref()
                    .get_tx_for_region(
                        &var_g.accession(),
                        &self.config.alt_aln_method,
                        r.start,
                        r.end,
                    )?
                    .into_iter()
                    .map(|rec| rec.tx_ac)
                    .collect::<Vec<_>>())
            }
            _ => return Err(anyhow::anyhow!("Not a GenomeVariant: {}", &var_g)),
        }
    }

    /// Normalize variant if requested and ignore errors.  This is better than checking whether
    /// the variant is intronic because future UTAs will support LRG, which will enable checking
    /// intronic variants.
    fn maybe_normalize(&self, var: &HgvsVariant) -> Result<HgvsVariant, anyhow::Error> {
        if self.config.normalize {
            self.inner
                .normalizer()?
                .normalize(var)
                .or_else(|_| Ok(var.clone()))
        } else {
            Ok(var.clone())
        }
    }

    /// Return chromosomal accession for the given transcript accession and the assembly and
    /// aln_method from configuration.
    fn alt_ac_for_tx_ac(&self, tx_ac: &str) -> Result<String, anyhow::Error> {
        let mut alt_acs = Vec::new();
        for record in self.provider.get_tx_mapping_options(tx_ac)? {
            if record.alt_aln_method == self.config.alt_aln_method
                && self.asm_accessions.contains(&record.alt_ac)
            {
                alt_acs.push(record.alt_ac);
            }
        }

        if alt_acs.is_empty() {
            return Err(anyhow::anyhow!(
                "No alignments for {} in {:?} using {}",
                tx_ac,
                self.config.assembly,
                self.config.alt_aln_method
            ));
        }

        if alt_acs.len() == 1 {
            Ok(alt_acs.first().unwrap().to_owned())
        } else {
            // alt_acs.len() > 1
            // Perform sanity check on result if more than one contig is returned.
            let mut alts = alt_acs
                .iter()
                .map(|ac| self.asm_map.get(ac))
                .filter(|x| x.is_some())
                .map(|x| x.unwrap().clone())
                .collect::<Vec<_>>();
            alts.sort();
            if alts.join("") != "XY" {
                return Err(anyhow::anyhow!(
                    "Multiple chromosomal alignments for {} in {:?} using {} \
                    (non-pseudoautosomal region) [{:?}]",
                    tx_ac,
                    self.config.assembly,
                    self.config.alt_aln_method,
                    &alts
                ));
            }

            // Assume PAR
            if self.config.in_par_assume == InParAssume::None {
                return Err(anyhow::anyhow!(
                    "Multiple chromosomal alignments for {} in {:?} using {} \
                    (likely pseudoautosomal region)",
                    tx_ac,
                    self.config.assembly,
                    self.config.alt_aln_method,
                ));
            }

            let alt_acs = alt_acs
                .into_iter()
                .filter(|ac| ac == self.config.in_par_assume.chrom_name())
                .collect::<Vec<_>>();
            if alt_acs.len() != 1 {
                Err(anyhow::anyhow!(
                    "Multiple chromosomal alignments for {} in {:?} using {} \
                    (in_par_assume={:?} selected {} of them)",
                    tx_ac,
                    self.config.assembly,
                    self.config.alt_aln_method,
                    self.config.in_par_assume,
                    alt_acs.len()
                ))
            } else {
                Ok(alt_acs.first().unwrap().to_owned())
            }
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{data::uta_sr::test_helpers::build_provider, static_data::Assembly};

    use super::{Config, Mapper};

    fn build_mapper_38() -> Result<Mapper, anyhow::Error> {
        let provider = build_provider()?;
        let config = Config {
            assembly: Assembly::Grch38,
            normalize: false,
            ..Config::default()
        };
        Ok(Mapper::new(config, provider))
    }

    #[test]
    fn smoke() -> Result<(), anyhow::Error> {
        build_mapper_38()?;
        Ok(())
    }

    /// The following is a port of the `test_hgvs_variantmapper_near_discrepancies.py` (sic!)
    mod near_discrepancies {
        use std::str::FromStr;

        use crate::parser::{HgvsVariant, NoRef};

        use super::build_mapper_38;

        #[test]
        fn run() -> Result<(), anyhow::Error> {
            let mapper = build_mapper_38()?;
            let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'#'))
                .from_path("tests/data/mapper/proj-near-disc.tsv")?;

            for row in rdr.records() {
                let row = row?;
                let disc_type = row.get(0).unwrap();
                let loc_type = row.get(1).unwrap();
                let variant = row.get(2).unwrap();
                let expected = row.get(3).unwrap();

                let var_n = HgvsVariant::from_str(variant)?;
                let var_g = mapper.n_to_g(&var_n)?;

                let actual = format!("{}", &NoRef(&var_g)).replace(['(', ')'], "");

                assert_eq!(
                    actual, expected,
                    "loc_type={} disc_type={}",
                    loc_type, disc_type,
                )
            }

            Ok(())
        }
    }
}
