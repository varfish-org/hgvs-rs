//! Code for variant mapping.

use std::rc::Rc;

use crate::{
    data::interface::Provider,
    parser::HgvsVariant,
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
}

impl Default for Config {
    fn default() -> Self {
        Self {
            replace_reference: true,
            strict_validation: false,
            prevalidation_level: ValidationLevel::Extrinsic,
            add_gene_symbol: false,
        }
    }
}

pub struct Mapper {
    config: Config,
    provider: Rc<dyn Provider>,
    validator: Box<dyn Validator>,
}

impl Mapper {
    fn new(config: &Config, provider: Rc<dyn Provider>) -> Mapper {
        let validator = config
            .prevalidation_level
            .validator(config.strict_validation, provider.clone());
        Mapper {
            config: config.clone(),
            provider,
            validator,
        }
    }

    fn config(&self) -> &Config {
        &self.config
    }

    /// Convert from genome to transcription.
    ///
    /// # Args
    ///
    /// * `var_g` -- `HgvsVariant::GenomeVariant` to project
    /// * `tx_ac` -- accession of transcript to project to
    /// * `alt_al_method` -- alignment method, e.g., `splign`
    fn g_to_t(
        &self,
        var_g: HgvsVariant,
        _tx_ac: &str,
        _alt_aln_method: &str,
    ) -> Result<(), anyhow::Error> {
        if let HgvsVariant::GenomeVariant {
            accession: _,
            gene_symbol: _,
            loc_edit: _,
        } = &var_g
        {
            self.validator.validate(&var_g)?;
            let _var_g = var_g.fill_ref(self.provider.as_ref())?;

            Ok(())
        } else {
            Err(anyhow::anyhow!(
                "Expected a GenomeVariant but received {}",
                &var_g
            ))
        }
    }
}
