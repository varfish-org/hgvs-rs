//! Code for enabling UTA access where sequences from a SeqRepo.
//!
//! * https://github.com/biocommons/uta
//! * https://github.com/biocommons/biocommons.seqrepo
//! * https://github.com/bihealth/seqrepo-rs

use std::path::PathBuf;
use std::rc::Rc;

use crate::data::uta::{Config as UtaConfig, Provider as UtaProvider};
use crate::data::{
    interface::GeneInfoRecord, interface::Provider as ProviderInterface, interface::TxExonsRecord,
    interface::TxForRegionRecord, interface::TxIdentityInfo, interface::TxInfoRecord,
    interface::TxMappingOptionsRecord, interface::TxSimilarityRecord,
};
use seqrepo::{AliasOrSeqId, Interface as SeqRepoInterface, SeqRepo};

/// Configurationf or the `data::uta::Provider`.
#[derive(Debug, PartialEq, Clone)]
pub struct Config {
    /// URL with the connection string, e.g.
    /// `"postgresql://anonymous:anonymous@uta.biocommons.org/uta'"`.
    pub db_url: String,
    /// The databaser schema to use, corresponds to the data version, e.g., `uta_20210129`.
    pub db_schema: String,
    /// Path to the seqrepo directory, e.g., `/usr/local/share/seqrepo/latest`.  The last path
    /// component is the "instance" name.
    pub seqrepo_path: String,
}

/// This provider provides information from a UTA and a SeqRepo.
///
/// Transcripts from a UTA Postgres database, sequences comes from a SeqRepo.  This makes
/// genome contig information available in contrast to `data::uta::Provider`.
pub struct Provider {
    inner: UtaProvider,
    seqrepo: Rc<dyn SeqRepoInterface>,
}

impl Provider {
    /// Create a new provider that uses UTA and SeqRepo information from the given configuration.
    ///
    /// This uses `seqrepo::SeqRepo` for the sequence repository.  You can inject any
    /// `seqrepo::Interface` implementation using `Provider::with_seqrepo`.
    pub fn new(config: Config) -> Result<Provider, anyhow::Error> {
        let seqrepo = PathBuf::from(&config.seqrepo_path);
        let path = seqrepo
            .parent()
            .ok_or(anyhow::anyhow!(
                "Could not get parent from {}",
                &config.seqrepo_path
            ))?
            .to_str()
            .unwrap()
            .to_string();
        let instance = seqrepo
            .file_name()
            .ok_or(anyhow::anyhow!(
                "Could not get basename from {}",
                &config.seqrepo_path
            ))?
            .to_str()
            .unwrap()
            .to_string();

        Ok(Self {
            inner: UtaProvider::with_config(&UtaConfig {
                db_url: config.db_url.clone(),
                db_schema: config.db_schema,
            })?,
            seqrepo: Rc::new(SeqRepo::new(path, &instance)?),
        })
    }

    /// Create a new provider allowing to inject a seqrepo.
    pub fn with_seqrepo(
        config: Config,
        seqrepo: Rc<dyn SeqRepoInterface>,
    ) -> Result<Provider, anyhow::Error> {
        Ok(Self {
            inner: UtaProvider::with_config(&UtaConfig {
                db_url: config.db_url.clone(),
                db_schema: config.db_schema,
            })?,
            seqrepo,
        })
    }
}

impl ProviderInterface for Provider {
    fn data_version(&self) -> &str {
        self.inner.data_version()
    }

    fn schema_version(&self) -> &str {
        self.inner.schema_version()
    }

    fn get_assembly_map(
        &self,
        assembly: crate::static_data::Assembly,
    ) -> linked_hash_map::LinkedHashMap<String, String> {
        self.inner.get_assembly_map(assembly)
    }

    fn get_gene_info(&self, hgnc: &str) -> Result<GeneInfoRecord, anyhow::Error> {
        self.inner.get_gene_info(hgnc)
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, anyhow::Error> {
        self.inner.get_pro_ac_for_tx_ac(tx_ac)
    }

    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, anyhow::Error> {
        let aos = AliasOrSeqId::Alias {
            value: ac.to_owned(),
            namespace: None,
        };
        self.seqrepo.fetch_sequence_part(&aos, begin, end)
    }

    fn get_similar_transcripts(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<TxSimilarityRecord>, anyhow::Error> {
        self.inner.get_similar_transcripts(tx_ac)
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, anyhow::Error> {
        self.inner.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, anyhow::Error> {
        self.inner.get_tx_for_gene(gene)
    }

    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, anyhow::Error> {
        self.inner
            .get_tx_for_region(alt_ac, alt_aln_method, start_i, end_i)
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, anyhow::Error> {
        self.inner.get_tx_identity_info(tx_ac)
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, anyhow::Error> {
        self.inner.get_tx_info(tx_ac, alt_ac, alt_aln_method)
    }

    fn get_tx_mapping_options(
        &self,
        tax_ac: &str,
    ) -> Result<Vec<TxMappingOptionsRecord>, anyhow::Error> {
        self.inner.get_tx_mapping_options(tax_ac)
    }
}

/// Code for helping setup of UTA providers, e.g., for setting up caching of SeqRepo results.
pub mod test_helpers {
    use std::{path::PathBuf, rc::Rc};

    use seqrepo::{
        CacheReadingSeqRepo, CacheWritingSeqRepo, Interface as SeqRepoInterface, SeqRepo,
    };

    use super::{Config, Provider, ProviderInterface};

    /// Setup a UTA Provider with data source depending on environment variables.
    ///
    /// See README.md for information on environment variable setup.
    pub fn build_provider() -> Result<Rc<dyn ProviderInterface>, anyhow::Error> {
        let db_url = std::env::var("TEST_UTA_DATABASE_URL")
            .expect("Environment variable TEST_UTA_DATABASE_URL undefined!");
        let db_schema = std::env::var("TEST_UTA_DATABASE_SCHEMA")
            .expect("Environment variable TEST_UTA_DATABASE_SCHEMA undefined!");
        let sr_cache_mode = std::env::var("TEST_SEQREPO_CACHE_MODE")
            .expect("Environment variable TEST_SEQREPO_CACHE_MODE undefined!");
        let sr_cache_path = std::env::var("TEST_SEQREPO_CACHE_PATH")
            .expect("Environment variable TEST_SEQREPO_CACHE_PATH undefined!");

        let (seqrepo, seqrepo_path) = if sr_cache_mode == "read" {
            let seqrepo: Rc<dyn SeqRepoInterface> =
                Rc::new(CacheReadingSeqRepo::new(sr_cache_path)?);
            (seqrepo, "".to_string())
        } else if sr_cache_mode == "write" {
            build_writing_sr(sr_cache_path)?
        } else {
            panic!("Invalid cache mode {}", &sr_cache_mode);
        };

        Ok(Rc::new(Provider::with_seqrepo(
            Config {
                db_url,
                db_schema,
                seqrepo_path,
            },
            seqrepo,
        )?))
    }

    /// Helper that builds the cache writing SeqRepo with inner stock SeqRepo.
    fn build_writing_sr(
        sr_cache_path: String,
    ) -> Result<(Rc<dyn SeqRepoInterface>, String), anyhow::Error> {
        let seqrepo_path = std::env::var("TEST_SEQREPO_PATH")
            .expect("Environment variable TEST_SEQREPO_PATH undefined!");
        let path_buf = PathBuf::from(seqrepo_path.clone());
        let path = path_buf
            .parent()
            .ok_or(anyhow::anyhow!(
                "Could not get parent from {}",
                &seqrepo_path
            ))?
            .to_str()
            .unwrap()
            .to_string();
        let instance = path_buf
            .file_name()
            .ok_or(anyhow::anyhow!(
                "Could not get basename from {}",
                &seqrepo_path
            ))?
            .to_str()
            .unwrap()
            .to_string();
        let seqrepo: Rc<dyn SeqRepoInterface> = Rc::new(CacheWritingSeqRepo::new(
            SeqRepo::new(path, &instance)?,
            &sr_cache_path,
        )?);
        Ok((seqrepo, seqrepo_path))
    }
}
