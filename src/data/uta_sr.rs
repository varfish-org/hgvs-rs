//! Code for enabling UTA access where sequences from a SeqRepo.
//!
//! * https://github.com/biocommons/uta
//! * https://github.com/biocommons/biocommons.seqrepo
//! * https://github.com/bihealth/seqrepo-rs

use std::path::PathBuf;
use std::rc::Rc;

use crate::data::uta;
use crate::data::{
    error::Error, interface, interface::GeneInfoRecord, interface::TxExonsRecord,
    interface::TxForRegionRecord, interface::TxIdentityInfo, interface::TxInfoRecord,
    interface::TxMappingOptionsRecord, interface::TxSimilarityRecord,
};
use seqrepo::{self, AliasOrSeqId, SeqRepo};

/// Configuration for the `data::uta_sr::Provider`.
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
    inner: uta::Provider,
    seqrepo: Rc<dyn seqrepo::Interface>,
}

impl Provider {
    /// Create a new provider that uses UTA and SeqRepo information from the given configuration.
    ///
    /// This uses `seqrepo::SeqRepo` for the sequence repository.  You can inject any
    /// `seqrepo::Interface` implementation using `Provider::with_seqrepo`.
    pub fn new(config: Config) -> Result<Provider, Error> {
        let seqrepo = PathBuf::from(&config.seqrepo_path);
        let path = seqrepo
            .parent()
            .ok_or(Error::PathParent(config.seqrepo_path.clone()))?
            .to_str()
            .expect("problem with path to string conversion")
            .to_string();
        let instance = seqrepo
            .file_name()
            .ok_or(Error::PathBasename(config.seqrepo_path.clone()))?
            .to_str()
            .expect("problem with path to string conversion")
            .to_string();

        Ok(Self {
            inner: uta::Provider::with_config(&uta::Config {
                db_url: config.db_url.clone(),
                db_schema: config.db_schema,
            })?,
            seqrepo: Rc::new(SeqRepo::new(path, &instance)?),
        })
    }

    /// Create a new provider allowing to inject a seqrepo.
    pub fn with_seqrepo(
        config: Config,
        seqrepo: Rc<dyn seqrepo::Interface>,
    ) -> Result<Provider, Error> {
        Ok(Self {
            inner: uta::Provider::with_config(&uta::Config {
                db_url: config.db_url.clone(),
                db_schema: config.db_schema,
            })?,
            seqrepo,
        })
    }
}

impl interface::Provider for Provider {
    fn data_version(&self) -> &str {
        self.inner.data_version()
    }

    fn schema_version(&self) -> &str {
        self.inner.schema_version()
    }

    fn get_assembly_map(
        &self,
        assembly: crate::static_data::Assembly,
    ) -> indexmap::IndexMap<String, String> {
        self.inner.get_assembly_map(assembly)
    }

    fn get_gene_info(&self, hgnc: &str) -> Result<GeneInfoRecord, Error> {
        self.inner.get_gene_info(hgnc)
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, Error> {
        self.inner.get_pro_ac_for_tx_ac(tx_ac)
    }

    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, Error> {
        let aos = AliasOrSeqId::Alias {
            value: ac.to_owned(),
            namespace: None,
        };
        self.seqrepo
            .fetch_sequence_part(&aos, begin, end)
            .map_err(Error::SeqRepoError)
    }

    fn get_acs_for_protein_seq(&self, seq: &str) -> Result<Vec<String>, Error> {
        self.inner.get_acs_for_protein_seq(seq)
    }

    fn get_similar_transcripts(&self, tx_ac: &str) -> Result<Vec<TxSimilarityRecord>, Error> {
        self.inner.get_similar_transcripts(tx_ac)
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, Error> {
        self.inner.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, Error> {
        self.inner.get_tx_for_gene(gene)
    }

    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, Error> {
        self.inner
            .get_tx_for_region(alt_ac, alt_aln_method, start_i, end_i)
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, Error> {
        self.inner.get_tx_identity_info(tx_ac)
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, Error> {
        self.inner.get_tx_info(tx_ac, alt_ac, alt_aln_method)
    }

    fn get_tx_mapping_options(&self, tx_ac: &str) -> Result<Vec<TxMappingOptionsRecord>, Error> {
        self.inner.get_tx_mapping_options(tx_ac)
    }
}

/// Code for helping setup of UTA providers, e.g., for setting up caching of SeqRepo results.
#[cfg(test)]
pub mod test_helpers {
    use anyhow::Error;
    use seqrepo::{CacheReadingSeqRepo, CacheWritingSeqRepo, SeqRepo};
    use std::{path::PathBuf, rc::Rc};

    use crate::data::interface;

    use super::{Config, Provider};

    /// Setup a UTA Provider with data source depending on environment variables.
    ///
    /// See README.md for information on environment variable setup.
    pub fn build_provider() -> Result<Rc<dyn interface::Provider>, Error> {
        log::debug!("building provider...");
        let db_url = std::env::var("TEST_UTA_DATABASE_URL")
            .expect("Environment variable TEST_UTA_DATABASE_URL undefined!");
        let db_schema = std::env::var("TEST_UTA_DATABASE_SCHEMA")
            .expect("Environment variable TEST_UTA_DATABASE_SCHEMA undefined!");
        let sr_cache_mode = std::env::var("TEST_SEQREPO_CACHE_MODE")
            .expect("Environment variable TEST_SEQREPO_CACHE_MODE undefined!");
        let sr_cache_path = std::env::var("TEST_SEQREPO_CACHE_PATH")
            .expect("Environment variable TEST_SEQREPO_CACHE_PATH undefined!");

        let (seqrepo, seqrepo_path) = if sr_cache_mode == "read" {
            log::debug!("reading provider...");
            let seqrepo: Rc<dyn seqrepo::Interface> =
                Rc::new(CacheReadingSeqRepo::new(sr_cache_path)?);
            log::debug!("construction done...");
            (seqrepo, "".to_string())
        } else if sr_cache_mode == "write" {
            log::debug!("writing provider...");
            build_writing_sr(sr_cache_path)?
        } else {
            panic!("Invalid cache mode {}", &sr_cache_mode);
        };
        log::debug!("now returning provider...");

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
    pub fn build_writing_sr(
        sr_cache_path: String,
    ) -> Result<(Rc<dyn seqrepo::Interface>, String), Error> {
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
            .expect("problem with path to string conversion")
            .to_string();
        let instance = path_buf
            .file_name()
            .ok_or(anyhow::anyhow!(
                "Could not get basename from {}",
                &seqrepo_path
            ))?
            .to_str()
            .expect("problem with path to string conversion")
            .to_string();
        let seqrepo: Rc<dyn seqrepo::Interface> = Rc::new(CacheWritingSeqRepo::new(
            SeqRepo::new(path, &instance)?,
            sr_cache_path,
        )?);
        Ok((seqrepo, seqrepo_path))
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
