//! Code for enabling UTA access.
//!
//! C.f. https://github.com/biocommons/uta

use linked_hash_map::LinkedHashMap;
use postgres::{Client, NoTls};
use std::fmt::Debug;

use crate::static_data::{Assembly, ASSEMBLY_INFOS};

use super::Interface;

/// Configurationf or the `data::uta::Provider`.
#[derive(Debug, PartialEq, Clone)]
pub struct Config {
    /// URL with the connection string, e.g.
    /// `"postgresql://anonymous:anonymous@uta.biocommons.org/uta'"`.
    pub db_url: String,
    /// The databaser schema to use, corresponds to the data version, e.g.,
    /// `uta_20210129`.
    pub db_schema: String,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            db_url: "postgresql://anonymous:anonymous@uta.biocommons.org:5432\
            /uta"
                .to_string(),
            db_schema: "uta_20210129".to_string(),
        }
    }
}

pub struct Provider {
    /// Configuration for the access.
    config: Config,
    /// Connection to the postgres database.
    conn: Client,
    /// The schema version, set on creation.
    schema_version: String,
}

impl Debug for Provider {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Provider")
            .field("config", &self.config)
            .field("conn", &"...")
            .finish()
    }
}

impl Provider {
    pub fn with_config(config: &Config) -> Result<Self, anyhow::Error> {
        let config = config.clone();
        let mut conn = Client::connect(&config.db_url, NoTls)?;
        let schema_version = Self::fetch_schema_version(&mut conn, &config.db_schema)?;
        Ok(Self {
            config,
            conn,
            schema_version,
        })
    }

    fn fetch_schema_version(conn: &mut Client, db_schema: &str) -> Result<String, anyhow::Error> {
        let sql = format!(
            "select key, value from {}.meta where key = 'schema_version'",
            db_schema
        );
        let row = conn.query_one(&sql, &[])?;
        Ok(row.get("value"))
    }
}

impl Interface for Provider {
    fn data_version(&self) -> &str {
        &self.config.db_schema
    }

    fn schema_version(&self) -> &str {
        &self.schema_version
    }

    fn get_assembly_map(&self, assembly: Assembly) -> LinkedHashMap<String, String> {
        LinkedHashMap::from_iter(
            ASSEMBLY_INFOS[assembly]
                .sequences
                .iter()
                .map(|record| (record.refseq_ac.clone(), record.name.clone())),
        )
    }

    fn get_gene_info(&self, gene: &str) -> Result<super::GeneInfo, anyhow::Error> {
        todo!()
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, anyhow::Error> {
        todo!()
    }

    fn get_seq(&self, ac: &str) -> String {
        todo!()
    }

    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, anyhow::Error> {
        todo!()
    }

    fn get_similar_transcripts(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<super::TxSimilarityRecord>, anyhow::Error> {
        todo!()
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<super::TxExonsRecord>, anyhow::Error> {
        todo!()
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<super::TxForGeneRecord>, anyhow::Error> {
        todo!()
    }

    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<super::TxForRegionRecord>, anyhow::Error> {
        todo!()
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<super::TxIdentityInfo, anyhow::Error> {
        todo!()
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<super::TxInfoRecord, anyhow::Error> {
        todo!()
    }

    fn get_tx_mapping_options(
        &self,
        tax_ac: &str,
    ) -> Result<Vec<super::TxMappingOptionsRecord>, anyhow::Error> {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use crate::{data::Interface, static_data::Assembly};

    use super::{Config, Provider};

    fn get_config() -> Config {
        Config {
            db_url: std::env::var("TEST_UTA_DATABASE_URL")
                .expect("Environment variable TEST_UTA_DATABASE_URL undefined!"),
            db_schema: std::env::var("TEST_UTA_DATABASE_SCHEMA")
                .expect("Environment variable TEST_UTA_DATABASE_SCHEMA undefined!"),
        }
    }

    #[test]
    fn construction() -> Result<(), anyhow::Error> {
        let config = get_config();
        let provider = Provider::with_config(&config)?;

        assert_eq!(provider.data_version(), config.db_schema);
        assert_eq!(provider.schema_version(), "1.1");

        Ok(())
    }

    #[test]
    fn get_assembly_map() -> Result<(), anyhow::Error> {
        let provider = Provider::with_config(&get_config())?;

        let am37 = provider.get_assembly_map(Assembly::Grch37);
        assert_eq!(am37.len(), 92);
        assert_eq!(am37.get("NC_000001.10"), Some(&"1".to_string()));

        let am38 = provider.get_assembly_map(Assembly::Grch38);
        assert_eq!(am38.len(), 455);
        assert_eq!(am38.get("NC_000001.11"), Some(&"1".to_string()));

        Ok(())
    }
}
