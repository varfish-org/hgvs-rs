//! Code for enabling UTA access.
//!
//! C.f. https://github.com/biocommons/uta

use linked_hash_map::LinkedHashMap;
use postgres::{Client, NoTls, Row};
use std::fmt::Debug;

use crate::static_data::{Assembly, ASSEMBLY_INFOS};

use super::{
    GeneInfoRecord, Interface, TxExonsRecord, TxForRegionRecord, TxIdentityInfo, TxInfoRecord,
    TxMappingOptionsRecord, TxSimilarityRecord,
};

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

impl TryFrom<Row> for GeneInfoRecord {
    type Error = anyhow::Error;

    fn try_from(row: Row) -> Result<Self, Self::Error> {
        let aliases: String = row.try_get("aliases")?;
        let aliases = aliases.split(',').map(|s| s.to_owned()).collect::<Vec<_>>();
        Ok(Self {
            hgnc: row.try_get("hgnc")?,
            maploc: row.try_get("maploc")?,
            descr: row.try_get("descr")?,
            summary: row.try_get("summary")?,
            aliases,
            added: row.try_get("added")?,
        })
    }
}

impl TryFrom<Row> for TxSimilarityRecord {
    type Error = anyhow::Error;

    fn try_from(row: Row) -> Result<Self, Self::Error> {
        Ok(Self {
            tx_ac1: row.try_get("tx_ac1")?,
            tx_ac2: row.try_get("tx_ac2")?,
            hgnc_eq: row.try_get("hgnc_eq").unwrap_or(false),
            cds_eq: row.try_get("cds_eq").unwrap_or(false),
            es_fp_eq: row.try_get("es_fp_eq").unwrap_or(false),
            cds_es_fp_eq: row.try_get("cds_es_fp_eq").unwrap_or(false),
            cds_exon_lengths_fp_eq: row.try_get("cds_exon_lengths_fp_eq").unwrap_or(false),
        })
    }
}

impl TryFrom<Row> for TxExonsRecord {
    type Error = anyhow::Error;

    fn try_from(row: Row) -> Result<Self, Self::Error> {
        Ok(Self {
            hgnc: row.try_get("hgnc")?,
            tx_ac: row.try_get("tx_ac")?,
            alt_ac: row.try_get("alt_ac")?,
            alt_aln_method: row.try_get("alt_aln_method")?,
            alt_strand: row.try_get("alt_strand")?,
            ord: row.try_get("ord")?,
            tx_start_i: row.try_get("tx_start_i")?,
            tx_end_i: row.try_get("tx_end_i")?,
            alt_start_i: row.try_get("alt_start_i")?,
            alt_end_i: row.try_get("alt_end_i")?,
            cigar: row.try_get("cigar")?,
            tx_aseq: row.try_get("tx_aseq")?,
            alt_aseq: row.try_get("alt_aseq")?,
            tx_exon_set_id: row.try_get("tx_exon_set_id")?,
            alt_exon_set_id: row.try_get("alt_exon_set_id")?,
            tx_exon_id: row.try_get("tx_exon_id")?,
            alt_exon_id: row.try_get("alt_exon_id")?,
            exon_aln_id: row.try_get("exon_aln_id")?,
        })
    }
}

impl TryFrom<Row> for TxForRegionRecord {
    type Error = anyhow::Error;

    fn try_from(row: Row) -> Result<Self, Self::Error> {
        Ok(Self {
            tx_ac: row.try_get("tx_ac")?,
            alt_ac: row.try_get("alt_ac")?,
            alt_strand: row.try_get("alt_strand")?,
            alt_aln_method: row.try_get("alt_aln_method")?,
            start_i: row.try_get("start_i")?,
            end_i: row.try_get("end_i")?,
        })
    }
}

impl TryFrom<Row> for TxIdentityInfo {
    type Error = anyhow::Error;

    fn try_from(row: Row) -> Result<Self, Self::Error> {
        Ok(Self {
            tx_ac: row.try_get("tx_ac")?,
            alt_ac: row.try_get("alt_ac")?,
            alt_aln_method: row.try_get("alt_aln_method")?,
            cds_start_i: row.try_get("cds_start_i")?,
            cds_end_i: row.try_get("cds_end_i")?,
            lengths: row.try_get("lengths")?,
            hgnc: row.try_get("hgnc")?,
        })
    }
}

impl TryFrom<Row> for TxInfoRecord {
    type Error = anyhow::Error;

    fn try_from(row: Row) -> Result<Self, Self::Error> {
        Ok(Self {
            hgnc: row.try_get("hgnc")?,
            cds_start_i: row.try_get("cds_start_i")?,
            cds_end_i: row.try_get("cds_end_i")?,
            tx_ac: row.try_get("tx_ac")?,
            alt_ac: row.try_get("alt_ac")?,
            alt_aln_method: row.try_get("alt_aln_method")?,
        })
    }
}

impl TryFrom<Row> for TxMappingOptionsRecord {
    type Error = anyhow::Error;

    fn try_from(row: Row) -> Result<Self, Self::Error> {
        Ok(Self {
            tx_ac: row.try_get("tx_ac")?,
            alt_ac: row.try_get("alt_ac")?,
            alt_aln_method: row.try_get("alt_aln_method")?,
        })
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
        let sql = format!("select key, value from {db_schema}.meta where key = 'schema_version'");
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

    fn get_gene_info(&mut self, hgnc: &str) -> Result<GeneInfoRecord, anyhow::Error> {
        let sql = format!(
            "SELECT * FROM {}.gene WHERE hgnc = $1",
            self.config.db_schema
        );
        self.conn.query_one(&sql, &[&hgnc])?.try_into()
    }

    fn get_pro_ac_for_tx_ac(&mut self, tx_ac: &str) -> Result<Option<String>, anyhow::Error> {
        let sql = format!(
            "SELECT pro_ac FROM {}.associated_accessions \
            WHERE tx_ac = $1 ORDER BY pro_ac DESC",
            self.config.db_schema
        );
        if let Some(row) = (self.conn.query(&sql, &[&tx_ac])?).into_iter().next() {
            return Ok(Some(row.try_get("pro_ac")?));
        }
        Ok(None)
    }

    fn get_seq(&mut self, ac: &str) -> Result<String, anyhow::Error> {
        self.get_seq_part(ac, None, None)
    }

    fn get_seq_part(
        &mut self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, anyhow::Error> {
        let sql = format!(
            "SELECT seq_id FROM {}.seq_anno WHERE ac = $1",
            self.config.db_schema
        );
        let seq_id: String = self.conn.query_one(&sql, &[&ac])?.try_get("seq_id")?;

        let sql = format!(
            "SELECT seq FROM {}.seq WHERE seq_id = $1",
            self.config.db_schema
        );
        let seq: String = self.conn.query_one(&sql, &[&seq_id])?.try_get("seq")?;

        let begin = begin.unwrap_or_default();
        let end = end
            .map(|end| std::cmp::min(end, seq.len()))
            .unwrap_or(seq.len());
        Ok(seq[begin..end].into())
    }

    fn get_similar_transcripts(
        &mut self,
        tx_ac: &str,
    ) -> Result<Vec<super::TxSimilarityRecord>, anyhow::Error> {
        let sql = format!(
            "SELECT * FROM {}.tx_similarity_v WHERE tx_ac1 = $1",
            self.config.db_schema
        );
        let mut result = Vec::new();
        for row in self.conn.query(&sql, &[&tx_ac])? {
            result.push(row.try_into()?);
        }
        Ok(result)
    }

    fn get_tx_exons(
        &mut self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<super::TxExonsRecord>, anyhow::Error> {
        let sql = format!(
            "SELECT * FROM {}.tx_exon_aln_v \
            WHERE tx_ac = $1 AND alt_ac = $2 and alt_aln_method = $3 \
            ORDER BY tx_start_i, tx_end_i, alt_start_i, alt_end_i, \
                     tx_exon_set_id, alt_exon_set_id",
            self.config.db_schema
        );
        let mut result = Vec::new();
        for row in self.conn.query(&sql, &[&tx_ac, &alt_ac, &alt_aln_method])? {
            result.push(row.try_into()?);
        }
        if result.is_empty() {
            Err(anyhow::anyhow!(
                "No tx_exons for tx_ac={}, alt_ac={}, alt_aln_method={}",
                &tx_ac,
                &alt_ac,
                &alt_aln_method
            ))
        } else {
            Ok(result)
        }
    }

    fn get_tx_for_gene(&mut self, gene: &str) -> Result<Vec<TxInfoRecord>, anyhow::Error> {
        let sql = format!(
            "SELECT hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method \
            FROM {}.transcript T \
            JOIN {}.exon_set es ON T.ac=ES.tx_ac WHERE alt_aln_method != 'transcript' \
            AND hgnc = $1
            ORDER BY hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method",
            self.config.db_schema, self.config.db_schema,
        );
        let mut result = Vec::new();
        for row in self.conn.query(&sql, &[&gene])? {
            result.push(row.try_into()?);
        }
        Ok(result)
    }

    fn get_tx_for_region(
        &mut self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, anyhow::Error> {
        let sql = format!(
            "SELECT tx_ac, alt_ac, alt_strand, alt_aln_method, \
                    min(start_i) AS start_i, max(end_i) AS end_i \
            FROM {}.exon_set es \
            JOIN {}.exon e ON es.exon_set_id = e.exon_set_id \
            WHERE alt_ac = $1
            GROUP BY tx_ac, alt_ac, alt_strand, alt_aln_method
            HAVING MIN(start_i) < $2 AND $3 <= MAX(end_i)
            ORDER BY tx_ac, alt_ac, alt_strand, alt_aln_method, start_i, end_i",
            self.config.db_schema, self.config.db_schema,
        );
        let mut result = Vec::new();
        for row in self.conn.query(&sql, &[&alt_ac, &start_i, &end_i])? {
            let record: TxForRegionRecord = row.try_into()?;
            // NB: The original Python code did not use alt_aln_method in the query either.
            if record.alt_aln_method == alt_aln_method {
                result.push(record);
            }
        }
        Ok(result)
    }

    fn get_tx_identity_info(&mut self, tx_ac: &str) -> Result<TxIdentityInfo, anyhow::Error> {
        let sql = format!(
            "SELECT DISTINCT(tx_ac), alt_ac, alt_aln_method, cds_start_i, \
                    cds_end_i, lengths, hgnc \
            FROM {}.tx_def_summary_v \
            WHERE tx_ac = $1 \
            ORDER BY tx_ac, alt_ac, alt_aln_method, cds_start_i, cds_end_i, lengths, hgnc",
            self.config.db_schema
        );
        self.conn.query_one(&sql, &[&tx_ac])?.try_into()
    }

    fn get_tx_info(
        &mut self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, anyhow::Error> {
        let sql = format!(
            "SELECT hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method
            FROM {}.transcript t \
            JOIN {}.exon_set es ON t.ac=es.tx_ac \
            WHERE tx_ac = $1 AND alt_ac = $2 AND alt_aln_method = $3 \
            ORDER BY hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method",
            self.config.db_schema, self.config.db_schema,
        );
        self.conn
            .query_one(&sql, &[&tx_ac, &alt_ac, &alt_aln_method])?
            .try_into()
    }

    fn get_tx_mapping_options(
        &mut self,
        tx_ac: &str,
    ) -> Result<Vec<TxMappingOptionsRecord>, anyhow::Error> {
        let sql = format!(
            "SELECT DISTINCT tx_ac, alt_ac, alt_aln_method \
            FROM {}.tx_exon_aln_v \
            WHERE tx_ac = $1 AND exon_aln_id IS NOT NULL \
            ORDER BY tx_ac, alt_ac, alt_aln_method",
            self.config.db_schema
        );
        let mut result = Vec::new();
        for row in self.conn.query(&sql, &[&tx_ac])? {
            result.push(row.try_into()?);
        }
        Ok(result)
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

    #[test]
    fn get_gene_info() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        assert_eq!(
            format!("{:?}", provider.get_gene_info("OMA1")?),
            "GeneInfoRecord { hgnc: \"OMA1\", maploc: \"1p32.2-p32.1\", \
            descr: \"OMA1 zinc metallopeptidase\", summary: \"OMA1 zinc metallopeptidase\", \
            aliases: [\"{2010001O09Rik\", \"DAB1\", \"MPRP-1\", \"MPRP1\", \"YKR087C\", \
            \"ZMPOMA1\", \"peptidase}\"], added: 2014-02-10T22:59:21.153414 }"
        );

        Ok(())
    }

    #[test]
    fn get_pro_ac_for_tx_ac() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        assert_eq!(
            provider.get_pro_ac_for_tx_ac("NM_130831.2")?,
            Some("NP_570844.1".to_string())
        );
        assert_eq!(provider.get_pro_ac_for_tx_ac("NM_130831.x")?, None);

        Ok(())
    }

    #[test]
    fn get_seq() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        assert_eq!(provider.get_seq("NM_001354664.1")?.len(), 6386);

        Ok(())
    }

    #[test]
    fn get_seq_part() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        assert_eq!(
            provider
                .get_seq_part("NM_001354664.1", Some(10), Some(100))?
                .len(),
            90
        );

        Ok(())
    }

    #[test]
    fn get_similar_transcripts() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        let records = provider.get_similar_transcripts("NM_001354664.1")?;

        assert_eq!(records.len(), 38,);
        assert_eq!(
            format!("{:?}", &records[0]),
            "TxSimilarityRecord { tx_ac1: \"NM_001354664.1\", tx_ac2: \
            \"ENST00000361150\", hgnc_eq: true, cds_eq: false, es_fp_eq: \
            false, cds_es_fp_eq: false, cds_exon_lengths_fp_eq: false }",
        );

        Ok(())
    }

    #[test]
    fn get_tx_exons() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        let records = provider.get_tx_exons("NM_001354664.1", "NC_000003.11", "splign")?;

        assert_eq!(records.len(), 30,);
        assert_eq!(
            format!("{:?}", &records[0]),
            "TxExonsRecord { hgnc: \"OPA1\", tx_ac: \"NM_001354664.1\", \
            alt_ac: \"NC_000003.11\", alt_aln_method: \"splign\", alt_strand: 1, \
            ord: 0, tx_start_i: 0, tx_end_i: 266, alt_start_i: 193310932, alt_end_i: \
            193311198, cigar: \"266=\", tx_aseq: None, alt_aseq: None, \
            tx_exon_set_id: 837345, alt_exon_set_id: 840099, tx_exon_id: 7236723, \
            alt_exon_id: 7270903, exon_aln_id: 4626713 }",
        );

        Ok(())
    }

    #[test]
    fn get_tx_for_gene() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        let records = provider.get_tx_for_gene("OMA1")?;

        assert_eq!(records.len(), 21,);
        assert_eq!(
            format!("{:?}", &records[0]),
            "TxInfoRecord { hgnc: \"OMA1\", cds_start_i: Some(0), cds_end_i: \
            Some(985), tx_ac: \"ENST00000421528\", alt_ac: \"NC_000001.10\", \
            alt_aln_method: \"genebuild\" }",
        );

        Ok(())
    }

    #[test]
    fn get_tx_for_region() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        let records = provider.get_tx_for_region("NC_000001.10", "splign", 58946391, 59012446)?;

        assert_eq!(records.len(), 2,);
        assert_eq!(
            format!("{:?}", &records[0]),
            "TxForRegionRecord { tx_ac: \"NM_145243.3\", alt_ac: \"NC_000001.10\", \
            alt_strand: -1, alt_aln_method: \"splign\", start_i: 58946390, \
            end_i: 59012446 }",
        );

        Ok(())
    }

    #[test]
    fn get_tx_identity_info() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        let record = provider.get_tx_identity_info("ENST00000421528")?;

        assert_eq!(
            format!("{:?}", &record),
            "TxIdentityInfo { tx_ac: \"ENST00000421528\", alt_ac: \"ENST00000421528\", \
            alt_aln_method: \"transcript\", cds_start_i: 0, cds_end_i: 985, lengths: \
            [24, 229, 174, 108, 129, 75, 150, 143, 1073], hgnc: \"OMA1\" }",
        );

        Ok(())
    }

    #[test]
    fn get_tx_info() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        let record = provider.get_tx_info("ENST00000421528", "NC_000001.10", "genebuild")?;

        assert_eq!(
            format!("{:?}", &record),
            "TxInfoRecord { hgnc: \"OMA1\", cds_start_i: Some(0), cds_end_i: \
            Some(985), tx_ac: \"ENST00000421528\", alt_ac: \"NC_000001.10\", \
            alt_aln_method: \"genebuild\" }",
        );

        Ok(())
    }

    #[test]
    fn get_tx_mapping_options() -> Result<(), anyhow::Error> {
        let mut provider = Provider::with_config(&get_config())?;

        let records = provider.get_tx_mapping_options("ENST00000421528")?;

        assert_eq!(records.len(), 1);
        assert_eq!(
            format!("{:?}", &records[0]),
            "TxMappingOptionsRecord { tx_ac: \"ENST00000421528\", alt_ac: \"NC_000001.10\", \
            alt_aln_method: \"genebuild\" }",
        );

        Ok(())
    }
}
