//! Access to `cdot` transcripts from local JSON files with sequences from a `seqrepo`.
//!
//! https://github.com/SACGF/cdot

use std::{collections::HashMap, path::PathBuf, rc::Rc, time::Instant};

use crate::{
    data::interface::{
        GeneInfoRecord, Provider as ProviderInterface, TxExonsRecord, TxForRegionRecord,
        TxIdentityInfo, TxInfoRecord, TxMappingOptionsRecord, TxSimilarityRecord,
    },
    static_data::{Assembly, ASSEMBLY_INFOS},
};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use chrono::NaiveDateTime;
use linked_hash_map::LinkedHashMap;
use seqrepo::{Interface as SeqRepoInterface, SeqRepo};

/// Configurationf or the `data::cdot::Provider`.
#[derive(Debug, PartialEq, Clone)]
pub struct Config {
    /// Paths to the gzip-ed JSON files to load
    pub json_paths: Vec<String>,
    /// Path to the seqrepo directory, e.g., `/usr/local/share/seqrepo/latest`.  The last path
    /// component is the "instance" name.
    pub seqrepo_path: String,
}

/// This provider provides information from a cDOT JSON and a SeqRepo.
///
/// Transcripts from a cDOT JSON Postgres database, sequences comes from a SeqRepo.
/// This makes genome contig information available.
///
/// # Remarks
///
/// The method are not implemented.
///
/// The method `get_tx_exons()` returns `None` for record entries `tx_aseq`, and `alt_aseq`
/// and `std::i32::MAX` for `tx_exon-set_id`, `alt_exon_set_id`, `tx_exon_id`, `alt_exon_id`,
/// `exon_aln_id`.
pub struct Provider {
    inner: TxProvider,
    seqrepo: Rc<dyn SeqRepoInterface>,
}

impl Provider {
    pub fn new(config: Config) -> Result<Self, anyhow::Error> {
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
            inner: TxProvider::with_config(
                config
                    .json_paths
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<&str>>()
                    .as_ref(),
            )?,
            seqrepo: Rc::new(SeqRepo::new(path, &instance)?),
        })
    }

    /// Create a new provider allowing to inject a seqrepo.
    pub fn with_seqrepo(
        config: Config,
        seqrepo: Rc<dyn SeqRepoInterface>,
    ) -> Result<Provider, anyhow::Error> {
        Ok(Self {
            inner: TxProvider::with_config(
                config
                    .json_paths
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<&str>>()
                    .as_ref(),
            )?,
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
        self.seqrepo.fetch_sequence_part(
            &seqrepo::AliasOrSeqId::Alias {
                value: ac.to_string(),
                namespace: None,
            },
            begin,
            end,
        )
    }

    fn get_acs_for_protein_seq(&self, seq: &str) -> Result<Vec<String>, anyhow::Error> {
        self.inner.get_acs_for_protein_seq(seq)
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
        tx_ac: &str,
    ) -> Result<Vec<TxMappingOptionsRecord>, anyhow::Error> {
        self.inner.get_tx_mapping_options(tx_ac)
    }
}

/// Data structures used for deserializing from cdot.
pub mod models {
    use linked_hash_map::LinkedHashMap;
    use serde::{Deserialize, Deserializer};

    /// Container for a cDot data file.
    #[derive(Deserialize, Debug, Clone)]
    pub struct Container {
        pub transcripts: LinkedHashMap<String, Transcript>,
        pub cdot_version: String,
        pub genome_builds: Vec<String>,
        pub genes: LinkedHashMap<String, Gene>,
    }

    /// Enum for representing the tags for transcripts.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Tag {
        Basic,
        EnsemblCanonical,
        ManeSelect,
        ManePlusClinical,
        RefSeqSelect,
    }

    #[derive(Deserialize, Debug, Clone)]
    pub struct Transcript {
        /// Transcript biotype, e.g., `vec![BioType::ProteinCoding]` for BRCA1.
        #[serde(default)]
        pub biotype: Option<Vec<BioType>>,
        /// Gene name, e.g., `"BRCA1"` for BRCA1.
        #[serde(default)]
        pub gene_name: Option<String>,
        /// Gene version, NCBI gene ID for NCBI, e.g., `"672"` for BRCA1, ENSEMBL gene id, e.g.,
        /// `"ENSG00000012048"` for BRCA1.
        pub gene_version: String,
        /// Alignments to the different genome builds.
        pub genome_builds: LinkedHashMap<String, GenomeAlignment>,
        /// HGNC identifier, e.g. `"1100"` for BRCA1 which is `HGNC:1100`.
        #[serde(default)]
        pub hgnc: Option<String>,
        /// Identifier of the transcript, same as key in `transcripts`, e.g., `"NM_007294.3"` for BRCA1.
        pub id: String,
        /// Identifier of corresponding protein, e.g., `"NP_009225.1"` for `"NM_007294.3"` of BRCA1.
        #[serde(default)]
        pub protein: Option<String>,
        /// Start codon of transcript, e.g., `232` for `"NM_007294.3"` of BRCA1.
        #[serde(default)]
        pub start_codon: Option<i32>,
        /// Stop codon of transcript, e.g., `232` for `"NM_007294.3"` of BRCA1.
        #[serde(default)]
        pub stop_codon: Option<i32>,
        /// Tags of the transcript.
        #[serde(default)]
        #[serde(deserialize_with = "deserialize_tag")]
        pub tag: Option<Vec<Tag>>,
    }

    /// Representation of the strand.
    #[derive(Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Strand {
        #[serde(rename = "+")]
        Plus,
        #[serde(rename = "-")]
        Minus,
    }

    /// Representation of an exon in the JSON.
    #[derive(Deserialize, Debug, Clone)]
    struct ExonHelper(i32, i32, i32, i32, i32, Option<String>);

    /// Representation of an exon after loading.
    #[derive(Debug, Clone)]
    pub struct Exon {
        /// Start position on reference.
        pub alt_start_i: i32,
        /// End position on reference.
        pub alt_end_i: i32,
        /// Exon number.
        pub ord: i32,
        /// CDS start coordinate.
        pub alt_cds_start_i: i32,
        /// CDS end coordinate.
        pub alt_cds_end_i: i32,
        /// Alignment of an exon in CIGAR format.
        pub cigar: String,
    }

    /// Representation of `transcripts.*.genome_builds` value.
    #[derive(Deserialize, Debug, Clone)]
    pub struct GenomeAlignment {
        /// CDS end position.
        #[serde(default)]
        pub cds_end: Option<i32>,
        /// CDS start position.
        #[serde(default)]
        pub cds_start: Option<i32>,
        /// ID of the contig.
        pub contig: String,
        /// List of exons.
        pub exons: Vec<Exon>,
        /// The strand.
        pub strand: Strand,
    }

    /// Enum for representing the biotypes.
    #[derive(Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
    pub enum BioType {
        #[serde(rename = "3prime_overlapping_ncrna")]
        ThreePrimeOverlappingNcRna,
        #[serde(rename = "antisense")]
        Antisense,
        #[serde(rename = "artifact")]
        Artifact,
        #[serde(rename = "IG_C_gene")]
        IgCGene,
        #[serde(rename = "IG_C_pseudogene")]
        IgCPseudogene,
        #[serde(rename = "IG_D_gene")]
        IgDGene,
        #[serde(rename = "IG_J_gene")]
        IdJGene,
        #[serde(rename = "IG_J_pseudogene")]
        IdJPseudogene,
        #[serde(rename = "IG_pseudogene")]
        IdPseudogene,
        #[serde(rename = "IG_V_gene")]
        IgVGene,
        #[serde(rename = "IG_V_pseudogene")]
        IgVPseudogene,
        #[serde(rename = "Mt_rRNA")]
        MtRna,
        #[serde(rename = "non_coding")]
        NonCoding,
        #[serde(rename = "polymorphic_pseudogene")]
        PolymorphicPseudogene,
        #[serde(rename = "processed_pseudogene")]
        ProcessedPseudogene,
        #[serde(rename = "protein_coding")]
        ProteinCoding,
        #[serde(rename = "pseudogene")]
        Pseudogene,
        #[serde(rename = "rRNA_pseudogene")]
        RnaPseudogene,
        #[serde(rename = "sense_intronic")]
        SenseIntronic,
        #[serde(rename = "sense_overlapping")]
        SenseOverlapping,
        #[serde(rename = "sRNA")]
        SRna,
        #[serde(rename = "TEC")]
        Tec,
        #[serde(rename = "transcribed_processed_pseudogene")]
        TranscriptProcessedPseudogene,
        #[serde(rename = "transcribed_unitary_pseudogene")]
        TranscriptUnitaryPseudogene,
        #[serde(rename = "transcribed_unprocessed_pseudogene")]
        TranscriptUnprocessedPseudogene,
        #[serde(rename = "translated_processed_pseudogene")]
        TranslatedProcessedPseudogene,
        #[serde(rename = "translated_unprocessed_pseudogene")]
        TranslatedUnprocessedPseudogene,
        #[serde(rename = "TR_C_gene")]
        TrCGene,
        #[serde(rename = "TR_D_gene")]
        TrDGene,
        #[serde(rename = "TR_J_gene")]
        TrJGene,
        #[serde(rename = "TR_J_pseudogene")]
        TrJPseudogene,
        #[serde(rename = "TR_V_gene")]
        TrVGene,
        #[serde(rename = "TR_V_pseudogene")]
        TrVPseudogene,
        #[serde(rename = "unitary_pseudogene")]
        UnitaryPseudogene,
        #[serde(rename = "unprocessed_pseudogene")]
        UnprocessedPseudogene,
    }

    #[derive(Deserialize, Debug, Clone)]
    pub struct Gene {
        #[serde(default)]
        #[serde(deserialize_with = "deserialize_gene_aliases")]
        pub aliases: Option<Vec<String>>,
        #[serde(default)]
        #[serde(deserialize_with = "deserialize_gene_biotype")]
        pub biotype: Option<Vec<BioType>>,
        #[serde(default)]
        pub description: Option<String>,
        pub gene_symbol: Option<String>,
        #[serde(default)]
        pub hgnc: Option<String>,
        #[serde(default)]
        pub map_location: Option<String>,
        #[serde(default)]
        pub summary: Option<String>,
        pub url: String,
    }

    /// Convert cdot gap to CIGAR string.
    ///
    /// ```
    // gap = 'M196 I1 M61 I1 M181'
    // CIGAR = '194=1D60=1D184='
    /// ```
    pub fn gap_to_cigar(gap: &str) -> String {
        let mut result = String::new();

        for gap_op in gap.split(' ') {
            result.push_str(&gap_op[1..]);
            let op = &gap_op[0..1];
            match op {
                "M" => result.push('='),
                "I" => result.push('D'),
                "D" => result.push('I'),
                _ => panic!("unknown"),
            }
        }

        result
    }

    impl<'de> Deserialize<'de> for Exon {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            Deserialize::deserialize(deserializer).map(
                |ExonHelper(alt_start_i, alt_end_i, ord, alt_cds_start_i, alt_cds_end_i, gap)| {
                    Exon {
                        alt_start_i,
                        alt_end_i,
                        ord,
                        alt_cds_start_i,
                        alt_cds_end_i,
                        cigar: if let Some(gap) = gap {
                            gap_to_cigar(&gap)
                        } else {
                            format!("{}M", alt_end_i - alt_start_i)
                        },
                    }
                },
            )
        }
    }

    fn str_to_tag(s: &str) -> Tag {
        match s {
            "basic" => Tag::Basic,
            "Ensembl_canonical" => Tag::EnsemblCanonical,
            "MANE_Plus_Clinical" => Tag::ManePlusClinical,
            "MANE_Select" => Tag::ManeSelect,
            "MANE Select" => Tag::ManeSelect,
            "RefSeq Select" => Tag::RefSeqSelect,
            _ => panic!("Invalid transcript tag {s:?}"),
        }
    }

    #[derive(Debug, Deserialize)]
    struct WrappedString(String);

    fn deserialize_tag<'de, D>(deserializer: D) -> Result<Option<Vec<Tag>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Option::<WrappedString>::deserialize(deserializer).map(|opt_wrapped| {
            opt_wrapped.map(|wrapped| wrapped.0.split(',').into_iter().map(str_to_tag).collect())
        })
    }

    fn deserialize_gene_aliases<'de, D>(deserializer: D) -> Result<Option<Vec<String>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let buf = Option::<String>::deserialize(deserializer)?;

        Ok(buf.map(|s| s.split(", ").into_iter().map(|s| s.to_string()).collect()))
    }

    fn str_to_biotype(s: &str) -> BioType {
        match s {
            "3prime_overlapping_ncrna" => BioType::ThreePrimeOverlappingNcRna,
            "antisense" => BioType::Antisense,
            "artifact" => BioType::Artifact,
            "IG_C_gene" => BioType::IgCGene,
            "IG_C_pseudogene" => BioType::IgCPseudogene,
            "IG_D_gene" => BioType::IgDGene,
            "IG_J_gene" => BioType::IdJGene,
            "IG_J_pseudogene" => BioType::IdJPseudogene,
            "IG_pseudogene" => BioType::IdPseudogene,
            "IG_V_gene" => BioType::IgVGene,
            "IG_V_pseudogene" => BioType::IgVPseudogene,
            "Mt_rRNA" => BioType::MtRna,
            "non_coding" => BioType::NonCoding,
            "polymorphic_pseudogene" => BioType::PolymorphicPseudogene,
            "processed_pseudogene" => BioType::ProcessedPseudogene,
            "protein_coding" => BioType::ProteinCoding,
            "pseudogene" => BioType::Pseudogene,
            "rRNA_pseudogene" => BioType::RnaPseudogene,
            "sense_intronic" => BioType::SenseIntronic,
            "sense_overlapping" => BioType::SenseOverlapping,
            "sRNA" => BioType::SRna,
            "TEC" => BioType::Tec,
            "transcribed_processed_pseudogene" => BioType::TranscriptProcessedPseudogene,
            "transcribed_unitary_pseudogene" => BioType::TranscriptUnitaryPseudogene,
            "transcribed_unprocessed_pseudogene" => BioType::TranscriptUnprocessedPseudogene,
            "translated_processed_pseudogene" => BioType::TranslatedProcessedPseudogene,
            "translated_unprocessed_pseudogene" => BioType::TranslatedUnprocessedPseudogene,
            "TR_C_gene" => BioType::TrCGene,
            "TR_D_gene" => BioType::TrDGene,
            "TR_J_gene" => BioType::TrJGene,
            "TR_J_pseudogene" => BioType::TrJPseudogene,
            "TR_V_gene" => BioType::TrVGene,
            "TR_V_pseudogene" => BioType::TrVPseudogene,
            "unitary_pseudogene" => BioType::UnitaryPseudogene,
            "unprocessed_pseudogene" => BioType::UnprocessedPseudogene,
            _ => panic!("Unknown biotype {s:?}"),
        }
    }

    fn deserialize_gene_biotype<'de, D>(deserializer: D) -> Result<Option<Vec<BioType>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let buf = Option::<String>::deserialize(deserializer)?;
        Ok(buf.map(|buf| {
            if buf.is_empty() {
                Vec::new()
            } else {
                buf.split(',').into_iter().map(str_to_biotype).collect()
            }
        }))
    }
}

/// Type alias for interval trees.
type IntervalTree = ArrayBackedIntervalTree<i32, String>;

/// Internal implementation of the transcript provider.
struct TxProvider {
    /// Genes by HGNC symbol.
    genes: HashMap<String, models::Gene>,
    /// Transcripts by ID.
    transcripts: HashMap<String, models::Transcript>,
    /// Transcript identifiers for each gene HGNC symbol.
    transcript_ids_for_gene: HashMap<String, Vec<String>>,
    /// Interval tree for each alternative contig.
    interval_trees: HashMap<String, IntervalTree>,
}

/// "Normal" associated functions and methods.
impl TxProvider {
    fn with_config(json_paths: &[&str]) -> Result<Self, anyhow::Error> {
        let mut genes = HashMap::new();
        let mut transcripts = HashMap::new();
        let mut transcript_ids_for_gene = HashMap::new();

        for json_path in json_paths {
            Self::load_and_extract(
                json_path,
                &mut transcript_ids_for_gene,
                &mut genes,
                &mut transcripts,
            )?;
        }

        log::debug!(
            "json::TxProvider -- #genes = {}, #transcripts = {}, #transcript_ids_for_gene = {}",
            genes.len(),
            transcripts.len(),
            transcript_ids_for_gene.len()
        );

        let interval_trees = Self::build_interval_trees(&transcripts);

        Ok(Self {
            genes,
            transcripts,
            transcript_ids_for_gene,
            interval_trees,
        })
    }

    fn load_and_extract(
        json_path: &&str,
        transcript_ids_for_gene: &mut HashMap<String, Vec<String>>,
        genes: &mut HashMap<String, models::Gene>,
        transcripts: &mut HashMap<String, models::Transcript>,
    ) -> Result<(), anyhow::Error> {
        log::debug!("Loading cdot transcripts from {:?}", json_path);
        let start = Instant::now();
        let models::Container {
            genes: c_genes,
            transcripts: c_txs,
            ..
        } = if json_path.ends_with(".gz") {
            serde_json::from_reader(flate2::bufread::GzDecoder::new(std::io::BufReader::new(
                std::fs::File::open(json_path)?,
            )))?
        } else {
            serde_json::from_reader(std::io::BufReader::new(std::fs::File::open(json_path)?))?
        };
        log::debug!(
            "loading / deserializing {} genes and {} transcripts from cdot took {:?}",
            c_genes.len(),
            c_txs.len(),
            start.elapsed()
        );

        let start = Instant::now();
        c_genes
            .values()
            .into_iter()
            .filter(|gene| {
                gene.gene_symbol.is_some()
                    && !gene.gene_symbol.as_ref().unwrap().is_empty()
                    && gene.map_location.is_some()
                    && !gene.map_location.as_ref().unwrap().is_empty()
            })
            .for_each(|gene| {
                let gene_symbol = gene.gene_symbol.as_ref().unwrap().clone();
                transcript_ids_for_gene
                    .entry(gene_symbol.clone())
                    .or_insert(Vec::new());
                genes.insert(gene_symbol, gene.clone());
            });
        c_txs
            .values()
            .into_iter()
            .filter(|tx| tx.gene_name.is_some() && !tx.gene_name.as_ref().unwrap().is_empty())
            .for_each(|tx| {
                let gene_name = tx.gene_name.as_ref().unwrap();
                transcript_ids_for_gene
                    .get_mut(gene_name)
                    .unwrap_or_else(|| panic!("tx {:?} for unknown gene {:?}", tx.id, gene_name))
                    .push(tx.id.clone());
                transcripts.insert(tx.id.clone(), tx.clone());
            });
        log::debug!("extracting datastructures took {:?}", start.elapsed());
        Ok(())
    }

    fn build_interval_trees(
        transcripts: &HashMap<String, models::Transcript>,
    ) -> HashMap<String, IntervalTree> {
        let start = Instant::now();
        log::debug!("Building interval trees...");

        let mut result = HashMap::new();
        for transcript in transcripts.values() {
            for genome_alignment in transcript.genome_builds.values() {
                let contig = &genome_alignment.contig;
                result.entry(contig.clone()).or_insert(IntervalTree::new());
                let tree = result.get_mut(contig).unwrap();
                let alt_start_i = genome_alignment
                    .exons
                    .iter()
                    .map(|exon| exon.alt_start_i)
                    .min()
                    .unwrap()
                    - 1;
                let alt_end_i = genome_alignment
                    .exons
                    .iter()
                    .map(|exon| exon.alt_end_i)
                    .max()
                    .unwrap();
                tree.insert(alt_start_i..alt_end_i, transcript.id.clone());
            }
        }

        for tree in result.values_mut() {
            tree.index();
        }

        log::debug!("Built interval trees in {:?}", start.elapsed());
        result
    }
}

/// What is returned by `Provider::schema_version()` and `Provider::data_version()`.
pub static REQUIRED_VERSION: &str = "1.1";

/// The alignment method returned for all cdot transcripts.
pub static NCBI_ALN_METHOD: &str = "splign";

/// Implementation for `ProviderInterface`.
impl TxProvider {
    fn data_version(&self) -> &str {
        REQUIRED_VERSION
    }

    fn schema_version(&self) -> &str {
        REQUIRED_VERSION
    }

    fn get_assembly_map(&self, assembly: Assembly) -> LinkedHashMap<String, String> {
        LinkedHashMap::from_iter(
            ASSEMBLY_INFOS[assembly]
                .sequences
                .iter()
                .map(|record| (record.refseq_ac.clone(), record.name.clone())),
        )
    }

    fn get_gene_info(&self, hgnc: &str) -> Result<GeneInfoRecord, anyhow::Error> {
        let gene = self
            .genes
            .get(hgnc)
            .ok_or(anyhow::anyhow!("No gene found for {:?}", hgnc))?;

        Ok(GeneInfoRecord {
            hgnc: gene
                .gene_symbol
                .as_ref()
                .expect("cannot happen; genes without symbol are not imported")
                .clone(),
            maploc: gene
                .map_location
                .as_ref()
                .expect("cannot happen; genes without map_location are not imported")
                .clone(),
            descr: gene
                .description
                .as_ref()
                .map(Clone::clone)
                .unwrap_or_default(),
            summary: gene.summary.as_ref().map(Clone::clone).unwrap_or_default(),
            aliases: gene.aliases.as_ref().map(Clone::clone).unwrap_or_default(),
            added: NaiveDateTime::default(),
        })
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, anyhow::Error> {
        let transcript = self
            .transcripts
            .get(tx_ac)
            .ok_or(anyhow::anyhow!("No transcript found for {:?}", &tx_ac))?;
        Ok(transcript.protein.clone())
    }

    /// Note from the original cdot Python code.
    ///
    /// This is not implemented. The only caller has comment: 'TODO: drop get_acs_for_protein_seq'
    /// And is only ever called as a backup when get_pro_ac_for_tx_ac fails
    #[allow(dead_code)]
    fn get_acs_for_protein_seq(&self, _seq: &str) -> Result<Vec<String>, anyhow::Error> {
        log::warn!(
            "cdot::data::json::TxProvider::get_acs_for_protein_seq() \
            This has not been implemented"
        );
        Ok(vec![])
    }

    /// Note from the original cdot Python code.
    ///
    /// UTA specific functionality that uses tx_similarity_v table
    /// This is not used by the HGVS library
    #[allow(dead_code)]
    fn get_similar_transcripts(
        &self,
        _tx_ac: &str,
    ) -> Result<Vec<TxSimilarityRecord>, anyhow::Error> {
        log::warn!(
            "cdot::data::json::TxProvider::get_similar_transcripts() \
            This has not been implemented"
        );
        Ok(vec![])
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, anyhow::Error> {
        let tx = self.transcripts.get(tx_ac).ok_or(anyhow::anyhow!(
            "Could not find transcript for {:?}",
            &tx_ac
        ))?;

        let genome_alignment = tx
            .genome_builds
            .values()
            .find(|genome_alignment| genome_alignment.contig == alt_ac);
        if let Some(genome_alignment) = genome_alignment {
            Ok(genome_alignment
                .exons
                .iter()
                .map(|exon| TxExonsRecord {
                    hgnc: tx
                        .gene_name
                        .as_ref()
                        .expect("cannot happen; transcript without gene name not imported")
                        .clone(),
                    tx_ac: tx_ac.to_string(),
                    alt_ac: genome_alignment.contig.clone(),
                    alt_aln_method: alt_aln_method.to_string(),
                    alt_strand: match genome_alignment.strand {
                        models::Strand::Plus => 1,
                        models::Strand::Minus => -1,
                    },
                    ord: exon.ord,
                    tx_start_i: exon.alt_cds_start_i - 1,
                    tx_end_i: exon.alt_cds_end_i,
                    alt_start_i: exon.alt_start_i,
                    alt_end_i: exon.alt_end_i,
                    cigar: exon.cigar.clone(),
                    tx_aseq: None,
                    alt_aseq: None,
                    tx_exon_set_id: std::i32::MAX,
                    alt_exon_set_id: std::i32::MAX,
                    tx_exon_id: std::i32::MAX,
                    alt_exon_id: std::i32::MAX,
                    exon_aln_id: std::i32::MAX,
                })
                .collect())
        } else {
            Err(anyhow::anyhow!(
                "Could not find alignment of {:?} to {:?}",
                tx_ac,
                alt_ac
            ))
        }
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, anyhow::Error> {
        if let Some(tx_acs) = self.transcript_ids_for_gene.get(gene) {
            let mut tmp = Vec::new();
            for tx_ac in tx_acs {
                let tx = self.transcripts.get(tx_ac).unwrap();
                let cds_start_i = tx.start_codon;
                let cds_end_i = tx.stop_codon;
                for genome_alignment in tx.genome_builds.values() {
                    let contig = &genome_alignment.contig;
                    let tx_start = genome_alignment.exons.first().unwrap().alt_start_i;
                    let tx_end = genome_alignment.exons.last().unwrap().alt_end_i;
                    let length = tx_end - tx_start;
                    let rec = TxInfoRecord {
                        hgnc: gene.to_string(),
                        cds_start_i,
                        cds_end_i,
                        tx_ac: tx_ac.clone(),
                        alt_ac: contig.clone(),
                        alt_aln_method: NCBI_ALN_METHOD.to_string(),
                    };
                    tmp.push((length, rec));
                }
            }

            // Sorted by length in descending order.
            tmp.sort_by(|a, b| b.0.cmp(&a.0));

            Ok(tmp.into_iter().map(|x| x.1).collect())
        } else {
            Ok(Vec::new())
        }
    }

    // NB: This implementation uses hgvs semantics, not cdot.
    //
    // cf. https://github.com/SACGF/cdot/issues/38
    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, anyhow::Error> {
        if alt_aln_method != NCBI_ALN_METHOD {
            return Ok(Vec::new());
        }

        if let Some(contig_itv_tree) = self.interval_trees.get(alt_ac) {
            let mut tmp = Vec::new();

            for entry in contig_itv_tree.find((start_i - 1)..(end_i)) {
                let tx_ac = entry.data();
                let tx = self.transcripts.get(tx_ac).unwrap();
                let genome_alignment = tx
                    .genome_builds
                    .values()
                    .into_iter()
                    .find(|genome_alignment| genome_alignment.contig == alt_ac);
                if let Some(genome_alignment) = genome_alignment {
                    let tx_start = genome_alignment.exons.first().unwrap().alt_start_i;
                    let tx_end = genome_alignment.exons.last().unwrap().alt_end_i;
                    let length = tx_end - tx_start;

                    let rec = TxForRegionRecord {
                        tx_ac: tx_ac.clone(),
                        alt_ac: alt_ac.to_string(),
                        alt_strand: match genome_alignment.strand {
                            models::Strand::Plus => 1,
                            models::Strand::Minus => -1,
                        },
                        alt_aln_method: alt_aln_method.to_string(),
                        start_i: tx_start,
                        end_i: tx_end,
                    };
                    tmp.push(((length, tx_ac.clone()), rec));
                }
            }

            // Sorted by length in descending order, break tie by tx accession.
            tmp.sort_by(|a, b| b.0.cmp(&a.0));

            Ok(tmp.into_iter().map(|x| x.1).collect())
        } else {
            Ok(Vec::new())
        }
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, anyhow::Error> {
        let tx = self.transcripts.get(tx_ac).ok_or(anyhow::anyhow!(
            "Could not find transcript for {:?}",
            &tx_ac
        ))?;

        let hgnc = tx
            .gene_name
            .as_ref()
            .expect("cannot happen; transcripts without gene_name not imported")
            .clone();

        let mut tmp = tx
            .genome_builds
            .values()
            .next()
            .unwrap()
            .exons
            .iter()
            .map(|exon| (exon.ord, exon.alt_cds_end_i + 1 - exon.alt_cds_start_i))
            .collect::<Vec<(i32, i32)>>();
        tmp.sort();

        let lengths = tmp.into_iter().map(|(_, length)| length).collect();
        Ok(TxIdentityInfo {
            tx_ac: tx_ac.to_string(),
            alt_ac: tx_ac.to_string(), // sic(!)
            alt_aln_method: String::from("transcript"),
            cds_start_i: tx.start_codon.unwrap_or_default(),
            cds_end_i: tx.stop_codon.unwrap_or_default(),
            lengths,
            hgnc,
        })
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, anyhow::Error> {
        let tx = self.transcripts.get(tx_ac).ok_or(anyhow::anyhow!(
            "Could not find transcript for {:?}",
            &tx_ac
        ))?;

        Ok(TxInfoRecord {
            hgnc: tx
                .gene_name
                .as_ref()
                .expect("cannot happen; transcript without gene name not imported")
                .clone(),
            cds_start_i: tx.start_codon,
            cds_end_i: tx.stop_codon,
            tx_ac: tx_ac.to_string(),
            alt_ac: alt_ac.to_string(),
            alt_aln_method: alt_aln_method.to_string(),
        })
    }

    fn get_tx_mapping_options(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<TxMappingOptionsRecord>, anyhow::Error> {
        let tx = self.transcripts.get(tx_ac).ok_or(anyhow::anyhow!(
            "Could not find transcript for {:?}",
            &tx_ac
        ))?;

        Ok(tx
            .genome_builds
            .values()
            .into_iter()
            .map(|genome_alignment| TxMappingOptionsRecord {
                tx_ac: tx_ac.to_string(),
                alt_ac: genome_alignment.contig.clone(),
                alt_aln_method: NCBI_ALN_METHOD.to_string(),
            })
            .collect())
    }
}

pub mod test_helpers {
    use std::rc::Rc;

    use crate::data::uta_sr::test_helpers::build_writing_sr;

    use super::{Config, Provider};
    use seqrepo::{CacheReadingSeqRepo, Interface as SeqRepoInterface};

    pub fn build_provider() -> Result<Provider, anyhow::Error> {
        let sr_cache_mode = std::env::var("TEST_SEQREPO_CACHE_MODE")
            .expect("Environment variable TEST_SEQREPO_CACHE_MODE undefined!");
        let sr_cache_path = std::env::var("TEST_SEQREPO_CACHE_PATH")
            .expect("Environment variable TEST_SEQREPO_CACHE_PATH undefined!");

        log::debug!("building provider...");
        let seqrepo = if sr_cache_mode == "read" {
            log::debug!("reading provider...");
            let seqrepo: Rc<dyn SeqRepoInterface> =
                Rc::new(CacheReadingSeqRepo::new(sr_cache_path)?);
            log::debug!("construction done...");
            seqrepo
        } else if sr_cache_mode == "write" {
            log::debug!("writing provider...");
            build_writing_sr(sr_cache_path)?.0
        } else {
            panic!("Invalid cache mode {}", &sr_cache_mode);
        };
        log::debug!("now returning provider...");

        let config = Config {
            json_paths: vec![String::from(
                "tests/data/data/cdot/cdot-0.2.12.refseq.grch37_grch38.brca1.json",
            )],
            seqrepo_path: String::from("nonexisting"),
        };
        // Note that we don't actually use the seqrepo instance except for initializing the provider.
        Provider::with_seqrepo(config, seqrepo)
    }
}

#[cfg(test)]
pub mod tests {
    use std::rc::Rc;
    use std::str::FromStr;

    use chrono::NaiveDateTime;
    use pretty_assertions::assert_eq;
    use test_log::test;

    use super::models::{gap_to_cigar, Container};
    use super::test_helpers::build_provider;
    use crate::data::interface::{
        GeneInfoRecord, Provider as ProviderInterface, TxExonsRecord, TxForRegionRecord,
        TxInfoRecord, TxMappingOptionsRecord, TxSimilarityRecord,
    };
    use crate::mapper::assembly::{Config as AssemblyMapperConfig, Mapper};
    use crate::parser::HgvsVariant;
    use crate::static_data::Assembly;

    #[test]
    fn deserialize_brca1() -> Result<(), anyhow::Error> {
        let json = std::fs::read_to_string(
            "tests/data/data/cdot/cdot-0.2.12.refseq.grch37_grch38.brca1.json",
        )?;

        let c: Container = serde_json::from_str(&json)?;
        assert_eq!(c.cdot_version, "0.2.12");

        let c_dump = format!("{:#?}\n", &c);
        let c_dump_expected = std::fs::read_to_string(
            "tests/data/data/cdot/cdot-0.2.12.refseq.grch37_grch38.brca1.txt",
        )?;
        assert_eq!(c_dump, c_dump_expected);

        Ok(())
    }

    #[test]
    fn gap_to_cigar_smoke() {
        assert_eq!(gap_to_cigar("M196 I1 M61 I1 M181"), "196=1D61=1D181=");
    }

    /// Deserialization of the big cdot files for benchmarking.
    #[cfg(deserialization_tests)]
    #[test]
    fn deserialize_big_files() -> Result<(), anyhow::Error> {
        let before = std::time::Instant::now();
        println!("ensembl...");
        let _ensembl: Container =
            serde_json::from_reader(flate2::bufread::GzDecoder::new(std::io::BufReader::new(
                std::fs::File::open("cdot-0.2.12.ensembl.grch37_grch38.json.gz")?,
            )))?;
        println!("Reading ensembl: {:?}", before.elapsed());

        let before = std::time::Instant::now();
        println!("refseq...");
        let _refseq: Container =
            serde_json::from_reader(flate2::bufread::GzDecoder::new(std::io::BufReader::new(
                std::fs::File::open("cdot-0.2.12.refseq.grch37_grch38.json.gz")?,
            )))?;
        println!("Reading refseq: {:?}", before.elapsed());

        Ok(())
    }

    #[test]
    fn provider_brca1_smoke() -> Result<(), anyhow::Error> {
        build_provider()?;
        Ok(())
    }

    #[test]
    fn provider_versions() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        assert_eq!(provider.data_version(), "1.1");
        assert_eq!(provider.schema_version(), "1.1");

        Ok(())
    }

    #[test]
    fn provider_get_assembly_map() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;
        assert_eq!(provider.get_assembly_map(Assembly::Grch37p10).len(), 275);
        assert_eq!(provider.get_assembly_map(Assembly::Grch38).len(), 455);

        Ok(())
    }

    #[test]
    fn provider_get_gene_info() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        assert!(provider.get_gene_info("BRCA2").is_err());

        let record = provider.get_gene_info("BRCA1")?;

        let expected = GeneInfoRecord {
            hgnc: String::from("BRCA1"),
            maploc: String::from("17q21.31"),
            descr: String::from("BRCA1 DNA repair associated"),
            summary: String::from("This gene encodes a 190 kD nuclear phosphoprotein that plays a role in maintaining genomic stability, and it also acts as a tumor suppressor. The BRCA1 gene contains 22 exons spanning about 110 kb of DNA. The encoded protein combines with other tumor suppressors, DNA damage sensors, and signal transducers to form a large multi-subunit protein complex known as the BRCA1-associated genome surveillance complex (BASC). This gene product associates with RNA polymerase II, and through the C-terminal domain, also interacts with histone deacetylase complexes. This protein thus plays a role in transcription, DNA repair of double-stranded breaks, and recombination. Mutations in this gene are responsible for approximately 40% of inherited breast cancers and more than 80% of inherited breast and ovarian cancers. Alternative splicing plays a role in modulating the subcellular localization and physiological function of this gene. Many alternatively spliced transcript variants, some of which are disease-associated mutations, have been described for this gene, but the full-length natures of only some of these variants has been described. A related pseudogene, which is also located on chromosome 17, has been identified. [provided by RefSeq, May 2020]"),
            aliases: vec![
                String::from("BRCAI"),
                String::from("BRCC1"),
                String::from("BROVCA1"),
                String::from("FANCS"),
                String::from("IRIS"),
                String::from("PNCA4"),
                String::from("PPP1R53"),
                String::from("PSCP"),
                String::from("RNF53"),
            ],
            added: NaiveDateTime::default(),
        };
        assert_eq!(expected, record);

        Ok(())
    }

    #[test]
    fn provider_get_pro_ac_for_tx_ac() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        assert!(provider.get_pro_ac_for_tx_ac("NM_007294.0").is_err());

        assert_eq!(
            provider.get_pro_ac_for_tx_ac("NM_007294.3")?,
            Some(String::from("NP_009225.1"))
        );

        Ok(())
    }

    #[test]
    fn provider_get_acs_for_protein_seq() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        assert_eq!(
            provider.get_acs_for_protein_seq("XXX")?,
            Vec::<String>::new()
        );

        Ok(())
    }

    #[test]
    fn provider_get_similar_transcripts() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        assert_eq!(
            provider.get_similar_transcripts("XXX")?,
            Vec::<TxSimilarityRecord>::new()
        );

        Ok(())
    }

    #[test]
    fn provider_get_tx_exons() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_exons("NM_007294.3", "NC_000017.10", "splign")?;

        let expected = vec![
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 22,
                tx_start_i: 5699,
                tx_end_i: 7207,
                alt_start_i: 41196311,
                alt_end_i: 41197819,
                cigar: String::from("1508M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 21,
                tx_start_i: 5638,
                tx_end_i: 5699,
                alt_start_i: 41199659,
                alt_end_i: 41199720,
                cigar: String::from("61M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 20,
                tx_start_i: 5564,
                tx_end_i: 5638,
                alt_start_i: 41201137,
                alt_end_i: 41201211,
                cigar: String::from("74M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 19,
                tx_start_i: 5509,
                tx_end_i: 5564,
                alt_start_i: 41203079,
                alt_end_i: 41203134,
                cigar: String::from("55M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 18,
                tx_start_i: 5425,
                tx_end_i: 5509,
                alt_start_i: 41209068,
                alt_end_i: 41209152,
                cigar: String::from("84M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 17,
                tx_start_i: 5384,
                tx_end_i: 5425,
                alt_start_i: 41215349,
                alt_end_i: 41215390,
                cigar: String::from("41M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 16,
                tx_start_i: 5306,
                tx_end_i: 5384,
                alt_start_i: 41215890,
                alt_end_i: 41215968,
                cigar: String::from("78M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 15,
                tx_start_i: 5218,
                tx_end_i: 5306,
                alt_start_i: 41219624,
                alt_end_i: 41219712,
                cigar: String::from("88M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 14,
                tx_start_i: 4907,
                tx_end_i: 5218,
                alt_start_i: 41222944,
                alt_end_i: 41223255,
                cigar: String::from("311M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 13,
                tx_start_i: 4716,
                tx_end_i: 4907,
                alt_start_i: 41226347,
                alt_end_i: 41226538,
                cigar: String::from("191M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 12,
                tx_start_i: 4589,
                tx_end_i: 4716,
                alt_start_i: 41228504,
                alt_end_i: 41228631,
                cigar: String::from("127M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 11,
                tx_start_i: 4417,
                tx_end_i: 4589,
                alt_start_i: 41234420,
                alt_end_i: 41234592,
                cigar: String::from("172M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 10,
                tx_start_i: 4328,
                tx_end_i: 4417,
                alt_start_i: 41242960,
                alt_end_i: 41243049,
                cigar: String::from("89M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 9,
                tx_start_i: 902,
                tx_end_i: 4328,
                alt_start_i: 41243451,
                alt_end_i: 41246877,
                cigar: String::from("3426M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 8,
                tx_start_i: 825,
                tx_end_i: 902,
                alt_start_i: 41247862,
                alt_end_i: 41247939,
                cigar: String::from("77M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 7,
                tx_start_i: 779,
                tx_end_i: 825,
                alt_start_i: 41249260,
                alt_end_i: 41249306,
                cigar: String::from("46M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 6,
                tx_start_i: 673,
                tx_end_i: 779,
                alt_start_i: 41251791,
                alt_end_i: 41251897,
                cigar: String::from("106M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 5,
                tx_start_i: 533,
                tx_end_i: 673,
                alt_start_i: 41256138,
                alt_end_i: 41256278,
                cigar: String::from("140M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 4,
                tx_start_i: 444,
                tx_end_i: 533,
                alt_start_i: 41256884,
                alt_end_i: 41256973,
                cigar: String::from("89M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 3,
                tx_start_i: 366,
                tx_end_i: 444,
                alt_start_i: 41258472,
                alt_end_i: 41258550,
                cigar: String::from("78M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 2,
                tx_start_i: 312,
                tx_end_i: 366,
                alt_start_i: 41267742,
                alt_end_i: 41267796,
                cigar: String::from("54M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 1,
                tx_start_i: 213,
                tx_end_i: 312,
                alt_start_i: 41276033,
                alt_end_i: 41276132,
                cigar: String::from("99M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
            TxExonsRecord {
                hgnc: String::from("BRCA1"),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
                alt_strand: -1,
                ord: 0,
                tx_start_i: 0,
                tx_end_i: 213,
                alt_start_i: 41277287,
                alt_end_i: 41277500,
                cigar: String::from("213M"),
                tx_aseq: None,
                alt_aseq: None,
                tx_exon_set_id: 2147483647,
                alt_exon_set_id: 2147483647,
                tx_exon_id: 2147483647,
                alt_exon_id: 2147483647,
                exon_aln_id: 2147483647,
            },
        ];

        assert_eq!(result, expected);

        Ok(())
    }

    #[test]
    fn provider_get_tx_for_gene() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_for_gene("BRCA1")?;

        let expected = vec![
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(232),
                cds_end_i: Some(5824),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(232),
                cds_end_i: Some(5824),
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(232),
                cds_end_i: Some(5887),
                tx_ac: String::from("NM_007300.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(232),
                cds_end_i: Some(5887),
                tx_ac: String::from("NM_007300.3"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(281),
                cds_end_i: Some(5732),
                tx_ac: String::from("NM_007297.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(281),
                cds_end_i: Some(5732),
                tx_ac: String::from("NM_007297.3"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(194),
                cds_end_i: Some(2294),
                tx_ac: String::from("NM_007299.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(194),
                cds_end_i: Some(2294),
                tx_ac: String::from("NM_007299.3"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(193),
                cds_end_i: Some(5785),
                tx_ac: String::from("XM_006722029.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(717),
                cds_end_i: Some(6309),
                tx_ac: String::from("XM_006722031.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(364),
                cds_end_i: Some(5746),
                tx_ac: String::from("XM_006722035.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(173),
                cds_end_i: Some(2453),
                tx_ac: String::from("XM_006722039.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(145),
                cds_end_i: Some(5734),
                tx_ac: String::from("XM_006722032.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(145),
                cds_end_i: Some(5734),
                tx_ac: String::from("XM_006722033.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(145),
                cds_end_i: Some(5734),
                tx_ac: String::from("XM_006722034.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(145),
                cds_end_i: Some(2428),
                tx_ac: String::from("XM_006722037.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(145),
                cds_end_i: Some(2425),
                tx_ac: String::from("XM_006722038.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(145),
                cds_end_i: Some(2425),
                tx_ac: String::from("XM_006722040.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(145),
                cds_end_i: Some(2302),
                tx_ac: String::from("XM_006722041.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(113),
                cds_end_i: Some(5705),
                tx_ac: String::from("NM_007294.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(113),
                cds_end_i: Some(5705),
                tx_ac: String::from("NM_007294.4"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(194),
                cds_end_i: Some(5645),
                tx_ac: String::from("NM_007297.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(194),
                cds_end_i: Some(5645),
                tx_ac: String::from("NM_007297.4"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(107),
                cds_end_i: Some(2207),
                tx_ac: String::from("NM_007299.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(107),
                cds_end_i: Some(2207),
                tx_ac: String::from("NM_007299.4"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(113),
                cds_end_i: Some(5768),
                tx_ac: String::from("NM_007300.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(113),
                cds_end_i: Some(5768),
                tx_ac: String::from("NM_007300.4"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: None,
                cds_end_i: None,
                tx_ac: String::from("NR_027676.2"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: None,
                cds_end_i: None,
                tx_ac: String::from("NR_027676.2"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: None,
                cds_end_i: None,
                tx_ac: String::from("NR_027676.1"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: None,
                cds_end_i: None,
                tx_ac: String::from("NR_027676.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(161),
                cds_end_i: Some(5753),
                tx_ac: String::from("XM_006722030.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(19),
                cds_end_i: Some(2299),
                tx_ac: String::from("NM_007298.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(19),
                cds_end_i: Some(2299),
                tx_ac: String::from("NM_007298.3"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
            TxInfoRecord {
                hgnc: String::from("BRCA1"),
                cds_start_i: Some(114),
                cds_end_i: Some(5325),
                tx_ac: String::from("XM_006722036.1"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
        ];

        assert_eq!(result, expected);

        Ok(())
    }

    #[test]
    fn provider_get_tx_for_region_empty() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        let result =
            provider.get_tx_for_region("NC_000017.10", "splign", 50_000_000, 50_000_001)?;

        let expected = vec![];

        assert_eq!(result, expected);

        Ok(())
    }

    #[test]
    fn provider_get_tx_for_region_brca1() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_for_region("NC_000017.10", "splign", 41196311, 41197819)?;

        let expected = vec![
            TxForRegionRecord {
                tx_ac: String::from("NM_007300.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277500,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277500,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007299.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277468,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007297.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277468,
            },
            TxForRegionRecord {
                tx_ac: String::from("NR_027676.2"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277381,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007300.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277381,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007299.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277381,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007297.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277381,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007294.4"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277381,
            },
            TxForRegionRecord {
                tx_ac: String::from("NR_027676.1"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41277340,
            },
            TxForRegionRecord {
                tx_ac: String::from("NM_007298.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_strand: -1,
                alt_aln_method: String::from("splign"),
                start_i: 41196311,
                end_i: 41276132,
            },
        ];

        assert_eq!(result, expected);

        Ok(())
    }

    #[test]
    fn provider_get_tx_info() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_info("NM_007294.3", "NC_000017.10", "splign")?;

        let expected = TxInfoRecord {
            hgnc: String::from("BRCA1"),
            cds_start_i: Some(232),
            cds_end_i: Some(5824),
            tx_ac: String::from("NM_007294.3"),
            alt_ac: String::from("NC_000017.10"),
            alt_aln_method: String::from("splign"),
        };

        assert_eq!(result, expected);

        Ok(())
    }

    #[test]
    fn provider_get_tx_mapping_options() -> Result<(), anyhow::Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_mapping_options("NM_007294.3")?;

        let expected = vec![
            TxMappingOptionsRecord {
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.10"),
                alt_aln_method: String::from("splign"),
            },
            TxMappingOptionsRecord {
                tx_ac: String::from("NM_007294.3"),
                alt_ac: String::from("NC_000017.11"),
                alt_aln_method: String::from("splign"),
            },
        ];

        assert_eq!(result, expected);

        Ok(())
    }

    fn build_mapper_37(normalize: bool) -> Result<Mapper, anyhow::Error> {
        let provider = Rc::new(build_provider()?);
        let config = AssemblyMapperConfig {
            assembly: Assembly::Grch37,
            normalize,
            ..Default::default()
        };
        Ok(Mapper::new(config, provider))
    }

    #[test]
    fn mapper_brca1_g_c() -> Result<(), anyhow::Error> {
        let mapper = build_mapper_37(false)?;
        let hgvs_g = "NC_000017.10:g.41197701G>C";
        let hgvs_c = "NM_007294.4:c.5586C>G";
        let hgvs_p = "NP_009225.1:p.His1862Gln";

        let var_g = HgvsVariant::from_str(hgvs_g)?;
        let var_c = mapper.g_to_c(&var_g, "NM_007294.4")?;
        let var_p = mapper.c_to_p(&var_c)?;

        assert_eq!(format!("{var_c}"), hgvs_c);
        assert_eq!(format!("{var_p}"), hgvs_p);

        Ok(())
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

// <LICENSE>
// MIT License
//
// Copyright (c) 2022 SACGF
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// </LICENSE>
