//! Access to `cdot` transcripts from local JSON files with sequences from a `seqrepo`.
//!
//! https://github.com/SACGF/cdot

use std::{collections::HashMap, path::PathBuf, sync::Arc, time::Instant};

use crate::{data::error::Error, data::interface::{
    self, GeneInfoRecord, TxExonsRecord, TxForRegionRecord, TxIdentityInfo, TxInfoRecord,
    TxMappingOptionsRecord, TxSimilarityRecord,
}, sequences::TranslationTable, Sequence};
use biocommons_bioutils::assemblies::{Assembly, ASSEMBLY_INFOS};

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use chrono::NaiveDateTime;
use indexmap::IndexMap;
use seqrepo::{self, SeqRepo};

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
/// and `i32::MAX` for `tx_exon-set_id`, `alt_exon_set_id`, `tx_exon_id`, `alt_exon_id`,
/// `exon_aln_id`.
pub struct Provider {
    inner: TxProvider,
    seqrepo: Arc<dyn seqrepo::Interface + Sync + Send>,
}

impl Provider {
    pub fn new(config: Config) -> Result<Self, Error> {
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
            inner: TxProvider::with_config(
                config
                    .json_paths
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<&str>>()
                    .as_ref(),
            )?,
            seqrepo: Arc::new(SeqRepo::new(path, &instance)?),
        })
    }

    /// Create a new provider allowing to inject a seqrepo.
    pub fn with_seqrepo(
        config: Config,
        seqrepo: Arc<dyn seqrepo::Interface + Sync + Send>,
    ) -> Result<Provider, Error> {
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

impl interface::Provider for Provider {
    fn data_version(&self) -> &str {
        self.inner.data_version()
    }

    fn schema_version(&self) -> &str {
        self.inner.schema_version()
    }

    fn get_assembly_map(
        &self,
        assembly: biocommons_bioutils::assemblies::Assembly,
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
    ) -> Result<Sequence, Error> {
        self.seqrepo
            .fetch_sequence_part(
                &seqrepo::AliasOrSeqId::Alias {
                    value: ac.to_string(),
                    namespace: None,
                },
                begin,
                end,
            )
            .map_err(Error::SeqRepoError)
    }

    fn get_acs_for_protein_seq(&self, seq: &[u8]) -> Result<Vec<String>, Error> {
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

/// Data structures used for deserializing from cdot.
pub mod models {
    use indexmap::IndexMap;
    use serde::{Deserialize, Deserializer, Serialize};

    /// Container for a cDot data file.
    #[derive(Deserialize, Serialize, Debug, Clone)]
    pub struct Container {
        pub transcripts: IndexMap<String, Transcript>,
        pub cdot_version: String,
        pub genome_builds: Vec<String>,
        pub genes: IndexMap<String, Gene>,
    }

    /// Enum for representing the tags for transcripts.
    #[derive(Deserialize, Serialize, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
    pub enum Tag {
        Basic,
        EnsemblCanonical,
        ManeSelect,
        ManePlusClinical,
        RefSeqSelect,
        GencodePrimary,
        Other,
    }

    #[derive(Deserialize, Serialize, Debug, Clone)]
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
        pub genome_builds: IndexMap<String, GenomeAlignment>,
        /// HGNC identifier, e.g. `"1100"` for BRCA1 which is `HGNC:1100`.
        #[serde(default)]
        pub hgnc: Option<String>,
        /// Identifier of the transcript, same as key in `transcripts`, e.g., `"NM_007294.3"` for BRCA1.
        pub id: String,
        /// When the transcript does not match the primary assembly (`"NC_"`) perfectly / the assembly itself is wrong/partial
        #[serde(skip_serializing_if = "Option::is_none")]
        #[serde(default)]
        pub partial: Option<u8>,
        /// Identifier of corresponding protein, e.g., `"NP_009225.1"` for `"NM_007294.3"` of BRCA1.
        #[serde(default)]
        pub protein: Option<String>,
        /// 0-based position of codon of transcript, e.g., `232` for `"NM_007294.3"` of BRCA1.
        #[serde(default)]
        pub start_codon: Option<i32>,
        /// End position of stop codon of transcript, e.g., `5824` for `"NM_007294.3"` of BRCA1.
        #[serde(default)]
        pub stop_codon: Option<i32>,
    }

    /// Representation of the strand.
    #[derive(Deserialize, Serialize, Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Strand {
        #[serde(rename = "+")]
        Plus,
        #[serde(rename = "-")]
        Minus,
    }

    /// Representation of an exon in the JSON.
    #[derive(Deserialize, Serialize, Debug, Clone)]
    struct ExonHelper(i32, i32, i32, i32, i32, Option<String>);

    /// Representation of an exon after loading.
    #[derive(Serialize, Debug, Clone)]
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
    #[derive(Deserialize, Serialize, Debug, Clone)]
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
        /// Tags of the transcript.
        #[serde(default)]
        #[serde(deserialize_with = "deserialize_tag")]
        pub tag: Option<Vec<Tag>>,
        /// Any notes for the transcript.
        #[serde(default)]
        pub note: Option<String>,
    }

    /// Enum for representing the biotypes.
    #[derive(Deserialize, Serialize, Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum BioType {
        #[serde(rename = "3prime_overlapping_ncrna")]
        ThreePrimeOverlappingNcRna,
        #[serde(rename = "aberrant_processed_transcript")]
        AberrantProcessedTranscript,
        #[serde(rename = "antisense")]
        Antisense,
        #[serde(rename = "antisense_RNA")]
        AntisenseRna,
        #[serde(rename = "artifact")]
        Artifact,
        #[serde(rename = "C_gene_segment")]
        CGeneSegment,
        #[serde(rename = "C_region")]
        CRegion,
        #[serde(rename = "C_region_pseudogene")]
        CRegionPseudogene,
        #[serde(rename = "D_gene_segment")]
        DGeneSegment,
        #[serde(rename = "D_segment")]
        DSegment,
        #[serde(rename = "guide_RNA")]
        GuideRna,
        #[serde(rename = "IG_C_gene")]
        IgCGene,
        #[serde(rename = "IG_C_pseudogene")]
        IgCPseudogene,
        #[serde(rename = "IG_D_gene")]
        IgDGene,
        #[serde(rename = "IG_J_gene")]
        IgJGene,
        #[serde(rename = "IG_J_pseudogene")]
        IgJPseudogene,
        #[serde(rename = "IG_pseudogene")]
        IgPseudogene,
        #[serde(rename = "IG_V_gene")]
        IgVGene,
        #[serde(rename = "IG_V_pseudogene")]
        IgVPseudogene,
        #[serde(rename = "J_gene_segment")]
        JGeneSegment,
        #[serde(rename = "J_segment")]
        JSegment,
        #[serde(rename = "J_segment_pseudogene")]
        JSegmentPseudogene,
        #[serde(rename = "lincRNA")]
        LincRna,
        #[serde(rename = "lnc_RNA", alias = "lncRNA")]
        LncRna,
        #[serde(rename = "miRNA")]
        MiRna,
        #[serde(rename = "misc_RNA")]
        MiscRna,
        #[serde(rename = "mRNA")]
        MRna,
        #[serde(rename = "Mt_rRNA")]
        MtRRna,
        #[serde(rename = "Mt_tRNA")]
        MtTRna,
        #[serde(rename = "nc_primary_transcript")]
        NcPrimaryTranscript,
        #[serde(rename = "ncRNA")]
        NcRna,
        #[serde(rename = "ncRNA_pseudogene")]
        NcRnaPseudogene,
        #[serde(rename = "NMD_transcript_variant")]
        NmdTranscriptVariant,
        #[serde(rename = "non_coding")]
        NonCoding,
        #[serde(rename = "other")]
        Other,
        #[serde(rename = "polymorphic_pseudogene")]
        PolymorphicPseudogene,
        #[serde(rename = "primary_transcript")]
        PrimaryTranscript,
        #[serde(rename = "processed_pseudogene")]
        ProcessedPseudogene,
        #[serde(rename = "processed_transcript")]
        ProcessedTranscript,
        #[serde(rename = "protein_coding")]
        ProteinCoding,
        #[serde(rename = "pseudogenic_transcript")]
        PseudogenicTranscript,
        #[serde(rename = "pseudogene")]
        Pseudogene,
        #[serde(rename = "RNase_MRP_RNA")]
        RnaseMrpRna,
        #[serde(rename = "RNase_P_RNA")]
        RnasePRna,
        #[serde(rename = "scaRNA")]
        ScnaRna,
        #[serde(rename = "scRNA")]
        ScnRna,
        #[serde(rename = "ribozyme")]
        Ribozyme,
        #[serde(rename = "rRNA")]
        RRna,
        #[serde(rename = "rRNA_pseudogene")]
        RnaPseudogene,
        #[serde(rename = "selenoprotein")]
        Selenoprotein,
        #[serde(rename = "sense_intronic")]
        SenseIntronic,
        #[serde(rename = "sense_overlapping")]
        SenseOverlapping,
        #[serde(rename = "snoRNA")]
        SnoRna,
        #[serde(rename = "snRNA")]
        SnRna,
        #[serde(rename = "sRNA")]
        SRna,
        #[serde(rename = "TEC")]
        Tec,
        #[serde(rename = "tRNA")]
        TRna,
        #[serde(rename = "telomerase_RNA")]
        TelomeraseRna,
        #[serde(rename = "transcribed_pseudogene")]
        TranscribedPseudogene,
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
        #[serde(rename = "unconfirmed_transcript")]
        UnconfirmedTranscript,
        #[serde(rename = "unprocessed_pseudogene")]
        UnprocessedPseudogene,
        #[serde(rename = "vault_RNA")]
        VaultRna,
        #[serde(rename = "vaultRNA_primary_transcript")]
        VaultRnaPrimaryTranscript,
        #[serde(rename = "V_segment")]
        VSegment,
        #[serde(rename = "V_gene_segment")]
        VGeneSegment,
        #[serde(rename = "V_segment_pseudogene")]
        VSegmentPseudogene,
        #[serde(rename = "Y_RNA")]
        YRna,
    }

    #[derive(Deserialize, Serialize, Debug, Clone)]
    pub struct Gene {
        #[serde(default)]
        #[serde(deserialize_with = "deserialize_gene_aliases")]
        pub aliases: Option<Vec<String>>,
        #[serde(default)]
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

    pub fn str_to_tag(s: &str) -> Tag {
        match s {
            "basic" => Tag::Basic,
            "Ensembl_canonical" => Tag::EnsemblCanonical,
            "MANE_Plus_Clinical" | "MANE Plus Clinical" => Tag::ManePlusClinical,
            "MANE_Select" | "MANE Select" => Tag::ManeSelect,
            "RefSeq Select" => Tag::RefSeqSelect,
            "GENCODE Primary" => Tag::GencodePrimary,
            _ => {
                log::trace!("unknown tag: {}", s);
                Tag::Other
            }
        }
    }

    #[derive(Debug, Deserialize, Serialize)]
    struct WrappedString(String);

    fn deserialize_tag<'de, D>(deserializer: D) -> Result<Option<Vec<Tag>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Option::<WrappedString>::deserialize(deserializer).map(|opt_wrapped| {
            opt_wrapped.map(|wrapped| wrapped.0.split(',').map(str_to_tag).collect())
        })
    }

    fn deserialize_gene_aliases<'de, D>(deserializer: D) -> Result<Option<Vec<String>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let buf = Option::<String>::deserialize(deserializer)?;

        Ok(buf.map(|s| s.split(", ").map(|s| s.to_string()).collect()))
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
    fn with_config(json_paths: &[&str]) -> Result<Self, Error> {
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
    ) -> Result<(), Error> {
        log::debug!("Loading cdot transcripts from {:?}", json_path);
        let start = Instant::now();
        let models::Container {
            genes: c_genes,
            transcripts: c_txs,
            ..
        } = if json_path.ends_with(".gz") {
            serde_json::from_reader(flate2::bufread::GzDecoder::new(std::io::BufReader::new(
                std::fs::File::open(json_path)
                    .map_err(|_e| Error::CdotJsonOpen(json_path.to_string()))?,
            )))
            .map_err(|_e| Error::CdotJsonParse(json_path.to_string()))?
        } else {
            serde_json::from_reader(std::io::BufReader::new(
                std::fs::File::open(json_path)
                    .map_err(|_e| Error::CdotJsonOpen(json_path.to_string()))?,
            ))
            .map_err(|_e| Error::CdotJsonParse(json_path.to_string()))?
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
            .filter(|gene| {
                gene.gene_symbol.is_some()
                    && !gene
                        .gene_symbol
                        .as_ref()
                        .expect("should not happen; empty gene symbol filtered out")
                        .is_empty()
                    && gene.map_location.is_some()
                    && !gene
                        .map_location
                        .as_ref()
                        .expect("should not happen; empty location filtered out")
                        .is_empty()
            })
            .for_each(|gene| {
                let gene_symbol = gene
                    .gene_symbol
                    .as_ref()
                    .expect("should not happen; empty gene symbol filtered out")
                    .clone();
                transcript_ids_for_gene
                    .entry(gene_symbol.clone())
                    .or_default();
                genes.insert(gene_symbol, gene.clone());
            });
        c_txs
            .values()
            .filter(|tx| {
                tx.gene_name.is_some()
                    && !tx
                        .gene_name
                        .as_ref()
                        .expect("should not happen; empty gene name filtered out")
                        .is_empty()
            })
            .for_each(|tx| {
                let gene_name = tx
                    .gene_name
                    .as_ref()
                    .expect("should not happen; empty gene name filtered out");
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
                let tree = result
                    .get_mut(contig)
                    .expect("should not happen; just inserted");
                let alt_start_i = genome_alignment
                    .exons
                    .iter()
                    .map(|exon| exon.alt_start_i)
                    .min()
                    .expect("should not happen; must have at least one exon")
                    - 1;
                let alt_end_i = genome_alignment
                    .exons
                    .iter()
                    .map(|exon| exon.alt_end_i)
                    .max()
                    .expect("should not happen; must have at least one exon");
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

/// Implementation for `interface::Provider`.
impl TxProvider {
    fn data_version(&self) -> &str {
        REQUIRED_VERSION
    }

    fn schema_version(&self) -> &str {
        REQUIRED_VERSION
    }

    fn get_assembly_map(&self, assembly: Assembly) -> IndexMap<String, String> {
        IndexMap::from_iter(
            ASSEMBLY_INFOS[assembly]
                .sequences
                .iter()
                .map(|record| (record.refseq_ac.clone(), record.name.clone())),
        )
    }

    fn get_gene_info(&self, hgnc: &str) -> Result<GeneInfoRecord, Error> {
        let gene = self
            .genes
            .get(hgnc)
            .ok_or(Error::NoGeneFound(hgnc.to_string()))?;

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
            descr: gene.description.clone().unwrap_or_default(),
            summary: gene.summary.clone().unwrap_or_default(),
            aliases: gene.aliases.clone().unwrap_or_default(),
            added: NaiveDateTime::default(),
        })
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, Error> {
        let transcript = self
            .transcripts
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;
        Ok(transcript.protein.clone())
    }

    /// Note from the original cdot Python code.
    ///
    /// This is not implemented. The only caller has comment: 'TODO: drop get_acs_for_protein_seq'
    /// And is only ever called as a backup when get_pro_ac_for_tx_ac fails
    #[allow(dead_code)]
    fn get_acs_for_protein_seq(&self, _seq: &[u8]) -> Result<Vec<String>, Error> {
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
    fn get_similar_transcripts(&self, _tx_ac: &str) -> Result<Vec<TxSimilarityRecord>, Error> {
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
    ) -> Result<Vec<TxExonsRecord>, Error> {
        let tx = self
            .transcripts
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;

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
                    tx_exon_set_id: i32::MAX,
                    alt_exon_set_id: i32::MAX,
                    tx_exon_id: i32::MAX,
                    alt_exon_id: i32::MAX,
                    exon_aln_id: i32::MAX,
                })
                .collect())
        } else {
            Err(Error::NoAlignmentFound(
                tx_ac.to_string(),
                alt_ac.to_string(),
            ))
        }
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, Error> {
        if let Some(tx_acs) = self.transcript_ids_for_gene.get(gene) {
            let mut tmp = Vec::new();
            for tx_ac in tx_acs {
                let tx = self
                    .transcripts
                    .get(tx_ac)
                    .expect("should not happen by construction");
                let cds_start_i = tx.start_codon;
                let cds_end_i = tx.stop_codon;
                for genome_alignment in tx.genome_builds.values() {
                    let contig = &genome_alignment.contig;
                    let tx_start = genome_alignment
                        .exons
                        .first()
                        .expect("should not happen; must have at least one exon")
                        .alt_start_i;
                    let tx_end = genome_alignment
                        .exons
                        .last()
                        .expect("should not happen; must have at least one exon")
                        .alt_end_i;
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
    ) -> Result<Vec<TxForRegionRecord>, Error> {
        if alt_aln_method != NCBI_ALN_METHOD {
            return Ok(Vec::new());
        }

        if let Some(contig_itv_tree) = self.interval_trees.get(alt_ac) {
            let mut tmp = Vec::new();

            for entry in contig_itv_tree.find((start_i - 1)..(end_i)) {
                let tx_ac = entry.data();
                let tx = self
                    .transcripts
                    .get(tx_ac)
                    .expect("should not happen by construction");
                let genome_alignment = tx
                    .genome_builds
                    .values()
                    .find(|genome_alignment| genome_alignment.contig == alt_ac);
                if let Some(genome_alignment) = genome_alignment {
                    let tx_start = genome_alignment
                        .exons
                        .first()
                        .expect("should not happen; must have at least one exon")
                        .alt_start_i;
                    let tx_end = genome_alignment
                        .exons
                        .last()
                        .expect("should not happen; must have at least one exon")
                        .alt_end_i;
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

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, Error> {
        let tx = self
            .transcripts
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;

        let needle = "UGA stop codon recoded as selenocysteine";
        let is_selenoprotein = tx.genome_builds.iter().any(|(_, genome_alignment)| {
            genome_alignment
                .note
                .clone()
                .unwrap_or_default()
                .contains(needle)
        });
        let is_selenoprotein = is_selenoprotein
            || tx.biotype.as_ref().map_or(false, |bt| {
                bt.contains(&crate::data::cdot::json::models::BioType::Selenoprotein)
            });

        let hgnc = tx
            .gene_name
            .as_ref()
            .expect("cannot happen; transcripts without gene_name not imported")
            .clone();

        let mut tmp = tx
            .genome_builds
            .values()
            .next()
            .expect("should not happen; must have at least one alignment remaining")
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
            translation_table: if is_selenoprotein {
                TranslationTable::Selenocysteine
            } else {
                TranslationTable::Standard
            },
        })
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, Error> {
        let tx = self
            .transcripts
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;

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

    fn get_tx_mapping_options(&self, tx_ac: &str) -> Result<Vec<TxMappingOptionsRecord>, Error> {
        let tx = self
            .transcripts
            .get(tx_ac)
            .ok_or(Error::NoTranscriptFound(tx_ac.to_string()))?;

        Ok(tx
            .genome_builds
            .values()
            .map(|genome_alignment| TxMappingOptionsRecord {
                tx_ac: tx_ac.to_string(),
                alt_ac: genome_alignment.contig.clone(),
                alt_aln_method: NCBI_ALN_METHOD.to_string(),
            })
            .collect())
    }
}

#[cfg(test)]
pub mod test_helpers {
    use anyhow::Error;
    use std::sync::Arc;

    use crate::data::uta_sr::test_helpers::build_writing_sr;

    use super::{Config, Provider};
    use seqrepo;

    pub fn build_provider() -> Result<Provider, Error> {
        let sr_cache_mode = std::env::var("TEST_SEQREPO_CACHE_MODE")
            .expect("Environment variable TEST_SEQREPO_CACHE_MODE undefined!");
        let sr_cache_path = std::env::var("TEST_SEQREPO_CACHE_PATH")
            .expect("Environment variable TEST_SEQREPO_CACHE_PATH undefined!");

        log::debug!("building provider...");
        let seqrepo = if sr_cache_mode == "read" {
            log::debug!("reading provider...");
            let seqrepo: Arc<dyn seqrepo::Interface + Send + Sync> =
                Arc::new(seqrepo::CacheReadingSeqRepo::new(sr_cache_path)?);
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
                "tests/data/data/cdot/cdot-0.2.21.refseq.grch37_grch38.brca1.json",
            )],
            seqrepo_path: String::from("nonexisting"),
        };
        // Note that we don't actually use the seqrepo instance except for initializing the provider.
        Provider::with_seqrepo(config, seqrepo).map_err(|e| anyhow::anyhow!(e))
    }
}

#[cfg(test)]
pub mod tests {
    use anyhow::Error;

    use std::str::FromStr;
    use std::sync::Arc;

    use pretty_assertions::assert_eq;
    use test_log::test;

    use super::models::{gap_to_cigar, Container};
    use super::test_helpers::build_provider;
    use crate::data::interface::{Provider, TxSimilarityRecord};
    use crate::mapper::assembly::{self, Mapper};
    use crate::parser::HgvsVariant;
    use biocommons_bioutils::assemblies::Assembly;

    #[test]
    fn test_sync() {
        fn is_sync<T: Sync>() {}
        is_sync::<super::Provider>();
    }

    #[test]
    fn deserialize_brca1() -> Result<(), Error> {
        let json = std::fs::read_to_string(
            "tests/data/data/cdot/cdot-0.2.21.refseq.grch37_grch38.brca1.json",
        )?;

        let c: Container = serde_json::from_str(&json)?;
        assert_eq!(c.cdot_version, "0.2.21");

        insta::assert_debug_snapshot!(&c);

        Ok(())
    }

    #[test]
    fn gap_to_cigar_smoke() {
        assert_eq!(gap_to_cigar("M196 I1 M61 I1 M181"), "196=1D61=1D181=");
    }

    /// Deserialization of the big cdot files for benchmarking.
    #[cfg(deserialization_tests)]
    #[test]
    fn deserialize_big_files() -> Result<(), Error> {
        let before = std::time::Instant::now();
        println!("ensembl...");
        let _ensembl: Container = serde_json::from_reader(flate2::bufread::GzDecoder::new(
            std::io::BufReader::new(std::fs::File::open("cdot-0.2.21.ensembl.grch37.json.gz")?),
        ))?;
        println!("Reading ensembl: {:?}", before.elapsed());

        let before = std::time::Instant::now();
        println!("refseq...");
        let _refseq: Container = serde_json::from_reader(flate2::bufread::GzDecoder::new(
            std::io::BufReader::new(std::fs::File::open("cdot-0.2.21.refseq.grch37.json.gz")?),
        ))?;
        println!("Reading refseq: {:?}", before.elapsed());

        Ok(())
    }

    #[test]
    fn provider_brca1_smoke() -> Result<(), Error> {
        build_provider()?;
        Ok(())
    }

    #[test]
    fn provider_versions() -> Result<(), Error> {
        let provider = build_provider()?;

        assert_eq!(provider.data_version(), "1.1");
        assert_eq!(provider.schema_version(), "1.1");

        Ok(())
    }

    #[test]
    fn provider_get_assembly_map() -> Result<(), Error> {
        let provider = build_provider()?;
        assert_eq!(provider.get_assembly_map(Assembly::Grch37p10).len(), 275);
        assert_eq!(provider.get_assembly_map(Assembly::Grch38).len(), 455);

        Ok(())
    }

    #[test]
    fn provider_get_gene_info() -> Result<(), Error> {
        let provider = build_provider()?;

        assert!(provider.get_gene_info("BRCA2").is_err());

        let record = provider.get_gene_info("BRCA1")?;

        insta::assert_debug_snapshot!(&record);

        Ok(())
    }

    #[test]
    fn provider_get_pro_ac_for_tx_ac() -> Result<(), Error> {
        let provider = build_provider()?;

        assert!(provider.get_pro_ac_for_tx_ac("NM_007294.0").is_err());

        assert_eq!(
            provider.get_pro_ac_for_tx_ac("NM_007294.3")?,
            Some(String::from("NP_009225.1"))
        );

        Ok(())
    }

    #[test]
    fn provider_get_acs_for_protein_seq() -> Result<(), Error> {
        let provider = build_provider()?;

        assert_eq!(
            provider.get_acs_for_protein_seq(b"XXX")?,
            Vec::<String>::new()
        );

        Ok(())
    }

    #[test]
    fn provider_get_similar_transcripts() -> Result<(), Error> {
        let provider = build_provider()?;

        assert_eq!(
            provider.get_similar_transcripts("XXX")?,
            Vec::<TxSimilarityRecord>::new()
        );

        Ok(())
    }

    #[test]
    fn provider_get_tx_exons() -> Result<(), Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_exons("NM_007294.3", "NC_000017.10", "splign")?;

        insta::assert_debug_snapshot!(&result);

        Ok(())
    }

    #[test]
    fn provider_get_tx_for_gene() -> Result<(), Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_for_gene("BRCA1")?;

        insta::assert_debug_snapshot!(&result);

        Ok(())
    }

    #[test]
    fn provider_get_tx_for_region_empty() -> Result<(), Error> {
        let provider = build_provider()?;

        let result =
            provider.get_tx_for_region("NC_000017.10", "splign", 50_000_000, 50_000_001)?;

        insta::assert_debug_snapshot!(&result);

        Ok(())
    }

    #[test]
    fn provider_get_tx_for_region_brca1() -> Result<(), Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_for_region("NC_000017.10", "splign", 41196311, 41197819)?;
        insta::assert_debug_snapshot!(&result);

        Ok(())
    }

    #[test]
    fn provider_get_tx_info() -> Result<(), Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_info("NM_007294.3", "NC_000017.10", "splign")?;

        insta::assert_debug_snapshot!(&result);

        Ok(())
    }

    #[test]
    fn provider_get_tx_mapping_options() -> Result<(), Error> {
        let provider = build_provider()?;

        let result = provider.get_tx_mapping_options("NM_007294.3")?;

        insta::assert_debug_snapshot!(&result);

        Ok(())
    }

    fn build_mapper_37(normalize: bool) -> Result<Mapper, Error> {
        let provider = Arc::new(build_provider()?);
        let config = assembly::Config {
            assembly: Assembly::Grch37,
            normalize,
            ..Default::default()
        };
        Ok(Mapper::new(config, provider))
    }

    #[test]
    fn mapper_brca1_g_c() -> Result<(), Error> {
        let mapper = build_mapper_37(false)?;
        let hgvs_g = "NC_000017.10:g.41197701G>C";
        let hgvs_n = "NM_007294.4:n.5699C>G";
        let hgvs_c = "NM_007294.4:c.5586C>G";
        let hgvs_p = "NP_009225.1:p.His1862Gln";

        let var_g = HgvsVariant::from_str(hgvs_g)?;
        let var_n = mapper.g_to_n(&var_g, "NM_007294.4")?;
        let var_c = mapper.g_to_c(&var_g, "NM_007294.4")?;
        let var_p = mapper.c_to_p(&var_c)?;

        assert_eq!(format!("{var_n}"), hgvs_n);
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
