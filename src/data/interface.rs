//! Definition of the interface for accessing the transcript database.

use chrono::NaiveDateTime;
use linked_hash_map::LinkedHashMap;

use crate::static_data::Assembly;

/// Information about a gene.
///
/// ```text
/// hgnc    | ATM
/// maploc  | 11q22-q23
/// descr   | ataxia telangiectasia mutated
/// summary | The protein encoded by this gene belongs to the PI3/PI4-kinase family. This...
/// aliases | AT1,ATA,ATC,ATD,ATE,ATDC,TEL1,TELO1
/// added   | 2014-02-04 21:39:32.57125
/// ```
#[derive(Debug, PartialEq)]
pub struct GeneInfoRecord {
    pub hgnc: String,
    pub maploc: String,
    pub descr: String,
    pub summary: String,
    pub aliases: Vec<String>,
    pub added: NaiveDateTime,
}

/// Information about similar transcripts.
///
/// ```text
/// tx_ac1                 | NM_001285829.1
/// tx_ac2                 | ENST00000341255
/// hgnc_eq                | f
/// cds_eq                 | f
/// es_fp_eq               | f
/// cds_es_fp_eq           | f
/// cds_exon_lengths_fp_eq | t
/// ```
///
/// Hint: "es" = "exon set", "fp" = "fingerprint", "eq" = "equal"
///
/// "Exon structure" refers to the start and end coordinates on a
/// specified reference sequence. Thus, having the same exon
/// structure means that the transcripts are defined on the same
/// reference sequence and have the same exon spans on that
/// sequence.
#[derive(Debug, PartialEq)]
pub struct TxSimilarityRecord {
    /// Accession of first transcript.
    pub tx_ac1: String,
    /// Accession of second transcript.
    pub tx_ac2: String,
    pub hgnc_eq: bool,
    /// Whether CDS sequences are identical.
    pub cds_eq: bool,
    /// Whether the full exon structures are identical (i.e., incl. UTR).
    pub es_fp_eq: bool,
    /// Whether the cds-clipped portions of the exon structures are identical
    /// (i.e., ecluding. UTR).
    pub cds_es_fp_eq: bool,
    pub cds_exon_lengths_fp_eq: bool,
}

///```text
/// hgnc            | TGDS
/// tx_ac           | NM_001304430.1
/// alt_ac          | NC_000013.10
/// alt_aln_method  | blat
/// alt_strand      | -1
/// ord             | 0
/// tx_start_i      | 0
/// tx_end_i        | 301
/// alt_start_i     | 95248228
/// alt_end_i       | 95248529
/// cigar           | 301=
/// tx_aseq         |
/// alt_aseq        |
/// tx_exon_set_id  | 348239
/// alt_exon_set_id | 722624
/// tx_exon_id      | 3518579
/// alt_exon_id     | 6063334
/// exon_aln_id     | 3461425
///```
#[derive(Debug, PartialEq)]
pub struct TxExonsRecord {
    pub hgnc: String,
    pub tx_ac: String,
    pub alt_ac: String,
    pub alt_aln_method: String,
    pub alt_strand: i16,
    pub ord: i32,
    pub tx_start_i: i32,
    pub tx_end_i: i32,
    pub alt_start_i: i32,
    pub alt_end_i: i32,
    pub cigar: String,
    pub tx_aseq: Option<String>,
    pub alt_aseq: Option<String>,
    pub tx_exon_set_id: i32,
    pub alt_exon_set_id: i32,
    pub tx_exon_id: i32,
    pub alt_exon_id: i32,
    pub exon_aln_id: i32,
}

/// ```text
/// tx_ac          | NM_001304430.2
/// alt_ac         | NC_000013.10
/// alt_strand     | -1
/// alt_aln_method | splign
/// start_i        | 95226307
/// end_i          | 95248406
/// ```
#[derive(Debug, PartialEq)]
pub struct TxForRegionRecord {
    pub tx_ac: String,
    pub alt_ac: String,
    pub alt_strand: i16,
    pub alt_aln_method: String,
    pub start_i: i32,
    pub end_i: i32,
}

/// ```text
/// tx_ac          | NM_199425.2
/// alt_ac         | NM_199425.2
/// alt_aln_method | transcript
/// cds_start_i    | 283
/// cds_end_i      | 1003
/// lengths        | {707,79,410}
/// hgnc           | VSX1
/// ```
#[derive(Debug, PartialEq)]
pub struct TxIdentityInfo {
    pub tx_ac: String,
    pub alt_ac: String,
    pub alt_aln_method: String,
    pub cds_start_i: i32,
    pub cds_end_i: i32,
    pub lengths: Vec<i32>,
    pub hgnc: String,
}

/// ```text
/// hgnc           | ATM
/// cds_start_i    | 385
/// cds_end_i      | 9556
/// tx_ac          | NM_000051.3
/// alt_ac         | AC_000143.1
/// alt_aln_method | splign
/// ```
#[derive(Debug, PartialEq)]
pub struct TxInfoRecord {
    pub hgnc: String,
    pub cds_start_i: Option<i32>,
    pub cds_end_i: Option<i32>,
    pub tx_ac: String,
    pub alt_ac: String,
    pub alt_aln_method: String,
}

/// ```text
/// -[ RECORD 1 ]--+----------------
/// tx_ac          | ENST00000000233
/// alt_ac         | NC_000007.13
/// alt_aln_method | genebuild
/// -[ RECORD 2 ]--+----------------
/// tx_ac          | ENST00000000412
/// alt_ac         | NC_000012.11
/// alt_aln_method | genebuild
/// ```
#[derive(Debug, PartialEq)]
pub struct TxMappingOptionsRecord {
    pub tx_ac: String,
    pub alt_ac: String,
    pub alt_aln_method: String,
}

/// Interface for data providers.
pub trait Provider {
    /// Return the data version, e.g., `uta_20180821`.
    fn data_version(&self) -> &str;

    /// Return the schema version, e.g., `"1.1"`.
    fn schema_version(&self) -> &str;

    /// Return a map from accession to chromosome name for the given assembly
    ///
    /// For example, when `assembly_name = "GRCh38.p5"`, the value for `"NC_000001.11"`
    /// would be `"1"`.
    ///
    /// # Arguments
    ///
    /// * `assembly` - The assembly to build the map for.
    fn get_assembly_map(&self, assembly: Assembly) -> LinkedHashMap<String, String>;

    /// Returns the basic information about the gene.
    ///
    /// # Arguments
    ///
    /// * `hgnc` - HGNC gene name
    fn get_gene_info(&self, hgnc: &str) -> Result<GeneInfoRecord, anyhow::Error>;

    /// Return the (single) associated protein accession for a given transcript accession,
    /// or None if not found.
    ///
    /// # Arguments
    ///
    /// * `tx_ac` -- transcript accession with version (e.g., 'NM_000051.3')
    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, anyhow::Error>;

    /// Return full sequence for the given accession.
    ///
    /// # Arguments
    ///
    /// * `ac` -- accession
    fn get_seq(&self, ac: &str) -> Result<String, anyhow::Error>;

    /// Return sequence part for the given accession.
    ///
    /// # Arguments
    ///
    /// * `ac` -- accession
    /// * `start` -- start position (0-based, start of sequence if missing)
    /// * `end` -- end position (0-based, end of sequence if missing)
    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, anyhow::Error>;

    /// Return a list of transcripts that are similar to the given transcript, with relevant
    /// similarity criteria.
    ///
    /// # Arguments
    ///
    /// * `tx_ac` -- transcript accession with version (e.g., 'NM_000051.3')
    fn get_similar_transcripts(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<TxSimilarityRecord>, anyhow::Error>;

    /// Return transcript exon info for supplied accession (tx_ac, alt_ac, alt_aln_method),
    /// or empty `Vec` if not found.
    ///
    /// # Arguments
    ///
    /// * `tx_ac` -- transcript accession with version (e.g., 'NM_000051.3')
    /// * `alt_ac` -- specific genomic sequence (e.g., NC_000011.4)
    /// * `alt_aln_method` -- sequence alignment method (e.g., splign, blat)
    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, anyhow::Error>;

    /// Return transcript info records for supplied gene, in order of decreasing length.
    ///
    /// # Arguments
    ///
    /// * `gene` - HGNC gene name
    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, anyhow::Error>;

    /// Return transcripts that overlap given region.
    ///
    /// # Arguments
    ///
    // * `alt_ac` -- reference sequence (e.g., NC_000007.13)
    // * `alt_aln_method` -- alignment method (e.g., splign)
    // * `start_i` -- 5' bound of region
    // * `end_i` -- 3' bound of region
    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, anyhow::Error>;

    /// Return features associated with a single transcript.
    ///
    /// # Arguments
    ///
    /// * `tx_ac` -- transcript accession with version (e.g., 'NM_199425.2')
    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, anyhow::Error>;

    /// Return a single transcript info for supplied accession (tx_ac, alt_ac, alt_aln_method), or None if not found.
    ///
    /// # Arguments
    ///
    /// * `tx_ac` -- transcript accession with version (e.g., 'NM_000051.3')
    /// * `alt_ac -- specific genomic sequence (e.g., NC_000011.4)
    /// * `alt_aln_method` -- sequence alignment method (e.g., splign, blat)
    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, anyhow::Error>;

    /// Return all transcript alignment sets for a given transcript accession (tx_ac).
    ///
    /// Returns empty list if transcript does not exist.  Use this method to discovery
    /// possible mapping options supported in the database.
    ///
    /// # Arguments
    ///
    /// * `tx_ac` -- transcript accession with version (e.g., 'NM_000051.3')
    fn get_tx_mapping_options(
        &self,
        tax_ac: &str,
    ) -> Result<Vec<TxMappingOptionsRecord>, anyhow::Error>;
}
