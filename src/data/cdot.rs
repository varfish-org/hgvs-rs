//! Access to `cdot` transcripts with sequences from a `seqrepo`.
//!
//! https://github.com/SACGF/cdot

use std::rc::Rc;

use seqrepo::Interface as SeqRepoInterface;

/// Configurationf or the `data::cdot::Provider`.
#[derive(Debug, PartialEq, Clone)]
pub struct Config {
    /// Paths to the gzip-ed JSON files to load
    pub json_paths: Vec<String>,
    /// Path to the seqrepo directory, e.g., `/usr/local/share/seqrepo/latest`.  The last path
    /// component is the "instance" name.
    pub seqrepo_path: String,
}

/// This provider provides information from a UTA and a SeqRepo.
///
/// Transcripts from a UTA Postgres database, sequences comes from a SeqRepo.  This makes
/// genome contig information available in contrast to `data::uta::Provider`.
pub struct Provider {
    inner: TxProvider,
    seqrepo: Rc<dyn SeqRepoInterface>,
}

/// Data structures used for deserializing from cdot.
pub mod models {
    use linked_hash_map::LinkedHashMap;
    use serde::{Deserialize, Deserializer};

    /// Container for a cDot data file.
    #[derive(Deserialize, Debug)]
    pub struct Container {
        pub transcripts: LinkedHashMap<String, Transcript>,
        pub cdot_version: String,
        pub genome_builds: Vec<String>,
        pub genes: LinkedHashMap<String, Gene>,
    }

    /// Enum for representing the tags for transcripts.
    #[derive(Debug)]
    pub enum Tag {
        Basic,
        EnsemblCanonical,
        ManeSelect,
        ManePlusClinical,
        RefSeqSelect,
    }

    #[derive(Deserialize, Debug)]
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
    #[derive(Deserialize, Debug)]
    pub enum Strand {
        #[serde(rename = "+")]
        Plus,
        #[serde(rename = "-")]
        Minus,
    }

    /// Representation of an exon in the JSON.
    #[derive(Deserialize, Debug)]
    struct ExonHelper(i32, i32, i32, i32, i32, Option<String>);

    /// Representation of an exon after loading.
    #[derive(Debug)]
    pub struct Exon {
        /// Start position on reference.
        pub alt_start_i: i32,
        /// End position on reference.
        pub alt_end_i: i32,
        /// CDS end coordinates.
        pub alt_exon_id: i32,
        /// CDS start coordinate.
        pub alt_cds_start_i: i32,
        /// CDS end coordinate.
        pub alt_cds_end_i: i32,
        /// Alignment of an exon in CIGAR format.
        pub cigar: String,
    }

    /// Representation of `transcripts.*.genome_builds` value.
    #[derive(Deserialize, Debug)]
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
    #[derive(Deserialize, Debug)]
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

    #[derive(Deserialize, Debug)]
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

        for gap_op in gap.split(" ") {
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
                |ExonHelper(
                    alt_start_i,
                    alt_end_i,
                    alt_exon_id,
                    alt_cds_start_i,
                    alt_cds_end_i,
                    gap,
                )| Exon {
                    alt_start_i,
                    alt_end_i,
                    alt_exon_id,
                    alt_cds_start_i,
                    alt_cds_end_i,
                    cigar: if let Some(gap) = gap {
                        gap_to_cigar(&gap)
                    } else {
                        format!("{}M", alt_end_i - alt_start_i + 1)
                    },
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
            _ => panic!("Invalid transcript tag {:?}", s),
        }
    }

    #[derive(Debug, Deserialize)]
    struct WrappedString(String);

    fn deserialize_tag<'de, D>(deserializer: D) -> Result<Option<Vec<Tag>>, D::Error>
    where
        D: Deserializer<'de>,
    {
        Option::<WrappedString>::deserialize(deserializer).map(|opt_wrapped| {
            opt_wrapped.map(|wrapped| {
                wrapped
                    .0
                    .split(",")
                    .into_iter()
                    .map(|s| str_to_tag(s))
                    .collect()
            })
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
            _ => panic!("Unknown biotype {:?}", s),
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
                buf.split(",")
                    .into_iter()
                    .map(|s| str_to_biotype(s))
                    .collect()
            }
        }))
    }
}

/// Internal implementation of the transcript provider.
struct TxProvider {}

#[cfg(test)]
pub mod tests {
    use pretty_assertions::assert_eq;

    use super::models::{gap_to_cigar, Container};

    #[test]
    fn test_deserialize_brca1() -> Result<(), anyhow::Error> {
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
    fn test_gap_to_cigar() {
        assert_eq!(gap_to_cigar("M196 I1 M61 I1 M181"), "196=1D61=1D181=");
    }

    #[cfg(skip_tests)]
    #[test]
    fn test_deserialize_big_files() -> Result<(), anyhow::Error> {
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
}
