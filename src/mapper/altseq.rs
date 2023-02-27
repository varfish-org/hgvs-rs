//! Code for building alternative sequence and convertion to HGVS.p.

use std::rc::Rc;

use crate::{
    data::interface::Provider,
    parser::HgvsVariant,
    sequences::{translate_cds, TranslationTable},
};

pub struct RefTranscriptData {
    /// Transcript nucleotide sequence.
    pub transcript_sequence: String,
    /// Translated amino acid sequence.
    pub aa_sequence: String,
    /// 1-based CDS start position on transcript.
    pub cds_start: i32,
    /// 1-based CDS end position on transcript.
    pub cds_stop: i32,
    /// Accession of the protein or `MD5_${md5sum}`.
    pub protein_accession: String,
}

impl RefTranscriptData {
    /// Construct new instance fetching data from the provider..
    ///
    /// # Args
    ///
    /// * `provider` -- Data provider to query.
    /// * `tx_ac` -- Transcript accession.
    /// * `pro_ac` -- Protein accession.
    pub fn new(
        provider: Rc<dyn Provider>,
        tx_ac: &str,
        pro_ac: Option<&str>,
    ) -> Result<Self, anyhow::Error> {
        let tx_info = provider.as_ref().get_tx_identity_info(tx_ac)?;
        let transcript_sequence = provider.as_ref().get_seq(tx_ac)?;

        // Use 1-based HGVS coordinates.
        let cds_start = tx_info.cds_start_i + 1;
        let cds_stop = tx_info.cds_end_i;

        // Coding sequences taht are not divisable by 3 are not yet supported.
        let tx_seq_to_translate =
            &transcript_sequence[((cds_start - 1) as usize)..(cds_stop as usize)];
        if tx_seq_to_translate.len() % 3 != 0 {
            return Err(anyhow::anyhow!(
                "Transcript {} is not supported because its sequence length of {} is not divible by 3.",
                tx_ac,
                tx_seq_to_translate.len()));
        }

        let aa_sequence =
            translate_cds(tx_seq_to_translate, true, "*", TranslationTable::Standard)?;
        let protein_accession = if let Some(pro_ac) = pro_ac {
            pro_ac.to_owned()
        } else if let Some(pro_ac) = provider.as_ref().get_pro_ac_for_tx_ac(tx_ac)? {
            pro_ac
        } else {
            // get_acs_for_protein_seq() will always return at least the MD5_ accession.
            //
            // NB: the following comment is from the original Python code.
            //
            // TODO: drop get_acs_for_protein_seq; use known mapping or digest (wo/pro ac inference)
            provider
                .as_ref()
                .get_acs_for_protein_seq(&aa_sequence)?
                .into_iter()
                .next()
                .unwrap()
        };

        Ok(Self {
            transcript_sequence,
            aa_sequence,
            cds_start,
            cds_stop,
            protein_accession,
        })
    }
}

/// Utility to insert an hgvs variant into a transcript sequence.
///
/// Generates a record corresponding to the modified transcript sequence, along with annotations
/// for use in conversion to an hgvsp tag.  Used in hgvsc to hgvsp conversion.
pub struct AltSeqBuilder<'a> {
    pub reference_data: &'a RefTranscriptData,
}

pub struct AltData {}

impl<'a> AltSeqBuilder<'a> {
    pub fn new(reference_data: &'a RefTranscriptData) -> Self {
        Self { reference_data }
    }

    pub fn build_altseq(&self) -> Result<Vec<AltData>, anyhow::Error> {
        todo!()
    }
}

pub struct AltSeqToHgvsp<'a> {
    pub var_c: &'a HgvsVariant,
    pub reference_data: &'a RefTranscriptData,
}

impl<'a> AltSeqToHgvsp<'a> {
    pub fn new(var_c: &'a HgvsVariant, reference_data: &'a RefTranscriptData) -> Self {
        Self {
            var_c,
            reference_data,
        }
    }

    pub fn build_hgvsp(&self) -> Result<HgvsVariant, anyhow::Error> {
        todo!()
    }
}
