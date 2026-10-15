//! Code for building alternative sequence and convertion to HGVS.p.

use std::{cmp::Ordering, sync::Arc};

use cached::proc_macro::cached;
use cached::LruCache;

use crate::{
    data::interface::{Provider, TranslationException},
    mapper::error::Error,
    parser::{
        Accession, CdsFrom, CdsPos, HgvsVariant, Mu, NaEdit, ProtInterval, ProtLocEdit, ProtPos,
        ProteinEdit, UncertainLengthChange,
    },
    sequences::{revcomp, translate_cds_with_exceptions, TranslationTable},
};

/// How the CDS of a transcript ends.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CdsEnd {
    /// Ends with a stop codon.
    StopCodon,
    /// Complete codons, no stop codon.
    NoStopCodon,
    /// Trailing bases of an incomplete codon.
    IncompleteCodon(String),
}

#[derive(Debug, Clone)]
pub struct RefTranscriptData {
    /// Transcript nucleotide sequence.
    pub transcript_sequence: Arc<str>,
    /// Translated amino acid sequence.
    pub aa_sequence: Arc<str>,
    /// 1-based CDS start position on transcript.
    pub cds_start: i32,
    /// 1-based CDS end position on transcript.
    pub cds_stop: i32,
    /// Accession of the protein or `MD5_${md5sum}`.
    pub protein_accession: Arc<str>,
    /// The translation table to use.
    pub translation_table: TranslationTable,
    /// Amino acids that the annotation assigns to single codons, whatever `translation_table`
    /// says.
    pub translation_exceptions: Vec<TranslationException>,
    /// c. base of the first base of each exon.  As in biocommons/hgvs, this is exact within
    /// the CDS only: it skips c.0 for the first exon only and ignores the `*` of the 3' UTR.
    pub exon_start_positions: Vec<i32>,
    /// c. base of the last base of each exon, see `exon_start_positions`.
    pub exon_end_positions: Vec<i32>,
    /// How the CDS ends.  `aa_sequence` holds its complete codons only.
    pub cds_end: CdsEnd,
}

/// The 0-based codon indices and the amino acids of the translation `exceptions`.
fn codon_exceptions(exceptions: &[TranslationException]) -> Vec<(usize, char)> {
    exceptions
        .iter()
        .filter_map(|exception| {
            let codon = usize::try_from(exception.position).ok()?.checked_sub(1)?;
            Some((codon, exception.amino_acid))
        })
        .collect()
}

/// The 0-based index in `ref_seq` after the last base that differs from `alt_seq`, or `None`
/// if both are equal.  The differing bases are placed as far 3' as possible: the common prefix
/// is removed first, then the common suffix of the rest.
fn changed_end(ref_seq: &str, alt_seq: &str) -> Option<usize> {
    if ref_seq == alt_seq {
        return None;
    }
    let (ref_seq, alt_seq) = (ref_seq.as_bytes(), alt_seq.as_bytes());
    let prefix = ref_seq
        .iter()
        .zip(alt_seq)
        .take_while(|(r, a)| r == a)
        .count();
    let suffix = ref_seq[prefix..]
        .iter()
        .rev()
        .zip(alt_seq[prefix..].iter().rev())
        .take_while(|(r, a)| r == a)
        .count();
    Some(ref_seq.len() - suffix)
}

/// The amino acids of `aa_seq` up to and including its first stop codon, or all of them if it
/// has none.
fn up_to_first_stop(aa_seq: &str) -> &str {
    aa_seq.find('*').map_or(aa_seq, |pos| &aa_seq[..=pos])
}

#[cached(
    ty = "LruCache<(String, String, String, Option<String>), Result<RefTranscriptData, Error>>",
    create = "{ LruCache::builder().max_size(1000).build().expect(\"valid cache size\") }",
    cache_err = true,
    convert = r#"{ (
        provider.data_version().to_string(),
        provider.schema_version().to_string(),
        tx_ac.to_string(),
        pro_ac.map(|s| s.to_string())
    ) }"#
)]
pub fn ref_transcript_data_cached(
    provider: Arc<dyn Provider + Send + Sync>,
    tx_ac: &str,
    pro_ac: Option<&str>,
) -> Result<RefTranscriptData, Error> {
    RefTranscriptData::new(provider, tx_ac, pro_ac)
}

impl RefTranscriptData {
    /// Construct new instance fetching data from the provider..
    ///
    /// If the length of the CDS is not a multiple of 3, the protein ends with the last
    /// complete codon, and `cds_end` holds the bases after it.  Unlike biocommons/hgvs, which
    /// rejects such a CDS, this returns an error only if a stop codon is in frame.
    ///
    /// # Args
    ///
    /// * `provider` -- Data provider to query.
    /// * `tx_ac` -- Transcript accession.
    /// * `pro_ac` -- Protein accession.
    pub fn new(
        provider: Arc<dyn Provider + Send + Sync>,
        tx_ac: &str,
        pro_ac: Option<&str>,
    ) -> Result<Self, Error> {
        let tx_info = provider.as_ref().get_tx_identity_info(tx_ac)?;
        let transcript_sequence: Arc<str> = provider.as_ref().get_seq(tx_ac)?.into();

        // Non-coding transcripts (e.g., NR_*) don't have a CDS defined.
        let (cds_start_i, cds_end_i) = match (tx_info.cds_start_i, tx_info.cds_end_i) {
            (Some(start), Some(end)) => (start, end),
            _ => return Err(Error::CdsUndefined(tx_ac.to_string())),
        };

        // Use 1-based HGVS coordinates.
        let cds_start = cds_start_i + 1;
        let cds_stop = cds_end_i;

        // Translate up to the last complete codon.  A CDS with an incomplete 3' end, e.g.,
        // GENCODE `cds_end_NF`, can end with an incomplete codon.
        let tx_seq_to_translate =
            &transcript_sequence[((cds_start - 1) as usize)..(cds_stop as usize)];
        let full_codons_len = tx_seq_to_translate.len() - tx_seq_to_translate.len() % 3;

        let translation_exceptions = provider.as_ref().get_tx_translation_exceptions(tx_ac)?;
        let aa_sequence: Arc<str> = translate_cds_with_exceptions(
            &tx_seq_to_translate[..full_codons_len],
            true,
            "*",
            tx_info.translation_table,
            &codon_exceptions(&translation_exceptions),
        )?
        .into();
        let cds_end = match &tx_seq_to_translate[full_codons_len..] {
            "" if aa_sequence.ends_with('*') => CdsEnd::StopCodon,
            "" => CdsEnd::NoStopCodon,
            // A stop codon in frame means that the CDS changes the frame instead, e.g., by
            // ribosomal slippage as in OAZ1.  A translation in one frame would be wrong.
            _ if aa_sequence.contains('*') => {
                return Err(Error::TranscriptLengthInvalid(
                    tx_ac.to_string(),
                    tx_seq_to_translate.len(),
                ))
            }
            bases => CdsEnd::IncompleteCodon(bases.to_string()),
        };
        let protein_accession: Arc<str> = if let Some(pro_ac) = pro_ac {
            pro_ac.into()
        } else if let Some(pro_ac) = provider.as_ref().get_pro_ac_for_tx_ac(tx_ac)? {
            pro_ac.into()
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
                .expect(
                    "get_acs_for_protein_seq() should always return at least the MD5_ accession",
                )
                .into()
        };

        let mut exon_start_positions = Vec::with_capacity(tx_info.lengths.len());
        let mut exon_end_positions: Vec<i32> = Vec::with_capacity(tx_info.lengths.len());
        for length in &tx_info.lengths {
            let (start, end) = match exon_end_positions.last() {
                Some(last_end) => (last_end + 1, last_end + length),
                None => (-cds_start_i, -cds_start_i + length),
            };
            exon_start_positions.push(start);
            exon_end_positions.push(end);
        }

        Ok(Self {
            transcript_sequence,
            aa_sequence,
            cds_start,
            cds_stop,
            protein_accession,
            translation_table: tx_info.translation_table,
            translation_exceptions,
            exon_start_positions,
            exon_end_positions,
            cds_end,
        })
    }
}

#[derive(Debug, Clone)]
pub struct AltTranscriptData {
    /// Transcript nucleotide sequence.
    #[allow(dead_code)]
    pub transcript_sequence: String,
    /// 1-letter amino acid sequence.
    pub aa_sequence: String,
    /// 1-based CDS start position.
    #[allow(dead_code)]
    pub cds_start: i32,
    /// 1-based CDS stop position.
    #[allow(dead_code)]
    pub cds_stop: i32,
    /// Protein accession number, e.g., `"NP_999999.2"`.
    #[allow(dead_code)]
    pub protein_accession: String,
    /// Whether this is a frameshift variant.
    pub is_frameshift: bool,
    /// 1-based AA start index for this variant.
    pub variant_start_aa: Option<i32>,
    /// Starting position (AA ref index) of the last framewshift which affects the rest of the
    /// sequence, ie.e., not offset by subsequent frameshifts.
    pub frameshift_start: Option<i32>,
    /// Whether this is a substitution AA variant.
    pub is_substitution: bool,
    /// Whether variant is "?".
    pub is_ambiguous: bool,
}

impl AltTranscriptData {
    /// Create a variant sequence using inputs from `VariantInserter`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        seq: &str,
        cds_start: i32,
        cds_stop: i32,
        is_frameshift: bool,
        variant_start_aa: Option<i32>,
        protein_accession: &str,
        ref_aa_sequence: &str,
        is_substitution: bool,
        is_ambiguous: bool,
        translation_table: TranslationTable,
    ) -> Result<Self, Error> {
        Self::new_owned(
            seq.to_owned(),
            cds_start,
            cds_stop,
            is_frameshift,
            variant_start_aa,
            protein_accession,
            ref_aa_sequence,
            is_substitution,
            is_ambiguous,
            translation_table,
            &[],
        )
    }

    /// Like `new`, but takes ownership of `seq`, and sets the amino acids of
    /// `translation_exceptions` at their 0-based codons of the CDS.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new_owned(
        seq: String,
        cds_start: i32,
        cds_stop: i32,
        is_frameshift: bool,
        variant_start_aa: Option<i32>,
        protein_accession: &str,
        ref_aa_sequence: &str,
        is_substitution: bool,
        is_ambiguous: bool,
        translation_table: TranslationTable,
        translation_exceptions: &[(usize, char)],
    ) -> Result<Self, Error> {
        let transcript_sequence = seq;
        let aa_sequence = if !transcript_sequence.is_empty() {
            // In case of SEPHS2 / HGNC:19686, the last amino acid is both a selenocysteine
            // and a stop codon.
            // We handle this by explicitly truncating the sequence at the Sec + stop codon.
            // This heuristic may not always be correct;
            // alternatively/additionally, we could check `protein_accession` for known cases.
            let seq_cds = if translation_table == TranslationTable::Selenocysteine
                && ref_aa_sequence.ends_with('U')
            {
                &transcript_sequence[((cds_start - 1) as usize)..cds_stop as usize]
            } else {
                &transcript_sequence[((cds_start - 1) as usize)..]
            };

            let seq_aa = if variant_start_aa.is_some() {
                translate_cds_with_exceptions(
                    seq_cds,
                    false,
                    "X",
                    translation_table,
                    translation_exceptions,
                )?
            } else {
                ref_aa_sequence.to_owned()
            };
            // The alt protein ends at `cds_stop`, which the caller moves by the length change
            // of the variant.  So after an in-frame change, the alt protein lines up with
            // `ref_aa_sequence`, and `AltSeqToHgvsp` can find the changed amino acids.  With a
            // stop codon, the cut is after the last stop codon up to `cds_stop`.  The alt
            // protein thus keeps a new stop codon in front of the old one.  With more than one
            // stop codon in the reference, it keeps the ones after the first (#96).  If no stop
            // codon is up to `cds_stop`, the alt protein goes on to the next stop codon or to
            // the end of the transcript.  Without a stop codon, e.g., in ENST00000420031.2, the
            // amino acids after `cds_stop` are unknown, so the alt protein ends at `cds_stop`.
            let orig_aa_len = (cds_stop - cds_start + 1) as usize / 3;
            let orig_aa_len = std::cmp::min(orig_aa_len, seq_aa.len());
            if !ref_aa_sequence.contains('*') {
                seq_aa[..orig_aa_len].to_owned()
            } else if let Some(stop_pos) = seq_aa[..orig_aa_len]
                .rfind('*')
                .or_else(|| seq_aa.find('*'))
            {
                seq_aa[..(stop_pos + 1)].to_owned()
            } else {
                seq_aa
            }
        } else {
            "".to_owned()
        };

        Ok(Self {
            transcript_sequence,
            aa_sequence,
            cds_start,
            cds_stop,
            protein_accession: protein_accession.to_owned(),
            is_frameshift,
            variant_start_aa,
            frameshift_start: None,
            is_substitution,
            is_ambiguous,
        })
    }
}

/// Utility enum for locating variant in a transcript.
#[derive(Debug, PartialEq)]
enum VariantLocation {
    Exon,
    Intron,
    FivePrimeUtr,
    ThreePrimeUtr,
    WholeGene,
}

/// Utility enum for edit type.
enum EditType {
    NaRefAlt,
    Dup,
    Inv,
    NotCds,
    WholeGeneDeleted,
}

/// Utility to insert an hgvs variant into a transcript sequence.
///
/// Generates a record corresponding to the modified transcript sequence, along with annotations
/// for use in conversion to an hgvsp tag.  Used in hgvsc to hgvsp conversion.
pub struct AltSeqBuilder<'a> {
    /// HgvsVariant::CdsVariant to build the alternative sequence for.
    pub var_c: HgvsVariant,
    /// Information about the transcript reference.
    pub reference_data: &'a RefTranscriptData,
    /// Whether the transcript reference amino acid sequence has multiple stop codons.
    pub ref_has_multiple_stops: bool,
    /// 0-based index of the first stop codon in the reference amino acid sequence.  It is
    /// difficult to handle transcripts with multiple stop codons, but we can make a reasonable
    /// prediction up to the first stop codon.  This rescues changes to the `p.Met1`.
    pub first_stop_pos: Option<usize>,
    /// See `variant::Config::ins_at_boundary_is_intronic`; `true` by default.
    pub ins_at_boundary_is_intronic: bool,
}

impl<'a> AltSeqBuilder<'a> {
    pub fn new(var_c: HgvsVariant, reference_data: &'a RefTranscriptData) -> Self {
        let ref_has_multiple_stops = reference_data.aa_sequence.matches('*').count() > 1;
        let first_stop_pos = reference_data.aa_sequence.find('*');

        Self {
            var_c,
            reference_data,
            ref_has_multiple_stops,
            first_stop_pos,
            ins_at_boundary_is_intronic: true,
        }
    }

    /// Given a variant and a sequence, incorporate the variant and return the new sequence.
    ///
    /// Data structure returned is analogous to the data structure used to return the variant
    /// sequence, but with an additional parameter denoting the start of a frameshift that should
    /// affect all bases downstream.
    ///
    /// # Returns
    ///
    /// Variant sequence data, or `Error::ExpectedCdsVariant` if `var_c` is not a
    /// `HgvsVariant::CdsVariant`.
    pub fn build_altseq(&self) -> Result<Vec<AltTranscriptData>, Error> {
        // NB: the following comment is from the original Python code.
        // Should loop over each allele rather than assume only 1 variant; return a list for now.

        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return Err(Error::ExpectedCdsVariant(self.var_c.to_string()));
        };
        let na_edit = loc_edit.edit.inner();
        let edit_type = match self.get_variant_region() {
            // Without a stop codon, the amino acids after the CDS are unknown.
            VariantLocation::Exon if self.after_cds_without_stop() => EditType::NotCds,
            VariantLocation::Exon => match na_edit {
                NaEdit::RefAlt { .. }
                | NaEdit::NumAlt { .. }
                | NaEdit::DelRef { .. }
                | NaEdit::DelNum { .. }
                | NaEdit::Ins { .. } => EditType::NaRefAlt,
                NaEdit::Dup { .. } => EditType::Dup,
                NaEdit::InvRef { .. } | NaEdit::InvNum { .. } => EditType::Inv,
            },
            VariantLocation::Intron
            // NB: the following comment is from the original Python code
            // TODO: handle case where variatn introduces a `Met` (new start)
            | VariantLocation::FivePrimeUtr
            | VariantLocation::ThreePrimeUtr => EditType::NotCds,
            VariantLocation::WholeGene => match na_edit {
                NaEdit::DelRef { .. } |
                NaEdit::DelNum { .. } => EditType::WholeGeneDeleted,
                NaEdit::Dup { .. } => {
                    tracing::warn!("Whole-gene duplication; consequence assumed to not affect protein product");
                    EditType::NotCds
                }
                NaEdit::InvRef { .. } |
                NaEdit::InvNum { .. } => {
                    tracing::warn!("Whole-gene inversion; consequence assumed to not affected protein product");
                    EditType::NotCds
                }
                NaEdit::RefAlt { .. } => {
                    tracing::warn!("The whole-gene variant {} is not a clean deletion. Assuming whole gene deletion.", self.var_c);
                    EditType::WholeGeneDeleted
                }
                _ => panic!("Invalid combination of whole gene variant location and NaEdit {na_edit:?}"),
            },
        };

        // Get the start of the "terminal" frameshift (i.e. one never "cancelled out").
        let alt_data = match edit_type {
            EditType::NaRefAlt => self.incorporate_delins(),
            EditType::Dup => self.incorporate_dup(),
            EditType::Inv => self.incorporate_inv(),
            EditType::NotCds => self.create_alt_equals_ref_noncds(),
            EditType::WholeGeneDeleted => self.create_no_protein(),
        }?;

        let alt_data = self.get_frameshift_start(alt_data);

        Ok(vec![alt_data])
    }

    /// Categorize variant by location in transcript.
    fn get_variant_region(&self) -> VariantLocation {
        match &self.var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => {
                let loc = loc_edit.loc.inner();
                let edit = loc_edit.edit.inner();
                let start_offset = loc.start.offset.unwrap_or(0);
                let end_offset = loc.end.offset.unwrap_or(0);
                let ref_data = self.reference_data;
                let at_exon_edge = ref_data.exon_start_positions.contains(&loc.start.base)
                    || ref_data.exon_end_positions.contains(&loc.end.base);
                // Without a stop codon, the CDS goes on after its last base, so there is no 3'
                // UTR, see `outside_cds()`.
                let has_stop_codon = ref_data.cds_end == CdsEnd::StopCodon;
                #[allow(clippy::nonminimal_bool)]
                if has_stop_codon
                    && (loc.start.cds_from == CdsFrom::End && loc.end.cds_from == CdsFrom::End
                        || ((edit.is_ins() || edit.is_dup()) && loc.end.cds_from == CdsFrom::End))
                {
                    // 3' UTR if both ends are 3' UTR or if the variant is an insertion/duplication
                    // and the end is 3' UTR.
                    VariantLocation::ThreePrimeUtr
                } else if loc.start.base < 0 && loc.end.base < 0 {
                    VariantLocation::FivePrimeUtr
                } else if loc.start.base < 0 && loc.end.cds_from == CdsFrom::End && !edit.is_dup() {
                    // Without a stop codon, a dup gets here, too.  Its copy can go before c.1.
                    VariantLocation::WholeGene
                } else if let Some(utr) = self.outside_cds() {
                    utr
                } else if self.ins_at_boundary_is_intronic && edit.is_dup() && at_exon_edge {
                    // The copy may sit at the exon boundary, like an insertion there.
                    VariantLocation::Intron
                } else if !self.ins_at_boundary_is_intronic
                    && ((edit.is_ins() && matches!((start_offset, end_offset), (-1, 0) | (0, 1)))
                        || (edit.is_dup() && (end_offset == -1 || start_offset == 1)))
                {
                    // The splice site stays intact, so the inserted bases extend the exon.
                    VariantLocation::Exon
                } else if loc.start.offset.is_some() || loc.end.offset.is_some() {
                    // Leave out anything intronic for now.
                    VariantLocation::Intron
                } else {
                    // Anything else that contains an exon.
                    VariantLocation::Exon
                }
            }
            _ => panic!("Must be CDS variant"),
        }
    }

    /// The UTR of a variant that leaves the CDS intact, or `None`.
    ///
    /// An ins, dup, del, delins, or inv leaves the CDS intact if one of its equivalent
    /// placements lies before c.1 or after the last base of the CDS.  This depends on the
    /// sequence, not on the description.  For example, `c.1dup` is `c.-1_1insA` if c.-1 is not
    /// an A, `c.1del` is `c.-1del` if c.-1 is an A, and `c.(N-2)_Ndup` puts its copy after c.N,
    /// the last base of the CDS.  The inserted bases can still form an upstream ATG, so the
    /// effect on the protein is unknown, not none.
    ///
    /// Without a stop codon, the CDS goes on after its last base, so there is no 3' UTR, see
    /// `after_cds_without_stop()`.
    fn outside_cds(&self) -> Option<VariantLocation> {
        let (before_cds, after_cds) = self.placements_outside_cds()?;
        if before_cds {
            Some(VariantLocation::FivePrimeUtr)
        } else if after_cds && self.reference_data.cds_end == CdsEnd::StopCodon {
            Some(VariantLocation::ThreePrimeUtr)
        } else {
            None
        }
    }

    /// Whether the CDS has no stop codon and the variant changes only bases after its last base
    /// c.N, e.g., `c.*1G>A`, or `c.Ndup` as its copy can go after c.N.  The CDS goes on after
    /// c.N, but its amino acids there are unknown.
    fn after_cds_without_stop(&self) -> bool {
        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return false;
        };
        self.reference_data.cds_end != CdsEnd::StopCodon
            && (loc_edit.loc.inner().start.cds_from == CdsFrom::End
                || self
                    .placements_outside_cds()
                    .is_some_and(|(_, after_cds)| after_cds))
    }

    /// Whether an equivalent placement of the variant lies before c.1, and whether one lies
    /// after the last base of the CDS, see `outside_cds()`.  `None` for the edits and positions
    /// that it does not check, e.g., a delins with an intron position.
    ///
    /// With `ins_at_boundary_is_intronic`, an ins or dup with an intron position stays
    /// intronic.  The other variants count only with exon positions.
    fn placements_outside_cds(&self) -> Option<(bool, bool)> {
        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return None;
        };
        let loc = loc_edit.loc.inner();
        let ref_data = self.reference_data;
        let seq = ref_data.transcript_sequence.as_bytes();
        let cds_start = usize::try_from(ref_data.cds_start - 1).ok()?;
        let cds_end = usize::try_from(ref_data.cds_stop).ok()?;
        // The 0-based index in `seq` of the base c.`pos.base`, ignoring the intron offset.
        let index = |pos: &CdsPos| {
            let base = i64::from(pos.base);
            let one_based = match pos.cds_from {
                // There is no c.0, c.-1 is the base before c.1.
                CdsFrom::Start if base < 0 => i64::from(ref_data.cds_start) + base,
                CdsFrom::Start => i64::from(ref_data.cds_start) + base - 1,
                CdsFrom::End => i64::from(ref_data.cds_stop) + base,
            };
            usize::try_from(one_based - 1)
                .ok()
                .filter(|&index| index < seq.len())
        };
        let offsets = (loc.start.offset.unwrap_or(0), loc.end.offset.unwrap_or(0));
        let at_boundary_is_exonic = !self.ins_at_boundary_is_intronic;
        // The variant replaces the bases `start..end` of `seq` with `alt_bases`.
        let inverted: String;
        let (start, end, alt_bases) = match (loc_edit.edit.inner(), offsets) {
            (NaEdit::Ins { alternative }, (0, 0)) => {
                let gap = index(&loc.start)? + 1;
                (gap, gap, alternative.as_str())
            }
            // At an exon boundary, the inserted bases go into the exon, as in
            // `get_variant_region()`.
            (NaEdit::Ins { alternative }, (-1, 0)) if at_boundary_is_exonic => {
                let gap = index(&loc.end)?;
                (gap, gap, alternative.as_str())
            }
            (NaEdit::Ins { alternative }, (0, 1)) if at_boundary_is_exonic => {
                let gap = index(&loc.start)? + 1;
                (gap, gap, alternative.as_str())
            }
            // The copy can go before the duplicated bases.  If these start before c.1, the copy
            // leaves the CDS intact, even if they start before the transcript.
            (NaEdit::Dup { .. }, (0, 0)) if loc.start.base < 0 => return Some((true, false)),
            (NaEdit::Dup { .. }, (0, 0)) => {
                let start = index(&loc.start)?;
                let end = index(&loc.end)?;
                (start, start, ref_data.transcript_sequence.get(start..=end)?)
            }
            // The copy of intron bases can go to the exon boundary.  Only the variant can tell
            // these bases.
            (NaEdit::Dup { reference }, (_, -1)) if at_boundary_is_exonic => {
                let gap = index(&loc.end)?;
                (gap, gap, reference.as_str())
            }
            (NaEdit::Dup { reference }, (1, _)) if at_boundary_is_exonic => {
                let gap = index(&loc.start)? + 1;
                (gap, gap, reference.as_str())
            }
            (NaEdit::DelRef { .. } | NaEdit::DelNum { .. }, (0, 0)) => {
                (index(&loc.start)?, index(&loc.end)? + 1, "")
            }
            // An inversion is a delins of the reverse complement.
            (NaEdit::InvRef { .. } | NaEdit::InvNum { .. }, (0, 0)) => {
                let start = index(&loc.start)?;
                let end = index(&loc.end)? + 1;
                inverted = revcomp(ref_data.transcript_sequence.get(start..end)?);
                (start, end, inverted.as_str())
            }
            // An identity, e.g., `c.1A=`, changes no base.
            (
                NaEdit::RefAlt {
                    reference,
                    alternative,
                },
                _,
            ) if reference == alternative => return None,
            (NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. }, (0, 0)) => (
                index(&loc.start)?,
                index(&loc.end)? + 1,
                alternative.as_str(),
            ),
            _ => return None,
        };
        if start == end && alt_bases.is_empty() {
            // Without the inserted bases, only this placement is known.
            Some((start <= cds_start, start >= cds_end))
        } else {
            // `alt` ends with the reference from c.1 on if and only if a placement in front of
            // c.1 gives `alt`.  Likewise, `alt` starts with the reference up to the last base
            // of the CDS if and only if a placement after that base gives `alt`.
            let alt = [seq.get(..start)?, alt_bases.as_bytes(), seq.get(end..)?].concat();
            if alt == seq {
                // Changes no base, e.g., the inversion of `GC`.
                return None;
            }
            Some((
                alt.ends_with(seq.get(cds_start..)?),
                alt.starts_with(seq.get(..cds_end)?),
            ))
        }
    }

    /// Helper to setup incorporate function.
    ///
    /// # Returns
    ///
    /// `(transcript sequence, cds start [1-based], cds stop [1-based], cds start index in
    /// seq [inc, 0-based], cds end index in seq [excl, 0-based])`
    fn setup_incorporate(&self) -> (String, i32, i32, usize, usize) {
        let mut start_end = Vec::new();

        match &self.var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => {
                for pos in &[&loc_edit.loc.inner().start, &loc_edit.loc.inner().end] {
                    match pos.cds_from {
                        CdsFrom::Start => {
                            if pos.base < 0 {
                                // 5' UTR
                                start_end.push(self.reference_data.cds_start as usize - 1);
                            } else if pos.offset.unwrap_or(0) < 0 {
                                // Intron before `pos.base`: the exon base before it.
                                start_end.push(
                                    self.reference_data.cds_start as usize - 1 + pos.base as usize
                                        - 2,
                                );
                            } else {
                                start_end.push(
                                    self.reference_data.cds_start as usize - 1 + pos.base as usize
                                        - 1,
                                );
                            }
                        }
                        CdsFrom::End => {
                            // 3' UTR
                            start_end.push(
                                self.reference_data.cds_stop as usize + pos.base as usize - 1,
                            );
                        }
                    }
                }
            }
            _ => panic!("invalid variant"),
        }

        let seq = self.reference_data.transcript_sequence.to_string();
        let start = std::cmp::min(start_end[0], seq.len());
        let end = std::cmp::min(start_end[1] + 1, seq.len());

        (
            seq,
            self.reference_data.cds_start,
            self.reference_data.cds_stop,
            start,
            end,
        )
    }

    /// The translation exceptions of the reference whose codons `alt_seq` keeps, as pairs of a
    /// 0-based codon of the CDS of `alt_seq` and the amino acid.
    ///
    /// The bases before the first and after the last difference to the reference are
    /// unchanged.  A codon of a translation exception in them keeps its amino acid if it is
    /// still in frame.
    fn alt_translation_exceptions(&self, alt_seq: &str) -> Vec<(usize, char)> {
        let exceptions = &self.reference_data.translation_exceptions;
        if exceptions.is_empty() {
            return Vec::new();
        }
        let Ok(cds_start) = usize::try_from(self.reference_data.cds_start - 1) else {
            return Vec::new();
        };
        let ref_seq = self.reference_data.transcript_sequence.as_bytes();
        let alt_seq = alt_seq.as_bytes();
        let prefix = ref_seq
            .iter()
            .zip(alt_seq)
            .take_while(|(r, a)| r == a)
            .count();
        let suffix = ref_seq
            .iter()
            .rev()
            .zip(alt_seq.iter().rev())
            .take(ref_seq.len().min(alt_seq.len()) - prefix)
            .take_while(|(r, a)| r == a)
            .count();

        codon_exceptions(exceptions)
            .into_iter()
            .filter_map(|(codon, amino_acid)| {
                let start = cds_start.checked_add(codon.checked_mul(3)?)?;
                let alt_start = if start.checked_add(3)? <= prefix {
                    start
                } else if start >= ref_seq.len() - suffix {
                    // `start + alt_seq.len() >= ref_seq.len()`, as `suffix <= alt_seq.len()`.
                    start + alt_seq.len() - ref_seq.len()
                } else {
                    return None;
                };
                let offset = alt_start.checked_sub(cds_start)?;
                (offset % 3 == 0).then_some((offset / 3, amino_acid))
            })
            .collect()
    }

    /// Whether the effect of a variant is unknown.  The effect is unknown if the reference has
    /// more than one stop codon and the variant changes the first one or a base after it.
    /// biocommons/hgvs gives `p.?` for every variant on such a transcript.
    ///
    /// The check takes the bases in which `alt_seq` differs from the reference transcript, at
    /// their 3'-most placement, so every description of a change gives the same result.  If no
    /// base differs, it takes `end`, the 0-based index after the last base of the variant.
    fn is_ambiguous(&self, alt_seq: &str, end: usize) -> bool {
        let Some(first_stop_pos) = self.first_stop_pos.filter(|_| self.ref_has_multiple_stops)
        else {
            return false;
        };
        let first_stop_start = (self.reference_data.cds_start - 1) as usize + 3 * first_stop_pos;
        let end = changed_end(&self.reference_data.transcript_sequence, alt_seq).unwrap_or(end);
        first_stop_start < end
    }

    /// Get starting position (AA ref index) of the last frameshift which affects the rest of
    /// the sequence, i.e. not offset by subsequent frameshifts.
    fn get_frameshift_start(&self, variant_data: AltTranscriptData) -> AltTranscriptData {
        AltTranscriptData {
            frameshift_start: if variant_data.is_frameshift {
                variant_data.variant_start_aa
            } else {
                variant_data.frameshift_start
            },
            ..variant_data
        }
    }

    fn incorporate_delins(&self) -> Result<AltTranscriptData, Error> {
        let (mut seq, cds_start, cds_stop, start, end) = self.setup_incorporate();
        let loc_range_start;

        let (reference, alternative) = match &self.var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => {
                loc_range_start = loc_edit.loc.inner().start.base;
                match loc_edit.edit.inner() {
                    NaEdit::RefAlt {
                        reference,
                        alternative,
                    } => (Some(reference.clone()), Some(alternative.clone())),
                    NaEdit::DelRef { reference } => (Some(reference.to_owned()), None),
                    NaEdit::Ins { alternative } => (None, Some(alternative.to_owned())),
                    _ => panic!("Can only work with concrete ref/alt"),
                }
            }
            _ => panic!("Can only work on CDS variants"),
        };
        let ref_len = reference.as_ref().map(|_| end - start).unwrap_or(0) as i32;
        let alt_len = alternative.as_ref().map(|alt| alt.len()).unwrap_or(0) as i32;
        let net_base_change = alt_len - ref_len;
        let cds_stop = cds_stop + net_base_change;

        // Incorporate the variant into the sequence (depending on the type).
        let mut is_substitution = false;
        let range = if end > seq.len() && reference.is_some() {
            tracing::warn!(
                    "Altered sequence range {:?} is incompatible with sequence length {:?}, clamping. Variant description is {}",
                    start..end,
                    seq.len(),
                    &self.var_c
                );
            start..seq.len()
        } else {
            start..end
        };
        match (reference, alternative) {
            (Some(reference), Some(alternative)) => {
                // delins or SNP
                seq.replace_range(range, &alternative);
                if reference.len() == 1 && alternative.len() == 1 {
                    is_substitution = true;
                }
            }
            (Some(_reference), None) => {
                // deletion
                seq.replace_range(range, "");
            }
            (None, Some(alternative)) => {
                // insertion
                seq.insert_str(start + 1, &alternative);
            }
            _ => panic!("This should not happen"),
        }

        let is_frameshift = net_base_change % 3 != 0;
        let is_ambiguous = self.is_ambiguous(&seq, end);

        // Use max. of mod 3 value and 1 (in the event that the indel starts in the 5' UTR range).
        let variant_start_aa = std::cmp::max((loc_range_start as f64 / 3.0).ceil() as i32, 1);
        let translation_exceptions = self.alt_translation_exceptions(&seq);

        AltTranscriptData::new_owned(
            seq,
            cds_start,
            cds_stop,
            is_frameshift,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            is_substitution,
            is_ambiguous,
            self.reference_data.translation_table,
            &translation_exceptions,
        )
    }

    fn incorporate_dup(&self) -> Result<AltTranscriptData, Error> {
        let (seq, cds_start, cds_stop, start, end) = self.setup_incorporate();

        // Intron bases are not in the transcript sequence, so take them from the variant.
        let dup_seq = if self.var_c.spans_intron() {
            match self.var_c.na_edit() {
                Some(NaEdit::Dup { reference }) if !reference.is_empty() => reference.as_str(),
                _ => return Err(Error::DuplicationReferenceMissing(self.var_c.to_string())),
            }
        } else {
            &seq[start..end]
        };
        let seq = format!("{}{}{}", &seq[..end], dup_seq, &seq[end..]);
        // `cds_stop` moves with the copy, as in `incorporate_delins()`, so that a dup gives the
        // same alt protein as the insertion of its copy.  The alt protein then keeps all of the
        // CDS: the stop codons after the copy (#96) and, without a stop codon, its last codons.
        // biocommons/hgvs keeps `cds_stop`, and biocommons/hgvs#896 proposes the same move.
        let cds_stop = cds_stop + dup_seq.len() as i32;

        let is_frameshift = dup_seq.len() % 3 != 0;

        let loc_end = match &self.var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => loc_edit.loc.inner().end.base,
            _ => panic!("can only work on CDS variants"),
        };
        let variant_start_aa = ((loc_end + 1) as f64 / 3.0).ceil() as i32;
        let translation_exceptions = self.alt_translation_exceptions(&seq);
        let is_ambiguous = self.is_ambiguous(&seq, end);

        AltTranscriptData::new_owned(
            seq,
            cds_start,
            cds_stop,
            is_frameshift,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            false,
            is_ambiguous,
            self.reference_data.translation_table,
            &translation_exceptions,
        )
    }

    /// Incorporate inversion into sequence.
    fn incorporate_inv(&self) -> Result<AltTranscriptData, Error> {
        let (seq, cds_start, cds_stop, start, end) = self.setup_incorporate();

        let seq = format!(
            "{}{}{}",
            &seq[..start],
            revcomp(&seq[start..end]),
            &seq[end..]
        );

        let loc_start = match &self.var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => loc_edit.loc.inner().start.base,
            _ => panic!("can only work on CDS variants"),
        };

        let variant_start_aa = std::cmp::max(((loc_start as f64) / 3.0).ceil() as i32, 1);
        let translation_exceptions = self.alt_translation_exceptions(&seq);
        let is_ambiguous = self.is_ambiguous(&seq, end);

        AltTranscriptData::new_owned(
            seq,
            cds_start,
            cds_stop,
            false,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            false,
            is_ambiguous,
            self.reference_data.translation_table,
            &translation_exceptions,
        )
    }

    /// Create an alt seq that matches the reference (for non-CDS variants).
    fn create_alt_equals_ref_noncds(&self) -> Result<AltTranscriptData, Error> {
        AltTranscriptData::new_owned(
            self.reference_data.transcript_sequence.to_string(),
            self.reference_data.cds_start,
            self.reference_data.cds_stop,
            false,
            None,
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            false,
            true,
            self.reference_data.translation_table,
            &[],
        )
    }

    /// Create a no-protein result.
    fn create_no_protein(&self) -> Result<AltTranscriptData, Error> {
        AltTranscriptData::new_owned(
            String::new(),
            -1,
            -1,
            false,
            None,
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            false,
            false,
            self.reference_data.translation_table,
            &[],
        )
    }
}

/// Build `HgvsVariant::ProtVariant` from information about change to transcript.
pub struct AltSeqToHgvsp<'a> {
    pub ref_data: &'a RefTranscriptData,
    pub alt_data: AltTranscriptData,
}

#[derive(Debug)]
struct AdHocRecord {
    start: i32,
    ins: String,
    del: String,
    is_frameshift: bool,
}

impl Default for AdHocRecord {
    fn default() -> Self {
        Self {
            start: -1,
            ins: "".to_owned(),
            del: "".to_owned(),
            is_frameshift: false,
        }
    }
}

impl<'a> AltSeqToHgvsp<'a> {
    pub fn new(ref_data: &'a RefTranscriptData, alt_data: AltTranscriptData) -> Self {
        Self { ref_data, alt_data }
    }

    /// Compare two amino acid sequences and generate a HGVS tag from the output.
    pub fn build_hgvsp(&self) -> Result<HgvsVariant, Error> {
        let mut records = Vec::new();

        if !self.alt_data.is_ambiguous && !self.alt_seq().is_empty() {
            let mut do_delins = true;

            let ref_seq = self.ref_seq();
            let alt_seq = self.alt_seq();
            let ref_len = ref_seq.len();
            let alt_len = alt_seq.len();

            // Helper to clamp indices safely and log a warning if out-of-bounds
            let clamp_index = |idx: usize, seq_len: usize, seq_name: &str| -> usize {
                let clamped = idx.min(seq_len);
                if clamped != idx {
                    log::warn!(
                        "Variant index {} exceeds {} sequence bounds (len {}). Clamping to {}.",
                        idx,
                        seq_name,
                        seq_len,
                        clamped
                    );
                }
                clamped
            };

            // Translation ends at the stop codon, so the protein is the same if `alt_seq` only
            // adds amino acids after the stop codon at the end of `ref_protein()`.  E.g., the
            // copy of `c.(N-4)_(N-2)dup` can form a new stop codon in front of the old one.
            let ref_protein = self.ref_protein();
            if ref_seq == alt_seq || ref_protein.ends_with('*') && alt_seq.starts_with(ref_protein)
            {
                // Silent p. variant.
                if let Some(start) = self.silent_start() {
                    let start_idx = start.checked_sub(1).and_then(|v| usize::try_from(v).ok());
                    let del = start_idx
                        .and_then(|i| ref_seq.as_bytes().get(i))
                        .map(|aa| (*aa as char).to_string())
                        .unwrap_or_default();
                    records.push(AdHocRecord {
                        start,
                        ins: del.clone(),
                        del,
                        is_frameshift: self.alt_data.is_frameshift,
                    });
                }
                do_delins = false;
            } else if self.is_substitution() && ref_len == alt_len {
                let mut diff_records = ref_seq
                    .as_bytes()
                    .iter()
                    .zip(alt_seq.as_bytes().iter())
                    .enumerate()
                    .filter(|(_i, (r, a))| r != a)
                    .map(|(i, (r, a))| AdHocRecord {
                        start: i as i32 + 1,
                        del: (*r as char).to_string(),
                        ins: (*a as char).to_string(),
                        is_frameshift: self.alt_data.is_frameshift,
                    })
                    .collect::<Vec<_>>();

                if diff_records.len() == 1 {
                    records.push(diff_records.drain(..).next().unwrap());
                    do_delins = false;
                }
            }

            if do_delins {
                let initial_start = self
                    .alt_data
                    .variant_start_aa
                    .expect("should not happen; must have start AA set")
                    .checked_sub(1)
                    .and_then(|v| usize::try_from(v).ok())
                    .unwrap_or(0);

                let safe_start = clamp_index(initial_start, ref_len.min(alt_len), "shared prefix");

                let matching_len = &ref_seq.as_bytes()[safe_start..]
                    .iter()
                    .zip(&alt_seq.as_bytes()[safe_start..])
                    .take_while(|(r, a)| r == a)
                    .count();

                let start = safe_start + matching_len;

                if self.alt_data.is_frameshift {
                    // Case: frameshifting delins or dup.
                    let ref_start = clamp_index(start, ref_len, "reference");
                    let alt_start = clamp_index(start, alt_len, "alt");

                    records.push(AdHocRecord {
                        start: start as i32 + 1,
                        ins: alt_seq[alt_start..].to_owned(),
                        del: ref_seq[ref_start..].to_owned(),
                        is_frameshift: self.alt_data.is_frameshift,
                    })
                } else {
                    // Case: non-frameshifting delins or dup.
                    let delta = alt_len as isize - ref_len as isize;
                    let offset = start + delta.unsigned_abs();

                    let ref_start = clamp_index(start, ref_len, "reference");
                    let ref_offset = clamp_index(offset, ref_len, "reference");
                    let alt_start = clamp_index(start, alt_len, "alt");
                    let alt_offset = clamp_index(offset, alt_len, "alt");

                    let (insertion, deletion, ref_sub, alt_sub) = match delta.cmp(&0) {
                        // if delta > 0 {
                        Ordering::Greater => (
                            // net insertion
                            alt_seq[alt_start..alt_offset].to_owned(),
                            "".to_string(),
                            ref_seq[ref_start..].to_owned(),
                            alt_seq[alt_offset..].to_owned(),
                        ),
                        Ordering::Less => (
                            // net deletion
                            "".to_string(),
                            ref_seq[ref_start..ref_offset].to_owned(),
                            ref_seq[ref_offset..].to_owned(),
                            alt_seq[alt_start..].to_owned(),
                        ),
                        Ordering::Equal => (
                            // size remains the same
                            "".to_string(),
                            "".to_string(),
                            ref_seq[ref_start..].to_owned(),
                            alt_seq[alt_start..].to_owned(),
                        ),
                    };

                    let diff_indices = ref_sub
                        .as_bytes()
                        .iter()
                        .zip(alt_sub.as_bytes())
                        .enumerate()
                        .filter(|(_i, (r, a))| r != a)
                        .map(|(i, _)| i)
                        .collect::<Vec<_>>();

                    let diff_indices = if diff_indices.is_empty()
                        && deletion.is_empty()
                        && insertion.starts_with('*')
                    {
                        vec![0]
                    } else {
                        diff_indices
                    };

                    let (deletion, insertion) = if !diff_indices.is_empty() {
                        let max_diff = diff_indices.last().unwrap() + 1;

                        let safe_ref_diff = max_diff.min(ref_sub.len());
                        let safe_alt_diff = max_diff.min(alt_sub.len());

                        (
                            format!("{}{}", deletion, &ref_sub[..safe_ref_diff]),
                            format!("{}{}", insertion, &alt_sub[..safe_alt_diff]),
                        )
                    } else {
                        (deletion, insertion)
                    };

                    records.push(AdHocRecord {
                        start: start as i32 + 1,
                        ins: insertion,
                        del: deletion,
                        is_frameshift: self.alt_data.is_frameshift,
                    });
                }
            }
        }

        if self.alt_data.is_ambiguous {
            Ok(self.create_variant(
                None,
                None,
                "",
                "",
                UncertainLengthChange::None,
                false,
                self.protein_accession(),
                self.alt_data.is_ambiguous,
                false,
                false,
                false,
                false,
                false,
            )?)
        } else if self.alt_seq().is_empty() {
            Ok(self.create_variant(
                None,
                None,
                "",
                "",
                UncertainLengthChange::None,
                false,
                self.protein_accession(),
                self.alt_data.is_ambiguous,
                false,
                false,
                true,
                false,
                false,
            )?)
        } else if let Some(var) = records.drain(..).next() {
            Ok(self.convert_to_hgvs_variant(var, self.protein_accession())?)
        } else {
            Err(Error::MultipleAAVariants)
        }
    }

    /// The 1-based position of the amino acid of a silent change.  This is the codon of the
    /// first base that differs between the reference and the alt transcript.  By the 3' rule,
    /// the variant changes this base, whatever its description.  E.g., if the CDS ends with
    /// `AAA TAG`, `c.25_26insAAT` and `c.27_29dup` both change the stop codon to `TAA`.  If the
    /// transcripts are the same in the CDS, e.g., for `c.30=`, this is `variant_start_aa`.
    fn silent_start(&self) -> Option<i32> {
        let ref_seq = self.ref_data.transcript_sequence.as_bytes();
        let alt_seq = self.alt_data.transcript_sequence.as_bytes();
        let cds_start = usize::try_from(self.ref_data.cds_start - 1).ok()?;
        let first_diff = ref_seq
            .iter()
            .zip(alt_seq)
            .take_while(|(r, a)| r == a)
            .count();
        if ref_seq != alt_seq && first_diff >= cds_start {
            i32::try_from((first_diff - cds_start) / 3 + 1).ok()
        } else {
            self.alt_data.variant_start_aa
        }
    }

    fn protein_accession(&self) -> &str {
        &self.ref_data.protein_accession
    }

    fn ref_seq(&self) -> &str {
        &self.ref_data.aa_sequence
    }

    /// The part of `ref_seq()` that the prediction covers.  With more than one stop codon, the
    /// prediction ends with the first one (#96), see `AltSeqBuilder::is_ambiguous()`.
    fn ref_protein(&self) -> &str {
        let ref_seq = self.ref_seq();
        if ref_seq.matches('*').count() > 1 {
            up_to_first_stop(ref_seq)
        } else {
            ref_seq
        }
    }

    fn alt_seq(&self) -> &str {
        &self.alt_data.aa_sequence
    }

    fn is_substitution(&self) -> bool {
        self.alt_data.is_substitution
    }

    fn convert_to_hgvs_variant(
        &self,
        record: AdHocRecord,
        protein_accession: &str,
    ) -> Result<HgvsVariant, Error> {
        let AdHocRecord {
            start,
            ins: insertion,
            del: deletion,
            is_frameshift,
        } = &record;

        let ref_protein = self.ref_protein();
        // Without a stop codon, the last amino acid of the reference protein is an ordinary
        // one, and the amino acids after it are unknown.
        let has_stop_codon = ref_protein.ends_with('*');

        // Handle the case of the variant being after the end of the protein sequence.  This can
        // happen when the variant is 5'-to-3' shifted beyond the stop codon.  In this case, we
        // can simply return `ProtLocEdit::NoChange` which will display as `p.=`.  Without a
        // stop codon, the effect is unknown, i.e., `p.?`.
        if *start as usize > ref_protein.len() {
            return Ok(HgvsVariant::ProtVariant {
                accession: Accession::new(protein_accession),
                gene_symbol: None,
                loc_edit: if has_stop_codon {
                    ProtLocEdit::NoChange
                } else {
                    ProtLocEdit::Unknown
                },
            });
        }

        // defaults
        let mut is_dup = false; // assume no dup
        let mut fsext_len = UncertainLengthChange::default(); // fs or ext length
        let mut is_sub = false;
        let mut is_ext = false;
        let mut is_init_met = false;
        let mut is_ambiguous = self.alt_data.is_ambiguous;
        let aa_start;
        let aa_end;
        let mut reference = String::new();
        let mut alternative = String::new();

        if *start == 1 {
            // initial methionine is modified
            // TODO: aa_start/aa_end was in Python code, not needed?
            // aa_start = ProtPos {
            //     aa: "M".to_owned(),
            //     number: 1,
            // };
            // aa_end = aa_start.clone();
            is_init_met = true;
            is_ambiguous = true;
        }

        if insertion.starts_with('*') {
            // stop codon at variant position
            aa_start = Some(ProtPos {
                aa: deletion
                    .chars()
                    .next()
                    .ok_or(Error::DeletionSequenceEmpty)?
                    .to_string(),
                number: *start,
            });
            aa_end = aa_start.clone();
            reference = "".to_string();
            alternative = "*".to_string();
            is_sub = true;
        } else if has_stop_codon && *start as usize == ref_protein.len() {
            // extension
            fsext_len = if self.alt_seq().ends_with('*') {
                // With more than one stop codon, both can go on after the stop codon at which
                // the prediction ends.
                let (insertion, deletion) = if ref_protein.len() < self.ref_seq().len() {
                    (up_to_first_stop(insertion), up_to_first_stop(deletion))
                } else {
                    (insertion.as_str(), deletion.as_str())
                };
                UncertainLengthChange::Known(insertion.len() as i32 - deletion.len() as i32)
            } else {
                UncertainLengthChange::Unknown
            };

            aa_start = Some(ProtPos {
                aa: "*".to_owned(),
                number: *start,
            });
            aa_end = aa_start.clone();

            "".clone_into(&mut reference);
            alternative = insertion
                .chars()
                .next()
                .map(|c| c.to_string())
                .unwrap_or_default();
            is_ext = true;
        } else if *is_frameshift {
            // frameshift; without a stop codon, the CDS can end before the new frame changes
            // a complete codon, and then the effect is unknown
            is_ambiguous |= !has_stop_codon && (insertion.is_empty() || insertion == deletion);
            aa_start = Some(ProtPos {
                aa: deletion
                    .chars()
                    .next()
                    .ok_or(Error::DeletionSequenceEmpty)?
                    .to_string(),
                number: *start,
            });
            aa_end = aa_start.clone();

            "".clone_into(&mut reference);
            alternative = insertion
                .chars()
                .next()
                .map(|c| c.to_string())
                .unwrap_or_default();

            fsext_len = insertion
                .find('*')
                .map(|pos| UncertainLengthChange::Known(pos as i32 + 1))
                .unwrap_or(UncertainLengthChange::Unknown);

            // ALL CASES BELOW HERE: no frameshift - sub/delins/dup
        } else if insertion == deletion {
            // silent
            aa_start = if *start == -1 {
                None
            } else {
                Some(ProtPos {
                    aa: deletion.clone(),
                    number: *start,
                })
            };
            aa_end = aa_start.clone();
        } else if insertion.len() == 1 && deletion.len() == 1 {
            // substitution
            aa_start = Some(ProtPos {
                aa: deletion.clone(),
                number: *start,
            });
            aa_end = aa_start.clone();
            "".clone_into(&mut reference);
            alternative.clone_from(insertion);
            is_sub = true;
        } else if !deletion.is_empty() {
            // delins OR deletion OR stop codon at variant position
            reference.clone_from(deletion);
            let end = start + deletion.len() as i32 - 1;

            aa_start = Some(ProtPos {
                aa: deletion
                    .chars()
                    .next()
                    .ok_or(Error::DeletionSequenceEmpty)?
                    .to_string(),
                number: *start,
            });
            if !insertion.is_empty() {
                // delins
                aa_end = if end > *start {
                    Some(ProtPos {
                        aa: deletion
                            .chars()
                            .last()
                            .ok_or(Error::DeletionSequenceEmpty)?
                            .to_string(),
                        number: end,
                    })
                } else {
                    aa_start.clone()
                };
                alternative.clone_from(insertion);
            } else {
                // deletion OR stop codon at variant position
                if has_stop_codon && deletion.len() as i32 + start == ref_protein.len() as i32 {
                    // stop codon at variant position
                    aa_end = aa_start.clone();
                    reference = "".to_string();
                    alternative = "*".to_string();
                    is_sub = true;
                } else {
                    // deletion
                    aa_end = if end > *start {
                        Some(ProtPos {
                            aa: deletion
                                .chars()
                                .last()
                                .ok_or(Error::DeletionSequenceEmpty)?
                                .to_string(),
                            number: end,
                        })
                    } else {
                        aa_start.clone()
                    };
                    alternative = "".to_string()
                }
            }
        } else if deletion.is_empty() {
            // insertion OR duplication OR extension
            let dup_start;
            (is_dup, dup_start) = self.check_if_ins_is_dup(*start, insertion);

            if is_dup {
                // is duplication
                let dup_end = dup_start + insertion.len() as i32 - 1;
                aa_start = Some(ProtPos {
                    aa: insertion
                        .chars()
                        .next()
                        .ok_or(Error::InsertionSequenceEmpty)?
                        .to_string(),
                    number: dup_start,
                });
                aa_end = Some(ProtPos {
                    aa: insertion
                        .chars()
                        .last()
                        .ok_or(Error::InsertionSequenceEmpty)?
                        .to_string(),
                    number: dup_end,
                });
                reference = "".to_string();
                alternative.clone_from(&reference);
            } else {
                // is non-dup insertion
                let start = std::cmp::max(2, *start as usize) - 1;
                let end = start + 1;

                aa_start = Some(ProtPos {
                    aa: self.ref_seq()[(start - 1)..start].to_owned(),
                    number: start as i32,
                });
                aa_end = Some(ProtPos {
                    aa: self.ref_seq()[(end - 1)..end].to_owned(),
                    number: end as i32,
                });
                reference = "".to_string();
                alternative.clone_from(insertion);
            }
        } else {
            panic!("Unexpected variant: {:?}", record);
        }

        self.create_variant(
            aa_start,
            aa_end,
            &reference,
            &alternative,
            fsext_len,
            is_dup,
            protein_accession,
            is_ambiguous,
            is_sub,
            is_ext,
            false,
            is_init_met,
            *is_frameshift,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn create_variant(
        &self,
        start: Option<ProtPos>,
        end: Option<ProtPos>,
        reference: &str,
        alternative: &str,
        fsext_len: UncertainLengthChange,
        is_dup: bool,
        acc: &str,
        is_ambiguous: bool,
        is_sub: bool,
        is_ext: bool,
        is_no_protein: bool,
        is_init_met: bool,
        is_frameshift: bool,
    ) -> Result<HgvsVariant, Error> {
        assert_eq!(start.is_some(), end.is_some());

        // If the `alternative` contains a stop codon (`*`/`X`) then we have to truncate
        // after it.
        let alternative = if let Some(pos) = alternative.find('*').or_else(|| alternative.find('X'))
        {
            &alternative[..=pos]
        } else {
            alternative
        };

        let loc_edit = if is_init_met {
            ProtLocEdit::InitiationUncertain
        } else if is_ambiguous {
            ProtLocEdit::Unknown
        } else if start.is_none() && !is_no_protein {
            ProtLocEdit::NoChange
        } else if is_no_protein {
            ProtLocEdit::NoProteinUncertain
        } else {
            let interval = ProtInterval {
                start: start.expect("must provide start"),
                end: end.expect("must provide end"),
            };
            // NB: order matters.
            // NB: in the original Python module, it is configurable whether the result
            // of protein prediction is certain or uncertain by a global configuration
            // variable.
            ProtLocEdit::Ordinary {
                edit: Mu::Certain(if is_sub {
                    // cases like Ter525Ter should be Ter525=
                    if reference == alternative
                        || alternative.len() == 1
                            && interval.start.aa == alternative
                            && interval.start.number == interval.end.number
                    {
                        ProteinEdit::Ident
                    } else {
                        ProteinEdit::Subst {
                            alternative: alternative.to_string(),
                        }
                    }
                } else if is_ext {
                    ProteinEdit::Ext {
                        aa_ext: Some(alternative.to_string()),
                        ext_aa: Some("*".to_string()),
                        change: fsext_len,
                    }
                } else if is_frameshift {
                    ProteinEdit::Fs {
                        alternative: Some(alternative.to_string()),
                        terminal: Some("*".to_owned()),
                        length: fsext_len,
                    }
                } else if is_dup {
                    ProteinEdit::Dup
                } else if reference.is_empty() == alternative.is_empty() {
                    if reference.len() > 1 || alternative.len() > 1 {
                        ProteinEdit::DelIns {
                            alternative: alternative.to_string(),
                        }
                    } else {
                        ProteinEdit::Subst {
                            alternative: alternative.to_string(),
                        }
                    }
                } else if alternative.is_empty() {
                    ProteinEdit::Del
                } else {
                    ProteinEdit::Ins {
                        alternative: alternative.to_string(),
                    }
                }),
                loc: Mu::Certain(interval),
            }
        };

        Ok(HgvsVariant::ProtVariant {
            accession: Accession::new(acc),
            gene_symbol: None,
            loc_edit,
        })
    }

    /// Helper to identity an insertion as a duplicate.
    fn check_if_ins_is_dup(&self, start: i32, insertion: &str) -> (bool, i32) {
        if insertion.len() + 1 >= start as usize {
            return (false, -1);
        }
        let dup_cand_start = start as usize - insertion.len() - 1;
        let dup_cand = &self.ref_seq()[dup_cand_start..dup_cand_start + insertion.len()];
        if insertion == dup_cand {
            (true, dup_cand_start as i32 + 1)
        } else {
            (false, -1)
        }
    }
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use crate::{mapper::Error, parser::HgvsVariant, sequences::TranslationTable};

    use super::{AltSeqBuilder, CdsEnd, RefTranscriptData};

    #[test]
    fn build_altseq_fails_for_non_cds_variant() -> Result<(), anyhow::Error> {
        let reference_data = RefTranscriptData {
            transcript_sequence: "AAAATCAAAATGAAAGCGAAAGCGTTTCGCGCGAAATAGGGG".into(),
            aa_sequence: "MKAKAFRAK*".into(),
            cds_start: 10,
            cds_stop: 39,
            protein_accession: "MOCK".into(),
            translation_table: TranslationTable::Standard,
            translation_exceptions: Vec::new(),
            exon_start_positions: Vec::new(),
            exon_end_positions: Vec::new(),
            cds_end: CdsEnd::StopCodon,
        };

        for hgvs in [
            "NM_999999.1:n.15A>G",
            "NC_000001.11:g.15A>G",
            "MOCK:p.Lys2Arg",
        ] {
            let builder = AltSeqBuilder::new(HgvsVariant::from_str(hgvs)?, &reference_data);
            assert!(
                matches!(builder.build_altseq(), Err(Error::ExpectedCdsVariant(_))),
                "{hgvs}"
            );
        }

        Ok(())
    }

    /// Without a stop codon, the CDS goes on after its last base c.27, so there is no 3' UTR.
    /// The transcript is `NM_999987.1` in the `c_to_p` tests of `variant.rs`: the CDS `ATG AAA
    /// GCG AAA GCG TTT CGC GCG AAA` (`MKAKAFRAK`) is followed by `GGGGGCCCCC`.
    #[test]
    fn variant_region_without_stop_codon() -> Result<(), anyhow::Error> {
        use super::VariantLocation;

        let reference_data = RefTranscriptData {
            transcript_sequence: "AAAATCAAACATGAAAGCGAAAGCGTTTCGCGCGAAAGGGGGCCCCC".into(),
            aa_sequence: "MKAKAFRAK".into(),
            cds_start: 11,
            cds_stop: 37,
            protein_accession: "MOCK".into(),
            translation_table: TranslationTable::Standard,
            translation_exceptions: Vec::new(),
            exon_start_positions: vec![-10],
            exon_end_positions: vec![37],
            cds_end: CdsEnd::NoStopCodon,
        };

        for ins_at_boundary_is_intronic in [false, true] {
            for (hgvs_c, expected) in [
                // The same insertion of an A after c.27.
                ("c.27_*1insA", VariantLocation::Exon),
                ("c.27dup", VariantLocation::Exon),
                // The same insertion of a G into the run of five Gs after c.27.
                ("c.27_*1insG", VariantLocation::Exon),
                ("c.*1dup", VariantLocation::Exon),
                ("c.*4_*5insG", VariantLocation::Exon),
                ("c.*1G>A", VariantLocation::Exon),
                ("c.*1_*3del", VariantLocation::Exon),
                ("c.*1_*3inv", VariantLocation::Exon),
                // The copy can go before c.-5, which leaves the CDS intact.
                ("c.-5_*2dup", VariantLocation::FivePrimeUtr),
                ("c.-5_*2del", VariantLocation::WholeGene),
            ] {
                let builder = AltSeqBuilder {
                    ins_at_boundary_is_intronic,
                    ..AltSeqBuilder::new(
                        HgvsVariant::from_str(&format!("NM_999987.1:{hgvs_c}"))?,
                        &reference_data,
                    )
                };
                assert_eq!(
                    builder.get_variant_region(),
                    expected,
                    "{hgvs_c} with ins_at_boundary_is_intronic = {ins_at_boundary_is_intronic}"
                );
            }
        }

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
