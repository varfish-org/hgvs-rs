//! Code for building alternative sequence and convertion to HGVS.p.

use std::{cmp::Ordering, rc::Rc};

use crate::{
    data::interface::Provider,
    parser::{
        Accession, CdsFrom, HgvsVariant, Mu, NaEdit, ProtInterval, ProtLocEdit, ProtPos,
        ProteinEdit, UncertainLengthChange,
    },
    sequences::{revcomp, translate_cds, TranslationTable},
};

#[derive(Debug, Clone)]
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

        // Coding sequences that are not divisable by 3 are not yet supported.
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

pub struct AltTranscriptData {
    /// Transcript nucleotide sequence.
    #[allow(dead_code)]
    transcript_sequence: String,
    /// 1-letter amino acid sequence.
    aa_sequence: String,
    /// 1-based CDS start position.
    #[allow(dead_code)]
    cds_start: i32,
    /// 1-based CDS stop position.
    #[allow(dead_code)]
    cds_stop: i32,
    /// Protein accession number, e.g., `"NP_999999.2"`.
    #[allow(dead_code)]
    protein_accession: String,
    /// Whether this is a frameshift variant.
    is_frameshift: bool,
    /// 1-based AA start index for this variant.
    variant_start_aa: Option<i32>,
    /// Starting position (AA ref index) of the last framewshift which affects the rest of the
    /// sequence, ie.e., not offset by subsequent frameshifts.
    frameshift_start: Option<i32>,
    /// Whether this is a substitution AA variant.
    is_substitution: bool,
    /// Whether variant is "?".
    is_ambiguous: bool,
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
        is_substitution: bool,
        is_ambiguous: bool,
    ) -> Result<Self, anyhow::Error> {
        let transcript_sequence = seq.to_owned();
        let aa_sequence = if !seq.is_empty() {
            let seq_cds = &transcript_sequence[((cds_start - 1) as usize)..];
            let seq_aa = translate_cds(seq_cds, false, "X", TranslationTable::Standard)?;
            let stop_pos = seq_aa[..((cds_stop - cds_start + 1) as usize / 3)]
                .rfind('*')
                .or_else(|| seq_aa.find('*'));
            if let Some(stop_pos) = stop_pos {
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
pub struct AltSeqBuilder {
    /// HgvsVariant::CdsVariant to build the alternative sequence for.
    pub var_c: HgvsVariant,
    /// Information about the transcript reference.
    pub reference_data: RefTranscriptData,
    /// Whether the transcript reference amino acid sequence has multiple stop codons.
    pub ref_has_multiple_stops: bool,
}

impl AltSeqBuilder {
    pub fn new(var_c: HgvsVariant, reference_data: RefTranscriptData) -> Self {
        if !matches!(&var_c, HgvsVariant::CdsVariant { .. }) {
            panic!("Must initialize with HgvsVariant::CdsVariant");
        }
        let ref_has_multiple_stops = reference_data.aa_sequence.matches('*').count() > 1;

        Self {
            var_c,
            reference_data,
            ref_has_multiple_stops,
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
    /// Variant sequence data.
    pub fn build_altseq(&self) -> Result<Vec<AltTranscriptData>, anyhow::Error> {
        // NB: the following common is from the original Python code.
        // Should loop over each allele rather than assume only 1 variant; return a list for now.

        let na_edit = self.var_c.na_edit().expect("Invalid CdsVariant");
        let edit_type = match self.get_variant_region() {
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
                NaEdit::DelRef {..} |
                NaEdit::DelNum { .. } => EditType::WholeGeneDeleted,
                NaEdit::Dup { .. } => {
                    log::warn!("Whole-gene duplication; consequence assumed to not affect protein product");
                    EditType::NotCds
                },
                NaEdit::InvRef { .. } |
                NaEdit::InvNum { .. } => {
                    log::warn!("Whole-gene inversion; consequence assumed to not affected protein product");
                    EditType::NotCds
                },
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
                if loc.start.cds_from == CdsFrom::End && loc.end.cds_from == CdsFrom::End {
                    VariantLocation::ThreePrimeUtr
                } else if loc.start.base < 0 && loc.end.base < 0 {
                    VariantLocation::FivePrimeUtr
                } else if loc.start.base < 0 && loc.end.cds_from == CdsFrom::End {
                    VariantLocation::WholeGene
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
                            } else if pos.offset.unwrap_or(0) <= 0 {
                                start_end.push(
                                    self.reference_data.cds_start as usize - 1 + pos.base as usize
                                        - 1,
                                );
                            } else {
                                start_end.push(
                                    self.reference_data.cds_start as usize - 1 + pos.base as usize,
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

        let seq = self.reference_data.transcript_sequence.to_owned();
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

    fn incorporate_delins(&self) -> Result<AltTranscriptData, anyhow::Error> {
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
        match (reference, alternative) {
            (Some(reference), Some(alternative)) => {
                // delins or SNP
                seq.replace_range(start..end, &alternative);
                if reference.len() == 1 && alternative.len() == 1 {
                    is_substitution = true;
                }
            }
            (Some(_reference), None) => {
                // deletion
                seq.replace_range(start..end, "");
            }
            (None, Some(alternative)) => {
                // insertion
                seq.insert_str(start + 1, &alternative);
            }
            _ => panic!("This should not happen"),
        }

        let is_frameshift = net_base_change % 3 != 0;

        // Use max. of mod 3 value and 1 (in the event that the indel starts in the 5' UTR range).
        let variant_start_aa = std::cmp::max((loc_range_start as f64 / 3.0).ceil() as i32, 1);

        AltTranscriptData::new(
            &seq,
            cds_start,
            cds_stop,
            is_frameshift,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            is_substitution,
            self.ref_has_multiple_stops,
        )
    }

    fn incorporate_dup(&self) -> Result<AltTranscriptData, anyhow::Error> {
        let (seq, cds_start, cds_stop, start, end) = self.setup_incorporate();

        let seq = format!(
            "{}{}{}{}",
            &seq[..start],
            &seq[start..end],
            &seq[start..end],
            &seq[end..]
        );

        let is_frameshift = (end - start) % 3 != 0;

        let loc_range = self
            .var_c
            .loc_range()
            .expect("Could not determine insertion location");
        let variant_start_aa = ((loc_range.end + 1) as f64 / 3.0).ceil() as i32;

        AltTranscriptData::new(
            &seq,
            cds_start,
            cds_stop,
            is_frameshift,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            false,
            self.ref_has_multiple_stops,
        )
    }

    /// Incorporate inversion into sequence.
    fn incorporate_inv(&self) -> Result<AltTranscriptData, anyhow::Error> {
        let (seq, cds_start, cds_stop, start, end) = self.setup_incorporate();

        let seq = format!(
            "{}{}{}",
            &seq[..start],
            revcomp(&seq[start..end]),
            &seq[end..]
        );

        let loc_range = self
            .var_c
            .loc_range()
            .expect("Could not determine insertion location");
        let variant_start_aa =
            std::cmp::max((((loc_range.start + 1) as f64) / 3.0).ceil() as i32, 1);

        AltTranscriptData::new(
            &seq,
            cds_start,
            cds_stop,
            false,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            false,
            self.ref_has_multiple_stops,
        )
    }

    /// Create an alt seq that matches the reference (for non-CDS variants).
    fn create_alt_equals_ref_noncds(&self) -> Result<AltTranscriptData, anyhow::Error> {
        AltTranscriptData::new(
            &self.reference_data.transcript_sequence,
            self.reference_data.cds_start,
            self.reference_data.cds_stop,
            false,
            None,
            &self.reference_data.protein_accession,
            false,
            true,
        )
    }

    /// Create a no-protein result.
    fn create_no_protein(&self) -> Result<AltTranscriptData, anyhow::Error> {
        AltTranscriptData::new(
            "",
            -1,
            -1,
            false,
            None,
            &self.reference_data.protein_accession,
            false,
            false,
        )
    }
}

/// Build `HgvsVariant::ProtVariant` from information about change to transcript.
pub struct AltSeqToHgvsp {
    pub ref_data: RefTranscriptData,
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

impl AltSeqToHgvsp {
    pub fn new(ref_data: RefTranscriptData, alt_data: AltTranscriptData) -> Self {
        Self { ref_data, alt_data }
    }

    /// Compare two amino acid sequences and generate a HGVS tag from the output.
    pub fn build_hgvsp(&self) -> Result<HgvsVariant, anyhow::Error> {
        let mut records = Vec::new();

        if !self.alt_data.is_ambiguous && !self.alt_seq().is_empty() {
            let mut do_delins = true;
            if self.ref_seq() == self.alt_seq() {
                // Silent p. variant.
                let start = self.alt_data.variant_start_aa;
                let record = if let Some(start) = start {
                    if start - 1 < self.ref_seq().len() as i32 {
                        let del = &self
                            .ref_seq()
                            .chars()
                            .nth(start as usize - 1)
                            .unwrap()
                            .to_string();
                        AdHocRecord {
                            start,
                            ins: del.clone(),
                            del: del.clone(),
                            is_frameshift: self.alt_data.is_frameshift,
                        }
                    } else {
                        Default::default()
                    }
                } else {
                    Default::default()
                };
                records.push(record);
                do_delins = false;
            } else if self.is_substitution() && self.ref_seq().len() == self.alt_seq().len() {
                let r = self.ref_seq().chars();
                let a = self.alt_seq().chars();
                let e = r.zip(a).enumerate();
                let mut diff_records = e
                    .filter(|(_i, (r, a))| r != a)
                    .map(|(i, (r, a))| AdHocRecord {
                        start: i as i32 + 1,
                        del: r.to_string(),
                        ins: a.to_string(),
                        is_frameshift: self.alt_data.is_frameshift,
                    })
                    .collect::<Vec<_>>();
                if diff_records.len() == 1 {
                    records.push(diff_records.drain(..).next().unwrap());
                    do_delins = false;
                }
            }

            if do_delins {
                let mut start = self.alt_data.variant_start_aa.unwrap() as usize - 1;
                while self.ref_seq().chars().nth(start) == self.alt_seq().chars().nth(start) {
                    start += 1;
                }
                if self.alt_data.is_frameshift {
                    // Case: frameshifting delins or dup.
                    let deletion = &self.ref_seq()[start..];
                    let insertion = &self.alt_seq()[start..];
                    records.push(AdHocRecord {
                        start: start as i32 + 1,
                        ins: insertion.to_owned(),
                        del: deletion.to_owned(),
                        is_frameshift: self.alt_data.is_frameshift,
                    })
                } else {
                    // Case: non-frameshifting delins or dup.
                    //
                    // Get size diff from diff in ref/alt lengths.
                    let delta = self.alt_seq().len() as isize - self.ref_seq().len() as isize;
                    let offset = start + delta.unsigned_abs();

                    let (insertion, deletion, ref_sub, alt_sub) = match delta.cmp(&0) {
                        // if delta > 0 {
                        Ordering::Greater => (
                            // net insertion
                            self.alt_seq()[start..offset].to_owned(),
                            "".to_string(),
                            self.ref_seq()[start..].to_owned(),
                            self.alt_seq()[offset..].to_owned(),
                        ),
                        Ordering::Less => (
                            // net deletion
                            "".to_string(),
                            self.ref_seq()[start..offset].to_owned(),
                            self.ref_seq()[offset..].to_owned(),
                            self.alt_seq()[start..].to_owned(),
                        ),
                        Ordering::Equal => (
                            // size remains the same
                            "".to_string(),
                            "".to_string(),
                            self.ref_seq()[start..].to_owned(),
                            self.alt_seq()[start..].to_owned(),
                        ),
                    };

                    // From start, get del/ins out to last difference.
                    let r = ref_sub.chars();
                    let a = alt_sub.chars();
                    let diff_indices = r
                        .zip(a)
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
                        (
                            format!("{}{}", deletion, &ref_sub[..max_diff]),
                            format!("{}{}", insertion, &alt_sub[..max_diff]),
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
            Err(anyhow::anyhow!("Got multiple AA variants - not supported!"))
        }
    }

    fn protein_accession(&self) -> &str {
        &self.ref_data.protein_accession
    }

    fn ref_seq(&self) -> &str {
        &self.ref_data.aa_sequence
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
    ) -> Result<HgvsVariant, anyhow::Error> {
        let AdHocRecord {
            start,
            ins: insertion,
            del: deletion,
            is_frameshift,
        } = &record;

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
                aa: deletion.chars().next().unwrap().to_string(),
                number: *start,
            });
            aa_end = aa_start.clone();
            reference = "".to_string();
            alternative = "*".to_string();
            is_sub = true;
        } else if *start as usize == self.ref_seq().len() {
            // extension
            fsext_len = if self.alt_seq().ends_with('*') {
                UncertainLengthChange::Known(insertion.len() as i32 - deletion.len() as i32)
            } else {
                UncertainLengthChange::Unknown
            };

            aa_start = Some(ProtPos {
                aa: "*".to_owned(),
                number: *start,
            });
            aa_end = aa_start.clone();

            reference = "".to_owned();
            alternative = insertion
                .chars()
                .next()
                .map(|c| c.to_string())
                .unwrap_or_default();
            is_ext = true;
        } else if *is_frameshift {
            // frameshift
            aa_start = Some(ProtPos {
                aa: deletion.chars().next().unwrap().to_string(),
                number: *start,
            });
            aa_end = aa_start.clone();

            reference = "".to_owned();
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
            reference = "".to_owned();
            alternative = insertion.clone();
            is_sub = true;
        } else if !deletion.is_empty() {
            // delins OR deletion OR stop codon at variant position
            reference = deletion.clone();
            let end = start + deletion.len() as i32 - 1;

            aa_start = Some(ProtPos {
                aa: deletion.chars().next().unwrap().to_string(),
                number: *start,
            });
            if !insertion.is_empty() {
                // delins
                aa_end = if end > *start {
                    Some(ProtPos {
                        aa: deletion.chars().last().unwrap().to_string(),
                        number: end,
                    })
                } else {
                    aa_start.clone()
                };
                alternative = insertion.clone();
            } else {
                // deletion OR stop codon at variant position
                if deletion.len() as i32 + start == self.ref_seq().len() as i32 {
                    // stop codon at variant position
                    aa_end = aa_start.clone();
                    reference = "".to_string();
                    alternative = "*".to_string();
                    is_sub = true;
                } else {
                    // deletion
                    aa_end = if end > *start {
                        Some(ProtPos {
                            aa: deletion.chars().last().unwrap().to_string(),
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
                    aa: insertion.chars().next().unwrap().to_string(),
                    number: dup_start,
                });
                aa_end = Some(ProtPos {
                    aa: insertion.chars().last().unwrap().to_string(),
                    number: dup_end,
                });
                reference = "".to_string();
                alternative = reference.clone();
            } else {
                // is non-dup insertion
                let start = start - 1;
                let end = start + 1;

                aa_start = Some(ProtPos {
                    aa: self
                        .ref_seq()
                        .chars()
                        .nth(start as usize - 1)
                        .unwrap()
                        .to_string(),
                    number: start,
                });
                aa_end = Some(ProtPos {
                    aa: self
                        .ref_seq()
                        .chars()
                        .nth(end as usize - 1)
                        .unwrap()
                        .to_string(),
                    number: end,
                });
                reference = "".to_string();
                alternative = insertion.clone();
            }
        } else {
            panic!("Unexpected variant: {:?}", &record);
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
    ) -> Result<HgvsVariant, anyhow::Error> {
        assert!(start.is_some() == end.is_some());

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
                loc: Mu::Certain(interval),
                edit: Mu::Certain(if is_sub {
                    ProteinEdit::Subst {
                        alternative: alternative.to_string(),
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
        let dup_cand_start = start as usize - insertion.len() - 1;
        let dup_cand = &self.ref_seq()[dup_cand_start..dup_cand_start + insertion.len()];
        if insertion == dup_cand {
            (true, dup_cand_start as i32 + 1)
        } else {
            (false, -1)
        }
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
