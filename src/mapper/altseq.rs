//! Code for building alternative sequence and convertion to HGVS.p.

use std::{cell::RefCell, cmp::Ordering, sync::Arc};

use quick_cache::{unsync::Cache, Equivalent};

use crate::{
    data::interface::{Provider, TranslationException},
    mapper::error::Error,
    parser::{
        Accession, CdsFrom, CdsPos, HgvsVariant, Mu, NaEdit, ProtInterval, ProtLocEdit, ProtPos,
        ProteinEdit, UncertainLengthChange,
    },
    sequences::{revcomp, translate_cds, translate_cds_with_exceptions, TranslationTable},
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

/// The `A`s that the poly(A) tail adds to complete the stop codon of `cds`, as biocommons/hgvs
/// adds them: `AA` after a final `T` and `A` after a final `TA`, with the vertebrate
/// mitochondrial translation table only.
fn poly_a(translation_table: TranslationTable, cds: &str) -> &'static str {
    match translation_table {
        TranslationTable::VertebrateMitochondrial if cds.len() % 3 == 1 && cds.ends_with('T') => {
            "AA"
        }
        TranslationTable::VertebrateMitochondrial if cds.len() % 3 == 2 && cds.ends_with("TA") => {
            "A"
        }
        _ => "",
    }
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

/// Cache key of `ref_transcript_data_shared`.
#[derive(Hash, PartialEq, Eq)]
struct RefTranscriptDataKey {
    data_version: String,
    schema_version: String,
    tx_ac: String,
    pro_ac: Option<String>,
}

/// Borrowed form of `RefTranscriptDataKey`, so that lookups need no allocation.
///
/// It must hash like `RefTranscriptDataKey`, hence the same fields in the same order.
#[derive(Hash)]
struct RefTranscriptDataKeyRef<'a> {
    data_version: &'a str,
    schema_version: &'a str,
    tx_ac: &'a str,
    pro_ac: Option<&'a str>,
}

impl Equivalent<RefTranscriptDataKey> for RefTranscriptDataKeyRef<'_> {
    fn equivalent(&self, key: &RefTranscriptDataKey) -> bool {
        self.data_version == key.data_version
            && self.schema_version == key.schema_version
            && self.tx_ac == key.tx_ac
            && self.pro_ac == key.pro_ac.as_deref()
    }
}

impl RefTranscriptDataKeyRef<'_> {
    fn to_owned(&self) -> RefTranscriptDataKey {
        RefTranscriptDataKey {
            data_version: self.data_version.to_string(),
            schema_version: self.schema_version.to_string(),
            tx_ac: self.tx_ac.to_string(),
            pro_ac: self.pro_ac.map(|s| s.to_string()),
        }
    }
}

thread_local! {
    /// Per thread, see `build_alignment_mapper_cached` in `variant.rs` for why.
    static REF_TRANSCRIPT_DATA_CACHE: RefCell<Cache<RefTranscriptDataKey, Result<Arc<RefTranscriptData>, Error>>> =
        RefCell::new(Cache::new(1000));
}

/// A cached version of `RefTranscriptData::new`.
pub fn ref_transcript_data_cached(
    provider: Arc<dyn Provider + Send + Sync>,
    tx_ac: &str,
    pro_ac: Option<&str>,
) -> Result<RefTranscriptData, Error> {
    ref_transcript_data_shared(&provider, tx_ac, pro_ac).map(|data| data.as_ref().clone())
}

/// Like `ref_transcript_data_cached`, but returns the cached value instead of a copy.
pub(crate) fn ref_transcript_data_shared(
    provider: &Arc<dyn Provider + Send + Sync>,
    tx_ac: &str,
    pro_ac: Option<&str>,
) -> Result<Arc<RefTranscriptData>, Error> {
    let key = RefTranscriptDataKeyRef {
        data_version: provider.data_version(),
        schema_version: provider.schema_version(),
        tx_ac,
        pro_ac,
    };
    if let Some(hit) = REF_TRANSCRIPT_DATA_CACHE.with_borrow(|c| c.get(&key).cloned()) {
        return hit;
    }
    let key = key.to_owned();
    let result = RefTranscriptData::new(provider.clone(), tx_ac, pro_ac).map(Arc::new);
    REF_TRANSCRIPT_DATA_CACHE.with_borrow_mut(|c| c.insert(key, result.clone()));
    result
}

impl RefTranscriptData {
    /// The `A`s that complete the stop codon of a mitochondrial CDS, see `new()`.
    pub(crate) fn poly_a(&self) -> &'static str {
        self.transcript_sequence
            .get((self.cds_start - 1) as usize..self.cds_stop as usize)
            .map_or("", |cds| poly_a(self.translation_table, cds))
    }

    /// Construct new instance fetching data from the provider..
    ///
    /// If the length of the CDS is not a multiple of 3, the protein ends with the last
    /// complete codon, and `cds_end` holds the bases after it.  Unlike biocommons/hgvs, which
    /// rejects such a CDS, this returns an error only if a stop codon is in frame.
    ///
    /// As in biocommons/hgvs, the poly(A) tail completes the stop codon of a CDS with the
    /// vertebrate mitochondrial translation table that ends with `T` or `TA`, see `poly_a()`.
    /// The missing `A`s go into the translated bases only, so `transcript_sequence`,
    /// `cds_stop` and the c. positions stay as they are.
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
        let cds = &transcript_sequence[((cds_start - 1) as usize)..(cds_stop as usize)];
        let tx_seq_to_translate = &format!("{cds}{}", poly_a(tx_info.translation_table, cds));
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
            None,
        )
    }

    /// Like `new`, but takes ownership of `seq`, and sets the amino acids of
    /// `translation_exceptions` at their 0-based codons of the CDS.  `poly_a` holds the `A`s
    /// that complete a mitochondrial stop codon and their 0-based index in `seq`.  They go into
    /// the translated bases only.
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
        poly_a: Option<(usize, &str)>,
    ) -> Result<Self, Error> {
        let transcript_sequence = seq;
        let aa_sequence = if !transcript_sequence.is_empty() {
            let cds_index = (cds_start - 1) as usize;
            let with_poly_a;
            // In case of SEPHS2 / HGNC:19686, the last amino acid is both a selenocysteine
            // and a stop codon.
            // We handle this by explicitly truncating the sequence at the Sec + stop codon.
            // This heuristic may not always be correct;
            // alternatively/additionally, we could check `protein_accession` for known cases.
            let seq_cds = if let Some((index, bases)) =
                poly_a.filter(|(index, _)| (cds_index..=transcript_sequence.len()).contains(index))
            {
                with_poly_a = format!(
                    "{}{bases}{}",
                    &transcript_sequence[cds_index..index],
                    &transcript_sequence[index..]
                );
                with_poly_a.as_str()
            } else if translation_table == TranslationTable::Selenocysteine
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
            // With a stop codon, the alt protein ends after the last stop codon up to
            // `cds_stop`, or else at the next one after it.  Without a stop codon, e.g., in
            // ENST00000420031.2, the amino acids after `cds_stop` are unknown, so the alt
            // protein ends at `cds_stop`, which the caller moves by the length change of the
            // variant.  So it lines up with `ref_aa_sequence`, also if it has a new stop codon
            // or if a stop codon follows `cds_stop`.
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

/// The region of a transcript that a c. variant changes, see
/// `AltSeqBuilder::variant_location()`.
///
/// As in biocommons/hgvs, the positions in the description decide the region, with one
/// exception: an insertion, dup, deletion, delins or inversion whose alt transcript sequence
/// keeps the start or stop codon in place is a UTR variant, see `FivePrimeUtr` and
/// `ThreePrimeUtr`.  `c_to_p` builds the altered protein only for `Cds`.  It gives `p.?` for a UTR or an intron, and `p.0?` for
/// the deletion of the whole gene.  It also gives `p.?` for a `Cds` variant after a CDS without
/// a stop codon, see there.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum VariantLocation {
    /// The variant changes bases of the CDS and has no intron position, e.g., `c.4A>G`.  It may
    /// also change bases of a UTR, e.g., `c.-2_2del`.
    ///
    /// Without `ins_at_boundary_is_intronic`, this also includes the insertions and
    /// duplications next to an exon boundary within the CDS, e.g., `c.16-1_16insT` and
    /// `c.16dup` if c.16 starts an exon.  The splice site stays intact, so the inserted bases
    /// extend the exon.
    ///
    /// Without a stop codon, the CDS goes on after its last base `c.N`, so there is no 3' UTR.
    /// This then includes the variants that end after `c.N`, e.g., `c.*1G>A` or
    /// `c.(N-1)_*1del`, and those whose alt sequence keeps the CDS up to `c.N`, e.g., `c.Ndup`.
    /// The amino acids after `c.N` are unknown, so `c_to_p` gives `p.?`.
    Cds,
    /// The variant has an intron position, e.g., `c.15+5del`, and lies neither in a UTR nor
    /// across the whole CDS.  This includes variants that also change bases of the CDS, e.g.,
    /// `c.15_15+3del`.
    ///
    /// With `ins_at_boundary_is_intronic`, this includes the insertions and duplications at an
    /// exon boundary within the CDS, e.g., `c.16-1_16insT` and `c.16dup` if c.16 starts an
    /// exon.  The inserted bases may go into the intron.
    Intron,
    /// Both ends of the variant lie before c.1, e.g., `c.-5C>G` or `c.-1+5del`.
    ///
    /// This also includes an insertion, dup, deletion, delins or inversion whose alt transcript
    /// sequence ends with the reference from c.1 on.  Then the same change fits in front of c.1,
    /// and the start codon stays where it is.  Examples are `c.-1_1insGCC`, `c.1dup` if c.-1 is
    /// not an A, `c.1del` if c.-1 is an A, and every dup from before c.1 into the CDS, such as
    /// `c.-2_2dup`, also if it starts before the transcript.  The inserted bases can still form
    /// an upstream ATG, so the effect on the protein is unknown.
    FivePrimeUtr,
    /// Both ends of the variant lie after the last base `c.N` of the CDS, e.g., `c.*3G>T`, or
    /// it is an insertion or duplication whose end does, e.g., `c.N_*1insA` or `c.-5_*2dup`.
    ///
    /// This also includes an insertion, dup, deletion, delins or inversion whose alt transcript
    /// sequence starts with the reference up to `c.N`.  Then the same change fits after `c.N`,
    /// and the stop codon stays where it is.  Examples are `c.Ndup`, `c.(N-2)_Ndup`, and `c.Ndel`
    /// if c.*1 is the same base.  Only a CDS with a stop codon has a 3' UTR, see `Cds`.
    ThreePrimeUtr,
    /// The variant spans the CDS from the 5' UTR to the 3' UTR, e.g., `c.-10_*10del`.  A
    /// duplication is `ThreePrimeUtr` instead with a stop codon, and `FivePrimeUtr` without one,
    /// as its copy can go before its first base.
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
    /// Whether the transcript reference amino acid sequence has multiple stop codons.  As in
    /// biocommons/hgvs, `c_to_p` then gives `p.?` for every variant.
    pub ref_has_multiple_stops: bool,
    /// 0-based index of the first stop codon in the reference amino acid sequence.
    #[deprecated(note = "unused, `c_to_p` gives `p.?` for a reference with several stop codons")]
    pub first_stop_pos: Option<usize>,
    /// See `variant::Config::ins_at_boundary_is_intronic`; `true` by default.
    pub ins_at_boundary_is_intronic: bool,
}

impl<'a> AltSeqBuilder<'a> {
    /// A builder for `var_c` on `reference_data`, with the defaults of `variant::Config`.
    #[allow(deprecated)] // first_stop_pos
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
        let edit_type = match self.variant_location()? {
            // Without a stop codon, the amino acids after the CDS are unknown.
            VariantLocation::Cds if self.after_cds_without_stop() => EditType::NotCds,
            VariantLocation::Cds => match na_edit {
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

    /// The region of the transcript that `var_c` changes.
    ///
    /// `build_altseq()` and thus `c_to_p` use this region.  It depends on
    /// `ins_at_boundary_is_intronic`, so it agrees with `c_to_p` only if this field has the
    /// value of `variant::Config::ins_at_boundary_is_intronic` of the mapper.
    /// `AltSeqBuilder::new()` sets it to `true`, the default of the `Config`.  With
    /// `variant::Config::shift_over_boundary`, `c_to_p` can build the protein from a shift of
    /// `var_c` into another region, which this method does not apply.  Without
    /// `ins_at_boundary_is_intronic`, the region of a dup of intron bases can also depend on
    /// these bases, which `c_to_p` fills in first.
    ///
    /// # Errors
    ///
    /// `Error::ExpectedCdsVariant` if `var_c` is not a `HgvsVariant::CdsVariant`.
    pub fn variant_location(&self) -> Result<VariantLocation, Error> {
        let location = match &self.var_c {
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
                    VariantLocation::Cds
                } else if loc.start.offset.is_some() || loc.end.offset.is_some() {
                    // Leave out anything intronic for now.
                    VariantLocation::Intron
                } else {
                    // Anything else that contains an exon.
                    VariantLocation::Cds
                }
            }
            _ => return Err(Error::ExpectedCdsVariant(self.var_c.to_string())),
        };
        Ok(location)
    }

    /// The UTR of a variant that leaves the CDS intact, or `None`.
    ///
    /// An ins, dup, del, delins, or inv leaves the CDS intact if the same change fits before
    /// c.1 or after the last base c.N of the CDS, i.e., if its alt transcript sequence keeps the
    /// reference from c.1 on or up to c.N.  This depends on the sequence, not on the
    /// description.  For example, if c.-1 is not an A, inserting the copy of `c.1dup` in front
    /// of c.1 gives the same sequence.  If c.-1 is an A, deleting it gives the same sequence as
    /// `c.1del`.  The copy of `c.(N-2)_Ndup` can go after c.N.  The inserted bases can still
    /// form an upstream ATG, so the effect on the protein is unknown, not none.
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

    /// Whether the CDS has no stop codon and the variant reaches after its last base c.N, i.e.,
    /// it ends with a c.* position, or it can go there, e.g., `c.Ndup`.  The CDS goes on after
    /// c.N, but its amino acids there are unknown.
    fn after_cds_without_stop(&self) -> bool {
        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return false;
        };
        self.reference_data.cds_end != CdsEnd::StopCodon
            && (loc_edit.loc.inner().end.cds_from == CdsFrom::End
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
            // `variant_location()`.
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

    /// The 0-based index in the alt transcript right after the first stop codon that the
    /// inserted bases `insert_seq` form in the reading frame, or `None`.  `insert_seq_idx` is the
    /// index of their first base.  Translation then ends within the inserted bases, so the
    /// variant is not a frameshift, and the bases after the stop codon do not count.
    ///
    /// Port of `_insert_stop_seq_end()` of biocommons/hgvs (biocommons/hgvs#830).  As there, a
    /// codon that starts before the inserted bases does not count.
    fn insert_stop_seq_end(
        &self,
        insert_seq: &str,
        insert_seq_idx: usize,
    ) -> Result<Option<usize>, Error> {
        let Some(frame_offset) = usize::try_from(self.reference_data.cds_start - 1)
            .ok()
            .and_then(|cds_start| insert_seq_idx.checked_sub(cds_start))
            .map(|offset| offset % 3)
        else {
            return Ok(None);
        };
        let framed = format!("{}{insert_seq}", "N".repeat(frame_offset));
        let ins_aa = translate_cds(&framed, false, "X", self.reference_data.translation_table)?;
        Ok(ins_aa
            .find('*')
            .map(|stop_aa_idx| insert_seq_idx + (stop_aa_idx + 1) * 3 - frame_offset))
    }

    /// The `A`s that complete a mitochondrial stop codon and their 0-based index in `alt_seq`,
    /// see `RefTranscriptData::new()`.  The change of `alt_seq` ends at the 0-based index `end`
    /// of the reference.  If it ends right after the CDS, e.g., an insertion there, the `A`s
    /// follow it.  `None` if there are none, or if the change reaches into the 3' UTR, e.g.,
    /// `c.16_*1del`, which removes the base that the `A`s follow.
    fn poly_a_in(&self, alt_seq: &str, end: usize) -> Option<(usize, &'static str)> {
        let bases = self.reference_data.poly_a();
        let cds_end = self.reference_data.cds_stop as usize;
        if bases.is_empty() || end > cds_end {
            return None;
        }
        // The change moves the end of the CDS by its length change.
        let index =
            (cds_end + alt_seq.len()).checked_sub(self.reference_data.transcript_sequence.len())?;
        Some((index, bases))
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

    /// Incorporate a substitution, deletion, insertion or delins into the transcript sequence.
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
        // An insertion changes no base and goes in at `start + 1`.
        let change_end = if reference.is_some() { end } else { start + 1 };
        let alt_len = alternative.as_ref().map(|alt| alt.len()).unwrap_or(0) as i32;
        let net_base_change = alt_len - ref_len;
        let mut cds_stop = cds_stop + net_base_change;
        let insert_stop = match (&reference, &alternative) {
            (None, Some(alternative)) => self
                .insert_stop_seq_end(alternative, start + 1)?
                .map(|stop_end| (stop_end, start + 1 + alternative.len())),
            _ => None,
        };

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

        let mut is_frameshift = net_base_change % 3 != 0;
        if let Some((stop_end, alt_end)) = insert_stop {
            is_frameshift = false;
            seq.replace_range(stop_end..alt_end, "");
            cds_stop -= (alt_end - stop_end) as i32;
        }

        // Use max. of mod 3 value and 1 (in the event that the indel starts in the 5' UTR range).
        let variant_start_aa = std::cmp::max((loc_range_start as f64 / 3.0).ceil() as i32, 1);
        let translation_exceptions = self.alt_translation_exceptions(&seq);
        let poly_a = self.poly_a_in(&seq, change_end);

        AltTranscriptData::new_owned(
            seq,
            cds_start,
            cds_stop,
            is_frameshift,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            is_substitution,
            self.ref_has_multiple_stops,
            self.reference_data.translation_table,
            &translation_exceptions,
            poly_a,
        )
    }

    /// Incorporate a dup into the transcript sequence: its copy goes right after the duplicated
    /// bases.
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
        let mut seq = format!("{}{}{}", &seq[..end], dup_seq, &seq[end..]);

        let mut is_frameshift = dup_seq.len() % 3 != 0;
        // Without a stop codon, the alt protein ends at `cds_stop`, see `AltTranscriptData`.  So
        // `cds_stop` moves with the copy, as in `incorporate_delins()`, and the alt protein
        // keeps the last codons of the CDS.  With a stop codon, biocommons/hgvs keeps `cds_stop`.
        let mut cds_stop = if self.reference_data.aa_sequence.contains('*') {
            cds_stop
        } else {
            cds_stop + dup_seq.len() as i32
        };
        if let Some(stop_end) = self.insert_stop_seq_end(dup_seq, end)? {
            is_frameshift = false;
            let dup_end = end + dup_seq.len();
            seq.replace_range(stop_end..dup_end, "");
            cds_stop -= (dup_end - stop_end) as i32;
        }

        let loc_end = match &self.var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => loc_edit.loc.inner().end.base,
            _ => panic!("can only work on CDS variants"),
        };
        let variant_start_aa = ((loc_end + 1) as f64 / 3.0).ceil() as i32;
        let translation_exceptions = self.alt_translation_exceptions(&seq);
        // The copy goes in at `end`.
        let poly_a = self.poly_a_in(&seq, end);

        AltTranscriptData::new_owned(
            seq,
            cds_start,
            cds_stop,
            is_frameshift,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            false,
            self.ref_has_multiple_stops,
            self.reference_data.translation_table,
            &translation_exceptions,
            poly_a,
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
        let poly_a = self.poly_a_in(&seq, end);

        AltTranscriptData::new_owned(
            seq,
            cds_start,
            cds_stop,
            false,
            Some(variant_start_aa),
            &self.reference_data.protein_accession,
            &self.reference_data.aa_sequence,
            false,
            self.ref_has_multiple_stops,
            self.reference_data.translation_table,
            &translation_exceptions,
            poly_a,
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
            None,
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
            None,
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

            if ref_seq == alt_seq {
                // Silent p. variant.
                if let Some(start) = self.alt_data.variant_start_aa {
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
    ) -> Result<HgvsVariant, Error> {
        let AdHocRecord {
            start,
            ins: insertion,
            del: deletion,
            is_frameshift,
        } = &record;

        // Without a stop codon, the last amino acid of the reference protein is an ordinary
        // one, and the amino acids after it are unknown.
        let has_stop_codon = self.ref_data.cds_end == CdsEnd::StopCodon;

        // Handle the case of the variant being after the end of the protein sequence.  This can
        // happen when the variant is 5'-to-3' shifted beyond the stop codon.  In this case, we
        // can simply return `ProtLocEdit::NoChange` which will display as `p.=`.  Without a
        // stop codon, the effect is unknown, i.e., `p.?`.
        if *start as usize > self.ref_seq().len() {
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
        } else if has_stop_codon && *start as usize == self.ref_seq().len() {
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
                if has_stop_codon && deletion.len() as i32 + start == self.ref_seq().len() as i32 {
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

    /// Cache lookups only hit if the borrowed key hashes like the owned one.
    #[test]
    fn ref_transcript_data_key_ref_matches_key() {
        use quick_cache::Equivalent;
        use std::hash::{BuildHasher, RandomState};

        let state = RandomState::new();
        for pro_ac in [None, Some("NP_001628.1")] {
            let key_ref = super::RefTranscriptDataKeyRef {
                data_version: "uta_20210129",
                schema_version: "1.1",
                tx_ac: "NM_001637.3",
                pro_ac,
            };
            let key = key_ref.to_owned();

            assert_eq!(state.hash_one(&key_ref), state.hash_one(&key));
            assert!(key_ref.equivalent(&key));
        }
    }

    #[test]
    fn non_cds_variant_fails() -> Result<(), anyhow::Error> {
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
            assert!(
                matches!(
                    builder.variant_location(),
                    Err(Error::ExpectedCdsVariant(_))
                ),
                "{hgvs}"
            );
        }

        Ok(())
    }

    /// Without a stop codon, the CDS goes on after its last base c.27, so there is no 3' UTR.
    /// The CDS is `ATG AAA GCG AAA GCG TTT CGC GCG AAA` (`MKAKAFRAK`), followed by
    /// `GGGGGCCCCC`.
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
                ("c.27_*1insA", VariantLocation::Cds),
                ("c.27dup", VariantLocation::Cds),
                // The same insertion of a G into the run of five Gs after c.27.
                ("c.27_*1insG", VariantLocation::Cds),
                ("c.*1dup", VariantLocation::Cds),
                ("c.*4_*5insG", VariantLocation::Cds),
                ("c.*1G>A", VariantLocation::Cds),
                ("c.*1_*3del", VariantLocation::Cds),
                ("c.*1_*3inv", VariantLocation::Cds),
                ("c.26_*1del", VariantLocation::Cds),
                ("c.27_*1delinsTC", VariantLocation::Cds),
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
                    builder.variant_location()?,
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
