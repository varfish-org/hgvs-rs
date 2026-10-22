//! Code for building alternative sequence and convertion to HGVS.p.

use std::{cell::RefCell, cmp::Ordering, sync::Arc};

use quick_cache::{unsync::Cache, Equivalent};

use crate::{
    data::interface::{Provider, TranslationException},
    mapper::error::Error,
    parser::{
        Accession, CdsFrom, CdsInterval, CdsPos, HgvsVariant, Mu, NaEdit, ProtInterval,
        ProtLocEdit, ProtPos, ProteinEdit, UncertainLengthChange,
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

/// The translation of `seq` by `translate_cds_with_exceptions(seq, false, "X", ..)` of its first
/// `codons` codons, the original length, and, if they have no stop codon, of the codons up to
/// the first one behind them.  With a stop codon in the reference, the alt protein of
/// `AltTranscriptData::new_owned()` ends with the last stop codon among the first `codons`
/// codons, or else with the first one behind them, so it needs no more.  This spares most of
/// the 3' UTR.  Without a stop codon in the reference, the alt protein ends after `codons`
/// codons, but this still goes on to the first stop codon behind them.
fn translate_until_stop(
    seq: &str,
    codons: usize,
    translation_table: TranslationTable,
    exceptions: &[(usize, char)],
) -> Result<String, Error> {
    // The codons to translate at a time behind the first `codons`.
    const CHUNK: usize = 64;
    let mut end = codons.saturating_mul(3).min(seq.len());
    let mut aa =
        translate_cds_with_exceptions(&seq[..end], false, "X", translation_table, exceptions)?;
    let mut stop = aa.contains('*');
    while !stop && end < seq.len() {
        let start = end;
        end = (start + 3 * CHUNK).min(seq.len());
        let first_codon = start / 3;
        let exceptions: Vec<_> = (exceptions.iter())
            .filter_map(|&(codon, amino_acid)| Some((codon.checked_sub(first_codon)?, amino_acid)))
            .collect();
        let part = translate_cds_with_exceptions(
            &seq[start..end],
            false,
            "X",
            translation_table,
            &exceptions,
        )?;
        stop = part.contains('*');
        aa.push_str(&part);
    }
    Ok(aa)
}

/// The length of the common prefix of `a` and `b`.
fn common_prefix_len(a: &[u8], b: &[u8]) -> usize {
    // Comparing chunks is much faster than comparing single bytes.
    let chunks: usize = (a.chunks(64).zip(b.chunks(64)))
        .take_while(|(a, b)| a == b)
        .map(|(a, _)| a.len())
        .sum();
    chunks
        + (a[chunks..].iter().zip(&b[chunks..]))
            .take_while(|(a, b)| a == b)
            .count()
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

/// Like `ref_transcript_data_cached`, but returns the cached value instead of a copy.  As there,
/// the cache is per thread.
pub fn ref_transcript_data_shared(
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
            false,
        )
    }

    /// Like `new`, but takes ownership of `seq`, and sets the amino acids of
    /// `translation_exceptions` at their 0-based codons of the CDS.  `poly_a` holds the `A`s
    /// that complete a mitochondrial stop codon and their 0-based index in `seq`.  They go into
    /// the translated bases only.  With `until_stop`, it translates only the codons up to the
    /// original length and, if they have no stop codon, those up to the first one, see
    /// `translate_until_stop()`.  This gives the same `aa_sequence`, but an invalid base behind
    /// the translated codons no longer fails the translation.  Without a stop codon in the
    /// reference, `aa_sequence` can end before the translated codons, and an invalid base in
    /// between still fails the translation.
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
        until_stop: bool,
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

            let orig_aa_len = (cds_stop - cds_start + 1) as usize / 3;
            let seq_aa = if variant_start_aa.is_none() {
                ref_aa_sequence.to_owned()
            } else if until_stop {
                translate_until_stop(
                    seq_cds,
                    orig_aa_len,
                    translation_table,
                    translation_exceptions,
                )?
            } else {
                translate_cds_with_exceptions(
                    seq_cds,
                    false,
                    "X",
                    translation_table,
                    translation_exceptions,
                )?
            };
            // With a stop codon, the alt protein ends after the last stop codon up to
            // `cds_stop`, or else at the next one after it.  Without a stop codon, e.g., in
            // ENST00000420031.2, the amino acids after `cds_stop` are unknown, so the alt
            // protein ends at `cds_stop`, which the caller moves by the length change of the
            // variant.  So it lines up with `ref_aa_sequence`, also if it has a new stop codon
            // or if a stop codon follows `cds_stop`.
            let orig_aa_len = std::cmp::min(orig_aa_len, seq_aa.len());
            // The stop codon is usually last, which spares the scan.
            if !(ref_aa_sequence.ends_with('*') || ref_aa_sequence.contains('*')) {
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
/// `ThreePrimeUtr`.  The option `ins_at_boundary_by_sequence` adds a second exception, see
/// `Cds`.  `c_to_p` builds the altered protein only for `Cds`.  It gives `p.?` for a
/// UTR or an intron, and `p.0?` for the deletion of the whole gene.  It also gives `p.?` for a
/// `Cds` variant after a CDS without a stop codon, see there.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum VariantLocation {
    /// The variant changes bases of the CDS and has no intron position, e.g., `c.4A>G`.  It may
    /// also change bases of a UTR, e.g., `c.-2_2del`.
    ///
    /// It also includes the insertions and duplications next to an exon boundary within the
    /// CDS whose bases join the exon, e.g., `c.16-1_16insT` and `c.16dup` if c.16 starts an
    /// exon.  `ins_at_boundary_is_intronic` decides this for all of them.  With
    /// `ins_at_boundary_by_sequence`, the alt sequence decides it for each one, if the intron
    /// bases next to the boundary are known: the bases join the exon if the splice-site
    /// dinucleotide moves away from the exon edge.
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
    /// It also includes the insertions and duplications at an exon boundary within the CDS
    /// whose bases go into the intron, see `Cds`.
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

/// Whether the alt sequence that replaces `seq[start..end]` with `alt` ends with
/// `seq[cds_start..]`, and whether it starts with `seq[..cds_end]`, see
/// `AltSeqBuilder::placements_outside_cds()`.  `None` if it changes no base, e.g., the inversion
/// of `GC`, or if an index lies outside `seq`.
///
/// If `start == end` and `alt` is empty, there are no bases to compare, and the result comes
/// from `start` alone, without a bounds check: whether it lies at or before `cds_start`, and
/// whether it lies at or behind `cds_end`.
///
/// The alt sequence is `seq[..start]`, `alt`, and `seq[end..]`, so it need not be built: the
/// comparison with `seq[cds_start..]` skips `seq[end..]`, and the one with `seq[..cds_end]`
/// skips `seq[..start]`.
fn keeps_ref_suffix_and_prefix(
    seq: &[u8],
    start: usize,
    end: usize,
    alt: &[u8],
    cds_start: usize,
    cds_end: usize,
) -> Option<(bool, bool)> {
    if start == end && alt.is_empty() {
        // No bases to compare, so the gap at `start` decides.
        return Some((start <= cds_start, start >= cds_end));
    }
    if start > seq.len() || end > seq.len() {
        return None;
    }
    let same_length = start + alt.len() == end;
    if same_length && alt == &seq[start..end] {
        return None;
    }
    if cds_start > seq.len() || cds_end > seq.len() {
        return None;
    }
    // The alt sequence from `lo` on must be `seq[cds_start..]`.  Behind `alt`, it is
    // `seq[end..]`, which lines up, so only the bases up to the end of `alt` count.
    let keeps_suffix = if end <= cds_start {
        true
    } else if let Some(lo) = (cds_start + start + alt.len()).checked_sub(end) {
        let cds_part = &seq[cds_start..end];
        match lo.checked_sub(start) {
            Some(skip) => alt[skip..] == *cds_part,
            None => {
                let (head, tail) = cds_part.split_at(start - lo);
                (same_length || seq[lo..start] == *head) && alt == tail
            }
        }
    } else {
        // The alt sequence is shorter than `seq[cds_start..]`.
        false
    };
    // The alt sequence up to `start` is `seq[..start]`.
    let keeps_prefix = if cds_end <= start {
        true
    } else if let Some(behind_alt) = (cds_end - start).checked_sub(alt.len()) {
        let alt_end = start + alt.len();
        alt == &seq[start..alt_end]
            && (same_length || seq.get(end..end + behind_alt) == Some(&seq[alt_end..cds_end]))
    } else {
        alt[..cds_end - start] == seq[start..cds_end]
    };
    Some((keeps_suffix, keeps_prefix))
}

/// Utility enum for edit type.
enum EditType {
    /// A substitution or delins: the interval is replaced by the given bases.
    Delins {
        alternative: String,
    },
    /// A deletion of the interval.
    Del,
    /// An insertion of the given bases after the start of the interval.
    Ins {
        alternative: String,
    },
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
    /// See `variant::Config::ins_at_boundary_by_sequence`; `false` by default.
    pub ins_at_boundary_by_sequence: bool,
    /// For an insertion or dup at an exon boundary, the two intron bases next to it in the
    /// orientation of the transcript: `+1` and `+2` after the last exon base, or `-2` and `-1`
    /// before the first one.  `ins_at_boundary_by_sequence` needs them, unless a dup of intron
    /// bases tells them.  `None` by default; `assembly::Mapper::alt_seq_builder()` sets them.
    pub boundary_intron_bases: Option<String>,
}

/// The side of the exon at whose boundary an insertion puts its bases.
enum ExonEdge {
    /// After the last exon base, in front of the donor site.
    End,
    /// Before the first exon base, behind the acceptor site.
    Start,
}

/// Key that orders c. positions along the transcript: 5' UTR and CDS positions by their base
/// (there is no c.0), then the 3' UTR, and within a base the intron before it, the base, the
/// intron after it.
fn cds_pos_order_key(pos: &CdsPos) -> (bool, i32, i32) {
    (
        pos.cds_from == CdsFrom::End,
        pos.base,
        pos.offset.unwrap_or(0),
    )
}

impl<'a> AltSeqBuilder<'a> {
    /// A builder for `var_c` on `reference_data`, with the defaults of `variant::Config`.
    #[allow(deprecated)] // first_stop_pos
    pub fn new(var_c: HgvsVariant, reference_data: &'a RefTranscriptData) -> Self {
        // One scan only: `rfind` finds the usual single stop codon at the end right away.
        let first_stop_pos = reference_data.aa_sequence.find('*');
        let ref_has_multiple_stops = first_stop_pos
            .is_some_and(|first| reference_data.aa_sequence.rfind('*') != Some(first));

        Self {
            var_c,
            reference_data,
            ref_has_multiple_stops,
            first_stop_pos,
            ins_at_boundary_is_intronic: true,
            ins_at_boundary_by_sequence: false,
            boundary_intron_bases: None,
        }
    }

    /// Given a variant and a sequence, incorporate the variant and return the new sequence.
    ///
    /// Data structure returned is analogous to the data structure used to return the variant
    /// sequence, but with an additional parameter denoting the start of a frameshift that should
    /// affect all bases downstream.
    ///
    /// # Errors
    ///
    /// - `Error::ExpectedCdsVariant` if `var_c` is not a `HgvsVariant::CdsVariant`.
    /// - `Error::ValidationFailed` with `StartMustBeLessThanEnd` if the start of the
    ///   interval lies after its end, e.g. `c.6_4del`.
    /// - `Error::CoordinateOutOfBounds` for a position that does not exist, i.e. `c.0` or
    ///   `c.*0`, or that lies after the CDS, e.g. `c.5000_5001insA` on a shorter transcript.  A
    ///   position within the CDS or a `c.*N` position after the end of the transcript sequence
    ///   is clamped to its last base instead, as in biocommons/hgvs.
    pub fn build_altseq(&self) -> Result<Vec<AltTranscriptData>, Error> {
        self.build_altseq_in(self.variant_location()?, false)
    }

    /// Like `build_altseq()`, for the `region` that `variant_location()` gives, so that a
    /// caller that has it need not compute it again.  With `until_stop`, it translates less of
    /// the alt transcript, see `AltTranscriptData::new_owned()`.
    pub(crate) fn build_altseq_in(
        &self,
        region: VariantLocation,
        until_stop: bool,
    ) -> Result<Vec<AltTranscriptData>, Error> {
        // NB: the following comment is from the original Python code.
        // Should loop over each allele rather than assume only 1 variant; return a list for now.

        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return Err(Error::ExpectedCdsVariant(self.var_c.to_string()));
        };
        let loc = loc_edit.loc.inner();
        // There is no c.0 and no c.*0.
        for pos in [&loc.start, &loc.end] {
            if pos.base == 0 {
                return Err(Error::CoordinateOutOfBounds(pos.to_string()));
            }
        }
        if cds_pos_order_key(&loc.start) > cds_pos_order_key(&loc.end) {
            return Err(
                crate::validator::Error::StartMustBeLessThanEnd(self.var_c.to_string()).into(),
            );
        }
        let na_edit = loc_edit.edit.inner();
        let edit_type = match region {
            // Without a stop codon, the amino acids after the CDS are unknown.
            VariantLocation::Cds if self.after_cds_without_stop() => EditType::NotCds,
            VariantLocation::Cds => match na_edit {
                // The deleted bases are those of the interval, whether the edit gives them
                // as a sequence or as a count, as in biocommons/hgvs.
                NaEdit::RefAlt { alternative, .. } | NaEdit::NumAlt { alternative, .. } => {
                    EditType::Delins {
                        alternative: alternative.clone(),
                    }
                }
                NaEdit::DelRef { .. } | NaEdit::DelNum { .. } => EditType::Del,
                NaEdit::Ins { alternative } => EditType::Ins {
                    alternative: alternative.clone(),
                },
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
                NaEdit::RefAlt { .. } | NaEdit::NumAlt { .. } => {
                    tracing::warn!("The whole-gene variant {} is not a clean deletion. Assuming whole gene deletion.", self.var_c);
                    EditType::WholeGeneDeleted
                }
                NaEdit::Ins { .. } => {
                    tracing::warn!("Whole-gene insertion; consequence assumed to not affect protein product");
                    EditType::NotCds
                }
            },
        };

        // Get the start of the "terminal" frameshift (i.e. one never "cancelled out").
        let alt_data = match edit_type {
            EditType::Delins { alternative } => {
                self.incorporate_delins(true, Some(alternative), until_stop)
            }
            EditType::Del => self.incorporate_delins(true, None, until_stop),
            EditType::Ins { alternative } => {
                self.incorporate_delins(false, Some(alternative), until_stop)
            }
            EditType::Dup => self.incorporate_dup(until_stop),
            EditType::Inv => self.incorporate_inv(until_stop),
            EditType::NotCds => self.create_alt_equals_ref_noncds(),
            EditType::WholeGeneDeleted => self.create_no_protein(),
        }?;

        let alt_data = self.get_frameshift_start(alt_data);

        Ok(vec![alt_data])
    }

    /// The region of the transcript that `var_c` changes.
    ///
    /// `build_altseq()` and thus `c_to_p` use this region.  It depends on the fields
    /// `ins_at_boundary_is_intronic`, `ins_at_boundary_by_sequence` and `boundary_intron_bases`.
    /// So it agrees with `c_to_p` only if they hold what `c_to_p` uses: the settings of
    /// `variant::Config` and, with `ins_at_boundary_by_sequence`, the intron bases from the
    /// genome.  `AltSeqBuilder::new()` sets the defaults of `variant::Config` and no intron
    /// bases.  `assembly::Mapper::alt_seq_builder()` returns the builder that its `c_to_p()`
    /// uses.  With `variant::Config::shift_over_boundary`, `c_to_p` can build the protein from a
    /// shift of `var_c` into another region, which this method does not apply.  The region of a
    /// dup of intron bases can also depend on these bases, which `c_to_p` fills in first.
    /// `c_to_p_with_region()` of either mapper returns the region that `c_to_p` uses.
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
                let boundary_is_intronic = self.boundary_is_intronic();
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
                } else if boundary_is_intronic && edit.is_dup() && at_exon_edge {
                    // The copy may sit at the exon boundary, like an insertion there.
                    VariantLocation::Intron
                } else if !boundary_is_intronic
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

    /// Whether the bases of an insertion or dup at an exon boundary are intronic, by
    /// `ins_at_boundary_by_sequence` if it can decide, else by `ins_at_boundary_is_intronic`.
    fn boundary_is_intronic(&self) -> bool {
        self.ins_at_boundary_by_sequence
            .then(|| self.keeps_splice_site())
            .flatten()
            .unwrap_or(self.ins_at_boundary_is_intronic)
    }

    /// Whether the alt sequence of an insertion or dup at an exon boundary keeps the
    /// splice-site dinucleotide right at the exon edge, see
    /// `variant::Config::ins_at_boundary_by_sequence`.  `None` if the variant is no such
    /// insertion or dup, or if the intron bases are unknown.
    fn keeps_splice_site(&self) -> Option<bool> {
        let (edge, _, inserted, intron_dup) = self.boundary_insertion()?;
        let intron = match &self.boundary_intron_bases {
            Some(bases) if bases.len() == 2 => bases.clone(),
            // A dup of at least two intron bases next to the boundary starts or ends with them.
            _ if intron_dup && inserted.len() >= 2 => match edge {
                ExonEdge::End => inserted[..2].to_string(),
                ExonEdge::Start => inserted[inserted.len() - 2..].to_string(),
            },
            _ => return None,
        };
        Some(match edge {
            ExonEdge::End => format!("{inserted}{intron}").starts_with(&intron),
            ExonEdge::Start => format!("{intron}{inserted}").ends_with(&intron),
        })
    }

    /// The c. interval of the two intron bases next to the exon boundary at which an insertion
    /// or dup puts its bases, e.g., `c.26+1_26+2` or `c.27-2_27-1`, or `None` for other
    /// variants, see `boundary_intron_bases`.
    pub(crate) fn boundary_intron_interval(&self) -> Option<CdsInterval> {
        let (edge, base, ..) = self.boundary_insertion()?;
        let pos = |offset| CdsPos {
            base,
            offset: Some(offset),
            cds_from: CdsFrom::Start,
        };
        Some(match edge {
            ExonEdge::End => CdsInterval {
                start: pos(1),
                end: pos(2),
            },
            ExonEdge::Start => CdsInterval {
                start: pos(-2),
                end: pos(-1),
            },
        })
    }

    /// The exon edge at which an insertion or dup at an exon boundary puts its bases, the c.
    /// base of the exon next to the boundary, the bases, and whether they are a copy of the
    /// intron bases next to the boundary, or `None` for other variants.  A dup puts its copy
    /// next to the duplicated bases, so a dup that ends with the last exon base or with intron
    /// base `-1` puts it at the boundary, as does one that starts with the first exon base or
    /// with intron base `+1`.
    fn boundary_insertion(&self) -> Option<(ExonEdge, i32, String, bool)> {
        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return None;
        };
        let loc = loc_edit.loc.inner();
        if loc.start.cds_from != CdsFrom::Start || loc.end.cds_from != CdsFrom::Start {
            return None;
        }
        let ref_data = self.reference_data;
        let offsets = (loc.start.offset.unwrap_or(0), loc.end.offset.unwrap_or(0));
        match (loc_edit.edit.inner(), offsets) {
            (NaEdit::Ins { alternative }, (0, 1)) => {
                Some((ExonEdge::End, loc.start.base, alternative.clone(), false))
            }
            (NaEdit::Ins { alternative }, (-1, 0)) => {
                Some((ExonEdge::Start, loc.end.base, alternative.clone(), false))
            }
            (NaEdit::Dup { reference }, (1, _)) if !reference.is_empty() => {
                Some((ExonEdge::End, loc.start.base, reference.clone(), true))
            }
            (NaEdit::Dup { reference }, (_, -1)) if !reference.is_empty() => {
                Some((ExonEdge::Start, loc.end.base, reference.clone(), true))
            }
            (NaEdit::Dup { .. }, (0, 0)) if loc.start.base > 0 => {
                let (edge, base) = if ref_data.exon_end_positions.contains(&loc.end.base) {
                    (ExonEdge::End, loc.end.base)
                } else if ref_data.exon_start_positions.contains(&loc.start.base) {
                    (ExonEdge::Start, loc.start.base)
                } else {
                    return None;
                };
                let cds_start = usize::try_from(ref_data.cds_start).ok()?;
                let start =
                    cds_start.checked_add_signed(isize::try_from(loc.start.base).ok()?)? - 2;
                let end = cds_start.checked_add_signed(isize::try_from(loc.end.base).ok()?)? - 1;
                Some((
                    edge,
                    base,
                    ref_data.transcript_sequence.get(start..end)?.to_string(),
                    false,
                ))
            }
            _ => None,
        }
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
    /// `after_cds_without_stop()`.  A stop codon that the poly(A) tail completes ends after c.N,
    /// so a change after c.N can still change it, e.g., `c.Ndup`.
    fn outside_cds(&self) -> Option<VariantLocation> {
        let (before_cds, after_cds) = self.placements_outside_cds()?;
        if before_cds {
            Some(VariantLocation::FivePrimeUtr)
        } else if after_cds
            && self.reference_data.cds_end == CdsEnd::StopCodon
            && self.reference_data.poly_a().is_empty()
        {
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
        let at_boundary_is_exonic = !self.boundary_is_intronic();
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
        keeps_ref_suffix_and_prefix(seq, start, end, alt_bases.as_bytes(), cds_start, cds_end)
    }

    /// Helper to setup incorporate function.
    ///
    /// # Returns
    ///
    /// `(transcript sequence, cds start [1-based], cds stop [1-based], cds start index in
    /// seq [inc, 0-based], cds end index in seq [excl, 0-based])`
    ///
    /// # Errors
    ///
    /// `Error::CoordinateOutOfBounds` for `c.0` and `c.*0`, which do not exist, and for a c.
    /// position after the end of the CDS.  A position after the end of the transcript
    /// sequence that is still within the CDS or the 3' UTR is clamped to the end, as in
    /// biocommons/hgvs.
    fn setup_incorporate(&self) -> Result<(String, i32, i32, usize, usize), Error> {
        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return Err(Error::ExpectedCdsVariant(self.var_c.to_string()));
        };
        let seq = self.reference_data.transcript_sequence.to_string();
        let cds_start = usize::try_from(self.reference_data.cds_start)
            .map_err(|_| Error::CdsUndefined(self.var_c.to_string()))?;
        let cds_stop = usize::try_from(self.reference_data.cds_stop)
            .map_err(|_| Error::CdsUndefined(self.var_c.to_string()))?;

        let mut start_end = Vec::new();
        for pos in [&loc_edit.loc.inner().start, &loc_edit.loc.inner().end] {
            let out_of_bounds = || Error::CoordinateOutOfBounds(pos.to_string());
            let index = match pos.cds_from {
                CdsFrom::Start if pos.base < 0 => {
                    // 5' UTR
                    cds_start - 1
                }
                CdsFrom::Start => {
                    // There is no c.0.
                    let base = usize::try_from(pos.base)
                        .ok()
                        .filter(|&base| base > 0)
                        .ok_or_else(out_of_bounds)?;
                    let index = if pos.offset.unwrap_or(0) < 0 {
                        // Intron before `pos.base`: the exon base before it.
                        (cds_start - 1 + base).checked_sub(2)
                    } else {
                        Some(cds_start - 1 + base - 1)
                    }
                    .ok_or_else(out_of_bounds)?;
                    // A position of the CDS that is past the end of the transcript sequence
                    // is clamped to the end, as a c.*N position; one past the CDS is an
                    // error.
                    if index < seq.len() {
                        index
                    } else if index < cds_stop {
                        seq.len()
                    } else {
                        return Err(out_of_bounds());
                    }
                }
                CdsFrom::End => {
                    // There is no c.*0.
                    let base = usize::try_from(pos.base)
                        .ok()
                        .filter(|&base| base > 0)
                        .ok_or_else(out_of_bounds)?;
                    // 3' UTR, clamped to the end of the transcript.
                    std::cmp::min(cds_stop + base - 1, seq.len())
                }
            };
            start_end.push(index);
        }

        let start = start_end[0];
        let end = std::cmp::min(start_end[1] + 1, seq.len());

        Ok((
            seq,
            self.reference_data.cds_start,
            self.reference_data.cds_stop,
            start,
            end,
        ))
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

    /// Incorporate a substitution or delins (`deletes` with `alternative`), a deletion
    /// (`deletes` without `alternative`) or an insertion (`alternative` without `deletes`)
    /// into the transcript sequence.
    fn incorporate_delins(
        &self,
        deletes: bool,
        alternative: Option<String>,
        until_stop: bool,
    ) -> Result<AltTranscriptData, Error> {
        let (mut seq, cds_start, cds_stop, start, end) = self.setup_incorporate()?;
        let loc_range_start = match &self.var_c {
            HgvsVariant::CdsVariant { loc_edit, .. } => loc_edit.loc.inner().start.base,
            _ => return Err(Error::ExpectedCdsVariant(self.var_c.to_string())),
        };
        // NB: as in biocommons/hgvs, the length of the deleted reference comes from the
        // interval, not from the given reference or count.
        let ref_len = if deletes { end - start } else { 0 };
        // An insertion changes no base and goes in at `start + 1`.
        let change_end = if deletes { end } else { start + 1 };
        let alt_len = alternative.as_ref().map(|alt| alt.len()).unwrap_or(0);
        let net_base_change = alt_len as i32 - ref_len as i32;
        let mut cds_stop = cds_stop + net_base_change;
        let insert_stop = match (deletes, &alternative) {
            (false, Some(alternative)) => self
                .insert_stop_seq_end(alternative, start + 1)?
                .map(|stop_end| (stop_end, start + 1 + alternative.len())),
            _ => None,
        };

        // Incorporate the variant into the sequence (depending on the type).
        let mut is_substitution = false;
        match (deletes, alternative) {
            (true, Some(alternative)) => {
                // delins or SNP
                seq.replace_range(start..end, &alternative);
                if ref_len == 1 && alternative.len() == 1 {
                    is_substitution = true;
                }
            }
            (true, None) => {
                // deletion
                seq.replace_range(start..end, "");
            }
            (false, alternative) => {
                // insertion, after the start base
                seq.insert_str(start + 1, alternative.as_deref().unwrap_or(""));
            }
        }

        // A deletion or delins through the stop codon is a frameshift, see
        // `through_stop_codon()`.
        let mut is_frameshift = net_base_change % 3 != 0 || self.through_stop_codon();
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
            // The bases after the end of the transcript are unknown.
            self.ref_has_multiple_stops || self.ends_after_transcript(),
            self.reference_data.translation_table,
            &translation_exceptions,
            poly_a,
            until_stop,
        )
    }

    /// Whether the variant starts in the CDS and ends in the 3' UTR, or is a deletion that
    /// starts before the stop codon and ends with its last base.  HGVS describes such a
    /// deletion, which removes the stop codon, as a frameshift: "deletions starting 5' of and
    /// including the translation termination (stop) codon are described as frameshift".
    /// biocommons/hgvs does so for the first case only (biocommons/hgvs#775), see
    /// `DIFFERENCES.md`.  A delins that ends with the last base of the stop codon counts by its
    /// length, as in biocommons/hgvs, as its bases can put a stop codon in the place of the old
    /// one.
    fn through_stop_codon(&self) -> bool {
        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return false;
        };
        let loc = loc_edit.loc.inner();
        let cds_len = self.reference_data.cds_stop - self.reference_data.cds_start + 1;
        loc.start.cds_from == CdsFrom::Start
            && match loc.end.cds_from {
                CdsFrom::End => true,
                CdsFrom::Start => {
                    matches!(
                        loc_edit.edit.inner(),
                        NaEdit::DelRef { .. } | NaEdit::DelNum { .. }
                    ) && self.reference_data.cds_end == CdsEnd::StopCodon
                        && cds_len % 3 == 0
                        && loc.start.base <= cds_len - 3
                        && loc.end.base == cds_len
                        && loc.end.offset.unwrap_or(0) == 0
                }
            }
    }

    /// Whether the variant ends after the last base of the transcript sequence.  This can
    /// happen without `strict_bounds`.
    fn ends_after_transcript(&self) -> bool {
        let HgvsVariant::CdsVariant { loc_edit, .. } = &self.var_c else {
            return false;
        };
        let end = &loc_edit.loc.inner().end;
        end.cds_from == CdsFrom::End
            && i64::from(self.reference_data.cds_stop) + i64::from(end.base)
                > self.reference_data.transcript_sequence.len() as i64
    }

    /// Incorporate a dup into the transcript sequence: its copy goes right after the duplicated
    /// bases.
    fn incorporate_dup(&self, until_stop: bool) -> Result<AltTranscriptData, Error> {
        let (seq, cds_start, cds_stop, start, end) = self.setup_incorporate()?;

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
        // `cds_stop` moves with the copy, as in `incorporate_delins()`, so that a dup gives the
        // same alt protein as the insertion of its copy.  Without a stop codon, the alt protein
        // then keeps the last codons of the CDS.  biocommons/hgvs keeps `cds_stop`, see
        // `DIFFERENCES.md`.
        let mut cds_stop = cds_stop + dup_seq.len() as i32;
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
            until_stop,
        )
    }

    /// Incorporate inversion into sequence.
    fn incorporate_inv(&self, until_stop: bool) -> Result<AltTranscriptData, Error> {
        let (seq, cds_start, cds_stop, start, end) = self.setup_incorporate()?;

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
            until_stop,
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
            false,
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
            false,
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
            // adds amino acids after the stop codon at the end of `ref_seq`.  E.g., the copy of
            // `c.(N-4)_(N-2)dup` can form a new stop codon in front of the old one.
            if ref_seq == alt_seq || ref_seq.ends_with('*') && alt_seq.starts_with(ref_seq) {
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
        let first_diff = common_prefix_len(ref_seq, alt_seq);
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
            // `c_to_p` marks the result as uncertain, see
            // `variant::Config::inferred_p_is_uncertain`.
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

    use super::{AltSeqBuilder, AltTranscriptData, CdsEnd, RefTranscriptData};

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

    /// `common_prefix_len()` counts across and within its chunks of 64 bytes.
    #[test]
    fn common_prefix_len() {
        for len in [0, 1, 63, 64, 65, 128, 130] {
            let a = vec![b'A'; len];
            for other_len in [len.saturating_sub(1), len, len + 1] {
                for diff in 0..=other_len {
                    let mut b = vec![b'A'; other_len];
                    if let Some(base) = b.get_mut(diff) {
                        *base = b'C';
                    }
                    let expected = a.iter().zip(&b).take_while(|(a, b)| a == b).count();
                    assert_eq!(super::common_prefix_len(&a, &b), expected, "{len} {diff}");
                    assert_eq!(super::common_prefix_len(&b, &a), expected, "{len} {diff}");
                }
            }
        }
    }

    /// `keeps_ref_suffix_and_prefix()` agrees with building the whole alt sequence and
    /// comparing its ends with the reference, as `placements_outside_cds()` did before.  It
    /// checks every replacement of up to 8 bases, and of none, at every position of short
    /// transcripts, by every sequence of up to 3 bases and by every copy of up to 8 bases near
    /// it, as for dups and shifted insertions.  The transcripts have repeats across c.1 and
    /// c.N, no stop codon, no UTRs, and a CDS that ends beyond the transcript.
    #[test]
    fn keeps_ref_suffix_and_prefix_agrees_with_alt_sequence() {
        fn by_alt_sequence(
            seq: &[u8],
            start: usize,
            end: usize,
            alt: &[u8],
            cds_start: usize,
            cds_end: usize,
        ) -> Option<(bool, bool)> {
            if start == end && alt.is_empty() {
                return Some((start <= cds_start, start >= cds_end));
            }
            let alt = [seq.get(..start)?, alt, seq.get(end..)?].concat();
            if alt == seq {
                return None;
            }
            Some((
                alt.ends_with(seq.get(cds_start..)?),
                alt.starts_with(seq.get(..cds_end)?),
            ))
        }

        // The 5' UTR, the CDS and the 3' UTR.
        let parts = [
            (
                "AAAATCAAAC",
                "ATGAAAGCGAAAGCGTTTCGCGCGAAATAG",
                "GGGGGGCCCCC",
            ),
            ("AAAATCAAAC", "ATGAAAGCGAAAGCGTTTCGCGCGAAA", "GGGGGCCCCC"),
            ("GAAAAAA", "ATGAAAAAAAAATAA", "AAAAAAG"),
            ("CACACAC", "ATGCACACATGA", "CACACA"),
            ("", "ATGAAATAA", ""),
        ];
        // The sequence, and the 0-based indices of c.1 and behind c.N.
        let mut transcripts: Vec<_> = (parts.iter())
            .map(|(utr5, cds, utr3)| {
                let cds_start = utr5.len();
                (
                    format!("{utr5}{cds}{utr3}"),
                    cds_start,
                    cds_start + cds.len(),
                )
            })
            .collect();
        transcripts.push(("ATGAAATA".to_string(), 0, 9));

        let mut short_alts = vec![Vec::new()];
        for len in 1..=3 {
            let longer: Vec<Vec<u8>> = (short_alts.iter())
                .filter(|alt| alt.len() == len - 1)
                .flat_map(|alt| b"ACGT".map(|base| [alt.as_slice(), &[base]].concat()))
                .collect();
            short_alts.extend(longer);
        }

        let mut outcomes = std::collections::HashSet::new();
        for (seq, cds_start, cds_end) in &transcripts {
            let (seq, cds_start, cds_end) = (seq.as_bytes(), *cds_start, *cds_end);
            for start in 0..=seq.len() + 1 {
                let mut alts: Vec<&[u8]> = short_alts.iter().map(Vec::as_slice).collect();
                for len in 1..=8 {
                    for from in start.saturating_sub(8)..=start + 8 {
                        alts.extend(seq.get(from..from + len));
                    }
                }
                // An end before `start` checks reversed intervals, too.
                for end in start.saturating_sub(2)..=start + 8 {
                    for alt in &alts {
                        let expected = by_alt_sequence(seq, start, end, alt, cds_start, cds_end);
                        assert_eq!(
                            super::keeps_ref_suffix_and_prefix(
                                seq, start, end, alt, cds_start, cds_end
                            ),
                            expected,
                            "{start}..{end} by {:?} in {:?}",
                            std::str::from_utf8(alt),
                            std::str::from_utf8(seq),
                        );
                        outcomes.insert(expected);
                    }
                }
            }
        }
        for outcome in [
            None,
            Some((false, false)),
            Some((true, false)),
            Some((false, true)),
        ] {
            assert!(outcomes.contains(&outcome), "{outcome:?}");
        }
    }

    /// `ATG AAA GCG AAA GCG TTT CGC GCG AAA TAG` (`MKAKAFRAK*`) with a 5' UTR of 9 and a 3'
    /// UTR of 2 bases.
    fn mock_reference_data() -> RefTranscriptData {
        RefTranscriptData {
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

    /// `translate_until_stop()` gives a prefix of the full translation that holds the first
    /// `codons` codons and, if they have no stop codon, ends with the first one behind them,
    /// across chunks and with translation exceptions behind the first chunk.
    #[test]
    fn translate_until_stop_is_a_prefix_of_the_full_translation() -> Result<(), anyhow::Error> {
        use crate::sequences::translate_cds_with_exceptions;

        let gcc = |n: usize| "GCC".repeat(n);
        let seqs = [
            format!("ATG{}TAA{}", gcc(10), gcc(5)),
            format!("ATG{}TGA{}TAG{}", gcc(70), gcc(60), gcc(3)),
            format!("ATG{}NNNTAAGC", gcc(150)),
            format!("ATG{}GC", gcc(150)),
            String::new(),
        ];
        let exception_sets = [
            vec![],
            vec![(3, 'U')],
            vec![(71, 'U'), (100, 'M'), (500, 'W')],
        ];
        for seq in &seqs {
            for exceptions in &exception_sets {
                let full = translate_cds_with_exceptions(
                    seq,
                    false,
                    "X",
                    TranslationTable::Standard,
                    exceptions,
                )?;
                for codons in [0, 1, 5, 11, 64, 65, 71, 200, usize::MAX / 2] {
                    let part = super::translate_until_stop(
                        seq,
                        codons,
                        TranslationTable::Standard,
                        exceptions,
                    )?;
                    let head = codons.min(full.len());
                    let expected = match full[..head].find('*') {
                        Some(_) => head,
                        None => full[head..].find('*').map_or(full.len(), |i| head + i + 1),
                    };
                    assert!(full.starts_with(&part), "{seq} {exceptions:?} {codons}");
                    assert!(part.len() >= expected, "{seq} {exceptions:?} {codons}");
                    assert!(
                        part.len() == full.len() || part[..expected].contains('*'),
                        "{seq} {exceptions:?} {codons}"
                    );
                }
            }
        }
        Ok(())
    }

    /// `c_to_p` translates the alt transcript only as far as the alt protein reaches.  This
    /// gives the same alt data and `p.` as the full translation of `build_altseq()`, for
    /// substitutions, deletions, insertions, dups, delins and inversions at every CDS position
    /// and across c.N: with a 3' UTR with `N` bases, without a stop codon, with a 3' UTR without
    /// stop codons that frameshifts and extensions read into, and with an invalid base in the
    /// 3' UTR.  Only there the full translation fails, where the alt protein ends before it.
    #[test]
    fn translation_until_stop_gives_the_same_protein() -> Result<(), anyhow::Error> {
        use super::{AltSeqToHgvsp, VariantLocation};
        use crate::sequences::translate_cds;
        use std::collections::HashSet;

        let gcc = "GCC".repeat(80);
        // The 5' UTR, the CDS, the 3' UTR, and whether the 3' UTR has an invalid base.
        let transcripts = [
            (
                "AAAATCAAAC",
                "ATGAAAGCGAAAGCGTTTCGCGCGAAATAG",
                "GCCAAATAAGGNNNCCTGAGGGTAGCCCCC".to_string(),
                false,
            ),
            (
                "AAAATCAAAC",
                "ATGAAAGCGAAAGCGTTTCGCGCGAAA",
                "GGGGGCCCCC".to_string(),
                false,
            ),
            (
                "AAAATCAAAC",
                "ATGAAAGCGAAAGCGTTTCGCGCGAAATAG",
                format!("{gcc}NNN{gcc}"),
                false,
            ),
            (
                "AAAATCAAAC",
                "ATGAAAGCGAAAGCGTTTCGCGCGAAATAG",
                "GCCAAATAAGGZGGCCCCC".to_string(),
                true,
            ),
        ];
        for (utr5, cds, utr3, invalid) in transcripts {
            let transcript_sequence = format!("{utr5}{cds}{utr3}");
            let aa_sequence = translate_cds(cds, true, "*", TranslationTable::Standard)?;
            let cds_end = if aa_sequence.ends_with('*') {
                CdsEnd::StopCodon
            } else {
                CdsEnd::NoStopCodon
            };
            let reference_data = RefTranscriptData {
                transcript_sequence: transcript_sequence.as_str().into(),
                aa_sequence: aa_sequence.into(),
                cds_start: utr5.len() as i32 + 1,
                cds_stop: (utr5.len() + cds.len()) as i32,
                protein_accession: "MOCK".into(),
                translation_table: TranslationTable::Standard,
                translation_exceptions: Vec::new(),
                exon_start_positions: vec![-(utr5.len() as i32)],
                exon_end_positions: vec![(transcript_sequence.len() - utr5.len()) as i32],
                cds_end,
            };

            let n = cds.len();
            let mut hgvs_cs = vec![
                "c.-2_2del".to_string(),
                "c.-1_1insA".to_string(),
                format!("c.{n}_*1insTT"),
                format!("c.{n}_*1insA"),
                "c.*1_*2insA".to_string(),
                "c.*2del".to_string(),
            ];
            for p in (n - 3)..=n {
                hgvs_cs.push(format!("c.{p}_*1del"));
                hgvs_cs.push(format!("c.{p}_*2delinsT"));
            }
            for p in 1..=n {
                let base = &cds[p - 1..p];
                for alt in ["A", "C", "G", "T"].into_iter().filter(|alt| *alt != base) {
                    hgvs_cs.push(format!("c.{p}{base}>{alt}"));
                }
                hgvs_cs.push(format!("c.{p}del"));
                hgvs_cs.push(format!("c.{p}dup"));
                if p < n {
                    for ins in ["A", "AA", "AAA", "TAA", "TGA", "C", "GCC"] {
                        hgvs_cs.push(format!("c.{p}_{}ins{ins}", p + 1));
                    }
                    for delins in ["T", "TAA", "GG", "TGAC"] {
                        hgvs_cs.push(format!("c.{p}_{}delins{delins}", p + 1));
                    }
                    hgvs_cs.push(format!("c.{p}_{}del", p + 1));
                    hgvs_cs.push(format!("c.{p}_{}dup", p + 1));
                }
                if p + 2 <= n {
                    hgvs_cs.push(format!("c.{p}_{}del", p + 2));
                    hgvs_cs.push(format!("c.{p}_{}dup", p + 2));
                    hgvs_cs.push(format!("c.{p}_{}inv", p + 2));
                }
                if p + 3 <= n {
                    hgvs_cs.push(format!("c.{p}_{}del", p + 3));
                }
            }

            let (mut same, mut fixed, mut regions) = (0, 0, HashSet::new());
            for hgvs_c in &hgvs_cs {
                let var_c = HgvsVariant::from_str(&format!("NM_999999.1:{hgvs_c}"))?;
                let builder = AltSeqBuilder::new(var_c, &reference_data);
                let region = builder.variant_location()?;
                regions.insert(region);
                let full = builder.build_altseq_in(region, false);
                let until_stop = builder.build_altseq_in(region, true);
                match (full, until_stop) {
                    (Ok(full), Ok(until_stop)) => {
                        assert_eq!(
                            format!("{full:?}"),
                            format!("{until_stop:?}"),
                            "{hgvs_c} in {transcript_sequence}"
                        );
                        let p = |alt_data: Vec<_>| {
                            alt_data
                                .into_iter()
                                .map(|data| AltSeqToHgvsp::new(&reference_data, data).build_hgvsp())
                                .map(|var_p| format!("{var_p:?}"))
                                .collect::<Vec<_>>()
                        };
                        assert_eq!(p(full), p(until_stop), "{hgvs_c} in {transcript_sequence}");
                        same += 1;
                    }
                    (Err(_), Ok(_)) => {
                        assert!(invalid, "{hgvs_c} in {transcript_sequence}");
                        fixed += 1;
                    }
                    (full, until_stop) => assert_eq!(
                        format!("{:?}", full.map(|_| ())),
                        format!("{:?}", until_stop.map(|_| ())),
                        "{hgvs_c} in {transcript_sequence}"
                    ),
                }
            }
            assert!(same > 0, "{transcript_sequence}");
            assert_eq!(fixed > 0, invalid, "{transcript_sequence}");
            for region in [VariantLocation::Cds, VariantLocation::FivePrimeUtr] {
                assert!(
                    regions.contains(&region),
                    "{region:?} {transcript_sequence}"
                );
            }
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

    /// Positions of the CDS that are not on the transcript, and c.0, which does not exist.
    #[test]
    fn build_altseq_fails_for_positions_off_the_transcript() -> Result<(), anyhow::Error> {
        let reference_data = mock_reference_data();

        for (hgvs, position) in [
            ("MOCK:c.5000_5001insA", "5000"),
            ("MOCK:c.5000A>T", "5000"),
            ("MOCK:c.5000del", "5000"),
            ("MOCK:c.1_5000del", "5000"),
            ("MOCK:c.5000_5001dup", "5000"),
            ("MOCK:c.0A>T", "0"),
            ("MOCK:c.-1_0insA", "0"),
            ("MOCK:c.*0A>T", "*0"),
            ("MOCK:c.30_*0del", "*0"),
        ] {
            let builder = AltSeqBuilder::new(HgvsVariant::from_str(hgvs)?, &reference_data);
            let result = builder.build_altseq();
            assert!(
                matches!(&result, Err(Error::CoordinateOutOfBounds(pos)) if pos == position),
                "{hgvs} gave {result:?}"
            );
        }

        Ok(())
    }

    /// On a transcript whose sequence ends before its CDS, a CDS position past the sequence
    /// is clamped to the end, like a `c.*N` position; a position past the CDS is an error.
    /// (A variant that ends after such a CDS, e.g. `c.15_*1del`, is a non-CDS variant, see
    /// `after_cds_without_stop()`.)
    #[test]
    fn build_altseq_clamps_cds_positions_past_the_sequence() -> Result<(), anyhow::Error> {
        let reference_data = RefTranscriptData {
            transcript_sequence: "AAAATCAAAATGAAAGCGAAAGCGTTTCGC".into(),
            aa_sequence: "MKAKAFR".into(),
            cds_start: 10,
            cds_stop: 39,
            protein_accession: "MOCK".into(),
            translation_table: TranslationTable::Standard,
            translation_exceptions: Vec::new(),
            exon_start_positions: Vec::new(),
            exon_end_positions: Vec::new(),
            cds_end: CdsEnd::NoStopCodon,
        };

        for (hgvs, alt_seq) in [("MOCK:c.15_30del", "AAAATCAAAATGAAAGCGAAAGC")] {
            let builder = AltSeqBuilder::new(HgvsVariant::from_str(hgvs)?, &reference_data);
            let alt_data = builder.build_altseq()?.remove(0);
            assert_eq!(alt_data.transcript_sequence, alt_seq, "{hgvs}");
        }
        let builder =
            AltSeqBuilder::new(HgvsVariant::from_str("MOCK:c.15_31del")?, &reference_data);
        assert!(
            matches!(builder.build_altseq(), Err(Error::CoordinateOutOfBounds(pos)) if pos == "31")
        );

        Ok(())
    }

    /// Intervals whose start lies after their end, within and across the CDS.
    #[test]
    fn build_altseq_fails_for_inverted_intervals() -> Result<(), anyhow::Error> {
        let reference_data = mock_reference_data();

        for hgvs in [
            "MOCK:c.6_4del",
            "MOCK:c.6_4dup",
            "MOCK:c.6_4inv",
            "MOCK:c.6_4delinsGG",
            "MOCK:c.6_4insA",
            "MOCK:c.5_-3del",
            "MOCK:c.*2_5del",
            "MOCK:c.4+1_4del",
        ] {
            let builder = AltSeqBuilder::new(HgvsVariant::from_str(hgvs)?, &reference_data);
            let result = builder.build_altseq();
            assert!(
                matches!(
                    &result,
                    Err(Error::ValidationFailed(
                        crate::validator::Error::StartMustBeLessThanEnd(_)
                    ))
                ),
                "{hgvs} gave {result:?}"
            );
        }

        Ok(())
    }

    /// An edit with a count gives the same alt transcript as the one with the sequence, as
    /// only the length of the reference is needed.
    #[test]
    fn build_altseq_supports_counted_edits() -> Result<(), anyhow::Error> {
        let reference_data = mock_reference_data();

        for (counted, with_sequence) in [
            ("MOCK:c.4_6del3", "MOCK:c.4_6delAAA"),
            ("MOCK:c.4_6del3insGG", "MOCK:c.4_6delAAAinsGG"),
            ("MOCK:c.4del1insG", "MOCK:c.4delAinsG"),
            ("MOCK:c.4_6inv3", "MOCK:c.4_6invAAA"),
            ("MOCK:c.-3_*2del35insA", "MOCK:c.-3_*2delinsA"),
            ("MOCK:c.-3_*2del35", "MOCK:c.-3_*2del"),
        ] {
            let build = |hgvs: &str| -> Result<AltTranscriptData, anyhow::Error> {
                let builder = AltSeqBuilder::new(HgvsVariant::from_str(hgvs)?, &reference_data);
                Ok(builder.build_altseq()?.remove(0))
            };
            let (counted_data, sequence_data) = (build(counted)?, build(with_sequence)?);
            assert_eq!(
                counted_data.transcript_sequence, sequence_data.transcript_sequence,
                "{counted}"
            );
            assert_eq!(
                counted_data.aa_sequence, sequence_data.aa_sequence,
                "{counted}"
            );
            assert_eq!(counted_data.cds_stop, sequence_data.cds_stop, "{counted}");
            assert_eq!(
                counted_data.is_frameshift, sequence_data.is_frameshift,
                "{counted}"
            );
            assert_eq!(
                counted_data.variant_start_aa, sequence_data.variant_start_aa,
                "{counted}"
            );
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
