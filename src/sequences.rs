//! Utility code for working with sequences.
//!
//! Partially ported over from `bioutils.sequences`.

pub use crate::sequences::error::Error;
use crate::Sequence;
use ahash::AHashMap;
use md5::{Digest, Md5};
use std::sync::LazyLock;

include!(concat!(env!("OUT_DIR"), "/tables_gen.rs"));

mod error {
    /// Error type for normalization of HGVS expressins.
    #[derive(thiserror::Error, Debug, Clone)]
    pub enum Error {
        #[error("invalid 1-letter aminoacid: {0} at {1}")]
        InvalidOneLetterAminoAcid(String, String),
        #[error("invalid 3-letter aminoacid: {0} at {1}")]
        InvalidThreeLetterAminoAcid(String, String),
        #[error("3-letter amino acid sequence length is not multiple of three: {0}")]
        InvalidThreeLetterAminoAcidLength(usize),
        #[error("codon is undefined in codon table: {0}")]
        UndefinedCodon(String),
        #[error("can only translate DNA sequences whose length is multiple of 3, but is: {0}")]
        UntranslatableDnaLenth(usize),
        #[error("character is not alphabetic: {0}")]
        NotAlphabetic(char),
    }
}

pub fn trim_common_prefixes(reference: &[u8], alternative: &[u8]) -> (usize, Sequence, Sequence) {
    if reference.is_empty() || alternative.is_empty() {
        return (0, reference.to_vec(), alternative.to_vec());
    }

    let mut trim = 0;
    while trim < reference.len() && trim < alternative.len() {
        if reference.get(trim) != alternative.get(trim) {
            break;
        }

        trim += 1;
    }

    (trim, reference[trim..].into(), alternative[trim..].into())
}

pub fn trim_common_suffixes(reference: &[u8], alternative: &[u8]) -> (usize, Sequence, Sequence) {
    if reference.is_empty() || alternative.is_empty() {
        return (0, reference.into(), alternative.into());
    }

    let mut trim = 0;
    let mut i_r = reference.len();
    let mut i_a = alternative.len();
    let mut pad = 0;
    while trim < reference.len() && trim < alternative.len() {
        trim += 1;
        assert!(i_r > 0);
        assert!(i_a > 0);
        i_r -= 1;
        i_a -= 1;

        if reference.get(i_r) != alternative.get(i_a) {
            pad = 1;
            break;
        }
    }

    (
        trim - pad,
        reference[..(i_r + pad)].into(),
        alternative[..(i_a + pad)].into(),
    )
}

/// Reverse complementing shortcut.
pub fn revcomp(seq: &[u8]) -> Sequence {
    bio::alphabets::dna::revcomp(seq)
}

/// Allow selection of translation table.
#[derive(
    Debug,
    Default,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    serde::Serialize,
    serde::Deserialize,
)]
pub enum TranslationTable {
    #[default]
    Standard,
    Selenocysteine,
    VertebrateMitochondrial,
}

/// Coerces string of 1- or 3-letter amino acids to 1-letter representation.
///
/// Fails if the sequence is not of valid 3/1-letter amino acids.
///
/// # Args
///
/// * `seq` -- An amino acid sequence.
///
/// # Returns
///
/// The sequence as one of 1-letter amino acids.
#[allow(dead_code)]
pub fn aa_to_aa1(seq: &[u8]) -> Result<String, Error> {
    if looks_like_aa3_p(seq) {
        aa3_to_aa1(seq)
    } else {
        Ok(String::from_utf8_lossy(seq).to_string())
    }
}

/// Coerces string of 1- or 3-letter amino acids to 3-letter representation.
///
/// Fails if the sequence is not of valid 3/1-letter amino acids.
///
/// # Args
///
/// * `seq` -- An amino acid sequence.
///
/// # Returns
///
/// The sequence as one of 1-letter amino acids.
#[allow(dead_code)]
pub fn aa_to_aa3(seq: &[u8]) -> Result<String, Error> {
    if looks_like_aa3_p(seq) {
        Ok(String::from_utf8_lossy(seq).to_string())
    } else {
        aa1_to_aa3(seq)
    }
}

/// Converts string of 1-letter amino acids to 3-letter amino acids.
///
/// Fails if the sequence is not of 1-letter amino acids.
///
/// # Args
///
/// * `seq` -- An amino acid sequence as 1-letter amino acids.
///
/// # Returns
///
/// The sequence as 3-letter amino acids.
#[allow(dead_code)]
pub fn aa1_to_aa3(seq: &[u8]) -> Result<String, Error> {
    if seq.is_empty() {
        return Ok(String::new());
    }

    let mut result = String::with_capacity(seq.len() * 3);

    for (i, aa1) in seq.iter().enumerate() {
        let aa3 = AA1_TO_AA3_STR[*aa1 as usize].ok_or_else(|| {
            Error::InvalidOneLetterAminoAcid(format!("{:?}", aa1), format!("{}", i + 1))
        })?;
        result.push_str(aa3);
    }

    Ok(result)
}

/// Converts string of 3-letter amino acids to 1-letter amino acids.
///
/// Fails if the sequence is not of 3-letter amino acids.
///
/// # Args
///
/// * `seq` -- An amino acid sequence as 3-letter amino acids.
///
/// # Returns
///
/// The sequence as 1-letter amino acids.
#[allow(dead_code)]
pub fn aa3_to_aa1(seq: &[u8]) -> Result<String, Error> {
    if seq.len() % 3 != 0 {
        return Err(Error::InvalidThreeLetterAminoAcidLength(seq.len()));
    }

    let mut result = String::with_capacity(seq.len() / 3);

    for (i, aa3) in seq.chunks(3).enumerate() {
        let aa1 = _aa3_to_aa1(aa3).ok_or_else(|| {
            Error::InvalidThreeLetterAminoAcid(format!("{:?}", aa3), format!("{}", i + 1))
        })? as char;
        result.push(aa1);
    }

    Ok(result)
}

/// Indicates whether a string looks like a 3-letter AA string.
///
/// # Args
///
/// * `looks_like_aa3_p` -- A sequence
///
/// # Returns
///
/// Whether the string is of the format of a 3-letter AA string.
#[allow(dead_code)]
fn looks_like_aa3_p(seq: &[u8]) -> bool {
    seq.len() % 3 == 0 && seq.get(1).map(|c| c.is_ascii_lowercase()).unwrap_or(true)
}

type Codon = [u8; 3];

/// Allow translation of `&[u8]` DNA codons to `u8` amino acids.
///
/// We use separate structs here to encapsulate getting the lazy static global data.
struct CodonTranslator {
    /// Mapping for "normalizing" DNA ASCII character (to upper case and `U -> T`).
    dna_ascii_map: [u8; 256],

    /// Mapping from DNA ASCII to 2-bit representation.
    dna_ascii_to_2bit: [u8; 256],

    /// IUPAC ambiguity codes.
    iupac_ambiguity_codes: [u8; 13],

    /// Mapping from 2bit DNA codon to amino acid 1-letter ASCII.
    codon_2bit_to_aa1: [u8; 64],

    /// Mapping from DNA 2-bit to amino acid 1-letter ASCII including degenerate codons.
    full_dna_to_aa1: &'static AHashMap<Codon, u8>,

    /// Buffer.
    codon: Codon,
}

static DNA_TO_AA1_LUT: LazyLock<AHashMap<Codon, u8>> = LazyLock::new(|| {
    let mut m = AHashMap::default();
    for (dna, aa1) in DNA_TO_AA1_LUT_VEC {
        assert_eq!(dna.len(), 3);
        let d = dna.as_bytes();
        m.insert([d[0], d[1], d[2]], aa1.as_bytes()[0]);
    }
    m
});

static DNA_TO_AA1_SEC: LazyLock<AHashMap<Codon, u8>> = LazyLock::new(|| {
    let mut m = AHashMap::default();
    for (dna, aa1) in DNA_TO_AA1_SEC_VEC {
        assert_eq!(dna.len(), 3);
        let d = dna.as_bytes();
        m.insert([d[0], d[1], d[2]], aa1.as_bytes()[0]);
    }
    m
});

static DNA_TO_AA1_CHRMT_VERTEBRATE: LazyLock<AHashMap<Codon, u8>> = LazyLock::new(|| {
    let mut m = AHashMap::default();
    for (dna, aa1) in DNA_TO_AA1_CHRMT_VERTEBRATE_VEC {
        assert_eq!(dna.len(), 3);
        let d = dna.as_bytes();
        m.insert([d[0], d[1], d[2]], aa1.as_bytes()[0]);
    }
    m
});

impl CodonTranslator {
    /// Initialize the struct.
    pub fn new(table: TranslationTable) -> Self {
        Self {
            dna_ascii_map: DNA_ASCII_MAP,
            dna_ascii_to_2bit: DNA_ASCII_TO_2BIT,
            iupac_ambiguity_codes: IUPAC_AMBIGUITY_CODES,

            codon_2bit_to_aa1: match table {
                TranslationTable::Standard => CODON_2BIT_TO_AA1_LUT,
                TranslationTable::Selenocysteine => CODON_2BIT_TO_AA1_SEC,
                TranslationTable::VertebrateMitochondrial => CODON_2BIT_TO_AA1_CHRMT_VERTEBRATE,
            },
            full_dna_to_aa1: match table {
                TranslationTable::Standard => &DNA_TO_AA1_LUT,
                TranslationTable::Selenocysteine => &DNA_TO_AA1_SEC,
                TranslationTable::VertebrateMitochondrial => &DNA_TO_AA1_CHRMT_VERTEBRATE,
            },

            codon: [0; 3],
        }
    }

    /// Translate the given codon to an amino acid.
    ///
    /// # Args
    ///
    /// * `codon` -- A codon.
    ///
    /// # Returns
    ///
    /// The corresponding amino acid.
    pub fn translate(&mut self, codon: &[u8]) -> Result<u8, Error> {
        // Normalize (to upper case etc.) codon.
        self.normalize_codon(codon);

        let translation = self
            // Attempt fast translation of codon
            .codon_to_aa1(&self.codon)
            // Fast translation fails, but slower hash map succeeded.
            .or_else(|| self.full_dna_to_aa1.get(&self.codon).copied())
            // If this contains an ambiguous code, set aa to X, otherwise, throw error
            .or_else(|| {
                codon
                    .iter()
                    .any(|c| self.iupac_ambiguity_codes.contains(c))
                    .then_some(b'X')
            });
        translation.ok_or_else(|| {
            Error::UndefinedCodon(
                std::str::from_utf8(codon)
                    .expect("cannot decode UTF-8")
                    .to_owned(),
            )
        })
    }

    fn dna3_to_2bit(&self, c: &[u8]) -> Option<u8> {
        let mut result = 0;
        for i in &c[..3] {
            result <<= 2;
            let tmp = self.dna_ascii_to_2bit[*i as usize];
            if tmp == 255 {
                return None;
            }
            result |= tmp;
        }
        Some(result)
    }

    /// Helper function to extract normalized codon to `self.codon`.
    fn normalize_codon(&mut self, codon: &[u8]) {
        for (i, c) in codon[..3].iter().enumerate() {
            self.codon[i] = self.dna_ascii_map[*c as usize];
        }
    }

    fn codon_to_aa1(&self, codon: &[u8]) -> Option<u8> {
        if let Some(val) = self.dna3_to_2bit(codon) {
            let tmp = self.codon_2bit_to_aa1[val as usize];
            if tmp == 0 {
                None
            } else {
                Some(tmp)
            }
        } else {
            DNA_TO_AA1_LUT.get(codon).copied()
        }
    }
}

/// Translates a DNA or RNA sequence into a single-letter amino acid sequence.
///
/// # Args
///
/// * `seq` -- A nucleotide sequence.
/// * `full_codons` -- If `true`, forces sequence to have length that is a multiple of 3
///    and return an `Err` otherwise.  If `false`, `ter_symbol` will be added as the last
///    amino acid.  This corresponds to biopython's behavior of padding the last codon with
///    `N` characters.
/// * `ter_symbol` -- Placeholder for the last amino acid if sequence length is not divisible
///    by three and `full_codons` is `false`.
/// * `translation_table` -- Indicates which codon to amino acid translation table to use.
///
/// # Returns
///
/// The corresponding single letter amino acid sequence.
pub fn translate_cds(
    seq: &[u8],
    full_codons: bool,
    ter_symbol: u8,
    translation_table: TranslationTable,
) -> Result<Sequence, Error> {
    if seq.is_empty() {
        return Ok(b"".into());
    }

    if full_codons && seq.len() % 3 != 0 {
        return Err(Error::UntranslatableDnaLenth(seq.len()));
    }

    // Translate the codons from the input to result.
    let mut translator = CodonTranslator::new(translation_table);
    let mut result = Vec::with_capacity(seq.len() / 3);
    for chunk in seq.chunks_exact(3) {
        result.push(translator.translate(chunk)?);
    }

    // Check for trailing bases and add the ter symbol if required.
    if !full_codons && seq.len() % 3 != 0 {
        result.push(ter_symbol);
    }

    Ok(result)
}

/// Converts sequence to normalized representation for hashing.
///
/// Essentially, removes whitespace and asterisks, and uppercases the string.
///
/// # Args
///
/// * `seq` -- The sequence to be normalized.
///
/// # Returns
///
/// The sequence as a string of uppercase letters.
pub fn normalize_sequence(seq: &[u8]) -> Result<Sequence, Error> {
    let mut result = Sequence::new();

    for c in seq {
        if !c.is_ascii_whitespace() && *c != b'*' {
            let c = c.to_ascii_uppercase();
            if c.is_ascii_alphabetic() {
                result.push(c)
            } else {
                return Err(Error::NotAlphabetic(c as char));
            }
        }
    }

    Ok(result)
}

/// Convert sequence to unicode MD5 hex digest.
///
/// Fails if normalization is not possible.
///
/// # Args
///
/// * `seq` -- A sequence
/// * `normalize` -- Whether to normalize the sequence before conversion, i.e., to ensure
///   representation as uppercase ltters without whitespace or asterisks.
///
/// # Returns
///
/// Unicode MD5 hex digest representation of sequence.
pub fn seq_md5(seq: &[u8], normalize: bool) -> Result<String, Error> {
    let seq = if normalize {
        normalize_sequence(seq)?
    } else {
        seq.to_owned()
    };
    let mut hasher = Md5::new();
    hasher.update(seq);
    let hash = hasher.finalize();
    let mut buf = [0u8; 64];
    let checksum =
        base16ct::lower::encode_str(&hash, &mut buf).expect("cannot perform base16 encoding");
    Ok(checksum.to_owned())
}

#[cfg(test)]
mod test {
    use super::*;

    use pretty_assertions::assert_eq;

    #[test]
    fn suffix_trimming() {
        assert_eq!(
            trim_common_suffixes(b"", b""),
            (0, vec![], vec![])
        );
        assert_eq!(
            trim_common_suffixes(b"", b"C"),
            (0, vec![], b"C".to_vec())
        );
        assert_eq!(
            trim_common_suffixes(b"C", b""),
            (0, b"C".to_vec(), vec![])
        );
        assert_eq!(
            trim_common_suffixes(b"A", b"AA"),
            (1, vec![], b"A".to_vec())
        );
        assert_eq!(
            trim_common_suffixes(b"AT", b"AG"),
            (0, b"AT".to_vec(), b"AG".to_vec())
        );
        assert_eq!(
            trim_common_suffixes(b"ATCG", b"AGCG"),
            (2, b"AT".to_vec(), b"AG".to_vec())
        );
    }

    #[test]
    fn prefix_trimming() {
        assert_eq!(
            trim_common_prefixes(b"", b""),
            (0, vec![], vec![])
        );
        assert_eq!(
            trim_common_prefixes(b"", b"C"),
            (0, vec![], b"C".to_vec())
        );
        assert_eq!(
            trim_common_prefixes(b"C", b""),
            (0, b"C".to_vec(), vec![])
        );
        assert_eq!(
            trim_common_prefixes(b"TA", b"GA"),
            (0, b"TA".to_vec(), b"GA".to_vec())
        );
        assert_eq!(
            trim_common_prefixes(b"CGTA", b"CGGA"),
            (2, b"TA".to_vec(), b"GA".to_vec())
        );
    }

    #[test]
    fn revcomp_cases() {
        assert_eq!(revcomp(b""), b"");
        assert_eq!(revcomp(b"A"), b"T");
        assert_eq!(revcomp(b"AG"), b"CT");
        assert_eq!(revcomp(b"CGAG"), b"CTCG");
    }

    #[test]
    fn aa_to_aa1_examples() -> Result<(), Error> {
        assert_eq!(aa_to_aa1(b"")?, "");
        assert_eq!(aa_to_aa1(b"CATSARELAME")?, "CATSARELAME");
        assert_eq!(
            aa_to_aa1(b"CysAlaThrSerAlaArgGluLeuAlaMetGlu")?,
            "CATSARELAME"
        );

        Ok(())
    }

    #[test]
    fn aa1_to_aa3_examples() -> Result<(), Error> {
        assert_eq!(aa1_to_aa3(b"")?, "");
        assert_eq!(
            aa1_to_aa3(b"CATSARELAME")?,
            "CysAlaThrSerAlaArgGluLeuAlaMetGlu"
        );

        Ok(())
    }

    #[test]
    fn aa3_to_aa1_examples() -> Result<(), Error> {
        assert!(aa3_to_aa1(b"Te").is_err());
        assert_eq!(aa3_to_aa1(b"")?, "");
        assert_eq!(
            aa3_to_aa1(b"CysAlaThrSerAlaArgGluLeuAlaMetGlu")?,
            "CATSARELAME"
        );

        Ok(())
    }

    #[test]
    fn translate_cds_examples() -> Result<(), Error> {
        assert_eq!(
            translate_cds(b"ATGCGA", true, b'*', TranslationTable::Standard)?,
            b"MR".to_vec()
        );
        assert_eq!(
            translate_cds(b"AUGCGA", true, b'*', TranslationTable::Standard)?,
            b"MR".to_vec()
        );
        assert_eq!(
            translate_cds(b"", true, b'*', TranslationTable::Standard)?,
            b"".to_vec()
        );
        assert!(translate_cds(b"AUGCG", true, b'*', TranslationTable::Standard).is_err());
        assert_eq!(
            translate_cds(b"AUGCG", false, b'*', TranslationTable::Standard)?,
            b"M*".to_vec()
        );
        assert_eq!(
            translate_cds(b"ATGTAN", true, b'*', TranslationTable::Standard)?,
            b"MX".to_vec()
        );
        assert_eq!(
            translate_cds(b"CCN", true, b'*', TranslationTable::Standard)?,
            b"P".to_vec()
        );
        assert_eq!(
            translate_cds(b"TRA", true, b'*', TranslationTable::Standard)?,
            b"*".to_vec()
        );
        assert_eq!(
            translate_cds(b"TTNTA", false, b'*', TranslationTable::Standard)?,
            b"X*".to_vec()
        );
        assert_eq!(
            translate_cds(b"CTB", true, b'*', TranslationTable::Standard)?,
            b"L".to_vec()
        );
        assert_eq!(
            translate_cds(b"AGM", true, b'*', TranslationTable::Standard)?,
            b"X".to_vec()
        );
        assert_eq!(
            translate_cds(b"GAS", true, b'*', TranslationTable::Standard)?,
            b"X".to_vec()
        );
        assert_eq!(
            translate_cds(b"CUN", true, b'*', TranslationTable::Standard)?,
            b"L".to_vec()
        );
        assert!(translate_cds(b"AUGCGQ", true, b'*', TranslationTable::Standard).is_err());

        Ok(())
    }

    #[test]
    fn seq_md5_examples() -> Result<(), Error> {
        assert_eq!(seq_md5(b"", true)?, "d41d8cd98f00b204e9800998ecf8427e");
        assert_eq!(seq_md5(b"ACGT", true)?, "f1f8f4bf413b16ad135722aa4591043e");
        assert_eq!(seq_md5(b"ACGT*", true)?, "f1f8f4bf413b16ad135722aa4591043e");
        assert_eq!(
            seq_md5(b" A C G T ", true)?,
            "f1f8f4bf413b16ad135722aa4591043e"
        );
        assert_eq!(seq_md5(b"acgt", true)?, "f1f8f4bf413b16ad135722aa4591043e");
        assert_eq!(seq_md5(b"acgt", false)?, "db516c3913e179338b162b2476d1c23f");

        Ok(())
    }

    #[test]
    fn normalize_sequence_examples() -> Result<(), Error> {
        assert_eq!(normalize_sequence(b"ACGT")?, b"ACGT".to_vec());
        assert_eq!(normalize_sequence(b"  A C G T * ")?, b"ACGT".to_vec());
        assert!(normalize_sequence(b"ACGT1").is_err());

        Ok(())
    }

    #[test]
    fn exercise_lazy_ds() {
        assert_eq!(DNA_ASCII_MAP[0], b'\0');
        assert_eq!(DNA_ASCII_TO_2BIT[b'A' as usize], 0);
        assert_eq!(AA3_TO_AA1_VEC[0], ("Ala", "A"));
        assert_eq!(DNA_TO_AA1_LUT_VEC[0], ("AAA", "K"));
        assert_eq!(DNA_TO_AA1_SEC_VEC[0], ("AAA", "K"));
        assert_eq!(DNA_TO_AA1_CHRMT_VERTEBRATE_VEC[0], ("AAA", "K"));
    }

    #[test]
    fn codon_translator_standard() -> Result<(), Error> {
        let mut translator = CodonTranslator::new(TranslationTable::Standard);

        // Non-denenerate codon.
        assert_eq!(translator.translate(b"AAA")?, b'K');
        // Degenerate codon.
        assert_eq!(translator.translate(b"AAR")?, b'K');

        Ok(())
    }

    #[test]
    fn codon_translator_sec() -> Result<(), Error> {
        let mut translator = CodonTranslator::new(TranslationTable::Selenocysteine);

        // Non-denenerate codon.
        assert_eq!(translator.translate(b"AAA")?, b'K');
        // Degenerate codon.
        assert_eq!(translator.translate(b"AAR")?, b'K');

        Ok(())
    }

    #[test]
    fn codon_translator_chrmt_vertebrate() -> Result<(), Error> {
        let mut translator = CodonTranslator::new(TranslationTable::Selenocysteine);

        // Non-denenerate codon.
        assert_eq!(translator.translate(b"AAA")?, b'K');
        // Degenerate codon.
        assert_eq!(translator.translate(b"AAR")?, b'K');

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
