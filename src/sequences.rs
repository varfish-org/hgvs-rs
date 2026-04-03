//! Utility code for working with sequences.
//!
//! Partially ported over from `bioutils.sequences`.

use md5::{Digest, Md5};

pub use crate::sequences::error::Error;

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

pub fn trim_common_prefixes(reference: &str, alternative: &str) -> (usize, String, String) {
    let (trim, r_slice, a_slice) = trim_common_prefixes_slice(reference, alternative);
    (trim, r_slice.to_string(), a_slice.to_string())
}

pub fn trim_common_suffixes(reference: &str, alternative: &str) -> (usize, String, String) {
    let (trim, r_slice, a_slice) = trim_common_suffixes_slice(reference, alternative);
    (trim, r_slice.to_string(), a_slice.to_string())
}

pub fn trim_common_prefixes_slice<'a>(
    reference: &'a str,
    alternative: &'a str,
) -> (usize, &'a str, &'a str) {
    let trim = reference
        .as_bytes()
        .iter()
        .zip(alternative.as_bytes().iter())
        .take_while(|(r, a)| r == a)
        .count();

    (trim, &reference[trim..], &alternative[trim..])
}

pub fn trim_common_suffixes_slice<'a>(
    reference: &'a str,
    alternative: &'a str,
) -> (usize, &'a str, &'a str) {
    let trim = reference
        .as_bytes()
        .iter()
        .rev()
        .zip(alternative.as_bytes().iter().rev())
        .take_while(|(r, a)| r == a)
        .count();

    (
        trim,
        &reference[..reference.len() - trim],
        &alternative[..alternative.len() - trim],
    )
}

/// Reverse complementing shortcut.
pub fn revcomp(seq: &str) -> String {
    std::str::from_utf8(&bio::alphabets::dna::revcomp(seq.as_bytes()))
        .expect("invalid utf-8 encoding")
        .to_string()
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
pub fn aa_to_aa1(seq: &str) -> Result<String, Error> {
    if looks_like_aa3_p(seq) {
        aa3_to_aa1(seq)
    } else {
        Ok(seq.to_string())
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
pub fn aa_to_aa3(seq: &str) -> Result<String, Error> {
    if looks_like_aa3_p(seq) {
        Ok(seq.to_string())
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
pub fn aa1_to_aa3(seq: &str) -> Result<String, Error> {
    if seq.is_empty() {
        return Ok(String::new());
    }

    let mut result = String::with_capacity(seq.len() * 3);

    for (i, aa1) in seq.as_bytes().iter().enumerate() {
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
pub fn aa3_to_aa1(seq: &str) -> Result<String, Error> {
    if seq.len() % 3 != 0 {
        return Err(Error::InvalidThreeLetterAminoAcidLength(seq.len()));
    }

    let mut result = String::with_capacity(seq.len() / 3);

    for (i, aa3) in seq.as_bytes().chunks(3).enumerate() {
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
fn looks_like_aa3_p(seq: &str) -> bool {
    seq.len() % 3 == 0 && seq.chars().nth(1).map(|c| c.is_lowercase()).unwrap_or(true)
}

/// Translates a DNA or RNA sequence into a single-letter amino acid sequence.
///
/// # Args
///
/// * `seq` -- A nucleotide sequence.
/// * `full_codons` -- If `true`, forces sequence to have length that is a multiple of 3
///   and return an `Err` otherwise.  If `false`, `ter_symbol` will be added as the last
///   amino acid.  This corresponds to biopython's behavior of padding the last codon with
///   `N` characters.
/// * `ter_symbol` -- Placeholder for the last amino acid if sequence length is not divisible
///   by three and `full_codons` is `false`.
/// * `translation_table` -- Indicates which codon to amino acid translation table to use.
///
/// # Returns
///
/// The corresponding single letter amino acid sequence.
pub fn translate_cds(
    seq: &str,
    full_codons: bool,
    ter_symbol: &str,
    translation_table: TranslationTable,
) -> Result<String, Error> {
    if seq.is_empty() {
        return Ok(String::new());
    }

    if full_codons && seq.len() % 3 != 0 {
        return Err(Error::UntranslatableDnaLenth(seq.len()));
    }

    let lut = match translation_table {
        TranslationTable::Standard => &CODON_IUPAC_TO_AA1_LUT,
        TranslationTable::Selenocysteine => &CODON_IUPAC_TO_AA1_SEC,
        TranslationTable::VertebrateMitochondrial => &CODON_IUPAC_TO_AA1_CHRMT_VERTEBRATE,
    };

    let mut result = Vec::with_capacity(seq.len() / 3 + ter_symbol.len());

    for chunk in seq.as_bytes().chunks_exact(3) {
        // map raw ascii to 4-bit IUPAC (handles U to T and lower to upper implicitly)
        let v0 = DNA_ASCII_TO_IUPAC[chunk[0] as usize];
        let v1 = DNA_ASCII_TO_IUPAC[chunk[1] as usize];
        let v2 = DNA_ASCII_TO_IUPAC[chunk[2] as usize];

        // 255 means invalid character
        if v0 != 255 && v1 != 255 && v2 != 255 {
            // pack into 12-bit integer
            let idx = ((v0 as usize) << 8) | ((v1 as usize) << 4) | (v2 as usize);
            result.push(lut[idx]);
        } else {
            return Err(Error::UndefinedCodon(
                std::str::from_utf8(chunk).unwrap().to_owned(),
            ));
        }
    }

    if !full_codons && seq.len() % 3 != 0 {
        result.extend_from_slice(ter_symbol.as_bytes());
    }

    // SAFETY: The translation tables and ter_symbol only return valid ASCII bytes.
    Ok(unsafe { String::from_utf8_unchecked(result) })
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
pub fn normalize_sequence(seq: &str) -> Result<String, Error> {
    let mut result = String::new();

    for c in seq.chars() {
        if !c.is_whitespace() && c != '*' {
            let c = c.to_ascii_uppercase();
            if c.is_alphabetic() {
                result.push(c)
            } else {
                return Err(Error::NotAlphabetic(c));
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
pub fn seq_md5(seq: &str, normalize: bool) -> Result<String, Error> {
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
            trim_common_suffixes("", ""),
            (0, "".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_suffixes("", "C"),
            (0, "".to_string(), "C".to_string())
        );
        assert_eq!(
            trim_common_suffixes("C", ""),
            (0, "C".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_suffixes("A", "AA"),
            (1, "".to_string(), "A".to_string())
        );
        assert_eq!(
            trim_common_suffixes("AT", "AG"),
            (0, "AT".to_string(), "AG".to_string())
        );
        assert_eq!(
            trim_common_suffixes("ATCG", "AGCG"),
            (2, "AT".to_string(), "AG".to_string())
        );
    }

    #[test]
    fn prefix_trimming() {
        assert_eq!(
            trim_common_prefixes("", ""),
            (0, "".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_prefixes("", "C"),
            (0, "".to_string(), "C".to_string())
        );
        assert_eq!(
            trim_common_prefixes("C", ""),
            (0, "C".to_string(), "".to_string())
        );
        assert_eq!(
            trim_common_prefixes("TA", "GA"),
            (0, "TA".to_string(), "GA".to_string())
        );
        assert_eq!(
            trim_common_prefixes("CGTA", "CGGA"),
            (2, "TA".to_string(), "GA".to_string())
        );
    }

    #[test]
    fn revcomp_cases() {
        assert_eq!(revcomp(""), "");
        assert_eq!(revcomp("A"), "T");
        assert_eq!(revcomp("AG"), "CT");
        assert_eq!(revcomp("CGAG"), "CTCG");
    }

    #[test]
    fn aa_to_aa1_examples() -> Result<(), Error> {
        assert_eq!(aa_to_aa1("")?, "");
        assert_eq!(aa_to_aa1("CATSARELAME")?, "CATSARELAME");
        assert_eq!(
            aa_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")?,
            "CATSARELAME"
        );

        Ok(())
    }

    #[test]
    fn aa1_to_aa3_examples() -> Result<(), Error> {
        assert_eq!(aa1_to_aa3("")?, "");
        assert_eq!(
            aa1_to_aa3("CATSARELAME")?,
            "CysAlaThrSerAlaArgGluLeuAlaMetGlu"
        );

        Ok(())
    }

    #[test]
    fn aa3_to_aa1_examples() -> Result<(), Error> {
        assert!(aa3_to_aa1("Te").is_err());
        assert_eq!(aa3_to_aa1("")?, "");
        assert_eq!(
            aa3_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")?,
            "CATSARELAME"
        );

        Ok(())
    }

    #[test]
    fn translate_cds_examples() -> Result<(), Error> {
        assert_eq!(
            translate_cds("ATGCGA", true, "*", TranslationTable::Standard)?,
            "MR"
        );
        assert_eq!(
            translate_cds("AUGCGA", true, "*", TranslationTable::Standard)?,
            "MR"
        );
        assert_eq!(
            translate_cds("", true, "*", TranslationTable::Standard)?,
            ""
        );
        assert!(translate_cds("AUGCG", true, "*", TranslationTable::Standard).is_err());
        assert_eq!(
            translate_cds("AUGCG", false, "*", TranslationTable::Standard)?,
            "M*"
        );
        assert_eq!(
            translate_cds("ATGTAN", true, "*", TranslationTable::Standard)?,
            "MX"
        );
        assert_eq!(
            translate_cds("CCN", true, "*", TranslationTable::Standard)?,
            "P"
        );
        assert_eq!(
            translate_cds("TRA", true, "*", TranslationTable::Standard)?,
            "*"
        );
        assert_eq!(
            translate_cds("TTNTA", false, "*", TranslationTable::Standard)?,
            "X*"
        );
        assert_eq!(
            translate_cds("CTB", true, "*", TranslationTable::Standard)?,
            "L"
        );
        assert_eq!(
            translate_cds("AGM", true, "*", TranslationTable::Standard)?,
            "X"
        );
        assert_eq!(
            translate_cds("GAS", true, "*", TranslationTable::Standard)?,
            "X"
        );
        assert_eq!(
            translate_cds("CUN", true, "*", TranslationTable::Standard)?,
            "L"
        );
        assert!(translate_cds("AUGCGQ", true, "*", TranslationTable::Standard).is_err());

        Ok(())
    }

    #[test]
    fn seq_md5_examples() -> Result<(), Error> {
        assert_eq!(seq_md5("", true)?, "d41d8cd98f00b204e9800998ecf8427e");
        assert_eq!(seq_md5("ACGT", true)?, "f1f8f4bf413b16ad135722aa4591043e");
        assert_eq!(seq_md5("ACGT*", true)?, "f1f8f4bf413b16ad135722aa4591043e");
        assert_eq!(
            seq_md5(" A C G T ", true)?,
            "f1f8f4bf413b16ad135722aa4591043e"
        );
        assert_eq!(seq_md5("acgt", true)?, "f1f8f4bf413b16ad135722aa4591043e");
        assert_eq!(seq_md5("acgt", false)?, "db516c3913e179338b162b2476d1c23f");

        Ok(())
    }

    #[test]
    fn normalize_sequence_examples() -> Result<(), Error> {
        assert_eq!(normalize_sequence("ACGT")?, "ACGT");
        assert_eq!(normalize_sequence("  A C G T * ")?, "ACGT");
        assert!(normalize_sequence("ACGT1").is_err());

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
