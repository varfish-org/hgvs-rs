//! Utility code for working with sequences.
//!
//! Partially ported over from `bioutils.sequences`.

use ahash::AHashMap;
use md5::{Digest, Md5};

pub use crate::sequences::error::Error;

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
    if reference.is_empty() || alternative.is_empty() {
        return (0, reference.to_string(), alternative.to_string());
    }

    let mut trim = 0;
    while trim < reference.len() && trim < alternative.len() {
        if reference.chars().nth(trim) != alternative.chars().nth(trim) {
            break;
        }

        trim += 1;
    }

    (
        trim,
        reference[trim..].to_string(),
        alternative[trim..].to_string(),
    )
}

pub fn trim_common_suffixes(reference: &str, alternative: &str) -> (usize, String, String) {
    if reference.is_empty() || alternative.is_empty() {
        return (0, reference.to_string(), alternative.to_string());
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

        if reference.chars().nth(i_r) != alternative.chars().nth(i_a) {
            pad = 1;
            break;
        }
    }

    (
        trim - pad,
        reference[..(i_r + pad)].to_string(),
        alternative[..(i_a + pad)].to_string(),
    )
}

/// Reverse complementing shortcut.
pub fn revcomp(seq: &str) -> String {
    std::str::from_utf8(&bio::alphabets::dna::revcomp(seq.as_bytes()))
        .expect("invalid utf-8 encoding")
        .to_string()
}

/// Mapping for DNA characters for normalization.
/// Built via
/// ```rust,no_run
/// let mut result = [0; 256];
/// for c in 0..=255 {
///     if c == b'u' || c == b'U' {
///         result[c as usize] = b'T';
///     } else if c.is_ascii_lowercase() {
///         result[c as usize] = c.to_ascii_uppercase();
///     } else {
///         result[c as usize] = c;
///     }
/// }
/// ```
/// Could probably be done in build.rs
const DNA_ASCII_MAP: [u8; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
    50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
    74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 84, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 65,
    66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 84, 86, 87, 88, 89,
    90, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140,
    141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
    179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
    198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216,
    217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235,
    236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
    255,
];
const DNA_ASCII_TO_2BIT: [u8; 256] = {
    let mut result = [255; 256];

    result[b'A' as usize] = 0;
    result[b'a' as usize] = 0;

    result[b'C' as usize] = 1;
    result[b'c' as usize] = 1;

    result[b'G' as usize] = 2;
    result[b'g' as usize] = 2;

    result[b'T' as usize] = 3;
    result[b't' as usize] = 3;
    result[b'U' as usize] = 3;
    result[b'u' as usize] = 3;
    result
};

fn dna3_to_2bit(c: &[u8]) -> Option<u8> {
    let mut result = 0;
    for i in 0..3 {
        result <<= 2;
        let tmp = DNA_ASCII_TO_2BIT[c[i] as usize];
        if tmp == 255 {
            return None;
        }
        result |= tmp;
    }
    Some(result)
}

pub const AA3_TO_AA1_VEC: &[(&str, &str)] = &[
    ("Ala", "A"),
    ("Arg", "R"),
    ("Asn", "N"),
    ("Asp", "D"),
    ("Cys", "C"),
    ("Gln", "Q"),
    ("Glu", "E"),
    ("Gly", "G"),
    ("His", "H"),
    ("Ile", "I"),
    ("Leu", "L"),
    ("Lys", "K"),
    ("Met", "M"),
    ("Phe", "F"),
    ("Pro", "P"),
    ("Ser", "S"),
    ("Thr", "T"),
    ("Trp", "W"),
    ("Tyr", "Y"),
    ("Val", "V"),
    ("Xaa", "X"),
    ("Ter", "*"),
    ("Sec", "U"),
];

const DNA_TO_AA1_LUT_VEC: &[(&str, &str)] = &[
    ("AAA", "K"),
    ("AAC", "N"),
    ("AAG", "K"),
    ("AAT", "N"),
    ("ACA", "T"),
    ("ACC", "T"),
    ("ACG", "T"),
    ("ACT", "T"),
    ("AGA", "R"),
    ("AGC", "S"),
    ("AGG", "R"),
    ("AGT", "S"),
    ("ATA", "I"),
    ("ATC", "I"),
    ("ATG", "M"),
    ("ATT", "I"),
    ("CAA", "Q"),
    ("CAC", "H"),
    ("CAG", "Q"),
    ("CAT", "H"),
    ("CCA", "P"),
    ("CCC", "P"),
    ("CCG", "P"),
    ("CCT", "P"),
    ("CGA", "R"),
    ("CGC", "R"),
    ("CGG", "R"),
    ("CGT", "R"),
    ("CTA", "L"),
    ("CTC", "L"),
    ("CTG", "L"),
    ("CTT", "L"),
    ("GAA", "E"),
    ("GAC", "D"),
    ("GAG", "E"),
    ("GAT", "D"),
    ("GCA", "A"),
    ("GCC", "A"),
    ("GCG", "A"),
    ("GCT", "A"),
    ("GGA", "G"),
    ("GGC", "G"),
    ("GGG", "G"),
    ("GGT", "G"),
    ("GTA", "V"),
    ("GTC", "V"),
    ("GTG", "V"),
    ("GTT", "V"),
    ("TAA", "*"),
    ("TAC", "Y"),
    ("TAG", "*"),
    ("TAT", "Y"),
    ("TCA", "S"),
    ("TCC", "S"),
    ("TCG", "S"),
    ("TCT", "S"),
    // caveat lector
    ("TGA", "*"),
    ("TGC", "C"),
    ("TGG", "W"),
    ("TGT", "C"),
    ("TTA", "L"),
    ("TTC", "F"),
    ("TTG", "L"),
    ("TTT", "F"),
    // degenerate codons
    ("AAR", "K"),
    ("AAY", "N"),
    ("ACB", "T"),
    ("ACD", "T"),
    ("ACH", "T"),
    ("ACK", "T"),
    ("ACM", "T"),
    ("ACN", "T"),
    ("ACR", "T"),
    ("ACS", "T"),
    ("ACV", "T"),
    ("ACW", "T"),
    ("ACY", "T"),
    ("AGR", "R"),
    ("AGY", "S"),
    ("ATH", "I"),
    ("ATM", "I"),
    ("ATW", "I"),
    ("ATY", "I"),
    ("CAR", "Q"),
    ("CAY", "H"),
    ("CCB", "P"),
    ("CCD", "P"),
    ("CCH", "P"),
    ("CCK", "P"),
    ("CCM", "P"),
    ("CCN", "P"),
    ("CCR", "P"),
    ("CCS", "P"),
    ("CCV", "P"),
    ("CCW", "P"),
    ("CCY", "P"),
    ("CGB", "R"),
    ("CGD", "R"),
    ("CGH", "R"),
    ("CGK", "R"),
    ("CGM", "R"),
    ("CGN", "R"),
    ("CGR", "R"),
    ("CGS", "R"),
    ("CGV", "R"),
    ("CGW", "R"),
    ("CGY", "R"),
    ("CTB", "L"),
    ("CTD", "L"),
    ("CTH", "L"),
    ("CTK", "L"),
    ("CTM", "L"),
    ("CTN", "L"),
    ("CTR", "L"),
    ("CTS", "L"),
    ("CTV", "L"),
    ("CTW", "L"),
    ("CTY", "L"),
    ("GAR", "E"),
    ("GAY", "D"),
    ("GCB", "A"),
    ("GCD", "A"),
    ("GCH", "A"),
    ("GCK", "A"),
    ("GCM", "A"),
    ("GCN", "A"),
    ("GCR", "A"),
    ("GCS", "A"),
    ("GCV", "A"),
    ("GCW", "A"),
    ("GCY", "A"),
    ("GGB", "G"),
    ("GGD", "G"),
    ("GGH", "G"),
    ("GGK", "G"),
    ("GGM", "G"),
    ("GGN", "G"),
    ("GGR", "G"),
    ("GGS", "G"),
    ("GGV", "G"),
    ("GGW", "G"),
    ("GGY", "G"),
    ("GTB", "V"),
    ("GTD", "V"),
    ("GTH", "V"),
    ("GTK", "V"),
    ("GTM", "V"),
    ("GTN", "V"),
    ("GTR", "V"),
    ("GTS", "V"),
    ("GTV", "V"),
    ("GTW", "V"),
    ("GTY", "V"),
    ("MGA", "R"),
    ("MGG", "R"),
    ("MGR", "R"),
    ("TAR", "*"),
    ("TAY", "Y"),
    ("TCB", "S"),
    ("TCD", "S"),
    ("TCH", "S"),
    ("TCK", "S"),
    ("TCM", "S"),
    ("TCN", "S"),
    ("TCR", "S"),
    ("TCS", "S"),
    ("TCV", "S"),
    ("TCW", "S"),
    ("TCY", "S"),
    ("TGY", "C"),
    ("TRA", "*"),
    ("TTR", "L"),
    ("TTY", "F"),
    ("YTA", "L"),
    ("YTG", "L"),
    ("YTR", "L"),
];

/// Translation table for selenocysteine.
const DNA_TO_AA1_SEC_VEC: &[(&str, &str)] = &[
    ("AAA", "K"),
    ("AAC", "N"),
    ("AAG", "K"),
    ("AAT", "N"),
    ("ACA", "T"),
    ("ACC", "T"),
    ("ACG", "T"),
    ("ACT", "T"),
    ("AGA", "R"),
    ("AGC", "S"),
    ("AGG", "R"),
    ("AGT", "S"),
    ("ATA", "I"),
    ("ATC", "I"),
    ("ATG", "M"),
    ("ATT", "I"),
    ("CAA", "Q"),
    ("CAC", "H"),
    ("CAG", "Q"),
    ("CAT", "H"),
    ("CCA", "P"),
    ("CCC", "P"),
    ("CCG", "P"),
    ("CCT", "P"),
    ("CGA", "R"),
    ("CGC", "R"),
    ("CGG", "R"),
    ("CGT", "R"),
    ("CTA", "L"),
    ("CTC", "L"),
    ("CTG", "L"),
    ("CTT", "L"),
    ("GAA", "E"),
    ("GAC", "D"),
    ("GAG", "E"),
    ("GAT", "D"),
    ("GCA", "A"),
    ("GCC", "A"),
    ("GCG", "A"),
    ("GCT", "A"),
    ("GGA", "G"),
    ("GGC", "G"),
    ("GGG", "G"),
    ("GGT", "G"),
    ("GTA", "V"),
    ("GTC", "V"),
    ("GTG", "V"),
    ("GTT", "V"),
    ("TAA", "*"),
    ("TAC", "Y"),
    ("TAG", "*"),
    ("TAT", "Y"),
    ("TCA", "S"),
    ("TCC", "S"),
    ("TCG", "S"),
    ("TCT", "S"),
    // caveat lector
    ("TGA", "U"),
    ("TGC", "C"),
    ("TGG", "W"),
    ("TGT", "C"),
    ("TTA", "L"),
    ("TTC", "F"),
    ("TTG", "L"),
    ("TTT", "F"),
    // degenerate codons
    ("AAR", "K"),
    ("AAY", "N"),
    ("ACB", "T"),
    ("ACD", "T"),
    ("ACH", "T"),
    ("ACK", "T"),
    ("ACM", "T"),
    ("ACN", "T"),
    ("ACR", "T"),
    ("ACS", "T"),
    ("ACV", "T"),
    ("ACW", "T"),
    ("ACY", "T"),
    ("AGR", "R"),
    ("AGY", "S"),
    ("ATH", "I"),
    ("ATM", "I"),
    ("ATW", "I"),
    ("ATY", "I"),
    ("CAR", "Q"),
    ("CAY", "H"),
    ("CCB", "P"),
    ("CCD", "P"),
    ("CCH", "P"),
    ("CCK", "P"),
    ("CCM", "P"),
    ("CCN", "P"),
    ("CCR", "P"),
    ("CCS", "P"),
    ("CCV", "P"),
    ("CCW", "P"),
    ("CCY", "P"),
    ("CGB", "R"),
    ("CGD", "R"),
    ("CGH", "R"),
    ("CGK", "R"),
    ("CGM", "R"),
    ("CGN", "R"),
    ("CGR", "R"),
    ("CGS", "R"),
    ("CGV", "R"),
    ("CGW", "R"),
    ("CGY", "R"),
    ("CTB", "L"),
    ("CTD", "L"),
    ("CTH", "L"),
    ("CTK", "L"),
    ("CTM", "L"),
    ("CTN", "L"),
    ("CTR", "L"),
    ("CTS", "L"),
    ("CTV", "L"),
    ("CTW", "L"),
    ("CTY", "L"),
    ("GAR", "E"),
    ("GAY", "D"),
    ("GCB", "A"),
    ("GCD", "A"),
    ("GCH", "A"),
    ("GCK", "A"),
    ("GCM", "A"),
    ("GCN", "A"),
    ("GCR", "A"),
    ("GCS", "A"),
    ("GCV", "A"),
    ("GCW", "A"),
    ("GCY", "A"),
    ("GGB", "G"),
    ("GGD", "G"),
    ("GGH", "G"),
    ("GGK", "G"),
    ("GGM", "G"),
    ("GGN", "G"),
    ("GGR", "G"),
    ("GGS", "G"),
    ("GGV", "G"),
    ("GGW", "G"),
    ("GGY", "G"),
    ("GTB", "V"),
    ("GTD", "V"),
    ("GTH", "V"),
    ("GTK", "V"),
    ("GTM", "V"),
    ("GTN", "V"),
    ("GTR", "V"),
    ("GTS", "V"),
    ("GTV", "V"),
    ("GTW", "V"),
    ("GTY", "V"),
    ("MGA", "R"),
    ("MGG", "R"),
    ("MGR", "R"),
    ("TAR", "*"),
    ("TAY", "Y"),
    ("TCB", "S"),
    ("TCD", "S"),
    ("TCH", "S"),
    ("TCK", "S"),
    ("TCM", "S"),
    ("TCN", "S"),
    ("TCR", "S"),
    ("TCS", "S"),
    ("TCV", "S"),
    ("TCW", "S"),
    ("TCY", "S"),
    ("TGY", "C"),
    ("TRA", "*"),
    ("TTR", "L"),
    ("TTY", "F"),
    ("YTA", "L"),
    ("YTG", "L"),
    ("YTR", "L"),
];

/// Vertebrate mitochondrial code, cf. https://en.wikipedia.org/wiki/Vertebrate_mitochondrial_code
const DNA_TO_AA1_CHRMT_VERTEBRATE_VEC: &[(&str, &str)] = &[
    ("AAA", "K"),
    ("AAC", "N"),
    ("AAG", "K"),
    ("AAT", "N"),
    ("ACA", "T"),
    ("ACC", "T"),
    ("ACG", "T"),
    ("ACT", "T"),
    // caveat lector
    ("AGA", "*"),
    ("AGC", "S"),
    // caveat lector
    ("AGG", "*"),
    ("AGT", "S"),
    // caveat lector
    ("ATA", "M"),
    ("ATC", "I"),
    ("ATG", "M"),
    ("ATT", "I"),
    ("CAA", "Q"),
    ("CAC", "H"),
    ("CAG", "Q"),
    ("CAT", "H"),
    ("CCA", "P"),
    ("CCC", "P"),
    ("CCG", "P"),
    ("CCT", "P"),
    ("CGA", "R"),
    ("CGC", "R"),
    ("CGG", "R"),
    ("CGT", "R"),
    ("CTA", "L"),
    ("CTC", "L"),
    ("CTG", "L"),
    ("CTT", "L"),
    ("GAA", "E"),
    ("GAC", "D"),
    ("GAG", "E"),
    ("GAT", "D"),
    ("GCA", "A"),
    ("GCC", "A"),
    ("GCG", "A"),
    ("GCT", "A"),
    ("GGA", "G"),
    ("GGC", "G"),
    ("GGG", "G"),
    ("GGT", "G"),
    ("GTA", "V"),
    ("GTC", "V"),
    ("GTG", "V"),
    ("GTT", "V"),
    ("TAA", "*"),
    ("TAC", "Y"),
    ("TAG", "*"),
    ("TAT", "Y"),
    ("TCA", "S"),
    ("TCC", "S"),
    ("TCG", "S"),
    ("TCT", "S"),
    // caveat lector
    ("TGA", "W"),
    ("TGC", "C"),
    ("TGG", "W"),
    ("TGT", "C"),
    ("TTA", "L"),
    ("TTC", "F"),
    ("TTG", "L"),
    ("TTT", "F"),
    // degenerate codons
    ("AAR", "K"),
    ("AAY", "N"),
    ("ACB", "T"),
    ("ACD", "T"),
    ("ACH", "T"),
    ("ACK", "T"),
    ("ACM", "T"),
    ("ACN", "T"),
    ("ACR", "T"),
    ("ACS", "T"),
    ("ACV", "T"),
    ("ACW", "T"),
    ("ACY", "T"),
    ("AGR", "R"),
    ("AGY", "S"),
    ("ATH", "I"),
    ("ATM", "I"),
    ("ATW", "I"),
    ("ATY", "I"),
    ("CAR", "Q"),
    ("CAY", "H"),
    ("CCB", "P"),
    ("CCD", "P"),
    ("CCH", "P"),
    ("CCK", "P"),
    ("CCM", "P"),
    ("CCN", "P"),
    ("CCR", "P"),
    ("CCS", "P"),
    ("CCV", "P"),
    ("CCW", "P"),
    ("CCY", "P"),
    ("CGB", "R"),
    ("CGD", "R"),
    ("CGH", "R"),
    ("CGK", "R"),
    ("CGM", "R"),
    ("CGN", "R"),
    ("CGR", "R"),
    ("CGS", "R"),
    ("CGV", "R"),
    ("CGW", "R"),
    ("CGY", "R"),
    ("CTB", "L"),
    ("CTD", "L"),
    ("CTH", "L"),
    ("CTK", "L"),
    ("CTM", "L"),
    ("CTN", "L"),
    ("CTR", "L"),
    ("CTS", "L"),
    ("CTV", "L"),
    ("CTW", "L"),
    ("CTY", "L"),
    ("GAR", "E"),
    ("GAY", "D"),
    ("GCB", "A"),
    ("GCD", "A"),
    ("GCH", "A"),
    ("GCK", "A"),
    ("GCM", "A"),
    ("GCN", "A"),
    ("GCR", "A"),
    ("GCS", "A"),
    ("GCV", "A"),
    ("GCW", "A"),
    ("GCY", "A"),
    ("GGB", "G"),
    ("GGD", "G"),
    ("GGH", "G"),
    ("GGK", "G"),
    ("GGM", "G"),
    ("GGN", "G"),
    ("GGR", "G"),
    ("GGS", "G"),
    ("GGV", "G"),
    ("GGW", "G"),
    ("GGY", "G"),
    ("GTB", "V"),
    ("GTD", "V"),
    ("GTH", "V"),
    ("GTK", "V"),
    ("GTM", "V"),
    ("GTN", "V"),
    ("GTR", "V"),
    ("GTS", "V"),
    ("GTV", "V"),
    ("GTW", "V"),
    ("GTY", "V"),
    ("MGA", "R"),
    ("MGG", "R"),
    ("MGR", "R"),
    ("TAR", "*"),
    ("TAY", "Y"),
    ("TCB", "S"),
    ("TCD", "S"),
    ("TCH", "S"),
    ("TCK", "S"),
    ("TCM", "S"),
    ("TCN", "S"),
    ("TCR", "S"),
    ("TCS", "S"),
    ("TCV", "S"),
    ("TCW", "S"),
    ("TCY", "S"),
    ("TGY", "C"),
    ("TRA", "*"),
    ("TTR", "L"),
    ("TTY", "F"),
    ("YTA", "L"),
    ("YTG", "L"),
    ("YTR", "L"),
];

/// Generated via:
/// ```rust,no_run
/// const _: &str = stringify!{
///     let mut result = [0; 64];
///     for (i, (dna3, aa1)) in DNA_TO_AA1_LUT_VEC.iter().enumerate() {
///         if i > 63 {
///             break;  // skip degenerate codons
///         }
///         let dna3_2bit = dna3_to_2bit(dna3.as_bytes()).expect("invalid dna3");
///         result[dna3_2bit as usize] = aa1.as_bytes()[0];
///     }
/// };
/// ```
///
const CODON_2BIT_TO_AA1_LUT: [u8; 64] = [
    75, 78, 75, 78, 84, 84, 84, 84, 82, 83, 82, 83, 73, 73, 77, 73, 81, 72, 81, 72, 80, 80, 80, 80,
    82, 82, 82, 82, 76, 76, 76, 76, 69, 68, 69, 68, 65, 65, 65, 65, 71, 71, 71, 71, 86, 86, 86, 86,
    42, 89, 42, 89, 83, 83, 83, 83, 42, 67, 87, 67, 76, 70, 76, 70,
];

lazy_static::lazy_static! {
    static ref AA1_TO_AA3: AHashMap<&'static [u8], &'static str> = {
        let mut m = AHashMap::default();
        for (aa3, aa1) in AA3_TO_AA1_VEC.iter() {
            m.insert(aa1.as_bytes(), *aa3);
        }
        m
    };

    static ref AA3_TO_AA1: AHashMap<&'static [u8], &'static str> = {
        let mut m = AHashMap::default();
        for (aa3, aa1) in AA3_TO_AA1_VEC.iter() {
            m.insert(aa3.as_bytes(), *aa1);
        }
        m
    };

    static ref DNA_TO_AA1_LUT: AHashMap<Codon, u8> = {
        let mut m = AHashMap::default();
        for (dna, aa1) in DNA_TO_AA1_LUT_VEC {
            assert_eq!(dna.len(), 3);
            let d = dna.as_bytes();
            m.insert([d[0], d[1], d[2]], aa1.as_bytes()[0]);
        }
        m
    };

    static ref DNA_TO_AA1_SEC: AHashMap<Codon, u8> = {
        let mut m = AHashMap::default();
        for (dna, aa1) in DNA_TO_AA1_SEC_VEC {
            assert_eq!(dna.len(), 3);
            let d = dna.as_bytes();
            m.insert([d[0], d[1], d[2]], aa1.as_bytes()[0]);
        }
        m
    };

    static ref DNA_TO_AA1_CHRMT_VERTEBRATE: AHashMap<Codon, u8> = {
        let mut m = AHashMap::default();
        for (dna, aa1) in DNA_TO_AA1_CHRMT_VERTEBRATE_VEC {
            assert_eq!(dna.len(), 3);
            let d = dna.as_bytes();
            m.insert([d[0], d[1], d[2]], aa1.as_bytes()[0]);
        }
        m
    };

    static ref CODON_2BIT_TO_AA1_SEC: [u8; 64] = {
        let mut result = [0; 64];
        for (i, (dna3, aa1)) in DNA_TO_AA1_SEC_VEC.iter().enumerate() {
            if i > 63 {
                break;  // skip degenerate codons
            }
            let dna3_2bit = dna3_to_2bit(dna3.as_bytes()).expect("invalid dna3");
            result[dna3_2bit as usize] = aa1.as_bytes()[0];
        }
        result
    };

    static ref CODON_2BIT_TO_AA1_CHRMT_VERTEBRATE: [u8; 64] = {
        let mut result = [0; 64];
        for (i, (dna3, aa1)) in DNA_TO_AA1_CHRMT_VERTEBRATE_VEC.iter().enumerate() {
            if i > 63 {
                break;  // skip degenerate codons
            }
            let dna3_2bit = dna3_to_2bit(dna3.as_bytes()).expect("invalid dna3");
            result[dna3_2bit as usize] = aa1.as_bytes()[0];
        }
        result
    };
}

const IUPAC_AMBIGUITY_CODES: [u8; 13] = *b"BDHVNUWSMKRYZ";

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

    for (i, aa1) in seq.as_bytes().chunks(1).enumerate() {
        let aa3 = AA1_TO_AA3.get(aa1).ok_or_else(|| {
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
        let aa1 = AA3_TO_AA1.get(aa3).ok_or(Error::InvalidOneLetterAminoAcid(
            format!("{:?}", aa3),
            format!("{}", i + 1),
        ))?;
        result.push_str(aa1);
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

impl CodonTranslator {
    /// Initialize the struct.
    pub fn new(table: TranslationTable) -> Self {
        Self {
            dna_ascii_map: DNA_ASCII_MAP,
            dna_ascii_to_2bit: DNA_ASCII_TO_2BIT,
            iupac_ambiguity_codes: IUPAC_AMBIGUITY_CODES,

            codon_2bit_to_aa1: match table {
                TranslationTable::Standard => CODON_2BIT_TO_AA1_LUT,
                TranslationTable::Selenocysteine => *CODON_2BIT_TO_AA1_SEC,
                TranslationTable::VertebrateMitochondrial => *CODON_2BIT_TO_AA1_CHRMT_VERTEBRATE,
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
    seq: &str,
    full_codons: bool,
    ter_symbol: &str,
    translation_table: TranslationTable,
) -> Result<String, Error> {
    if seq.is_empty() {
        return Ok("".to_string());
    }

    if full_codons && seq.len() % 3 != 0 {
        return Err(Error::UntranslatableDnaLenth(seq.len()));
    }

    // Translate the codons from the input to result.
    let mut translator = CodonTranslator::new(translation_table);
    let mut result = String::with_capacity(seq.len() / 3);
    for chunk in seq.as_bytes().chunks_exact(3) {
        result.push(char::from(translator.translate(chunk)?));
    }

    // Check for trailing bases and add the ter symbol if required.
    if !full_codons && seq.len() % 3 != 0 {
        result.push_str(ter_symbol);
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

    #[test]
    fn exercise_lazy_ds() {
        assert!(DNA_ASCII_MAP[0] == b'\0');
        assert!(DNA_ASCII_TO_2BIT[b'A' as usize] == 0);
        assert!(AA3_TO_AA1_VEC[0] == ("Ala", "A"));
        assert!(DNA_TO_AA1_LUT_VEC[0] == ("AAA", "K"));
        assert!(DNA_TO_AA1_SEC_VEC[0] == ("AAA", "K"));
        assert!(DNA_TO_AA1_CHRMT_VERTEBRATE_VEC[0] == ("AAA", "K"));
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
