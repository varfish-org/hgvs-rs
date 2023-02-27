//! Utility code for working with sequences.
//!
//! Partially ported over from `bioutils.sequences`.

use md5::{Digest, Md5};
use phf::phf_map;

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
        .unwrap()
        .to_string()
}

static AA3_TO_AA1_LUT: phf::Map<&'static str, &'static str> = phf_map! {
    "Ala" => "A",
    "Arg" => "R",
    "Asn" => "N",
    "Asp" => "D",
    "Cys" => "C",
    "Gln" => "Q",
    "Glu" => "E",
    "Gly" => "G",
    "His" => "H",
    "Ile" => "I",
    "Leu" => "L",
    "Lys" => "K",
    "Met" => "M",
    "Phe" => "F",
    "Pro" => "P",
    "Ser" => "S",
    "Thr" => "T",
    "Trp" => "W",
    "Tyr" => "Y",
    "Val" => "V",
    "Xaa" => "X",
    "Ter" => "*",
    "Sec" => "U",
};

static AA1_TO_AA3_LUT: phf::Map<&'static str, &'static str> = phf_map! {
    "A" => "Ala",
    "R" => "Arg",
    "N" => "Asn",
    "D" => "Asp",
    "C" => "Cys",
    "Q" => "Gln",
    "E" => "Glu",
    "G" => "Gly",
    "H" => "His",
    "I" => "Ile",
    "L" => "Leu",
    "K" => "Lys",
    "M" => "Met",
    "F" => "Phe",
    "P" => "Pro",
    "S" => "Ser",
    "T" => "Thr",
    "W" => "Trp",
    "Y" => "Tyr",
    "V" => "Val",
    "X" => "Xaa",
    "*" => "Ter",
    "U" => "Sec",
};

/// NCBI standard translation table.
static DNA_TO_AA1_LUT: phf::Map<&'static str, &'static str> = phf_map! {
    "AAA" => "K",
    "AAC" => "N",
    "AAG" => "K",
    "AAT" => "N",
    "ACA" => "T",
    "ACC" => "T",
    "ACG" => "T",
    "ACT" => "T",
    "AGA" => "R",
    "AGC" => "S",
    "AGG" => "R",
    "AGT" => "S",
    "ATA" => "I",
    "ATC" => "I",
    "ATG" => "M",
    "ATT" => "I",
    "CAA" => "Q",
    "CAC" => "H",
    "CAG" => "Q",
    "CAT" => "H",
    "CCA" => "P",
    "CCC" => "P",
    "CCG" => "P",
    "CCT" => "P",
    "CGA" => "R",
    "CGC" => "R",
    "CGG" => "R",
    "CGT" => "R",
    "CTA" => "L",
    "CTC" => "L",
    "CTG" => "L",
    "CTT" => "L",
    "GAA" => "E",
    "GAC" => "D",
    "GAG" => "E",
    "GAT" => "D",
    "GCA" => "A",
    "GCC" => "A",
    "GCG" => "A",
    "GCT" => "A",
    "GGA" => "G",
    "GGC" => "G",
    "GGG" => "G",
    "GGT" => "G",
    "GTA" => "V",
    "GTC" => "V",
    "GTG" => "V",
    "GTT" => "V",
    "TAA" => "*",
    "TAC" => "Y",
    "TAG" => "*",
    "TAT" => "Y",
    "TCA" => "S",
    "TCC" => "S",
    "TCG" => "S",
    "TCT" => "S",
    // caveat lector
    "TGA" => "*",
    "TGC" => "C",
    "TGG" => "W",
    "TGT" => "C",
    "TTA" => "L",
    "TTC" => "F",
    "TTG" => "L",
    "TTT" => "F",
    // degenerate codons
    "AAR" => "K",
    "AAY" => "N",
    "ACB" => "T",
    "ACD" => "T",
    "ACH" => "T",
    "ACK" => "T",
    "ACM" => "T",
    "ACN" => "T",
    "ACR" => "T",
    "ACS" => "T",
    "ACV" => "T",
    "ACW" => "T",
    "ACY" => "T",
    "AGR" => "R",
    "AGY" => "S",
    "ATH" => "I",
    "ATM" => "I",
    "ATW" => "I",
    "ATY" => "I",
    "CAR" => "Q",
    "CAY" => "H",
    "CCB" => "P",
    "CCD" => "P",
    "CCH" => "P",
    "CCK" => "P",
    "CCM" => "P",
    "CCN" => "P",
    "CCR" => "P",
    "CCS" => "P",
    "CCV" => "P",
    "CCW" => "P",
    "CCY" => "P",
    "CGB" => "R",
    "CGD" => "R",
    "CGH" => "R",
    "CGK" => "R",
    "CGM" => "R",
    "CGN" => "R",
    "CGR" => "R",
    "CGS" => "R",
    "CGV" => "R",
    "CGW" => "R",
    "CGY" => "R",
    "CTB" => "L",
    "CTD" => "L",
    "CTH" => "L",
    "CTK" => "L",
    "CTM" => "L",
    "CTN" => "L",
    "CTR" => "L",
    "CTS" => "L",
    "CTV" => "L",
    "CTW" => "L",
    "CTY" => "L",
    "GAR" => "E",
    "GAY" => "D",
    "GCB" => "A",
    "GCD" => "A",
    "GCH" => "A",
    "GCK" => "A",
    "GCM" => "A",
    "GCN" => "A",
    "GCR" => "A",
    "GCS" => "A",
    "GCV" => "A",
    "GCW" => "A",
    "GCY" => "A",
    "GGB" => "G",
    "GGD" => "G",
    "GGH" => "G",
    "GGK" => "G",
    "GGM" => "G",
    "GGN" => "G",
    "GGR" => "G",
    "GGS" => "G",
    "GGV" => "G",
    "GGW" => "G",
    "GGY" => "G",
    "GTB" => "V",
    "GTD" => "V",
    "GTH" => "V",
    "GTK" => "V",
    "GTM" => "V",
    "GTN" => "V",
    "GTR" => "V",
    "GTS" => "V",
    "GTV" => "V",
    "GTW" => "V",
    "GTY" => "V",
    "MGA" => "R",
    "MGG" => "R",
    "MGR" => "R",
    "TAR" => "*",
    "TAY" => "Y",
    "TCB" => "S",
    "TCD" => "S",
    "TCH" => "S",
    "TCK" => "S",
    "TCM" => "S",
    "TCN" => "S",
    "TCR" => "S",
    "TCS" => "S",
    "TCV" => "S",
    "TCW" => "S",
    "TCY" => "S",
    "TGY" => "C",
    "TRA" => "*",
    "TTR" => "L",
    "TTY" => "F",
    "YTA" => "L",
    "YTG" => "L",
    "YTR" => "L",
};

/// Translation table for selenocysteine.
static DNA_TO_AA1_SEC: phf::Map<&'static str, &'static str> = phf_map! {
    "AAA" => "K",
    "AAC" => "N",
    "AAG" => "K",
    "AAT" => "N",
    "ACA" => "T",
    "ACC" => "T",
    "ACG" => "T",
    "ACT" => "T",
    "AGA" => "R",
    "AGC" => "S",
    "AGG" => "R",
    "AGT" => "S",
    "ATA" => "I",
    "ATC" => "I",
    "ATG" => "M",
    "ATT" => "I",
    "CAA" => "Q",
    "CAC" => "H",
    "CAG" => "Q",
    "CAT" => "H",
    "CCA" => "P",
    "CCC" => "P",
    "CCG" => "P",
    "CCT" => "P",
    "CGA" => "R",
    "CGC" => "R",
    "CGG" => "R",
    "CGT" => "R",
    "CTA" => "L",
    "CTC" => "L",
    "CTG" => "L",
    "CTT" => "L",
    "GAA" => "E",
    "GAC" => "D",
    "GAG" => "E",
    "GAT" => "D",
    "GCA" => "A",
    "GCC" => "A",
    "GCG" => "A",
    "GCT" => "A",
    "GGA" => "G",
    "GGC" => "G",
    "GGG" => "G",
    "GGT" => "G",
    "GTA" => "V",
    "GTC" => "V",
    "GTG" => "V",
    "GTT" => "V",
    "TAA" => "*",
    "TAC" => "Y",
    "TAG" => "*",
    "TAT" => "Y",
    "TCA" => "S",
    "TCC" => "S",
    "TCG" => "S",
    "TCT" => "S",
    // caveat lector
    "TGA" => "U",
    "TGC" => "C",
    "TGG" => "W",
    "TGT" => "C",
    "TTA" => "L",
    "TTC" => "F",
    "TTG" => "L",
    "TTT" => "F",
    // degenerate codons
    "AAR" => "K",
    "AAY" => "N",
    "ACB" => "T",
    "ACD" => "T",
    "ACH" => "T",
    "ACK" => "T",
    "ACM" => "T",
    "ACN" => "T",
    "ACR" => "T",
    "ACS" => "T",
    "ACV" => "T",
    "ACW" => "T",
    "ACY" => "T",
    "AGR" => "R",
    "AGY" => "S",
    "ATH" => "I",
    "ATM" => "I",
    "ATW" => "I",
    "ATY" => "I",
    "CAR" => "Q",
    "CAY" => "H",
    "CCB" => "P",
    "CCD" => "P",
    "CCH" => "P",
    "CCK" => "P",
    "CCM" => "P",
    "CCN" => "P",
    "CCR" => "P",
    "CCS" => "P",
    "CCV" => "P",
    "CCW" => "P",
    "CCY" => "P",
    "CGB" => "R",
    "CGD" => "R",
    "CGH" => "R",
    "CGK" => "R",
    "CGM" => "R",
    "CGN" => "R",
    "CGR" => "R",
    "CGS" => "R",
    "CGV" => "R",
    "CGW" => "R",
    "CGY" => "R",
    "CTB" => "L",
    "CTD" => "L",
    "CTH" => "L",
    "CTK" => "L",
    "CTM" => "L",
    "CTN" => "L",
    "CTR" => "L",
    "CTS" => "L",
    "CTV" => "L",
    "CTW" => "L",
    "CTY" => "L",
    "GAR" => "E",
    "GAY" => "D",
    "GCB" => "A",
    "GCD" => "A",
    "GCH" => "A",
    "GCK" => "A",
    "GCM" => "A",
    "GCN" => "A",
    "GCR" => "A",
    "GCS" => "A",
    "GCV" => "A",
    "GCW" => "A",
    "GCY" => "A",
    "GGB" => "G",
    "GGD" => "G",
    "GGH" => "G",
    "GGK" => "G",
    "GGM" => "G",
    "GGN" => "G",
    "GGR" => "G",
    "GGS" => "G",
    "GGV" => "G",
    "GGW" => "G",
    "GGY" => "G",
    "GTB" => "V",
    "GTD" => "V",
    "GTH" => "V",
    "GTK" => "V",
    "GTM" => "V",
    "GTN" => "V",
    "GTR" => "V",
    "GTS" => "V",
    "GTV" => "V",
    "GTW" => "V",
    "GTY" => "V",
    "MGA" => "R",
    "MGG" => "R",
    "MGR" => "R",
    "TAR" => "*",
    "TAY" => "Y",
    "TCB" => "S",
    "TCD" => "S",
    "TCH" => "S",
    "TCK" => "S",
    "TCM" => "S",
    "TCN" => "S",
    "TCR" => "S",
    "TCS" => "S",
    "TCV" => "S",
    "TCW" => "S",
    "TCY" => "S",
    "TGY" => "C",
    "TRA" => "*",
    "TTR" => "L",
    "TTY" => "F",
    "YTA" => "L",
    "YTG" => "L",
    "YTR" => "L",
};

static IUPAC_AMBIGUITY_CODES: &str = "BDHVNUWSMKRYZ";

/// Allow selection of translation table.
pub enum TranslationTable {
    Standard,
    #[allow(dead_code)]
    Selenocysteine,
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
pub fn aa_to_aa1(seq: &str) -> Result<String, anyhow::Error> {
    if looks_like_aa3_p(seq) {
        aa3_to_aa1(seq)
    } else {
        Ok(seq.to_string())
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
pub fn aa1_to_aa3(seq: &str) -> Result<String, anyhow::Error> {
    let mut result = String::new();

    for i in 0..seq.len() {
        let aa3 = AA1_TO_AA3_LUT.get(&seq[i..(i + 1)]).ok_or(anyhow::anyhow!(
            "Invalid 1-letter amino acid: {}",
            &seq[i..(i + 1)]
        ))?;
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
pub fn aa3_to_aa1(seq: &str) -> Result<String, anyhow::Error> {
    if seq.len() % 3 != 0 {
        return Err(anyhow::anyhow!(
            "3-letter amino acid sequence length is not multiple of three"
        ));
    }

    let mut result = String::new();

    for i in 0..(seq.len() / 3) {
        let aa3 = &seq[(i * 3)..((i + 1) * 3)];
        let aa1 = AA3_TO_AA1_LUT
            .get(aa3)
            .ok_or(anyhow::anyhow!("Invalid 3-letter amino acid: {}", &aa3,))?;
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
) -> Result<String, anyhow::Error> {
    if seq.is_empty() {
        return Ok("".to_string());
    }

    if full_codons && seq.len() % 3 != 0 {
        return Err(anyhow::anyhow!(
            "Sequence length must be a multiple of three"
        ));
    }

    let translation_table = match translation_table {
        TranslationTable::Standard => &DNA_TO_AA1_LUT,
        TranslationTable::Selenocysteine => &DNA_TO_AA1_SEC,
    };

    let seq = seq.replace('u', "t").replace('U', "T").to_uppercase();

    let mut result = String::new();
    for i in 0..(seq.len() / 3) {
        let codon = &seq[(i * 3)..((i + 1) * 3)];
        let aa = translation_table.get(codon);
        if let Some(aa) = aa {
            result.push_str(aa);
        } else {
            // If this contains an ambiguous code, set aa to X, otherwise, throw error
            let mut ok = false;
            for i in 0..codon.len() {
                if IUPAC_AMBIGUITY_CODES.contains(&codon[i..(i + 1)]) {
                    ok = true;
                    result.push('X');
                    break;
                }
            }
            if !ok {
                return Err(anyhow::anyhow!(
                    "Codon {} at position {}..{} is undefined in codon table",
                    codon,
                    i + 1,
                    i + 3
                ));
            }
        }
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
pub fn normalize_sequence(seq: &str) -> Result<String, anyhow::Error> {
    let mut result = String::new();

    for c in seq.chars() {
        if !c.is_whitespace() && c != '*' {
            let c = c.to_ascii_uppercase();
            if c.is_alphabetic() {
                result.push(c)
            } else {
                return Err(anyhow::anyhow!("Character {} is not alphetic", c));
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
pub fn seq_md5(seq: &str, normalize: bool) -> Result<String, anyhow::Error> {
    let seq = if normalize {
        normalize_sequence(seq)?
    } else {
        seq.to_owned()
    };
    let mut hasher = Md5::new();
    hasher.update(seq);
    let hash = hasher.finalize();
    let mut buf = [0u8; 64];
    let checksum = base16ct::lower::encode_str(&hash, &mut buf).unwrap();
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
    fn aa_to_aa1_examples() -> Result<(), anyhow::Error> {
        assert_eq!(aa_to_aa1("")?, "");
        assert_eq!(aa_to_aa1("CATSARELAME")?, "CATSARELAME");
        assert_eq!(
            aa_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")?,
            "CATSARELAME"
        );

        Ok(())
    }

    #[test]
    fn aa1_to_aa3_examples() -> Result<(), anyhow::Error> {
        assert_eq!(aa1_to_aa3("")?, "");
        assert_eq!(
            aa1_to_aa3("CATSARELAME")?,
            "CysAlaThrSerAlaArgGluLeuAlaMetGlu"
        );

        Ok(())
    }

    #[test]
    fn aa3_to_aa1_examples() -> Result<(), anyhow::Error> {
        assert!(aa3_to_aa1("Te").is_err());
        assert_eq!(aa3_to_aa1("")?, "");
        assert_eq!(
            aa3_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")?,
            "CATSARELAME"
        );

        Ok(())
    }

    #[test]
    fn translate_cds_examples() -> Result<(), anyhow::Error> {
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
    fn seq_md5_examples() -> Result<(), anyhow::Error> {
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
    fn normalize_sequence_examples() -> Result<(), anyhow::Error> {
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
