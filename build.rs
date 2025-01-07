use std::env;
use std::fs::File;
use std::io::{BufWriter, Result, Write};
use std::path::Path;

fn main() -> Result<()> {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("tables_gen.rs");
    let mut f = File::create(&dest_path).map(BufWriter::new)?;

    include_hardcoded_translation_tables(&mut f)?;
    generate_dna_ascii_map(&mut f)?;

    generate_codon_2bit_to_aa1_lut(&mut f)?;
    generate_codon_2bit_to_aa1_sec(&mut f)?;
    generate_codon_2bit_to_aa1_chrmt_vertebrate(&mut f)?;

    generate_aa1_to_aa3_str_lookup_function(&mut f)?;
    generate_aa3_to_aa1_lookup_function(&mut f)?;

    f.flush()?;
    println!("cargo::rerun-if-changed=build.rs");
    Ok(())
}

fn generate_dna_ascii_map(f: &mut BufWriter<File>) -> Result<()> {
    let mut result = [0; 256];
    for c in 0..=255 {
        if c == b'u' || c == b'U' {
            result[c as usize] = b'T';
        } else if c.is_ascii_lowercase() {
            result[c as usize] = c.to_ascii_uppercase();
        } else {
            result[c as usize] = c;
        }
    }

    writeln!(f, "/// Mapping for DNA characters for normalization.")?;
    write!(f, "const DNA_ASCII_MAP: [u8; 256] = [")?;
    for v in result {
        write!(f, "{}, ", v)?;
    }
    writeln!(f, "];")?;
    Ok(())
}
fn generate_codon_2bit_to_aa1_lut(f: &mut BufWriter<File>) -> Result<()> {
    let mut result = [0; 64];
    for (i, (dna3, aa1)) in DNA_TO_AA1_LUT_VEC.iter().enumerate() {
        if i > 63 {
            break; // skip degenerate codons
        }
        let dna3_2bit = dna3_to_2bit(dna3.as_bytes()).expect("invalid dna3");
        result[dna3_2bit as usize] = aa1.as_bytes()[0];
    }
    write!(f, "const CODON_2BIT_TO_AA1_LUT: [u8; 64] = [")?;
    for v in result {
        write!(f, "{}, ", v)?;
    }
    writeln!(f, "];")?;
    Ok(())
}

fn generate_codon_2bit_to_aa1_sec(f: &mut BufWriter<File>) -> Result<()> {
    let mut result = [0; 64];
    for (i, (dna3, aa1)) in DNA_TO_AA1_SEC_VEC.iter().enumerate() {
        if i > 63 {
            break; // skip degenerate codons
        }
        let dna3_2bit = dna3_to_2bit(dna3.as_bytes()).expect("invalid dna3");
        result[dna3_2bit as usize] = aa1.as_bytes()[0];
    }
    write!(f, "const CODON_2BIT_TO_AA1_SEC: [u8; 64] = [")?;
    for v in result {
        write!(f, "{}, ", v)?;
    }
    writeln!(f, "];")?;
    Ok(())
}

fn generate_codon_2bit_to_aa1_chrmt_vertebrate(f: &mut BufWriter<File>) -> Result<()> {
    let mut result = [0; 64];
    for (i, (dna3, aa1)) in DNA_TO_AA1_CHRMT_VERTEBRATE_VEC.iter().enumerate() {
        if i > 63 {
            break; // skip degenerate codons
        }
        let dna3_2bit = dna3_to_2bit(dna3.as_bytes()).expect("invalid dna3");
        result[dna3_2bit as usize] = aa1.as_bytes()[0];
    }
    write!(f, "const CODON_2BIT_TO_AA1_CHRMT_VERTEBRATE: [u8; 64] = [")?;
    for v in result {
        write!(f, "{}, ", v)?;
    }
    writeln!(f, "];")?;
    Ok(())
}

fn generate_aa1_to_aa3_str_lookup_function(f: &mut BufWriter<File>) -> Result<()> {
    writeln!(f, "fn _aa1_to_aa3_str(aa1: u8) -> Option<&'static str> {{")?;
    writeln!(f, "    match aa1 {{")?;
    for (aa3, aa1) in AA3_TO_AA1_VEC {
        writeln!(f, "        b'{}' => Some(\"{}\"),", aa1, aa3)?;
    }
    writeln!(f, r"        _ => None,")?;
    writeln!(f, "    }}")?;
    writeln!(f, "}}")?;
    Ok(())
}

fn generate_aa3_to_aa1_lookup_function(f: &mut BufWriter<File>) -> Result<()> {
    writeln!(f, "fn _aa3_to_aa1(aa3: &[u8]) -> Option<u8> {{")?;
    writeln!(f, "    match aa3 {{")?;
    for (aa3, aa1) in AA3_TO_AA1_VEC {
        writeln!(f, "        b\"{}\" => Some(b'{}'),", aa3, aa1)?;
    }
    writeln!(f, "        _ => None,")?;
    writeln!(f, "    }}")?;
    writeln!(f, "}}")?;
    Ok(())
}

fn include_hardcoded_translation_tables(f: &mut BufWriter<File>) -> Result<()> {
    let text = include_str!("tables.in");
    writeln!(f, "{}", text)?;
    Ok(())
}

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

// Hard-coded translation tables from src/tables.rs

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
