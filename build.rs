use std::env;
use std::fs::File;
use std::io::{BufWriter, Result, Write};
use std::path::Path;

include!("tables.in");

fn main() -> Result<()> {
    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("tables_gen.rs");
    let mut f = File::create(&dest_path).map(BufWriter::new)?;

    generate_dna_ascii_to_iupac(&mut f)?;
    generate_4096_lut("CODON_IUPAC_TO_AA1_LUT", DNA_TO_AA1_LUT_VEC, &mut f)?;
    generate_4096_lut("CODON_IUPAC_TO_AA1_SEC", DNA_TO_AA1_SEC_VEC, &mut f)?;
    generate_4096_lut(
        "CODON_IUPAC_TO_AA1_CHRMT_VERTEBRATE",
        DNA_TO_AA1_CHRMT_VERTEBRATE_VEC,
        &mut f,
    )?;

    generate_aa1_to_aa3_str_lookup_function(&mut f)?;
    generate_aa1_to_aa3_str_lookup_table(&mut f)?;
    generate_aa3_to_aa1_lookup_function(&mut f)?;

    f.flush()?;
    println!("cargo::rerun-if-changed=build.rs");
    println!("cargo::rerun-if-changed=tables.in");
    Ok(())
}

fn generate_dna_ascii_to_iupac(f: &mut BufWriter<File>) -> Result<()> {
    let mut result = [255u8; 256];
    let mappings = [
        (b'A', 0),
        (b'C', 1),
        (b'G', 2),
        (b'T', 3),
        (b'U', 3),
        (b'R', 4),
        (b'Y', 5),
        (b'S', 6),
        (b'W', 7),
        (b'K', 8),
        (b'M', 9),
        (b'B', 10),
        (b'D', 11),
        (b'H', 12),
        (b'V', 13),
        (b'N', 14),
    ];

    for &(c, val) in mappings.iter() {
        result[c as usize] = val;
        // Handle lowercase implicitly
        result[(c as char).to_ascii_lowercase() as usize] = val;
    }

    writeln!(
        f,
        "/// Maps ASCII to 4-bit IUPAC ID. Invalid characters map to 255."
    )?;
    write!(f, "const DNA_ASCII_TO_IUPAC: [u8; 256] = [")?;
    for v in result {
        write!(f, "{}, ", v)?;
    }
    writeln!(f, "];")?;
    Ok(())
}

fn generate_4096_lut(name: &str, vec: &[(&str, &str)], f: &mut BufWriter<File>) -> Result<()> {
    let mut result = [0u8; 4096];

    for i in 0..15 {
        for j in 0..15 {
            for k in 0..15 {
                result[(i << 8) | (j << 4) | k] = b'X';
            }
        }
    }

    let get_nibble = |c: char| -> usize {
        let mappings = [
            ('A', 0),
            ('C', 1),
            ('G', 2),
            ('T', 3),
            ('U', 3),
            ('R', 4),
            ('Y', 5),
            ('S', 6),
            ('W', 7),
            ('K', 8),
            ('M', 9),
            ('B', 10),
            ('D', 11),
            ('H', 12),
            ('V', 13),
            ('N', 14),
        ];
        mappings.iter().find(|(m, _)| *m == c).unwrap().1
    };

    for &(dna, aa) in vec {
        let chars: Vec<char> = dna.chars().collect();
        let idx = (get_nibble(chars[0]) << 8) | (get_nibble(chars[1]) << 4) | get_nibble(chars[2]);
        result[idx] = aa.as_bytes()[0];
    }

    write!(f, "const {}: [u8; 4096] = [", name)?;
    for v in result {
        write!(f, "{}, ", v)?;
    }
    writeln!(f, "];")?;
    Ok(())
}

fn generate_aa1_to_aa3_str_lookup_function(f: &mut BufWriter<File>) -> Result<()> {
    writeln!(
        f,
        "const fn _aa1_to_aa3_str(aa1: u8) -> Option<&'static str> {{"
    )?;
    writeln!(f, "    match aa1 {{")?;
    for (aa3, aa1) in AA3_TO_AA1_VEC {
        writeln!(f, "        b'{}' => Some(\"{}\"),", aa1, aa3)?;
    }
    writeln!(f, r"        _ => None,")?;
    writeln!(f, "    }}")?;
    writeln!(f, "}}")?;
    Ok(())
}

fn generate_aa1_to_aa3_str_lookup_table(f: &mut BufWriter<File>) -> Result<()> {
    let mut result = [""; 256];
    for (aa3, aa1) in AA3_TO_AA1_VEC {
        result[aa1.as_bytes()[0] as usize] = aa3;
    }
    write!(f, "const AA1_TO_AA3_STR: [Option<&str>; 256] = [")?;
    for v in result {
        if v.is_empty() {
            write!(f, "None, ")?;
        } else {
            write!(f, r##"Some("{}"), "##, v)?;
        }
    }
    writeln!(f, "];")?;
    Ok(())
}

fn generate_aa3_to_aa1_lookup_function(f: &mut BufWriter<File>) -> Result<()> {
    writeln!(f, "const fn _aa3_to_aa1(aa3: &[u8]) -> Option<u8> {{")?;
    writeln!(f, "    match aa3 {{")?;
    for (aa3, aa1) in AA3_TO_AA1_VEC {
        writeln!(f, "        b\"{}\" => Some(b'{}'),", aa3, aa1)?;
    }
    writeln!(f, "        _ => None,")?;
    writeln!(f, "    }}")?;
    writeln!(f, "}}")?;
    Ok(())
}
