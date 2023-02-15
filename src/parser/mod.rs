//! This module contains the code for HGVS variant descriptions.
//!
//! The parsing functionality is provided through `Type::parse()` functions.
//! The data structures also provide the `Display` trait for conversion to
//! strings etc.

mod ds;
mod impl_parse;
mod parse_funcs;

use std::str::FromStr;

pub use crate::parser::ds::*;
pub use crate::parser::impl_parse::*;

impl FromStr for HgvsVariant {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|e| anyhow::anyhow!("problem parsing HGVS string {:?}", &e))
            .map(|(_rest, variant)| variant)
    }
}

#[cfg(test)]
mod test {
    use std::{str::FromStr, fs::File, io::{BufReader, BufRead}};

    use crate::parser::{Accession, CdsFrom, CdsInterval, CdsLocEdit, CdsPos, Mu, NaEdit};

    use super::{HgvsVariant, Parseable};

    #[test]
    fn from_str_basic() -> Result<(), anyhow::Error> {
        assert_eq!(
            HgvsVariant::from_str("NM_01234.5:c.22+1A>T")?,
            HgvsVariant::CdsVariant {
                accession: Accession {
                    value: "NM_01234.5".to_string()
                },
                gene_symbol: None,
                loc_edit: CdsLocEdit {
                    loc: Mu::Certain(CdsInterval {
                        begin: CdsPos {
                            base: 22,
                            offset: Some(1),
                            cds_from: CdsFrom::Start
                        },
                        end: CdsPos {
                            base: 22,
                            offset: Some(1),
                            cds_from: CdsFrom::Start
                        }
                    }),
                    edit: Mu::Certain(NaEdit::RefAlt {
                        reference: "A".to_string(),
                        alternative: "T".to_string()
                    })
                }
            }
        );

        Ok(())
    }

    #[test]
    fn x() -> Result<(), anyhow::Error> {
        // HgvsVariant::from_str("AC_01234.5:c.76_78del")?;
        // CdsLocEdit::parse("76_78del")?;
        CdsInterval::parse("76_78")?;
        Ok(())
    }

    // This test uses the "gauntlet" file from the hgvs package.
    #[test]
    fn hgvs_gauntlet() -> Result<(), anyhow::Error> {
        let reader = BufReader::new(File::open("tests/data/gauntlet")?);

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if !line.starts_with("#") && !line.is_empty() {
                assert!(
                    HgvsVariant::from_str(line).is_ok(),
                    "line = {}", line
                )
            }
        }

        Ok(())
    }

    // This test uses the "reject" file from the hgvs package.
    #[test]
    fn hgvs_reject() -> Result<(), anyhow::Error> {
        let reader = BufReader::new(File::open("tests/data/reject")?);

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if !line.starts_with("#") && !line.is_empty() {
                assert!(
                    !HgvsVariant::from_str(line).is_ok(),
                    "line = {}", line
                )
            }
        }

        Ok(())
    }
}
