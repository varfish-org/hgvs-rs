//! This module contains the code for HGVS variant descriptions.
//!
//! The parsing functionality is provided through `Type::parse()` functions.
//! The data structures also provide the `Display` trait for conversion to
//! strings etc.

mod display;
mod ds;
mod impl_parse;
mod impl_validate;
mod parse_funcs;

use std::str::FromStr;

pub use crate::parser::display::*;
pub use crate::parser::ds::*;
use crate::parser::impl_parse::*;

impl FromStr for HgvsVariant {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|e| anyhow::anyhow!("problem parsing HGVS string {:?}", &e))
            .map(|(_rest, variant)| variant)
    }
}

impl FromStr for GenomeInterval {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|e| anyhow::anyhow!("problem with parsing HGVS genome interval: {:?}", &e))
            .map(|(_rest, g_interval)| g_interval)
    }
}

impl FromStr for TxInterval {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|e| anyhow::anyhow!("problem with parsing HGVS transcript interval: {:?}", &e))
            .map(|(_rest, g_interval)| g_interval)
    }
}

impl FromStr for CdsInterval {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|e| anyhow::anyhow!("problem with parsing HGVS CDS interval: {:?}", &e))
            .map(|(_rest, g_interval)| g_interval)
    }
}

#[cfg(test)]
mod test {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
        str::FromStr,
    };

    use crate::parser::{
        Accession, CdsFrom, CdsInterval, CdsLocEdit, CdsPos, GenomeInterval, Mu, NaEdit,
    };

    use super::HgvsVariant;

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
                        start: CdsPos {
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
    fn not_ok() -> Result<(), anyhow::Error> {
        assert!(HgvsVariant::from_str("x").is_err());

        Ok(())
    }

    // This test uses the "gauntlet" file from the hgvs package.
    #[test]
    fn hgvs_gauntlet() -> Result<(), anyhow::Error> {
        let reader = BufReader::new(File::open("tests/data/parser/gauntlet")?);

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if !line.starts_with('#') && !line.is_empty() {
                let result = HgvsVariant::from_str(line);
                assert!(result.is_ok(), "line = {}; result = {:?}", &line, &result);
            }
        }

        Ok(())
    }

    // This test uses the "reject" file from the hgvs package.
    #[test]
    fn hgvs_reject() -> Result<(), anyhow::Error> {
        let reader = BufReader::new(File::open("tests/data/parser/reject")?);

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if !line.starts_with('#') && !line.is_empty() {
                assert!(HgvsVariant::from_str(line).is_err(), "line = {line}")
            }
        }

        Ok(())
    }

    // Test genome interval parsing.
    #[test]
    fn genome_interval_from_str() -> Result<(), anyhow::Error> {
        assert!(GenomeInterval::from_str("x").is_err());
        assert_eq!(
            GenomeInterval::from_str("1")?,
            GenomeInterval {
                start: Some(1),
                end: Some(1)
            }
        );
        assert_eq!(
            GenomeInterval::from_str("1_1")?,
            GenomeInterval {
                start: Some(1),
                end: Some(1)
            }
        );
        assert_eq!(
            GenomeInterval::from_str("?_1")?,
            GenomeInterval {
                start: None,
                end: Some(1)
            }
        );
        assert_eq!(
            GenomeInterval::from_str("1_?")?,
            GenomeInterval {
                start: Some(1),
                end: None
            }
        );
        assert_eq!(
            GenomeInterval::from_str("?_?")?,
            GenomeInterval {
                start: None,
                end: None
            }
        );

        Ok(())
    }
}
