//! This module contains the code for HGVS variant descriptions.
//!
//! The parsing functionality is provided through `Type::parse()` functions.
//! The data structures also provide the `Display` trait for conversion to
//! strings etc.

mod display;
mod ds;
mod error;
mod impl_parse;
mod impl_validate;
mod parse_funcs;

use std::str::FromStr;

pub use crate::parser::display::*;
pub use crate::parser::ds::*;
pub use crate::parser::error::*;
use crate::parser::impl_parse::*;

impl FromStr for HgvsVariant {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|_e| Error::InvalidHgvsVariant(s.to_string()))
            .map(|(_rest, variant)| variant)
    }
}

impl FromStr for GenomeInterval {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|_e| Error::InvalidGenomeInterval(s.to_string()))
            .map(|(_rest, g_interval)| g_interval)
    }
}

impl FromStr for TxInterval {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|_e| Error::InvalidTxInterval(s.to_string()))
            .map(|(_rest, g_interval)| g_interval)
    }
}

impl FromStr for CdsInterval {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
            .map_err(|_e| Error::InvalidCdsInterval(s.to_string()))
            .map(|(_rest, g_interval)| g_interval)
    }
}

#[cfg(test)]
mod test {
    use anyhow::Error;
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
    fn from_str_basic() -> Result<(), Error> {
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
                        reference: b"A".to_vec(),
                        alternative: b"T".to_vec()
                    })
                }
            }
        );

        Ok(())
    }

    #[test]
    fn not_ok() -> Result<(), Error> {
        assert!(HgvsVariant::from_str("x").is_err());

        Ok(())
    }

    // This test uses the "gauntlet" file from the hgvs package.
    #[test]
    fn hgvs_gauntlet() -> Result<(), Error> {
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
    fn hgvs_reject() -> Result<(), Error> {
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
    fn genome_interval_from_str() -> Result<(), Error> {
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
