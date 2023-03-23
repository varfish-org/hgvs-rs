//! Provide implementation of validation to data structures.

use std::ops::Range;

use crate::validator::Validateable;

use super::{
    CdsInterval, CdsLocEdit, GenomeInterval, GenomeLocEdit, HgvsVariant, MtLocEdit, NaEdit,
    ProtLocEdit, RnaLocEdit, TxLocEdit,
};

impl Validateable for NaEdit {
    fn validate(&self) -> Result<(), anyhow::Error> {
        match &self {
            NaEdit::RefAlt {
                reference,
                alternative,
            } => {
                if reference.is_empty() && alternative.is_empty() {
                    Err(anyhow::anyhow!("Reference or alternative must be non-empty; use DelRef or Ins otherwise ({}/{})", reference, alternative))
                } else {
                    Ok(())
                }
            }
            NaEdit::NumAlt { count, alternative } => {
                if *count < 1 {
                    Err(anyhow::anyhow!(
                        "Number of deleted bases must be >=1 but was {}",
                        count
                    ))
                } else if alternative.is_empty() {
                    Err(anyhow::anyhow!(
                        "Number of alternative bases must non-empty",
                    ))
                } else {
                    Ok(())
                }
            }
            NaEdit::DelRef { reference: _ } => Ok(()),
            NaEdit::DelNum { count } => {
                if *count < 1 {
                    Err(anyhow::anyhow!(
                        "Number of deleted bases must be >=1 but was {}",
                        count
                    ))
                } else {
                    Ok(())
                }
            }
            NaEdit::Ins { alternative: _ } => Ok(()),
            NaEdit::Dup { reference: _ } => Ok(()),
            NaEdit::InvRef { reference: _ } => Ok(()),
            NaEdit::InvNum { count } => {
                if *count < 1 {
                    Err(anyhow::anyhow!(
                        "Number of inverted bases must be >=1 but was {}",
                        count
                    ))
                } else {
                    Ok(())
                }
            }
        }
    }
}

impl Validateable for HgvsVariant {
    fn validate(&self) -> Result<(), anyhow::Error> {
        // NB: we only need to validate `self.loc_edit`.  The cases that the Python library
        // considers are fended off by the Rust type system.
        match &self {
            HgvsVariant::CdsVariant { loc_edit, .. } => loc_edit.validate(),
            HgvsVariant::GenomeVariant { loc_edit, .. } => loc_edit.validate(),
            HgvsVariant::MtVariant { loc_edit, .. } => loc_edit.validate(),
            HgvsVariant::TxVariant { loc_edit, .. } => loc_edit.validate(),
            HgvsVariant::ProtVariant { loc_edit, .. } => loc_edit.validate(),
            HgvsVariant::RnaVariant { loc_edit, .. } => loc_edit.validate(),
        }
    }
}

impl Validateable for CdsLocEdit {
    fn validate(&self) -> Result<(), anyhow::Error> {
        let loc = self.loc.inner();
        loc.validate()?;

        let maybe_range: Result<Range<i32>, _> = loc.clone().try_into();
        let range = if let Ok(range) = maybe_range {
            range
        } else {
            log::trace!(
                "Skipping CDS location because loc cannot be converted to range: {:?}",
                loc
            );
            return Ok(());
        };

        match self.edit.inner() {
            NaEdit::RefAlt { .. }
            | NaEdit::DelRef { .. }
            | NaEdit::Dup { .. }
            | NaEdit::InvRef { .. } => {
                // We cannot make assumptions about reference length as we can have positon
                // offsets.
                Ok(())
            }
            NaEdit::Ins { .. } => {
                if range.len() != 2 {
                    Err(anyhow::anyhow!(
                        "Locus length must be 1 for insertions ({}, {:?})",
                        &self,
                        &range
                    ))
                } else {
                    Ok(())
                }
            }
            NaEdit::DelNum { count } | NaEdit::NumAlt { count, .. } | NaEdit::InvNum { count } => {
                if range.len() as i32 != *count {
                    Err(anyhow::anyhow!(
                        "Length implied by coordinates must equal count ({})",
                        &self
                    ))
                } else {
                    Ok(())
                }
            }
        }
    }
}

impl Validateable for CdsInterval {
    fn validate(&self) -> Result<(), anyhow::Error> {
        Ok(()) // TODO
    }
}

impl Validateable for GenomeLocEdit {
    fn validate(&self) -> Result<(), anyhow::Error> {
        self.loc.inner().validate()?;
        self.edit.inner().validate()
    }
}

impl Validateable for GenomeInterval {
    fn validate(&self) -> Result<(), anyhow::Error> {
        if let Some(start) = self.start {
            if start < 1 {
                return Err(anyhow::anyhow!("Start must be >= 1 but was {}", start));
            }
        }
        if let Some(end) = self.end {
            if end < 1 {
                return Err(anyhow::anyhow!("End must be >= 1 but was {}", end));
            }
        }
        if let (Some(start), Some(end)) = (self.start, self.end) {
            if start > end {
                return Err(anyhow::anyhow!(
                    "Start must be <= end position but was {} > {}",
                    start,
                    end
                ));
            }
        }

        Ok(())
    }
}

impl Validateable for MtLocEdit {
    fn validate(&self) -> Result<(), anyhow::Error> {
        Ok(()) // TODO
    }
}

impl Validateable for TxLocEdit {
    fn validate(&self) -> Result<(), anyhow::Error> {
        Ok(()) // TODO
    }
}

impl Validateable for RnaLocEdit {
    fn validate(&self) -> Result<(), anyhow::Error> {
        Ok(()) // TODO
    }
}

impl Validateable for ProtLocEdit {
    fn validate(&self) -> Result<(), anyhow::Error> {
        Ok(()) // TODO
    }
}

#[cfg(test)]
mod test {
    use crate::{parser::GenomeInterval, validator::Validateable};

    #[test]
    fn validate_genomeinterval() -> Result<(), anyhow::Error> {
        let g_interval = GenomeInterval {
            start: Some(1),
            end: Some(2),
        };
        assert!(g_interval.validate().is_ok());

        let g_interval = GenomeInterval {
            start: Some(-1),
            end: Some(2),
        };
        assert!(g_interval.validate().is_err());

        let g_interval = GenomeInterval {
            start: Some(1),
            end: Some(-2),
        };
        assert!(g_interval.validate().is_err());

        let g_interval = GenomeInterval {
            start: Some(2),
            end: Some(1),
        };
        assert!(g_interval.validate().is_err());

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
