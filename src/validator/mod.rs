//! Implementation of validation.

mod error;

use std::sync::Arc;

use log::{error, warn};

pub use crate::validator::error::Error;
use crate::{
    data::interface::Provider,
    mapper::{variant::Config, variant::Mapper},
    parser::HgvsVariant,
};

/// Trait for validating of variants, locations etc.
pub trait Validateable {
    fn validate(&self) -> Result<(), Error>;
}

/// Validation level specification.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum ValidationLevel {
    /// No validation.
    Null,
    /// Only inspect the variant description itself.
    Intrinsic,
    /// Full validation including checks based on sequence and intrinsics.
    Full,
}

impl ValidationLevel {
    pub fn validator(
        &self,
        strict: bool,
        provider: Arc<dyn Provider + Send + Sync>,
    ) -> Arc<dyn Validator + Send + Sync> {
        match self {
            ValidationLevel::Null => Arc::new(NullValidator::new()),
            ValidationLevel::Intrinsic => Arc::new(IntrinsicValidator::new(strict)),
            ValidationLevel::Full => Arc::new(FullValidator::new(strict, provider)),
        }
    }
}

/// Trait for validators.
pub trait Validator {
    /// Return whether validation is strict.
    ///
    /// Validation is strict if errors cause `Err` results rather than just logging a warning.
    fn is_strict(&self) -> bool;

    /// Validate the given variant.
    ///
    /// Depending on the configuration and implementation of the validator, an `Err` will be
    /// returned or only a warning will be logged.
    fn validate(&self, var: &HgvsVariant) -> Result<(), Error>;
}

/// A validator that performs no validation.
pub struct NullValidator {}

impl NullValidator {
    pub fn new() -> Self {
        Self {}
    }
}

impl Default for NullValidator {
    fn default() -> Self {
        Self::new()
    }
}

impl Validator for NullValidator {
    fn is_strict(&self) -> bool {
        false
    }

    fn validate(&self, _var: &HgvsVariant) -> Result<(), Error> {
        Ok(())
    }
}

/// A validator that only performs intrinsic validation.
///
/// This means that only the variant description itself is checked without considering the
/// actual sequence.
pub struct IntrinsicValidator {
    strict: bool,
}

impl IntrinsicValidator {
    pub fn new(strict: bool) -> Self {
        Self { strict }
    }
}

impl Validator for IntrinsicValidator {
    fn is_strict(&self) -> bool {
        self.strict
    }

    fn validate(&self, var: &HgvsVariant) -> Result<(), Error> {
        let res = var.validate();
        match (&res, self.is_strict()) {
            (Ok(_), _) => Ok(()),
            (Err(_), false) => {
                warn!("Validation of {} failed: {:?}", var, res);
                Ok(())
            }
            (Err(_), true) => {
                error!("Validation of {} failed: {:?}", var, res);
                res
            }
        }
    }
}

/// Attempts to determine if the HGVS name validates against external data sources
pub struct ExtrinsicValidator {
    strict: bool,
    #[allow(dead_code)]
    mapper: Mapper,
}

impl ExtrinsicValidator {
    pub fn new(strict: bool, provider: Arc<dyn Provider + Send + Sync>) -> Self {
        let config = Config {
            replace_reference: false,
            strict_validation: false,
            prevalidation_level: ValidationLevel::Null,
            add_gene_symbol: false,
            strict_bounds: true,
            renormalize_g: false,
            genome_seq_available: true,
        };
        Self {
            strict,
            mapper: Mapper::new(&config, provider),
        }
    }
}

impl Validator for ExtrinsicValidator {
    fn is_strict(&self) -> bool {
        self.strict
    }

    fn validate(&self, var: &HgvsVariant) -> Result<(), Error> {
        // Check transcripts bounds
        match var {
            HgvsVariant::CdsVariant { .. } | HgvsVariant::TxVariant { .. } => {
                let res = self.check_tx_bound(var);
                if res.is_err() {
                    if self.is_strict() {
                        error!("Validation of {} failed: {:?}", var, res);
                        return res;
                    } else {
                        warn!("Validation of {} failed: {:?}", var, res);
                    }
                }
            }
            _ => {}
        }

        // Check CDS bounds
        {
            let res = self.check_cds_bound(var);
            if res.is_err() {
                if self.is_strict() {
                    error!("Validation of {} failed: {:?}", var, res);
                    return res;
                } else {
                    warn!("Validation of {} failed: {:?}", var, res);
                }
            }
        }

        // Check reference.
        {
            let res = self.check_ref(var);
            if res.is_err() {
                if self.is_strict() {
                    error!("Validation of {} failed: {:?}", var, res);
                    return res;
                } else {
                    warn!("Validation of {} failed: {:?}", var, res);
                }
            }
        }

        Ok(())
    }
}

impl ExtrinsicValidator {
    fn check_tx_bound(&self, _var: &HgvsVariant) -> Result<(), Error> {
        Ok(()) // TODO
    }

    fn check_cds_bound(&self, _var: &HgvsVariant) -> Result<(), Error> {
        Ok(()) // TODO
    }

    fn check_ref(&self, _var: &HgvsVariant) -> Result<(), Error> {
        Ok(()) // TODO
    }
}

/// Full validator performing both intrinsic and extrinsic validation.
pub struct FullValidator {
    intrinsic: IntrinsicValidator,
    extrinsic: ExtrinsicValidator,
}

impl FullValidator {
    pub fn new(strict: bool, provider: Arc<dyn Provider + Send + Sync>) -> Self {
        Self {
            intrinsic: IntrinsicValidator::new(strict),
            extrinsic: ExtrinsicValidator::new(strict, provider),
        }
    }
}

impl Validator for FullValidator {
    fn is_strict(&self) -> bool {
        self.intrinsic.is_strict()
    }

    fn validate(&self, var: &HgvsVariant) -> Result<(), Error> {
        self.intrinsic.validate(var)?;
        self.extrinsic.validate(var)
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
