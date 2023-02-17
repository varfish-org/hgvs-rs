//! Implementation of validation.

use std::rc::Rc;

use crate::{data::interface::Provider, parser::HgvsVariant};

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum ValidationLevel {
    Null,
    Intrinsic,
    Extrinsic,
}

impl ValidationLevel {
    pub fn validator(&self, strict: bool, provider: Rc<dyn Provider>) -> Box<dyn Validator> {
        match self {
            ValidationLevel::Null => Box::new(NullValidator::new(strict)),
            ValidationLevel::Intrinsic => Box::new(IntrinsicValidator::new(strict)),
            ValidationLevel::Extrinsic => Box::new(ExtrinsicValidator::new(strict, provider)),
        }
    }
}

pub trait Validator {
    fn is_strict(&self) -> bool;

    fn validate(&mut self, var: &HgvsVariant) -> Result<(), anyhow::Error>;
}

pub struct NullValidator {
    strict: bool,
}

impl NullValidator {
    pub fn new(strict: bool) -> Self {
        Self { strict }
    }
}

impl Validator for NullValidator {
    fn is_strict(&self) -> bool {
        self.strict
    }

    fn validate(&mut self, _var: &HgvsVariant) -> Result<(), anyhow::Error> {
        todo!()
    }
}

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

    fn validate(&mut self, _var: &HgvsVariant) -> Result<(), anyhow::Error> {
        todo!()
    }
}

pub struct ExtrinsicValidator {
    strict: bool,
    provider: Rc<dyn Provider>,
}

impl ExtrinsicValidator {
    pub fn new(strict: bool, provider: Rc<dyn Provider>) -> Self {
        Self { strict, provider }
    }
}

impl Validator for ExtrinsicValidator {
    fn is_strict(&self) -> bool {
        self.strict
    }

    fn validate(&mut self, _var: &HgvsVariant) -> Result<(), anyhow::Error> {
        todo!()
    }
}
