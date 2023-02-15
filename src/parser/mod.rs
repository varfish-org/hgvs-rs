//! This module contains the code for HGVS variant descriptions.
//!
//! The parsing functionality is provided through `Type::parse()` functions.
//! The data structures also provide the `Display` trait for conversion to
//! strings etc.

mod ds;
mod impl_parse;
mod parse_funcs;

pub use crate::parser::ds::*;
pub use crate::parser::impl_parse::*;
