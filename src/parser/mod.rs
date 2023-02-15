//! This module contains the code for HGVS variant descriptions.
//!
//! The parsing functionality is provided through `Type::parse()` functions.
//! The data structures also provide the `Display` trait for conversion to
//! strings etc.

mod ds;
mod impls_parse;
mod funcs;

pub use crate::parser::ds::*;
pub use crate::parser::impls_parse::*;
