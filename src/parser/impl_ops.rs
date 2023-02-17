//! Implementation of operations on the data structures.

use crate::data::interface::Provider;

use super::HgvsVariant;

impl HgvsVariant {
    /// Fill reference bases using the given data provider.
    ///
    /// # Args
    ///
    /// * `provider` -- The `Provider` to use for fetching reference bases.
    pub fn fill_ref(&self, _provider: &dyn Provider) -> Self {
        todo!()
    }
}
