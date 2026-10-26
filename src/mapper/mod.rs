//! Code supporting mapping between coordinate systems.

pub mod alignment;
pub mod altseq;
pub mod assembly;
pub mod cigar;
mod error;
pub mod variant;

pub use error::Error;

use std::sync::{Arc, Weak};

use crate::data::interface::Provider;

/// The identity of `provider` in the cache keys of the mappers: its address.
///
/// Each cache entry holds a `Weak` reference to its provider, see `CacheEntry`.  This keeps the
/// address allocated, so that no other provider can have it while the entry exists.
pub(crate) fn provider_id(provider: &Arc<dyn Provider + Send + Sync>) -> usize {
    Arc::as_ptr(provider).cast::<()>().addr()
}

/// A cached value of the provider whose `provider_id()` is in the key of the entry.
pub(crate) struct CacheEntry<T> {
    _provider: Weak<dyn Provider + Send + Sync>,
    pub value: T,
}

impl<T: Clone> CacheEntry<T> {
    pub fn new(provider: &Arc<dyn Provider + Send + Sync>, value: &T) -> Self {
        Self {
            _provider: Arc::downgrade(provider),
            value: value.clone(),
        }
    }
}
