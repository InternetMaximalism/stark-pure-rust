use crate::hash::Digest;
use crate::utils::blake;

use serde::{Deserialize, Serialize};
use std::fmt::Debug;

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BlakeDigest(pub Vec<u8>);

impl Debug for BlakeDigest {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    f.write_fmt(format_args!("\"0x{}\"", hex::encode(&self.0)))
  }
}

impl AsRef<[u8]> for BlakeDigest {
  fn as_ref(&self) -> &[u8] {
    &(self.0)
  }
}

impl Default for BlakeDigest {
  fn default() -> Self {
    Self(vec![])
  }
}

impl Digest for BlakeDigest {
  fn hash(message: &[u8]) -> Self {
    Self(blake(message))
  }
}
