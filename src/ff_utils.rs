use hex::FromHexError;

pub trait ToBytes {
  fn to_bytes_be(&self) -> Result<Vec<u8>, FromHexError>;
  fn to_bytes_le(&self) -> Result<Vec<u8>, FromHexError>;
}
