use hex::FromHexError;
use std::ops::{Add, Mul, Sub};

pub trait ScalarOps:
  Sized
  + Add<Self, Output = Self>
  + Sub<Self, Output = Self>
  + Mul<Self, Output = Self>
  + for<'r> Add<&'r Self, Output = Self>
  + for<'r> Sub<&'r Self, Output = Self>
  + for<'r> Mul<&'r Self, Output = Self>
{
}

pub trait ToBytes {
  fn to_bytes_be(&self) -> Result<Vec<u8>, FromHexError>;
  fn to_bytes_le(&self) -> Result<Vec<u8>, FromHexError>;
}

pub trait FromBytes
where
  Self: Sized,
{
  fn from_bytes_be(value: &[u8]) -> Option<Self>;
  fn from_bytes_le(value: &[u8]) -> Option<Self>;
}
