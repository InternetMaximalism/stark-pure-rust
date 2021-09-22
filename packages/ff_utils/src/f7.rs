use crate::ff_utils::{FromBytes, ToBytes};
use core::iter::FromIterator;
use ff::PrimeField;
use hex::{FromHexError, ToHex};
use num::bigint::BigUint;

#[derive(PrimeField)]
#[PrimeFieldModulus = "7"]
#[PrimeFieldGenerator = "3"]
#[PrimeFieldReprEndianness = "little"]
pub struct F7([u64; 1]);

impl ToHex for F7 {
  // Parse a Fp value to a hex string with 0x-prefix.
  fn encode_hex<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.to_repr());
    T::from_iter(repr[2..].chars())
  }

  fn encode_hex_upper<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.to_repr());
    T::from_iter(format!("{}", &repr.to_uppercase()[2..]).chars())
  }
}

#[test]
fn test_f7_to_hex() {
  let input = 1;
  let x = F7::from(input);
  assert_eq!(x.encode_hex::<String>(), format!("{:016x}", input));
  assert_eq!(x.encode_hex_upper::<String>(), format!("{:016X}", input));
}

impl ToBytes for F7 {
  fn to_bytes_be(&self) -> Result<Vec<u8>, FromHexError> {
    hex::decode(&self.encode_hex::<String>())
  }
  fn to_bytes_le(&self) -> Result<Vec<u8>, FromHexError> {
    let mut res = self.to_bytes_be()?;
    res.reverse();
    Ok(res)
  }
}

#[test]
fn test_f7_to_bytes() {
  let input = 1;
  let x = F7::from(input);
  let x_bytes_be = x.to_bytes_be().unwrap();
  let x_bytes_le = x.to_bytes_le().unwrap();
  assert_eq!(x_bytes_be, vec![0, 0, 0, 0, 0, 0, 0, 1]);
  assert_eq!(x_bytes_be.len(), 8);
  assert_eq!(x_bytes_le, vec![1, 0, 0, 0, 0, 0, 0, 0]);
  assert_eq!(x_bytes_le.len(), 8);
}

impl FromBytes for F7 {
  fn from_bytes_be(value: &[u8]) -> Option<Self> {
    Self::from_str(&BigUint::from_bytes_be(value.as_ref()).to_str_radix(10))
  }
  fn from_bytes_le(value: &[u8]) -> Option<Self> {
    Self::from_str(&BigUint::from_bytes_le(value.as_ref()).to_str_radix(10))
  }
}
