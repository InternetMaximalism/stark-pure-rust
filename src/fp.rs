use core::iter::FromIterator;
use hex::{ToHex, FromHexError};
use ff::PrimeField;
use crate::ff_utils::ToBytes;

#[derive(PrimeField)]
#[PrimeFieldModulus = "115792089237316195423570985008687907853269984665640564039457584006405596119041"]
#[PrimeFieldGenerator = "7"]
#[PrimeFieldReprEndianness = "little"]
pub struct Fp([u64; 5]);

impl ToHex for Fp {
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
fn test_fp_to_hex() {
  let input = 31;
  let x = Fp::from(input);
  assert_eq!(x.encode_hex::<String>(), format!("{:064x}", input));
  assert_eq!(x.encode_hex_upper::<String>(), format!("{:064X}", input));
}

impl ToBytes for Fp {
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
fn test_fp_to_bytes() {
  let input = 31;
  let x = Fp::from(input);
  let x_bytes_be = x.to_bytes_be().unwrap();
  let x_bytes_le = x.to_bytes_le().unwrap();
  assert_eq!(x_bytes_be, vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 31]);
  assert_eq!(x_bytes_be.len(), 32);
  assert_eq!(x_bytes_le, vec![31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
  assert_eq!(x_bytes_le.len(), 32);
}
