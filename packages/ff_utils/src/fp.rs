use crate::ff_utils::{FromBytes, ScalarOps, ToBytes};
use core::iter::FromIterator;
use ff::{Field, PrimeField, PrimeFieldRepr};
use hex::{FromHexError, ToHex};
use num::bigint::BigUint;
use std::ops::{Add, Mul, Sub};

// #[PrimeFieldModulus = "115792089237316195423570985008687907853269984665640564039457584006405596119041"]
#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
pub struct Fp(U256);

impl ToHex for Fp {
  // Parse a Fp value to a hex string with 0x-prefix.
  fn encode_hex<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.into_repr());
    T::from_iter(repr[2..].chars())
  }

  fn encode_hex_upper<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.into_repr());
    T::from_iter(repr.to_uppercase()[2..].chars())
  }
}

#[test]
fn test_fp_to_hex() {
  let input = 31;
  let x = Fp::from_repr(U256::from(input)).unwrap();
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
  let x = Fp::from_repr(U256::from(input)).unwrap();
  let x_bytes_be = x.to_bytes_be().unwrap();
  let x_bytes_le = x.to_bytes_le().unwrap();
  assert_eq!(
    x_bytes_be,
    vec![
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      31
    ]
  );
  assert_eq!(x_bytes_be.len(), 32);
  assert_eq!(
    x_bytes_le,
    vec![
      31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0
    ]
  );
  assert_eq!(x_bytes_le.len(), 32);
}

impl FromBytes for Fp {
  fn from_bytes_be(value: &[u8]) -> Option<Self> {
    Self::from_str(&BigUint::from_bytes_be(value.as_ref()).to_str_radix(10))
  }
  fn from_bytes_le(value: &[u8]) -> Option<Self> {
    Self::from_str(&BigUint::from_bytes_le(value.as_ref()).to_str_radix(10))
  }
}

#[test]
fn test_fp_from_bytes() {
  let input_be = vec![
    0u8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    31,
  ];
  let input_le = vec![
    31u8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,
  ];
  let answer = 31;
  let answer_fp = Fp::from_repr(U256::from(answer)).unwrap();
  let input_bytes_be = Fp::from_bytes_be(&input_be).unwrap();
  let input_bytes_le = Fp::from_bytes_le(&input_le).unwrap();
  assert_eq!(input_bytes_be, answer_fp);
  assert_eq!(input_bytes_le, answer_fp);
}

#[test]
fn test_fp_from_bytes_big_number() {
  let input_be = [
    1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 0,
  ];
  let _ = Fp::from_bytes_be(&input_be).unwrap();

  let input_le = [
    1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 0,
  ];
  let _ = Fp::from_bytes_le(&input_le).unwrap();
}

impl Add<Self> for Fp {
  type Output = Self;

  fn add(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.add_assign(&rhs);
    a
  }
}

impl<'r> Add<&'r Self> for Fp {
  type Output = Self;

  fn add(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.add_assign(rhs);
    a
  }
}

impl Sub<Self> for Fp {
  type Output = Self;

  fn sub(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.sub_assign(&rhs);
    a
  }
}

impl<'r> Sub<&'r Self> for Fp {
  type Output = Self;

  fn sub(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.sub_assign(rhs);
    a
  }
}

impl Mul<Self> for Fp {
  type Output = Self;

  fn mul(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.mul_assign(&rhs);
    a
  }
}

impl<'r> Mul<&'r Self> for Fp {
  type Output = Self;

  fn mul(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.mul_assign(rhs);
    a
  }
}

impl ScalarOps for Fp {}
