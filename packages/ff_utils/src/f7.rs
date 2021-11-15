use crate::ff_utils::{FromBytes, ScalarOps, ToBytes};
use ff::{Field, PrimeField, PrimeFieldRepr};
// use bellman::PrimeField;
use core::iter::FromIterator;
use hex::{FromHexError, ToHex};
use num::bigint::BigUint;
use std::ops::{Add, Mul, Sub};

#[derive(PrimeField)]
#[PrimeFieldModulus = "7"]
#[PrimeFieldGenerator = "3"]
pub struct F7(U8);

impl ToHex for F7 {
  // Parse a Fp value to a hex string with 0x-prefix.
  fn encode_hex<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.into_repr());
    T::from_iter(repr[2..].chars())
  }

  fn encode_hex_upper<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.into_repr());
    T::from_iter(format!("{}", &repr.to_uppercase()[2..]).chars())
  }
}

#[test]
fn test_f7_to_hex() {
  let input = 1;
  let x = F7::from_repr(U8::from(input)).unwrap();
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
  let x = F7::from_repr(U8::from(input)).unwrap();
  let x_bytes_be = x.to_bytes_be().unwrap();
  let x_bytes_le = x.to_bytes_le().unwrap();
  assert_eq!(x_bytes_be, vec![0, 0, 0, 0, 0, 0, 0, 1]);
  assert_eq!(x_bytes_be.len(), 8);
  assert_eq!(x_bytes_le, vec![1, 0, 0, 0, 0, 0, 0, 0]);
  assert_eq!(x_bytes_le.len(), 8);
}

// override from_str method because F7::from_repr(U8::from(10)) cause panic.
impl F7 {
  fn from_str(s: &str) -> Option<Self> {
    if s.is_empty() {
      return None;
    }

    if s == "0" {
      return Some(Self::zero());
    }

    let mut res = Self::zero();

    let mut first_digit = true;

    for c in s.chars() {
      match c.to_digit(10) {
        Some(c) => {
          if first_digit {
            if c == 0 {
              return None;
            }

            first_digit = false;
          } else {
            res.mul_assign(&Self::from_repr(U8::from(10)).unwrap());
          }

          res.add_assign(&Self::from_repr(U8::from(c as u64)).unwrap());
        }
        None => {
          return None;
        }
      }
    }

    Some(res)
  }
}

impl FromBytes for F7 {
  fn from_bytes_be(value: &[u8]) -> Option<Self> {
    Self::from_str(&BigUint::from_bytes_be(value.as_ref()).to_str_radix(10))
  }
  fn from_bytes_le(value: &[u8]) -> Option<Self> {
    Self::from_str(&BigUint::from_bytes_le(value.as_ref()).to_str_radix(10))
  }
}

#[test]
fn test_fp_from_bytes() {
  let input = vec![6];
  let answer = 6;
  let answer_fp = F7::from_repr(U8::from(answer)).unwrap();
  let input_bytes = F7::from_bytes_be(&input).unwrap();
  println!("input_bytes: {:?}", input_bytes);
  assert_eq!(input_bytes, answer_fp);
}

#[test]
#[should_panic]
fn test_panic_fp_from_bytes() {
  let input = vec![31];
  let answer = 31;
  let answer_fp = F7::from_repr(U8::from(answer)).unwrap();
  let input_bytes = F7::from_bytes_be(&input).unwrap();
  assert_eq!(input_bytes, answer_fp);
}

impl Add<Self> for F7 {
  type Output = Self;

  fn add(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.add_assign(&rhs);
    a
  }
}

impl<'r> Add<&'r Self> for F7 {
  type Output = Self;

  fn add(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.add_assign(rhs);
    a
  }
}

impl Sub<Self> for F7 {
  type Output = Self;

  fn sub(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.sub_assign(&rhs);
    a
  }
}

impl<'r> Sub<&'r Self> for F7 {
  type Output = Self;

  fn sub(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.sub_assign(rhs);
    a
  }
}

impl Mul<Self> for F7 {
  type Output = Self;

  fn mul(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.mul_assign(&rhs);
    a
  }
}

impl<'r> Mul<&'r Self> for F7 {
  type Output = Self;

  fn mul(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.mul_assign(rhs);
    a
  }
}

impl ScalarOps for F7 {}
