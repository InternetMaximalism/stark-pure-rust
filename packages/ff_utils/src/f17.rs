use crate::ff_utils::{FromBytes, ScalarOps, ToBytes};
use ff::{Field, PrimeField, PrimeFieldRepr};
// use bellman::PrimeField;
use core::iter::FromIterator;
use hex::{FromHexError, ToHex};
use num::bigint::BigUint;
use std::ops::{Add, Mul, Sub};

#[derive(PrimeField)]
#[PrimeFieldModulus = "17"]
#[PrimeFieldGenerator = "3"]
pub struct F17(U8);

impl ToHex for F17 {
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
  let x = F17::from_repr(U8::from(input)).unwrap();
  assert_eq!(x.encode_hex::<String>(), format!("{:016x}", input));
  assert_eq!(x.encode_hex_upper::<String>(), format!("{:016X}", input));
}

impl ToBytes for F17 {
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
  let x = F17::from_repr(U8::from(input)).unwrap();
  let x_bytes_be = x.to_bytes_be().unwrap();
  let x_bytes_le = x.to_bytes_le().unwrap();
  assert_eq!(x_bytes_be, vec![0, 0, 0, 0, 0, 0, 0, 1]);
  assert_eq!(x_bytes_be.len(), 8);
  assert_eq!(x_bytes_le, vec![1, 0, 0, 0, 0, 0, 0, 0]);
  assert_eq!(x_bytes_le.len(), 8);
}

impl FromBytes for F17 {
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
  let answer_fp = F17::from_repr(U8::from(answer)).unwrap();
  let input_bytes = F17::from_bytes_be(&input).unwrap();
  println!("input_bytes: {:?}", input_bytes);
  assert_eq!(input_bytes, answer_fp);
}

#[test]
#[should_panic]
fn test_panic_fp_from_bytes() {
  let input = vec![31];
  let answer = 31;
  let answer_fp = F17::from_repr(U8::from(answer)).unwrap();
  let input_bytes = F17::from_bytes_be(&input).unwrap();
  assert_eq!(input_bytes, answer_fp);
}

impl Add<Self> for F17 {
  type Output = Self;

  fn add(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.add_assign(&rhs);
    a
  }
}

impl<'r> Add<&'r Self> for F17 {
  type Output = Self;

  fn add(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.add_assign(rhs);
    a
  }
}

impl Sub<Self> for F17 {
  type Output = Self;

  fn sub(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.sub_assign(&rhs);
    a
  }
}

impl<'r> Sub<&'r Self> for F17 {
  type Output = Self;

  fn sub(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.sub_assign(rhs);
    a
  }
}

impl Mul<Self> for F17 {
  type Output = Self;

  fn mul(self, rhs: Self) -> Self::Output {
    let mut a = self;
    a.mul_assign(&rhs);
    a
  }
}

impl<'r> Mul<&'r Self> for F17 {
  type Output = Self;

  fn mul(self, rhs: &Self) -> Self::Output {
    let mut a = self;
    a.mul_assign(rhs);
    a
  }
}

impl ScalarOps for F17 {}

// #[test]
// fn test_f17() {
//   use bellman::plonk::polynomials::Polynomial;
//   use bellman::worker::Worker;
//   // let values = [0u8, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7]
//   let values = [12u8, 9, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
//     .iter()
//     .map(|v| F17::from_bytes_le(&[*v]).unwrap())
//     .collect::<Vec<_>>();
//   // let worker = Worker::new();
//   // let coeffs = Polynomial::from_values(values)
//   //   .unwrap()
//   //   .ifft(&worker)
//   //   .into_coeffs();
//   let powers = expand_root_of_unity(F17::multiplicative_generator());
//   println!("powers: {:?}", powers);

//   let coeffs = simple_ft(&values, &powers);
//   println!("coeffs: {:?}", coeffs);
//   let dividend = [8u8, 1]
//     .iter()
//     .map(|v| F17::from_bytes_le(&[*v]).unwrap())
//     .collect::<Vec<_>>();
//   let quotient = div_polys(&coeffs, &dividend);
//   println!("quotient: {:?}", quotient);
//   let coeffs2 = mul_polys(&quotient, &dividend);
//   println!("coeffs2: {:?}", coeffs2);
//   assert_eq!(coeffs, coeffs2);
// }

// fn simple_ft<T: PrimeField + ScalarOps>(values: &[T], roots_of_unity: &[T]) -> Vec<T> {
//   let m = roots_of_unity.len();
//   let n = values.len();
//   let mut values = values.to_vec();
//   if m > n {
//     for _ in 0..(m - n) {
//       values.push(T::zero());
//     }
//   }

//   let mut o = vec![T::zero(); m];
//   for i in 0..m {
//     for j in 0..m {
//       let v = values[j] * roots_of_unity[(i * j) % m];
//       o[i].add_assign(&v);
//     }
//   }

//   o
// }

// pub fn expand_root_of_unity<T: PrimeField + ScalarOps>(root_of_unity: T) -> Vec<T> {
//   let mut output = vec![T::one()];
//   let mut current_root = root_of_unity;
//   while current_root != T::one() {
//     output.push(current_root);
//     current_root.mul_assign(&root_of_unity);
//   }

//   output
// }

// pub fn mul_polys<T: PrimeField + ScalarOps>(a: &[T], b: &[T]) -> Vec<T> {
//   let mut o = vec![T::zero(); a.len() + b.len() - 1];
//   for i in 0..a.len() {
//     for j in 0..b.len() {
//       o[i + j] = o[i + j] + a[i] * b[j];
//     }
//   }

//   o
// }

// pub fn div_polys<T: PrimeField + ScalarOps>(a: &[T], b: &[T]) -> Vec<T> {
//   let zeros_len = b
//     .iter()
//     .rev()
//     .take_while(|&&coeff| coeff == T::zero())
//     .collect::<Vec<&T>>()
//     .len();
//   let b = b[..(b.len() - zeros_len)].to_vec();

//   assert!(a.len() >= b.len());
//   let mut c = a.to_vec();
//   let mut o = vec![];
//   let mut apos = a.len() - 1;
//   let bpos = b.len() - 1;
//   let diff = apos - bpos;
//   for d in (0..(diff + 1)).rev() {
//     let quot = c[apos] * b[bpos].inverse().unwrap();
//     o.push(quot);
//     for i in (0..(bpos + 1)).rev() {
//       let v = b[i] * quot;
//       c[d + i].sub_assign(&v);
//     }

//     apos -= 1;
//   }

//   o.reverse();
//   o
// }
