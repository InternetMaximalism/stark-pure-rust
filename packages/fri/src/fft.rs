use ff::PrimeField;
use std::convert::TryInto;

pub fn expand_root_of_unity<T: PrimeField>(root_of_unity: T) -> Vec<T> {
  let mut output = vec![T::one()];
  let mut current_root = root_of_unity;
  while current_root != T::one() {
    output.push(current_root);
    current_root *= root_of_unity;
  }

  output
}

#[test]
fn test_expand_root_of_unity() {
  use crate::f7::F7;

  let root_of_unity = F7::multiplicative_generator();
  let roots_of_unity: Vec<F7> = [1, 3, 2, 6, 4, 5]
    .iter()
    .map(|x| F7::from(*x as u64))
    .collect();
  let res = expand_root_of_unity(root_of_unity);
  assert_eq!(res, roots_of_unity);

  use crate::ff_utils::ToBytes;
  use crate::fp::Fp;
  use crate::utils::parse_bytes_to_u64_vec;
  use ff::Field;
  use num::bigint::BigUint;

  let precision = 65536usize;
  let times_nmr = BigUint::from_bytes_be(&(Fp::zero() - Fp::one()).to_bytes_be().unwrap());
  let times_dnm = BigUint::from_bytes_be(&precision.to_be_bytes());
  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) /precision
  let root_of_unity = Fp::multiplicative_generator().pow_vartime(&times);
  let res = expand_root_of_unity(root_of_unity);
  assert_eq!(res.len(), precision);
}

fn roots_of_unity_rev<T: PrimeField>(roots: &[T]) -> Vec<T> {
  let mut roots_rev = roots.to_vec();
  roots_rev.push(T::one());
  roots_rev.reverse();
  roots_rev.pop();
  roots_rev
}

#[test]
fn test_rev_expand_root_of_unity() {
  use crate::f7::F7;

  let roots: Vec<F7> = [1, 3, 2, 6, 4, 5]
    .iter()
    .map(|x| F7::from(*x as u64))
    .collect();
  let roots_rev = roots_of_unity_rev(&roots);
  let answer: Vec<F7> = [1, 5, 4, 6, 2, 3].iter().map(|x| F7::from(*x)).collect();
  assert_eq!(roots_rev, answer);
}

fn _simple_ft<T: PrimeField>(values: &[T], roots_of_unity: &[T]) -> Vec<T> {
  let m = roots_of_unity.len();
  let n = values.len();
  let mut values = values.to_vec();
  if m > n {
    for _ in 0..(m - n) {
      values.push(T::zero());
    }
  }

  let mut o = vec![T::zero(); m];
  for i in 0..m {
    for j in 0..m {
      o[i] += values[j] * roots_of_unity[(i * j) % m]
    }
  }

  o
}

#[test]
fn test_simple_ft() {
  use crate::f7::F7;

  let values: Vec<F7> = [1, 2, 0].iter().map(|x| F7::from(*x)).collect();
  let roots_of_unity: Vec<F7> = [1, 2, 4].iter().map(|x| F7::from(*x)).collect();
  let answer: Vec<F7> = [3, 5, 2].iter().map(|x| F7::from(*x)).collect();
  let res = _simple_ft(&values, &roots_of_unity);
  assert_eq!(res, answer);

  let values: Vec<F7> = [0, 1, 1, 0].iter().map(|x| F7::from(*x)).collect();
  let roots_of_unity: Vec<F7> = [1, 2, 4].iter().map(|x| F7::from(*x)).collect();
  let answer: Vec<F7> = [2, 6, 6].iter().map(|x| F7::from(*x)).collect();
  let res = _simple_ft(&values, &roots_of_unity);
  assert_eq!(res, answer);
}

fn _fft<T: PrimeField>(values: &[T], roots_of_unity: &[T]) -> Vec<T> {
  if values.len() <= 4 {
    return _simple_ft(values, roots_of_unity);
  }

  let m = roots_of_unity.len();
  let mut values = values.to_vec();
  if values.len() % 2 == 1 {
    values.push(T::zero());
  }
  let n = values.len();

  let mut values_from0_step2 = Vec::with_capacity(n / 2);
  let mut values_from1_step2 = Vec::with_capacity(n / 2);
  for i in 0..n {
    if i % 2 == 0 {
      values_from0_step2.push(values[i]);
    } else {
      values_from1_step2.push(values[i]);
    }
  }

  let mut roots_of_unity_step2 = Vec::with_capacity(m / 2);
  for i in 0..m {
    if i % 2 == 0 {
      roots_of_unity_step2.push(roots_of_unity[i]);
    }
  }

  let l = _fft(&values_from0_step2, &roots_of_unity_step2);
  let r = _fft(&values_from1_step2, &roots_of_unity_step2);
  let mut o = vec![T::zero(); l.len() * 2];
  for i in 0..(l.len()) {
    let x = l[i];
    let y = r[i];
    let y_times_out = y * roots_of_unity[i % m];
    o[i] = x + y_times_out;
    o[i + l.len()] = x - y_times_out;
  }

  o
}

pub fn fft<T: PrimeField>(values: &[T], root_of_unity: T) -> Vec<T> {
  let roots = expand_root_of_unity(root_of_unity);
  assert!(values.len() <= roots.len());
  _fft(values, &roots)
}

#[test]
fn test_fft() {
  use crate::f7::F7;

  let values: Vec<F7> = [1, 0, 0, 0, 0, 0].iter().map(|x| F7::from(*x)).collect();
  let res = fft(&values, F7::multiplicative_generator());
  let answer: Vec<F7> = [1, 1, 1, 1, 1, 1].iter().map(|x| F7::from(*x)).collect();
  assert_eq!(res, answer);

  let values: Vec<F7> = [1, 0, 2, 1, 0, 1].iter().map(|x| F7::from(*x)).collect();
  let res = fft(&values, F7::multiplicative_generator());
  let answer: Vec<F7> = [5, 2, 0, 1, 1, 4].iter().map(|x| F7::from(*x)).collect();
  assert_eq!(res, answer);
}

fn _inv_fft<T: PrimeField>(values: &[T], roots_rev: &[T]) -> Vec<T> {
  let m: T = T::from(roots_rev.len().try_into().unwrap());
  let inv_len = m.invert().unwrap();
  _fft(values, &roots_rev)
    .iter()
    .map(|x| *x * inv_len)
    .collect()
}

pub fn inv_fft<T: PrimeField>(values: &[T], root_of_unity: T) -> Vec<T> {
  let roots = expand_root_of_unity(root_of_unity);
  let roots_rev = roots_of_unity_rev(&roots);
  assert!(values.len() <= roots_rev.len());
  _inv_fft(values, &roots_rev)
}

#[test]
fn test_inv_fft() {
  use crate::f7::F7;

  let g2 = F7::multiplicative_generator();
  let values: Vec<F7> = [1, 1, 1, 1, 1, 1].iter().map(|x| F7::from(*x)).collect();
  let answer: Vec<F7> = [1, 0, 0, 0, 0, 0].iter().map(|x| F7::from(*x)).collect();
  let res = inv_fft(&values, g2);
  assert_eq!(res, answer);

  let values: Vec<F7> = [5, 2, 0, 1, 1, 4].iter().map(|x| F7::from(*x)).collect();
  let answer: Vec<F7> = [1, 0, 2, 1, 0, 1].iter().map(|x| F7::from(*x)).collect();
  let res = inv_fft(&values, g2);
  assert_eq!(res, answer);
}

pub fn mul_polys<T: PrimeField>(a: &[T], b: &[T], root_of_unity: T) -> Vec<T> {
  let roots = expand_root_of_unity(root_of_unity);
  let x1 = _fft(a, &roots);
  let x2 = _fft(b, &roots);

  let n = x1.len();
  let mut mul_values = Vec::with_capacity(n);
  for i in 0..n {
    mul_values.push(x1[i] * x2[i]);
  }

  let roots_rev = roots_of_unity_rev(&roots);
  _inv_fft(&mul_values, &roots_rev)
}

#[test]
fn test_mul_polys_via_fft() {
  use crate::f7::F7;

  let g2 = F7::multiplicative_generator();
  let a: Vec<F7> = [5, 0, 1].iter().map(|x| F7::from(*x)).collect();
  let b: Vec<F7> = [2, 3, 1].iter().map(|x| F7::from(*x)).collect();
  let answer: Vec<F7> = [3, 1, 0, 3, 1, 0].iter().map(|x| F7::from(*x)).collect();
  let res = mul_polys(&a, &b, g2);
  assert_eq!(res, answer);

  let a: Vec<F7> = [5, 0, 0, 1, 1].iter().map(|x| F7::from(*x)).collect();
  let b: Vec<F7> = [3, 0, 0, 2].iter().map(|x| F7::from(*x)).collect();
  let answer: Vec<F7> = [3, 2, 0, 6, 3, 0].iter().map(|x| F7::from(*x)).collect();
  let res = mul_polys(&a, &b, g2);
  // comb(a, b) = [15, 0, 0, 6, 3, 0, 2, 2]
  assert_eq!(res, answer);
}

use crate::poly_utils::multi_inv;

pub fn div_polys<T: PrimeField>(a: &[T], b: &[T], root_of_unity: T) -> Vec<T> {
  let roots = expand_root_of_unity(root_of_unity);
  let x1 = _fft(a, &roots);
  let x2 = _fft(b, &roots);
  let inv_x2: Vec<T> = multi_inv(&x2);

  let n = x1.len();
  let mut mul_values = Vec::with_capacity(n);
  for i in 0..n {
    mul_values.push(x1[i] * inv_x2[i]);
  }

  let roots_rev = roots_of_unity_rev(&roots);
  println!("{:?}", _inv_fft(&inv_x2, &roots_rev));
  _inv_fft(&mul_values, &roots_rev)
}

#[test]
fn test_div_polys_via_fft() {
  use crate::f7::F7;

  let g2 = F7::multiplicative_generator();
  let a: Vec<F7> = [3, 1, 0, 3, 1, 0].iter().map(|x| F7::from(*x)).collect();
  let b: Vec<F7> = [2, 3, 1].iter().map(|x| F7::from(*x)).collect();
  let answer: Vec<F7> = [5, 0, 1].iter().map(|x| F7::from(*x)).collect();
  let res = div_polys(&a, &b, g2);
  assert_eq!(res, answer);
}
