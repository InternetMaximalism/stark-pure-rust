use std::convert::TryInto;
use blake2::{Blake2s, Digest};
use ff::PrimeField;

pub fn blake(message: &[u8]) -> String {
  let mut hasher = Blake2s::new();
  hasher.update(message);
  let res = hasher.finalize();
  format!("{:x}", res)
}

#[test]
pub fn test_blake() {
  let message = b"hello world";
  let res = blake(message);
  let answer = "9aec6806794561107e594b1f6a8a6b0c92a0cba9acf5e5e93cca06f781813b0b";
  assert_eq!(res, answer.to_string());
}

pub fn is_a_power_of_2(x: i32) -> bool {
  if x == 1 {
    true
  } else if x % 2 == 1 {
    false
  } else {
    is_a_power_of_2(x / 2)
  }
}

pub fn get_pseudorandom_indices(
  seed: String,
  modulus: i32,
  count: usize,
  exclude_multiples_of: i32,
) -> Vec<usize> {
  assert!(modulus < 2i32.pow(24));
  let mut data = seed.clone();
  while data.len() < 4 * count {
    data += &blake(&data[(data.len() - 32)..].as_ref());
  }
  println!("{:?}", data);
  if exclude_multiples_of == 0 {
    (0..(count * 4))
      .step_by(4)
      .map(|i| i32::from_str_radix(&data[i..(i + 4)], 16).unwrap() % modulus)
      .map(|i| i.try_into().unwrap())
      .collect()
  } else {
    let real_modulus = modulus * (exclude_multiples_of - 1) / exclude_multiples_of;
    (0..(count * 4))
      .step_by(4)
      .map(|i| i32::from_str_radix(&data[i..(i + 4)], 16).unwrap() % real_modulus)
      .map(|i| i + 1 + i / (exclude_multiples_of - 1))
      .map(|i| i.try_into().unwrap())
      .collect()
  }
}

#[test]
fn test_get_pseudorandom_indices() {
  let res = get_pseudorandom_indices(blake(b"hello world"), 7, 5, 0);
  let answer = [5, 2, 0, 5, 5];
  assert_eq!(res, answer);

  let res = get_pseudorandom_indices(blake(b"hello another world"), 7, 20, 0);
  let answer = [6, 5, 4, 6, 6, 4, 6, 6, 3, 5, 0, 1, 4, 3, 3, 3, 5, 1, 1, 5];
  assert_eq!(res, answer);
}
