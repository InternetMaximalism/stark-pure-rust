use blake2::{Blake2s, Digest};
use std::cmp::min;
use std::convert::TryInto;

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

pub fn is_a_power_of_2(x: usize) -> bool {
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

fn get_power_sequence(x: u64, steps: usize) -> Vec<u64> {
  let mut powers = vec![1u64];
  for _ in 1..steps {
    powers.push(powers.last().unwrap() * x);
  }
  powers
}

pub fn parse_bytes_to_u64_vec(mut xs: &[u8]) -> Vec<u64> {
  let mut output = vec![];
  while xs.len() > 0 {
    let mut y = 0u64;
    let sub_xs: Vec<u64> = xs
      .iter()
      .take(8)
      .map(|x| (*x).try_into().unwrap())
      .collect();
    for (power, x) in get_power_sequence(256, 8).iter().zip(sub_xs) {
      y += power * x;
    }
    output.push(y);
    xs = &xs[min(8, xs.len())..];
  }
  output
}

#[test]
fn test_parse_bytes_to_u64_vec() {
  let xs = &[1u8, 1, 0, 0, 0, 0, 0, 0, 255, 0];
  let res = parse_bytes_to_u64_vec(xs);
  let answer = &[257u64, 255u64];
  assert_eq!(&res, answer);
}
