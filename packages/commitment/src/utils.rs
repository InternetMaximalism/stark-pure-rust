use blake2::*;
use std::cmp::min;
use std::convert::TryInto;

pub fn blake(message: &[u8]) -> Vec<u8> {
  let mut hasher = Blake2s::new();
  hasher.update(message);
  let res = hasher.finalize();
  hex::decode(format!("{:x}", res)).unwrap()
}

#[test]
pub fn test_blake() {
  let message = b"hello world";
  let res = blake(message);
  let answer =
    hex::decode("9aec6806794561107e594b1f6a8a6b0c92a0cba9acf5e5e93cca06f781813b0b").unwrap();
  assert_eq!(res, answer);

  let res2 = blake(&answer);
  let answer2 =
    hex::decode("8ea974646c2be3c16f9f52a2e5ebb3d2df7ba184a6440e47fc6fcce6e9d9bdc4").unwrap();
  assert_eq!(res2, answer2);
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

pub fn log2_ceil(value: usize) -> usize {
  if value == 0 {
    panic!("The first argument must be a positive number.")
  }

  if value == 1 {
    return 0;
  }

  let mut log_value = 1;
  let mut tmp_value = value - 1;
  while tmp_value > 1 {
    tmp_value /= 2;
    log_value += 1;
  }

  log_value
}

#[test]
#[should_panic]
pub fn test_error_log2_ceil() {
  log2_ceil(0);
}

#[test]
pub fn test_log2_ceil() {
  let res0 = log2_ceil(1);
  assert_eq!(res0, 0);

  let res1 = log2_ceil(2);
  assert_eq!(res1, 1);

  let res2 = log2_ceil(3);
  assert_eq!(res2, 2);

  let res3 = log2_ceil(4);
  assert_eq!(res3, 2);

  let res4 = log2_ceil(5);
  assert_eq!(res4, 3);

  let res5 = log2_ceil(127);
  assert_eq!(res5, 7);
}

pub fn get_pseudorandom_indices(
  seed: &[u8],
  modulus: u32,
  count: usize,
  exclude_multiples_of: u32,
) -> Vec<u32> {
  assert!(modulus < 2u32.pow(24));
  let mut data = seed.to_vec();
  while data.len() < 4 * count {
    let new_data = blake(&data[(data.len() - 32)..]);
    data.extend(&new_data);
  }
  if exclude_multiples_of == 0 {
    return (0..(count * 4))
      .step_by(4)
      .map(|i| u32::from_be_bytes([data[i], data[i + 1], data[i + 2], data[i + 3]]) % modulus)
      .collect();
  }

  let real_modulus = modulus * (exclude_multiples_of - 1) / exclude_multiples_of;
  return (0..(count * 4))
    .step_by(4)
    .map(|i| u32::from_be_bytes([data[i], data[i + 1], data[i + 2], data[i + 3]]) % real_modulus)
    .map(|i| i + 1 + i / (exclude_multiples_of - 1))
    .collect();
}

#[test]
fn test_get_pseudorandom_indices() {
  let res = get_pseudorandom_indices(&blake(b"hello world"), 7, 5, 0);
  let answer = [5, 5, 5, 3, 5];
  assert_eq!(res, answer);

  let res = get_pseudorandom_indices(&blake(b"hello another world"), 7, 20, 0);
  let answer = [3, 0, 2, 4, 4, 1, 4, 2, 5, 1, 3, 2, 1, 0, 0, 1, 6, 5, 2, 3];
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
