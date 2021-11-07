use commitment::hash::Digest;
use commitment::merkle_tree::{MerkleTree, Proof};
use commitment::pallarel_merkle_tree::ParallelMerkleTree;
use ff::PrimeField;
use ff_utils::ff_utils::{FromBytes, ToBytes};
use fri::fri::FriProof;
use fri::poly_utils::{lagrange_interp, multi_inv, sparse};
use fri::utils::{blake, get_pseudorandom_indices};
#[allow(unused_imports)]
use log::{debug, info, warn};
use num::bigint::BigUint;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::convert::TryInto;

pub fn log2_ceil(value: usize) -> u32 {
  let mut log_value = 1;
  let mut tmp = value;
  while tmp > 1 {
    tmp /= 2;
    log_value += 1;
  }

  log_value
}

pub fn parse_be_bytes_to_decimal(value: &[u8]) -> String {
  BigUint::from_bytes_be(value).to_str_radix(10)
}

pub fn u32_be_bytes_to_u8_be_bytes(values: [u32; 8]) -> [u8; 32] {
  let mut output = [0u8; 32];
  for (i, value) in values.iter().enumerate() {
    for (j, &v) in value.to_be_bytes().iter().enumerate() {
      output[4 * i + j] = v;
    }
  }

  output
}

// pub fn u32_le_bytes_to_u8_be_bytes(values: [u32; 8]) -> [u8; 32] {
//   let mut output = [0u8; 32];
//   for (i, value) in values.iter().rev().enumerate() {
//     for (j, &v) in value.to_be_bytes().iter().enumerate() {
//       output[4 * i + j] = v;
//     }
//   }

//   output
// }

pub fn mk_seed(messages: &[Vec<u8>]) -> String {
  let mut message: Vec<u8> = vec![];
  for m in messages {
    message.extend(m);
  }
  parse_be_bytes_to_decimal(&blake(&message))
}

pub fn r1cs_computational_trace<T: PrimeField>(coefficients: &[T], witness: &[T]) -> Vec<T> {
  let n_constraints = coefficients.len() / witness.len() / 3;
  let n_wires = witness.len();
  assert_eq!(coefficients.len(), 3 * n_constraints * witness.len());

  let mut computational_trace: Vec<T> = witness
    .iter()
    .cycle()
    .take(3 * n_constraints * witness.len())
    .map(|&x| x)
    .collect();
  debug_assert_eq!(computational_trace.len(), 3 * n_constraints * witness.len());

  for (i, &coeff) in coefficients.iter().enumerate() {
    if i % n_wires == 0 {
      computational_trace.push(coeff);
    } else {
      computational_trace.push(coeff * witness[i % n_wires] + computational_trace.last().unwrap());
    }
  }

  computational_trace
}

// pub fn r1cs_computational_trace<T: PrimeField + FromBytes>(
//   r1cs: &R1csContents,
//   witness: &[T],
// ) -> Vec<Vec<(u32, T)>> {
//   let n_constraints = r1cs.header.n_constraints;
//   let n_wires = r1cs.header.n_wires;

//   let R1csContents {
//     version: _,
//     header: _,
//     constraints: Constraints(constraints),
//   } = r1cs;
//   let computational_trace = vec![];
//   for Constraint { factors } in constraints {
//     let new_factors = vec![];
//     for j in 0..3 {
//       let Factor {
//         n_coefficient,
//         coefficients,
//       } = factors[j];
//       let mut acc = T::zero();
//       let accumulations = vec![];
//       for Coefficient { wire_id, value } in coefficients {
//         let value = value.iter().fold(vec![], |acc, x| {
//           acc.extend(x.to_le_bytes());
//           acc
//         });
//         acc = acc + T::from_bytes_le(value).unwrap();
//         accumulations.push((wire_id, acc))
//       }

//       computational_trace.push(accumulations);
//     }
//   }

//   computational_trace
// }

#[derive(Serialize, Deserialize, Debug)]
pub struct StarkProof<H: Digest> {
  pub m_root: H,
  pub l_root: H,
  pub a_root: H,
  pub main_branches: Vec<Proof<Vec<u8>, H>>,
  pub linear_comb_branches: Vec<Proof<Vec<u8>, H>>,
  pub fri_proof: Vec<FriProof<H>>,
}

// const modulus = 2**256 - 2**32 * 351 + 1;
// const non_residue = 7;
pub const LOG_EXTENSION_FACTOR: usize = 3usize;
pub const EXTENSION_FACTOR: usize = 8usize; // >= 4 (for low-degree proof)
pub const SPOT_CHECK_SECURITY_FACTOR: usize = 80usize;

pub fn calc_max_log_precision<T: PrimeField + ToBytes>() -> u32 {
  let mut ff_order_be = (T::zero() - T::one()).to_bytes_be().unwrap();
  let mut log_max_precision = 0;
  loop {
    match ff_order_be.pop() {
      Some(0) => {
        log_max_precision += 8;
      }
      Some(rem) => {
        debug_assert!(rem != 0);
        let mut rem2 = rem;
        while rem2 % 2 == 0 {
          log_max_precision += 1;
          rem2 /= 2;
        }
        break;
      }
      None => {
        break;
      }
    }
  }

  log_max_precision
}

pub fn convert_usize_iter_to_ff_vec<T: PrimeField + FromBytes, I: IntoIterator<Item = usize>>(
  iter: I,
) -> Vec<T> {
  iter
    .into_iter()
    .map(|v| T::from_bytes_le(v.to_le_bytes().as_ref()).unwrap())
    .collect::<Vec<_>>()
}

// Z(X) = X^steps - 1
pub fn calc_z_polynomial<T: PrimeField>(steps: usize) -> Vec<T> {
  let mut sparse_z1 = HashMap::new();
  sparse_z1.insert(0, -T::one());
  sparse_z1.insert(steps, T::one());
  sparse(sparse_z1)
}

// Q1(j) = F0(j) * (P(j) - F1(j) * P(j - 1) - K(j) * S(j))
pub fn calc_q1_evaluations<T: PrimeField>(
  s_evaluations: &[T],
  k_evaluations: &[T],
  p_evaluations: &[T],
  f0_evaluations: &[T],
  f1_evaluations: &[T],
  precision: usize,
  skips: usize,
) -> Vec<T> {
  let mut q1_evaluations = vec![];
  for j in 0..precision {
    let s_of_x = s_evaluations[j];
    let k_of_x = k_evaluations[j];
    let p_of_prev_x = p_evaluations[(j + precision - skips) % precision];
    let p_of_x = p_evaluations[j];
    let f0 = f0_evaluations[j];
    let f1 = f1_evaluations[j];

    let q1_of_x = f0 * (p_of_x - f1 * p_of_prev_x - k_of_x * s_of_x);
    q1_evaluations.push(q1_of_x);
    // if j % skips == 0 {
    //   println!(
    //     "{:?}, {:?}, {:?}, {:?}",
    //     f0,
    //     f1,
    //     p_of_x - k_of_x * s_of_x,
    //     p_of_prev_x
    //   );
    // }
  }

  q1_evaluations
}

// Q2(j) = F2(j) * (P(j + 2k) - P(j + k) * P(j))
// where k := original_steps / 3;
pub fn calc_q2_evaluations<T: PrimeField>(
  p_evaluations: &[T],
  f2_evaluations: &[T],
  precision: usize,
  skips: usize,
  original_steps: usize,
) -> Vec<T> {
  let mut q2_evaluations = vec![];
  for j in 0..precision {
    let j = j;
    let j2 = (j + original_steps / 3 * skips) % precision;
    let j3 = (j + original_steps / 3 * 2 * skips) % precision;
    let a_eval = p_evaluations[j];
    let b_eval = p_evaluations[j2];
    let c_eval = p_evaluations[j3];
    let f2 = f2_evaluations[j];

    let q2_of_x = f2 * (c_eval - a_eval * b_eval);
    q2_evaluations.push(q2_of_x);
    // if j % skips == 0 {
    //   println!(
    //     "{:?}, {:?}, {:?}, {:?}",
    //     c_eval - a_eval * b_eval,
    //     f0,
    //     f2,
    //     q2_of_x
    //   );
    // }
  }

  q2_evaluations
}

pub fn get_accumulator_tree_root<T: PrimeField + ToBytes, H: Digest>(
  permuted_indices: &[usize],
  witness_trace: &[T],
) -> H {
  let accumulator_str = permuted_indices
    .iter()
    .zip(witness_trace)
    .map(|(&p_val, &a_val)| {
      let mut res = vec![];
      res.extend(p_val.to_le_bytes());
      res.extend(a_val.to_bytes_le().unwrap());
      res
    })
    .collect::<Vec<_>>();

  let mut a_tree: ParallelMerkleTree<Vec<u8>, H> = ParallelMerkleTree::new(accumulator_str);
  a_tree.update();
  a_tree.gen_proofs(&[]);
  let a_root = a_tree.get_root().unwrap();

  a_root
}

pub fn get_random_ff_values<T: PrimeField + FromBytes>(
  seed: &[u8],
  modulus: u32,
  size: usize,
  exclude_multiples_of: u32,
) -> Vec<T> {
  let accumulator_randomness =
    get_pseudorandom_indices(seed, modulus, size * 8, exclude_multiples_of);
  // println!("accumulator_randomness: {:?}", accumulator_randomness);
  let random_values = accumulator_randomness
    .chunks(8)
    .map(|rand| T::from_bytes_le(&u32_be_bytes_to_u8_be_bytes(rand.try_into().unwrap())).unwrap())
    .collect::<Vec<_>>();
  debug_assert_eq!(random_values.len(), size);

  random_values
}

// A(g1^j) = A(g1^(j-1)) * val_nmr / val_dnm with A(g1^(-1)) = 1
pub fn calc_a_mini_evaluations<T: PrimeField>(
  witness_trace: &[T],
  ext_indices: &[T],
  ext_permuted_indices: &[T],
  random_values: &[T],
  steps: usize,
  skips: usize,
) -> Vec<T> {
  let r = random_values;
  let mut a_nmr_evaluations: Vec<T> = vec![];
  let mut a_dnm_evaluations: Vec<T> = vec![];
  let mut val_nmr_list = vec![];
  let mut val_dnm_list = vec![];
  for j in 0..steps {
    let last_acc_nmr = if j != 0 {
      a_nmr_evaluations.last().unwrap().clone()
    } else {
      T::one()
    };
    let last_acc_dnm = if j != 0 {
      a_dnm_evaluations.last().unwrap().clone()
    } else {
      T::one()
    };
    let val_nmr = r[0] + r[1] * ext_indices[j * skips] + r[2] * witness_trace[j];
    let val_dnm = r[0] + r[1] * ext_permuted_indices[j * skips] + r[2] * witness_trace[j];
    let acc_nmr = val_nmr * last_acc_nmr;
    let acc_dnm = val_dnm * last_acc_dnm;

    val_nmr_list.push(val_nmr);
    val_dnm_list.push(val_dnm);
    a_nmr_evaluations.push(acc_nmr);
    a_dnm_evaluations.push(acc_dnm);
    if val_nmr == T::zero() || val_dnm == T::zero() {
      println!("value is zero: {:?}", j);
    }
  }

  let inv_a_dnm_evaluations = multi_inv(&a_dnm_evaluations);
  let a_mini_evaluations: Vec<T> = a_nmr_evaluations
    .iter()
    .zip(&inv_a_dnm_evaluations)
    .map(|(&a_nmr, &inv_a_dnm)| a_nmr * inv_a_dnm)
    .collect();

  a_mini_evaluations
}

// Q3(g2^j) = A(g2^j) * val_dnm - A(g2^(j-skips)) * val_nmr
// where val_nmr = r0 + r1 *          ext_indices[j] + r2 * S(g2^j)
//       val_dnm = r0 + r1 * ext_permuted_indices[j] + r2 * S(g2^j)
pub fn calc_q3_evaluations<T: PrimeField>(
  s_evaluations: &[T],
  a_evaluations: &[T],
  ext_indices: &[T],
  ext_permuted_indices: &[T],
  random_values: &[T],
  precision: usize,
  skips: usize,
) -> Vec<T> {
  let r = random_values;
  assert!(r.len() >= 3);

  let mut q3_evaluations = vec![];
  for j in 0..precision {
    // A(j) * val_dnm = A(j - 1) * val_nmr
    let val_nmr = r[0] + r[1] * ext_indices[j] + r[2] * s_evaluations[j];
    let val_dnm = r[0] + r[1] * ext_permuted_indices[j] + r[2] * s_evaluations[j];
    let prev_j = (j + precision - skips) % precision;
    let q3_of_x = a_evaluations[j] * val_dnm - a_evaluations[prev_j] * val_nmr;
    q3_evaluations.push(q3_of_x);
    // if j % skips == 0 {
    //   println!(
    //     "a {:?}, {:?}, {:?}, {:?}",
    //     val_nmr * val_dnm.invert().unwrap(),
    //     a_evaluations[prev_j],
    //     a_evaluations[prev_j] * val_nmr * val_dnm.invert().unwrap(),
    //     a_evaluations[j],
    //   );
    // }
  }

  q3_evaluations
}

// Compute D1(x) = Q1(x) / Z(x)
pub fn calc_d1_polynomial<T: PrimeField>(q1_evaluations: &[T], inv_z_evaluations: &[T]) -> Vec<T> {
  for (pos, (&q, &z)) in q1_evaluations.iter().zip(inv_z_evaluations).enumerate() {
    if z == T::zero() {
      assert_eq!(q, T::zero(), "invalid D1: {:?} {:?} {:?}", pos, q, z);
    }
  }
  q1_evaluations
    .iter()
    .zip(inv_z_evaluations)
    .map(|(&q1, &z1i)| q1 * z1i)
    .collect()
}

// Compute D2(x) = Q2(x) / Z(x)
pub fn calc_d2_polynomial<T: PrimeField>(q2_evaluations: &[T], inv_z_evaluations: &[T]) -> Vec<T> {
  for (pos, (&q, &z)) in q2_evaluations.iter().zip(inv_z_evaluations).enumerate() {
    if z == T::zero() {
      assert_eq!(q, T::zero(), "invalid D2: {:?} {:?} {:?}", pos, q, z);
    }
  }
  q2_evaluations
    .iter()
    .zip(inv_z_evaluations)
    .map(|(&q2, &z2i)| q2 * z2i)
    .collect()
}

// Compute D3(x) = Q3(x) / Z(x)
pub fn calc_d3_polynomial<T: PrimeField>(q3_evaluations: &[T], inv_z_evaluations: &[T]) -> Vec<T> {
  for (pos, (&q, &z)) in q3_evaluations.iter().zip(inv_z_evaluations).enumerate() {
    if z == T::zero() {
      assert_eq!(q, T::zero(), "invalid D3: {:?} {:?} {:?}", pos, q, z);
    }
  }
  q3_evaluations
    .iter()
    .zip(inv_z_evaluations)
    .map(|(&q3, &z3i)| q3 * z3i)
    .collect()
}

// I2(X) where I2(g1^w) = public_wires[k] for all (k, w) in public_first_indices
pub fn calc_i2_polynomial<T: PrimeField>(
  public_first_indices: &[(usize, usize)],
  xs: &[T],
  public_wires: &[T],
  skips: usize,
) -> Vec<T> {
  let mut x_vals: Vec<T> = vec![];
  let mut y_vals: Vec<T> = vec![];
  for (k, w) in public_first_indices {
    x_vals.push(xs[skips * w]);
    y_vals.push(public_wires[*k]);
  }

  lagrange_interp(&x_vals, &y_vals)
}

// Zb2(X) where Zb2(g1^w) = 0 for all (k, w) in public_first_indices
pub fn calc_zb2_evaluations<T: PrimeField>(
  public_first_indices: &[(usize, usize)],
  xs: &[T],
  precision: usize,
  skips: usize,
) -> Vec<T> {
  let mut zb2_evaluations = vec![T::one(); precision];
  for (_, w) in public_first_indices {
    let j = w * skips;
    zb2_evaluations = zb2_evaluations
      .iter()
      .enumerate()
      .map(|(i, &val)| val * (xs[i] - xs[j]))
      .collect();
  }

  zb2_evaluations
}

// I3(X) where I3(g1^(-1)) = 1
pub fn calc_i3_polynomial<T: PrimeField>(xs: &[T], skips: usize) -> Vec<T> {
  let x_of_last_step = xs[xs.len() - skips];
  let x_vals = vec![x_of_last_step];
  let y_vals = vec![T::one()];
  lagrange_interp(&x_vals, &y_vals)
}

// Zb3(X) = X - g1^(-1)
pub fn calc_zb3_evaluations<T: PrimeField>(xs: &[T], precision: usize, skips: usize) -> Vec<T> {
  let x_of_last_step = xs[xs.len() - skips];
  let zb3_evaluations = vec![T::one(); precision];
  zb3_evaluations
    .iter()
    .enumerate()
    .map(|(i, &val)| val * (xs[i] - x_of_last_step))
    .collect()
}

// B2(x) = (S(x) - I2(x)) / Z_b2(x)
pub fn calc_b2_evaluations<T: PrimeField>(
  s_evaluations: &[T],
  i2_evaluations: &[T],
  inv_zb2_evaluations: &[T],
) -> Vec<T> {
  for (pos, ((&zb2, &s), &i2)) in inv_zb2_evaluations
    .iter()
    .zip(s_evaluations)
    .zip(i2_evaluations)
    .enumerate()
  {
    if zb2 == T::zero() {
      assert_eq!(s, i2, "invalid B2: {:?} {:?} {:?}", pos, s, i2);
    }
  }

  s_evaluations
    .iter()
    .zip(i2_evaluations)
    .zip(inv_zb2_evaluations)
    .map(|((&s, &i2), &inv_zb)| (s - i2) * inv_zb)
    .collect()
}

// B3(x) = (A(x) - I3(x)) / Z_b3(x)
pub fn calc_b3_evaluations<T: PrimeField>(
  a_evaluations: &[T],
  i3_evaluations: &[T],
  inv_zb3_evaluations: &[T],
) -> Vec<T> {
  for (pos, ((&zb3, &a), &i3)) in inv_zb3_evaluations
    .iter()
    .zip(a_evaluations)
    .zip(i3_evaluations)
    .enumerate()
  {
    if zb3 == T::zero() {
      assert_eq!(a, i3, "invalid B3: {:?} {:?} {:?}", pos, a, i3);
    }
  }

  a_evaluations
    .iter()
    .zip(i3_evaluations)
    .zip(inv_zb3_evaluations)
    .map(|((&a, &i3), &inv_zb)| (a - i3) * inv_zb)
    .collect()
}
