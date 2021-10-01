use std::collections::HashMap;

use crate::utils::*;
use ff::PrimeField;
use ff_utils::ff_utils::{FromBytes, ToBytes};
use fri::delayed::Delayed;
use fri::fft::{best_fft, expand_root_of_unity, inv_best_fft};
use fri::fri::prove_low_degree;
use fri::lazily;
use fri::merkle_tree3::{gen_multi_proofs_multi_core, BlakeDigest};
use fri::multicore::Worker;
use fri::poly_utils::{eval_poly_at, lagrange_interp, multi_inv, sparse};
use fri::utils::{get_pseudorandom_indices, parse_bytes_to_u64_vec};
use num::bigint::BigUint;
use std::convert::TryInto;

pub fn mk_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
  witness_trace: &[T],
  computational_trace: &[T],
  public_wires: &[T],
  public_first_indices: &[(usize, usize)],
  permuted_indices: &[usize],
  coefficients: &[T], // This argument may be good to use HashMap<usize, T> instead of Vec<T> because we can omit zero coefficients from it.
  flag0: &[T],
  flag1: &[T],
  flag2: &[T],
  n_constraints: usize,
  n_wires: usize,
) -> StarkProof {
  // println!("n_cons: {:?}", n_constraints);
  // println!("n_wires: {:?}", n_wires);
  // println!("n_public_wires: {:?}", public_wires.len());
  // println!("last_coeff_list: {:?}", last_coeff_list);
  let original_steps = coefficients.len();
  println!("original_steps: {:?}", original_steps);
  assert!(original_steps <= 3 * n_constraints * n_wires);
  assert!(original_steps % 3 == 0);
  assert_eq!(witness_trace.len(), original_steps);
  assert_eq!(computational_trace.len(), original_steps);

  let mut log_steps = 1;
  let mut tmp_steps = original_steps - 1;
  while tmp_steps > 1 {
    tmp_steps /= 2;
    log_steps += 1;
  }
  let mut steps = 2usize.pow(log_steps);
  if steps < 8 {
    steps = 8;
  }

  let precision = steps * EXTENSION_FACTOR;
  println!("precision: {:?}", precision);
  println!(
    "original_steps * skips {:?}",
    original_steps * EXTENSION_FACTOR
  );

  let log_precision = log_steps + LOG_EXTENSION_FACTOR as u32;
  {
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
    println!("max_precision: 2^{:?}", log_max_precision);
    assert!(precision <= 2usize.pow(log_max_precision));
  }

  let mut permuted_indices = permuted_indices.to_vec();
  // println!("permuted_indices.len(): {:?}", permuted_indices.len());
  permuted_indices.extend(original_steps..steps);
  // println!("permuted_indices.len(): {:?}", permuted_indices.len());

  let mut witness_trace = witness_trace.to_vec();
  // println!("witness_trace.len(): {:?}", witness_trace.len());
  witness_trace.extend(vec![T::zero(); steps - original_steps]);
  // println!("witness_trace.len(): {:?}", witness_trace.len());

  // println!("computational_trace: {:?}", computational_trace);
  let mut computational_trace = computational_trace.to_vec();
  computational_trace.extend(vec![T::zero(); steps - original_steps]);

  let mut coefficients = coefficients.to_vec();
  coefficients.extend(vec![T::zero(); steps - original_steps]);

  // Root of unity such that x^precision=1

  let ff_order = T::zero() - T::one();

  let times_nmr = BigUint::from_bytes_le(&ff_order.to_bytes_le().unwrap());
  let times_dnm = BigUint::from_bytes_le(&precision.to_le_bytes());
  // println!("{:?} {:?}", &times_nmr, &times_dnm);
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  {
    let times = parse_bytes_to_u64_vec(&times_nmr.to_bytes_le()); // modulus - 1
    let unity = T::multiplicative_generator().pow_vartime(&times);
    debug_assert_eq!(unity, T::one());
  }

  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus
  println!("expand g2");
  let xs = expand_root_of_unity(g2); // Powers of the higher-order root of unity
  let skips = precision / steps; // EXTENSION_FACTOR
  let g1 = xs[skips]; // Root of unity such that x^steps=1
  let log_order_of_g1 = log_steps as u32;
  let log_order_of_g2 = log_precision as u32;

  // Interpolate the computational trace into a polynomial P, with each step
  // along a successive power of g1
  println!("calculate expanding polynomials");

  let worker = Worker::new();

  let s_polynomial = inv_best_fft(witness_trace.clone(), &g1, &worker, log_order_of_g1); // S(X)
  let s_evaluations = best_fft(s_polynomial, &g2, &worker, log_order_of_g2);
  // println!("s_evaluations: {:?}", s_evaluations);
  println!("Converted witness trace into a polynomial and low-degree extended it");

  // println!("trace: {:?}", computational_trace);
  let p_polynomial = inv_best_fft(computational_trace, &g1, &worker, log_order_of_g1); // P(X)
  let p_evaluations = best_fft(p_polynomial, &g2, &worker, log_order_of_g2);
  // println!("p_evaluations: {:?}", p_evaluations);
  println!("Converted computational trace into a polynomial and low-degree extended it");

  // println!("coeff: {:?}", coefficients);
  let k_polynomial = inv_best_fft(coefficients, &g1, &worker, log_order_of_g1); // K(X)
  let k_evaluations = best_fft(k_polynomial, &g2, &worker, log_order_of_g2);
  // println!("k: {:?}", k_evaluations);
  println!("Converted coefficients into a polynomial and low-degree extended it");

  let f0_polynomial = inv_best_fft(flag0.to_vec(), &g1, &worker, log_order_of_g1);
  let f0_evaluations = best_fft(f0_polynomial.clone(), &g2, &worker, log_order_of_g2);
  // println!("f0_evaluations: {:?}", f0_evaluations);

  let f1_polynomial = inv_best_fft(flag1.to_vec(), &g1, &worker, log_order_of_g1);
  let f1_evaluations = best_fft(f1_polynomial, &g2, &worker, log_order_of_g2);
  // println!("f1_evaluations: {:?}", f1_evaluations);

  let f2_polynomial = inv_best_fft(flag2.to_vec(), &g1, &worker, log_order_of_g1);
  let f2_evaluations = best_fft(f2_polynomial, &g2, &worker, log_order_of_g2);
  // println!("f2_evaluations: {:?}", f2_evaluations);

  let z1_polynomial = {
    let mut sparse_z1 = HashMap::new();
    sparse_z1.insert(0, -T::one());
    sparse_z1.insert(steps, T::one());
    sparse(sparse_z1)
  };
  let z1_evaluations = best_fft(z1_polynomial, &g2, &worker, log_order_of_g2);
  // println!("z1_evaluations: {:?}", z1_evaluations);

  let mut q1_evaluations = vec![];
  for j in 0..precision {
    let s_of_x = s_evaluations[j];
    let k_of_x = k_evaluations[j];
    let p_of_prev_x = p_evaluations[(j + precision - skips) % precision];
    let p_of_x = p_evaluations[j];
    let f0 = f0_evaluations[j];
    let f1 = f1_evaluations[j];

    // Q1(j) = F0(j) * (P(j) - F1(j) * P(j - 1) - K(j) * S(j))
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
  // println!("q1_evaluations: {:?}", q1_evaluations);
  println!("Computed Q1 polynomial");

  let mut q2_evaluations = vec![];
  for j in 0..precision {
    let j = j;
    let j2 = (j + original_steps / 3 * skips) % precision;
    let j3 = (j + original_steps / 3 * 2 * skips) % precision;
    let a_eval = p_evaluations[j];
    let b_eval = p_evaluations[j2];
    let c_eval = p_evaluations[j3];
    let f2 = f2_evaluations[j];

    // Q2(j) = F0(j) * F2(j) * (P(j + 2k) - P(j + k) * P(j))
    // where k := original_steps / 3;
    let q2_of_x = f2 * (c_eval - a_eval * b_eval);
    // if j % skips == 0 {
    //   println!(
    //     "{:?}, {:?}, {:?}, {:?}",
    //     c_eval - a_eval * b_eval,
    //     f0,
    //     f2,
    //     q2_of_x
    //   );
    // }
    q2_evaluations.push(q2_of_x);
  }
  // println!("f0_evaluations: {:?}", f0_evaluations);
  // println!("f2_evaluations: {:?}", f2_evaluations);
  // println!("q2_evaluations: {:?}", q2_evaluations);
  println!("Computed Q2 polynomial");

  let converted_indices = (0..steps)
    .map(|v| T::from_bytes_le(v.to_le_bytes().as_ref()).unwrap())
    .collect::<Vec<_>>();
  let index_polynomial = inv_best_fft(converted_indices, &g1, &worker, log_order_of_g1);
  let ext_indices = best_fft(index_polynomial, &g2, &worker, log_order_of_g2);
  println!("Computed extended index polynomial");

  let converted_permuted_indices = permuted_indices
    .clone()
    .iter()
    .map(|v| T::from_bytes_le(v.to_le_bytes().as_ref()).unwrap())
    .collect::<Vec<_>>();
  let permuted_polynomial = inv_best_fft(converted_permuted_indices, &g1, &worker, log_order_of_g1);
  let ext_permuted_indices = best_fft(permuted_polynomial, &g2, &worker, log_order_of_g2);
  // println!("ext_permuted_indices: {:?}", ext_permuted_indices);
  println!("Computed extended permuted index polynomial");

  let accumulator_str = permuted_indices
    .iter()
    .zip(&witness_trace)
    .map(|(&p_val, &a_val)| {
      Delayed::new(move || {
        let mut res = vec![];
        res.extend(p_val.to_le_bytes());
        res.extend(a_val.to_bytes_le().unwrap());
        res
      })
    })
    .collect::<Vec<_>>();
  let (_, a_root) =
    gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(&accumulator_str, &[], &worker);
  let accumulator_randomness =
    get_pseudorandom_indices(a_root.as_ref(), precision as u32, 3 * 8, 0);
  println!("accumulator_randomness: {:?}", accumulator_randomness);

  let mut a_nmr_evaluations: Vec<T> = vec![];
  let mut a_dnm_evaluations: Vec<T> = vec![];
  let mut val_nmr_list = vec![];
  let mut val_dnm_list = vec![];
  let r1 = T::from_bytes_le(&u32_be_bytes_to_u8_be_bytes(
    accumulator_randomness[0..8].try_into().unwrap(),
  ))
  .unwrap();
  let r2 = T::from_bytes_le(&u32_be_bytes_to_u8_be_bytes(
    accumulator_randomness[8..16].try_into().unwrap(),
  ))
  .unwrap();
  let r3 = T::from_bytes_le(&u32_be_bytes_to_u8_be_bytes(
    accumulator_randomness[16..24].try_into().unwrap(),
  ))
  .unwrap();
  println!("r: {:?} {:?} {:?}", r1, r2, r3);

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
    let val_nmr = r1 + r2 * ext_indices[j * skips] + r3 * witness_trace[j];
    let val_dnm = r1 + r2 * ext_permuted_indices[j * skips] + r3 * witness_trace[j];
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
    .map(|(&a_nmr, &a_dnm)| a_nmr * a_dnm)
    .collect();

  let a_polynomials = inv_best_fft(a_mini_evaluations, &g1, &worker, log_order_of_g1);
  let a_evaluations = best_fft(a_polynomials, &g2, &worker, log_order_of_g2);

  let mut q3_evaluations = vec![];
  for j in 0..precision {
    // A(j) * val_dnm = A(j - 1) * val_nmr
    let val_nmr = r1 + r2 * ext_indices[j] + r3 * s_evaluations[j];
    let val_dnm = r1 + r2 * ext_permuted_indices[j] + r3 * s_evaluations[j];
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

  let z3_polynomial = {
    let mut sparse_z3_nmr = HashMap::new();
    sparse_z3_nmr.insert(0, -T::one());
    sparse_z3_nmr.insert(steps, T::one());
    sparse(sparse_z3_nmr)
  };
  let z3_evaluations = best_fft(z3_polynomial, &g2, &worker, log_order_of_g2);

  // Compute D(x) = Q(x) / Z(x)
  let inv_z1_evaluations = multi_inv(&z1_evaluations);
  // let inv_z2_evaluations = multi_inv(&z2_evaluations);
  let inv_z3_evaluations = multi_inv(&z3_evaluations);
  println!("Computed Q3 polynomial");

  let d1_evaluations: Vec<T> = q1_evaluations
    .iter()
    .zip(&inv_z1_evaluations)
    .map(|(&q1, &z1i)| q1 * z1i)
    .collect();
  // println!("d1_evaluations: {:?}", d1_evaluations);

  let d2_evaluations: Vec<T> = q2_evaluations
    .iter()
    .zip(&inv_z1_evaluations)
    .map(|(&q2, &z2i)| q2 * z2i)
    .collect();

  let d3_evaluations: Vec<T> = q3_evaluations
    .iter()
    .zip(&inv_z3_evaluations)
    .map(|(&q3, &z3i)| q3 * z3i)
    .collect();

  println!("Computed D polynomial");

  // let interpolant = {
  //   let mut x_vals = vec![];
  //   let mut y_vals = vec![];
  //   for (j, n_coeff) in last_coeff_list[0..public_wires.len()]
  //     .to_vec()
  //     .iter()
  //     .enumerate()
  //   {
  //     x_vals.push(xs[n_coeff * skips]);
  //     y_vals.push(public_wires[j]);
  //   }

  //   lagrange_interp(&x_vals, &y_vals)
  // };
  // let i_evaluations: Vec<T> = xs.iter().map(|&x| eval_poly_at(&interpolant, x)).collect();
  // OR
  // i_evaluations = best_fft(interpolant, modulus, g2)

  let interpolant2 = {
    let mut x_vals: Vec<T> = vec![];
    let mut y_vals: Vec<T> = vec![];
    for (k, w) in public_first_indices {
      x_vals.push(xs[skips * w]);
      y_vals.push(public_wires[*k]);
    }
    lagrange_interp(&x_vals, &y_vals)
  };

  let i2_evaluations: Vec<T> = xs.iter().map(|&x| eval_poly_at(&interpolant2, x)).collect();

  let x_of_last_step = xs[(steps - 1) * skips];
  let interpolant3 = {
    let x_vals = vec![x_of_last_step];
    let y_vals = vec![T::one()];
    lagrange_interp(&x_vals, &y_vals)
  };

  let i3_evaluations: Vec<T> = xs.iter().map(|&x| eval_poly_at(&interpolant3, x)).collect();

  // let mut zb_evaluations = vec![T::one(); precision];
  // for w in 0..n_wires {
  //   let j = w * skips;
  //   zb_evaluations = zb_evaluations
  //     .iter()
  //     .enumerate()
  //     .map(|(i, &val)| val * (xs[i] - xs[j]))
  //     .collect();
  // }
  // for k in 0..n_constraints {
  //   let j = last_coeff_list[k] * skips;
  //   zb_evaluations = zb_evaluations
  //     .iter()
  //     .enumerate()
  //     .map(|(i, &val)| val * (xs[i] - xs[j]))
  //     .collect();
  // }
  // let inv_zb_evaluations: Vec<T> = multi_inv(&zb_evaluations);

  let mut zb2_evaluations = vec![T::one(); precision];
  for (_, w) in public_first_indices {
    let j = w * skips;
    zb2_evaluations = zb2_evaluations
      .iter()
      .enumerate()
      .map(|(i, &val)| val * (xs[i] - xs[j]))
      .collect();
  }
  let inv_zb2_evaluations: Vec<T> = multi_inv(&zb2_evaluations);

  let mut zb3_evaluations = vec![T::one(); precision];
  zb3_evaluations = zb3_evaluations
    .iter()
    .enumerate()
    .map(|(i, &val)| val * (xs[i] - x_of_last_step))
    .collect();
  let inv_zb3_evaluations: Vec<T> = multi_inv(&zb3_evaluations);

  // B(x) = (P(x) - I(x)) / Z_b(x)
  // let b_evaluations: Vec<T> = p_evaluations
  //   .iter()
  //   .zip(&i_evaluations)
  //   .zip(&inv_zb_evaluations)
  //   .map(|((&p, &i), &inv_zb)| (p - i) * inv_zb)
  //   .collect();
  // B2(x) = (S(x) - I2(x)) / Z_b2(x)
  let b2_evaluations: Vec<T> = s_evaluations
    .iter()
    .zip(&i2_evaluations)
    .zip(&inv_zb2_evaluations)
    .map(|((&s, &i), &inv_zb)| (s - i) * inv_zb)
    .collect();
  // B3(x) = (A(x) - I3(x)) / Z_b3(x)
  for (pos, ((&zb, &a), &i)) in zb3_evaluations
    .iter()
    .zip(&a_evaluations)
    .zip(&i3_evaluations)
    .enumerate()
  {
    if zb == T::zero() {
      if a != i {
        println!("pos: {:?} {:?} {:?}", pos, a, i);
      }
    }
  }
  let b3_evaluations: Vec<T> = a_evaluations
    .iter()
    .zip(&i3_evaluations)
    .zip(&inv_zb3_evaluations)
    .map(|((&a, &i), &inv_zb)| (a - i) * inv_zb)
    .collect();
  println!("Computed B polynomial");

  // Compute their Merkle root
  let poly_evaluations_str = p_evaluations
    .iter()
    .zip(&a_evaluations)
    .zip(&s_evaluations)
    .zip(&d1_evaluations)
    .zip(&d2_evaluations)
    .zip(&d3_evaluations)
    .zip(&b2_evaluations)
    .zip(&b3_evaluations)
    .map(
      |(((((((&p_val, &a_val), &s_val), &d1_val), &d2_val), &d3_val), &b_val), &b3_val)| {
        lazily!(
          let mut res = vec![];
          res.extend(p_val.to_bytes_le().unwrap());
          res.extend(a_val.to_bytes_le().unwrap());
          res.extend(s_val.to_bytes_le().unwrap());
          res.extend(d1_val.to_bytes_le().unwrap());
          res.extend(d2_val.to_bytes_le().unwrap());
          res.extend(d3_val.to_bytes_le().unwrap());
          res.extend(b_val.to_bytes_le().unwrap());
          res.extend(b3_val.to_bytes_le().unwrap());
          res
        )
      },
    )
    .collect::<Vec<_>>();
  println!("Compute Merkle tree for the plain evaluations");
  // let start = std::time::Instant::now();

  let (_, m_root) =
    gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(&poly_evaluations_str, &[], &worker);
  // let mut m_tree: PermutedParallelMerkleTree<Vec<u8>, BlakeDigest> =
  //   PermutedParallelMerkleTree::new(&worker);
  // m_tree.update(poly_evaluations_str);
  // let m_root = m_tree.root().unwrap();
  // println!("m_root: {:?}", m_root);
  // let end: std::time::Duration = start.elapsed();
  // println!(
  //   "Computed Merkle root for the plain evaluations: {}.{:03}s",
  //   end.as_secs(),
  //   end.subsec_nanos() / 1_000_000
  // );
  println!("Computed Merkle root for the plain evaluations");

  // Based on the hashes of P, D and B, we select a random linear combination
  // of P * x^steps, P, B * x^steps, B and D, and prove the low-degreeness of that,
  // instead of proving the low-degreeness of P, B and D separately
  let k0 = T::one();
  let k1 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x01".to_vec()])).unwrap();
  let k2 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x02".to_vec()])).unwrap();
  let k3 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x03".to_vec()])).unwrap();
  let k4 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x04".to_vec()])).unwrap();
  let k5 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x05".to_vec()])).unwrap();
  let k6 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x06".to_vec()])).unwrap();
  let k7 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x07".to_vec()])).unwrap();
  let k8 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x08".to_vec()])).unwrap();
  let k9 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x09".to_vec()])).unwrap();
  let k10 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x0a".to_vec()])).unwrap();

  // Compute the linear combination. We don't even both calculating it in
  // coefficient form; we just compute the evaluations
  // let g2_to_the_steps = xs[steps];
  let g2_to_the_steps = g2.pow_vartime(&parse_bytes_to_u64_vec(&steps.to_le_bytes()));
  let mut powers = vec![T::one()];
  for _ in 1..precision {
    powers.push(g2_to_the_steps * powers.last().unwrap());
  }

  let mut l_evaluations: Vec<T> = vec![];
  for (
    (((((((&p_of_x, &a_of_x), &s_of_x), &d1_of_x), &d2_of_x), &d3_of_x), &b_of_x), &b3_of_x),
    &x_to_the_steps,
  ) in p_evaluations
    .iter()
    .zip(&a_evaluations)
    .zip(&s_evaluations)
    .zip(&d1_evaluations)
    .zip(&d2_evaluations)
    .zip(&d3_evaluations)
    .zip(&b2_evaluations)
    .zip(&b3_evaluations)
    .zip(&powers)
    .take(precision)
  {
    l_evaluations.push(
      k0 * d1_of_x
        + k1 * d2_of_x
        + k2 * d3_of_x
        + k3 * p_of_x
        + k4 * p_of_x * x_to_the_steps
        + k5 * b_of_x
        + k6 * b_of_x * x_to_the_steps
        + k7 * b3_of_x // TODO
        + k8 * b3_of_x * x_to_the_steps // TODO
        + k9 * a_of_x
        + k10 * s_of_x,
    );
  }

  let l_evaluations_str = l_evaluations
    .iter()
    .map(|x| lazily!(x.to_bytes_le().unwrap()))
    .collect::<Vec<_>>();

  println!("Compute Merkle tree for the linear combination of evaluations");
  // let start = std::time::Instant::now();
  let (_, l_root) =
    gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(&l_evaluations_str, &[], &worker);
  // let mut l_m_tree: PermutedParallelMerkleTree<Vec<u8>, BlakeDigest> =
  //   PermutedParallelMerkleTree::new(&worker);
  // l_m_tree.update(l_evaluations_str);
  // let l_root = l_m_tree.root().unwrap();
  // println!("l_m_root: {:?}", l_root);
  // let end: std::time::Duration = start.elapsed();
  // println!(
  //   "Computed Merkle root for the linear combination of evaluations: {}.{:03}s",
  //   end.as_secs(),
  //   end.subsec_nanos() / 1_000_000
  // );
  println!("Computed Merkle root for the linear combination of evaluations");

  // Do some spot checks of the Merkle tree at pseudo-random coordinates, excluding
  // multiples of `skips`
  let positions = get_pseudorandom_indices(
    l_root.as_ref(),
    precision as u32,
    SPOT_CHECK_SECURITY_FACTOR,
    skips as u32,
  )
  .iter()
  .map(|&i| i as usize)
  .collect::<Vec<usize>>();
  // println!("{:?}", positions);

  // let start = std::time::Instant::now();
  let (linear_comb_branches, _) =
    gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(&l_evaluations_str, &positions, &worker);
  // let linear_comb_branches = mk_multi_branch(&l_m_tree, &positions);
  // let end: std::time::Duration = start.elapsed();
  // println!(
  //   "Computed Merkle root for the linear combination of evaluations: {}.{:03}s",
  //   end.as_secs(),
  //   end.subsec_nanos() / 1_000_000
  // );
  println!("Computed Merkle root for the linear combination of evaluations");

  let mut augmented_positions = vec![];
  for &j in positions.iter().peekable() {
    augmented_positions.extend([
      j,
      (j + precision - skips) % precision,
      (j + original_steps / 3 * skips) % precision,
      (j + original_steps / 3 * 2 * skips) % precision,
    ]);
  }

  // let start = std::time::Instant::now();
  let (main_branches, _) = gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(
    &poly_evaluations_str,
    &augmented_positions,
    &worker,
  );
  // let end: std::time::Duration = start.elapsed();
  // println!(
  //   "Computed Merkle root for the plain evaluations: {}.{:03}s",
  //   end.as_secs(),
  //   end.subsec_nanos() / 1_000_000
  // );
  println!("Computed Merkle root for the plain evaluations");

  // let main_branches = mk_multi_branch(&m_tree, &augmented_positions);
  println!("Computed {} spot checks", SPOT_CHECK_SECURITY_FACTOR);

  // Return the Merkle roots of P and D, the spot check Merkle proofs,
  // and low-degree proofs of P and D
  let fri_proof = prove_low_degree::<T>(&l_evaluations, g2, precision / 4, skips as u32);

  StarkProof {
    m_root,
    l_root,
    a_root,
    main_branches,
    linear_comb_branches,
    fri_proof,
  }
  // println!("STARK computed in %.4f sec" % (time.time() - start_time))
}
