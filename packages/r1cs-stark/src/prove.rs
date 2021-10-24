use crate::utils::*;
use commitment::delayed::Delayed;
use commitment::lazily;
use commitment::merkle_tree::{gen_multi_proofs_multi_core, BlakeDigest};
use commitment::multicore::Worker;
use ff::PrimeField;
use ff_utils::ff_utils::{FromBytes, ToBytes};
use fri::fft::{best_fft, expand_root_of_unity, inv_best_fft};
use fri::fri::prove_low_degree;
use fri::poly_utils::{eval_poly_at, multi_inv};
use fri::utils::{get_pseudorandom_indices, parse_bytes_to_u64_vec};
use num::bigint::BigUint;
use std::io::Error;

pub fn mk_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
  witness_trace: &[T],
  computational_trace: &[T],
  public_wires: &[T],
  public_first_indices: &[(usize, usize)],
  permuted_indices: &[usize],
  coefficients: &[T],
  flag0: &[T],
  flag1: &[T],
  flag2: &[T],
  n_constraints: usize,
  n_wires: usize,
) -> Result<StarkProof, Error> {
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

  let log_steps = log2_ceil(original_steps - 1);
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
  let log_max_precision = calc_max_log_precision::<T>();
  println!("max_precision: 2^{:?}", log_max_precision);
  assert!(precision <= 2usize.pow(log_max_precision));

  let mut permuted_indices = permuted_indices.to_vec();
  permuted_indices.extend(original_steps..steps);
  // println!("permuted_indices.len(): {:?}", permuted_indices.len());

  let mut coefficients = coefficients.to_vec();
  coefficients.extend(vec![T::zero(); steps - original_steps]);

  let mut witness_trace = witness_trace.to_vec();
  witness_trace.extend(vec![T::zero(); steps - original_steps]);
  // println!("witness_trace.len(): {:?}", witness_trace.len());

  let mut computational_trace = computational_trace.to_vec();
  computational_trace.extend(vec![T::zero(); steps - original_steps]);
  // println!("computational_trace: {:?}", computational_trace);

  let ff_order = T::zero() - T::one();
  let times_nmr = BigUint::from_bytes_le(&ff_order.to_bytes_le().unwrap());
  let times_dnm = BigUint::from_bytes_le(&precision.to_le_bytes());
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  {
    let times = parse_bytes_to_u64_vec(&times_nmr.to_bytes_le()); // modulus - 1
    let unity = T::multiplicative_generator().pow_vartime(&times);
    debug_assert_eq!(unity, T::one()); // g0^(modulus - 1) = 1 mod modulus
  }

  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision = 1 mod modulus
  let xs = expand_root_of_unity(g2); // Powers of the higher-order root of unity
  let skips = precision / steps; // EXTENSION_FACTOR
  let g1 = xs[skips]; // root of unity x such that x^steps = 1
  let log_order_of_g1 = log_steps as u32;
  let log_order_of_g2 = log_precision as u32;

  // Interpolate the computational trace into a polynomial P, with each step
  // along a successive power of g1
  println!("calculate expanding polynomials");
  let worker = Worker::new();

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
  println!("Converted flags into a polynomial and low-degree extended it");

  let s_polynomial = inv_best_fft(witness_trace.clone(), &g1, &worker, log_order_of_g1); // S(X)
  let s_evaluations = best_fft(s_polynomial, &g2, &worker, log_order_of_g2);
  // println!("s_evaluations: {:?}", s_evaluations);
  println!("Converted witness trace into a polynomial and low-degree extended it");

  let p_polynomial = inv_best_fft(computational_trace, &g1, &worker, log_order_of_g1); // P(X)
  let p_evaluations = best_fft(p_polynomial, &g2, &worker, log_order_of_g2);
  // println!("p_evaluations: {:?}", p_evaluations);
  println!("Converted computational trace into a polynomial and low-degree extended it");

  let z_polynomial = calc_z_polynomial(steps);
  let z_evaluations = best_fft(z_polynomial, &g2, &worker, log_order_of_g2);
  // println!("z_evaluations: {:?}", z_evaluations);
  println!("Computed Z polynomial");

  let q1_evaluations = calc_q1_evaluations(
    &s_evaluations,
    &k_evaluations,
    &p_evaluations,
    &f0_evaluations,
    &f1_evaluations,
    precision,
    skips,
  );
  // println!("q1_evaluations: {:?}", q1_evaluations);
  println!("Computed Q1 polynomial");

  let q2_evaluations = calc_q2_evaluations(
    &p_evaluations,
    &f2_evaluations,
    precision,
    skips,
    original_steps,
  );
  // println!("f0_evaluations: {:?}", f0_evaluations);
  // println!("f2_evaluations: {:?}", f2_evaluations);
  // println!("q2_evaluations: {:?}", q2_evaluations);
  println!("Computed Q2 polynomial");

  // let converted_indices = (0..steps)
  //   .map(|v| T::from_bytes_le(v.to_le_bytes().as_ref()).unwrap())
  //   .collect::<Vec<_>>();
  let converted_indices = convert_usize_iter_to_ff_vec(0..steps);
  let index_polynomial = inv_best_fft(converted_indices, &g1, &worker, log_order_of_g1);
  let ext_indices = best_fft(index_polynomial, &g2, &worker, log_order_of_g2);
  println!("Computed extended index polynomial");

  let converted_permuted_indices = convert_usize_iter_to_ff_vec(permuted_indices.clone());
  let permuted_polynomial = inv_best_fft(converted_permuted_indices, &g1, &worker, log_order_of_g1);
  let ext_permuted_indices = best_fft(permuted_polynomial, &g2, &worker, log_order_of_g2);
  // println!("ext_permuted_indices: {:?}", ext_permuted_indices);
  println!("Computed extended permuted index polynomial");

  let a_root: BlakeDigest = get_accumulator_tree_root(&permuted_indices, &witness_trace, &worker);
  let r: Vec<T> = get_random_ff_values(a_root.as_ref(), precision as u32, 3, 0);
  // println!("r: {:?}", r);

  let a_mini_evaluations = calc_a_mini_evaluations(
    &witness_trace,
    &ext_indices,
    &ext_permuted_indices,
    &r,
    steps,
    skips,
  );
  let a_polynomials = inv_best_fft(a_mini_evaluations, &g1, &worker, log_order_of_g1);
  let a_evaluations = best_fft(a_polynomials, &g2, &worker, log_order_of_g2);
  println!("Computed A polynomial");

  // let z3_polynomial = calc_z_polynomial(steps);
  // let z3_evaluations = best_fft(z3_polynomial, &g2, &worker, log_order_of_g2);
  // println!("z3_evaluations: {:?}", z3_evaluations);
  // println!("Computed Z3 polynomial");

  let q3_evaluations = calc_q3_evaluations(
    &s_evaluations,
    &a_evaluations,
    &ext_indices,
    &ext_permuted_indices,
    &r,
    precision,
    skips,
  );
  println!("Computed Q3 polynomial");

  let inv_z_evaluations = multi_inv(&z_evaluations);
  // let inv_z2_evaluations = multi_inv(&z2_evaluations);
  // let inv_z3_evaluations = multi_inv(&z3_evaluations);

  let d1_evaluations = calc_d1_polynomial(&q1_evaluations, &inv_z_evaluations);
  println!("Computed D1 polynomial");

  let d2_evaluations = calc_d2_polynomial(&q2_evaluations, &inv_z_evaluations);
  println!("Computed D2 polynomial");

  let d3_evaluations: Vec<T> = calc_d3_polynomial(&q3_evaluations, &inv_z_evaluations);
  println!("Computed D3 polynomial");

  let interpolant2 = calc_i2_polynomial(public_first_indices, &xs, &public_wires, skips);
  let i2_evaluations: Vec<T> = xs.iter().map(|&x| eval_poly_at(&interpolant2, x)).collect();

  let interpolant3 = calc_i3_polynomial(&xs, skips);
  let i3_evaluations: Vec<T> = xs.iter().map(|&x| eval_poly_at(&interpolant3, x)).collect();
  println!("Computed I polynomial");

  let zb2_evaluations = calc_zb2_evaluations(&public_first_indices, &xs, precision, skips);
  let zb3_evaluations = calc_zb3_evaluations(&xs, precision, skips);
  println!("Computed Zb polynomial");

  let inv_zb2_evaluations: Vec<T> = multi_inv(&zb2_evaluations);
  let b2_evaluations = calc_b2_evaluations(&s_evaluations, &i2_evaluations, &inv_zb2_evaluations);

  let inv_zb3_evaluations: Vec<T> = multi_inv(&zb3_evaluations);
  let b3_evaluations = calc_b3_evaluations(&a_evaluations, &i3_evaluations, &inv_zb3_evaluations);
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

  let (_, m_root) =
    gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(&poly_evaluations_str, &[], &worker);
  println!("Computed Merkle root for the plain evaluations");

  // Based on the hashes of P, D and B, we select a random linear combination
  // of P * x^steps, P, B * x^steps, B and D, and prove that it is low-degree polynomial,
  // instead of proving that P, B and D are low degree polynomials separately.
  // let k0 = T::one();
  // let k1 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x01".to_vec()])).unwrap();
  // ...
  // let k10 = T::from_str(&mk_seed(&[m_root.as_ref().to_vec(), b"\x0a".to_vec()])).unwrap();
  let mut k = vec![T::one()];
  for i in 1u8..11 {
    k.push(
      T::from_str(&mk_seed(&[
        m_root.as_ref().to_vec(),
        i.to_be_bytes().to_vec(),
      ]))
      .unwrap(),
    );
  }

  // Compute the linear combination. We don't even both calculating it in
  // coefficient form; we just compute the evaluations
  let g2_to_the_steps = xs[steps];
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
      k[0] * d1_of_x
        + k[1] * d2_of_x
        + k[2] * d3_of_x
        + k[3] * p_of_x
        + k[4] * p_of_x * x_to_the_steps
        + k[5] * b_of_x
        + k[6] * b_of_x * x_to_the_steps
        + k[7] * b3_of_x
        + k[8] * b3_of_x * x_to_the_steps
        + k[9] * a_of_x
        + k[10] * s_of_x,
    );
  }

  let l_evaluations_str = l_evaluations
    .iter()
    .map(|x| lazily!(x.to_bytes_le().unwrap()))
    .collect::<Vec<_>>();

  let (_, l_root) =
    gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(&l_evaluations_str, &[], &worker);
  println!("Computed Merkle root for the linear combination of evaluations");

  // Do some spot checks of the Merkle tree at pseudo-random coordinates,
  // excluding multiples of `skips`
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

  let (linear_comb_branches, _) =
    gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(&l_evaluations_str, &positions, &worker);
  println!("Computed Merkle branch for the linear combination of evaluations");

  let mut augmented_positions = vec![];
  for &j in positions.iter().peekable() {
    augmented_positions.extend([
      j,
      (j + precision - skips) % precision,
      (j + original_steps / 3 * skips) % precision,
      (j + original_steps / 3 * 2 * skips) % precision,
    ]);
  }

  let (main_branches, _) = gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(
    &poly_evaluations_str,
    &augmented_positions,
    &worker,
  );
  println!("Computed Merkle branch for the plain evaluations");

  // Return the Merkle roots of P and D, the spot check Merkle proofs,
  // and low-degree proofs of P and D
  let fri_proof = prove_low_degree::<T>(&l_evaluations, g2, precision / 4, skips as u32);
  println!("Computed {} spot checks", SPOT_CHECK_SECURITY_FACTOR);

  Ok(StarkProof {
    m_root,
    l_root,
    a_root,
    main_branches,
    linear_comb_branches,
    fri_proof,
  })
}
