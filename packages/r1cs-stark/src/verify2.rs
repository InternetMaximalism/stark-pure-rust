use std::collections::HashMap;

use crate::utils::*;
use ff::PrimeField;
use ff_utils::ff_utils::{FromBytes, ToBytes};
use fri::fft::{best_fft, expand_root_of_unity, inv_best_fft};
use fri::fri::verify_low_degree_proof;
use fri::merkle_tree2::verify_multi_branch;
use fri::multicore::Worker;
use fri::poly_utils::{eval_poly_at, lagrange_interp, sparse};
use fri::utils::{get_pseudorandom_indices, parse_bytes_to_u64_vec};
use num::bigint::BigUint;

pub fn verify_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
  proof: StarkProof,
  public_wires: &[T],
  public_first_indices: &[(usize, usize)],
  permuted_indices: &[usize],
  coefficients: &[T],
  flag0: &[T],
  flag1: &[T],
  flag2: &[T],
  n_constraints: usize,
  n_wires: usize,
) -> Result<bool, String> {
  let original_steps = coefficients.len();
  println!("original_steps: {:?}", original_steps);
  assert!(original_steps <= 3 * n_constraints * n_wires);
  assert!(original_steps % 3 == 0);

  let mut log_steps = 1;
  let mut tmp_steps = original_steps - 1;
  while tmp_steps > 1 {
    tmp_steps /= 2;
    log_steps += 1;
  }
  let mut steps = 2usize.pow(log_steps);
  assert!(steps <= 2usize.pow(32));
  if steps < 8 {
    steps = 8;
  }

  let precision = steps * EXTENSION_FACTOR;
  let log_precision = log_steps + LOG_EXTENSION_FACTOR as u32;

  let mut permuted_indices = permuted_indices.to_vec();
  permuted_indices.extend(original_steps..steps);

  let mut coefficients = coefficients.to_vec();
  coefficients.extend(vec![T::zero(); steps - original_steps]);

  // start_time = time.time()
  let StarkProof {
    m_root,
    l_root,
    main_branches,
    linear_comb_branches,
    fri_proof,
  } = proof;

  // Get (steps)th root of unity
  let times_nmr = BigUint::from_bytes_le(&(T::zero() - T::one()).to_bytes_le().unwrap());
  let times_dnm = BigUint::from_bytes_le(&precision.to_le_bytes());
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus
  let xs = expand_root_of_unity(g2);
  let skips = precision / steps; // EXTENSION_FACTOR
  let g1 = xs[skips];
  let log_order_of_g1 = log_steps as u32;
  let log_order_of_g2 = log_precision as u32;

  let worker = Worker::new();

  // K(X)
  // println!("coeff: {:?}", coefficients);
  let k_polynomial = inv_best_fft(coefficients, &g1, &worker, log_order_of_g1);
  // let k_evaluations = best_fft(k_polynomial, &g2, &worker, log_order_of_g2);
  // println!("k: {:?}", k_evaluations);
  println!("Converted coefficients into a polynomial and low-degree extended it");

  let f0_polynomial = inv_best_fft(flag0.to_vec(), &g1, &worker, log_order_of_g1);
  // let mut f_evaluations = f_polynomial;
  // if f_evaluations.len() < precision {
  //   let mut padding = vec![T::zero(); (precision - f_evaluations.len()) as usize];
  //   f_evaluations.append(&mut padding);
  // }
  // let f_evaluations = best_fft(f_polynomial, &g2, &worker, log_order_of_g2);
  let f1_polynomial = inv_best_fft(flag1.to_vec(), &g1, &worker, log_order_of_g1);
  let f2_polynomial = inv_best_fft(flag2.to_vec(), &g1, &worker, log_order_of_g1);

  // Verifies the low-degree proofs
  assert!(
    verify_low_degree_proof(l_root.clone(), g2, &fri_proof, precision / 4, skips as u32).unwrap()
  );

  // Performs the spot checks
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
  let positions = get_pseudorandom_indices(
    l_root.as_ref(),
    precision as u32,
    SPOT_CHECK_SECURITY_FACTOR,
    skips as u32,
  );
  let mut augmented_positions = vec![];
  for &j in positions.iter().peekable() {
    // println!(
    //   "{:?} {:?} {:?} {:?}",
    //   j,
    //   (j + precision - skips) % precision,
    //   (j + original_steps / 3 * skips) % precision,
    //   (j + 2 * original_steps / 3 * skips) % precision,
    // );
    augmented_positions.extend([
      j,
      (j + precision - skips) % precision,
      (j + original_steps / 3 * skips) % precision,
      (j + 2 * original_steps / 3 * skips) % precision,
    ]);
  }

  // println!("{:?}", positions);
  let main_branches = main_branches;
  let main_branch_leaves =
    verify_multi_branch(&m_root, &augmented_positions, main_branches).unwrap();
  let linear_comb_branch_leaves =
    verify_multi_branch(&l_root, &positions, linear_comb_branches).unwrap();

  let z1_polynomial = {
    let mut sparse_z1 = HashMap::new();
    sparse_z1.insert(0, -T::one());
    sparse_z1.insert(steps, T::one());
    sparse(sparse_z1)
  };
  let z1_evaluations = best_fft(z1_polynomial, &g2, &worker, log_order_of_g2);

  // let mut z2_evaluations: Vec<T> = vec![T::one(); precision];
  // // for k in 0..n_constraints {
  // //   let j = ((3 * k + 2) * n_wires - 1) * skips;
  // //   z2_evaluations = z2_evaluations
  // //     .iter()
  // //     .enumerate()
  // //     .map(|(i, &val)| val * (xs[i] - xs[j]))
  // //     .collect();
  // // }

  // for j in 0..precision {
  //   if j < original_steps * skips && j % skips == 0 && last_coeff_list.contains(&(j / skips)) {
  //     z2_evaluations = z2_evaluations
  //       .iter()
  //       .enumerate()
  //       .map(|(i, &val)| val * (xs[i] - xs[j]))
  //       .collect();
  //   }
  // }

  let z3_polynomial = {
    let mut sparse_z3_nmr = HashMap::new();
    sparse_z3_nmr.insert(0, -T::one());
    sparse_z3_nmr.insert(steps, T::one());
    sparse(sparse_z3_nmr)
  };
  let z3_evaluations = best_fft(z3_polynomial, &g2, &worker, log_order_of_g2);

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

  // let interpolant = {
  //   let mut x_vals = vec![];
  //   let mut y_vals = vec![];
  //   for (j, n_coeff) in n_coeff_list[0..public_wires.len()]
  //     .to_vec()
  //     .iter()
  //     .enumerate()
  //   {
  //     x_vals.push(xs[n_coeff * skips]);
  //     y_vals.push(public_wires[j]);
  //   }

  //   lagrange_interp(&x_vals, &y_vals)
  // };

  let interpolant2 = {
    let mut x_vals: Vec<T> = vec![];
    let mut y_vals: Vec<T> = vec![];
    for (k, w) in public_first_indices {
      x_vals.push(xs[skips * w]);
      y_vals.push(public_wires[*k]);
    }
    lagrange_interp(&x_vals, &y_vals)
  };

  let x_of_last_step = xs[(steps - 1) * skips];
  let interpolant3 = {
    let x_vals = vec![x_of_last_step];
    let y_vals = vec![T::one()];
    lagrange_interp(&x_vals, &y_vals)
  };
  println!("Computed boundary polynomial");

  let r1 = T::zero(); // TODO: randomize
  let r2 = T::one(); // TODO: randomize
  let r3 = T::multiplicative_generator(); // TODO: randomize
  for (i, &pos) in positions.iter().enumerate() {
    let x = xs[pos]; // g2.pow_vartime(&parse_bytes_to_u64_vec(&pos.to_le_bytes()));

    let mut m_branch0 = main_branch_leaves[i * 4].clone();
    let mut m_branch10 = main_branch_leaves[i * 4 + 1].clone();
    let mut m_branch20 = main_branch_leaves[i * 4 + 2].clone();
    let mut m_branch30 = main_branch_leaves[i * 4 + 3].clone();
    let l_of_x = T::from_bytes_le(&linear_comb_branch_leaves[i]).unwrap();

    let mut m_branch1 = m_branch0.split_off(32); // m_branch0 = leaves[i * 5][..32]
    let mut m_branch2 = m_branch1.split_off(32); // m_branch1 = leaves[i * 5][32..64]
    let mut m_branch3 = m_branch2.split_off(32); // m_branch2 = leaves[i * 5][64..96]
    let mut m_branch4 = m_branch3.split_off(32); // m_branch3 = leaves[i * 5][96..128]
    let mut m_branch5 = m_branch4.split_off(32); // m_branch4 = leaves[i * 5][128..160]
    let mut m_branch6 = m_branch5.split_off(32); // m_branch5 = leaves[i * 5][160..192]
    let mut m_branch7 = m_branch6.split_off(32); // m_branch6 = leaves[i * 5][192..224]
    let _ = m_branch7.split_off(32); // m_branch7 = leaves[i * 5][224..]
    let mut m_branch11 = m_branch10.split_off(32); // m_branch10 = leaves[i * 5 + 1][..32]
    let _ = m_branch11.split_off(32); // m_branch11 = leaves[i * 5 + 1][32..64]
    let _ = m_branch20.split_off(32); // m_branch20 = leaves[i * 5 + 2][..32]
    let _ = m_branch30.split_off(32); // m_branch30 = leaves[i * 5 + 3][..32]
    let p_of_x = T::from_bytes_le(&m_branch0).unwrap();
    let p_of_prev_x = T::from_bytes_le(&m_branch10).unwrap();
    let a_of_x = T::from_bytes_le(&m_branch1).unwrap();
    let a_of_prev_x = T::from_bytes_le(&m_branch11).unwrap();
    let s_of_x = T::from_bytes_le(&m_branch2).unwrap();
    let d1_of_x = T::from_bytes_le(&m_branch3).unwrap();
    let d2_of_x = T::from_bytes_le(&m_branch4).unwrap();
    let d3_of_x = T::from_bytes_le(&m_branch5).unwrap();
    let b_of_x = T::from_bytes_le(&m_branch6).unwrap();
    let b3_of_x = T::from_bytes_le(&m_branch7).unwrap();
    let p_of_x_plus_w = T::from_bytes_le(&m_branch20).unwrap();
    let p_of_x_plus_2w = T::from_bytes_le(&m_branch30).unwrap();

    let z1_value = z1_evaluations[pos];
    // let z2_value = z2_evaluations[pos];
    let z3_value = z3_evaluations[pos];

    let k_of_x = eval_poly_at(&k_polynomial, x);
    let f0 = eval_poly_at(&f0_polynomial, x);
    let f1 = eval_poly_at(&f1_polynomial, x);
    let f2 = eval_poly_at(&f2_polynomial, x);

    // Check first transition constraints Q1(x) = Z1(x) * D1(x)
    assert_eq!(
      f0 * (p_of_x - f1 * p_of_prev_x - k_of_x * s_of_x),
      z1_value * d1_of_x
    );

    // Check second transition constraints Q2(x) = Z2(x) * D2(x)
    assert_eq!(
      f2 * (p_of_x_plus_2w - p_of_x * p_of_x_plus_w),
      // p_of_x_plus_2w - p_of_x * p_of_x_plus_w,
      z1_value * d2_of_x
    );

    let val_nmr = r1 + r2 * ext_indices[pos] + r3 * s_of_x;
    let val_dnm = r1 + r2 * ext_permuted_indices[pos] + r3 * s_of_x;

    // Check third transition constraints Q3(x) = Z3(x) * D3(x)
    assert_eq!(a_of_x * val_dnm - a_of_prev_x * val_nmr, z3_value * d3_of_x);

    // Check boundary constraints P(x) - I(x) = Zb(x) * B(x)

    // let mut zb_of_x = T::one();
    // for w in 0..n_wires {
    //   let j = w * skips;
    //   zb_of_x = zb_of_x * (x - xs[j]);
    // }
    // for k in 0..n_constraints {
    //   let j = n_coeff_list[k] * skips;
    //   zb_of_x = zb_of_x * (x - xs[j]);
    // }

    // assert_eq!(p_of_x - eval_poly_at(&interpolant, x), zb_of_x * b_of_x);

    // println!(
    //   "{:03} {:?} {:?}    {:?} {:?}",
    //   pos,
    //   p_of_x,
    //   eval_poly_at(&interpolant, x),
    //   zb_of_x,
    //   b_of_x
    // );

    let mut zb2_of_x = T::one();
    for (_, w) in public_first_indices {
      zb2_of_x *= x - xs[w * skips];
    }

    assert_eq!(s_of_x - eval_poly_at(&interpolant2, x), zb2_of_x * b_of_x);

    let zb3_of_x = x - x_of_last_step;
    assert_eq!(a_of_x - eval_poly_at(&interpolant3, x), zb3_of_x * b3_of_x);

    // Check correctness of the linear combination
    let x_to_the_steps = x.pow_vartime(&parse_bytes_to_u64_vec(&steps.to_le_bytes()));
    assert_eq!(
      l_of_x,
      k0 * d1_of_x
        + k1 * d2_of_x
        + k2 * d3_of_x
        + k3 * p_of_x
        + k4 * p_of_x * x_to_the_steps
        + k5 * b_of_x
        + k6 * b_of_x * x_to_the_steps
        + k7 * b3_of_x
        + k8 * b3_of_x * x_to_the_steps
        + k9 * a_of_x
        + k10 * s_of_x,
    );
  }

  println!("Verified {} consistency checks", SPOT_CHECK_SECURITY_FACTOR);
  // println!("Verified STARK in %.4f sec" % (time.time() - start_time));
  Ok(true)
}
