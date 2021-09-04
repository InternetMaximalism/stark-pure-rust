use std::collections::HashMap;

use crate::utils::*;
use ff::PrimeField;
use fri::ff_utils::{FromBytes, ToBytes};
use fri::fft::{expand_root_of_unity, fft, inv_fft};
use fri::fri::prove_low_degree;
use fri::permuted_tree::{get_root, merklize, mk_multi_branch};
use fri::poly_utils::{div_polys, eval_poly_at, lagrange_interp, multi_inv, sparse};
use fri::utils::{get_pseudorandom_indices, parse_bytes_to_u64_vec};
use num::bigint::BigUint;

pub fn mk_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
  witness_trace: &[T],
  computational_trace: &[T],
  public_wires: &[T],
  public_first_indices: &[(usize, usize)],
  permuted_indices: &[usize],
  last_coeff_list: &[usize],
  coefficients: &[T], // This argument may be good to use HashMap<usize, T> instead of Vec<T> because we can omit zero coefficients from it.
  n_constraints: usize,
  n_wires: usize,
) -> StarkProof {
  let original_steps = coefficients.len();
  assert!(original_steps <= 3 * n_constraints * n_wires);
  assert!(original_steps % 3 == 0);
  assert_eq!(witness_trace.len(), original_steps);
  assert_eq!(computational_trace.len(), original_steps);
  assert_eq!(coefficients.len(), original_steps);
  // assert_eq!(3 * (last_coeff_list.len() - 1), original_steps);

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
  let mut permuted_indices = permuted_indices.to_vec();
  permuted_indices.extend(original_steps..steps);

  let mut witness_trace = witness_trace.to_vec();
  witness_trace.extend(vec![T::zero(); steps - original_steps]);

  let mut permuted_witness_trace = vec![];
  for j in 0..permuted_indices.len() {
    permuted_witness_trace.push(witness_trace[permuted_indices[j]].clone());
  }

  let mut computational_trace = computational_trace.to_vec();
  computational_trace.extend(vec![T::zero(); steps - original_steps]);

  let mut coefficients = coefficients.to_vec();
  coefficients.extend(vec![T::zero(); steps - original_steps]);

  // Root of unity such that x^precision=1
  let times_nmr = BigUint::from_bytes_be(&(T::zero() - T::one()).to_bytes_be().unwrap());
  let times_dnm = BigUint::from_bytes_be(&precision.to_be_bytes());
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus
  let xs = expand_root_of_unity(g2); // Powers of the higher-order root of unity
  let skips = precision / steps; // EXTENSION_FACTOR
  let g1 = xs[skips]; // Root of unity such that x^steps=1

  // Interpolate the computational trace into a polynomial P, with each step
  // along a successive power of g1

  let s_polynomial = inv_fft(&witness_trace, g1); // S(X)
  let s_evaluations = fft(&s_polynomial, g2);
  println!("Converted witness trace into a polynomial and low-degree extended it");

  let p_polynomial = inv_fft(&computational_trace, g1); // P(X)
  let p_evaluations = fft(&p_polynomial, g2);
  // println!("{:?}", computational_trace);
  // println!("{:?}", p_polynomial);
  // println!("{:?}", p_evaluations);
  println!("Converted computational trace into a polynomial and low-degree extended it");

  let permuted_s_polynomial = inv_fft(&permuted_witness_trace, g1); // P(X)
  let permuted_s_evaluations = fft(&permuted_s_polynomial, g2);
  println!("Converted permuted computational trace into a polynomial and low-degree extended it");

  // println!("{:?}", coefficients);
  let k_polynomial = inv_fft(&coefficients, g1); // K(X)
                                                 // println!("{:?}", k_polynomial);
  let k_evaluations = fft(&k_polynomial, g2);
  // println!("{:?}", k_evaluations);
  println!("Converted coefficients into a polynomial and low-degree extended it");

  let mut ext_last_coeff_list: Vec<usize> = vec![];
  ext_last_coeff_list.append(&mut last_coeff_list.to_vec());
  ext_last_coeff_list.append(
    &mut last_coeff_list
      .iter()
      .map(|i| i + original_steps / 3)
      .collect(),
  );
  ext_last_coeff_list.append(
    &mut last_coeff_list
      .iter()
      .map(|i| i + original_steps / 3 * 2)
      .collect(),
  );
  println!("ext_last_coeff_list.len(): {:?}", ext_last_coeff_list.len());

  // Create the composed polynomial such that
  // Q(g1^j) = P(g1^(j-1)) + P(g1^(j % n_constraints))*K(g1^j) - P(g1^j)
  // let z1_nmr_evaluations: Vec<T> = (0..precision)
  //   .map(|i| xs[(i * steps) % precision] - T::one())
  //   .collect();
  // let z1_nmr_inv = multi_inv(&z1_nmr_evaluations);
  // let mut z1_dnm_evaluations: Vec<T> = vec![T::one(); precision];
  let mut z1_evaluations: Vec<T> = vec![T::one(); precision];
  let mut q1_evaluations = vec![];
  for j in 0..precision {
    let s_of_x = s_evaluations[j]; // S(g1^j)
    let k_of_x = k_evaluations[j]; // K(g1^j)
    let p_of_prev_x = p_evaluations[(j + precision - skips) % precision]; // P(g1^(j-1))
    let p_of_x = p_evaluations[j]; // P(g1^j)

    q1_evaluations.push(p_of_x - p_of_prev_x - k_of_x * s_of_x);

    // if j == 50 {
    //   println!(
    //     "{:?} {:?}    {:?} {:?}    {:?}",
    //     p_of_x,
    //     p_of_prev_x_plus_half,
    //     p_of_x_plus_half,
    //     k_of_x_plus_half,
    //     p_of_x_plus_half - p_of_prev_x_plus_half - k_of_x_plus_half * p_of_x
    //   );
    // }

    // if (p_of_x_plus_half - p_of_prev_x_plus_half - k_of_x_plus_half * p_of_x == T::zero()) {
    //   println!("valid {}", j);
    // }

    if j < original_steps * skips && j % skips == 0 && !ext_last_coeff_list.contains(&(j / skips)) {
      z1_evaluations = z1_evaluations
        .iter()
        .enumerate()
        .map(|(i, &val)| val * (xs[i] - xs[j]))
        .collect();
    }
  }

  println!("Computed Q1 polynomial");
  // println!("{:?}", q1_evaluations);
  // println!("{:?}", z1_evaluations);
  // let z1_polynomial = inv_fft(&z1_evaluations, g2);
  // println!("{:?}", z1_polynomial);

  // let mut z2_dnm_evaluations: Vec<T> = vec![T::one(); precision];
  let mut z2_evaluations: Vec<T> = vec![T::one(); precision];
  let mut q2_evaluations = vec![];
  for j in 0..precision {
    let j1 = j;
    let j2 = (j1 + last_coeff_list.len() * skips) % precision;
    let j3 = (j2 + last_coeff_list.len() * skips) % precision;
    let a_eval = p_evaluations[j1];
    let b_eval = p_evaluations[j2];
    let c_eval = p_evaluations[j3];
    q2_evaluations.push(c_eval - a_eval * b_eval);

    // if j % skips == 0
    // // c_eval == a_eval * b_eval
    // {
    //   println!("{:03} {:?} {:?}    {:?}", j, a_eval, b_eval, c_eval);
    // }

    if j < original_steps * skips && j % skips == 0 && last_coeff_list.contains(&(j / skips)) {
      // println!("{:?} {:?}", j, c_eval - a_eval * b_eval);
      z2_evaluations = z2_evaluations
        .iter()
        .enumerate()
        .map(|(i, &val)| val * (xs[i] - xs[j]))
        .collect();
    }
  }

  println!("Computed Q2 polynomial");

  // let mut z2_dnm_evaluations: Vec<T> = vec![T::one(); precision];

  let mut a_nmr_evaluations: Vec<T> = vec![];
  let mut a_dnm_evaluations: Vec<T> = vec![];
  let mut val_nmr_list = vec![];
  let mut val_dnm_list = vec![];
  let r1 = T::zero(); // TODO: randomize
  let r2 = T::one(); // TODO: randomize
  let r3 = T::multiplicative_generator(); // TODO: randomize
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
    let val_nmr = r1 + r2 * T::from(j as u64) + r3 * witness_trace[j];
    let val_dnm = r1 + r2 * T::from(j as u64) + r3 * witness_trace[permuted_indices[j]];
    let acc_nmr = val_nmr * last_acc_nmr;
    let acc_dnm = val_dnm * last_acc_dnm;

    val_nmr_list.push(val_nmr);
    val_dnm_list.push(val_dnm);
    a_nmr_evaluations.push(acc_nmr);
    a_dnm_evaluations.push(acc_dnm);
  }

  let inv_a_dnm_evaluations = multi_inv(&a_dnm_evaluations);
  let a_mini_evaluations: Vec<T> = a_nmr_evaluations
    .iter()
    .zip(&inv_a_dnm_evaluations)
    .map(|(&a_nmr, &a_dnm)| a_nmr * a_dnm)
    .collect();
  let a_polynomials = inv_fft(&a_mini_evaluations, g1);
  let a_evaluations = fft(&a_polynomials, g2);

  let inv_val_dnm_list = multi_inv(&val_dnm_list);
  let v_mini_evaluations: Vec<T> = val_nmr_list
    .iter()
    .zip(&inv_val_dnm_list)
    .map(|(&v_nmr, &v_dnm)| v_nmr * v_dnm)
    .collect();
  let v_polynomials = inv_fft(&v_mini_evaluations, g1);
  let v_evaluations = fft(&v_polynomials, g2);

  let mut q3_evaluations = vec![];
  for j in 0..precision {
    // A(j) * val_dnm = A(j - 1) * val_nmr
    q3_evaluations.push(
      a_evaluations[j] - a_evaluations[(j + precision - steps) % precision] * v_evaluations[j],
    );
  }

  let mut sparse_z3_nmr = HashMap::new();
  sparse_z3_nmr.insert(0, -T::one());
  sparse_z3_nmr.insert(steps, T::one());
  let z3_dnm = [-T::one(), T::one()];
  let z3_polynomial = div_polys(&sparse(sparse_z3_nmr), &z3_dnm);
  let z3_evaluations = fft(&z3_polynomial, g2);

  // Compute D(x) = Q(x) / Z(x)
  // let z_nmr_evaluations: Vec<T> = (0..precision)
  //   .map(|i| xs[(i * steps) % precision] - T::one())
  //   .collect();
  // let z_nmr_evaluations: Vec<T> = (0..precision)
  //   .map(|i| xs[(i * steps) % precision] - T::one())
  //   .collect();
  // let inv_z_nmr_evaluations: Vec<T> = multi_inv(&z_nmr_evaluations);
  let inv_z1_evaluations = multi_inv(&z1_evaluations);
  let inv_z2_evaluations = multi_inv(&z2_evaluations);
  let inv_z3_evaluations = multi_inv(&z3_evaluations);
  // let inv_z1_evaluations: Vec<T> = z1_dnm_evaluations
  //   .iter()
  //   .zip(&inv_z_nmr_evaluations)
  //   .map(|(&z1d, &zni)| z1d * zni)
  //   .collect();
  // let z1_evaluations: Vec<T> = inv_z1_dnm_evaluations
  //   .iter()
  //   .zip(&z_nmr_evaluations)
  //   .map(|(&z1d, &zni)| z1d * zni)
  //   .collect();
  // for (i, (&q1, &z1)) in q1_evaluations.iter().zip(&z1_evaluations).enumerate() {
  //   if z1 == T::zero() {
  //     if q1 != T::zero() {
  //       assert_eq!(i, 0);
  //     } else {
  //       println!("valid {:?}", i);
  //     }
  //   }
  // }
  // let z1_polynomial = inv_fft(&z1_evaluations, g2);
  // let q1_polynomial = reduction_poly(&inv_fft(&q1_evaluations, g2), precision);
  // let d1_polynomial = div_polys(&q1_polynomial, &z1_polynomial);
  // let q1_polynomial_copy = reduction_poly(&mul_polys(&z1_polynomial, &d1_polynomial), precision);
  // assert_eq!(q1_polynomial, q1_polynomial_copy);
  // println!("{:?}", d1_polynomial);

  println!("Computed Q3 polynomial");

  // let d1_evaluations = fft(&d1_polynomial, g2);
  // println!("{:?}", d1_evaluations);
  let d1_evaluations: Vec<T> = q1_evaluations
    .iter()
    .zip(&inv_z1_evaluations)
    .map(|(&q1, &z1i)| q1 * z1i)
    .collect();

  // println!("P1: {:?}", p_evaluations[59]);
  // println!("P1: {:?}", p_evaluations[51]);
  // println!("K1: {:?}", k_evaluations[59]);
  // println!("S: {:?}", s_evaluations[59]);
  // println!("Q1: {:?}", q1_evaluations[59]);
  // println!("Z1: {:?}", z1_evaluations[59]);
  // println!("D1: {:?}", d1_evaluations[59]);
  // println!("Z1_dnm: {:?}", z1_dnm_evaluations[393]);
  // println!("Z1_nmr^-1: {:?}", inv_z_nmr_evaluations[393]);
  // println!(
  //   "Z1^-1: {:?}",
  //   z1_dnm_evaluations[393] * inv_z_nmr_evaluations[393]
  // );
  // println!("Z1: {:?}", z1_evaluations[393]);
  // println!("D1*Z1: {:?}", d1_evaluations[393] * z1_evaluations[393]);
  // println!(
  //   "Q1/Z1: {:?}",
  //   q1_evaluations[393] * z1_evaluations[393].invert().unwrap()
  // );

  let d2_evaluations: Vec<T> = q2_evaluations
    .iter()
    .zip(&inv_z2_evaluations)
    .map(|(&q2, &z2i)| q2 * z2i)
    .collect();

  let d3_evaluations: Vec<T> = q3_evaluations
    .iter()
    .zip(&inv_z3_evaluations)
    .map(|(&q3, &z3i)| q3 * z3i)
    .collect();
  // println!("deg Q2: {:?}", q2_evaluations[510]);
  // println!("deg D2: {:?}", d2_evaluations[510]);
  // println!("Z2_dnm: {:?}", z2_dnm_evaluations[510]);
  // println!("Z2_nmr^-1: {:?}", inv_z_nmr_evaluations[510]);
  // println!(
  //   "Z2^-1: {:?}",
  //   z2_dnm_evaluations[510] * inv_z_nmr_evaluations[510]
  // );
  // println!(
  //   "Z2: {:?}",
  //   (z2_dnm_evaluations[510] * inv_z_nmr_evaluations[510])
  //     .invert()
  //     .unwrap()
  // );
  println!("Computed D polynomial");

  // let interpolant = {
  //   let mut x_vals = vec![];
  //   let mut y_vals = vec![];
  //   // println!("witness_trace: {:?}", witness_trace);
  //   // println!("coefficients: {:?}", coefficients);
  //   // println!("last_coeff_list: {:?}", last_coeff_list);
  //   // println!("public_wires: {:?}", public_wires);
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
  // i_evaluations = fft(interpolant, modulus, g2)

  let interpolant2 = {
    let mut x_vals: Vec<T> = vec![];
    let mut y_vals: Vec<T> = vec![];
    for (k, w) in public_first_indices {
      x_vals.push(xs[steps * w]);
      y_vals.push(public_wires[*k]);
    }
    lagrange_interp(&x_vals, &y_vals)
  };

  let i2_evaluations: Vec<T> = xs.iter().map(|&x| eval_poly_at(&interpolant2, x)).collect();

  let interpolant3 = {
    let x_vals = vec![T::one(), xs[(original_steps - 1) * skips]];
    let y_vals = vec![
      T::one(), // val_nmr_list[0] * val_dnm_list[0].invert().unwrap(),
      T::one(),
    ];
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
    .map(|(i, &val)| val * (xs[i] - T::one()) * (xs[i] - xs[(original_steps - 1) * skips]))
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
  let b3_evaluations: Vec<T> = a_evaluations
    .iter()
    .zip(&i3_evaluations)
    .zip(&inv_zb3_evaluations)
    .map(|((&a, &i), &inv_zb)| (a - i) * inv_zb)
    .collect();
  println!("Computed B polynomial");

  // Compute their Merkle root
  let poly_evaluations_str: Vec<Vec<u8>> = p_evaluations
    .iter()
    .zip(&a_evaluations)
    .zip(&s_evaluations)
    .zip(&permuted_s_evaluations)
    .zip(&d1_evaluations)
    .zip(&d2_evaluations)
    .zip(&d3_evaluations)
    .zip(&b2_evaluations)
    .zip(&b3_evaluations)
    .map(
      |(
        (((((((&p_val, &a_val), &s_val), &permuted_s_val), &d1_val), &d2_val), &d3_val), &b_val),
        &b3_val,
      )| {
        let mut res = vec![];
        res.extend(p_val.to_bytes_be().unwrap());
        res.extend(a_val.to_bytes_be().unwrap());
        res.extend(s_val.to_bytes_be().unwrap());
        res.extend(permuted_s_val.to_bytes_be().unwrap());
        res.extend(d1_val.to_bytes_be().unwrap());
        res.extend(d2_val.to_bytes_be().unwrap());
        res.extend(d3_val.to_bytes_be().unwrap());
        res.extend(b_val.to_bytes_be().unwrap());
        res.extend(b3_val.to_bytes_be().unwrap());
        res
      },
    )
    .collect();
  println!("Compute Merkle tree");

  let m_tree = merklize(&poly_evaluations_str);
  let m_root = get_root(&m_tree);
  println!("Computed hash root");

  // Based on the hashes of P, D and B, we select a random linear combination
  // of P * x^steps, P, B * x^steps, B and D, and prove the low-degreeness of that,
  // instead of proving the low-degreeness of P, B and D separately
  let k0 = T::one();
  let k1 = T::from_str(&mk_seed(&[m_root.clone(), b"\x01".to_vec()])).unwrap();
  let k2 = T::from_str(&mk_seed(&[m_root.clone(), b"\x02".to_vec()])).unwrap();
  let k3 = T::from_str(&mk_seed(&[m_root.clone(), b"\x03".to_vec()])).unwrap();
  let k4 = T::from_str(&mk_seed(&[m_root.clone(), b"\x04".to_vec()])).unwrap();
  let k5 = T::from_str(&mk_seed(&[m_root.clone(), b"\x05".to_vec()])).unwrap();
  let k6 = T::from_str(&mk_seed(&[m_root.clone(), b"\x06".to_vec()])).unwrap();
  let k7 = T::from_str(&mk_seed(&[m_root.clone(), b"\x07".to_vec()])).unwrap();
  let k8 = T::from_str(&mk_seed(&[m_root.clone(), b"\x08".to_vec()])).unwrap();
  let k9 = T::from_str(&mk_seed(&[m_root.clone(), b"\x09".to_vec()])).unwrap();
  let k10 = T::from_str(&mk_seed(&[m_root.clone(), b"\x0a".to_vec()])).unwrap();
  let k11 = T::from_str(&mk_seed(&[m_root.clone(), b"\x0b".to_vec()])).unwrap();

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
    (
      (
        ((((((&p_of_x, &a_of_x), &s_of_x), permuted_s_of_x), &d1_of_x), &d2_of_x), &d3_of_x),
        &b_of_x,
      ),
      &b3_of_x,
    ),
    &x_to_the_steps,
  ) in p_evaluations
    .iter()
    .zip(&a_evaluations)
    .zip(&s_evaluations)
    .zip(&permuted_s_evaluations)
    .zip(&d1_evaluations)
    .zip(&d2_evaluations)
    .zip(&d3_evaluations)
    .zip(&b2_evaluations)
    .zip(&b3_evaluations)
    .zip(&powers)
    .take(precision)
  {
    // if l_evaluations.len() == 50 || l_evaluations.len() == 59 {
    //   println!(
    //     "{:?}  {:?}  {:?}  {:?}  {:?}  {:?}  {:?}",
    //     l_evaluations.len(),
    //     p_of_x,
    //     s_of_x,
    //     permuted_s_of_x,
    //     d1_of_x,
    //     d2_of_x,
    //     d3_of_x
    //   );
    // }
    // println!("{:?}  {:?}  {:?}", l_evaluations.len(), p_of_x, s_of_x,);
    l_evaluations.push(
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
        + k10 * s_of_x
        + k11 * permuted_s_of_x,
    );
  }

  let l_evaluations_str: Vec<Vec<u8>> = l_evaluations
    .iter()
    .map(|x| x.to_bytes_be().unwrap())
    .collect();
  let l_m_tree = merklize(&l_evaluations_str);
  println!("Computed random linear combination");

  // Do some spot checks of the Merkle tree at pseudo-random coordinates, excluding
  // multiples of `skips`
  let positions = get_pseudorandom_indices(
    &get_root(&l_m_tree),
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
    //   (j + last_coeff_list.len() * skips) % precision,
    //   (j + 2 * last_coeff_list.len() * skips) % precision,
    // );
    augmented_positions.extend([
      j,
      (j + precision - skips) % precision,
      (j + last_coeff_list.len() * skips) % precision,
      (j + 2 * last_coeff_list.len() * skips) % precision,
    ]);
  }
  // let branches: Vec<T> = vec![];
  // for pos in positions:
  //     branches.append(mk_branch(m_tree, pos))
  //     branches.append(mk_branch(m_tree, (pos + skips) % precision))
  //     branches.append(mk_branch(l_m_tree, pos))
  println!("Computed {} spot checks", SPOT_CHECK_SECURITY_FACTOR);

  // Return the Merkle roots of P and D, the spot check Merkle proofs,
  // and low-degree proofs of P and D
  let fri_proof = prove_low_degree::<T>(&l_evaluations, g2, precision / 4, skips as u32);

  // println!("{:?}", positions);

  StarkProof {
    m_root: get_root(&m_tree).clone(),
    l_root: get_root(&l_m_tree).clone(),
    main_branches: mk_multi_branch(&m_tree, &augmented_positions),
    linear_comb_branches: mk_multi_branch(&l_m_tree, &positions),
    fri_proof: fri_proof,
  }
  // println!("STARK computed in %.4f sec" % (time.time() - start_time))
}
