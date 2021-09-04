use crate::utils::*;
use ff::PrimeField;
use fri::ff_utils::{FromBytes, ToBytes};
use fri::fft::{expand_root_of_unity, fft, inv_fft};
use fri::fri::prove_low_degree;
use fri::permuted_tree::{get_root, merklize, mk_multi_branch};
use fri::poly_utils::{eval_poly_at, lagrange_interp, multi_inv};
use fri::utils::{get_pseudorandom_indices, parse_bytes_to_u64_vec};
use num::bigint::BigUint;

pub fn mk_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
  computational_trace: &[T],
  public_input: &[T],
  n_coeff_list: &[usize],
  coefficients: &[T], // This argument may be good to use HashMap<usize, T> instead of Vec<T> because we can omit zero coefficients from it.
  n_constraints: usize,
  n_wires: usize,
) -> StarkProof {
  let mut constants = coefficients.to_vec();
  constants.extend(coefficients.to_vec());
  let original_steps = constants.len();
  assert!(original_steps <= 6 * n_constraints * n_wires);
  assert_eq!(computational_trace.len(), original_steps);

  let mut log_steps = 1;
  let mut tmp_steps = original_steps - 1;
  while tmp_steps > 1 {
    tmp_steps /= 2;
    log_steps += 1;
  }
  let steps = 2usize.pow(log_steps);

  assert!(steps <= 2usize.pow(32));
  let precision = steps * EXTENSION_FACTOR;

  let mut computational_trace = computational_trace.to_vec();
  computational_trace.extend(vec![T::zero(); steps - original_steps]);

  let mut constants = constants.to_vec();
  constants.extend(vec![T::zero(); steps - original_steps]);

  // Root of unity such that x^precision=1
  let times_nmr = BigUint::from_bytes_be(&(T::zero() - T::one()).to_bytes_be().unwrap());
  let times_dnm = BigUint::from_bytes_be(&precision.to_be_bytes());
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus

  // Root of unity such that x^steps=1
  let skips = precision / steps; // EXTENSION_FACTOR

  // Powers of the higher-order root of unity
  let xs = expand_root_of_unity(g2);
  let g1 = xs[skips];

  // Interpolate the computational trace into a polynomial P, with each step
  // along a successive power of g1
  let p_polynomial = inv_fft(&computational_trace, g1);
  let p_evaluations = fft(&p_polynomial, g2);
  println!("Converted computational steps into a polynomial and low-degree extended it");

  let k_polynomial = inv_fft(&constants, g1);
  let k_evaluations = fft(&k_polynomial, g2);
  println!("Converted constants into a polynomial and low-degree extended it");

  // Create the composed polynomial such that
  // Q(g1^j) = P(g1^(j-1)) + P(g1^(j % n_constraints))*K(g1^j) - P(g1^j)
  // let z1_num_evaluations: Vec<T> = (0..precision)
  //   .map(|i| xs[(i * steps) % precision] - T::one())
  //   .collect();
  // let z1_num_inv = multi_inv(&z1_num_evaluations);
  // let mut z1_dnm_evaluations: Vec<T> = vec![T::one(); precision];
  let mut z1_evaluations: Vec<T> = vec![T::one(); precision];
  let mut q1_evaluations = vec![];
  for j in 0..precision {
    let half = 3 * n_constraints * n_wires * skips;
    let p_of_x_plus_half = p_evaluations[(j + half) % precision]; // P(g1^next_k)
    let p_of_prev_x_plus_half = p_evaluations[(j + half - skips) % precision]; // P(g1^(next_k-1))
    let p_of_x = p_evaluations[j]; // P(g1^(next_k % n_wires))
    let k_of_x_plus_half = k_evaluations[(j + half) % precision]; // K(g1^next_k)

    q1_evaluations.push(p_of_x_plus_half - p_of_prev_x_plus_half - k_of_x_plus_half * p_of_x);
    // if j == 393 {
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

    // if (j >= original_steps * skips / 2 || j % (n_wires * skips) == 0 || j % skips != 0)
    //   && (j % skips == 0)
    // {
    //   z1_dnm_evaluations = z1_dnm_evaluations
    //     .iter()
    //     .enumerate()
    //     .map(|(i, &val)| val * (xs[i] - xs[j]))
    //     .collect();
    // }

    if j < original_steps * skips / 2 && j % skips == 0 && j % (n_wires * skips) != 0 {
      z1_evaluations = z1_evaluations
        .iter()
        .enumerate()
        .map(|(i, &val)| val * (xs[i] - xs[j]))
        .collect();
    }
  }

  // let mut z2_dnm_evaluations: Vec<T> = vec![T::one(); precision];
  let mut z2_evaluations: Vec<T> = vec![T::one(); precision];
  let mut q2_evaluations = vec![];
  for j in 0..precision {
    let j1 = j;
    let j2 = (j1 + n_coeff_list.len() * skips) % precision;
    let j3 = (j2 + n_coeff_list.len() * skips) % precision;
    let a_eval = p_evaluations[j1];
    let b_eval = p_evaluations[j2];
    let c_eval = p_evaluations[j3];
    q2_evaluations.push(c_eval - a_eval * b_eval);

    // if j % skips == 0
    // // c_eval == a_eval * b_eval
    // {
    //   println!("{:03} {:?} {:?}    {:?}", j, a_eval, b_eval, c_eval);
    // }

    // if j < original_steps * skips
    //   && j >= (2 * n_wires - 1) * skips
    //   && (j - (2 * n_wires - 1) * skips) % (3 * n_wires * skips) == 0
    // if (j >= original_steps * skips / 2
    //   || j < (n_wires - 1) * skips
    //   || (j - (n_wires - 1) * skips) % (3 * n_wires * skips) != 0)
    //   && j % skips == 0
    // {
    //   z2_dnm_evaluations = z2_dnm_evaluations
    //     .iter()
    //     .enumerate()
    //     .map(|(i, &val)| val * (xs[i] - xs[j]))
    //     .collect();
    // }

    if j < original_steps * skips
      && j >= original_steps * skips / 2
      && (original_steps * skips / 2 - j) % skips == 0
      && (n_coeff_list.contains(&(&(original_steps * skips / 2 - j) / skips)))
    {
      z2_evaluations = z2_evaluations
        .iter()
        .enumerate()
        .map(|(i, &val)| val * (xs[i] - xs[j]))
        .collect();
    }
  }
  println!("Computed Q polynomial");

  // Compute D(x) = Q(x) / Z(x)
  // let z_nmr_evaluations: Vec<T> = (0..precision)
  //   .map(|i| xs[(i * steps) % precision] - T::one())
  //   .collect();
  // let z_nmr_evaluations: Vec<T> = (0..precision)
  //   .map(|i| xs[(i * steps) % precision] - T::one())
  //   .collect();
  // let inv_z_nmr_evaluations: Vec<T> = multi_inv(&z_nmr_evaluations);
  let inv_z1_evaluations: Vec<T> = multi_inv(&z1_evaluations);
  let inv_z2_evaluations: Vec<T> = multi_inv(&z2_evaluations);
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

  // let d1_evaluations = fft(&d1_polynomial, g2);
  // println!("{:?}", d1_evaluations);
  let d1_evaluations: Vec<T> = q1_evaluations
    .iter()
    .zip(&inv_z1_evaluations)
    .map(|(&q1, &z1i)| q1 * z1i)
    .collect();
  // println!("P1: {:?}", p_evaluations[(393 + 240 - 8) % 512]);
  // println!("P1: {:?}", p_evaluations[(393 + 240) % 512]);
  // println!("P1: {:?}", p_evaluations[393]);
  // println!("K1: {:?}", k_evaluations[(393 + 240) % 512]);
  // println!("Q1: {:?}", q1_evaluations[393]);
  // println!("D1: {:?}", d1_evaluations[393]);
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
  println!("Computed D polynomial");
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

  // {
  //   let ys = &d1_evaluations;
  //   // let ys: Vec<T> = inv_z1_dnm_evaluations
  //   //   .iter()
  //   //   .map(|&z1| {
  //   //     if z1 == T::zero() {
  //   //       T::zero()
  //   //     } else {
  //   //       z1.invert().unwrap()
  //   //     }
  //   //   })
  //   //   .collect();
  //   // let ys: Vec<T> = inv_z1_dnm_evaluations
  //   //   .iter()
  //   //   .zip(&z_nmr_evaluations)
  //   //   .map(|(&z1di, &z1n)| (z1di * z1n))
  //   //   .collect();
  //   // let ys = multi_inv(&inv_z1_dnm_evaluations);
  //   // let ys: Vec<T> = ys
  //   //   .iter()
  //   //   .zip(&q1_evaluations)
  //   //   .map(|(&z1i, &q1)| (z1i * q1))
  //   //   .collect();
  //   let mut pts: Vec<usize> = (0..ys.len())
  //     .filter(|&x| x % (skips as usize) != 0)
  //     .collect();
  //   for &pos in pts.iter().peekable() {
  //     assert_eq!(
  //       q1_evaluations[pos],
  //       d1_evaluations[pos] * z1_evaluations[pos]
  //     );
  //   }
  //   let rest = pts.split_off(7 * ys.len() / 8 - 1); // 3 * precision / 4 + 1); // pts[max_deg_plus_1..]
  //   let x_vals: Vec<T> = pts.iter().map(|&pos| xs[pos]).collect();
  //   let y_vals: Vec<T> = pts.iter().map(|&pos| ys[pos]).collect();
  //   let poly = lagrange_interp(&x_vals, &y_vals);

  //   // println!("{:?}", poly);
  //   // println!("{:?}", poly.len());
  //   for (_, &pos) in rest.iter().enumerate() {
  //     assert_eq!(eval_poly_at(&poly, xs[pos]), ys[pos]);
  //   }
  // }

  let interpolant = {
    let mut x_vals = vec![];
    let mut y_vals = vec![];
    for w in 0..n_wires {
      x_vals.push(xs[w * skips]);
      y_vals.push(public_input[w]);
    }

    for k in 0..n_constraints {
      x_vals.push(xs[n_coeff_list[k] * skips]);
      y_vals.push(public_input[n_coeff_list[k]]);
    }

    lagrange_interp(&x_vals, &y_vals)
  };
  let i_evaluations: Vec<T> = xs.iter().map(|&x| eval_poly_at(&interpolant, x)).collect();
  // OR
  // i_evaluations = fft(interpolant, modulus, g2)

  let mut zb_evaluations = vec![T::one(); precision];
  for w in 0..n_wires {
    let j = w * skips;
    zb_evaluations = zb_evaluations
      .iter()
      .enumerate()
      .map(|(i, &val)| val * (xs[i] - xs[j]))
      .collect();
  }
  for k in 0..n_constraints {
    let j = n_coeff_list[k] * skips;
    zb_evaluations = zb_evaluations
      .iter()
      .enumerate()
      .map(|(i, &val)| val * (xs[i] - xs[j]))
      .collect();
  }
  let inv_zb_evaluations: Vec<T> = multi_inv(&zb_evaluations);

  // B(x) = (P(x) - I(x)) / Z_b(x)
  let b_evaluations: Vec<T> = p_evaluations
    .iter()
    .zip(&i_evaluations)
    .zip(&inv_zb_evaluations)
    .map(|((&p, &i), &inv_zb)| (p - i) * inv_zb)
    .collect();
  println!("Computed B polynomial");

  // Compute their Merkle root
  let poly_evaluations_str: Vec<Vec<u8>> = p_evaluations
    .iter()
    .zip(&d1_evaluations)
    .zip(&d2_evaluations)
    .zip(&b_evaluations)
    .map(|(((&p_val, &d1_val), &d2_val), &b_val)| {
      let mut res = vec![];
      res.extend(p_val.to_bytes_be().unwrap());
      res.extend(d1_val.to_bytes_be().unwrap());
      res.extend(d2_val.to_bytes_be().unwrap());
      res.extend(b_val.to_bytes_be().unwrap());
      res
    })
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

  // Compute the linear combination. We don't even both calculating it in
  // coefficient form; we just compute the evaluations
  // let g2_to_the_steps = xs[steps];
  let g2_to_the_steps = g2.pow_vartime(&parse_bytes_to_u64_vec(&steps.to_le_bytes()));
  let mut powers = vec![T::one()];
  for _ in 1..precision {
    powers.push(g2_to_the_steps * powers.last().unwrap());
  }

  let mut l_evaluations: Vec<T> = vec![];
  for ((((&p_of_x, &d1_of_x), &d2_of_x), &b_of_x), &x_to_the_steps) in p_evaluations
    .iter()
    .zip(&d1_evaluations)
    .zip(&d2_evaluations)
    .zip(&b_evaluations)
    .zip(&powers)
    .take(precision)
  {
    l_evaluations.push(
      k0 * d1_of_x
        + k1 * d2_of_x
        + k2 * p_of_x
        + k3 * p_of_x * x_to_the_steps
        + k4 * b_of_x
        + k5 * b_of_x * x_to_the_steps,
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
    let half = 3 * n_constraints * n_wires * skips;
    augmented_positions.extend([
      j,
      (j + half - skips) % precision,
      (j + half) % precision,
      (j + n_wires * public_input.len() * skips) % precision,
      (j + 2 * n_wires * public_input.len() * skips) % precision,
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

  StarkProof {
    m_root: get_root(&m_tree).clone(),
    l_root: get_root(&l_m_tree).clone(),
    main_branches: mk_multi_branch(&m_tree, &augmented_positions),
    linear_comb_branches: mk_multi_branch(&l_m_tree, &positions),
    fri_proof: fri_proof,
  }
  // println!("STARK computed in %.4f sec" % (time.time() - start_time))
}
