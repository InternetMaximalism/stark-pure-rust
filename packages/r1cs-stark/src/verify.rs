// use crate::utils::*;
// use ff::PrimeField;
// use fri::ff_utils::{FromBytes, ToBytes};
// use fri::fft::{expand_root_of_unity, inv_fft};
// use fri::fri::verify_low_degree_proof;
// use fri::permuted_tree::verify_multi_branch;
// use fri::poly_utils::{eval_poly_at, lagrange_interp};
// use fri::utils::{get_pseudorandom_indices, parse_bytes_to_u64_vec};
// use num::bigint::BigUint;

// pub fn verify_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
//   proof: StarkProof,
//   public_input: &[T],
//   n_coeff_list: &[usize],
//   coefficients: &[T],
//   n_constraints: usize,
//   n_wires: usize,
// ) -> Result<bool, String> {
//   // assert_eq!(public_input.len(), n_wires);

//   let mut constants = coefficients.to_vec();
//   constants.extend(coefficients.to_vec());
//   let original_steps = constants.len();
//   assert!(original_steps <= 6 * n_constraints * n_wires);

//   let mut log_steps = 1;
//   let mut tmp_steps = original_steps - 1;
//   while tmp_steps > 1 {
//     tmp_steps /= 2;
//     log_steps += 1;
//   }
//   let steps = 2usize.pow(log_steps);

//   // start_time = time.time()
//   assert!(steps <= 2usize.pow(32));

//   let StarkProof {
//     m_root,
//     l_root,
//     main_branches,
//     linear_comb_branches,
//     fri_proof,
//   } = proof;

//   let precision = steps * EXTENSION_FACTOR;
//   let skips = precision / steps; // EXTENSION_FACTOR

//   // Get (steps)th root of unity
//   let times_nmr = BigUint::from_bytes_be(&(T::zero() - T::one()).to_bytes_be().unwrap());
//   let times_dnm = BigUint::from_bytes_be(&precision.to_be_bytes());
//   assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
//   let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
//   let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus
//   let xs = expand_root_of_unity(g2);
//   let g1 = xs[skips];

//   // Gets the polynomial representing the constants
//   let k_polynomial = inv_fft(&constants, g1);

//   // Verifies the low-degree proofs
//   assert!(
//     verify_low_degree_proof(l_root.clone(), g2, &fri_proof, precision / 4, skips as u32).unwrap()
//   );

//   // Performs the spot checks
//   let k0 = T::one();
//   let k1 = T::from_str(&mk_seed(&[m_root.clone(), b"\x01".to_vec()])).unwrap();
//   let k2 = T::from_str(&mk_seed(&[m_root.clone(), b"\x02".to_vec()])).unwrap();
//   let k3 = T::from_str(&mk_seed(&[m_root.clone(), b"\x03".to_vec()])).unwrap();
//   let k4 = T::from_str(&mk_seed(&[m_root.clone(), b"\x04".to_vec()])).unwrap();
//   let k5 = T::from_str(&mk_seed(&[m_root.clone(), b"\x05".to_vec()])).unwrap();
//   let positions = get_pseudorandom_indices(
//     &l_root.clone(),
//     precision as u32,
//     SPOT_CHECK_SECURITY_FACTOR,
//     skips as u32,
//   );
//   let mut augmented_positions = vec![];
//   for &j in positions.iter().peekable() {
//     let half = 3 * n_constraints * n_wires * skips;
//     augmented_positions.extend([
//       j,
//       (j + half - skips) % precision,
//       (j + half) % precision,
//       (j + n_wires * skips) % precision,
//       (j + 2 * n_wires * skips) % precision,
//     ]);
//   }

//   let main_branch_leaves = verify_multi_branch(&m_root, &augmented_positions, &main_branches);
//   let linear_comb_branch_leaves = verify_multi_branch(&l_root, &positions, &linear_comb_branches);

//   // let z1_num_evaluations: Vec<T> = (0..precision)
//   //   .map(|i| xs[(i * steps) % precision] - T::one())
//   //   .collect();
//   // let mut z1_den_evaluations: Vec<T> = vec![T::one(); precision];
//   // for j in 0..precision {
//   //   if j >= original_steps * skips || (j >= skips && (j - skips) % (n_wires * skips) == 0) {
//   //     z1_den_evaluations = z1_den_evaluations
//   //       .iter()
//   //       .enumerate()
//   //       .map(|(i, &val)| val * (xs[i] - xs[j]))
//   //       .collect();
//   //   }
//   // }
//   // let z1_den_inv = multi_inv(&z1_den_evaluations);
//   let mut z1_evaluations: Vec<T> = vec![T::one(); precision];
//   for j in 0..precision {
//     if j < original_steps * skips / 2 && j % skips == 0 && j % (n_wires * skips) != 0 {
//       z1_evaluations = z1_evaluations
//         .iter()
//         .enumerate()
//         .map(|(i, &val)| val * (xs[i] - xs[j]))
//         .collect();
//     }
//   }

//   let mut z2_evaluations: Vec<T> = vec![T::one(); precision];
//   // for k in 0..n_constraints {
//   //   let j = ((3 * k + 2) * n_wires - 1) * skips;
//   //   z2_evaluations = z2_evaluations
//   //     .iter()
//   //     .enumerate()
//   //     .map(|(i, &val)| val * (xs[i] - xs[j]))
//   //     .collect();
//   // }

//   for j in 0..precision {
//     if j < original_steps * skips
//       && j >= original_steps * skips / 2
//       && (original_steps * skips / 2 - j) % skips == 0
//       && (n_coeff_list.contains(&(&(original_steps * skips / 2 - j) / skips)))
//     {
//       z2_evaluations = z2_evaluations
//         .iter()
//         .enumerate()
//         .map(|(i, &val)| val * (xs[i] - xs[j]))
//         .collect();
//     }
//   }

//   for (i, &pos) in positions.iter().enumerate() {
//     let mut m_branch0 = main_branch_leaves[i * 5].clone();
//     let mut m_branch1 = main_branch_leaves[i * 5 + 1].clone();
//     let mut m_branch2 = main_branch_leaves[i * 5 + 2].clone();
//     let mut m_branch10 = main_branch_leaves[i * 5 + 3].clone();
//     let mut m_branch11 = main_branch_leaves[i * 5 + 4].clone();
//     let l_of_x = T::from_bytes_be(linear_comb_branch_leaves[i].clone()).unwrap();

//     let mut m_branch3 = m_branch0.split_off(32); // m_branch0 = leaves[i * 5][..32]
//     let mut m_branch4 = m_branch3.split_off(32); // m_branch3 = leaves[i * 5][32..64]
//     let mut m_branch5 = m_branch4.split_off(32); // m_branch4 = leaves[i * 5][64..96]
//     let _ = m_branch5.split_off(32); // m_branch5 = leaves[i * 5][96..]
//     let _ = m_branch1.split_off(32); // m_branch1 = leaves[i * 5 + 1][..32]
//     let _ = m_branch2.split_off(32); // m_branch2 = leaves[i * 5 + 2][..32]
//     let _ = m_branch10.split_off(32); // m_branch10 = leaves[i * 5 + 3][..32]
//     let _ = m_branch11.split_off(32); // m_branch11 = leaves[i * 5 + 4][..32]
//     let p_of_x = T::from_bytes_be(m_branch0).unwrap();
//     let p_of_prev_x_plus_half = T::from_bytes_be(m_branch1).unwrap();
//     let p_of_x_plus_half = T::from_bytes_be(m_branch2).unwrap();
//     let d1_of_x = T::from_bytes_be(m_branch3).unwrap();
//     let d2_of_x = T::from_bytes_be(m_branch4).unwrap();
//     let b_of_x = T::from_bytes_be(m_branch5).unwrap();
//     let p_of_x_plus_w = T::from_bytes_be(m_branch10).unwrap();
//     let p_of_x_plus_2w = T::from_bytes_be(m_branch11).unwrap();

//     let z1_value = z1_evaluations[pos];
//     let z2_value = z2_evaluations[pos];

//     let x = xs[pos]; // g2.pow_vartime(&parse_bytes_to_u64_vec(&pos.to_le_bytes()));
//     let k_of_x_plus_half = eval_poly_at(
//       &k_polynomial,
//       xs[(pos + 3 * n_constraints * n_wires * skips) % precision],
//     );
//     // Check first transition constraints Q1(x) = Z1(x) * D1(x)
//     // println!(
//     //   "{:03} {:?} {:?}    {:?} {:?}    {:?} {:?}",
//     //   pos, p_of_x, p_of_prev_x_plus_half, p_of_x_plus_half, k_of_x_plus_half, z1_value, d1_of_x
//     // );
//     assert_eq!(
//       p_of_x_plus_half - p_of_prev_x_plus_half - k_of_x_plus_half * p_of_x,
//       z1_value * d1_of_x
//     );

//     // Check second transition constraints Q2(x) = Z2(x) * D2(x)
//     // println!(
//     //   "{:03} {:?} {:?}    {:?} {:?}    {:?} {:?}",
//     //   pos,
//     //   p_of_x,
//     //   p_of_x_plus_w,
//     //   p_of_x_plus_2w,
//     //   z2_value * d2_of_x,
//     //   z2_value,
//     //   d2_of_x
//     // );
//     assert_eq!(p_of_x_plus_2w - p_of_x * p_of_x_plus_w, z2_value * d2_of_x);

//     // Check boundary constraints P(x) - I(x) = Zb(x) * B(x)
//     let interpolant = {
//       let mut x_vals = vec![];
//       let mut y_vals = vec![];
//       for w in 0..n_wires {
//         x_vals.push(xs[n_coeff_list[w] * skips]);
//         y_vals.push(public_input[n_coeff_list[w]]);
//       }

//       for k in 0..n_constraints {
//         x_vals.push(xs[n_coeff_list[k] * skips]);
//         y_vals.push(public_input[n_coeff_list[k]]);
//       }

//       lagrange_interp(&x_vals, &y_vals)
//     };

//     let mut zb_of_x = T::one();
//     for w in 0..n_wires {
//       let j = w * skips;
//       zb_of_x = zb_of_x * (x - xs[j]);
//     }
//     for k in 0..n_constraints {
//       let j = n_coeff_list[k] * skips;
//       zb_of_x = zb_of_x * (x - xs[j]);
//     }
//     // println!(
//     //   "{:03} {:?} {:?}    {:?} {:?}",
//     //   pos,
//     //   p_of_x,
//     //   eval_poly_at(&interpolant, x),
//     //   zb_of_x,
//     //   b_of_x
//     // );
//     assert_eq!(p_of_x - eval_poly_at(&interpolant, x), zb_of_x * b_of_x);

//     // Check correctness of the linear combination
//     let x_to_the_steps = x.pow_vartime(&parse_bytes_to_u64_vec(&steps.to_le_bytes()));
//     assert_eq!(
//       l_of_x,
//       k0 * d1_of_x
//         + k1 * d2_of_x
//         + k2 * p_of_x
//         + k3 * p_of_x * x_to_the_steps
//         + k4 * b_of_x
//         + k5 * b_of_x * x_to_the_steps,
//     );
//   }

//   println!("Verified {} consistency checks", SPOT_CHECK_SECURITY_FACTOR);
//   // println!("Verified STARK in %.4f sec" % (time.time() - start_time));
//   Ok(true)
// }
