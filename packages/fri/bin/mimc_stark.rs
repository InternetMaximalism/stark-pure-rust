use ff::PrimeField;
use fri::ff_utils::{FromBytes, ToBytes};
use fri::fft::{expand_root_of_unity, fft, inv_fft};
use fri::fri::{prove_low_degree, verify_low_degree_proof, FriProof};
use fri::permuted_tree::{get_root, merklize, mk_multi_branch, verify_multi_branch};
use fri::poly_utils::{eval_poly_at, lagrange_interp, mul_polys, multi_inv};
use fri::utils::{blake, get_pseudorandom_indices, is_a_power_of_2, parse_bytes_to_u64_vec};
use num::bigint::BigUint;
use std::convert::TryInto;

fn parse_hex_to_decimal(value: &[u8]) -> String {
  BigUint::from_bytes_be(value).to_str_radix(10)
}

fn mk_seed(messages: &[Vec<u8>]) -> String {
  let mut message: Vec<u8> = vec![];
  for m in messages {
    message.extend(m);
  }
  parse_hex_to_decimal(&blake(&message))
}

pub fn mimc_computational_trace<T: PrimeField>(
  inp: T,
  steps: usize,
  round_constants: &[T],
) -> Vec<T> {
  let mut computational_trace = vec![inp];
  for i in 0..(steps - 1) {
    let last = *computational_trace.last().unwrap();
    let last_pow3 = last * last * last;
    computational_trace.push(last_pow3 + round_constants[i % round_constants.len()]);
  }
  computational_trace
}

// pub fn mimc<T: PrimeField>(inp: T, steps: usize, round_constants: &[T]) -> T {
//   let computational_trace = mimc_computational_trace(inp, steps, round_constants);
//   *computational_trace.last().unwrap()
// }

pub struct StarkProof {
  m_root: Vec<u8>,
  l_root: Vec<u8>,
  main_branches: Vec<Vec<Vec<u8>>>,
  linear_comb_branches: Vec<Vec<Vec<u8>>>,
  fri_proof: Vec<FriProof>,
}

// const modulus = 2**256 - 2**32 * 351 + 1;
// const non_residue = 7;
const EXTENSION_FACTOR: usize = 8usize;
const SPOT_CHECK_SECURITY_FACTOR: usize = 80usize;

// Generate a STARK for a MIMC calculation
pub fn mk_mimc_proof<T: PrimeField + FromBytes + ToBytes>(
  inp: T,
  steps: usize,
  round_constants: &[T],
) -> StarkProof {
  // start_time = time.time()
  // Some constraints to make our job easier
  assert!(steps <= 2usize.pow(32)); // EXTENSION_FACTOR
  assert!(is_a_power_of_2(steps));
  assert!(is_a_power_of_2(round_constants.len()));
  assert!(round_constants.len() < steps);

  let precision = steps * EXTENSION_FACTOR;

  // Root of unity such that x^precision=1
  let times_nmr = BigUint::from_bytes_be(&(T::zero() - T::one()).to_bytes_be().unwrap());
  let times_dnm = BigUint::from_bytes_be(&precision.to_be_bytes());
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus

  // Root of unity such that x^steps=1
  let skips = precision / steps;
  let g1 = g2.pow_vartime(&parse_bytes_to_u64_vec(&skips.to_le_bytes()));

  // Powers of the higher-order root of unity
  let xs = expand_root_of_unity(g2);
  let last_step_position = xs[(steps - 1) * EXTENSION_FACTOR];

  // Generate the computational trace
  let computational_trace = mimc_computational_trace(inp, steps, round_constants);
  let output = *computational_trace.last().unwrap();
  println!("Done generating computational trace");

  // Interpolate the computational trace into a polynomial P, with each step
  // along a successive power of g1
  let computational_trace_polynomial = inv_fft(&computational_trace, g1);
  let p_evaluations = fft(&computational_trace_polynomial, g2);
  println!("Converted computational steps into a polynomial and low-degree extended it");

  let skips2 = steps / round_constants.len();
  let constants_mini_polynomial = inv_fft(
    round_constants,
    g1.pow_vartime(&parse_bytes_to_u64_vec(&skips2.to_le_bytes())),
  );
  // constants_polynomial = [0 if i % skips2 else constants_mini_polynomial[i//skips2] for i in range(steps)]
  let constants_mini_extension = fft(
    &constants_mini_polynomial,
    g2.pow_vartime(&parse_bytes_to_u64_vec(&skips2.to_le_bytes())),
  );
  println!("Converted round constants into a polynomial and low-degree extended it");

  // Create the composed polynomial such that
  // C(P(x), P(g1*x), K(x)) = P(g1*x) - P(x)**3 - K(x)
  let mut c_of_p_evaluations = vec![];
  for i in 0..precision {
    let p_eval_pow3 = p_evaluations[i] * p_evaluations[i] * p_evaluations[i];
    c_of_p_evaluations.push(
      p_evaluations[(i + EXTENSION_FACTOR) % precision]
        - p_eval_pow3
        - constants_mini_extension[i % constants_mini_extension.len()],
    );
  }
  println!("Computed C(P, K) polynomial");

  // Compute D(x) = C(P(x), P(g1*x), K(x)) / Z(x)
  // Z(x) = (x^steps - 1) / (x - x_at_last_step)
  let z_num_evaluations: Vec<T> = (0..precision)
    .map(|i| xs[(i * steps) % precision] - T::one())
    .collect();
  let z_num_inv = multi_inv(&z_num_evaluations);
  let z_den_evaluations: Vec<T> = (0..precision).map(|i| xs[i] - last_step_position).collect();
  let d_evaluations: Vec<T> = c_of_p_evaluations
    .iter()
    .zip(z_den_evaluations)
    .zip(z_num_inv)
    .map(|((cp, zd), zni)| *cp * zd * zni)
    .collect();
  println!("Computed D polynomial");

  // Compute interpolant of ((1, input), (x_at_last_step, output))
  let interpolant = lagrange_interp(&[T::one(), last_step_position], &[inp, output]); // lagrange_interp_2
  let i_evaluations: Vec<T> = xs.iter().map(|x| eval_poly_at(&interpolant, *x)).collect();
  // OR
  // i_evaluations = fft(interpolant, modulus, g2)

  let zero_poly2 = mul_polys(
    &[T::zero() - T::one(), T::one()],
    &[T::zero() - last_step_position, T::one()],
  );
  let z2_evaluations: Vec<T> = xs.iter().map(|x| eval_poly_at(&zero_poly2, *x)).collect();
  let inv_z2_evaluations: Vec<T> = multi_inv(&z2_evaluations);

  let b_evaluations: Vec<T> = p_evaluations
    .iter()
    .zip(&i_evaluations)
    .zip(&inv_z2_evaluations)
    .map(|((&p, &i), &inv_z2)| (p - i) * inv_z2)
    .collect();
  println!("Computed B polynomial");

  // Compute their Merkle root
  let p_evaluations_str: Vec<Vec<u8>> = p_evaluations
    .iter()
    .zip(&d_evaluations)
    .zip(&b_evaluations)
    .map(|((&p_val, &d_val), &b_val)| {
      let mut res = vec![];
      res.extend(p_val.to_bytes_be().unwrap());
      res.extend(d_val.to_bytes_be().unwrap());
      res.extend(b_val.to_bytes_be().unwrap());
      res
    })
    .collect();
  let m_tree = merklize(&p_evaluations_str);
  let m_root = get_root(&m_tree);
  println!("Computed hash root");

  // Based on the hashes of P, D and B, we select a random linear combination
  // of P * x^steps, P, B * x^steps, B and D, and prove the low-degreeness of that,
  // instead of proving the low-degreeness of P, B and D separately
  let k1 = T::from_str(&mk_seed(&[m_root.clone(), b"\x01".to_vec()])).unwrap();
  let k2 = T::from_str(&mk_seed(&[m_root.clone(), b"\x02".to_vec()])).unwrap();
  let k3 = T::from_str(&mk_seed(&[m_root.clone(), b"\x03".to_vec()])).unwrap();
  let k4 = T::from_str(&mk_seed(&[m_root.clone(), b"\x04".to_vec()])).unwrap();

  // Compute the linear combination. We don"t even both calculating it in
  // coefficient form; we just compute the evaluations
  let g2_to_the_steps = xs[steps];
  let mut powers = vec![T::one()];
  for _ in 1..precision {
    powers.push(*powers.last().unwrap() * g2_to_the_steps);
  }

  let mut l_evaluations: Vec<T> = vec![];
  for (((&d_eval, &p_eval), &b_eval), &power) in d_evaluations
    .iter()
    .zip(&p_evaluations)
    .zip(&b_evaluations)
    .zip(&powers)
    .take(precision)
  {
    l_evaluations
      .push(d_eval + p_eval * k1 + p_eval * k2 * power + b_eval * k3 + b_eval * power * k4);
  }

  let l_evaluations_str: Vec<Vec<u8>> = l_evaluations
    .iter()
    .map(|x| x.to_bytes_be().unwrap())
    .collect();
  let l_m_tree = merklize(&l_evaluations_str);
  println!("Computed random linear combination");

  // Do some spot checks of the Merkle tree at pseudo-random coordinates, excluding
  // multiples of `EXTENSION_FACTOR`
  // let branches: Vec<T> = vec![];
  let positions = get_pseudorandom_indices(
    &get_root(&l_m_tree),
    precision.try_into().unwrap(),
    SPOT_CHECK_SECURITY_FACTOR,
    EXTENSION_FACTOR.try_into().unwrap(),
  );
  let mut augmented_positions = vec![];
  for &x in positions.iter().peekable() {
    augmented_positions.extend([x, (x + skips) % precision]);
  }
  // for pos in positions:
  //     branches.append(mk_branch(m_tree, pos))
  //     branches.append(mk_branch(m_tree, (pos + skips) % precision))
  //     branches.append(mk_branch(l_m_tree, pos))
  println!("Computed {} spot checks", SPOT_CHECK_SECURITY_FACTOR);

  // Return the Merkle roots of P and D, the spot check Merkle proofs,
  // and low-degree proofs of P and D
  let fri_proof = prove_low_degree::<T>(
    &l_evaluations,
    g2,
    steps * 2,
    EXTENSION_FACTOR.try_into().unwrap(),
  );

  StarkProof {
    m_root: get_root(&m_tree).clone(),
    l_root: get_root(&l_m_tree).clone(),
    main_branches: mk_multi_branch(&m_tree, &augmented_positions),
    linear_comb_branches: mk_multi_branch(&l_m_tree, &positions),
    fri_proof: fri_proof,
  }
  // println!("STARK computed in %.4f sec" % (time.time() - start_time))
}

// Verifies a STARK
pub fn verify_mimc_proof<T: PrimeField + FromBytes + ToBytes>(
  inp: T,
  steps: usize,
  round_constants: &[T],
  output: T,
  proof: StarkProof,
) -> Result<bool, &str> {
  let StarkProof {
    m_root,
    l_root,
    main_branches,
    linear_comb_branches,
    fri_proof,
  } = proof;

  // start_time = time.time()
  assert!(steps <= 2usize.pow(32)); // EXTENSION_FACTOR
  assert!(is_a_power_of_2(steps));
  assert!(is_a_power_of_2(round_constants.len()));
  assert!(round_constants.len() < steps);

  let precision = steps * EXTENSION_FACTOR;

  // Get (steps)th root of unity
  let times_nmr = BigUint::from_bytes_be(&(T::zero() - T::one()).to_bytes_be().unwrap());
  let times_dnm = BigUint::from_bytes_be(&precision.to_be_bytes());
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus
  let skips = precision / steps;

  // Gets the polynomial representing the round constants
  let skips2 = steps / round_constants.len();
  let constants_mini_polynomial = inv_fft(
    round_constants,
    g2.pow_vartime(&parse_bytes_to_u64_vec(
      &(EXTENSION_FACTOR * skips2).to_le_bytes(),
    )),
  );

  // Verifies the low-degree proofs
  assert!(verify_low_degree_proof(
    l_root.clone(),
    g2,
    &fri_proof,
    steps * 2,
    EXTENSION_FACTOR.try_into().unwrap()
  )
  .unwrap());

  // Performs the spot checks
  let k1 = T::from_str(&mk_seed(&[m_root.clone(), b"\x01".to_vec()])).unwrap();
  let k2 = T::from_str(&mk_seed(&[m_root.clone(), b"\x02".to_vec()])).unwrap();
  let k3 = T::from_str(&mk_seed(&[m_root.clone(), b"\x03".to_vec()])).unwrap();
  let k4 = T::from_str(&mk_seed(&[m_root.clone(), b"\x04".to_vec()])).unwrap();
  let positions = get_pseudorandom_indices(
    &l_root.clone(),
    precision.try_into().unwrap(),
    SPOT_CHECK_SECURITY_FACTOR,
    EXTENSION_FACTOR.try_into().unwrap(),
  );
  let mut augmented_positions = vec![];
  for &x in positions.iter().peekable() {
    augmented_positions.extend([x, (x + skips) % precision]);
  }
  let last_step_position = g2.pow_vartime(&parse_bytes_to_u64_vec(
    &((steps - 1) * skips).to_le_bytes(),
  ));

  let main_branch_leaves = verify_multi_branch(&m_root, &augmented_positions, &main_branches);
  let linear_comb_branch_leaves = verify_multi_branch(&l_root, &positions, &linear_comb_branches);

  for i in 0..positions.len() {
    let pos = positions[i];
    let x = g2.pow_vartime(&parse_bytes_to_u64_vec(&pos.to_le_bytes()));
    let x_to_the_steps = x.pow_vartime(&parse_bytes_to_u64_vec(&steps.to_le_bytes()));
    let mut m_branch1 = main_branch_leaves[i * 2].clone();
    let mut m_branch2 = main_branch_leaves[i * 2 + 1].clone();
    let l_of_x = T::from_bytes_be(linear_comb_branch_leaves[i].clone()).unwrap();

    let mut m_branch3 = m_branch1.split_off(32); // m_branch1 = leaves[i * 2][..32]
    let m_branch4 = m_branch3.split_off(32); // m_branch3 = leaves[i * 2][32..64], m_branch4 = leaves[i * 2][64..]
    let _ = m_branch2.split_off(32); // m_branch2 = main_branch_leaves[i * 2 + 1][..32]
    let p_of_x = T::from_bytes_be(m_branch1).unwrap();
    let p_of_g1x = T::from_bytes_be(m_branch2).unwrap();
    let d_of_x = T::from_bytes_be(m_branch3).unwrap();
    let b_of_x = T::from_bytes_be(m_branch4).unwrap();

    let z_value = (x.pow_vartime(&parse_bytes_to_u64_vec(&steps.to_le_bytes())) - T::one())
      * (x - last_step_position).invert().unwrap();
    let k_of_x = eval_poly_at(
      &constants_mini_polynomial,
      x.pow_vartime(&parse_bytes_to_u64_vec(&skips2.to_le_bytes())),
    );

    // Check transition constraints C(P(x)) = Z(x) * D(x)
    assert_eq!(
      p_of_g1x - p_of_x * p_of_x * p_of_x - k_of_x - z_value * d_of_x,
      T::zero()
    );

    // Check boundary constraints B(x) * Q(x) + I(x) = P(x)
    let interpolant = lagrange_interp(&[T::one(), last_step_position], &[inp, output]); // lagrange_interp_2
    let zero_poly2 = mul_polys(
      &[T::zero() - T::one(), T::one()],
      &[T::zero() - last_step_position, T::one()],
    );
    assert_eq!(
      p_of_x - b_of_x * eval_poly_at(&zero_poly2, x) - eval_poly_at(&interpolant, x),
      T::zero()
    );

    // Check correctness of the linear combination
    assert_eq!(
      l_of_x
        - d_of_x
        - k1 * p_of_x
        - k2 * p_of_x * x_to_the_steps
        - k3 * b_of_x
        - k4 * b_of_x * x_to_the_steps,
      T::zero()
    );
  }

  println!("Verified {} consistency checks", SPOT_CHECK_SECURITY_FACTOR);
  // println!("Verified STARK in %.4f sec" % (time.time() - start_time));
  Ok(true)
}

fn main() {
  // use crate::{mk_mimc_proof, verify_mimc_proof};
  use stark_pure_rust::fp::Fp;
  use stark_pure_rust::fri::fri_proof_bin_length;
  use stark_pure_rust::merkle_tree::bin_length;

  type TargetFiniteField = Fp;

  let inp = TargetFiniteField::from(3u64);
  let log_steps: u32 = match std::env::args().nth(1) {
    Some(x) => x
      .parse()
      .expect("The first command line argument was not a u32 value."),
    _ => 13,
  };
  let steps = 2usize.pow(log_steps);
  // Full STARK test
  // constants = [random.randrange(modulus) for i in range(64)]
  let constants: Vec<TargetFiniteField> = (0u64..64)
    .map(|i| TargetFiniteField::from(i.pow(7) ^ 42))
    .collect();
  let computational_trace = mimc_computational_trace(inp, steps, &constants);
  let proof = mk_mimc_proof(inp, steps, &constants);

  let StarkProof {
    m_root: _,
    l_root: _,
    main_branches,
    linear_comb_branches,
    fri_proof,
  } = &proof;
  let len1 = bin_length(main_branches) + bin_length(linear_comb_branches);
  let len2 = fri_proof_bin_length(fri_proof);
  println!(
    "Approx proof length: {} (branches), {} (FRI proof), {} (total)",
    len1,
    len2,
    len1 + len2
  );

  let answer = *computational_trace.last().unwrap(); // mimc(inp, steps, &constants)
  assert!(verify_mimc_proof(inp, steps, &constants, answer, proof).unwrap());
  println!("Done proof verification");
}
