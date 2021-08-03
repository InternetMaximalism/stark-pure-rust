use ff::PrimeField;
use num::bigint::BigUint;
use stark_pure_rust::ff_utils::{FromBytes, ToBytes};
use stark_pure_rust::fft::{expand_root_of_unity, fft, inv_fft};
use stark_pure_rust::fri::{prove_low_degree, verify_low_degree_proof, FriProof};
use stark_pure_rust::permuted_tree::{get_root, merklize, mk_multi_branch, verify_multi_branch};
use stark_pure_rust::poly_utils::{eval_poly_at, lagrange_interp, multi_inv};
use stark_pure_rust::utils::{blake, get_pseudorandom_indices, parse_bytes_to_u64_vec};

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

#[derive(Debug, PartialEq)]
pub struct Coefficient {
  pub wire_id: u32,
  pub value: [u32; 8],
}

#[derive(Debug, PartialEq)]
pub struct Factor {
  pub n_coefficient: u32,
  pub coefficients: Vec<Coefficient>,
}

#[derive(Debug, PartialEq)]
pub struct Constraint {
  pub factors: Vec<Factor>,
}

#[derive(Debug, PartialEq)]
pub struct Constraints(pub Vec<Constraint>);

#[derive(Debug, PartialEq)]
pub struct Version(pub u32);

#[derive(Debug, PartialEq)]
pub struct Header {
  pub field_size: u32,
  pub prime_number: [u128; 2],
  pub n_wires: u32,
  pub n_public_outputs: u32,
  pub n_public_inputs: u32,
  pub n_private_inputs: u32,
  pub n_labels: u64,
  pub n_constraints: u32,
}

pub struct R1csContents {
  pub version: Version,
  pub header: Header,
  pub constraints: Constraints,
}

pub trait VerifyForm {
  fn verify_form(&self) -> bool;
}

impl VerifyForm for Factor {
  fn verify_form(&self) -> bool {
    self.n_coefficient as usize == self.coefficients.len()
  }
}

impl VerifyForm for R1csContents {
  fn verify_form(&self) -> bool {
    let is_valid_constraints = self.header.n_constraints as usize == self.constraints.0.len();
    //let is_valid_labels = self.header.n_labels as usize == self.wire2_labelid_map.0.len();
    let mut is_valid_factors = true;
    for constraint in self.constraints.0.iter() {
      if constraint.factors.len() != 3
        || !constraint.factors[0].verify_form()
        || !constraint.factors[1].verify_form()
        || !constraint.factors[2].verify_form()
      {
        is_valid_factors = false;
        break;
      }
    }
    return is_valid_constraints && /*is_valid_labels &&*/ is_valid_factors;
  }
}

pub fn r1cs_computational_trace<T: PrimeField>(witness: &[T], coefficients: &[T]) -> Vec<T> {
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

#[test]
fn test() {
  R1csContents {
    version: Version(1),
    header: Header {
      field_size: 32,
      prime_number: [
        53438638232309528389504892708671455233,
        64323764613183177041862057485226039389,
      ],
      n_wires: 5,
      n_public_outputs: 1,
      n_public_inputs: 2,
      n_private_inputs: 0,
      n_labels: 6,
      n_constraints: 2,
    },
    constraints: Constraints(vec![
      Constraint {
        factors: vec![
          Factor {
            n_coefficient: 1,
            coefficients: vec![Coefficient {
              wire_id: 3,
              value: [
                4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
                811880050,
              ],
            }],
          },
          Factor {
            n_coefficient: 1,
            coefficients: vec![Coefficient {
              wire_id: 3,
              value: [1, 0, 0, 0, 0, 0, 0, 0],
            }],
          },
          Factor {
            n_coefficient: 1,
            coefficients: vec![Coefficient {
              wire_id: 4,
              value: [
                4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
                811880050,
              ],
            }],
          },
        ],
      },
      Constraint {
        factors: vec![
          Factor {
            n_coefficient: 1,
            coefficients: vec![Coefficient {
              wire_id: 2,
              value: [
                4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
                811880050,
              ],
            }],
          },
          Factor {
            n_coefficient: 1,
            coefficients: vec![Coefficient {
              wire_id: 4,
              value: [1, 0, 0, 0, 0, 0, 0, 0],
            }],
          },
          Factor {
            n_coefficient: 1,
            coefficients: vec![Coefficient {
              wire_id: 1,
              value: [
                4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
                811880050,
              ],
            }],
          },
        ],
      },
    ]),
  };
}

pub struct StarkProof {
  m_root: Vec<u8>,
  l_root: Vec<u8>,
  main_branches: Vec<Vec<Vec<u8>>>,
  linear_comb_branches: Vec<Vec<Vec<u8>>>,
  fri_proof: Vec<FriProof>,
}

// const modulus = 2**256 - 2**32 * 351 + 1;
// const non_residue = 7;
const EXTENSION_FACTOR: usize = 8usize; // >= 4 (for low-degree proof)
const SPOT_CHECK_SECURITY_FACTOR: usize = 80usize;

// Generate a STARK for a MIMC calculation
pub fn mk_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
  computational_trace: &[T],
  coefficients: &[T], // This argument may be good to use HashMap<usize, T> instead of Vec<T> because we can omit zero coefficients from it.
  n_constraints: usize,
  n_wires: usize,
) -> StarkProof {
  let original_steps = 6 * n_constraints * n_wires;
  assert_eq!(computational_trace.len(), original_steps);

  let mut constants = coefficients.to_vec();
  constants.extend(coefficients.to_vec());
  assert_eq!(constants.len(), original_steps);

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
    let j2 = (j1 + n_wires * skips) % precision;
    let j3 = (j2 + n_wires * skips) % precision;
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
      && (j - (n_wires - 1) * skips) % (3 * n_wires * skips) == 0
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

  // Compute interpolant of ((1, input), (x_at_last_step, output))
  let interpolant = {
    let mut x_vals = vec![];
    let mut y_vals = vec![];
    for w in 0..n_wires {
      x_vals.push(xs[w * skips]);
      y_vals.push(computational_trace[w]);
    }

    for k in 0..n_constraints {
      let w = (3 * k + 1) * n_wires;
      x_vals.push(xs[w * skips]);
      y_vals.push(computational_trace[w]);
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
    let j = (3 * k + 1) * n_wires * skips;
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
      (j + n_wires * skips) % precision,
      (j + 2 * n_wires * skips) % precision,
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

// Verifies a STARK
pub fn verify_r1cs_proof<T: PrimeField + FromBytes + ToBytes>(
  proof: StarkProof,
  computational_trace: &[T],
  coefficients: &[T],
  n_constraints: usize,
  n_wires: usize,
) -> Result<bool, String> {
  let original_steps = 6 * n_constraints * n_wires;
  assert_eq!(computational_trace.len(), original_steps);

  let mut constants = coefficients.to_vec();
  constants.extend(coefficients.to_vec());
  assert_eq!(constants.len(), original_steps);

  let mut log_steps = 1;
  let mut tmp_steps = original_steps - 1;
  while tmp_steps > 1 {
    tmp_steps /= 2;
    log_steps += 1;
  }
  let steps = 2usize.pow(log_steps);
  let mut computational_trace = computational_trace.to_vec();
  computational_trace.extend(vec![T::zero(); steps - original_steps]);

  // start_time = time.time()
  assert!(steps <= 2usize.pow(32));

  let StarkProof {
    m_root,
    l_root,
    main_branches,
    linear_comb_branches,
    fri_proof,
  } = proof;

  let precision = steps * EXTENSION_FACTOR;
  let skips = precision / steps; // EXTENSION_FACTOR

  // Get (steps)th root of unity
  let times_nmr = BigUint::from_bytes_be(&(T::zero() - T::one()).to_bytes_be().unwrap());
  let times_dnm = BigUint::from_bytes_be(&precision.to_be_bytes());
  assert!(&times_nmr % &times_dnm == BigUint::from(0u8));
  let times = parse_bytes_to_u64_vec(&(times_nmr / times_dnm).to_bytes_le()); // (modulus - 1) / precision
  let g2 = T::multiplicative_generator().pow_vartime(&times); // g2^precision == 1 mod modulus
  let xs = expand_root_of_unity(g2);
  let g1 = xs[skips];

  // Gets the polynomial representing the constants
  let k_polynomial = inv_fft(&constants, g1);

  // Verifies the low-degree proofs
  assert!(
    verify_low_degree_proof(l_root.clone(), g2, &fri_proof, precision / 4, skips as u32).unwrap()
  );

  // Performs the spot checks
  let k0 = T::one();
  let k1 = T::from_str(&mk_seed(&[m_root.clone(), b"\x01".to_vec()])).unwrap();
  let k2 = T::from_str(&mk_seed(&[m_root.clone(), b"\x02".to_vec()])).unwrap();
  let k3 = T::from_str(&mk_seed(&[m_root.clone(), b"\x03".to_vec()])).unwrap();
  let k4 = T::from_str(&mk_seed(&[m_root.clone(), b"\x04".to_vec()])).unwrap();
  let k5 = T::from_str(&mk_seed(&[m_root.clone(), b"\x05".to_vec()])).unwrap();
  let positions = get_pseudorandom_indices(
    &l_root.clone(),
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
      (j + n_wires * skips) % precision,
      (j + 2 * n_wires * skips) % precision,
    ]);
  }

  let main_branch_leaves = verify_multi_branch(&m_root, &augmented_positions, &main_branches);
  let linear_comb_branch_leaves = verify_multi_branch(&l_root, &positions, &linear_comb_branches);

  // let z1_num_evaluations: Vec<T> = (0..precision)
  //   .map(|i| xs[(i * steps) % precision] - T::one())
  //   .collect();
  // let mut z1_den_evaluations: Vec<T> = vec![T::one(); precision];
  // for j in 0..precision {
  //   if j >= original_steps * skips || (j >= skips && (j - skips) % (n_wires * skips) == 0) {
  //     z1_den_evaluations = z1_den_evaluations
  //       .iter()
  //       .enumerate()
  //       .map(|(i, &val)| val * (xs[i] - xs[j]))
  //       .collect();
  //   }
  // }
  // let z1_den_inv = multi_inv(&z1_den_evaluations);
  let mut z1_evaluations: Vec<T> = vec![T::one(); precision];
  for j in 0..precision {
    if j < original_steps * skips / 2 && j % skips == 0 && j % (n_wires * skips) != 0 {
      z1_evaluations = z1_evaluations
        .iter()
        .enumerate()
        .map(|(i, &val)| val * (xs[i] - xs[j]))
        .collect();
    }
  }

  let mut z2_evaluations: Vec<T> = vec![T::one(); precision];
  // for k in 0..n_constraints {
  //   let j = ((3 * k + 2) * n_wires - 1) * skips;
  //   z2_evaluations = z2_evaluations
  //     .iter()
  //     .enumerate()
  //     .map(|(i, &val)| val * (xs[i] - xs[j]))
  //     .collect();
  // }

  for j in 0..precision {
    if j < original_steps * skips
      && j >= original_steps * skips / 2
      && (j - (n_wires - 1) * skips) % (3 * n_wires * skips) == 0
    {
      z2_evaluations = z2_evaluations
        .iter()
        .enumerate()
        .map(|(i, &val)| val * (xs[i] - xs[j]))
        .collect();
    }
  }

  for (i, &pos) in positions.iter().enumerate() {
    let mut m_branch0 = main_branch_leaves[i * 5].clone();
    let mut m_branch1 = main_branch_leaves[i * 5 + 1].clone();
    let mut m_branch2 = main_branch_leaves[i * 5 + 2].clone();
    let mut m_branch10 = main_branch_leaves[i * 5 + 3].clone();
    let mut m_branch11 = main_branch_leaves[i * 5 + 4].clone();
    let l_of_x = T::from_bytes_be(linear_comb_branch_leaves[i].clone()).unwrap();

    let mut m_branch3 = m_branch0.split_off(32); // m_branch0 = leaves[i * 5][..32]
    let mut m_branch4 = m_branch3.split_off(32); // m_branch3 = leaves[i * 5][32..64]
    let mut m_branch5 = m_branch4.split_off(32); // m_branch4 = leaves[i * 5][64..96]
    let _ = m_branch5.split_off(32); // m_branch5 = leaves[i * 5][96..]
    let _ = m_branch1.split_off(32); // m_branch1 = leaves[i * 5 + 1][..32]
    let _ = m_branch2.split_off(32); // m_branch2 = leaves[i * 5 + 2][..32]
    let _ = m_branch10.split_off(32); // m_branch10 = leaves[i * 5 + 3][..32]
    let _ = m_branch11.split_off(32); // m_branch11 = leaves[i * 5 + 4][..32]
    let p_of_x = T::from_bytes_be(m_branch0).unwrap();
    let p_of_prev_x_plus_half = T::from_bytes_be(m_branch1).unwrap();
    let p_of_x_plus_half = T::from_bytes_be(m_branch2).unwrap();
    let d1_of_x = T::from_bytes_be(m_branch3).unwrap();
    let d2_of_x = T::from_bytes_be(m_branch4).unwrap();
    let b_of_x = T::from_bytes_be(m_branch5).unwrap();
    let p_of_x_plus_w = T::from_bytes_be(m_branch10).unwrap();
    let p_of_x_plus_2w = T::from_bytes_be(m_branch11).unwrap();

    let z1_value = z1_evaluations[pos];
    let z2_value = z2_evaluations[pos];

    let x = xs[pos]; // g2.pow_vartime(&parse_bytes_to_u64_vec(&pos.to_le_bytes()));
    let k_of_x_plus_half = eval_poly_at(
      &k_polynomial,
      xs[(pos + 3 * n_constraints * n_wires * skips) % precision],
    );
    // Check first transition constraints Q1(x) = Z1(x) * D1(x)
    // println!(
    //   "{:03} {:?} {:?}    {:?} {:?}    {:?} {:?}",
    //   pos, p_of_x, p_of_prev_x_plus_half, p_of_x_plus_half, k_of_x_plus_half, z1_value, d1_of_x
    // );
    assert_eq!(
      p_of_x_plus_half - p_of_prev_x_plus_half - k_of_x_plus_half * p_of_x,
      z1_value * d1_of_x
    );

    // Check second transition constraints Q2(x) = Z2(x) * D2(x)
    // println!(
    //   "{:03} {:?} {:?}    {:?} {:?}    {:?} {:?}",
    //   pos,
    //   p_of_x,
    //   p_of_x_plus_w,
    //   p_of_x_plus_2w,
    //   z2_value * d2_of_x,
    //   z2_value,
    //   d2_of_x
    // );
    assert_eq!(p_of_x_plus_2w - p_of_x * p_of_x_plus_w, z2_value * d2_of_x);

    // Check boundary constraints P(x) - I(x) = Zb(x) * B(x)
    let interpolant = {
      let mut x_vals = vec![];
      let mut y_vals = vec![];
      for w in 0..n_wires {
        x_vals.push(xs[w * skips]);
        y_vals.push(computational_trace[w]);
      }

      for k in 0..n_constraints {
        let j = (3 * k + 1) * n_wires;
        x_vals.push(xs[j * skips]);
        y_vals.push(computational_trace[j]);
      }
      lagrange_interp(&x_vals, &y_vals)
    };

    let mut zb_of_x = T::one();
    for w in 0..n_wires {
      let j = w * skips;
      zb_of_x = zb_of_x * (x - xs[j]);
    }
    for k in 0..n_constraints {
      let j = (3 * k + 1) * n_wires * skips;
      zb_of_x = zb_of_x * (x - xs[j]);
    }
    // println!(
    //   "{:03} {:?} {:?}    {:?} {:?}",
    //   pos,
    //   p_of_x,
    //   eval_poly_at(&interpolant, x),
    //   zb_of_x,
    //   b_of_x
    // );
    assert_eq!(p_of_x - eval_poly_at(&interpolant, x), zb_of_x * b_of_x);

    // Check correctness of the linear combination
    let x_to_the_steps = x.pow_vartime(&parse_bytes_to_u64_vec(&steps.to_le_bytes()));
    assert_eq!(
      l_of_x,
      k0 * d1_of_x
        + k1 * d2_of_x
        + k2 * p_of_x
        + k3 * p_of_x * x_to_the_steps
        + k4 * b_of_x
        + k5 * b_of_x * x_to_the_steps,
    );
  }

  println!("Verified {} consistency checks", SPOT_CHECK_SECURITY_FACTOR);
  // println!("Verified STARK in %.4f sec" % (time.time() - start_time));
  Ok(true)
}

fn main() {
  // use crate::{mk_mimc_proof, verify_mimc_proof};
  // use stark_pure_rust::ff::PrimeFiled;
  use stark_pure_rust::fp::Fp;
  use stark_pure_rust::fri::fri_proof_bin_length;
  use stark_pure_rust::merkle_tree::bin_length;

  type TargetFF = Fp;

  // let inp = TargetFF::from(3u64);
  // Full STARK test
  let a_coeff: Vec<Vec<TargetFF>> = vec![vec![0u64, 0, 0, 1, 0], vec![0, 0, 1, 0, 0]]
    .iter()
    .map(|xs| xs.iter().map(|x| TargetFF::from(*x)).collect())
    .collect();
  let b_coeff: Vec<Vec<TargetFF>> = vec![vec![0u64, 0, 0, 1, 0], vec![0, 0, 0, 0, 1]]
    .iter()
    .map(|xs| xs.iter().map(|x| TargetFF::from(*x)).collect())
    .collect();
  let c_coeff: Vec<Vec<TargetFF>> = vec![vec![0u64, 0, 0, 0, 1], vec![0, 1, 0, 0, 0]]
    .iter()
    .map(|xs| xs.iter().map(|x| TargetFF::from(*x)).collect())
    .collect();
  let coefficients: Vec<TargetFF> =
    a_coeff
      .iter()
      .zip(&b_coeff)
      .zip(&c_coeff)
      .fold(vec![], |mut acc, ((a, b), c)| {
        acc.extend(a);
        acc.extend(b);
        acc.extend(c);
        acc
      });
  let witness: Vec<TargetFF> = vec![1u64, 12, 3, 2, 4]
    .iter()
    .map(|x| TargetFF::from(*x))
    .collect();

  // Generate the computational trace
  let computational_trace = r1cs_computational_trace(&witness, &coefficients);
  // println!("Done generating computational trace");
  let proof = mk_r1cs_proof(&computational_trace, &coefficients, 2, 5);

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

  assert!(verify_r1cs_proof(proof, &computational_trace, &coefficients, 2, 5).unwrap());
  println!("Done proof verification");
}
