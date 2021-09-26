use ff::PrimeField;
use fri::merkle_tree2::Proof;
use fri::utils::blake;
use fri::{fri::FriProof, merkle_tree2::BlakeDigest};
use num::bigint::BigUint;
use serde::{Deserialize, Serialize};

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
pub struct StarkProof {
  pub m_root: BlakeDigest,
  pub l_root: BlakeDigest,
  pub a_root: BlakeDigest,
  pub main_branches: Vec<Proof<Vec<u8>, BlakeDigest>>,
  pub linear_comb_branches: Vec<Proof<Vec<u8>, BlakeDigest>>,
  pub fri_proof: Vec<FriProof>,
}

// const modulus = 2**256 - 2**32 * 351 + 1;
// const non_residue = 7;
pub const LOG_EXTENSION_FACTOR: usize = 3usize;
pub const EXTENSION_FACTOR: usize = 8usize; // >= 4 (for low-degree proof)
pub const SPOT_CHECK_SECURITY_FACTOR: usize = 80usize;

// fn get_r1cs_contents() -> R1csContents {
//   R1csContents {
//     version: Version(1),
//     header: Header {
//       field_size: 32,
//       prime_number: [
//         53438638232309528389504892708671455233,
//         64323764613183177041862057485226039389,
//       ],
//       n_wires: 5,
//       n_public_outputs: 1,
//       n_public_inputs: 2,
//       n_private_inputs: 0,
//       n_labels: 6,
//       n_constraints: 2,
//     },
//     constraints: Constraints(vec![
//       Constraint {
//         factors: vec![
//           Factor {
//             n_coefficient: 1,
//             coefficients: vec![Coefficient {
//               wire_id: 3,
//               value: [
//                 4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
//                 811880050,
//               ],
//             }],
//           },
//           Factor {
//             n_coefficient: 1,
//             coefficients: vec![Coefficient {
//               wire_id: 3,
//               value: [1, 0, 0, 0, 0, 0, 0, 0],
//             }],
//           },
//           Factor {
//             n_coefficient: 1,
//             coefficients: vec![Coefficient {
//               wire_id: 4,
//               value: [
//                 4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
//                 811880050,
//               ],
//             }],
//           },
//         ],
//       },
//       Constraint {
//         factors: vec![
//           Factor {
//             n_coefficient: 1,
//             coefficients: vec![Coefficient {
//               wire_id: 2,
//               value: [
//                 4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
//                 811880050,
//               ],
//             }],
//           },
//           Factor {
//             n_coefficient: 1,
//             coefficients: vec![Coefficient {
//               wire_id: 4,
//               value: [1, 0, 0, 0, 0, 0, 0, 0],
//             }],
//           },
//           Factor {
//             n_coefficient: 1,
//             coefficients: vec![Coefficient {
//               wire_id: 1,
//               value: [
//                 4026531840, 1138881939, 2042196113, 674490440, 2172737629, 3092268470, 3778125865,
//                 811880050,
//               ],
//             }],
//           },
//         ],
//       },
//     ]),
//   }
// }
