use ff::PrimeField;
use fri::fri::FriProof;
use fri::utils::blake;
use num::bigint::BigUint;
use serde::{Serialize, Deserialize};

pub fn parse_hex_to_decimal(value: &[u8]) -> String {
  BigUint::from_bytes_be(value).to_str_radix(10)
}

pub fn mk_seed(messages: &[Vec<u8>]) -> String {
  let mut message: Vec<u8> = vec![];
  for m in messages {
    message.extend(m);
  }
  parse_hex_to_decimal(&blake(&message))
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

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct StarkProof {
  pub m_root: Vec<u8>,
  pub l_root: Vec<u8>,
  pub main_branches: Vec<Vec<Vec<u8>>>,
  pub linear_comb_branches: Vec<Vec<Vec<u8>>>,
  pub fri_proof: Vec<FriProof>,
}

// const modulus = 2**256 - 2**32 * 351 + 1;
// const non_residue = 7;
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
