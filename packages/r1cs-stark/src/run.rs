use crate::prove::mk_r1cs_proof;
use crate::utils::*;
use crate::verify::verify_r1cs_proof;
use circom2bellman_core::{read_bytes, Constraints, Factor, Header, R1csContents};
use ff::Field;
use fri::ff_utils::FromBytes;
use fri::fp::Fp;
use fri::fri::fri_proof_bin_length;
use fri::merkle_tree::bin_length;
use num::bigint::BigUint;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;

pub fn prove_with_witness(r1cs: &R1csContents, witness: &[Vec<u8>]) -> StarkProof {
  let Header {
    field_size: _,
    prime_number: _,
    n_public_outputs: _,
    n_public_inputs: _,
    n_private_inputs: _,
    n_labels: _,
    n_constraints,
    n_wires,
  } = r1cs.header;
  let Constraints(constraints) = &r1cs.constraints;

  type TargetFF = Fp; // TODO: Use field_size.

  let coefficients: Vec<TargetFF> = constraints
    .iter()
    .map(|constraint| {
      let Factor {
        n_coefficient,
        coefficients,
      } = &constraint.factors[0];
      let mut a_coeff = vec![TargetFF::zero(); n_wires as usize];
      for x in coefficients {
        a_coeff[x.wire_id as usize] =
          TargetFF::from_bytes_be(n_coefficient.to_be_bytes().to_vec()).unwrap();
      }

      let Factor {
        n_coefficient,
        coefficients,
      } = &constraint.factors[1];
      let mut b_coeff = vec![TargetFF::zero(); n_wires as usize];
      for x in coefficients {
        b_coeff[x.wire_id as usize] =
          TargetFF::from_bytes_be(n_coefficient.to_be_bytes().to_vec()).unwrap();
      }

      let Factor {
        n_coefficient,
        coefficients,
      } = &constraint.factors[2];
      let mut c_coeff = vec![TargetFF::zero(); n_wires as usize];
      for x in coefficients {
        c_coeff[x.wire_id as usize] =
          TargetFF::from_bytes_be(n_coefficient.to_be_bytes().to_vec()).unwrap();
      }

      let mut acc = vec![];
      acc.extend(a_coeff);
      acc.extend(b_coeff);
      acc.extend(c_coeff);
      acc
    })
    .fold(vec![], |mut acc, val| {
      acc.extend(val);
      acc
    });

  let witness: Vec<TargetFF> = witness
    .iter()
    .map(|x| TargetFF::from_bytes_be(x.to_vec()).unwrap())
    .collect();

  // Generate the computational trace
  let computational_trace = r1cs_computational_trace(&coefficients, &witness);
  // println!("Done generating computational trace");

  mk_r1cs_proof(
    &computational_trace,
    &coefficients,
    n_constraints as usize,
    n_wires as usize,
  )
}

// TODO: Input boundary conditions instead of witness.
fn verify_with_witness(r1cs: &R1csContents, witness: &[Vec<u8>], proof: StarkProof) {
  let Header {
    field_size: _,
    prime_number: _,
    n_public_outputs: _,
    n_public_inputs: _,
    n_private_inputs: _,
    n_labels: _,
    n_constraints,
    n_wires,
  } = r1cs.header;
  let Constraints(constraints) = &r1cs.constraints;

  type TargetFF = Fp; // TODO: Use field_size.

  let coefficients: Vec<TargetFF> = constraints
    .iter()
    .map(|constraint| {
      let Factor {
        n_coefficient,
        coefficients,
      } = &constraint.factors[0];
      let mut a_coeff = vec![TargetFF::zero(); n_wires as usize];
      for x in coefficients {
        a_coeff[x.wire_id as usize] =
          TargetFF::from_bytes_be(n_coefficient.to_be_bytes().to_vec()).unwrap();
      }

      let Factor {
        n_coefficient,
        coefficients,
      } = &constraint.factors[1];
      let mut b_coeff = vec![TargetFF::zero(); n_wires as usize];
      for x in coefficients {
        b_coeff[x.wire_id as usize] =
          TargetFF::from_bytes_be(n_coefficient.to_be_bytes().to_vec()).unwrap();
      }

      let Factor {
        n_coefficient,
        coefficients,
      } = &constraint.factors[2];
      let mut c_coeff = vec![TargetFF::zero(); n_wires as usize];
      for x in coefficients {
        c_coeff[x.wire_id as usize] =
          TargetFF::from_bytes_be(n_coefficient.to_be_bytes().to_vec()).unwrap();
      }

      let mut acc = vec![];
      acc.extend(a_coeff);
      acc.extend(b_coeff);
      acc.extend(c_coeff);
      acc
    })
    .fold(vec![], |mut acc, val| {
      acc.extend(val);
      acc
    });

  let witness: Vec<TargetFF> = witness
    .iter()
    .map(|x| TargetFF::from_bytes_be(x.to_vec()).unwrap())
    .collect();

  assert!(verify_r1cs_proof(
    proof,
    &witness,
    &coefficients,
    n_constraints as usize,
    n_wires as usize
  )
  .unwrap());
  println!("Done proof verification");
}

pub fn prove_with_file_path<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>>(
  r1cs_file_path: P,
  witness_json_path: Q,
  proof_json_path: R,
) -> Result<(), std::io::Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_bytes(&raw_r1cs);

  let witness_file = File::open(witness_json_path)?;
  let witness_reader = BufReader::new(witness_file);
  let witness: Vec<String> = serde_json::from_reader(witness_reader)?;
  let witness: Vec<Vec<u8>> = witness
    .iter()
    .map(|x| x.parse::<BigUint>().unwrap().to_bytes_be())
    .collect();

  let proof = prove_with_witness(&r1cs, &witness);
  let serialized_proof = serde_json::to_string(&proof)?;
  let mut file = File::create(proof_json_path)?;
  write!(file, "{}", serialized_proof)?;

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

  Ok(())
}

pub fn verify_with_file_path<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>>(
  r1cs_file_path: P,
  witness_json_path: Q,
  proof_json_path: R,
) -> Result<(), std::io::Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_bytes(&raw_r1cs);

  let witness_file = File::open(witness_json_path)?;
  let witness_reader = BufReader::new(witness_file);
  let witness: Vec<String> = serde_json::from_reader(witness_reader)?;
  let witness: Vec<Vec<u8>> = witness
    .iter()
    .map(|x| x.parse::<BigUint>().unwrap().to_bytes_be())
    .collect();

  let proof_file = File::open(proof_json_path)?;
  let proof_reader = BufReader::new(proof_file);
  let proof: StarkProof = serde_json::from_reader(proof_reader)?;

  verify_with_witness(&r1cs, &witness, proof);

  Ok(())
}

pub fn run_with_file_path<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>>(
  r1cs_file_path: P,
  witness_json_path: Q,
  proof_json_path: R,
) -> Result<(), std::io::Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file)
    .read_to_end(&mut raw_r1cs)
    ?;
  let r1cs = read_bytes(&raw_r1cs);

  let witness_file = File::open(witness_json_path)?;
  let witness_reader = BufReader::new(witness_file);
  let witness: Vec<String> = serde_json::from_reader(witness_reader)?;
  let witness: Vec<Vec<u8>> = witness
    .iter()
    .map(|x| x.parse::<BigUint>().unwrap().to_bytes_be())
    .collect();

  let proof = prove_with_witness(&r1cs, &witness);
  let serialized_proof = serde_json::to_string(&proof)?;
  let mut file = File::create(proof_json_path)?;
  write!(file, "{}", serialized_proof)?;
  verify_with_witness(&r1cs, &witness, proof);

  Ok(())
}

#[test]
fn test_run_with_file_path() {
  let r1cs_file_path = "./tests/mul_bn128.r1cs";
  let witness_json_path = "./tests/mul_bn128_wtns_valid.json";
  let proof_json_path = "./tests/mul_bn128_proof.json";
  prove_with_file_path(r1cs_file_path, witness_json_path, proof_json_path).unwrap();
  verify_with_file_path(r1cs_file_path, witness_json_path, proof_json_path).unwrap();
}
