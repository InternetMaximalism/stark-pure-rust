use crate::prove::mk_r1cs_proof;
use crate::reader::{
  read_r1cs, read_witness, Coefficient, Constraints, Factor, Header, R1csContents,
};
use crate::utils::*;
use crate::verify::verify_r1cs_proof;

use circom2bellman_core::Constraint;
use commitment::hash::Digest;
use ff::{Field, PrimeField};
use ff_utils::ff_utils::FromBytes;
use ff_utils::fp::Fp;
#[allow(unused_imports)]
use log::{debug, info};
use num::BigUint;
use serde::{Deserialize, Serialize};
use std::cmp::max;
use std::fs::File;
use std::io::Error;
use std::io::{BufReader, Read, Write};
use std::path::Path;

fn calc_coefficients<T: PrimeField + FromBytes>(
  constraints: &[Constraint],
  n_wires: usize,
) -> (Vec<T>, Vec<Vec<(u8, usize)>>, Vec<usize>) {
  let mut a_coeff_list: Vec<T> = vec![];
  let mut b_coeff_list: Vec<T> = vec![];
  let mut c_coeff_list: Vec<T> = vec![];

  let mut wire_using_list: Vec<Vec<(u8, usize)>> = vec![vec![]; n_wires as usize];
  // let mut wire_prev_list: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
  let mut acc_n_coeff = 0usize;
  let mut last_coeff_list = vec![];
  for constraint in constraints {
    let Factor {
      n_coefficient: n_a_coeff,
      coefficients: a_coefficients,
    } = &constraint.factors[0];
    let Factor {
      n_coefficient: n_b_coeff,
      coefficients: b_coefficients,
    } = &constraint.factors[1];
    let Factor {
      n_coefficient: n_c_coeff,
      coefficients: c_coefficients,
    } = &constraint.factors[2];
    let n_coeff = *max(max(n_a_coeff, n_b_coeff), n_c_coeff);

    let n_wires = n_wires as usize;

    for i in 0..n_coeff {
      if i < *n_a_coeff {
        let Coefficient { wire_id, value } = &a_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let c = T::from_bytes_le(value).unwrap();
        wire_using_list[wire_id].push((0, a_coeff_list.len()));
        a_coeff_list.push(c);
      } else {
        let wire_id = n_wires - 1;
        let c = T::zero();
        wire_using_list[wire_id].push((0, a_coeff_list.len()));
        a_coeff_list.push(c);
      }
    }

    for i in 0..n_coeff {
      if i < *n_b_coeff {
        let Coefficient { wire_id, value } = &b_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let c = T::from_bytes_le(value).unwrap();
        wire_using_list[wire_id].push((1, b_coeff_list.len()));
        b_coeff_list.push(c);
      } else {
        let wire_id = n_wires - 1;
        let c = T::zero();
        wire_using_list[wire_id].push((1, b_coeff_list.len()));
        b_coeff_list.push(c);
      }
    }

    for i in 0..n_coeff {
      if i < *n_c_coeff {
        let Coefficient { wire_id, value } = &c_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let c = T::from_bytes_le(value).unwrap();
        wire_using_list[wire_id].push((2, c_coeff_list.len()));
        c_coeff_list.push(c);
      } else {
        let wire_id = n_wires - 1;
        let c = T::zero();
        wire_using_list[wire_id].push((2, c_coeff_list.len()));
        c_coeff_list.push(c);
      }
    }

    acc_n_coeff += n_coeff as usize;
    last_coeff_list.push(acc_n_coeff - 1);
  }

  let mut coefficients = vec![];
  coefficients.extend(a_coeff_list);
  coefficients.extend(b_coeff_list);
  coefficients.extend(c_coeff_list);

  // println!("coefficients: {:?}", coefficients);

  (coefficients, wire_using_list, last_coeff_list)
}

fn calc_coefficients_and_witness<T: PrimeField + FromBytes>(
  constraints: &[Constraint],
  witness: &[T],
  n_wires: usize,
) -> (Vec<T>, Vec<T>, Vec<T>, Vec<Vec<(u8, usize)>>, Vec<usize>) {
  let mut a_wit_list: Vec<T> = vec![];
  let mut b_wit_list: Vec<T> = vec![];
  let mut c_wit_list: Vec<T> = vec![];
  let mut a_trace: Vec<T> = vec![];
  let mut b_trace: Vec<T> = vec![];
  let mut c_trace: Vec<T> = vec![];
  let mut a_coeff_list: Vec<T> = vec![];
  let mut b_coeff_list: Vec<T> = vec![];
  let mut c_coeff_list: Vec<T> = vec![];

  let mut wire_using_list: Vec<Vec<(u8, usize)>> = vec![vec![]; n_wires as usize]; // Vec<Vec<Position>>
  let mut acc_n_coeff = 0usize;
  let mut last_coeff_list = vec![];
  for (_index, constraint) in constraints.iter().enumerate() {
    let Factor {
      n_coefficient: n_a_coeff,
      coefficients: a_coefficients,
    } = &constraint.factors[0];
    let Factor {
      n_coefficient: n_b_coeff,
      coefficients: b_coefficients,
    } = &constraint.factors[1];
    let Factor {
      n_coefficient: n_c_coeff,
      coefficients: c_coefficients,
    } = &constraint.factors[2];
    let n_coeff = *max(max(n_a_coeff, n_b_coeff), n_c_coeff);

    // let target_index = |v| v == 0;

    let n_wires = n_wires as usize;

    let mut t = T::zero();
    for i in 0..n_coeff {
      if i < *n_a_coeff {
        let Coefficient { wire_id, value } = &a_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        // if target_index(index) {
        //   println!("A: {:?} {:?}", w, a_coefficients[i as usize]);
        // }
        let c = T::from_bytes_le(value).unwrap();
        // println!("w: {:?}", w);
        // println!("c: {:?}", c);
        // println!("t: {:?}", t);
        t = t + c * w;
        wire_using_list[wire_id].push((0, a_coeff_list.len()));
        a_wit_list.push(w);
        a_coeff_list.push(c);
        a_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = T::zero();
        // println!("w': {:?}", w);
        // println!("c': {:?}", c);
        // println!("t: {:?}", t);
        wire_using_list[wire_id].push((0, a_coeff_list.len()));
        a_wit_list.push(w);
        a_coeff_list.push(c);
        a_trace.push(t);
      }
    }
    // if target_index(index) {
    //   println!("t_A: {:?}", t);
    // }

    let mut t = T::zero();
    for i in 0..n_coeff {
      if i < *n_b_coeff {
        let Coefficient { wire_id, value } = &b_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        // if target_index(index) {
        //   println!("B: {:?} {:?}", w, b_coefficients[i as usize]);
        // }
        let c = T::from_bytes_le(value).unwrap();
        // println!("w: {:?}", w);
        // println!("c: {:?}", c);
        // println!("t: {:?}", t);
        t = t + c * w;
        wire_using_list[wire_id].push((1, b_coeff_list.len()));
        b_wit_list.push(w);
        b_coeff_list.push(c);
        b_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = T::zero();
        // println!("w': {:?}", w);
        // println!("c': {:?}", c);
        // println!("t: {:?}", t);
        wire_using_list[wire_id].push((1, b_coeff_list.len()));
        b_wit_list.push(w);
        b_coeff_list.push(c);
        b_trace.push(t);
      }
    }
    // if target_index(index) {
    //   println!("t_B: {:?}", t);
    // }

    let mut t = T::zero();
    for i in 0..n_coeff {
      if i < *n_c_coeff {
        let Coefficient { wire_id, value } = &c_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        // if target_index(index) {
        //   println!("C: {:?} {:?}", w, c_coefficients[i as usize]);
        // }
        let c = T::from_bytes_le(value).unwrap();
        // println!("w: {:?}", w);
        // println!("c: {:?}", c);
        // println!("t: {:?}", t);
        t = t + c * w;
        wire_using_list[wire_id].push((2, c_coeff_list.len()));
        c_wit_list.push(w);
        c_coeff_list.push(c);
        c_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = T::zero();
        // println!("w': {:?}", w);
        // println!("c': {:?}", c);
        // println!("t: {:?}", t);
        wire_using_list[wire_id].push((2, c_coeff_list.len()));
        c_wit_list.push(w);
        c_coeff_list.push(c);
        c_trace.push(t);
      }
    }
    // if target_index(index) {
    //   println!("t_C: {:?}", t);
    // }

    acc_n_coeff += n_coeff as usize;
    last_coeff_list.push(acc_n_coeff - 1);
  }

  let mut witness_trace = vec![];
  witness_trace.extend(a_wit_list);
  witness_trace.extend(b_wit_list);
  witness_trace.extend(c_wit_list);

  let mut computational_trace = vec![];
  computational_trace.extend(a_trace);
  computational_trace.extend(b_trace);
  computational_trace.extend(c_trace);

  let mut coefficients = vec![];
  coefficients.extend(a_coeff_list);
  coefficients.extend(b_coeff_list);
  coefficients.extend(c_coeff_list);

  // println!("coefficients: {:?}", coefficients);
  debug_assert_eq!(witness_trace.len(), coefficients.len());
  debug_assert_eq!(computational_trace.len(), coefficients.len());

  (
    witness_trace,
    computational_trace,
    coefficients,
    wire_using_list,
    last_coeff_list,
  )
}

fn calc_flags<T: PrimeField>(
  last_coeff_list: &[usize],
  coefficients_len: usize,
) -> (Vec<T>, Vec<T>, Vec<T>) {
  assert_eq!(coefficients_len % 3, 0);
  let a_trace_len = coefficients_len / 3;

  let flag0 = vec![T::one(); coefficients_len];
  let mut flag1 = vec![T::one(); coefficients_len];
  let first_coeff_list = last_coeff_list
    .iter()
    .map(|v| (v + 1) % a_trace_len)
    .collect::<Vec<_>>();
  for k in first_coeff_list {
    flag1[k] = T::zero();
    flag1[k + a_trace_len] = T::zero();
    flag1[k + a_trace_len * 2] = T::zero();
  }

  let mut flag2 = vec![T::zero(); coefficients_len];
  for k in last_coeff_list {
    flag2[*k] = T::from(1u64);
  }

  (flag0, flag1, flag2)
}

pub fn prove_with_witness<H: Digest>(
  r1cs: &R1csContents,
  witness: &[Vec<u8>],
) -> Result<StarkProof<H>, Error> {
  let Header {
    field_size: _,
    prime_number,
    n_public_inputs,
    n_public_outputs,
    n_private_inputs: _,
    n_labels,
    n_constraints,
    n_wires,
  } = r1cs.header;
  let Constraints(constraints) = &r1cs.constraints;

  let n_constraints = n_constraints as usize;
  let n_wires = n_wires as usize;
  println!("n_constraints: {:?}", n_constraints);
  println!("n_wires: {:?}", n_wires);
  println!(
    "n_public_wires: {:?} + {:?}",
    n_public_inputs, n_public_outputs
  );
  println!("n_labels: {:?}", n_labels);

  // println!("prime_number: {:?}", prime_number);
  // assert_eq!(
  //   prime_number,
  //   [
  //     53438638232309528389504892708671455233,
  //     64323764613183177041862057485226039389
  //   ]
  // );
  assert_eq!(
    prime_number,
    [
      1, 0, 0, 240, 147, 245, 225, 67, 145, 112, 185, 121, 72, 232, 51, 40, 93, 88, 129, 129, 182,
      69, 80, 184, 41, 160, 49, 225, 114, 78, 100, 48
    ]
  );
  type TargetFF = Fp; // TODO: Use r1cs.header.field_size.

  // assert_eq!(witness.len(), n_wires as usize);
  let witness: Vec<TargetFF> = witness
    .iter()
    .map(|x| TargetFF::from_bytes_le(x).unwrap())
    .collect();
  assert_eq!(witness[0], TargetFF::one());
  let public_wires = witness[..(1 + n_public_inputs as usize + n_public_outputs as usize)].to_vec();

  // Generate the computational trace
  println!("Generate coefficients");
  let start = std::time::Instant::now();
  let (witness_trace, computational_trace, coefficients, wire_using_list, last_coeff_list) =
    calc_coefficients_and_witness(constraints, &witness, n_wires as usize);
  let end: std::time::Duration = start.elapsed();
  println!(
    "Generated coefficients: {}.{:03}s",
    end.as_secs(),
    end.subsec_nanos() / 1_000_000
  );
  // println!("last_coeff_list: {:?}", last_coeff_list);

  let start = std::time::Instant::now();

  let (flag0, flag1, flag2) = calc_flags(&last_coeff_list, coefficients.len());
  let end: std::time::Duration = start.elapsed();
  println!(
    "Generated flags: {}.{:03}s",
    end.as_secs(),
    end.subsec_nanos() / 1_000_000
  );

  let a_trace_len = coefficients.len() / 3;
  println!("a_trace_len: {}", a_trace_len);

  println!("permuted_indices");
  let start = std::time::Instant::now();
  let mut permuted_indices = vec![0usize; computational_trace.len()];
  for vs in wire_using_list.iter() {
    if vs.len() == 0 {
      continue;
    }
    let mut old_w = a_trace_len * vs.last().unwrap().0 as usize + vs.last().unwrap().1;
    for &(k, v) in vs.iter() {
      let w = a_trace_len * k as usize + v;
      permuted_indices[w] = old_w;
      old_w = w;
    }
  }
  // println!("permuted_indices: {:?}", permuted_indices);
  // println!("wire_using_list: {:?}", wire_using_list);
  let end: std::time::Duration = start.elapsed();
  println!(
    "Generated permuted_indices: {}.{:03}s",
    end.as_secs(),
    end.subsec_nanos() / 1_000_000
  );

  println!("public_first_indices");
  let start = std::time::Instant::now();
  let mut public_first_indices = vec![];
  for w in 0..public_wires.len() {
    if wire_using_list[w].len() > 0 {
      let (k, v) = *wire_using_list[w].first().unwrap();
      public_first_indices.push((w, a_trace_len * k as usize + v));
    }
  }
  // println!("public_first_indices: {:?}", public_first_indices);
  let end: std::time::Duration = start.elapsed();
  println!(
    "Generated public_first_indices: {}.{:03}s",
    end.as_secs(),
    end.subsec_nanos() / 1_000_000
  );

  // test_mk_d2_proof(
  //   &witness_trace,
  //   &computational_trace,
  //   &permuted_indices,
  //   &coefficients,
  //   &flag2,
  //   n_constraints,
  //   n_wires,
  // )
  // .unwrap();

  mk_r1cs_proof::<TargetFF, H>(
    &witness_trace,
    &computational_trace,
    &public_wires,
    &public_first_indices,
    &permuted_indices,
    &coefficients,
    &flag0,
    &flag1,
    &flag2,
    n_constraints,
    n_wires,
  )
}

fn verify_with_witness<H: Digest>(
  r1cs: &R1csContents,
  public_wires: &[Vec<u8>],
  proof: StarkProof<H>,
) -> Result<bool, Error> {
  let Header {
    field_size: _,
    prime_number: _,
    n_public_inputs: _,
    n_public_outputs: _,
    n_private_inputs: _,
    n_labels: _,
    n_constraints,
    n_wires,
  } = r1cs.header;
  let Constraints(constraints) = &r1cs.constraints;

  let n_constraints = n_constraints as usize;
  let n_wires = n_wires as usize;

  type TargetFF = Fp; // TODO: Use r1cs.header.field_size.

  let public_wires: Vec<TargetFF> = public_wires
    .iter()
    .map(|x| TargetFF::from_bytes_le(x).unwrap())
    .collect();
  assert_eq!(public_wires[0], TargetFF::one());

  println!("Generate coefficients");
  let (coefficients, wire_using_list, last_coeff_list) =
    calc_coefficients(constraints, n_wires as usize);

  let a_trace_len = coefficients.len() / 3;

  println!("Generate flags");
  let (flag0, flag1, flag2) = calc_flags(&last_coeff_list, coefficients.len());

  println!("permuted_indices");
  let mut permuted_indices = vec![0usize; coefficients.len()];
  for vs in wire_using_list.iter() {
    if vs.len() == 0 {
      continue;
    }
    let mut old_w = a_trace_len * vs.last().unwrap().0 as usize + vs.last().unwrap().1;
    for &(k, v) in vs.iter() {
      let w = a_trace_len * k as usize + v;
      permuted_indices[w] = old_w;
      old_w = w;
    }
  }

  println!("public_first_indices");
  let mut public_first_indices = vec![];
  for w in 0..public_wires.len() {
    if wire_using_list[w].len() > 0 {
      let (k, v) = *wire_using_list[w].first().unwrap();
      public_first_indices.push((w, a_trace_len * k as usize + v));
    }
  }

  verify_r1cs_proof(
    proof,
    &public_wires,
    &public_first_indices,
    &permuted_indices,
    &coefficients,
    &flag0,
    &flag1,
    &flag2,
    n_constraints,
    n_wires,
  )
}

pub fn prove_with_file_path<
  P: AsRef<Path>,
  Q: AsRef<Path>,
  R: AsRef<Path>,
  H: Digest + Serialize,
>(
  r1cs_file_path: P,
  witness_file_path: Q,
  proof_json_path: R,
) -> Result<(), Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_r1cs(&raw_r1cs);

  let witness_file = File::open(witness_file_path)?;
  let mut raw_witness = vec![];
  BufReader::new(witness_file).read_to_end(&mut raw_witness)?;
  let witness = read_witness(&raw_witness);

  let proof = prove_with_witness::<H>(&r1cs, &witness)?;
  let serialized_proof = serde_json::to_string(&proof)?;
  let mut file = File::create(proof_json_path)?;
  write!(file, "{}", serialized_proof)?;

  Ok(())
}

pub fn verify_with_file_path<
  P: AsRef<Path>,
  Q: AsRef<Path>,
  R: AsRef<Path>,
  H: Digest + for<'a> Deserialize<'a>,
>(
  r1cs_file_path: P,
  witness_file_path: Q,
  proof_json_path: R,
) -> Result<(), Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_r1cs(&raw_r1cs);

  let witness_file = File::open(witness_file_path)?;
  let mut raw_witness = vec![];
  BufReader::new(witness_file).read_to_end(&mut raw_witness)?;
  let witness = read_witness(&raw_witness);

  let proof_file = File::open(proof_json_path)?;
  let proof_reader = BufReader::new(proof_file);
  let proof: StarkProof<H> = serde_json::from_reader(proof_reader)?;

  // TODO: Use public_wires instead of witness.
  let public_wires = witness
    [..(1 + r1cs.header.n_public_inputs as usize + r1cs.header.n_public_outputs as usize)]
    .to_vec();
  assert!(verify_with_witness(&r1cs, &public_wires, proof)?);
  println!("Done proof verification");

  Ok(())
}

pub fn run_with_file_path<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>, H: Digest + Serialize>(
  r1cs_file_path: P,
  witness_file_path: Q,
  proof_json_path: R,
) -> Result<(), Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_r1cs(&raw_r1cs);

  let witness_file = File::open(witness_file_path)?;
  let mut raw_witness = vec![];
  BufReader::new(witness_file).read_to_end(&mut raw_witness)?;
  let witness = read_witness(&raw_witness);

  let witness_ff = witness
    .iter()
    .map(|v| BigUint::from_bytes_le(v))
    .collect::<Vec<_>>();
  let wtns_json_file_path = "./tests/sha256_2_test.wtns.json";
  let mut wtns_json_file = File::create(wtns_json_file_path)?;
  write!(wtns_json_file, "{:?}", witness_ff)?;

  let proof = prove_with_witness::<H>(&r1cs, &witness)?;
  let serialized_proof = serde_json::to_string(&proof)?;
  let mut file = File::create(proof_json_path)?;
  write!(file, "{}", serialized_proof)?;

  let public_wires = witness
    [..(1 + r1cs.header.n_public_inputs as usize + r1cs.header.n_public_outputs as usize)]
    .to_vec();
  assert!(verify_with_witness(&r1cs, &public_wires, proof)?);
  println!("Done proof verification");

  Ok(())
}

#[test]
fn test_prove_and_verify_with_file_path() {
  use commitment::blake::BlakeDigest;

  let r1cs_file_path = "./tests/compute.r1cs";
  let witness_file_path = "./tests/compute.wtns";
  let proof_json_path = "./tests/compute_proof.json";
  prove_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path)
    .unwrap();
  verify_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path)
    .unwrap();
}

#[test]
fn test_run_compute_with_file_path() {
  use commitment::blake::BlakeDigest;

  let r1cs_file_path = "./tests/compute.r1cs";
  let witness_file_path = "./tests/compute.wtns";
  let proof_json_path = "./tests/compute_proof.json";
  run_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path)
    .unwrap();
}

#[test]
fn test_run_pedersen_with_file_path() {
  use commitment::blake::BlakeDigest;

  let r1cs_file_path = "./tests/pedersen_test.r1cs";
  let witness_file_path = "./tests/pedersen_test.wtns";
  let proof_json_path = "./tests/pedersen_test_proof.json";
  run_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path)
    .unwrap();
}

#[test]
fn test_run_poseidon_with_file_path() {
  use commitment::blake::BlakeDigest;

  let r1cs_file_path = "./tests/poseidon3_test.r1cs";
  let witness_file_path = "./tests/poseidon3_test.wtns";
  let proof_json_path = "./tests/poseidon3_test_proof.json";
  run_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path)
    .unwrap();
}
