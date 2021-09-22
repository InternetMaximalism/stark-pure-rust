use crate::prove2::mk_r1cs_proof;
use crate::utils::*;
use crate::verify2::verify_r1cs_proof;
use circom2bellman_core::{read_bytes, Coefficient, Constraints, Factor, Header, R1csContents};
use ff::Field;
use ff_utils::ff_utils::FromBytes;
use ff_utils::fp::Fp;
// use fri::fri::fri_proof_bin_length;
// use fri::merkle_tree::bin_length;
use num::bigint::BigUint;
use std::cmp::max;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;

pub fn prove_with_witness(r1cs: &R1csContents, witness: &[Vec<u8>]) -> StarkProof {
  let Header {
    field_size: _,
    prime_number,
    n_public_inputs,
    n_public_outputs,
    n_private_inputs: _,
    n_labels: _,
    n_constraints,
    n_wires,
  } = r1cs.header;
  let Constraints(constraints) = &r1cs.constraints;

  let n_constraints = n_constraints as usize;
  let n_wires = n_wires as usize;
  println!("n_constraints: {:?}", n_constraints);
  println!("n_wires: {:?}", n_wires);

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

  // assert_eq!(witness.len(), n_wires as usize);
  let witness: Vec<TargetFF> = witness
    .iter()
    .map(|x| TargetFF::from_bytes_le(x).unwrap())
    .collect();
  assert_eq!(witness[0], TargetFF::one());
  let public_wires = witness[..(1 + n_public_inputs as usize + n_public_outputs as usize)].to_vec();

  type TargetFF = Fp; // TODO: Use r1cs.header.field_size.

  println!("Generate coefficients");
  let mut a_wit_list: Vec<TargetFF> = vec![];
  let mut b_wit_list: Vec<TargetFF> = vec![];
  let mut c_wit_list: Vec<TargetFF> = vec![];
  let mut a_trace: Vec<TargetFF> = vec![];
  let mut b_trace: Vec<TargetFF> = vec![];
  let mut c_trace: Vec<TargetFF> = vec![];
  let mut a_coeff_list: Vec<TargetFF> = vec![];
  let mut b_coeff_list: Vec<TargetFF> = vec![];
  let mut c_coeff_list: Vec<TargetFF> = vec![];

  let mut wire_using_list: Vec<Vec<(u8, usize)>> = vec![vec![]; n_wires as usize]; // Vec<Vec<Position>>
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

    let mut t = TargetFF::zero();
    for i in 0..n_coeff {
      if i < *n_a_coeff {
        let Coefficient { wire_id, value } = &a_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        let c = TargetFF::from_bytes_le(value).unwrap();
        // println!("w: {:?}", w);
        // println!("c: {:?}", c);
        // println!("t: {:?}", t);
        t = t + c * w;
        wire_using_list[wire_id].push((0, a_trace.len()));
        a_wit_list.push(w);
        a_coeff_list.push(c);
        a_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = TargetFF::zero();
        // println!("w': {:?}", w);
        // println!("c': {:?}", c);
        // println!("t: {:?}", t);
        wire_using_list[wire_id].push((0, a_trace.len()));
        a_wit_list.push(w);
        a_coeff_list.push(c);
        a_trace.push(t);
      }
    }

    let mut t = TargetFF::zero();
    for i in 0..n_coeff {
      if i < *n_b_coeff {
        let Coefficient { wire_id, value } = &b_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        let c = TargetFF::from_bytes_le(value).unwrap();
        // println!("w: {:?}", w);
        // println!("c: {:?}", c);
        // println!("t: {:?}", t);
        t = t + c * w;
        wire_using_list[wire_id].push((1, b_trace.len()));
        b_wit_list.push(w);
        b_coeff_list.push(c);
        b_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = TargetFF::zero();
        // println!("w': {:?}", w);
        // println!("c': {:?}", c);
        // println!("t: {:?}", t);
        wire_using_list[wire_id].push((1, b_trace.len()));
        b_wit_list.push(w);
        b_coeff_list.push(c);
        b_trace.push(t);
      }
    }

    let mut t = TargetFF::zero();
    for i in 0..n_coeff {
      if i < *n_c_coeff {
        let Coefficient { wire_id, value } = &c_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        let c = TargetFF::from_bytes_le(value).unwrap();
        // println!("w: {:?}", w);
        // println!("c: {:?}", c);
        // println!("t: {:?}", t);
        t = t + c * w;
        wire_using_list[wire_id].push((2, c_trace.len()));
        c_wit_list.push(w);
        c_coeff_list.push(c);
        c_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = TargetFF::zero();
        // println!("w': {:?}", w);
        // println!("c': {:?}", c);
        // println!("t: {:?}", t);
        wire_using_list[wire_id].push((2, c_trace.len()));
        c_wit_list.push(w);
        c_coeff_list.push(c);
        c_trace.push(t);
      }
    }

    acc_n_coeff += n_coeff as usize;
    last_coeff_list.push(acc_n_coeff - 1);
  }

  let a_trace_len = a_trace.len();

  let mut flag0 = vec![];
  let mut flag1 = vec![];
  let mut flag2 = vec![];
  for k in 0..(a_coeff_list.len() + b_coeff_list.len() + c_coeff_list.len()) {
    let f0 = 1u64;
    let f1: u64 = if !last_coeff_list.contains(&((k + a_trace_len - 1) % a_trace_len)) {
      1
    } else {
      0
    };
    let f2: u64 = if last_coeff_list.contains(&k) { 1 } else { 0 };
    flag0.push(TargetFF::from(f0));
    flag1.push(TargetFF::from(f1));
    flag2.push(TargetFF::from(f0 * f2));
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
  debug_assert_eq!(coefficients.len(), witness_trace.len());
  debug_assert_eq!(computational_trace.len(), witness_trace.len());

  // Generate the computational trace
  println!("Done generating computational trace");

  println!("permuted_indices");
  // TODO: take a long time
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

  println!("public_first_indices");
  let mut public_first_indices = vec![];
  for w in 0..public_wires.len() {
    if wire_using_list[w].len() > 0 {
      let (k, v) = *wire_using_list[w].first().unwrap();
      public_first_indices.push((w, a_trace_len * k as usize + v));
    }
  }
  // println!("public_first_indices: {:?}", public_first_indices);

  mk_r1cs_proof(
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

fn verify_with_witness(r1cs: &R1csContents, witness: &[Vec<u8>], proof: StarkProof) {
  let Header {
    field_size: _,
    prime_number: _,
    n_public_inputs,
    n_public_outputs,
    n_private_inputs: _,
    n_labels: _,
    n_constraints,
    n_wires,
  } = r1cs.header;
  let Constraints(constraints) = &r1cs.constraints;

  let n_constraints = n_constraints as usize;
  let n_wires = n_wires as usize;

  type TargetFF = Fp; // TODO: Use r1cs.header.field_size.

  let witness: Vec<TargetFF> = witness
    .iter()
    .map(|x| TargetFF::from_bytes_le(x).unwrap())
    .collect();
  assert_eq!(witness[0], TargetFF::one());
  let public_wires = witness[..(1 + n_public_inputs as usize + n_public_outputs as usize)].to_vec();

  let mut a_wit_list: Vec<TargetFF> = vec![];
  let mut b_wit_list: Vec<TargetFF> = vec![];
  let mut c_wit_list: Vec<TargetFF> = vec![];
  let mut a_trace: Vec<TargetFF> = vec![];
  let mut b_trace: Vec<TargetFF> = vec![];
  let mut c_trace: Vec<TargetFF> = vec![];
  let mut a_coeff_list: Vec<TargetFF> = vec![];
  let mut b_coeff_list: Vec<TargetFF> = vec![];
  let mut c_coeff_list: Vec<TargetFF> = vec![];

  let mut wire_using_list: Vec<Vec<(usize, usize)>> = vec![vec![]; n_wires as usize];
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
      let mut t = TargetFF::zero();
      if i < *n_a_coeff {
        let Coefficient { wire_id, value } = &a_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        let c = TargetFF::from_bytes_le(value).unwrap();
        t = t + c * w;
        wire_using_list[wire_id].push((0, a_trace.len()));
        a_wit_list.push(w);
        a_coeff_list.push(c);
        a_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = TargetFF::zero();
        wire_using_list[wire_id].push((0, a_trace.len()));
        a_wit_list.push(w);
        a_coeff_list.push(c);
        a_trace.push(t);
      }
    }

    for i in 0..n_coeff {
      let mut t = TargetFF::zero();
      if i < *n_b_coeff {
        let Coefficient { wire_id, value } = &b_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        let c = TargetFF::from_bytes_le(value).unwrap();
        t = t + c * w;
        wire_using_list[wire_id].push((1, b_trace.len()));
        b_wit_list.push(w);
        b_coeff_list.push(c);
        b_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = TargetFF::zero();
        wire_using_list[wire_id].push((1, b_trace.len()));
        b_wit_list.push(w);
        b_coeff_list.push(c);
        b_trace.push(t);
      }
    }

    for i in 0..n_coeff {
      let mut t = TargetFF::zero();
      if i < *n_c_coeff {
        let Coefficient { wire_id, value } = &c_coefficients[i as usize];
        let wire_id = *wire_id as usize;
        let w = witness[wire_id];
        let c = TargetFF::from_bytes_le(value).unwrap();
        t = t + c * w;
        wire_using_list[wire_id].push((2, c_trace.len()));
        c_wit_list.push(w);
        c_coeff_list.push(c);
        c_trace.push(t);
      } else {
        let wire_id = n_wires - 1;
        let w = witness[wire_id];
        let c = TargetFF::zero();
        wire_using_list[wire_id].push((2, c_trace.len()));
        c_wit_list.push(w);
        c_coeff_list.push(c);
        c_trace.push(t);
      }
    }

    acc_n_coeff += n_coeff as usize;
    last_coeff_list.push(acc_n_coeff - 1);
  }

  let a_trace_len = a_trace.len();

  let mut flag0 = vec![];
  let mut flag1 = vec![];
  let mut flag2 = vec![];
  for k in 0..(a_coeff_list.len() + b_coeff_list.len() + c_coeff_list.len()) {
    let f0 = 1u64;
    let f1: u64 = if !last_coeff_list.contains(&((k + a_trace_len - 1) % a_trace_len)) {
      1
    } else {
      0
    };
    let f2: u64 = if last_coeff_list.contains(&k) { 1 } else { 0 };
    flag0.push(TargetFF::from(f0));
    flag1.push(TargetFF::from(f1));
    flag2.push(TargetFF::from(f0 * f2));
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
  debug_assert_eq!(coefficients.len(), witness_trace.len());
  debug_assert_eq!(computational_trace.len(), witness_trace.len());

  println!("permuted_indices");
  // TODO: take a long time
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

  // let wire_prev_list: Vec<(usize, usize)> = wire_prev_list
  //   .iter()
  //   .map(|((k, v), (prev_k, prev_v))| (k * a_trace.len() + v, prev_k * a_trace.len() + prev_v))
  //   .collect();
  // let mut permuted_indices = vec![];
  // for w in 0..wire_prev_list.len() {
  //   permuted_indices.push(
  //     wire_prev_list
  //       .iter()
  //       .filter(|&&v| v.0 == w)
  //       .collect::<Vec<_>>()[0]
  //       .1,
  //   );
  // }
  // let permuted_indices = (0..witness_trace.len()).collect::<Vec<_>>();
  // println!("permuted_indices: {:?}", permuted_indices);

  println!("public_first_indices");
  let mut public_first_indices = vec![];
  for w in 0..public_wires.len() {
    if wire_using_list[w].len() > 0 {
      let (k, v) = *wire_using_list[w].first().unwrap();
      public_first_indices.push((w, a_trace_len * k as usize + v));
    }
  }

  assert!(verify_r1cs_proof(
    proof,
    &public_wires,
    &public_first_indices,
    &permuted_indices,
    &coefficients,
    &flag0,
    &flag1,
    &flag2,
    n_constraints,
    n_wires
  )
  .unwrap());
  println!("Done proof verification");
}

pub fn prove_with_file_path<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>>(
  r1cs_file_path: P,
  witness_file_path: Q,
  proof_json_path: R,
) -> Result<(), std::io::Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_bytes(&raw_r1cs);

  let witness_file = File::open(witness_file_path)?;
  let witness_reader = BufReader::new(witness_file);
  let witness: Vec<String> = serde_json::from_reader(witness_reader)?;
  let witness: Vec<Vec<u8>> = witness
    .iter()
    .map(|x| x.parse::<BigUint>().unwrap().to_bytes_le())
    .collect();

  let proof = prove_with_witness(&r1cs, &witness);
  let serialized_proof = serde_json::to_string(&proof)?;
  let mut file = File::create(proof_json_path)?;
  write!(file, "{}", serialized_proof)?;

  // let StarkProof {
  //   m_root: _,
  //   l_root: _,
  //   main_branches,
  //   linear_comb_branches,
  //   fri_proof,
  // } = &proof;
  // let len1 = bin_length(main_branches) + bin_length(linear_comb_branches);
  // let len2 = fri_proof_bin_length(fri_proof);
  // println!(
  //   "Approx proof length: {} (branches), {} (FRI proof), {} (total)",
  //   len1,
  //   len2,
  //   len1 + len2
  // );

  Ok(())
}

pub fn verify_with_file_path<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>>(
  r1cs_file_path: P,
  witness_file_path: Q,
  proof_json_path: R,
) -> Result<(), std::io::Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_bytes(&raw_r1cs);

  let witness_file = File::open(witness_file_path)?;
  let witness_reader = BufReader::new(witness_file);
  let witness: Vec<String> = serde_json::from_reader(witness_reader)?;
  let witness: Vec<Vec<u8>> = witness
    .iter()
    .map(|x| x.parse::<BigUint>().unwrap().to_bytes_le())
    .collect();

  let proof_file = File::open(proof_json_path)?;
  let proof_reader = BufReader::new(proof_file);
  let proof: StarkProof = serde_json::from_reader(proof_reader)?;

  verify_with_witness(&r1cs, &witness, proof);

  Ok(())
}

pub fn run_with_file_path<P: AsRef<Path>, Q: AsRef<Path>, R: AsRef<Path>>(
  r1cs_file_path: P,
  witness_file_path: Q,
  proof_json_path: R,
) -> Result<(), std::io::Error> {
  let r1cs_file = File::open(r1cs_file_path)?;
  let mut raw_r1cs = vec![];
  BufReader::new(r1cs_file).read_to_end(&mut raw_r1cs)?;
  let r1cs = read_bytes(&raw_r1cs);

  // let witness_file = File::open(witness_file_path)?;
  // let witness_reader = BufReader::new(witness_file);
  // let witness: Vec<String> = serde_json::from_reader(witness_reader)?;
  // let witness: Vec<Vec<u8>> = witness
  //   .iter()
  //   .map(|x| x.parse::<BigUint>().unwrap().to_bytes_le())
  //   .collect();

  use bytes::Buf;

  let witness_file = File::open(witness_file_path)?;
  let mut raw_witness = vec![];
  BufReader::new(witness_file).read_to_end(&mut raw_witness)?;
  fn read_witness(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut p = &bytes[..];
    let mut witness = vec![];
    let magic = p.get_u32_le();
    assert_eq!(magic, 1936618615); // wtns
    for _ in 0..5 {
      p.get_u32_le();
    }

    let field_size = p.get_u32_le();
    let mut field_order = BigUint::from(0u32);
    let mut power = BigUint::from(1u32);
    for _ in 0..(field_size / 4) {
      field_order += p.get_u32_le() * &power;
      power *= BigUint::from(1u64 << 32);
    }
    println!("field size: {}", field_order);

    let n_wires = p.get_u32_le();
    p.get_u32_le(); // n_constraints
    p.get_u32_le();
    p.get_u32_le();

    for _ in 0..n_wires {
      let mut field_order = BigUint::from(0u32);
      let mut power = BigUint::from(1u32);
      for _ in 0..(field_size / 4) {
        field_order += p.get_u32_le() * &power;
        power *= BigUint::from(1u64 << 32);
      }

      witness.push(field_order.to_bytes_le().to_vec());
    }

    witness
  }

  let witness = read_witness(&raw_witness);
  // println!("{:?}", witness);

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
  let witness_file_path = "./tests/mul_bn128.wtns";
  let proof_json_path = "./tests/mul_bn128_proof.json";
  prove_with_file_path(r1cs_file_path, witness_file_path, proof_json_path).unwrap();
  verify_with_file_path(r1cs_file_path, witness_file_path, proof_json_path).unwrap();
}
