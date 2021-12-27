use commitment::blake::BlakeDigest;
use r1cs_stark::run::run_with_file_path;

fn main() {
  let mut args = std::env::args().skip(1);
  let r1cs_file_path = args.next().unwrap();
  let witness_file_path = args.next().unwrap();
  let proof_json_path = args.next().unwrap();
  run_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path)
    .unwrap();
}

#[test]
fn test_main_func_to_run_compute() {
  use commitment::blake::BlakeDigest;

  let r1cs_file_path = "./packages/r1cs-stark/tests/compute.r1cs";
  let witness_file_path = "./packages/r1cs-stark//tests/compute.wtns";
  let proof_json_path = "./packages/r1cs-stark//tests/compute_proof.json";
  run_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path)
    .unwrap();
}
