use r1cs_stark::run::run_with_file_path;

fn main() {
  let mut args = std::env::args().skip(1);
  let r1cs_file_path = args.next().unwrap();
  let witness_json_path = args.next().unwrap();
  let proof_json_path = args.next().unwrap();
  run_with_file_path(r1cs_file_path, witness_json_path, proof_json_path).unwrap();
}
