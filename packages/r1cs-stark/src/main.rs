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
