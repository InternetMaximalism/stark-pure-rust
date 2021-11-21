use std::time::Duration;

use commitment::blake::BlakeDigest;
use r1cs_stark::run::run_with_file_path;
use self_meter::Meter;

fn main() {
  let mut args = std::env::args().skip(1);
  let r1cs_file_path = args.next().unwrap();
  let witness_file_path = args.next().unwrap();
  let proof_json_path = args.next().unwrap();

  let mut meter = Meter::new(Duration::new(1, 0)).unwrap();
  meter.track_current_thread("main");
  meter.scan().unwrap();
  let result =
    run_with_file_path::<_, _, _, BlakeDigest>(r1cs_file_path, witness_file_path, proof_json_path);
  meter.scan().unwrap();
  let reports = meter.report().unwrap();
  println!("reports: {:#?}", reports);

  result.unwrap();
}
