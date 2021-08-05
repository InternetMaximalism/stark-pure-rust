const {
  prove: _prove,
  verify: _verify,
  prove_with_file_path: _prove_with_file_path,
  verify_with_file_path: _verify_with_file_path,
} = require("../index.node");

function prove(
  r1cs_file_path: string,
  witness_json_path: string,
  proof_json_path: string
): void {
  return _prove(r1cs_file_path, witness_json_path, proof_json_path);
}

function verify(
  r1cs_file_path: string,
  witness_json_path: string,
  proof_json_path: string
) {
  return _verify(r1cs_file_path, witness_json_path, proof_json_path);
}

function prove_with_file_path(
  r1cs_file_path: string,
  witness_json_path: string,
  proof_json_path: string
): void {
  return _prove_with_file_path(
    r1cs_file_path,
    witness_json_path,
    proof_json_path
  );
}

function verify_with_file_path(
  r1cs_file_path: string,
  witness_json_path: string,
  proof_json_path: string
) {
  return _verify_with_file_path(
    r1cs_file_path,
    witness_json_path,
    proof_json_path
  );
}

export { prove, verify, prove_with_file_path, verify_with_file_path };
