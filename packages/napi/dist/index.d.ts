declare function prove(r1cs_file_path: string, witness_json_path: string, proof_json_path: string): void;
declare function verify(r1cs_file_path: string, witness_json_path: string, proof_json_path: string): any;
declare function prove_with_file_path(r1cs_file_path: string, witness_json_path: string, proof_json_path: string): void;
declare function verify_with_file_path(r1cs_file_path: string, witness_json_path: string, proof_json_path: string): any;
export { prove, verify, prove_with_file_path, verify_with_file_path };
