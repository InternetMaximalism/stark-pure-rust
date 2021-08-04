# R1CS STARK

## How to Use

```sh
cargo run <r1cs_file_path> <witness_json_path> <proof_json_path>
```

For example,

```sh
cargo run ./tests/mul_bn128.r1cs ./tests/mul_bn128_wtns_valid.json ./tests/mul_bn128_proof.json
```

## API

```rust
fn sample_code() -> Result<(), std::io::Error> {
  use r1cs_stark::run::{prove_with_file_path, verify_with_file_path};

  let r1cs_file_path = "./tests/mul_bn128.r1cs";
  let witness_json_path = "./tests/mul_bn128_wtns_valid.json";
  let proof_json_path = "./tests/mul_bn128_proof.json";

  // r1cs + witness -> proof
  prove_with_file_path(r1cs_file_path, witness_json_path, proof_json_path)?;

  // r1cs + public input/output + proof -> verification result
  verify_with_file_path(r1cs_file_path, witness_json_path, proof_json_path)?;
}
```
