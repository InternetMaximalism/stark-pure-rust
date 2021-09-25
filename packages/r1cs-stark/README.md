# R1CS STARK

## How to Use

```sh
cargo run <r1cs_file_path> <witness_file_path> <proof_json_path>
```

For example,

```sh
cargo run ./tests/compute.r1cs ./tests/compute.wtns ./tests/compute_proof.json
```

## API

```rust
fn sample_code() -> Result<(), std::io::Error> {
  use r1cs_stark::run::{prove_with_file_path, verify_with_file_path};

  let r1cs_file_path = "./tests/compute.r1cs";
  let witness_file_path = "./tests/compute.wtns";
  let proof_json_path = "./tests/compute_proof.json";

  // r1cs + witness -> proof
  prove_with_file_path(r1cs_file_path, witness_file_path, proof_json_path)?;

  // r1cs + public input/output + proof -> verification result
  verify_with_file_path(r1cs_file_path, witness_file_path, proof_json_path)?;
}
```

## Format

### R1CS

[iden3 R1CS format](https://github.com/iden3/r1csfile/blob/master/doc/r1cs_bin_format.md)

### Witness

array of decimal string which is less than field size

```json
["1", "12", "3", "2", "4"]
```

## Road Map

- introduce private input
