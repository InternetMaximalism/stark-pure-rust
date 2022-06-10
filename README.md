# zk-STARKs Library

[original implementation](https://github.com/ethereum/research/tree/master/mimc_stark)

## How to Use

First, Install [Rust](https://www.rust-lang.org/tools/install) (cargo version >= 1.56.0).

```sh
cd packages/r1cs-stark
rustup override set nightly
```

```sh
CIRCUIT_NAME="compute"
RUST_BACKTRACE=1 cargo run ./tests/${CIRCUIT_NAME}.r1cs ./tests/${CIRCUIT_NAME}.wtns ./tests/${CIRCUIT_NAME}_proof.json
```

`CIRCUIT_NAME` allows following values.

- `compute`
- `bits`
- `pedersen_test`
- `poseidon3_test`
- `sha256_2_test`

## Test with Docker

If your machine is not Linux, we can use this method to measure memory loads.

### build

```sh
docker build -t r1cs-stark .
```

### run

```sh
docker run -it -v $PWD:/root r1cs-stark \
  sh -c 'cd packages/r1cs-stark; CIRCUIT_NAME="compute"; RUST_BACKTRACE=1 cargo run ./tests/${CIRCUIT_NAME}.r1cs ./tests/${CIRCUIT_NAME}.wtns ./tests/${CIRCUIT_NAME}_proof.json > ./tests/log_${CIRCUIT_NAME}.txt'
```

`CIRCUIT_NAME` allows following values.

- `compute`
- `bits`
- `pedersen_test`
- `poseidon3_test`
- <s>`sha256_2_test`</s> (Docker process may be killed due to high memory loads)

## Todo

- Calculate public wires for verifier from input.json .

- Introduce faster algorithm for inv FFT and FFT.

- Reduce the number of constraints of zkSTARK.
