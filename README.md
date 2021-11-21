# zk-STARKs Library

[original implementation](https://github.com/ethereum/research/tree/master/mimc_stark)

## How to Use

- [r1cs-stark](packages/r1cs-stark/README.md)

- [r1cs-stark-node-api](packages/napi/README.md)

```sh
cd packages/r1cs-stark
CIRCUIT_NAME="compute"
RUST_BACKTRACE=1 cargo run ./tests/${CIRCUIT_NAME}.r1cs ./tests/${CIRCUIT_NAME}.wtns ./tests/${CIRCUIT_NAME}_proof.json
```

The value of `CIRCUIT_NAME` allows `compute`, `bits`, `pedersen_test` or `sha256_2_test`.

## Test with Docker

If your machine is not Linux, we recommend to use this method to measure memory load.

### build

```sh
docker build -t r1cs-stark .
```

### run

```sh
docker run -it -v $PWD:/root r1cs-stark \
  sh -c 'cd packages/r1cs-stark; CIRCUIT_NAME="compute"; RUST_BACKTRACE=1 cargo run ./tests/${CIRCUIT_NAME}.r1cs ./tests/${CIRCUIT_NAME}.wtns ./tests/${CIRCUIT_NAME}_proof.json > ./tests/log_${CIRCUIT_NAME}.txt'
```

The value of `CIRCUIT_NAME` allows `compute`, `bits` or `pedersen_test`.
