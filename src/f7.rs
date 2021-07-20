use ff::{PrimeField};

#[derive(PrimeField)]
#[PrimeFieldModulus = "7"]
#[PrimeFieldGenerator = "3"]
#[PrimeFieldReprEndianness = "little"]
pub struct F7([u64; 1]);
