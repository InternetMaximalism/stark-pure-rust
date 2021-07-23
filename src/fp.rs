use core::iter::FromIterator;
use hex::ToHex;
use ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
#[PrimeFieldReprEndianness = "little"]
pub struct Fp([u64; 4]);

impl ToHex for Fp {
  // Parse a Fp value to a hex string with 0x-prefix.
  fn encode_hex<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.to_repr());
    T::from_iter(repr.chars())
  }

  fn encode_hex_upper<T: FromIterator<char>>(&self) -> T {
    let repr = format!("{:?}", self.to_repr());
    T::from_iter(format!("0x{}", &repr.to_uppercase()[2..]).chars())
  }
}

#[test]
fn test_fp_to_hex() {
  let input = 31;
  let x = Fp::from(input);
  assert_eq!(x.encode_hex::<String>(), format!("0x{:064x}", input));
  assert_eq!(x.encode_hex_upper::<String>(), format!("0x{:064X}", input));
}
