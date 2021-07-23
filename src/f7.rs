use core::iter::FromIterator;
use hex::ToHex;
use ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "7"]
#[PrimeFieldGenerator = "3"]
#[PrimeFieldReprEndianness = "little"]
pub struct F7([u64; 1]);

impl ToHex for F7 {
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
fn test_f7_to_hex() {
  let input = 1;
  let x = F7::from(input);
  assert_eq!(x.encode_hex::<String>(), format!("0x{:016x}", input));
  assert_eq!(x.encode_hex_upper::<String>(), format!("0x{:016X}", input));
}
