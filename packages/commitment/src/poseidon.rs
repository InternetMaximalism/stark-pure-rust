use crate::hash::Digest;
use blstrs::Scalar as Fr;
use generic_array::typenum::U2;
use neptune::poseidon::{HashMode, PoseidonConstants};
use neptune::*;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PoseidonDigest(pub Vec<u8>);

impl Debug for PoseidonDigest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("\"0x{}\"", hex::encode(&self.0)))
    }
}

impl AsRef<[u8]> for PoseidonDigest {
    fn as_ref(&self) -> &[u8] {
        &(self.0)
    }
}

impl Default for PoseidonDigest {
    fn default() -> Self {
        Self(vec![])
    }
}

impl Digest for PoseidonDigest {
    fn hash(message: &[u8]) -> Self {
        let message_len = message.len();
        assert!(message_len <= 64);
        let mut padded_message = message.to_vec();
        padded_message.extend(vec![0u8; (((message_len - 1) / 32) + 1) * 32 - message_len]);
        println!("padded_message: {:?}", padded_message);
        let s = padded_message
            .chunks(32)
            .map(|b| {
                Fr::from_bytes_le(&{
                    let mut res = [0u8; 32];
                    for i in 0..b.len() {
                        res[i] = b[i];
                    }
                    res
                })
                .unwrap()
            })
            .collect::<Vec<_>>();
        let constants: PoseidonConstants<Fr, U2> =
            PoseidonConstants::new_with_strength(Strength::Standard);
        let mut h: Poseidon<Fr> = Poseidon::new(&constants);
        h.reset();
        println!("s: {:?}", s);
        for scalar in s {
            h.input(scalar).unwrap();
        }

        let digest = h.hash_in_mode(HashMode::Correct).to_bytes_le();
        println!("digest: {:?}", digest);
        Self(digest.to_vec())
    }
}

#[test]
fn test_poseidon_3() {
    let mut message = (0..3).collect::<Vec<_>>();
    message.extend(vec![0u8; 64 - message.len()]);
    let res = PoseidonDigest::hash(&message);
    assert_eq!(
        hex::encode(res.as_ref()),
        "b3a1a3cfaebc3a557d52dd3e25076f7f7b51f2bf46f5289d66c389b51477ec25"
    );
}

#[test]
fn test_poseidon_32() {
    let mut message = (0..32).collect::<Vec<_>>();
    message.extend(vec![0u8; 64 - message.len()]);
    let res = PoseidonDigest::hash(&message);
    assert_eq!(
        hex::encode(res.as_ref()),
        "0e67a788ec648e60632957f8d10b71f12fba0050a7688bdad9de2e78dbf5495b"
    );
}

#[test]
fn test_poseidon_63() {
    let mut message = (0..63).collect::<Vec<_>>();
    message.extend(vec![0u8; 64 - message.len()]);
    let res = PoseidonDigest::hash(&message);
    assert_eq!(
        hex::encode(res.as_ref()),
        "ddae0004ffee05d6da43777af82faa1f0c6ac08d7048f9a4ddf6d2b259f7075e"
    );
}

#[test]
fn test_poseidon_64() {
    let message = (0..64).collect::<Vec<_>>();
    let res = PoseidonDigest::hash(&message);
    assert_eq!(
        hex::encode(res.as_ref()),
        "93bde2916aec7310f6e07faa70f14ed0c173832adcc03aeaed230f94540f0632"
    );
}

#[test]
#[should_panic]
fn test_poseidon_65() {
    let message = (0..65).collect::<Vec<_>>();
    let _ = PoseidonDigest::hash(&message);
}
