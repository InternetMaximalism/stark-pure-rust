use bytes::Buf;
pub use circom2bellman_core::{
    read_bytes as read_r1cs, Coefficient, Constraints, Factor, Header, R1csContents,
};
use num::bigint::BigUint;

pub fn read_witness(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut p = bytes;
    let mut witness = vec![];
    let magic = p.get_u32_le();
    assert_eq!(magic, 1936618615); // wtns
    for _ in 0..5 {
        p.get_u32_le();
    }

    let field_size = p.get_u32_le();
    let mut field_order = BigUint::from(0u32);
    let mut power = BigUint::from(1u32);
    for _ in 0..(field_size / 4) {
        field_order += p.get_u32_le() * &power;
        power *= BigUint::from(1u64 << 32);
    }
    // println!("field size: {}", field_order);

    let n_wires = p.get_u32_le();
    p.get_u32_le(); // n_constraints
    p.get_u32_le();
    p.get_u32_le();

    for _ in 0..n_wires {
        let mut field_order = BigUint::from(0u32);
        let mut power = BigUint::from(1u32);
        for _ in 0..(field_size / 4) {
            field_order += p.get_u32_le() * &power;
            power *= BigUint::from(1u64 << 32);
        }

        witness.push(field_order.to_bytes_le().to_vec());
    }

    witness
}

#[test]
fn test_read_r1cs() {
    use std::fs::File;
    use std::io::{BufReader, Read};

    let r1cs_file_path = "./tests/compute.r1cs";
    let r1cs_json_path = "./tests/compute.r1cs.json";
    let r1cs_file = File::open(r1cs_file_path).unwrap();
    let mut raw_r1cs = vec![];
    let mut r1cs_reader = BufReader::new(r1cs_file);
    r1cs_reader.read_to_end(&mut raw_r1cs).unwrap();
    let r1cs = read_r1cs(&raw_r1cs);

    let r1cs_json_file = File::open(r1cs_json_path).unwrap();
    let r1cs_json_reader = BufReader::new(r1cs_json_file);
    let answer_r1cs: R1csContents = serde_json::from_reader(r1cs_json_reader).unwrap();

    assert_eq!(r1cs, answer_r1cs);
}

#[test]
fn test_read_witness() {
    use std::fs::File;
    use std::io::{BufReader, Read};

    let witness_file_path = "./tests/compute.wtns";
    let witness_file = File::open(witness_file_path).unwrap();
    let mut raw_witness = vec![];
    BufReader::new(witness_file)
        .read_to_end(&mut raw_witness)
        .unwrap();
    let witness = read_witness(&raw_witness);
    assert_eq!(
        witness,
        vec![
            vec![1],
            vec![
                135, 136, 135, 103, 17, 74, 207, 218, 212, 163, 232, 164, 38, 238, 216, 34, 56,
                221, 180, 135, 36, 249, 144, 247, 19, 79, 126, 26, 164, 114, 177, 5,
            ],
            vec![17],
            vec![33, 1],
            vec![49, 19],
        ],
    );
}
