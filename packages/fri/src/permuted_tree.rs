use crate::merkle_tree::{
  bin_length as _bin_length, get_root as _get_root, merklize as _merklize, mk_branch as _mk_branch,
  mk_multi_branch as _mk_multi_branch, verify_branch as _verify_branch,
  verify_multi_branch as _verify_multi_branch,
};
use std::convert::TryInto;

pub fn permute4_values(values: &[Vec<u8>]) -> Vec<Vec<u8>> {
  let mut o: Vec<Vec<u8>> = vec![];
  let ld4 = values.len() / 4;
  for i in 0..ld4 {
    o.extend([
      values[i].clone(),
      values[i + ld4].clone(),
      values[i + ld4 * 2].clone(),
      values[i + ld4 * 3].clone(),
    ]);
  }
  o
}

pub fn permute4_index(index: usize, length: usize) -> usize {
  let ld4 = length / 4;
  index / ld4 + (index % ld4) * 4
}

pub fn permute4_indices(indices: &[usize], length: usize) -> Vec<usize> {
  let ld4 = length / 4;
  indices
    .iter()
    .map(|index| index / ld4 + (index % ld4) * 4)
    .collect()
}

pub fn merklize(values: &[Vec<u8>]) -> Vec<Vec<u8>> {
  _merklize(&permute4_values(values))
}

pub fn get_root(tree: &[Vec<u8>]) -> Vec<u8> {
  _get_root(tree)
}

pub fn mk_branch(tree: &[Vec<u8>], index: usize) -> Vec<Vec<u8>> {
  let num_of_leaves = (tree.len() + 1) / 2;
  _mk_branch(tree, permute4_index(index, num_of_leaves))
}

pub fn verify_branch(root: &[u8], index: usize, proof: &[Vec<u8>]) -> Vec<u8> {
  let num_of_leaves = 2usize.pow(proof.len().try_into().unwrap()) / 2;
  _verify_branch(root, permute4_index(index, num_of_leaves), proof)
}

#[test]
fn test_single_permuted_proof() {
  let index = 2;
  let leaves: Vec<Vec<u8>> = vec![
    hex::decode("7fffffff").unwrap(),
    hex::decode("80000000").unwrap(),
    hex::decode("00000003").unwrap(),
    hex::decode("00000000").unwrap(),
  ];
  let merkle_tree = merklize(&leaves);
  let merkle_root = get_root(&merkle_tree);
  let proof = mk_branch(&merkle_tree, index);
  assert_eq!(
    proof
      .iter()
      .map(|node| hex::encode(node))
      .collect::<Vec<String>>(),
    [
      "00000003",
      "00000000",
      "f086026887af5fd609b58ecc4fec9ad514dba2c6fed57078d1f40ba0b2ecc4ca",
    ]
  );

  let res = verify_branch(&merkle_root, index, &proof);
  assert_eq!(res, leaves[index]);
}

pub fn mk_multi_branch(tree: &[Vec<u8>], indices: &[usize]) -> Vec<Vec<Vec<u8>>> {
  _mk_multi_branch(tree, &permute4_indices(indices, (tree.len() + 1) / 2))
}

pub fn verify_multi_branch(
  root: &[u8],
  indices: &[usize],
  proofs: &[Vec<Vec<u8>>],
) -> Vec<Vec<u8>> {
  let num_of_leaves = 2usize.pow(proofs.first().unwrap().len().try_into().unwrap()) / 2;
  _verify_multi_branch(root, &permute4_indices(indices, num_of_leaves), proofs)
}

pub fn bin_length(proof: &[Vec<Vec<u8>>]) -> usize {
  _bin_length(proof)
}
