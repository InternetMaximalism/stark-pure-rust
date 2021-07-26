use crate::merkle_tree::{
  bin_length as _bin_length, merklize as _merklize, mk_branch as _mk_branch,
  mk_multi_branch as _mk_multi_branch, verify_branch as _verify_branch,
  verify_multi_branch as _verify_multi_branch,
};
use std::convert::TryInto;

pub fn permute4_values(values: &[String]) -> Vec<String> {
  let mut o = vec![];
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

pub fn permute4_index(x: usize, value: usize) -> usize {
  let ld4 = value / 4;
  x / ld4 + (x % ld4) * 4
}

pub fn permute4_indices(xs: &[usize], value: usize) -> Vec<usize> {
  let ld4 = value / 4;
  xs.iter().map(|x| x / ld4 + (x % ld4) * 4).collect()
}

pub fn merklize(values: &[String]) -> Vec<String> {
  _merklize(permute4_values(values))
}

pub fn mk_branch(tree: &[String], index: usize) -> Vec<String> {
  let n = permute4_index(index, tree.len() / 2);
  _mk_branch(tree, n)
}

pub fn verify_branch(root: &String, index: usize, proof: &[String]) -> String {
  let proof_size: usize = (2u32.pow(proof.len().try_into().unwrap()) / 2)
    .try_into()
    .unwrap();
  _verify_branch(root, permute4_index(index, proof_size), proof)
}

pub fn mk_multi_branch(tree: &[String], indices: &[usize]) -> Vec<Vec<String>> {
  _mk_multi_branch(tree, &permute4_indices(indices, tree.len() / 2))
}

pub fn verify_multi_branch(
  root: &String,
  indices: &[usize],
  proofs: &[Vec<String>],
) -> Vec<String> {
  let p: u32 = proofs[0].len().try_into().unwrap();
  _verify_multi_branch(root, &permute4_indices(indices, 2usize.pow(p) / 2), proofs)
}

pub fn bin_length(proof: &[Vec<String>]) -> usize {
  _bin_length(proof)
}
