use ff::PrimeField;
use std::convert::TryInto;
use crate::utils::{blake, is_a_power_of_2};

pub fn merklize(nodes: Vec<String>) -> Vec<String> {
  // let mut nodes: Vec<String> = leaves;
  let n = nodes.len();
  assert!(is_a_power_of_2(n.try_into().unwrap()));

  let mut nodes_clone = nodes.clone();

  for i in 0..(n - 1) {
    let message = nodes_clone[2 * i].to_string() + &nodes_clone[2 * i + 1];
    let internal_hash = blake(&message.as_ref());
    nodes_clone.push(internal_hash);
  }
  nodes_clone
}

#[test]
fn test_merklize() {
  use crate::fp::Fp;
  let three = Fp::from(3);
  let three_repr = three.to_repr();
  println!("{:?}", three_repr);
  let leaves = vec![i32::max_value(), i32::min_value(), 3, 0];
  let merkle_tree = merklize(leaves.iter().map(|x| format!("{:08x}", x)).collect());
  assert_eq!(
    merkle_tree,
    [
      "7fffffff",
      "80000000",
      "00000003",
      "00000000",
      "bf873d9fd14913779d20856a11379a595536cfa8f447c6cd36d6e68c827e0547",
      "9d5e04de8033aaeb39a13b8f7f8b435844b3abfb3beffaa1f6eaf8e5df35e87a",
      "649cf0b90a7a1f7eadf5e63ce52e47dc8ec06dedcacd47a89015a32f95780d50"
    ]
  );

  let leaves = vec![1];
  let merkle_tree = merklize(leaves.iter().map(|x| format!("{:08x}", x)).collect());
  assert_eq!(merkle_tree, ["00000001"]);
}

pub fn get_root(tree: &[String]) -> String {
  tree[tree.len() - 1].clone()
}

pub fn mk_branch(tree: &[String], index: usize) -> Vec<String> {
  let mut tree_clone = tree.to_vec();
  tree_clone.push(format!("{:08x}", 0));
  tree_clone.reverse();
  let mut index_clone = (tree.len() + 1) - 1 - index.clone();
  let mut o = vec![tree_clone[index_clone].clone()];
  while index_clone > 1 {
    o.push(tree_clone[index_clone ^ 1].clone());
    index_clone /= 2;
  }
  o
}

pub fn verify_branch(root: &String, index: usize, proof: Vec<String>) -> String {
  let proof_size: usize = 2u32
    .pow(proof.len().try_into().unwrap())
    .try_into()
    .unwrap();
  let mut index_clone: usize = index.clone() + proof_size;
  let mut v = proof[0].clone();
  for i in 1..(proof.len()) {
    let node = &proof[i];
    let message = if index_clone % 2 == 1 {
      node.to_string() + &v
    } else {
      v.to_string() + &node
    };
    v = blake(&message.as_ref());
    index_clone /= 2;
  }
  assert!(v.to_string() == *root);
  proof[0].clone()
}

#[test]
fn test_single_proof() {
  let index = 2;
  let leaves = vec![i32::max_value(), i32::min_value(), 3, 0];
  let merkle_tree = merklize(leaves.iter().map(|x| format!("{:08x}", x)).collect());
  let merkle_root = get_root(&merkle_tree);
  let proof = mk_branch(&merkle_tree, index);
  assert_eq!(
    proof,
    [
      "00000003",
      "00000000",
      "bf873d9fd14913779d20856a11379a595536cfa8f447c6cd36d6e68c827e0547"
    ]
  );

  let res = verify_branch(&merkle_root, index, proof);
  assert_eq!(res, format!("{:08x}", leaves[index]));
}

pub fn mk_multi_branch(tree: &[String], indices: &[usize]) -> Vec<Vec<String>> {
  let mut proofs = vec![];
  for i in indices {
    let proof = mk_branch(tree, *i);
    proofs.push(proof);
  }

  proofs
}

pub fn verify_multi_branch(
  root: &String,
  indices: &[usize],
  proofs: Vec<Vec<String>>,
) -> Vec<String> {
  (0..indices.len())
    .map(|key| {
      let i = indices[key];
      let b = proofs[key].clone();
      verify_branch(root, i, b)
    })
    .collect()
}

#[test]
fn test_multi_proof() {
  let indices = [1, 2];
  let inputs = vec![i32::max_value(), i32::min_value(), 3, 0];
  let leaves: Vec<String> = inputs.iter().map(|x| format!("{:08x}", x)).collect();
  let merkle_tree = merklize(leaves);
  let merkle_root = get_root(&merkle_tree);
  let proofs = mk_multi_branch(&merkle_tree, &indices);
  assert_eq!(
    proofs,
    [
      [
        &leaves[2],
        &leaves[3],
        "bf873d9fd14913779d20856a11379a595536cfa8f447c6cd36d6e68c827e0547"
      ],
      [&leaves[1], &leaves[0], ""]
    ]
  );

  let res = verify_multi_branch(&merkle_root, &indices, proofs);
  let answer: Vec<String> = indices.iter().map(|index| leaves[*index]).collect();
  assert_eq!(res, answer);
}

fn bin_length(proof: Vec<Vec<String>>) -> usize {
  proof.len() * 2
    + proof
      .iter()
      .map(|xs| xs.iter().fold(0, |acc, x| acc + x.len()) + xs.len() / 8)
      .fold(0, |acc, val| acc + val)
}
