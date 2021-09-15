use crate::multicore::Worker;
use crate::utils::{blake, is_a_power_of_2};
use std::fmt::Debug;
use std::fmt::Error;

pub trait Element: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {}
pub trait Digest: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {
  fn hash(message: &[u8]) -> Self;
}

pub trait MerkleTree<E: Element, H: Digest> {
  fn update<I: IntoIterator<Item = E>>(&mut self, into: I);
  fn gen_proof(&self, index: usize) -> Proof<E, H>;
}

pub struct SerialMerkleTree<E: Element, H: Digest> {
  leaves: Vec<E>,
  nodes: Vec<Vec<H>>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Blake2s(Vec<u8>);

impl AsRef<[u8]> for Blake2s {
  fn as_ref(&self) -> &[u8] {
    &(self.0)
  }
}

impl Default for Blake2s {
  fn default() -> Self {
    Blake2s(vec![])
  }
}

impl Digest for Blake2s {
  fn hash(message: &[u8]) -> Self {
    Blake2s(blake(message))
  }
}

impl<E: Element, H: Digest> SerialMerkleTree<E, H> {
  fn new() -> Self {
    let leaves = vec![];
    let nodes = vec![];
    SerialMerkleTree { leaves, nodes }
  }
}

impl<E: Element, H: Digest> MerkleTree<E, H> for SerialMerkleTree<E, H> {
  fn update<I: IntoIterator<Item = E>>(&mut self, into: I) {
    let leaves: Vec<E> = into.into_iter().collect();
    assert!(is_a_power_of_2(leaves.len()));
    let mut nodes: Vec<Vec<H>> = vec![];
    nodes.push(
      leaves
        .iter()
        .map(|message| H::hash(message.as_ref()))
        .collect(),
    );
    let mut current_nodes = nodes.last().unwrap();
    while current_nodes.len() < 2 {
      let next_nodes: Vec<H> = current_nodes
        .chunks(2)
        .map(|node| {
          let mut message: Vec<u8> = vec![];
          message.extend(node[0].as_ref());
          message.extend(node[1].as_ref());
          H::hash(&message)
        })
        .collect();
      nodes.push(next_nodes);
      current_nodes = nodes.last().unwrap();
    }
    self.leaves = leaves;
    self.nodes = nodes;
  }

  // Returns the nodes that need to verify Merkle proof without root and leaf.
  fn gen_proof(&self, index: usize) -> Proof<E, H> {
    let mut proof_nodes = vec![];
    let mut tmp = index;
    for row in self.nodes.iter().take(self.nodes.len() - 1) {
      proof_nodes.push(row[tmp].clone());
      tmp /= 2;
    }
    let proof_leaf = self.leaves[index].clone();
    Proof {
      leaf: proof_leaf,
      nodes: proof_nodes,
    }
  }
}

pub struct ParallelMerkleTree<E: Element, H: Digest> {
  worker: Worker,
  leaves: Vec<E>,
  nodes: Vec<Vec<H>>,
}

impl<E: Element, H: Digest> ParallelMerkleTree<E, H> {
  fn new(worker: Worker) -> Self {
    let leaves = vec![];
    let nodes = vec![];
    ParallelMerkleTree {
      worker,
      leaves,
      nodes,
    }
  }
}

impl<E: Element, H: Digest> MerkleTree<E, H> for ParallelMerkleTree<E, H> {
  fn update<I: IntoIterator<Item = E>>(&mut self, into: I) {
    let leaves: Vec<E> = into.into_iter().collect();
    assert!(is_a_power_of_2(leaves.len()));

    let mut nodes: Vec<Vec<H>> = vec![];
    let mut split_nodes: Vec<Vec<Vec<H>>> = vec![];
    let mut last_nodes_len = nodes.last().unwrap().len();
    while last_nodes_len <= self.worker.cpus {
      last_nodes_len /= 2;
      let next_nodes: Vec<H> = vec![H::default(); last_nodes_len];
      nodes.push(next_nodes);
    }

    self.worker.scope(leaves.len(), |scope, chunk| {
      for (idx, (sub_leaves, wrapped_nodes)) in leaves
        .chunks(chunk)
        .zip(split_nodes.chunks_mut(chunk))
        .enumerate()
      {
        scope.spawn(move |_| {
          let mut sub_tree = SerialMerkleTree::<E, H>::new();
          sub_tree.update(sub_leaves.to_vec());
          for (i, sub_nodes) in sub_tree.nodes.iter().enumerate() {
            wrapped_nodes[i][idx] = sub_nodes.to_vec();
          }
        });
      }
    });

    let mut current_nodes = nodes.last().unwrap();
    while current_nodes.len() <= 1 {
      let next_nodes: Vec<H> = current_nodes
        .chunks(2)
        .map(|node| {
          let mut message: Vec<u8> = vec![];
          message.extend(node[0].as_ref());
          message.extend(node[1].as_ref());
          H::hash(&message)
        })
        .collect();
      nodes.push(next_nodes);
      current_nodes = nodes.last().unwrap();
    }

    self.leaves = leaves;
    self.nodes = nodes;
  }

  // Returns the nodes that need to verify Merkle proof without root and leaf.
  fn gen_proof(&self, index: usize) -> Proof<E, H> {
    let mut proof_nodes = vec![];
    let mut tmp = index;
    for row in self.nodes.iter().take(self.nodes.len() - 1) {
      proof_nodes.push(row[tmp].clone());
      tmp /= 2;
    }
    let proof_leaf = self.leaves[index].clone();
    Proof {
      leaf: proof_leaf,
      nodes: proof_nodes,
    }
  }
}

pub struct Proof<E: Element, H: Digest> {
  leaf: E,
  nodes: Vec<H>,
}

impl<E: Element, H: Digest> Proof<E, H> {
  fn validate(&self, root: H, index: usize) -> Result<E, Error> {
    let mut tmp = index;
    let mut current_node = H::hash(self.leaf.as_ref());
    for node in self.nodes.iter() {
      let mut message: Vec<u8> = vec![];
      if tmp % 2 == 0 {
        message.extend(current_node.as_ref());
        message.extend(node.as_ref());
      } else {
        message.extend(node.as_ref());
        message.extend(current_node.as_ref());
      }
      current_node = H::hash(&message);
      tmp /= 2;
    }
    assert_eq!(current_node, root);
    Ok(self.leaf.clone())
  }
}

pub fn mk_multi_branch<E: Element, H: Digest, T: MerkleTree<E, H>>(
  tree: &T,
  indices: &[usize],
) -> Vec<Proof<E, H>> {
  let mut proofs = vec![];
  for index in indices {
    let proof = tree.gen_proof(*index);
    proofs.push(proof);
  }

  proofs
}

pub fn verify_multi_branch<E: Element, H: Digest>(
  root: H,
  indices: &[usize],
  proofs: &[Proof<E, H>],
) -> Result<Vec<E>, Error> {
  Ok(
    indices
      .iter()
      .zip(proofs)
      .map(|(index, proof)| proof.validate(root.clone(), *index).expect(""))
      .collect(),
  )
}

#[test]
fn test_multi_proof() {
  let indices = [1, 2];
  let leaves: Vec<Vec<u8>> = vec![
    hex::decode("7fffffff").unwrap(),
    hex::decode("80000000").unwrap(),
    hex::decode("00000003").unwrap(),
    hex::decode("00000000").unwrap(),
  ];
  let merkle_tree = merklize(&leaves);
  let merkle_root = get_root(&merkle_tree);
  let proofs = mk_multi_branch(&merkle_tree, &indices);
  assert_eq!(
    proofs
      .iter()
      .map(|proof| proof
        .iter()
        .map(|node| hex::encode(node))
        .collect::<Vec<String>>())
      .collect::<Vec<Vec<String>>>(),
    [
      [
        "80000000",
        "7fffffff",
        "bfc3f121b61735fb3ac096a730b5f52909dc6f76c681f74fa2d59cf54cc4c74d"
      ],
      [
        "00000003",
        "00000000",
        "f086026887af5fd609b58ecc4fec9ad514dba2c6fed57078d1f40ba0b2ecc4ca",
      ],
    ]
  );

  let res = verify_multi_branch(&merkle_root, &indices, &proofs);
  let answer: Vec<Vec<u8>> = indices.iter().map(|index| leaves[*index].clone()).collect();
  assert_eq!(res, answer);
}

pub fn bin_length(proof: &[Vec<Vec<u8>>]) -> usize {
  proof.len() * 2
    + proof
      .iter()
      .map(|xs| xs.iter().fold(0, |acc, x| acc + x.len()) + xs.len() / 8)
      .fold(0, |acc, val| acc + val)
}
