use crate::multicore::Worker;
use crate::utils::{blake, is_a_power_of_2};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::fmt::Error;

pub trait Element: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {}
pub trait Digest: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {
  fn hash(message: &[u8]) -> Self;
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Proof<E: Element, H: Digest> {
  pub leaf: E,
  pub nodes: Vec<H>,
}

pub trait MerkleTree<E: Element, H: Digest> {
  fn width(&self) -> usize;
  fn root(&self) -> H;
  fn update<I: IntoIterator<Item = E>>(&mut self, into: I);
  fn gen_proof(&self, index: usize) -> Proof<E, H>;
}

pub struct SerialMerkleTree<E: Element, H: Digest> {
  leaves: Vec<E>,
  nodes: Vec<Vec<H>>,
}

pub struct PermutedSerialMerkleTree<E: Element, H: Digest> {
  pub leaves: Vec<E>,
  pub nodes: Vec<Vec<H>>,
}

#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct BlakeDigest(Vec<u8>);

impl AsRef<[u8]> for BlakeDigest {
  fn as_ref(&self) -> &[u8] {
    &(self.0)
  }
}

impl Default for BlakeDigest {
  fn default() -> Self {
    BlakeDigest(vec![])
  }
}

impl Digest for BlakeDigest {
  fn hash(message: &[u8]) -> Self {
    BlakeDigest(blake(message))
  }
}

impl Element for [u8; 32] {}

impl Element for Vec<u8> {}

impl<E: Element, H: Digest> SerialMerkleTree<E, H> {
  pub fn new() -> Self {
    let leaves = vec![];
    let nodes = vec![];
    SerialMerkleTree { leaves, nodes }
  }
}

impl<E: Element, H: Digest> MerkleTree<E, H> for SerialMerkleTree<E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn root(&self) -> H {
    self.nodes.last().unwrap()[0].clone()
  }

  fn update<I: IntoIterator<Item = E>>(&mut self, leaves: I) {
    let leaves: Vec<E> = leaves.into_iter().collect::<Vec<E>>();
    // println!("leaves: {:?}", leaves);
    assert!(is_a_power_of_2(leaves.len()));

    let mut nodes: Vec<Vec<H>> = vec![];
    nodes.push(
      leaves
        .iter()
        .map(|message| H::hash(message.as_ref()))
        .collect(),
    );
    let mut current_nodes = nodes.last().unwrap();

    while current_nodes.len() >= 2 {
      let next_nodes: Vec<H> = current_nodes
        .chunks(2)
        .map(|node| {
          let mut message: Vec<u8> = vec![];
          message.extend(node[0].as_ref());
          message.extend(node[1].as_ref());
          let hash = H::hash(&message);
          // println!("node: {:?}", node[0]);
          // println!("node: {:?}", node[1]);
          // println!("hash: {:?}", hash);
          hash
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
    let index = permute4_index(index, self.width());
    let mut tmp = index;
    let mut proof_nodes = vec![];

    for row in self.nodes.iter().take(self.nodes.len() - 1) {
      proof_nodes.push(row[tmp ^ 1].clone());
      tmp /= 2;
    }
    let proof_leaf = self.leaves[index].clone();
    Proof {
      leaf: proof_leaf,
      nodes: proof_nodes,
    }
  }
}

impl<E: Element, H: Digest> PermutedSerialMerkleTree<E, H> {
  pub fn new() -> Self {
    let leaves = vec![];
    let nodes = vec![];
    PermutedSerialMerkleTree { leaves, nodes }
  }
}

impl<E: Element, H: Digest> MerkleTree<E, H> for PermutedSerialMerkleTree<E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn root(&self) -> H {
    self.nodes.last().unwrap()[0].clone()
  }

  fn update<I: IntoIterator<Item = E>>(&mut self, leaves: I) {
    let leaves: Vec<E> = permute4_values(&leaves.into_iter().collect::<Vec<E>>());
    // println!("leaves: {:?}", leaves);
    assert!(is_a_power_of_2(leaves.len()));

    let mut nodes: Vec<Vec<H>> = vec![];
    nodes.push(
      leaves
        .iter()
        .map(|message| H::hash(message.as_ref()))
        .collect(),
    );
    let mut current_nodes = nodes.last().unwrap();

    while current_nodes.len() >= 2 {
      let next_nodes: Vec<H> = current_nodes
        .chunks(2)
        .map(|node| {
          let mut message: Vec<u8> = vec![];
          message.extend(node[0].as_ref());
          message.extend(node[1].as_ref());
          let hash = H::hash(&message);
          // println!("node: {:?}", node[0]);
          // println!("node: {:?}", node[1]);
          // println!("hash: {:?}", hash);
          hash
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
    let index = permute4_index(index, self.width());
    let mut tmp = index;
    let mut proof_nodes = vec![];

    for row in self.nodes.iter().take(self.nodes.len() - 1) {
      proof_nodes.push(row[tmp ^ 1].clone());
      tmp /= 2;
    }
    let proof_leaf = self.leaves[index].clone();
    Proof {
      leaf: proof_leaf,
      nodes: proof_nodes,
    }
  }
}

#[test]
fn test_serial_single_proof() {
  let index = 2;
  let leaves: Vec<Vec<u8>> = vec![
    hex::decode("7fffffff").unwrap(),
    hex::decode("80000000").unwrap(),
    hex::decode("00000003").unwrap(),
    hex::decode("00000000").unwrap(),
    hex::decode("7ffffffe").unwrap(),
    hex::decode("80000001").unwrap(),
    hex::decode("00000004").unwrap(),
    hex::decode("00000001").unwrap(),
    hex::decode("7ffffffd").unwrap(),
    hex::decode("80000002").unwrap(),
    hex::decode("00000005").unwrap(),
    hex::decode("00000002").unwrap(),
    hex::decode("7ffffffc").unwrap(),
    hex::decode("80000003").unwrap(),
    hex::decode("00000006").unwrap(),
    hex::decode("00000003").unwrap(),
  ];

  let mut merkle_tree: PermutedSerialMerkleTree<Vec<u8>, BlakeDigest> =
    PermutedSerialMerkleTree::new();
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.root();
  assert_eq!(
    merkle_root,
    BlakeDigest(vec![
      40, 219, 9, 138, 84, 182, 133, 203, 147, 236, 201, 244, 138, 114, 92, 232, 66, 52, 230, 57,
      190, 63, 198, 159, 103, 7, 104, 125, 173, 114, 29, 194
    ])
  );
  let proof = merkle_tree.gen_proof(index);
  assert_eq!(proof.leaf, [0, 0, 0, 3]);
  assert_eq!(
    proof.nodes,
    [
      BlakeDigest(vec![
        166, 235, 119, 119, 238, 181, 180, 76, 197, 98, 65, 89, 4, 181, 184, 44, 184, 141, 150,
        152, 150, 19, 15, 61, 138, 104, 19, 45, 150, 128, 112, 202
      ]),
      BlakeDigest(vec![
        23, 222, 231, 110, 6, 40, 238, 204, 106, 124, 246, 240, 22, 170, 28, 53, 124, 233, 113, 14,
        47, 251, 9, 126, 0, 2, 107, 102, 119, 12, 163, 143
      ]),
      BlakeDigest(vec![
        94, 199, 39, 51, 110, 252, 14, 162, 17, 167, 228, 147, 45, 179, 116, 52, 119, 69, 179, 80,
        23, 228, 196, 8, 177, 135, 16, 167, 114, 173, 247, 170
      ]),
      BlakeDigest(vec![
        193, 66, 115, 166, 190, 216, 172, 122, 147, 67, 172, 142, 20, 70, 83, 77, 69, 185, 123, 86,
        41, 116, 123, 29, 34, 212, 80, 20, 18, 22, 223, 185
      ])
    ]
  );

  assert_eq!(proof.validate(&merkle_root, index), Ok(vec![0, 0, 0, 3]));
}

pub struct PermutedParallelMerkleTree<'a, E: Element, H: Digest> {
  worker: &'a Worker,
  pub leaves: Vec<E>,
  pub nodes: Vec<Vec<H>>,
}

impl<'a, E: Element, H: Digest> PermutedParallelMerkleTree<'a, E, H> {
  pub fn new(worker: &'a Worker) -> Self {
    let leaves = vec![];
    let nodes = vec![];
    PermutedParallelMerkleTree {
      worker,
      leaves,
      nodes,
    }
  }
}

impl<'a, E: Element, H: Digest> MerkleTree<E, H> for PermutedParallelMerkleTree<'a, E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn root(&self) -> H {
    self.nodes.last().unwrap()[0].clone()
  }

  fn update<I: IntoIterator<Item = E>>(&mut self, leaves: I) {
    let leaves: Vec<E> = permute4_values(&leaves.into_iter().collect::<Vec<E>>());
    // println!("leaves: {:?}", leaves);
    assert!(is_a_power_of_2(leaves.len()));

    let mut nodes: Vec<Vec<H>> = vec![];

    let mut split_nodes: Vec<Vec<Vec<H>>> = vec![vec![]; self.worker.cpus];
    let mut last_nodes_len = leaves.len();
    while last_nodes_len >= self.worker.cpus {
      last_nodes_len /= 2;
      for i in 0..self.worker.cpus {
        split_nodes[i].push(vec![]);
      }
      nodes.push(vec![]);
    }

    let leaves_chunk_size = leaves.len() / self.worker.cpus;
    let split_nodes_chunk_size = split_nodes.len() / self.worker.cpus;
    self.worker.scope(self.worker.cpus, |scope, _| {
      for (_, (sub_leaves, wrapped_nodes)) in leaves
        .chunks(leaves_chunk_size)
        .zip(split_nodes.chunks_mut(split_nodes_chunk_size))
        .enumerate()
      {
        scope.spawn(move |_| {
          let mut sub_tree = SerialMerkleTree::<E, H>::new();
          sub_tree.update(sub_leaves.to_vec());
          for (i, sub_nodes) in sub_tree.nodes.iter().enumerate() {
            wrapped_nodes[0][i] = sub_nodes.to_vec();
          }
        });
      }
    });

    for wrapped_nodes in split_nodes {
      // flatten
      for (i, v) in wrapped_nodes.iter().enumerate() {
        // println!("v.clone(): {} {:?}", i, v.clone().len());
        nodes[i].extend(v.clone());
      }
    }

    let mut current_nodes = nodes.last().unwrap();
    while current_nodes.len() >= 2 {
      let next_nodes: Vec<H> = current_nodes
        .chunks(2)
        .map(|node| {
          let mut message: Vec<u8> = vec![];
          message.extend(node[0].as_ref());
          message.extend(node[1].as_ref());
          let hash = H::hash(&message);
          // println!("node: {:?}", node[0]);
          // println!("node: {:?}", node[1]);
          // println!("hash: {:?}", hash);
          hash
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
    let index = permute4_index(index, self.width());
    let mut tmp = index;
    let mut proof_nodes = vec![];
    for row in self.nodes.iter().take(self.nodes.len() - 1) {
      proof_nodes.push(row[tmp ^ 1].clone());
      tmp /= 2;
    }
    let proof_leaf = self.leaves[index].clone();
    Proof {
      leaf: proof_leaf,
      nodes: proof_nodes,
    }
  }
}

#[test]
fn test_parallel_single_proof() {
  let index = 2;
  let leaves: Vec<Vec<u8>> = vec![
    hex::decode("7fffffff").unwrap(),
    hex::decode("80000000").unwrap(),
    hex::decode("00000003").unwrap(),
    hex::decode("00000000").unwrap(),
    hex::decode("7ffffffe").unwrap(),
    hex::decode("80000001").unwrap(),
    hex::decode("00000004").unwrap(),
    hex::decode("00000001").unwrap(),
    hex::decode("7ffffffd").unwrap(),
    hex::decode("80000002").unwrap(),
    hex::decode("00000005").unwrap(),
    hex::decode("00000002").unwrap(),
    hex::decode("7ffffffc").unwrap(),
    hex::decode("80000003").unwrap(),
    hex::decode("00000006").unwrap(),
    hex::decode("00000003").unwrap(),
  ];

  let worker = Worker::new();
  let mut merkle_tree = PermutedParallelMerkleTree::<'_, Vec<u8>, BlakeDigest>::new(&worker);
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.root();
  assert_eq!(
    merkle_root,
    BlakeDigest(vec![
      40, 219, 9, 138, 84, 182, 133, 203, 147, 236, 201, 244, 138, 114, 92, 232, 66, 52, 230, 57,
      190, 63, 198, 159, 103, 7, 104, 125, 173, 114, 29, 194
    ])
  );
  let proof = merkle_tree.gen_proof(index);
  assert_eq!(proof.leaf, [0, 0, 0, 3]);
  assert_eq!(
    proof.nodes,
    [
      BlakeDigest(vec![
        166, 235, 119, 119, 238, 181, 180, 76, 197, 98, 65, 89, 4, 181, 184, 44, 184, 141, 150,
        152, 150, 19, 15, 61, 138, 104, 19, 45, 150, 128, 112, 202
      ]),
      BlakeDigest(vec![
        23, 222, 231, 110, 6, 40, 238, 204, 106, 124, 246, 240, 22, 170, 28, 53, 124, 233, 113, 14,
        47, 251, 9, 126, 0, 2, 107, 102, 119, 12, 163, 143
      ]),
      BlakeDigest(vec![
        94, 199, 39, 51, 110, 252, 14, 162, 17, 167, 228, 147, 45, 179, 116, 52, 119, 69, 179, 80,
        23, 228, 196, 8, 177, 135, 16, 167, 114, 173, 247, 170
      ]),
      BlakeDigest(vec![
        193, 66, 115, 166, 190, 216, 172, 122, 147, 67, 172, 142, 20, 70, 83, 77, 69, 185, 123, 86,
        41, 116, 123, 29, 34, 212, 80, 20, 18, 22, 223, 185
      ])
    ]
  );

  assert_eq!(proof.validate(&merkle_root, index), Ok(vec![0, 0, 0, 3]));
}

impl<E: Element, H: Digest> Proof<E, H> {
  fn height(&self) -> usize {
    self.nodes.len()
  }

  pub fn validate(&self, root: &H, index: usize) -> Result<E, Error> {
    let num_of_leaves = 2usize.pow(self.height() as u32);
    let index = permute4_index(index, num_of_leaves);
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
    assert_eq!(current_node, root.clone());
    Ok(self.leaf.clone())
  }
}

pub fn mk_multi_branch<E: Element, H: Digest, T: MerkleTree<E, H>>(
  tree: &T,
  indices: &[usize],
) -> Vec<Proof<E, H>> {
  // println!("indices: {:?}", indices);
  let mut proofs = vec![];
  for &index in indices {
    let proof = tree.gen_proof(index);
    proofs.push(proof);
  }

  proofs
}

pub fn verify_multi_branch<E: Element, H: Digest>(
  root: &H,
  indices: &[usize],
  proofs: Vec<Proof<E, H>>,
) -> Result<Vec<E>, Error> {
  // println!("indices: {:?}", indices);
  Ok(
    indices
      .iter()
      .zip(proofs)
      .map(|(index, proof)| proof.validate(&root, *index).expect(""))
      .collect(),
  )
}

#[test]
fn test_serial_multi_proof() {
  let indices = [2, 7, 13];
  let leaves: Vec<Vec<u8>> = vec![hex::decode("7fffffff").unwrap(); 1 << 20];
  println!("generate Merkle tree");

  let mut merkle_tree: PermutedSerialMerkleTree<Vec<u8>, BlakeDigest> =
    PermutedSerialMerkleTree::new();
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.root();
  assert_eq!(
    merkle_root,
    BlakeDigest(vec![
      174, 144, 67, 84, 246, 248, 226, 220, 207, 80, 125, 76, 79, 127, 103, 100, 40, 133, 31, 137,
      14, 116, 212, 75, 41, 10, 138, 231, 20, 15, 60, 218
    ])
  );

  println!("generate Merkle proof");
  let proofs = mk_multi_branch(&merkle_tree, &indices);
  assert_eq!(proofs[0].leaf, [127, 255, 255, 255]);
  assert_eq!(
    proofs[0].nodes[0],
    BlakeDigest(vec![
      183, 43, 83, 113, 206, 255, 164, 224, 26, 161, 132, 156, 219, 135, 5, 64, 110, 20, 121, 29,
      179, 89, 248, 38, 188, 1, 163, 146, 237, 38, 182, 185
    ])
  );

  verify_multi_branch(&merkle_root, &indices, proofs).unwrap();
  // assert_eq!(Ok(vec![0, 0, 0, 3]));
}

#[test]
fn test_parallel_multi_proof() {
  let indices = [2, 7, 13];
  let leaves: Vec<Vec<u8>> = vec![hex::decode("7fffffff").unwrap(); 1 << 20];
  println!("generate Merkle tree");

  let worker = Worker::new();
  let mut merkle_tree: PermutedParallelMerkleTree<Vec<u8>, BlakeDigest> =
    PermutedParallelMerkleTree::new(&worker);
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.root();
  assert_eq!(
    merkle_root,
    BlakeDigest(vec![
      174, 144, 67, 84, 246, 248, 226, 220, 207, 80, 125, 76, 79, 127, 103, 100, 40, 133, 31, 137,
      14, 116, 212, 75, 41, 10, 138, 231, 20, 15, 60, 218
    ])
  );

  println!("generate Merkle proof");
  let proofs = mk_multi_branch(&merkle_tree, &indices);
  assert_eq!(proofs[0].leaf, [127, 255, 255, 255]);
  assert_eq!(
    proofs[0].nodes[0],
    BlakeDigest(vec![
      183, 43, 83, 113, 206, 255, 164, 224, 26, 161, 132, 156, 219, 135, 5, 64, 110, 20, 121, 29,
      179, 89, 248, 38, 188, 1, 163, 146, 237, 38, 182, 185
    ])
  );

  verify_multi_branch(&merkle_root, &indices, proofs).unwrap();
  // assert_eq!(Ok(vec![0, 0, 0, 3]));
}

pub fn permute4_values<E: Element>(values: &[E]) -> Vec<E> {
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

pub fn bin_length(proof: &[Vec<Vec<u8>>]) -> usize {
  proof.len() * 2
    + proof
      .iter()
      .map(|xs| xs.iter().fold(0, |acc, x| acc + x.len()) + xs.len() / 8)
      .fold(0, |acc, val| acc + val)
}
