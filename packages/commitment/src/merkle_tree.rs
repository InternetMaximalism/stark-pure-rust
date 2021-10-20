use crate::delayed::Delayed;
use crate::multicore::Worker;
use crate::utils::{blake, is_a_power_of_2, log2_ceil};
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::fmt::Error;
use std::ops::Deref;

pub trait Element: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {}
pub trait Digest: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {
  fn hash(message: &[u8]) -> Self;
}

// `nodes` is arranged in the direction from the leaves to the roots of the tree structure.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Proof<E: Element, H: Digest> {
  pub leaf: E,
  pub nodes: Vec<H>,
}

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BlakeDigest(Vec<u8>);

impl Debug for BlakeDigest {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    f.write_fmt(format_args!("\"0x{}\"", hex::encode(&self.0)))
  }
}

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

impl<E: Element, H: Digest> Proof<E, H> {
  pub fn height(&self) -> usize {
    self.nodes.len()
  }

  pub fn validate(&self, root: &H, index: usize) -> Result<E, Error> {
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

// This function is the single thread version of gen_multi_proofs_multi_core.
// This implementation is only for testing usage.
pub fn gen_multi_proofs_in_place<E: Element, H: Digest>(
  current_layer: &mut [H],
  indices: &[usize],
  steps: usize,
) -> Vec<Proof<E, H>> {
  // If the overlapping of proof is omitted, the better practice is to permute leaves.
  // Because we pack unused nodes in one place, it is possible to increase the number of overlapping nodes.
  // let leaves: Vec<E> = permute4_values(&leaves);
  // let indices = permute4_indices(indices, leaves.len());
  // println!("leaves: {:?}", leaves);
  // println!("indices: {:?}", indices);
  assert!(is_a_power_of_2(current_layer.len()));

  let mut proofs = vec![
    Proof {
      leaf: E::default(),
      nodes: vec![]
    };
    indices.len()
  ];

  let mut log_nodes_interval = log2_ceil(steps);
  assert_eq!(1 << log_nodes_interval, steps);
  let current_nodes_len = current_layer.len();
  while (1 << log_nodes_interval) < current_nodes_len {
    for (i, &tmp_index) in indices.iter().enumerate() {
      let twin_index = ((tmp_index >> log_nodes_interval) ^ 1) << log_nodes_interval;
      proofs[i].nodes.push(current_layer[twin_index].clone());
    }

    let nodes_interval = 1 << log_nodes_interval;
    let chunks_num = nodes_interval << 1;
    for node in current_layer.chunks_mut(chunks_num) {
      let mut message: Vec<u8> = vec![];
      message.extend(node[0].as_ref());
      message.extend(node[nodes_interval].as_ref());
      let hash = H::hash(&message);
      // println!("node: {:?}", node[0]);
      // println!("node: {:?}", node[1]);
      // println!("hash: {:?}", hash);
      node[0] = hash;
    }

    log_nodes_interval += 1;
  }

  proofs
}

// This function returns the Merkle proofs at multiple indices.
// In making a Merkle tree, this function forget unnecessary nodes
// that are not relevant to the proof.
pub fn gen_multi_proofs_multi_core<E: Element, H: Digest>(
  leaves: &[Delayed<'_, E>],
  indices: &[usize],
  worker: &Worker,
) -> (Vec<Proof<E, H>>, H) {
  // println!("leaves: {:?}", leaves);
  // println!("indices: {:?}", indices);
  assert!(is_a_power_of_2(leaves.len()));

  let start = std::time::Instant::now();
  let mut sorted_indices = indices.iter().enumerate().collect::<Vec<_>>().clone();
  sorted_indices.sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
  let end: std::time::Duration = start.elapsed();
  println!(
    "Sort leaves: {}.{:03}s",
    end.as_secs(),
    end.subsec_nanos() / 1_000_000
  );

  let start = std::time::Instant::now();

  let proof_leaves = sorted_indices
    .iter()
    .map(|(_, &index)| leaves[index].deref().clone())
    .collect::<Vec<_>>();

  let mut current_nodes = leaves
    .iter()
    .map(|message| H::hash(message.as_ref()))
    .collect::<Vec<_>>();

  let end: std::time::Duration = start.elapsed();
  println!(
    "Collect leaves: {}.{:03}s",
    end.as_secs(),
    end.subsec_nanos() / 1_000_000
  );

  let start = std::time::Instant::now();

  let chunks_num = 1 << worker.log_num_cpus();
  let chunk_size = current_nodes.len() / chunks_num;
  let mut next_indices = vec![None; chunks_num];
  let mut sub_proofs = current_nodes
    .chunks_mut(chunk_size)
    .zip(next_indices.chunks_mut(1))
    .enumerate()
    .map(|(chunk_index, (sub_current_nodes, next_index))| {
      // Filter index from chunk_index * chunk_size to (chunk_index + 1) * chunk_size
      // and subtract chunk_index * chunk_size.
      let sub_indices = sorted_indices
        .iter()
        .filter_map(|(_, &index)| {
          if index >= chunk_index * chunk_size && index < (chunk_index + 1) * chunk_size {
            Some(index - chunk_index * chunk_size)
          } else {
            None
          }
        })
        .collect::<Vec<_>>();
      // println!("sub_indices: {:?}", sub_indices);
      if sub_indices.len() != 0 {
        next_index[0] = Some(chunk_index);
      }

      gen_multi_proofs_in_place::<E, H>(sub_current_nodes, &sub_indices, 1)
    })
    .flatten()
    .collect::<Vec<_>>();

  // let next_indices = next_indices.iter().filter_map(|&v| v).collect::<Vec<_>>();
  // println!("next_indices: {:?}", next_indices);

  // The n-th sub-tree root exists the (chunk_size * n)-th elements of current_nodes.
  // println!("chunk_size: {:?}", chunk_size);
  let next_proofs = gen_multi_proofs_in_place::<E, H>(
    &mut current_nodes,
    &(0..chunks_num).map(|v| v * chunk_size).collect::<Vec<_>>(),
    chunk_size,
  );

  let end: std::time::Duration = start.elapsed();
  println!(
    "Calculated tree: {}.{:03}s",
    end.as_secs(),
    end.subsec_nanos() / 1_000_000
  );

  let merkle_root = current_nodes[0].clone();

  for (i, ((_, &index), leaf)) in sorted_indices.iter().zip(proof_leaves).enumerate() {
    let next_index = index / chunk_size;
    sub_proofs[i].leaf = leaf;
    sub_proofs[i]
      .nodes
      .extend(next_proofs[next_index].nodes.clone());
  }

  let mut enhanced_sub_proofs = sub_proofs.iter().zip(sorted_indices).collect::<Vec<_>>();
  enhanced_sub_proofs.sort_by(|(_, (a, _)), (_, (b, _))| a.partial_cmp(b).unwrap());
  let proofs = enhanced_sub_proofs
    .iter()
    .map(|(v, _)| v.clone().clone())
    .collect::<Vec<_>>();
  (proofs, merkle_root)
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
fn test_gen_multi_proofs_multi_core() {
  use crate::lazily;

  // let leaves = vec![
  //   [0, 0, 0, 0],
  //   [0, 0, 0, 1],
  //   ...,
  //   [0, 0, 0, 15],
  // ];
  let leaves: Vec<Vec<u8>> = (0..16)
    .map(|i| hex::decode(format!("{:08x}", i)).unwrap().to_vec())
    .collect::<Vec<_>>();
  // println!("leaves: {:?}", leaves);

  let worker = Worker::new_with_cpus(4);

  let indices = [10, 4, 6, 3, 6, 8];
  let (merkle_proofs, merkle_root) = gen_multi_proofs_multi_core::<Vec<u8>, BlakeDigest>(
    &leaves
      .iter()
      .map(|v| lazily!(v.to_vec()))
      .collect::<Vec<_>>(),
    &indices,
    &worker,
  );
  // println!("{:?} {:?}", merkle_proofs, merkle_root);

  let mut merkle_tree: SerialMerkleTree<Vec<u8>, BlakeDigest> = SerialMerkleTree::new();
  merkle_tree.update(leaves.clone());
  let merkle_root2 = merkle_tree.root().unwrap();
  assert_eq!(merkle_root, merkle_root2, "invalid Merkle root");

  let mut merkle_proofs2 = vec![];
  merkle_proofs2.push(merkle_tree.gen_proof(indices[0]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[1]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[2]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[3]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[4]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[5]));
  // println!("merkle_proofs1: {:?}", merkle_proofs);
  // println!("merkle_proofs2: {:?}", merkle_proofs2);
  for i in 0..indices.len() {
    assert_eq!(
      merkle_proofs[i].leaf, merkle_proofs2[i].leaf,
      "invalid proof leaf: {:?}",
      i
    );
    assert_eq!(
      merkle_proofs[i].nodes, merkle_proofs2[i].nodes,
      "invalid proof nodes: {:?}",
      i
    );
  }

  verify_multi_branch(&merkle_root, &indices, merkle_proofs).unwrap();
}

pub trait MerkleTree<E: Element, H: Digest> {
  // This method returns the number of leaves.
  fn width(&self) -> usize;

  // This method returns the Merkle root.
  // The Merkle root is the first node of last layer.
  fn root(&self) -> Option<H>;

  // This method makes the Merkle tree.
  fn update<I: IntoIterator<Item = E>>(&mut self, into: I);

  // This method returns the nodes that need to verify Merkle proof without root and leaf.
  fn gen_proof(&self, index: usize) -> Proof<E, H>;
}

// `SerialMerkleTree` describe a simple Merkle tree structure.
// This implement is only for testing usage.
// A layer is the set of the tree nodes which are the same depth.
// The first layer is the hash of each leaves
// The last layer is the tree root.
pub struct SerialMerkleTree<E: Element, H: Digest> {
  leaves: Vec<E>,
  layers: Vec<Vec<H>>,
}

impl<E: Element, H: Digest> SerialMerkleTree<E, H> {
  pub fn new() -> Self {
    let leaves = vec![];
    let layers = vec![];
    SerialMerkleTree { leaves, layers }
  }
}

impl<E: Element, H: Digest> MerkleTree<E, H> for SerialMerkleTree<E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn root(&self) -> Option<H> {
    if let Some(n) = self.layers.last() {
      Some(n[0].clone())
    } else {
      None
    }
  }

  fn update<I: IntoIterator<Item = E>>(&mut self, leaves: I) {
    let leaves: Vec<E> = leaves.into_iter().collect::<Vec<E>>();
    // println!("leaves: {:?}", leaves);
    // assert!(is_a_power_of_2(leaves.len()));

    let mut layers: Vec<Vec<H>> = vec![];
    layers.push(
      leaves
        .iter()
        .map(|message| H::hash(message.as_ref()))
        .collect(),
    );
    let mut current_layer = layers.last().unwrap();

    while current_layer.len() >= 2 {
      // the next layer's i-th node is the hash of the the current layer's 2i-th and (2i+1)-th nodes.
      let next_layer: Vec<H> = current_layer
        .chunks(2)
        .map(|node_pair| {
          let mut message: Vec<u8> = vec![];
          message.extend(node_pair[0].as_ref());
          if node_pair.len() > 1 {
            message.extend(node_pair[1].as_ref());
          }
          // println!("node0: {:?}", node_pair[0]);
          // println!("node1: {:?}", node_pair[1]);

          // the hash of node_pair[0] and node_pair[1]
          let hash = H::hash(&message);
          // println!("hash: {:?}", hash);

          hash
        })
        .collect();
      layers.push(next_layer);
      current_layer = layers.last().unwrap();
    }
    self.leaves = leaves;
    self.layers = layers;
  }

  fn gen_proof(&self, index: usize) -> Proof<E, H> {
    let mut tmp = index;
    let mut proof_nodes = vec![];

    for row in self.layers.iter().take(self.layers.len() - 1) {
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

pub struct ParallelMerkleTree<'a, E: Element, H: Digest> {
  worker: &'a Worker,
  pub leaves: Vec<E>,
  pub layers: Vec<Vec<H>>,
}

impl<'a, E: Element, H: Digest> ParallelMerkleTree<'a, E, H> {
  pub fn new(worker: &'a Worker) -> Self {
    let leaves = vec![];
    let layers = vec![];
    Self {
      worker,
      leaves,
      layers,
    }
  }
}

impl<'a, E: Element, H: Digest> MerkleTree<E, H> for ParallelMerkleTree<'a, E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn root(&self) -> Option<H> {
    if let Some(n) = self.layers.last() {
      Some(n[0].clone())
    } else {
      None
    }
  }

  fn update<I: IntoIterator<Item = E>>(&mut self, leaves: I) {
    let leaves: Vec<E> = leaves.into_iter().collect::<Vec<E>>();
    // println!("leaves: {:?}", leaves);
    // assert!(is_a_power_of_2(leaves.len()));

    let mut layers: Vec<Vec<H>> = vec![];

    let mut split_layers: Vec<Vec<Vec<H>>> = vec![vec![]; self.worker.cpus];
    let mut last_layer_len = leaves.len();
    while last_layer_len >= self.worker.cpus {
      last_layer_len /= 2;
      for i in 0..self.worker.cpus {
        split_layers[i].push(vec![]);
      }
      layers.push(vec![]);
    }

    let leaves_chunk_size = leaves.len() / self.worker.cpus;
    let split_layers_chunk_size = split_layers.len() / self.worker.cpus;
    self.worker.scope(self.worker.cpus, |scope, _| {
      for (_, (sub_leaves, wrapped_layer)) in leaves
        .chunks(leaves_chunk_size)
        .zip(split_layers.chunks_mut(split_layers_chunk_size))
        .enumerate()
      {
        scope.spawn(move |_| {
          let mut sub_tree = SerialMerkleTree::<E, H>::new();
          sub_tree.update(sub_leaves.to_vec());
          for (i, sub_layer) in sub_tree.layers.iter().enumerate() {
            wrapped_layer[0][i] = sub_layer.to_vec();
          }
        });
      }
    });

    for wrapped_layer in split_layers {
      // flatten
      for (i, v) in wrapped_layer.iter().enumerate() {
        // println!("v.clone(): {} {:?}", i, v.clone().len());
        layers[i].extend(v.clone());
      }
    }

    let mut current_layer = layers.last().unwrap();
    while current_layer.len() >= 2 {
      let next_layer: Vec<H> = current_layer
        .chunks(2)
        .map(|node_pair| {
          let mut message: Vec<u8> = vec![];
          message.extend(node_pair[0].as_ref());
          if node_pair.len() > 1 {
            message.extend(node_pair[1].as_ref());
          }
          // println!("node0: {:?}", node_pair[0]);
          // println!("node1: {:?}", node_pair[1]);

          let hash = H::hash(&message);
          // println!("hash: {:?}", hash);
          hash
        })
        .collect();
      layers.push(next_layer);
      current_layer = layers.last().unwrap();
    }

    self.leaves = leaves;
    self.layers = layers;
  }

  fn gen_proof(&self, index: usize) -> Proof<E, H> {
    let mut tmp = index;
    let mut proof_nodes = vec![];
    for row in self.layers.iter().take(self.layers.len() - 1) {
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
  let mut merkle_tree = ParallelMerkleTree::<'_, Vec<u8>, BlakeDigest>::new(&worker);
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.root().unwrap();
  assert_eq!(
    merkle_root,
    BlakeDigest(
      hex::decode("9f04496db6a8c505e88a7db289161a540a0cb953ef81c9b86103f0d6d12e8e15").unwrap()
    )
  );
  let proof = merkle_tree.gen_proof(index);
  assert_eq!(proof.leaf, hex::decode("00000003").unwrap());
  assert_eq!(
    proof.nodes,
    [
      BlakeDigest(
        hex::decode("4cd90cc0d54239ee5b3fd9989b4ef4cbebbbdd08410758cbd2d291fa364c82d5").unwrap()
      ),
      BlakeDigest(
        hex::decode("2e3d3579213e0a992d60b503f1d8fe331b8bd548e227e8dbd741ca1752077b84").unwrap()
      ),
      BlakeDigest(
        hex::decode("9a8c87bb98f1b2e0f7036a27a343dc8fd649bedc737093c2080a34c6b9f6f375").unwrap()
      ),
      BlakeDigest(
        hex::decode("ef459d75e20ce2f3fc4378ff20fe2d594fbcf16cccd986c2e0d3df41bd3bbe44").unwrap()
      )
    ]
  );

  assert_eq!(
    proof.validate(&merkle_root, index),
    Ok(hex::decode("00000003").unwrap())
  );
}

#[test]
fn test_parallel_multi_proof() {
  let indices = [2, 7, 13];
  let leaves: Vec<Vec<u8>> = vec![hex::decode("7fffffff").unwrap(); 1 << 12];
  println!("generate Merkle tree");

  let worker = Worker::new();
  let mut merkle_tree: ParallelMerkleTree<Vec<u8>, BlakeDigest> = ParallelMerkleTree::new(&worker);
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.root().unwrap();
  assert_eq!(
    merkle_root,
    BlakeDigest(
      hex::decode("a0d91c3115f9e4d9f142e7cb2f413c10f0f2f9f65d9f918b80f852f9ebc06ebc").unwrap()
    )
  );

  println!("generate Merkle proof");
  let proofs = mk_multi_branch(&merkle_tree, &indices);
  assert_eq!(proofs[0].leaf, hex::decode("7fffffff").unwrap());
  assert_eq!(
    proofs[0].nodes[0],
    BlakeDigest(
      hex::decode("b72b5371ceffa4e01aa1849cdb8705406e14791db359f826bc01a392ed26b6b9").unwrap()
    )
  );

  verify_multi_branch(&merkle_root, &indices, proofs).unwrap();
  // assert_eq!(Ok(vec![0, 0, 0, 3]));
}
