use crate::delayed::Delayed;
use crate::hash::Digest;
use crate::lazily;
use crate::merkle_tree::{Element, MerkleTree, Proof};
use crate::multicore::Worker;
use crate::utils::{is_a_power_of_2, log2_ceil};
use std::ops::Deref;

pub struct MerkleProofInPlace<'a, E: Element, H: Digest> {
  worker: Worker,
  pub leaves: Vec<Delayed<'a, E>>,
  root: Option<H>,
}

impl<'a, E: Element, H: Digest> MerkleProofInPlace<'a, E, H> {
  pub fn new() -> Self {
    let worker = Worker::new();
    let leaves = vec![];
    let root = Some(H::default());
    Self {
      worker,
      leaves,
      root,
    }
  }
}

impl<'a, E: Element, H: Digest> MerkleTree<E, H> for MerkleProofInPlace<'a, E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn get_root(&self) -> Option<H> {
    if let Some(root) = self.root.clone() {
      Some(root)
    } else {
      // self.gen_proofs(&[]);
      // Some(self.root.unwrap().clone())
      None
    }
  }

  fn update<I: IntoIterator<Item = E>>(&mut self, leaves: I) {
    self.leaves = leaves
      .into_iter()
      .map(|x| lazily!(x.clone()))
      .collect::<Vec<_>>();
  }

  fn gen_proofs(&mut self, indices: &[usize]) -> Vec<Proof<E, H>> {
    let (merkle_proof, merkle_root) =
      gen_multi_proofs_multi_core(&self.leaves, indices, &self.worker);
    self.root = Some(merkle_root);
    merkle_proof
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

#[test]
fn test_gen_multi_proofs_multi_core() {
  use crate::blake::BlakeDigest;
  use crate::lazily;
  use crate::merkle_tree::verify_multi_branch;
  use crate::serial_merkle_tree::SerialMerkleTree;

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
  let merkle_root2 = merkle_tree.get_root().unwrap();
  assert_eq!(merkle_root, merkle_root2, "invalid Merkle root");

  let merkle_proofs2 = merkle_tree.gen_proofs(&indices);
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
