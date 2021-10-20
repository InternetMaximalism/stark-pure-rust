use crate::delayed::Delayed;
#[warn(unused_imports)]
use crate::lazily;
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
    // let num_of_leaves = 2usize.pow(self.height() as u32);
    // let index = permute4_index(index, num_of_leaves);
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

pub fn gen_multi_proofs_in_place<E: Element, H: Digest>(
  current_nodes: &mut [H],
  indices: &[usize],
  steps: usize,
) -> Vec<Proof<E, H>> {
  // let leaves: Vec<E> = permute4_values(&leaves);
  // let indices = permute4_indices(indices, leaves.len());
  // println!("leaves: {:?}", leaves);
  // println!("indices: {:?}", indices);
  assert!(is_a_power_of_2(current_nodes.len()));

  let mut proofs = vec![
    Proof {
      leaf: E::default(),
      nodes: vec![]
    };
    indices.len()
  ];

  let mut log_nodes_interval = log2_ceil(steps);
  assert_eq!(1 << log_nodes_interval, steps);
  let current_nodes_len = current_nodes.len();
  while (1 << log_nodes_interval) < current_nodes_len {
    for (i, &tmp_index) in indices.iter().enumerate() {
      let twin_index = ((tmp_index >> log_nodes_interval) ^ 1) << log_nodes_interval;
      proofs[i].nodes.push(current_nodes[twin_index].clone());
    }

    let nodes_interval = 1 << log_nodes_interval;
    let chunks_num = nodes_interval << 1;
    for node in current_nodes.chunks_mut(chunks_num) {
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

pub fn gen_multi_proofs_multi_core<E: Element, H: Digest>(
  leaves: &[Delayed<'_, E>],
  indices: &[usize],
  worker: &Worker,
) -> (Vec<Proof<E, H>>, H) {
  // let leaves: Vec<E> = permute4_values(&leaves);
  // let indices = permute4_indices(indices, leaves.len());
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
fn test_serial_internal_roots() {
  // let leaves: Vec<Vec<u8>> = vec![
  //   hex::decode("7fffffff").unwrap(),
  //   hex::decode("80000000").unwrap(),
  //   hex::decode("00000003").unwrap(),
  //   hex::decode("00000000").unwrap(),
  //   hex::decode("7ffffffe").unwrap(),
  //   hex::decode("80000001").unwrap(),
  //   hex::decode("00000004").unwrap(),
  //   hex::decode("00000001").unwrap(),
  //   hex::decode("7ffffffd").unwrap(),
  //   hex::decode("80000002").unwrap(),
  //   hex::decode("00000005").unwrap(),
  //   hex::decode("00000002").unwrap(),
  //   hex::decode("7ffffffc").unwrap(),
  //   hex::decode("80000003").unwrap(),
  //   hex::decode("00000006").unwrap(),
  //   hex::decode("00000003").unwrap(),
  // ];

  let leaves: Vec<Vec<u8>> = vec![
    hex::decode("00000000").unwrap(),
    hex::decode("00000001").unwrap(),
    hex::decode("00000002").unwrap(),
    hex::decode("00000003").unwrap(),
    hex::decode("00000004").unwrap(),
    hex::decode("00000005").unwrap(),
    hex::decode("00000006").unwrap(),
    hex::decode("00000007").unwrap(),
    hex::decode("00000008").unwrap(),
    hex::decode("00000009").unwrap(),
    hex::decode("0000000a").unwrap(),
    hex::decode("0000000b").unwrap(),
    hex::decode("0000000c").unwrap(),
    hex::decode("0000000d").unwrap(),
    hex::decode("0000000e").unwrap(),
    hex::decode("0000000f").unwrap(),
  ];

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
  let mut merkle_proofs2 = vec![];
  merkle_proofs2.push(merkle_tree.gen_proof(indices[0]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[1]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[2]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[3]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[4]));
  merkle_proofs2.push(merkle_tree.gen_proof(indices[5]));
  // println!("merkle_proofs2: {:?}", merkle_proofs2);
  // println!("merkle_proofs1: {:?}", merkle_proofs);
  // println!(
  //   "{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}",
  //   merkle_proofs2[0].nodes,
  //   merkle_proofs2[1].nodes,
  //   merkle_proofs2[2].nodes,
  //   merkle_proofs2[3].nodes,
  //   merkle_proofs2[4].nodes,
  //   merkle_proofs2[5].nodes,
  //   merkle_proofs2[6].nodes,
  //   merkle_proofs2[7].nodes,
  // );
  // println!(
  //   "{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}",
  //   merkle_proofs[0].nodes,
  //   merkle_proofs[1].nodes,
  //   merkle_proofs[2].nodes,
  //   merkle_proofs[3].nodes,
  //   merkle_proofs[4].nodes,
  //   merkle_proofs[5].nodes,
  //   merkle_proofs[6].nodes,
  //   merkle_proofs[7].nodes
  // );

  // for i in 0..indices.len() {
  //   assert_eq!(merkle_proofs[i].leaf, merkle_proofs2[i].leaf);
  //   assert_eq!(merkle_proofs[i].nodes, merkle_proofs2[i].nodes);
  // }
  assert_eq!(merkle_root, merkle_root2);

  verify_multi_branch(&merkle_root, &indices, merkle_proofs).unwrap();
}

impl<E: Element, H: Digest> SerialMerkleTree<E, H> {
  pub fn new() -> Self {
    let leaves = vec![];
    let nodes = vec![];
    SerialMerkleTree { leaves, nodes }
  }

  pub fn gen_internal_roots<I: IntoIterator<Item = E>>(&mut self, leaves: I, height: u32) -> H {
    let mut leaves = leaves.into_iter();
    // let leaves: Vec<E> = leaves.into_iter().collect::<Vec<E>>();
    // println!("leaves: {:?}", leaves);

    // assert!(is_a_power_of_2(leaves.len()));
    let chunk_size = 2usize.pow(height);

    let mut internal_roots = vec![H::default(); chunk_size];
    for (chunk_index, internal_root) in internal_roots.chunks_mut(1).enumerate() {
      let beginning_index = chunk_index * chunk_size;
      let mut sub_nodes: Vec<H> = vec![];
      println!("{}", chunk_index);

      for i in beginning_index..(beginning_index + chunk_size) {
        let digest = H::hash(leaves.nth(i).get_or_insert(E::default()).as_ref());
        sub_nodes.push(digest);
      }

      while sub_nodes.len() >= 2 {
        sub_nodes = sub_nodes
          .chunks(2)
          .map(|node| {
            let mut message: Vec<u8> = vec![];
            message.extend(node[0].as_ref());
            if node.len() > 1 {
              message.extend(node[1].as_ref());
            }
            H::hash(&message)
          })
          .collect::<Vec<_>>();
      }

      internal_root[0] = sub_nodes[0].clone();
    }

    while internal_roots.len() >= 2 {
      internal_roots = internal_roots
        .chunks(2)
        .map(|node| {
          let mut message: Vec<u8> = vec![];
          message.extend(node[0].as_ref());
          if node.len() > 1 {
            message.extend(node[1].as_ref());
          }
          H::hash(&message)
        })
        .collect::<Vec<_>>();
    }

    internal_roots[0].clone()
  }
}

pub trait MerkleTree<E: Element, H: Digest> {
  fn width(&self) -> usize;
  fn root(&self) -> Option<H>;
  fn update<I: IntoIterator<Item = E>>(&mut self, into: I);
  fn gen_proof(&self, index: usize) -> Proof<E, H>;
}

pub struct SerialMerkleTree<E: Element, H: Digest> {
  leaves: Vec<E>,
  nodes: Vec<Vec<H>>,
}

impl<E: Element, H: Digest> MerkleTree<E, H> for SerialMerkleTree<E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn root(&self) -> Option<H> {
    if let Some(n) = self.nodes.last() {
      Some(n[0].clone())
    } else {
      None
    }
  }

  fn update<I: IntoIterator<Item = E>>(&mut self, leaves: I) {
    let leaves: Vec<E> = leaves.into_iter().collect::<Vec<E>>();
    // println!("leaves: {:?}", leaves);

    // assert!(is_a_power_of_2(leaves.len()));

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
          if node.len() > 1 {
            message.extend(node[1].as_ref());
          }
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
