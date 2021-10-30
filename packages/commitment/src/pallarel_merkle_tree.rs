use crate::{
  hash::Digest,
  merkle_tree::{Element, MerkleTree, Proof},
  multicore::Worker,
  serial_merkle_tree::SerialMerkleTree,
};

pub struct ParallelMerkleTree<E: Element, H: Digest> {
  worker: Worker,
  pub leaves: Vec<E>,
  pub layers: Vec<Vec<H>>,
}

impl<E: Element, H: Digest> ParallelMerkleTree<E, H> {
  pub fn new() -> Self {
    let worker = Worker::new();
    let leaves = vec![];
    let layers = vec![];
    Self {
      worker,
      leaves,
      layers,
    }
  }
}

impl<E: Element, H: Digest> ParallelMerkleTree<E, H> {
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

impl<E: Element, H: Digest> MerkleTree<E, H> for ParallelMerkleTree<E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn get_root(&self) -> Option<H> {
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

  fn gen_proofs(&mut self, indices: &[usize]) -> Vec<Proof<E, H>> {
    // println!("indices: {:?}", indices);
    let mut proofs = vec![];
    for &index in indices {
      let proof = self.gen_proof(index);
      proofs.push(proof);
    }

    proofs
  }
}

#[test]
fn test_parallel_single_proof() {
  use crate::blake::BlakeDigest;

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

  let mut merkle_tree = ParallelMerkleTree::<Vec<u8>, BlakeDigest>::new();
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.get_root().unwrap();
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
  use crate::blake::BlakeDigest;
  use crate::merkle_tree::verify_multi_branch;

  let indices = [2, 7, 13];
  let leaves: Vec<Vec<u8>> = vec![hex::decode("7fffffff").unwrap(); 1 << 12];
  println!("generate Merkle tree");

  let mut merkle_tree: ParallelMerkleTree<Vec<u8>, BlakeDigest> = ParallelMerkleTree::new();
  merkle_tree.update(leaves);
  let merkle_root = merkle_tree.get_root().unwrap();
  assert_eq!(
    merkle_root,
    BlakeDigest(
      hex::decode("a0d91c3115f9e4d9f142e7cb2f413c10f0f2f9f65d9f918b80f852f9ebc06ebc").unwrap()
    )
  );

  println!("generate Merkle proof");
  let proofs = merkle_tree.gen_proofs(&indices);
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
