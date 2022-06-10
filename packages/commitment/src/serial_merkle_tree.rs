use crate::{
  hash::Digest,
  merkle_tree::{Element, MerkleTree, Proof},
};

// `SerialMerkleTree` describe a simple Merkle tree structure.
// This implement is only for testing usage.
// A layer is the set of the tree nodes which are the same depth.
// The first layer is the hash of each leaves
// The last layer is the tree root.
pub struct SerialMerkleTree<E: Element, H: Digest> {
  pub leaves: Vec<E>,
  pub layers: Vec<Vec<H>>,
}

impl<E: Element, H: Digest> SerialMerkleTree<E, H> {
  pub fn new() -> Self {
    let leaves = vec![];
    let layers = vec![];
    Self { leaves, layers }
  }
}

impl<E: Element, H: Digest> SerialMerkleTree<E, H> {
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

impl<E: Element, H: Digest> MerkleTree<E, H> for SerialMerkleTree<E, H> {
  fn width(&self) -> usize {
    self.leaves.len()
  }

  fn get_root(&self) -> Option<H> {
    self.layers.last().map(|n| n[0].clone())
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
          
          // println!("hash: {:?}", hash);

          H::hash(&message)
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
