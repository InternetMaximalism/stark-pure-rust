use crate::hash::Digest;

use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::fmt::Error;

pub trait Element: AsRef<[u8]> + Clone + Default + Sync + Send + Eq + Debug {}

impl Element for [u8; 32] {}

impl Element for Vec<u8> {}

// `nodes` is arranged in the direction from the leaves to the roots of the tree structure.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Proof<E: Element, H: Digest> {
    pub leaf: E,
    pub nodes: Vec<H>,
}

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

pub fn verify_multi_branch<E: Element, H: Digest>(
    root: &H,
    indices: &[usize],
    proofs: Vec<Proof<E, H>>,
) -> Result<Vec<E>, Error> {
    // println!("indices: {:?}", indices);

    Ok(indices
        .iter()
        .zip(proofs)
        .map(|(index, proof)| proof.validate(root, *index).expect(""))
        .collect())
}

pub trait MerkleTree<E: Element, H: Digest> {
    // This method returns the number of leaves.
    fn width(&self) -> usize;

    // This method returns the Merkle root.
    // The Merkle root is the first node of last layer.
    fn get_root(&self) -> Option<H>;

    // This method makes the Merkle tree.
    fn update<I: IntoIterator<Item = E>>(&mut self, into: I);

    // This method returns the nodes that need to verify Merkle proof without leaf.
    fn gen_proofs(&mut self, indices: &[usize]) -> Vec<Proof<E, H>>;
}
