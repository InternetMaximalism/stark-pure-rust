use crate::interpolation::{eval_quartic, interpolate, multi_interp_4};
use crate::utils::{expand_root_of_unity, get_pseudorandom_indices};

use bellman::plonk::polynomials::Polynomial;
use bellman::worker::Worker;
use commitment::hash::Digest;
use commitment::merkle_proof_in_place::MerkleProofInPlace;
use commitment::merkle_tree::{verify_multi_branch, MerkleTree, Proof};
use ff_utils::ff::PrimeField;
use ff_utils::ff_utils::{FromBytes, ScalarOps, ToBytes};
use serde::{Deserialize, Serialize};
use std::convert::TryInto;

const MIN_DEG_DIRECT_CHECKING: usize = 16;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum FriProof<H: Digest> {
  Last {
    last: Vec<Vec<u8>>,
  },
  Middle {
    root2: H,
    column_branches: Vec<Proof<Vec<u8>, H>>,
    poly_branches: Vec<Proof<Vec<u8>, H>>,
  },
}

pub fn prove_low_degree_directly<T: PrimeField + FromBytes + ToBytes + ScalarOps, H: Digest>(
  values: &[T],
  root_of_unity: T,
  max_deg_plus_1: usize,
  exclude_multiples_of: u32,
) -> Vec<FriProof<H>> {
  let worker = Worker::new();
  prove_low_degree_rec(
    vec![],
    worker,
    values,
    root_of_unity,
    max_deg_plus_1,
    max_deg_plus_1,
    exclude_multiples_of,
  )
}

pub fn prove_low_degree<T: PrimeField + FromBytes + ToBytes + ScalarOps, H: Digest>(
  values: &[T],
  root_of_unity: T,
  max_deg_plus_1: usize,
  exclude_multiples_of: u32,
) -> Vec<FriProof<H>> {
  let worker = Worker::new();
  prove_low_degree_rec(
    vec![],
    worker,
    values,
    root_of_unity,
    max_deg_plus_1,
    MIN_DEG_DIRECT_CHECKING,
    exclude_multiples_of,
  )
}

fn prove_low_degree_rec<T: PrimeField + FromBytes + ToBytes + ScalarOps, H: Digest>(
  mut acc: Vec<FriProof<H>>,
  worker: Worker,
  values: &[T],
  root_of_unity: T,
  max_deg_plus_1: usize,
  deg_direct_checking: usize,
  exclude_multiples_of: u32,
) -> Vec<FriProof<H>> {
  println!(
    "Proving {} values are degree <= {}",
    values.len(),
    max_deg_plus_1
  );
  assert!(
    deg_direct_checking >= MIN_DEG_DIRECT_CHECKING,
    "the degree of direct checking is too low"
  );

  // Calculate the set of x coordinates
  // let xs = expand_root_of_unity(root_of_unity);
  let xs = {
    let coeffs = vec![T::one(); values.len()];
    let mut polyn = Polynomial::from_coeffs(coeffs).unwrap();
    polyn.distribute_powers(&worker, root_of_unity);
    polyn.into_coeffs()
  };

  // If the degree we are checking for is less than or equal to 32,
  // use the polynomial directly as a proof
  if max_deg_plus_1 <= deg_direct_checking {
    println!("Produced FRI proof");
    let mut pts: Vec<usize> = if exclude_multiples_of != 0 {
      (0..values.len())
        .filter(|&x| x % (exclude_multiples_of as usize) != 0)
        .collect()
    } else {
      (0..values.len()).collect()
    };

    let rest = pts.split_off(max_deg_plus_1); // pts[max_deg_plus_1..]

    // let x_vals: Vec<T> = pts.iter().map(|&pos| xs[pos]).collect();
    // let y_vals: Vec<T> = pts.iter().map(|&pos| values[pos]).collect();
    // let poly = Polynomial::from_coeffs(interpolate(&xy_pairs).unwrap()).unwrap();
    let xy_pairs = pts
      .iter()
      .map(|&pos| (xs[pos], values[pos]))
      .collect::<Vec<_>>();
    let poly = Polynomial::from_coeffs(interpolate(&xy_pairs).unwrap()).unwrap();

    for (_, &pos) in rest.iter().enumerate() {
      // Fail to prove low degree if this error occurs.
      debug_assert_eq!(poly.evaluate_at(&worker, xs[pos]), values[pos]);
    }

    acc.push(FriProof::Last {
      last: values.iter().map(|x| x.to_bytes_le().unwrap()).collect(),
    });
    return acc;
  }

  // Calculate the set of x coordinates
  // let xs = expand_root_of_unity(root_of_unity);
  // assert_eq!(values.len(), xs.len());

  // Put the values into a Merkle tree. This is the root that the
  // proof will be checked against
  let encoded_values = values
    .iter()
    .map(|x| x.to_bytes_le().unwrap())
    .collect::<Vec<_>>();
  // let m = merklize(&encoded_values);
  // let mut m: PermutedParallelMerkleTree<Vec<u8>, H> =
  //   PermutedParallelMerkleTree::new(&worker);
  // m.update(encoded_values.collect::<Vec<Vec<u8>>>());
  let mut m_tree: MerkleProofInPlace<Vec<u8>, H> = MerkleProofInPlace::new();
  m_tree.update(encoded_values);
  m_tree.gen_proofs(&[]);
  let m_root = m_tree.get_root().unwrap();

  // Select a pseudo-random x coordinate
  // let m_root = m.root().unwrap();
  let special_x = T::from_bytes_le(m_root.as_ref()).unwrap();

  // Calculate the "column" at that x coordinate
  // (see https://vitalik.ca/general/2017/11/22/starks_part_2.html)
  // We calculate the column by Lagrange-interpolating each row, and not
  // directly from the polynomial, as this is more efficient
  let quarter_len = xs.len() / 4;
  let xsets: Vec<[T; 4]> = (0..quarter_len)
    .map(|i| {
      let mut o = [T::zero(); 4];
      for j in 0..4 {
        o[j] = xs[i + quarter_len * j];
      }
      o
    })
    .collect();
  let ysets: Vec<[T; 4]> = (0..quarter_len)
    .map(|i| {
      let mut o = [T::zero(); 4];
      for j in 0..4 {
        o[j] = values[i + quarter_len * j];
      }
      o
    })
    .collect();
  let x_polys = multi_interp_4(&xsets, &ysets);
  let column: Vec<T> = x_polys
    .iter()
    .map(|p| eval_quartic(*p, special_x))
    .collect();
  let encoded_column = column
    .iter()
    .map(|p| p.to_bytes_le().unwrap())
    .collect::<Vec<_>>();
  let mut m2_tree: MerkleProofInPlace<Vec<u8>, H> = MerkleProofInPlace::new();
  m2_tree.update(encoded_column);
  m2_tree.gen_proofs(&[]);
  let m2_root = m2_tree.get_root().unwrap();
  // let mut m2_tree: PermutedParallelMerkleTree<Vec<u8>, H> =
  //   PermutedParallelMerkleTree::new(&worker);
  // m2_tree.update(encoded_column.collect::<Vec<Vec<u8>>>());
  // let m2_root = m2_tree.root().unwrap();
  // let m2 = merklize(&encoded_column);
  // let m2_root = get_root(&m2).clone();

  // Pseudo-randomly select y indices to sample
  let ys = get_pseudorandom_indices(
    m2_root.as_ref(),
    column.len() as u32,
    40,
    exclude_multiples_of,
  )
  .iter()
  .map(|&i| i as usize)
  .collect::<Vec<_>>();
  let column_branches = m2_tree.gen_proofs(&ys);

  // Compute the positions for the values in the polynomial
  let poly_positions: Vec<usize> = ys
    .iter()
    .map(|y| {
      (0..4)
        .map(|j| y + (xs.len() / 4) * j)
        .collect::<Vec<usize>>()
    })
    // .fold(Vec::with_capacity(ys.len()), |acc, val| {
    //   acc.iter().chain(val.iter()).collect::<Vec<usize>>()
    // });
    .flatten()
    .collect();
  let poly_branches = m_tree.gen_proofs(&poly_positions);

  // This component of the proof, including Merkle branches
  acc.push(FriProof::Middle {
    root2: m2_root,
    column_branches,
    poly_branches,
  });

  // Recurse...
  prove_low_degree_rec(
    acc,
    worker,
    &column,
    root_of_unity.pow(&[4u64]),
    max_deg_plus_1 / 4,
    deg_direct_checking,
    exclude_multiples_of,
  )
}

pub fn verify_low_degree_proof<T: PrimeField + FromBytes + ToBytes + ScalarOps, H: Digest>(
  merkle_root: H,
  root_of_unity: T,
  proof: &[FriProof<H>],
  max_deg_plus_1: usize,
  exclude_multiples_of: u32,
) -> Result<bool, &str> {
  let worker = Worker::new();
  verify_low_degree_proof_rec(
    worker,
    merkle_root,
    root_of_unity,
    proof,
    max_deg_plus_1,
    exclude_multiples_of,
  )
}

pub fn verify_low_degree_proof_rec<T: PrimeField + FromBytes + ToBytes + ScalarOps, H: Digest>(
  worker: Worker,
  mut merkle_root: H,
  mut root_of_unity: T,
  proof: &[FriProof<H>],
  mut max_deg_plus_1: usize,
  exclude_multiples_of: u32,
) -> Result<bool, &str> {
  // Calculate which root of unity we're working with
  let mut test_val = root_of_unity;
  let mut rou_deg = 1u64;
  while test_val != T::one() {
    rou_deg *= 2;
    test_val.mul_assign(&test_val.clone());
  }
  debug_assert_eq!(root_of_unity.pow(&[rou_deg]), T::one());

  // Powers of the given root of unity 1, p, p**2, p**3 such that p**4 = 1
  let quartic_roots_of_unity = [
    T::one(),
    root_of_unity.pow(&[rou_deg / 4]),
    root_of_unity.pow(&[rou_deg / 2]),
    root_of_unity.pow(&[rou_deg * 3 / 4]),
  ];

  // Verify the recursive components of the proof
  for prf in proof.iter().take(proof.len() - 1) {
    let (root2, column_branches, poly_branches) = match prf.clone() {
      FriProof::Middle {
        root2,
        column_branches,
        poly_branches,
      } => (root2, column_branches, poly_branches),
      _ => {
        return Err("FRI proofs must consist of FriProof::Middle except the last element.");
      }
    };
    println!("Verifying degree <= {:?}", max_deg_plus_1);

    // Calculate the pseudo-random x coordinate
    let special_x = T::from_bytes_le(merkle_root.as_ref()).unwrap();

    // Calculate the pseudo-randomly sampled y indices
    let ys = get_pseudorandom_indices(
      root2.as_ref(),
      (rou_deg / 4).try_into().unwrap(),
      40,
      exclude_multiples_of,
    )
    .iter()
    .map(|&i| i as usize)
    .collect::<Vec<usize>>();

    // Compute the positions for the values in the polynomial
    let mut poly_positions: Vec<usize> = vec![];
    for y in ys.clone() {
      for j in 0usize..4usize {
        let rou_deg_divided_by4: usize = (rou_deg / 4).try_into().unwrap();
        poly_positions.push(j * rou_deg_divided_by4 + y);
      }
    }

    // Verify Merkle branches
    let column_values = verify_multi_branch(&root2, &ys, column_branches).unwrap();
    let poly_values = verify_multi_branch(&merkle_root, &poly_positions, poly_branches).unwrap();

    // For each y coordinate, get the x coordinates on the row, the values on
    // the row, and the value at that y from the column
    let mut x_coords: Vec<[T; 4]> = vec![];
    let mut rows: Vec<[T; 4]> = vec![];
    let mut column_vals = vec![];
    for i in 0..ys.len() {
      // The x coordinates from the polynomial
      let x1 = root_of_unity.pow(&[ys[i].try_into().unwrap()]);
      let mut x_coord = [T::zero(); 4];
      for j in 0..4 {
        let mut v = quartic_roots_of_unity[j];
        v.mul_assign(&x1);
        x_coord[j] = v;
      }
      x_coords.push(x_coord);

      // The values from the original polynomial
      let mut row = [T::zero(); 4];
      for j in 0..4 {
        row[j] = T::from_bytes_le(&poly_values[i * 4 + j]).unwrap();
      }
      rows.push(row);
      column_vals.push(T::from_bytes_le(&column_values[i]).unwrap());
    }

    // Verify for each selected y coordinate that the four points from the
    // polynomial and the one point from the column that are on that y
    // coordinate are on the same deg < 4 polynomial
    let polys = multi_interp_4(&x_coords, &rows);

    for key in 0..polys.len() {
      let p = polys[key];
      let c = column_vals[key];
      assert_eq!(eval_quartic(p, special_x), c);
    }

    // Update constants to check the next proof
    merkle_root = root2.clone();
    root_of_unity = root_of_unity.pow(&[4]);
    max_deg_plus_1 /= 4;
    rou_deg /= 4;
  }

  println!("Verifying degree <= {:?}", max_deg_plus_1);
  assert!(
    max_deg_plus_1 >= MIN_DEG_DIRECT_CHECKING / 2,
    "the degree of direct checking is too low"
  );

  let last_data = if let Some(FriProof::Last { last }) = proof.last() {
    last.iter().map(|v| v.to_vec()).collect::<Vec<_>>()
  } else {
    return Err("The last element of FRI proofs must be FriProof::Last.");
    // return Err(SimpleError::new("The last element of FRI proofs must be FriProof::Last."));
  };
  let last_data_len = last_data.len();
  assert!(last_data.len() > max_deg_plus_1);

  let decoded_last_data: Vec<T> = last_data
    .iter()
    .map(|x| T::from_bytes_le(x).unwrap())
    .collect();

  // Check the Merkle root matches up
  let mut m_tree: MerkleProofInPlace<Vec<u8>, H> = MerkleProofInPlace::new();
  m_tree.update(last_data);
  m_tree.gen_proofs(&[]);
  let m_root = m_tree.get_root().unwrap();
  // let mut m_tree: PermutedParallelMerkleTree<Vec<u8>, H> =
  //   PermutedParallelMerkleTree::new(&worker);
  // m_tree.update(last_data.clone());
  // let m_tree = merklize(last_data);
  assert_eq!(m_root, merkle_root);

  // Check the degree of the last_data
  let xs = expand_root_of_unity(root_of_unity);
  let mut pts: Vec<usize> = if exclude_multiples_of != 0 {
    (0..last_data_len)
      .filter(|&pos| pos % (exclude_multiples_of as usize) != 0)
      .collect()
  } else {
    (0..last_data_len).collect()
  };

  // Verify the direct components of the proof
  let rest = pts.split_off(max_deg_plus_1); // pts[max_deg_plus_1..]

  // let x_vals: Vec<T> = pts.iter().map(|&pos| xs[pos]).collect();
  // let y_vals: Vec<T> = pts.iter().map(|&pos| decoded_last_data[pos]).collect();
  // let poly = Polynomial::from_coeffs(lagrange_interp(&x_vals, &y_vals)).unwrap();
  let xy_pairs = pts
    .iter()
    .map(|&pos| (xs[pos], decoded_last_data[pos]))
    .collect::<Vec<_>>();
  let poly = Polynomial::from_coeffs(interpolate(&xy_pairs).unwrap()).unwrap();
  for (_, &pos) in rest.iter().enumerate() {
    assert_eq!(poly.evaluate_at(&worker, xs[pos]), decoded_last_data[pos]);
  }

  println!("FRI proof verified");
  Ok(true)
}

// pub fn fri_proof_bin_length(fri_proof: &[FriProof]) -> usize {
//   fri_proof
//     .iter()
//     .map(|x| match x {
//       FriProof::Middle {
//         root2: _,
//         column_branches,
//         poly_branches,
//       } => 32 + bin_length(column_branches) + bin_length(poly_branches),
//       FriProof::Last { last } => last.iter().map(|x| x.len()).sum(),
//     })
//     .sum()
// }
