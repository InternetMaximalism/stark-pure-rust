use ff::PrimeField;
use hex::ToHex;
use num::bigint::BigUint;
use std::convert::TryInto;
// use simple_error::SimpleError;
use crate::fft::expand_root_of_unity;
use crate::permuted_tree::{merklize, mk_multi_branch, verify_multi_branch, bin_length};
use crate::poly_utils::{eval_poly_at, eval_quartic, lagrange_interp, multi_interp_4};
use crate::utils::get_pseudorandom_indices;

pub fn parse_hex_to_decimal(value: String) -> String {
  BigUint::from_bytes_be(value.as_ref()).to_str_radix(10)
}

pub enum FriProof {
  Last {
    last: Vec<String>,
  },
  Middle {
    root2: String,
    column_branches: Vec<Vec<String>>,
    poly_branches: Vec<Vec<String>>,
  },
}

pub fn prove_low_degree<T: PrimeField + ToHex>(
  values: &[T],
  root_of_unity: T,
  max_deg_plus_1: usize,
  exclude_multiples_of: i32,
) -> Vec<FriProof> {
  prove_low_degree_rec(
    vec![],
    values,
    root_of_unity,
    max_deg_plus_1,
    exclude_multiples_of,
  )
}

fn prove_low_degree_rec<T: PrimeField + ToHex>(
  mut acc: Vec<FriProof>,
  values: &[T],
  root_of_unity: T,
  max_deg_plus_1: usize,
  exclude_multiples_of: i32,
) -> Vec<FriProof> {
  println!(
    "Proving {} values are degree <= {}",
    values.len(),
    max_deg_plus_1
  );

  // If the degree we are checking for is less than or equal to 32,
  // use the polynomial directly as a proof
  if max_deg_plus_1 <= 16 {
    println!("Produced FRI proof");
    acc.push(FriProof::Last {
      last: values.iter().map(|x| x.encode_hex()).collect(),
    });
    return acc;
  }

  // Calculate the set of x coordinates
  let xs = expand_root_of_unity(root_of_unity);
  assert_eq!(values.len(), xs.len());

  // Put the values into a Merkle tree. This is the root that the
  // proof will be checked against
  let encoded_values: Vec<String> = values.iter().map(|x| x.encode_hex::<String>()).collect();
  let m = merklize(&encoded_values);

  // Select a pseudo-random x coordinate
  let special_x = T::from_str(&parse_hex_to_decimal(m.last().unwrap().clone())).unwrap();

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
  let encoded_column: Vec<String> = column.iter().map(|p| p.encode_hex::<String>()).collect();
  let m2 = merklize(&encoded_column);
  let m2_root = m2[1].clone();

  // Pseudo-randomly select y indices to sample
  let ys = get_pseudorandom_indices(
    m2_root.clone(),
    column.len().try_into().unwrap(),
    40,
    exclude_multiples_of,
  );

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

  // This component of the proof, including Merkle branches
  acc.push(FriProof::Middle {
    root2: m2_root,
    column_branches: mk_multi_branch(&m2, &ys),
    poly_branches: mk_multi_branch(&m, &poly_positions),
  });

  // Recurse...
  prove_low_degree_rec(
    acc,
    &column,
    root_of_unity.pow_vartime(&[4u64]),
    max_deg_plus_1 / 4,
    exclude_multiples_of,
  )
}

pub fn verify_low_degree_proof<T: PrimeField>(
  merkle_root: String,
  root_of_unity: T,
  proof: &[FriProof],
  max_deg_plus_1: usize,
  exclude_multiples_of: i32,
) -> Result<bool, &str> {
  verify_low_degree_proof_rec(
    merkle_root,
    root_of_unity,
    proof,
    max_deg_plus_1,
    exclude_multiples_of,
  )
}

pub fn verify_low_degree_proof_rec<T: PrimeField>(
  mut merkle_root: String,
  mut root_of_unity: T,
  proof: &[FriProof],
  mut max_deg_plus_1: usize,
  exclude_multiples_of: i32,
) -> Result<bool, &str> {
  // Calculate which root of unity we're working with
  let mut test_val = root_of_unity;
  let mut rou_deg = 1u64;
  while test_val != T::one() {
    rou_deg *= 2;
    test_val = test_val * test_val;
  }

  // Powers of the given root of unity 1, p, p**2, p**3 such that p**4 = 1
  let quartic_roots_of_unity = [
    T::one(),
    root_of_unity.pow_vartime(&[rou_deg / 4]),
    root_of_unity.pow_vartime(&[rou_deg / 2]),
    root_of_unity.pow_vartime(&[rou_deg * 3 / 4]),
  ];

  // Verify the recursive components of the proof
  for prf in proof.iter().take(proof.len() - 1) {
    let (root2, column_branches, poly_branches) = match prf {
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
    let special_x = T::from_str(&parse_hex_to_decimal(merkle_root.clone())).unwrap();

    // Calculate the pseudo-randomly sampled y indices
    let ys = get_pseudorandom_indices(
      root2.clone(),
      (rou_deg / 4).try_into().unwrap(),
      40,
      exclude_multiples_of,
    );

    // Compute the positions for the values in the polynomial
    let mut poly_positions: Vec<usize> = vec![];
    for y in ys.clone() {
      for j in 0usize..4usize {
        let rou_deg_divided_by4: usize = (rou_deg / 4).try_into().unwrap();
        poly_positions.push(j * rou_deg_divided_by4 + y);
      }
    }

    // Verify Merkle branches
    let column_values = verify_multi_branch(root2, &ys, column_branches);
    let poly_values = verify_multi_branch(&merkle_root, &poly_positions, poly_branches);

    // For each y coordinate, get the x coordinates on the row, the values on
    // the row, and the value at that y from the column
    let mut x_coords: Vec<[T; 4]> = vec![];
    let mut rows: Vec<[T; 4]> = vec![];
    let mut column_vals = vec![];
    for i in 0..ys.len() {
      // The x coordinates from the polynomial
      let x1 = root_of_unity.pow_vartime(&[ys[i].try_into().unwrap()]);
      let mut x_coord = [T::zero(); 4];
      for j in 0..4 {
        x_coord[j] = quartic_roots_of_unity[j] * x1;
      }
      x_coords.push(x_coord);

      // The values from the original polynomial
      let mut row = [T::zero(); 4];
      for j in 0..4 {
        row[j] = T::from_str(&parse_hex_to_decimal(poly_values[i * 4 + j].clone())).unwrap();
      }
      rows.push(row);
      column_vals.push(T::from_str(&parse_hex_to_decimal(column_values[i].clone())).unwrap());
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
    root_of_unity = root_of_unity.pow_vartime(&[4]);
    max_deg_plus_1 /= 4;
    rou_deg /= 4;
  }

  // Verify the direct components of the proof
  println!("Verifying degree <= {:?}", max_deg_plus_1);
  assert!(max_deg_plus_1 <= 16);

  // Check the Merkle root matches up
  let last_data = if let Some(FriProof::Last { last }) = proof.last() {
    last
  } else {
    return Err("The last element of FRI proofs must be FriProof::Last.");
    // return Err(SimpleError::new("The last element of FRI proofs must be FriProof::Last."));
  };
  let m_tree = merklize(last_data);
  assert_eq!(m_tree[1], merkle_root);

  // Check the degree of the last_data
  let powers = expand_root_of_unity(root_of_unity);
  let mut pts: Vec<usize> = if exclude_multiples_of != 0 {
    (0..last_data.len())
      .filter(|x| {
        let exclude_multiples_of_usize: usize = exclude_multiples_of.try_into().unwrap();
        *x % exclude_multiples_of_usize != 0
      })
      .collect()
  } else {
    (0..last_data.len()).collect()
  };

  let decoded_last_data: Vec<T> = last_data.iter().map(|x| T::from_str(&parse_hex_to_decimal(x.clone())).unwrap()).collect();
  let rest = pts.split_off(max_deg_plus_1.try_into().unwrap());
  let xs: Vec<T> = pts.iter().map(|x| powers[*x]).collect();
  let ys: Vec<T> = pts.iter().map(|x| decoded_last_data[*x]).collect();
  let poly = lagrange_interp(&xs, &ys);
  for x in rest {
    assert!(eval_poly_at(&poly, powers[x]) == decoded_last_data[x]);
  }

  println!("FRI proof verified");
  Ok(true)
}

pub fn fri_proof_bin_length(fri_proof: &[FriProof]) -> usize {
  fri_proof
    .iter()
    .map(|x| match x {
      FriProof::Middle {
        root2: _,
        column_branches,
        poly_branches,
      } => 32 + bin_length(column_branches) + bin_length(poly_branches),
      FriProof::Last { last } => last.iter().map(|x| x.len()).sum(),
    })
    .sum()
}
