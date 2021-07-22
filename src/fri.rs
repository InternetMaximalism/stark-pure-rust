use ff::PrimeField;
use std::convert::TryInto;
// use crate::ff_utils::Fp;
use crate::fft::expand_root_of_unity;
use crate::permuted_tree::{merklize, mk_multi_branch, verify_multi_branch};
use crate::poly_utils::{eval_poly_at, eval_quartic, lagrange_interp, multi_interp_4};
use crate::utils::get_pseudorandom_indices;

pub struct FriProof {
  root2: String,
  column_branches: Vec<Vec<String>>,
  poly_branches: Vec<Vec<String>>,
}

pub fn prove_low_degree<T: PrimeField>(
  values: &[T],
  root_of_unity: T,
  max_deg_plus_1: i32,
  exclude_multiples_of: i32,
) -> Vec<FriProof> {
  println!(
    "Proving %{} values are degree <= %{}",
    values.len(),
    max_deg_plus_1
  );

  // If the degree we are checking for is less than or equal to 32,
  // use the polynomial directly as a proof
  if max_deg_plus_1 <= 16 {
    println!("Produced FRI proof");
    return vec![values
      .iter()
      .map(|x| format!("{:?}", x.to_repr()))
      .collect()];
  }

  // Calculate the set of x coordinates
  let xs = expand_root_of_unity(root_of_unity);
  assert_eq!(values.len(), xs.len());

  // Put the values into a Merkle tree. This is the root that the
  // proof will be checked against
  let m = merklize(values);

  // Select a pseudo-random x coordinate
  let special_x = T::from_str(&m[m.len() - 1]).unwrap();

  // Calculate the "column" at that x coordinate
  // (see https://vitalik.ca/general/2017/11/22/starks_part_2.html)
  // We calculate the column by Lagrange-interpolating each row, and not
  // directly from the polynomial, as this is more efficient
  let quarter_len = xs.len() / 4;
  let xsets: Vec<[T; 4]> = (0..quarter_len)
    .map(|i| {
      let o = [T::zero(); 4];
      for j in 0..4 {
        o[j] = xs[i + quarter_len * j];
      }
      o
    })
    .collect();
  let ysets: Vec<[T; 4]> = (0..quarter_len)
    .map(|i| {
      let o = [T::zero(); 4];
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
  let m2 = merklize(&column);

  // Pseudo-randomly select y indices to sample
  let ys = get_pseudorandom_indices(
    m2[1],
    column.len().try_into().unwrap(),
    40,
    exclude_multiples_of,
  );

  // Compute the positions for the values in the polynomial
  let poly_positions: Vec<i32> = ys
    .iter()
    .map(|y| (0..4).map(|j| y + (xs.len() / 4).try_into().unwrap() * j))
    .fold(Vec::with_capacity(ys.len()), |acc, val| {
      acc.extend(val);
      acc
    }); // .flat()

  // This component of the proof, including Merkle branches
  let mut o = vec![FriProof {
    root2: m2[1],
    column_branches: mk_multi_branch(m2, ys),
    poly_branches: mk_multi_branch(m, poly_positions),
  }];

  // Recurse...
  o.push(prove_low_degree(
    &column,
    root_of_unity.pow_vartime(4u64),
    max_deg_plus_1 / 4,
    modulus,
    exclude_multiples_of,
  ));
  o
}

fn verify_low_degree_proof<T: PrimeField>(
  merkle_root: String,
  root_of_unity: T,
  proof: &[FriProof<T>],
  maxdeg_plus_1: i32,
  exclude_multiples_of: i32,
) -> bool {
  // Calculate which root of unity we're working with
  let mut testval = root_of_unity;
  let mut roudeg = 1u64;
  while testval != T::one() {
    roudeg *= 2;
    testval = testval * testval;
  }

  // Powers of the given root of unity 1, p, p**2, p**3 such that p**4 = 1
  let quartic_roots_of_unity = [
    T::one(),
    root_of_unity.pow_vartime(roudeg / 4 as u64),
    root_of_unity.pow_vartime(roudeg / 2 as u64),
    root_of_unity.pow_vartime(roudeg * 3 / 4 as u64),
  ];

  // Verify the recursive components of the proof
  for prf in proof[0..(proof.len() - 1)] {
    let root2 = prf.root2;
    let column_branches = prf.column_branches;
    let poly_branches = prf.poly_branches;
    println!("Verifying degree <= {:?}", maxdeg_plus_1);

    // Calculate the pseudo-random x coordinate
    let special_x = T::from_str(merkle_root);

    // Calculate the pseudo-randomly sampled y indices
    let ys = get_pseudorandom_indices(
      root2,
      (roudeg / 4).try_into().unwrap(),
      40,
      exclude_multiples_of = exclude_multiples_of,
    );

    // Compute the positions for the values in the polynomial
    let poly_positions = ys
      .iter()
      .fold(vec![], |y| [y + 0..4.map(|j| (roudeg / 4) * j)]);

    // Verify Merkle branches
    let column_values = verify_multi_branch(root2, ys, column_branches);
    let poly_values = verify_multi_branch(merkle_root, poly_positions, poly_branches);

    // For each y coordinate, get the x coordinates on the row, the values on
    // the row, and the value at that y from the column
    let mut xcoords = vec![];
    let mut rows = vec![];
    let mut columnvals = vec![];
    for i in 0..ys.len() {
      let y = ys[i];
      // The x coordinates from the polynomial
      let x1 = root_of_unity.pow_vartime(y);
      xcoords.push((0..4).map(|j| quartic_roots_of_unity[j] * x1));

      // The values from the original polynomial
      let row = poly_values
        .drain(i * 4..i * 4 + 4)
        .iter()
        .map(|x| T::from_str(x));
      rows.push(row);

      columnvals.push(T::from_str(column_values[i]));
    }

    // Verify for each selected y coordinate that the four points from the
    // polynomial and the one point from the column that are on that y
    // coordinate are on the same deg < 4 polynomial
    let polys = multi_interp_4(xcoords, rows);

    for key in 0..polys.len() {
      let p = polys[key];
      let c = columnvals[key];
      assert_eq!(eval_quartic(p, special_x), c);
    }

    // Update constants to check the next proof
    merkle_root = root2;
    root_of_unity = root_of_unity.pow_vartime(4);
    maxdeg_plus_1 /= 4;
    roudeg /= 4;
  }

  // Verify the direct components of the proof
  let data = [
    T::from_str(proof[proof.len() - 1].root2),
    T::from_str(proof[proof.len() - 1].column_branches),
    T::from_str(proof[proof.len() - 1].poly_branches),
  ];
  println!("Verifying degree <= {:?}", maxdeg_plus_1);
  assert!(maxdeg_plus_1 <= 16);

  // Check the Merkle root matches up
  let mtree = merklize(data);
  assert_eq!(mtree[1], merkle_root);

  // Check the degree of the data
  let powers = expand_root_of_unity(root_of_unity);
  let mut pts: Vec<usize> = if exclude_multiples_of != 0 {
    (0..data.len())
      .filter(|x| *x % exclude_multiples_of.try_into().unwrap() != 0)
      .collect()
  } else {
    (0..data.len()).collect()
  };

  let rest = pts.split_off(maxdeg_plus_1.try_into().unwrap());
  let poly = lagrange_interp(pts.iter().map(|x| powers[*x]), pts.iter().map(|x| data[*x]));
  for x in rest {
    assert!(eval_poly_at(poly, powers[x]) == data[x]);
  }

  println!("FRI proof verified");
  true
}
