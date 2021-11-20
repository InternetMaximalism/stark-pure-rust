use ff_utils::{ff::PrimeField, ff_utils::ScalarOps};

// p[0] + p[1] * x + p[2] * x^2 + p[3] * x^3
pub fn eval_quartic<T: PrimeField>(p: [T; 4], x: T) -> T {
  let mut res = p[3]; // p[3]
  res.mul_assign(&x); // p[3] * x
  res.add_assign(&p[2]); // p[2] + p[3] * x
  res.mul_assign(&x); // (p[2] + p[3] * x) * x
  res.add_assign(&p[1]); // p[1] + (p[2] + p[3] * x) * x
  res.mul_assign(&x); // (p[1] + (p[2] + p[3] * x) * x) * x
  res.add_assign(&p[0]); // p[0] + (p[1] + (p[2] + p[3] * x) * x) * x
  res
}

#[test]
fn test_eval_quartic() {
  use bellman::Field;
  use ff_utils::f7::F7;
  use ff_utils::ff_utils::FromBytes;

  let mut poly = [F7::one(); 4];
  for (i, v) in [1u32, 2, 3, 4].iter().enumerate() {
    poly[i] = F7::from_bytes_be(&v.to_be_bytes()).unwrap();
  }
  let x = F7::from_bytes_be(&[5]).unwrap();
  let res = eval_quartic(poly, x);
  let answer = F7::from_bytes_be(&[5]).unwrap();
  assert_eq!(res, answer);
}

pub fn multi_inv<T: PrimeField + ScalarOps>(values: &[T]) -> Vec<T> {
  let mut partials = vec![T::one()];
  debug_assert!(partials.len() > 0);
  for i in 0..values.len() {
    partials.push(
      *partials.last().unwrap()
        * (if values[i] != T::zero() {
          values[i]
        } else {
          T::one()
        }),
    );
  }

  let mut inv = partials[partials.len() - 1].inverse().unwrap();

  let mut outputs = vec![T::zero(); values.len()];
  for i in (0..values.len()).rev() {
    outputs[i] = if values[i] != T::zero() {
      partials[i] * inv
    } else {
      T::zero()
    };

    let v = if values[i] != T::zero() {
      values[i]
    } else {
      T::one()
    };
    inv.mul_assign(&v);
  }

  outputs
}

/// Perform a Lagrange interpolation for a set of points
/// It's O(n^2) operations, so use with caution
pub fn interpolate<T: PrimeField>(points: &[(T, T)]) -> Option<Vec<T>> {
  let max_degree_plus_one = points.len();
  assert!(
    max_degree_plus_one >= 2,
    "should interpolate for degree >= 1"
  );
  let mut coeffs = vec![T::zero(); max_degree_plus_one];
  // external iterator
  for (k, p_k) in points.iter().enumerate() {
    let (x_k, y_k) = p_k;
    // coeffs from 0 to max_degree - 1
    let mut contribution = vec![T::zero(); max_degree_plus_one];
    let mut demoninator = T::one();
    let mut max_contribution_degree = 0;
    // internal iterator
    for (j, p_j) in points.iter().enumerate() {
      let (x_j, _) = p_j;
      if j == k {
        continue;
      }

      let mut diff = *x_k;
      diff.sub_assign(&x_j);
      demoninator.mul_assign(&diff);

      if max_contribution_degree == 0 {
        max_contribution_degree = 1;
        contribution
          .get_mut(0)
          .expect("must have enough coefficients")
          .sub_assign(&x_j);
        contribution
          .get_mut(1)
          .expect("must have enough coefficients")
          .add_assign(&T::one());
      } else {
        let mul_by_minus_x_j: Vec<T> = contribution
          .iter()
          .map(|el| {
            let mut tmp = *el;
            tmp.mul_assign(&x_j);
            tmp.negate();

            tmp
          })
          .collect();

        contribution.insert(0, T::zero());
        contribution.truncate(max_degree_plus_one);

        assert_eq!(mul_by_minus_x_j.len(), max_degree_plus_one);
        for (i, c) in contribution.iter_mut().enumerate() {
          let other = mul_by_minus_x_j
            .get(i)
            .expect("should have enough elements");
          c.add_assign(&other);
        }
      }
    }

    demoninator = demoninator.inverse().expect("denominator must be non-zero");
    for (i, this_contribution) in contribution.into_iter().enumerate() {
      let c = coeffs.get_mut(i).expect("should have enough coefficients");
      let mut tmp = this_contribution;
      tmp.mul_assign(&demoninator);
      tmp.mul_assign(&y_k);
      c.add_assign(&tmp);
    }
  }

  Some(coeffs)
}

pub fn evaluate_at_x<T: PrimeField>(coeffs: &[T], x: &T) -> T {
  let mut res = T::zero();
  let mut pow = T::one();
  for c in coeffs.iter() {
    let mut tmp = c.clone();
    tmp.mul_assign(&pow);
    res.add_assign(&tmp);

    pow.mul_assign(&x);
  }

  res
}

// #[test]
// fn test_interpolation_1() {
//   use bellman::pairing::bn256::{Bn256, Fr};
//   let points = vec![
//     (Fr::zero(), Fr::one()),
//     (Fr::one(), Fr::from_str("2").unwrap()),
//   ];
//   let interpolation_res =
//     interpolate::<Bn256>(&points[..]).expect("must interpolate a linear func");
//   assert_eq!(interpolation_res.len(), 2);
//   for (i, c) in interpolation_res.iter().enumerate() {
//     println!("Coeff {} = {}", i, c);
//   }

//   for (_i, p) in points.iter().enumerate() {
//     let (x, y) = p;
//     let val = evaluate_at_x::<Bn256>(&interpolation_res[..], &x);
//     assert_eq!(*y, val);
//     println!("Eval at {} = {}, original value = {}", x, val, y);
//   }
// }

// #[test]
// fn test_interpolation_powers_of_2() {
//   use bellman::pairing::bn256::{Bn256, Fr};
//   const MAX_POWER: u32 = Fr::CAPACITY;

//   let mut points: Vec<(Fr, Fr)> = vec![];
//   let mut power = Fr::one();
//   let two = Fr::from_str("2").unwrap();
//   for i in 0..MAX_POWER {
//     let x = Fr::from_str(&i.to_string()).unwrap();
//     let y = power.clone();
//     points.push((x, y));

//     power.mul_assign(&two);
//   }
//   let interpolation_res = interpolate::<Bn256>(&points[..]).expect("must interpolate");
//   assert_eq!(*interpolation_res.get(0).unwrap(), Fr::one());
//   assert_eq!(
//     interpolation_res.len(),
//     points.len(),
//     "array sized must match"
//   );
//   assert_eq!(
//     interpolation_res.len(),
//     MAX_POWER as usize,
//     "array size must be equal to the max power"
//   );

//   for (_i, p) in points.iter().enumerate() {
//     let (x, y) = p;
//     let val = evaluate_at_x::<Bn256>(&interpolation_res[..], &x);
//     // println!("Eval at {} = {}, original value = {}", x, val, y);
//     // assert!(*y == val, format!("must assert equality for x = {}", x) );
//     assert_eq!(*y, val);
//   }
// }

// Optimized version of the above restricted to deg-4 polynomials
pub fn multi_interp_4<T: PrimeField + ScalarOps>(
  xsets: &[[T; 4]],
  ysets: &[[T; 4]],
) -> Vec<[T; 4]> {
  let mut data = vec![];
  let mut inv_targets = vec![];
  for key in 0..xsets.len() {
    let xs = xsets[key];
    let ys = ysets[key];
    let x01 = xs[0] * xs[1];
    let x02 = xs[0] * xs[2];
    let x03 = xs[0] * xs[3];
    let x12 = xs[1] * xs[2];
    let x13 = xs[1] * xs[3];
    let x23 = xs[2] * xs[3];
    let eq0 = [
      T::zero() - x12 * xs[3],
      x12 + x13 + x23,
      T::zero() - xs[1] - xs[2] - xs[3],
      T::one(),
    ];
    let eq1 = [
      T::zero() - x02 * xs[3],
      x02 + x03 + x23,
      T::zero() - xs[0] - xs[2] - xs[3],
      T::one(),
    ];
    let eq2 = [
      T::zero() - x01 * xs[3],
      x01 + x03 + x13,
      T::zero() - xs[0] - xs[1] - xs[3],
      T::one(),
    ];
    let eq3 = [
      T::zero() - x01 * xs[2],
      x01 + x02 + x12,
      T::zero() - xs[0] - xs[1] - xs[2],
      T::one(),
    ];
    let e0 = eval_quartic(eq0, xs[0]);
    let e1 = eval_quartic(eq1, xs[1]);
    let e2 = eval_quartic(eq2, xs[2]);
    let e3 = eval_quartic(eq3, xs[3]);
    data.push([ys, eq0, eq1, eq2, eq3] as [[T; 4]; 5]);
    inv_targets.extend([e0, e1, e2, e3]);
  }

  let inv_alls = multi_inv(&inv_targets);
  let mut outputs = vec![];
  for i in 0..data.len() {
    let [ys, eq0, eq1, eq2, eq3] = data[i];
    let inv_y0 = ys[0] * inv_alls[i * 4 + 0];
    let inv_y1 = ys[1] * inv_alls[i * 4 + 1];
    let inv_y2 = ys[2] * inv_alls[i * 4 + 2];
    let inv_y3 = ys[3] * inv_alls[i * 4 + 3];
    let mut o = [T::zero(); 4];
    for i in 0..4 {
      o[i] = eq0[i] * inv_y0 + eq1[i] * inv_y1 + eq2[i] * inv_y2 + eq3[i] * inv_y3;
    }

    outputs.push(o);
  }

  // assert!(o, [self.lagrange_interp_4(xs, ys) for xs, ys in zip(xsets, ysets)]);
  outputs
}
