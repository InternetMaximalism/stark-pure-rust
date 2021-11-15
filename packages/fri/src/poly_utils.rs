// use bellman::plonk::polynomials::{Coefficients, Polynomial};
// // use core::ops::{
// //   Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
// // };
// // use std::io::Error;
// use ff_utils::ff::PrimeField;
// use ff_utils::ff_utils::ScalarOps;
// use std::cmp::max;
// use std::collections::HashMap;
// // use num_traits::{Num, NumAssign};

// // pub struct Polyn = Vec;

// // impl<T: Clone + Num> Add<Polyn<T>> for Polyn<T> {
// //   fn add(self, other: Self) -> Self {
// //     (0..max(self.len(), other.len()))
// //       .map(|i| if i < self.len() { self[i] } else { T::zero() } + if i < other.len() { other[i] } else { T::zero() })
// //       .collect()
// //   }
// // }

// // impl<T: Clone + Num> Add<T> for Polyn<T> {
// //   fn add(self, other: T) -> Self {
// //     (0..self.len())
// //       .map(|i| if i < self.len() { self[i] } else { T::zero() } + if i < 1 { other } else { T::zero() })
// //       .collect()
// //   }
// // }

// // impl<'a, T: Clone + Num> Add<&'a T> for Polyn<T> {}
// // impl<'a, 'b, T: Clone + Num> Add<&'a T> for &'b Polyn<T> {}
// // impl<'a, 'b, T: Clone + Num> Add<&'b Polyn<T>> for &'a Polyn<T> {}
// // impl<'a, T: Clone + Num> Add<Polyn<T>> for &'a Polyn<T> {}
// // impl<'a, T: Clone + Num> Add<T> for &'a Polyn<T> {}
// // impl<'a, T: Clone + NumAssign> AddAssign<&'a Polyn<T>> for Polyn<T> {}
// // impl<'a, T: Clone + NumAssign> AddAssign<&'a T> for Polyn<T> {}
// // impl<T: Clone + NumAssign> AddAssign<Polyn<T>> for Polyn<T> {}
// // impl<T: Clone + NumAssign> AddAssign<T> for Polyn<T> {}

// pub fn multi_inv<T: PrimeField>(values: &[T]) -> Vec<T> {
//   let mut partials = vec![T::one()];
//   debug_assert!(partials.len() > 0);
//   for i in 0..values.len() {
//     partials.push(
//       *partials.last().unwrap()
//         * (if values[i] != T::zero() {
//           values[i]
//         } else {
//           T::one()
//         }),
//     );
//   }

//   let mut inv = partials[partials.len() - 1].invert().unwrap();

//   let mut outputs = vec![T::zero(); values.len()];
//   for i in (0..values.len()).rev() {
//     outputs[i] = if values[i] != T::zero() {
//       partials[i] * inv
//     } else {
//       T::zero()
//     };

//     inv *= if values[i] != T::zero() {
//       values[i]
//     } else {
//       T::one()
//     };
//   }

//   outputs
// }

// #[test]
// fn test_multi_inv() {
//   use ff_utils::f7::F7;

//   let values: Vec<F7> = [1u64, 3, 2, 6, 4, 5].iter().map(|x| F7::from(*x)).collect();
//   let res = multi_inv(&values);
//   let answer: Vec<F7> = [1u64, 5, 4, 6, 2, 3].iter().map(|x| F7::from(*x)).collect();
//   assert_eq!(res, answer);

//   let values: Vec<F7> = [0u64, 1, 5, 4, 0, 6, 2, 3, 0]
//     .iter()
//     .map(|x| F7::from(*x))
//     .collect();
//   let res = multi_inv(&values);
//   let answer: Vec<F7> = [0u64, 1, 3, 2, 0, 6, 4, 5, 0]
//     .iter()
//     .map(|x| F7::from(*x))
//     .collect();
//   assert_eq!(res, answer);
// }

// // pub fn eval_poly_at<T: PrimeField>(polyn: &[T], x: T) -> T {
// //   let mut y = T::zero();
// //   let mut power_of_x = T::one();
// //   for polyn_coeff in polyn {
// //     let v = power_of_x * polyn_coeff;
// //     y.add_assign(&v);
// //     power_of_x.mul_assign(&x);
// //   }

// //   y
// // }

// // #[test]
// // fn test_eval_poly_at() {
// //   use ff_utils::f7::F7;

// //   // p(x) = 1 + 2x + x^3
// //   let p: Vec<F7> = [1u64, 2, 0, 1].iter().map(|c| F7::from(*c)).collect();
// //   let x = F7::from(2u64);
// //   // p(2) = 1 + 2*2 + 2^3 = 6
// //   let y = eval_poly_at(&p, x);
// //   let answer = F7::from(6u64);
// //   assert_eq!(y, answer);
// // }

// pub fn add_polys<T: PrimeField + ScalarOps>(a: &[T], b: &[T]) -> Vec<T> {
//   (0..max(a.len(), b.len()))
//     .map(
//       |i| if i < a.len() { a[i] } else { T::zero() } + if i < b.len() { b[i] } else { T::zero() },
//     )
//     .collect()
// }

// #[test]
// fn test_add_polys() {
//   use ff_utils::f7::F7;

//   // p1(x) = 4 + 2x + x^3
//   let p1: Vec<F7> = [4u64, 2, 0, 1].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p3(x) = 3 + 3x + 2x^2 + x^3
//   let p3 = add_polys(&p1, &p2);
//   let answer: Vec<F7> = [3u64, 3, 2, 1].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p3, answer);
// }

// pub fn sub_polys<T: PrimeField + ScalarOps>(a: &[T], b: &[T]) -> Vec<T> {
//   (0..max(a.len(), b.len()))
//     .map(
//       |i| if i < a.len() { a[i] } else { T::zero() } - if i < b.len() { b[i] } else { T::zero() },
//     )
//     .collect()
// }

// #[test]
// fn test_sub_polys() {
//   use ff_utils::f7::F7;

//   // p1(x) = 4 + 2x + x^3
//   let p1: Vec<F7> = [4u64, 2, 0, 1].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p3(x) = 5 + x + 5x^2 + x^3
//   let p3 = sub_polys(&p1, &p2);
//   let answer: Vec<F7> = [5u64, 1, 5, 1].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p3, answer);
// }

// pub fn mul_by_const<T: PrimeField>(a: &[T], scalar: T) -> Vec<T> {
//   a.iter().map(|x| *x * scalar).collect()
// }

// #[test]
// fn test_mul_by_const() {
//   use ff_utils::f7::F7;

//   // p(x) = 4 + 2x + x^3
//   let p: Vec<F7> = [4u64, 2, 0, 1].iter().map(|c| F7::from(*c)).collect();
//   let scalar = F7::from(5u64);
//   // q(x) = 5 * p(x)
//   let q = mul_by_const(&p, scalar);
//   let answer: Vec<F7> = [6u64, 3, 0, 5].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(q, answer);
// }

// pub fn reduction_poly<T: PrimeField>(a: &[T], n: usize) -> Vec<T> {
//   let mut o = vec![T::zero(); n];
//   for i in 0..a.len() {
//     o[i % n] += a[i];
//   }

//   o
// }

// #[test]
// fn test_reduction_poly() {
//   use ff_utils::f7::F7;

//   let p: Vec<F7> = [4u64, 2, 0, 1, 3, 2].iter().map(|c| F7::from(*c)).collect();
//   let q = reduction_poly(&p, 4);
//   let answer: Vec<F7> = [0u64, 4, 0, 1].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(q, answer);

//   let p: Vec<F7> = [4u64, 2].iter().map(|c| F7::from(*c)).collect();
//   let q = reduction_poly(&p, 4);
//   let answer: Vec<F7> = [4u64, 2, 0, 0].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(q, answer);
// }

// // recommend to use the DFT in the case of high degree
// pub fn mul_polys<T: PrimeField>(a: &[T], b: &[T]) -> Vec<T> {
//   let mut o = vec![T::zero(); a.len() + b.len() - 1];
//   for i in 0..a.len() {
//     for j in 0..b.len() {
//       o[i + j] = o[i + j] + a[i] * b[j];
//     }
//   }

//   o
// }

// #[test]
// fn test_mul_polys_low_degree() {
//   use ff_utils::f7::F7;

//   // p1(x) = 4 + 2x + x^3
//   let p1: Vec<F7> = [4u64, 2, 0, 1].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p3(x) = 3 + 2x + 3x^2 + 3x^3 + 1x^4 + 2x^5
//   let p3 = mul_polys(&p1, &p2);
//   let answer: Vec<F7> = [3u64, 2, 3, 3, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p3, answer);
// }

// pub fn poly_scale<T: PrimeField>(a: &[T], n: usize) -> Vec<T> {
//   let mut res = vec![T::zero(); n];
//   res.extend(a);
//   res
// }

// // recommend to use the DFT in the case of high degree
// pub fn div_polys<T: PrimeField>(a: &[T], b: &[T]) -> Vec<T> {
//   let zeros_len = b
//     .iter()
//     .rev()
//     .take_while(|&&coeff| coeff == T::zero())
//     .collect::<Vec<&T>>()
//     .len();
//   let b = b[..(b.len() - zeros_len)].to_vec();

//   assert!(a.len() >= b.len());
//   let mut c = a.to_vec();
//   let mut o = vec![];
//   let mut apos = a.len() - 1;
//   let bpos = b.len() - 1;
//   let diff = apos - bpos;
//   for d in (0..(diff + 1)).rev() {
//     let quot = c[apos] * b[bpos].invert().unwrap();
//     o.push(quot);
//     for i in (0..(bpos + 1)).rev() {
//       c[d + i] -= b[i] * quot;
//     }

//     apos -= 1;
//   }

//   o.reverse();
//   o
// }

// #[test]
// fn test_div_polys_low_degree() {
//   use ff_utils::f7::F7;

//   // p1(x) = 3 + 2x + 3x^2 + 3x^3 + 1x^4 + 2x^5
//   let p1: Vec<F7> = [3u64, 2, 3, 3, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p3(x) = 0
//   let p3 = div_polys(&p1, &p2);
//   let answer: Vec<F7> = [4u64, 2, 0, 1].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p3, answer);
// }

// #[test]
// #[should_panic(expected = "assertion failed: a.len() >= b.len()")]
// fn test_panic_div_polys_low_degree() {
//   use ff_utils::f7::F7;

//   // p1(x) = 3 + 2x + 3x^2 + 3x^3 + 1x^4 + 2x^5
//   let p1: Vec<F7> = [3u64, 2, 3, 3, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   div_polys(&p2, &p1);
// }

// // recommend to use the DFT in the case of high degree
// pub fn mod_polys<T: PrimeField>(a: &[T], b: &[T]) -> Vec<T> {
//   sub_polys(a, &mul_polys(b, &div_polys(a, b)))
//     .drain(0..(b.len() - 1))
//     .collect()
// }

// #[test]
// fn test_mod_polys_low_degree() {
//   use ff_utils::f7::F7;

//   // p1(x) = 5 + 4x + 3x^2 + 3x^3 + 1x^4 + 2x^5
//   let p1: Vec<F7> = [5u64, 4, 3, 3, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p3(x) = 2 + 2x
//   let p3 = mod_polys(&p1, &p2);
//   let answer: Vec<F7> = [2u64, 2].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p3, answer);

//   // p1(x) = 5 + 2x + 3x^2 + 3x^3 + 1x^4 + 2x^5
//   let p1: Vec<F7> = [5u64, 2, 3, 3, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p3(x) = 0
//   let p3 = mod_polys(&p1, &p2);
//   let answer: Vec<F7> = [2u64, 0].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p3, answer);

//   // p1(x) = 3 + 2x + 3x^2 + 3x^3 + 1x^4 + 2x^5
//   let p1: Vec<F7> = [3u64, 2, 3, 3, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p2(x) = 6 + x + 2x^2
//   let p2: Vec<F7> = [6u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   // p3(x) = 0
//   let p3 = mod_polys(&p1, &p2);
//   let answer: Vec<F7> = [0u64, 0].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p3, answer);
// }

// // Build a polynomial from a few coefficients
// pub fn sparse<T: PrimeField>(coeff_dict: HashMap<usize, T>) -> Vec<T> {
//   let mut o = vec![T::zero(); coeff_dict.keys().max().unwrap() + 1];
//   for k in coeff_dict.keys() {
//     o[*k] = coeff_dict[k];
//   }
//   o
// }

// #[test]
// fn test_sparse() {
//   use ff_utils::f7::F7;

//   let mut coeff: HashMap<usize, F7> = HashMap::new();
//   coeff.insert(1usize, F7::from(3u64));
//   coeff.insert(5usize, F7::from(1u64));
//   let p: Vec<F7> = sparse(coeff);
//   let answer: Vec<F7> = [0u64, 3, 0, 0, 0, 1].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(p, answer);

//   let mut coeff: HashMap<usize, F7> = HashMap::new();
//   coeff.insert(1usize, F7::from(3u64));
//   coeff.insert(5usize, F7::from(1u64));
//   coeff.insert(6usize, F7::from(0u64));
//   let p: Vec<F7> = sparse(coeff);
//   let answer: Vec<F7> = [0u64, 3, 0, 0, 0, 1, 0]
//     .iter()
//     .map(|c| F7::from(*c))
//     .collect();
//   assert_eq!(p, answer);
// }

// // Build a polynomial that returns 0 at all specified xs
// pub fn zpoly<T: PrimeField>(xs: &[T]) -> Vec<T> {
//   let mut root = vec![T::one()];
//   for i in 0..xs.len() {
//     root.push(T::zero());
//     for j in (0..(i + 1)).rev() {
//       root[j + 1] = root[j + 1] - root[j] * xs[i];
//     }
//   }

//   root.reverse();
//   root
// }

// #[test]
// fn test_zpoly() {
//   use ff_utils::f7::F7;

//   let xs: Vec<F7> = [0u64, 1, 2].iter().map(|c| F7::from(*c)).collect();
//   let ys = zpoly(&xs);
//   let answer: Vec<F7> = [0u64, 2, 4, 1].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(ys, answer);

//   let xs: Vec<F7> = [0u64, 3, 3].iter().map(|c| F7::from(*c)).collect();
//   let ys = zpoly(&xs);
//   let answer: Vec<F7> = [0u64, 2, 1, 1].iter().map(|c| F7::from(*c)).collect();
//   assert_eq!(ys, answer);

//   let xs: Vec<F7> = [1, 2, 3, 4, 5, 6].iter().map(|c| F7::from(*c)).collect();
//   let ys = zpoly(&xs);
//   let answer: Vec<F7> = [6u64, 0, 0, 0, 0, 0, 1]
//     .iter()
//     .map(|c| F7::from(*c))
//     .collect();
//   assert_eq!(ys, answer);

//   let xs: Vec<F7> = [1, 3, 2, 6, 5, 4].iter().map(|c| F7::from(*c)).collect();
//   let ys = zpoly(&xs);
//   assert_eq!(ys, answer);
// }

// // Given p+1 y values and x values with no errors, recovers the original
// // p+1 degree polynomial.
// // Lagrange interpolation works roughly in the following way.
// // 1. Suppose you have a set of points, eg. x = [1, 2, 3], y = [2, 5, 10]
// // 2. For each x, generate a polynomial which equals its corresponding
// //    y coordinate at that point and 0 at all other points provided.
// // 3. Add these polynomials together.
// pub fn lagrange_interp<T: PrimeField>(xs: &[T], ys: &[T]) -> Vec<T> {
//   // Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
//   let root = zpoly(xs);
//   assert_eq!(root.len(), ys.len() + 1);
//   // print(root)
//   // Generate per-value numerator polynomials, eg. for x=x2,
//   // (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
//   // polynomial back by each x coordinate
//   let nums: Vec<Polynomial<T, Coefficients>> = xs
//     .iter()
//     .map(|x| Polynomial::from_coeffs(div_polys(&root, &vec![T::zero() - x, T::one()])).unwrap())
//     .collect();
//   // Generate denominators by evaluating numerator polys at each x
//   let denoms: Vec<T> = (0..xs.len())
//     .map(|i| nums[i].evaluate_at(&worker, xs[i]))
//     .collect();
//   let inv_denoms = multi_inv(&denoms);
//   // Generate output polynomial, which is the sum of the per-value numerator
//   // polynomials rescaled to have the right y values
//   let mut b = vec![T::zero(); ys.len()];
//   for i in 0..xs.len() {
//     let yslice = ys[i] * inv_denoms[i];
//     for j in 0..ys.len() {
//       if nums[i][j] != T::zero() && ys[i] != T::zero() {
//         b[j] += nums[i][j] * yslice;
//       }
//     }
//   }

//   b
// }

// // Optimized poly evaluation for degree 4
// pub fn eval_quartic<T: PrimeField>(p: [T; 4], x: T) -> T {
//   let xsq = x * x;
//   let xcb = xsq * x;
//   p[0] + p[1] * x + p[2] * xsq + p[3] * xcb
// }

// // Optimized version of the above restricted to deg-4 polynomials
// pub fn multi_interp_4<T: PrimeField>(xsets: &[[T; 4]], ysets: &[[T; 4]]) -> Vec<[T; 4]> {
//   let mut data = vec![];
//   let mut inv_targets = vec![];
//   for key in 0..xsets.len() {
//     let xs = xsets[key];
//     let ys = ysets[key];
//     let x01 = xs[0] * xs[1];
//     let x02 = xs[0] * xs[2];
//     let x03 = xs[0] * xs[3];
//     let x12 = xs[1] * xs[2];
//     let x13 = xs[1] * xs[3];
//     let x23 = xs[2] * xs[3];
//     let eq0 = [
//       T::zero() - x12 * xs[3],
//       x12 + x13 + x23,
//       T::zero() - xs[1] - xs[2] - xs[3],
//       T::one(),
//     ];
//     let eq1 = [
//       T::zero() - x02 * xs[3],
//       x02 + x03 + x23,
//       T::zero() - xs[0] - xs[2] - xs[3],
//       T::one(),
//     ];
//     let eq2 = [
//       T::zero() - x01 * xs[3],
//       x01 + x03 + x13,
//       T::zero() - xs[0] - xs[1] - xs[3],
//       T::one(),
//     ];
//     let eq3 = [
//       T::zero() - x01 * xs[2],
//       x01 + x02 + x12,
//       T::zero() - xs[0] - xs[1] - xs[2],
//       T::one(),
//     ];
//     let e0 = eval_quartic(eq0, xs[0]);
//     let e1 = eval_quartic(eq1, xs[1]);
//     let e2 = eval_quartic(eq2, xs[2]);
//     let e3 = eval_quartic(eq3, xs[3]);
//     data.push([ys, eq0, eq1, eq2, eq3] as [[T; 4]; 5]);
//     inv_targets.extend([e0, e1, e2, e3]);
//   }

//   let inv_alls = multi_inv(&inv_targets);
//   let mut outputs = vec![];
//   for i in 0..data.len() {
//     let [ys, eq0, eq1, eq2, eq3] = data[i];
//     let inv_y0 = ys[0] * inv_alls[i * 4 + 0];
//     let inv_y1 = ys[1] * inv_alls[i * 4 + 1];
//     let inv_y2 = ys[2] * inv_alls[i * 4 + 2];
//     let inv_y3 = ys[3] * inv_alls[i * 4 + 3];
//     let mut o = [T::zero(); 4];
//     for i in 0..4 {
//       o[i] = eq0[i] * inv_y0 + eq1[i] * inv_y1 + eq2[i] * inv_y2 + eq3[i] * inv_y3;
//     }

//     outputs.push(o);
//   }

//   // assert!(o, [self.lagrange_interp_4(xs, ys) for xs, ys in zip(xsets, ysets)]);
//   outputs
// }
