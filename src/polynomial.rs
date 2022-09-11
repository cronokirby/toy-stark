use itertools::iterate;

use crate::field::{Field, ROOT_OF_UNITY_ORDER};
use crate::util;

/// Get a root of unity with a certain degree.
///
/// This degree must be low enough, otherwise our field doesn't contain that root.
fn root_of_unity(lg_degree: u32) -> Field {
    assert!(lg_degree <= ROOT_OF_UNITY_ORDER);
    // Every squaring halves the degree of the polynomial.
    // This means that to get to the right degree, we need to square a number
    // of times based on the difference between the order of our root and the desired order.
    let mut out = Field::root_of_unity();
    for _ in 0..(ROOT_OF_UNITY_ORDER - lg_degree) {
        out *= out;
    }
    out
}

/// A single iteration of our NTT algorithm.
fn vector_ntt_iter<const INVERSE: bool>(data: &mut [Field], power: Field, skip: usize) {
    let n = data.len();
    for start in (0..)
        .map(|start| 2 * start * skip)
        .take_while(|&start| start < n)
    {
        let mut w = Field::one();
        for i in start..start + skip {
            let j = i + skip;
            let (a, b) = (data[i], data[j]);
            // These two transformations are almost inverses, with w being the inverse
            // in one branch. The difference is that in the INVERSE branch,
            // we end up with twice the result we need. We correct for this later.
            if INVERSE {
                data[i] = a + b;
                data[j] = (a - b) * w;
            } else {
                // Hopefully the compiler does common subexpression elimination here.
                data[i] = a + w * b;
                data[j] = a - w * b;
            }
            w *= power;
        }
    }
}

/// Calculate the NTT in place.
///
/// Note that this will accept and return polynomials in reverse bit order,
/// and evaluations in normal bit order.
fn vector_ntt<const INVERSE: bool>(data: &mut [Field]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    assert!(n.is_power_of_two());
    let lg_n = util::lg(n);
    // This conversion is fine since we already need to check that lg_n is less than 32.
    let lg_n_usize = lg_n as usize;

    // See CLRS, Section 30.3 for the algorithm here, a bit tricky to explain otherwise.

    // Create powers of the root, up to 1
    let mut root = root_of_unity(lg_n);
    if INVERSE {
        root = root.inverse();
    }
    let w_powers: Vec<Field> = iterate(root, |x| x * x).take(lg_n_usize).collect();

    // We reverse the order we operate in when calculating the inverse transform.
    if INVERSE {
        let skips = iterate(1 << (lg_n_usize - 1), |x| x >> 1).take(lg_n_usize);
        w_powers
            .into_iter()
            .zip(skips)
            .for_each(|(power, skip)| vector_ntt_iter::<INVERSE>(data, power, skip));
    } else {
        let skips = iterate(1, |x| x << 1).take(lg_n_usize);
        w_powers
            .into_iter()
            .rev()
            .zip(skips)
            .for_each(|(power, skip)| vector_ntt_iter::<INVERSE>(data, power, skip));
    }

    // We also need to scale the result appropriately when doing the inverse transform.
    if INVERSE {
        let adjust = Field::from(1 << lg_n).inverse();
        for di in data.iter_mut() {
            *di *= adjust;
        }
    };
}

/// Represents a polynomial over our base field.
///
/// We require that the length of this polynomial be a power of two.
#[derive(Clone, Debug, PartialEq)]
struct Polynomial {
    // The length of this vector, N, is a power of two.
    //
    // This contains a polynomial of length N (degree N - 1).
    //
    // The coefficients are in a bit of an odd order. Given a polynomial
    //   a0 + a1 X + ... + a(N-1) X^(N-1)
    // we have coefficients[i] = a(reverse(i, lg N)), with the bits of the index,
    // considered as a lg N bit number, reversed.
    //
    // Another way of looking at it is that the first half contains the polynomial
    // whose first bit is 0, the second half those whose first bit is 1.
    // In each half, we split into two new halves, now based on the second bit,
    // etc.
    coefficients: Vec<Field>,
}

impl Polynomial {
    /// Calculate the number theoretic transform of this polynomial.
    ///
    /// This will return a vector containing the evaluations of the polynomial
    /// at roots of unity.
    pub fn ntt(&self) -> NTTPolynomial {
        let mut data = self.coefficients.clone();
        vector_ntt::<false>(&mut data);
        NTTPolynomial { evaluations: data }
    }

    /// Evaluate this polynomial at a given point.
    pub fn evaluate(&self, x: Field) -> Field {
        let n = self.coefficients.len();
        let lg_n = util::lg(n);
        let mut acc = Field::zero();
        // Use Horner's method, except that we need to reverse the bit order.
        for i in (0..n).rev().map(|i| util::reverse(i, lg_n)) {
            acc = acc * x + self.coefficients[i];
        }
        acc
    }
}

/// Represents the number theoretic transform of a polynomial.
///
/// The transform of a polynomial f is given by evaluating f at roots of unity:
///   f(w^0), f(w^1), f(w^2), f(w^3), ...
///
/// We require that the length of these evaluations be a power of two.
#[derive(Clone, Debug, PartialEq)]
struct NTTPolynomial {
    // The length of this vector, N, is a power of two.
    //
    // The vector contains f(w^i) for i in [N], with w an Nth root of unity.
    // The order of evaluations is f(w^0), f(w^1), f(w^2), ...
    evaluations: Vec<Field>,
}

impl NTTPolynomial {
    /// Create a polynomial from the evaluations over a domain.
    ///
    /// The number of evaluations, N, must be a power of two.
    ///
    /// Given a polynomial f, the list contains the evaluations:
    ///   f(w^0), f(w^1), f(w^2), ...
    /// for an Nth root of unity w.
    pub fn new(evaluations: &[Field]) -> Self {
        let n = evaluations.len();
        assert!(
            n.is_power_of_two(),
            "evaluation length ({n}) is not a power of two."
        );

        NTTPolynomial {
            evaluations: evaluations.to_owned(),
        }
    }

    /// Interpolate the evaluations, recovering the underlying polynomial.
    pub fn interpolate(&self) -> Polynomial {
        let mut data = self.evaluations.clone();
        vector_ntt::<true>(&mut data);
        Polynomial { coefficients: data }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::field::generators::arb_field;
    use proptest::{collection::vec, prelude::*};

    fn arb_polynomial(lg_degree: u32) -> impl Strategy<Value = Polynomial> {
        vec(arb_field(), 1 << lg_degree).prop_map(|coefficients| Polynomial { coefficients })
    }

    proptest! {
        #[test]
        fn test_double_ntt_is_identity(f in arb_polynomial(8)) {
            let f_hat = f.ntt();
            let f_prime = f_hat.interpolate();
            assert_eq!(f, f_prime);
        }
    }

    proptest! {
        #[test]
        fn test_ntt_vs_manual_evaluation(f in arb_polynomial(2)) {
            let f_hat = f.ntt();
            let w = root_of_unity(2);
            for (root, &f_root) in iterate(Field::one(), |x| x * w).zip(f_hat.evaluations.iter()) {
                assert_eq!(f.evaluate(root), f_root);
            }
        }
    }

    #[test]
    fn test_ntt_of_x() {
        let x = Polynomial {
            coefficients: vec![Field::zero(), Field::zero(), Field::one(), Field::zero()],
        };
        let w = root_of_unity(2);
        let expected_ntt = NTTPolynomial {
            evaluations: vec![Field::one(), w, -Field::one(), -w],
        };
        assert_eq!(x.ntt(), expected_ntt);
    }

    #[test]
    fn test_polynomial_evaluation_example() {
        let f = Polynomial {
            coefficients: vec![
                Field::from(1),
                Field::from(4),
                Field::from(2),
                Field::from(3),
            ],
        };
        let x = Field::from(7);
        let expected = f.coefficients[0]
            + f.coefficients[1] * x * x
            + f.coefficients[2] * x
            + f.coefficients[3] * x * x * x;
        assert_eq!(f.evaluate(x), expected);
    }
}
