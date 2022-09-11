use itertools::iterate;

use crate::field::{Field, ROOT_OF_UNITY_ORDER};
use crate::util;

fn root_of_unity(lg_degree: u32) -> Field {
    assert!(lg_degree <= ROOT_OF_UNITY_ORDER);
    let mut out = Field::root_of_unity();
    for _ in 0..(ROOT_OF_UNITY_ORDER - lg_degree) {
        out *= out;
    }
    out
}

fn vector_ntt_iter<const INVERSE: bool>(data: &mut [Field], power: Field, skip: usize) {
    let n = data.len();
    let mut w = Field::one();
    for start in (0..)
        .map(|start| 2 * start * skip)
        .take_while(|&start| start < n)
    {
        for i in start..start + skip {
            let j = i + skip;
            let (a, b) = (data[i], data[j]);
            if INVERSE {
                data[i] = a + b;
                data[j] = (a - b) * w;
            } else {
                data[i] = a + w * b;
                data[j] = a - w * b;
            }
        }
        w *= power;
    }
}

fn vector_ntt<const INVERSE: bool>(data: &mut [Field]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    assert!(n.is_power_of_two());
    let lg_n = util::lg(n);
    // This conversion is fine since we already need to check that lg_n is less than 32.
    let lg_n_usize = lg_n as usize;

    // Create powers of the root, up to 1
    let mut root = root_of_unity(lg_n);
    if INVERSE {
        root = root.inverse();
    }
    let w_powers: Vec<Field> = iterate(root, |x| x * x).take(lg_n_usize).collect();

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

    if INVERSE {
        let adjust = Field::from(1 << lg_n).inverse();
        for di in data.iter_mut() {
            *di *= adjust;
        }
    };
}

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
    pub fn ntt(&self) -> NTTPolynomial {
        let mut data = self.coefficients.clone();
        vector_ntt::<false>(&mut data);
        NTTPolynomial { evaluations: data }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct NTTPolynomial {
    // The length of this vector, N, is a power of two.
    //
    // The vector contains f(w^i) for i in [N], with w an Nth root of unity.
    // The order of the evaluations is such that evaluations[i] is w^reverse(i, lg N).
    //
    // See the comments on `Polynomial` for some discussion on the choice of this convention.
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

    fn arb_polynomial() -> impl Strategy<Value = Polynomial> {
        vec(arb_field(), 1 << 8).prop_map(|coefficients| Polynomial { coefficients })
    }

    proptest! {
        #[test]
        fn test_double_ntt_is_identity(f in arb_polynomial()) {
            dbg!(&f);
            let f_hat = f.ntt();
            dbg!(&f_hat);
            let f_prime = f_hat.interpolate();
            assert_eq!(f, f_prime);
        }
    }

    #[test]
    fn test_ntt_of_x() {
        let x = Polynomial {
            coefficients: vec![Field::zero(), Field::one()],
        };
        let expected_ntt = NTTPolynomial {
            evaluations: vec![Field::one(), -Field::one()],
        };
        assert_eq!(x.ntt(), expected_ntt);
    }
}
