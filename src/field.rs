use auto_ops::impl_op_ex;
use std::mem;

/// The modulus P := 2^64 - 2^32 + 1.
///
/// This is a prime number, and determines the base field we use for constraints.
const P: u64 = u64::wrapping_neg(1 << 32) + 1;

#[derive(Clone, Copy, Debug, PartialEq)]
struct Field(u64);

impl Field {
    /// Create the zero element of this field.
    ///
    /// This is the identity for addition.
    pub fn zero() -> Self {
        Field(0)
    }

    /// Create the one element of this field.
    ///
    /// This is the identity for multiplication.
    pub fn one() -> Self {
        Field(1)
    }

    /// A generator of this field.
    ///
    /// This is an element x such that every other non-zero element y = x^a,
    /// for some exponent a.
    ///
    /// Sometimes you need to ensure that two sets of element don't intersect,
    /// but you don't want to mess with the multiplicative structure of those elements:
    /// multiplying by the generator allows separating the two sets in this way.
    pub fn generator() -> Self {
        Field(7)
    }

    /// A 2^32th root of unity.
    ///
    /// This is a value x such that x^(2^32) = 1.
    ///
    /// This is very useful, because it allows us to perform Number Theoretic Transforms
    /// over this field, which is a very quick method to multiply polynomials.
    pub fn root_of_unity() -> Self {
        Field(20033703337)
    }

    /// Add another field element to this one.
    fn add_mut(&mut self, other: &Field) {
        // Because a, b are at most P - 1, the result of addition is at most
        // 2P - 2, which can be made < P with a single subtraction. Our strategy
        // Is thus to perform the addition, and then conditionally subtract P
        // if the result is greater than P.
        let (addition, overflow) = self.0.overflowing_add(other.0);
        let (with_p_subtracted, underflow) = addition.overflowing_sub(P);
        // If we overflowed the 64 bits with the addition, we know that we need to
        // subtract P for sure. Otherwise, we choose the subtraction if it didn't underflow.
        self.0 = if overflow || !underflow {
            with_p_subtracted
        } else {
            addition
        };
    }

    /// Return the result of adding this field element to another.
    fn add(&self, other: &Field) -> Field {
        // Simply reuse the mutable addition defined earlier.
        let mut out = *self;
        out.add_mut(other);
        out
    }

    /// Subtract another field element from this one.
    fn sub_mut(&mut self, other: &Field) {
        // The strategy here is to perform the subtraction, and then add back P
        // if an underflow happened. This will always result in a value < P.
        // If no underflow happened, this is clear, since both values were < P.
        // If an underflow happened, the largest result we can have is -1. Adding
        // P gives us P - 1, which is < P, so everything works.
        let (subtraction, underflow) = self.0.overflowing_sub(other.0);
        self.0 = if underflow {
            subtraction.wrapping_add(P)
        } else {
            subtraction
        };
    }

    /// Return the result of subtracting another field element from this one.
    fn sub(&self, other: &Field) -> Field {
        let mut out = *self;
        out.sub_mut(other);
        out
    }

    /// Reduce a 128 bit integer into a field element.
    fn reduce_128(x: u128) -> Field {
        // We exploit special properties of the field.
        //
        // First, 2^64 = 2^32 - 1 mod P.
        //
        // Second, 2^96 = 2^32(2^32 - 1) = 2^64 - 2^32 = -1 mod P.
        //
        // Thus, if we write a 128 bit integer x as:
        //     x = c 2^96 + b 2^64 + a
        // We have:
        //     x = b (2^32 - 1) + (a - c) mod P
        // And this expression will be our strategy for performing the reduction.
        let a = x as u64;
        let b = ((x >> 64) & 0xFF_FF_FF_FF) as u64;
        let c = (x >> 96) as u64;

        // While we lean on existing code, we need to be careful because some of
        // these types are partially reduced.
        //
        // First, if we look at a - c, the end result with our field code can
        // be any 64 bit value (consider c = 0). We can also make the same assumption
        // for (b << 32) - b. The question then becomes, is Field(x) + Field(y)
        // ok even if both x and y are arbitrary u64 values?
        //
        // Yes. Even if x and y have the maximum value, a single subtraction of P
        // would suffice to make their sum < P. Thus, our strategy for field addition
        // will always work.
        (Field(a) - Field(c)) + Field((b << 32) - b)
    }

    /// Return the multiplication of this field element with another.
    fn mul(&self, other: &Field) -> Field {
        // 128 bit multiplications don't optimize too badly actually.
        Field::reduce_128(u128::from(self.0) * u128::from(other.0))
    }

    /// Modify this element by multiplying it with another.
    fn mul_mut(&mut self, other: &Field) {
        *self = self.mul(other);
    }

    /// Return self / 2 mod P.
    ///
    /// This is a utility function, mainly used when inverting field elements.
    fn half_mod_p(&self) -> Field {
        if self.0 & 1 == 0 {
            Field(self.0 >> 1)
        } else {
            let (addition, carry) = self.0.overflowing_add(P);
            Field((addition >> 1) | (u64::from(carry) << 63))
        }
    }

    /// Return the inverse of this field element.
    ///
    /// This is an element x such that x * self = 1.
    pub fn inverse(&self) -> Field {
        let mut a = self.0;
        let mut u = Field(1u64);
        let mut b = P;
        let mut v = Field(0u64);

        while a != 0 {
            if a & 1 == 0 {
                a >>= 1;
                u = u.half_mod_p();
            } else {
                if a < b {
                    mem::swap(&mut a, &mut b);
                    mem::swap(&mut u, &mut v);
                }
                a = (a - b) >> 1;
                u = (u - v).half_mod_p();
            }
        }

        v
    }
}

// Now, use all of the functions we've defined inside of the struct to implement
// all of the associated operators.
//
// We use macros in order to generate implementations for both plain values and
// references, which is quite convenient.
impl_op_ex!(+ |a: &Field, b: &Field| -> Field { a.add(b) });
impl_op_ex!(+= |a: &mut Field, b: &Field| { a.add_mut(b) });
impl_op_ex!(-|a: &Field, b: &Field| -> Field { a.sub(b) });
impl_op_ex!(-= |a: &mut Field, b: &Field| { a.sub_mut(b) });
impl_op_ex!(*|a: &Field, b: &Field| -> Field { a.mul(b) });
impl_op_ex!(*= |a: &mut Field, b: &Field| { a.mul_mut(b) });
impl_op_ex!(/|a: &Field, b: &Field| -> Field { a.mul(&b.inverse()) });
impl_op_ex!(/= |a: &mut Field, b: &Field| { a.mul_mut(&b.inverse()) });

// We might want to create the field from u64s, for example, when deserializing.
impl From<u64> for Field {
    fn from(x: u64) -> Self {
        Field(x % P)
    }
}

// Also useful to be able to serialize fields as u64s.
impl From<Field> for u64 {
    fn from(x: Field) -> Self {
        x.0
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use proptest::prelude::*;

    /// A strategy to generate arbitrary field elements.
    fn arb_field() -> impl Strategy<Value = Field> {
        // While not cryptographically secure to reduce such a small element,
        // this is more than enough for our testing purposes.
        any::<u64>().prop_map(|x| Field(x % P))
    }

    /// A strategy to generate non-zero field elements.
    fn arb_field_non_zero() -> impl Strategy<Value = Field> {
        arb_field().prop_filter("field elements must be non-zero", |x| x != &Field::zero())
    }

    proptest! {
        #[test]
        fn test_addition_commutative(a in arb_field(), b in arb_field()) {
            assert_eq!(a + b, b + a);
        }
    }

    proptest! {
        #[test]
        fn test_addition_identity(a in arb_field()) {
            assert_eq!(a + Field::zero(), a);
        }
    }

    proptest! {
        #[test]
        fn test_subtraction_with_self_is_zero(a in arb_field()) {
            assert_eq!(a - a, Field::zero());
        }
    }

    proptest! {
        #[test]
        fn test_multiplication_commutative(a in arb_field(), b in arb_field()) {
            assert_eq!(a * b, b * a);
        }
    }

    proptest! {
        #[test]
        fn test_multiplication_identity(a in arb_field()) {
            assert_eq!(a * Field::one(), a);
        }
    }

    proptest! {
        #[test]
        fn test_division_by_self_is_one(a in arb_field_non_zero()) {
            assert_eq!(a / a, Field::one());
        }
    }

    #[test]
    fn test_one_plus_one_is_two() {
        assert_eq!(Field::from(1) + Field::from(1), Field::from(2));
    }

    #[test]
    fn test_subtraction_examples() {
        assert_eq!(Field::one() - Field::one(), Field::zero());
        assert_eq!(Field::zero() - Field::one(), Field::from(P - 1));
    }

    #[test]
    fn test_multiplication_examples() {
        assert_eq!(Field::from(2) * Field::from(2), Field::from(4));
        assert_eq!(Field::from(P - 1) * Field::from(P - 1), Field::one());
    }

    #[test]
    fn test_inverse_one_is_one() {
        assert_eq!(Field::one().inverse(), Field::one());
    }

    #[test]
    fn test_inverse_minus_one_is_one() {
        assert_eq!(Field::from(P - 1).inverse(), Field::from(P - 1));
    }

    #[test]
    fn test_root_of_unity() {
        let mut out = Field::root_of_unity();
        for _ in 0..31 {
            out = out * out;
        }
        assert_ne!(out, Field::one());
        out = out * out;
        assert_eq!(out, Field::one());
    }
}
