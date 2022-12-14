use auto_ops::impl_op_ex;
use std::mem;

/// The modulus P := 2^64 - 2^32 + 1.
///
/// This is a prime number, and determines the base field we use for constraints.
const P: u64 = u64::wrapping_neg(1 << 32) + 1;
/// The order of the root of unity in our field.
///
/// This is more specifically the logarithm of the order, which is a more useful
/// quantity.
pub const ROOT_OF_UNITY_ORDER: u32 = 32;

/// Our base field.
///
/// This is a field with ~64 bits, and with some nice properties, like having
/// a root of unity of degree 2^32.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Field(u64);

impl Field {
    /// The zero element of this field.
    ///
    /// This is the identity for addition.
    pub const fn zero() -> Self {
        Field(0)
    }

    /// The one element of this field.
    ///
    /// This is the identity for multiplication.
    pub const fn one() -> Self {
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
    pub const fn generator() -> Self {
        Field(7)
    }

    /// A 2^32th root of unity.
    ///
    /// This is a value x such that x^(2^32) = 1.
    ///
    /// This is very useful, because it allows us to perform Number Theoretic Transforms
    /// over this field, which is a very quick method to multiply polynomials.
    pub const fn root_of_unity() -> Self {
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
impl_op_ex!(-|a: &Field| -> Field { Field::zero().sub(a) });
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

/// Represents an element of an extension field.
///
/// This is an extension of our base field, used when we need a larger field.
/// This field has ~192 bits, which we can leverage to reach 128 bits of security
/// in many places.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExtensionField {
    // Represents a polynomial a0 + a1 X + a2 X^2, in that order.
    data: [Field; 3],
}

impl ExtensionField {
    /// Create an extension field element from a base field element.
    const fn from_field(x: Field) -> Self {
        Self {
            data: [x, Field::zero(), Field::zero()],
        }
    }

    /// The zero element in this field.
    ///
    /// This is the identity for addition, and the same as the base field
    pub const fn zero() -> Self {
        Self::from_field(Field::zero())
    }

    /// The one element in this field.
    ///
    /// This is the identity for multiplication, and the same as the base field
    pub const fn one() -> Self {
        Self::from_field(Field::one())
    }

    /// A generator of this field.
    ///
    /// An element x such that every other non-zero element can be written as y = x^a,
    /// for some exponent a.
    ///
    /// Note that this is not the same generator as the base field, since this
    /// extension field is larger.
    pub const fn generator() -> Self {
        Self {
            data: [Field(4), Field(1), Field::zero()],
        }
    }

    /// Modify this element by adding another.
    ///
    /// Having this is useful because extension field elements are somewhat large.
    fn add_mut(&mut self, other: &Self) {
        self.data
            .iter_mut()
            .zip(other.data.iter())
            .for_each(|(a, b)| *a += b);
    }

    /// Return the result of adding two elements together.
    fn add(&self, other: &Self) -> Self {
        let mut out = *self;
        out.add_mut(other);
        out
    }

    /// Modify this element by subtracting another.
    fn sub_mut(&mut self, other: &Self) {
        self.data
            .iter_mut()
            .zip(other.data.iter())
            .for_each(|(a, b)| *a -= b)
    }

    /// Returning the result of subtracting another element from this one.
    fn sub(&self, other: &Self) -> Self {
        let mut out = *self;
        out.sub_mut(other);
        out
    }

    // Modify this element by multiplying it with another.
    fn mul_mut(&mut self, other: &Self) {
        // First, define c(X) = (a0 + a1 X + a2 X^2)(b0 + b1 X + b2 X^2).
        //
        // This is a polynomial of degree 4, so we only have 5 coefficients.
        // Conveniently, by going from c4 to c0, we can reuse the space in self
        // for c0, c1, and c2.
        let c4 = self.data[2] * other.data[2];
        let c3 = self.data[2] * other.data[1] + self.data[1] * other.data[2];
        self.data[2] = self.data[2] * other.data[0]
            + self.data[1] * other.data[1]
            + self.data[0] * other.data[2];
        self.data[1] = self.data[1] * other.data[0] + self.data[0] * other.data[1];
        self.data[0] *= other.data[0];

        // The polynomial defining this field is f(X) = X^3 + 3. This means that
        // X^3 = -3 mod f. So, c(x) mod 4 = (c0 + c1 X + c2 X^2) - 3(c3 + c4 X).

        // It's faster to just subtract 3 times, rather than multiplying by -3.
        for _ in 0..3 {
            self.data[1] -= c4;
            self.data[0] -= c3;
        }
    }

    /// Return the result of multiplying this element with another.
    fn mul(&self, other: &Self) -> Self {
        let mut out = *self;
        out.mul_mut(other);
        out
    }

    /// Apply the map of frobenius to this element.
    ///
    /// This map sends x to x^P, with P the order of the base field.
    /// This is an essential component of our inversion algorithm.
    fn frobenius(&mut self) {
        // Notice that since P * x = 0, we have (a + b)^P = a^P + b^P.
        // (To see why, use the binomial expansion, and notice that all but the
        // first and last terms have P as a factor).
        //
        // This means that (a0 + a1 X + a2 X^2)^P = a0^P + a1^P X^P + a2^P X^2P.
        // Now, in the base field, because there are P - 1 nonzero elements,
        // we have a^P = P. This allows us to write:
        //  (a0 + a1 X + a2 X^2)^P = a0 + a1 X^P + a2 X^2P
        // now we just need to determine how to handle those powers of X.

        // Write e = 3 q + r, then since X^3 = -3, we have:
        //   X^e = X^(3q + r) = -3^q X^r
        // Thus, X^P = -3^(P // 3) X^(P % 3), and X^2P = (X^P)^2.
        //
        // We could do these calculations by hand, by computer, or just get
        // sage to do them for us.
        //
        // We have:
        //  X^P = 4294967295 * X
        //  X^2P = 18446744065119617025 * X^2
        self.data[1] *= Field(4294967295);
        self.data[2] *= Field(18446744065119617025);
    }

    /// Replace this field element with its inverse.
    ///
    /// This is an element y such that x * y = 1.
    pub fn invert(&mut self) {
        // The exponent r := (P^3 - 1) / (P - 1) is such that a^r is in the base field.
        // To see why, notice that (a^r)^(P - 1) = 1. This means that the element
        // a^r must belong to a multiplicative subgroup of size P - 1, which is
        // just the non-zero elements of the base field. This is because
        // (P^3 - 1) = (P - 1)(P^2 + P + 1).
        //
        // Our inversion strategy is thus to calculate:
        //
        //  (a^r)^-1 * a^(r - 1)
        //
        // With inversion done in the base field.
        //
        // Notice that r - 1 = P^2 + P, so our strategy for calculating a^(r - 1) is:
        //  a^(r - 1) = (a * a^P)^P.

        // First, put a^(r - 1) into self.
        let mut a = *self;
        self.frobenius();
        *self *= a;
        self.frobenius();
        // Now, put a^r into a.
        a *= &*self;
        // Calculate (a^r)^-1, knowing that we're in the base field.
        let a_inv = a.data[0].inverse();
        // Multiply this with a^(r - 1), to get our inverse.
        self.data.iter_mut().for_each(|x| *x *= a_inv);
    }

    /// Return the inverse of a field element.
    ///
    /// This is an element y such that x * y = 1.
    pub fn inverse(&self) -> ExtensionField {
        let mut out = *self;
        out.invert();
        out
    }
}

// Now, use all of the functions we've defined inside of the struct to implement
// all of the associated operators.
//
// We use macros in order to generate implementations for both plain values and
// references, which is quite convenient.
impl_op_ex!(+ |a: &ExtensionField, b: &ExtensionField| -> ExtensionField { a.add(b) });
impl_op_ex!(+= |a: &mut ExtensionField, b: &ExtensionField| { a.add_mut(b) });
impl_op_ex!(-|a: &ExtensionField| -> ExtensionField { ExtensionField::zero().sub(a) });
impl_op_ex!(-|a: &ExtensionField, b: &ExtensionField| -> ExtensionField { a.sub(b) });
impl_op_ex!(-= |a: &mut ExtensionField, b: &ExtensionField| { a.sub_mut(b) });
impl_op_ex!(*|a: &ExtensionField, b: &ExtensionField| -> ExtensionField { a.mul(b) });
impl_op_ex!(*= |a: &mut ExtensionField, b: &ExtensionField| { a.mul_mut(b) });
impl_op_ex!(/|a: &ExtensionField, b: &ExtensionField| -> ExtensionField { a * b.inverse()});
impl_op_ex!(/= |a: &mut ExtensionField, b: &ExtensionField| { *a *= b.inverse() });

// Very commonly, we want to convert elements of the base field into the extension field.
impl From<Field> for ExtensionField {
    fn from(x: Field) -> Self {
        Self::from_field(x)
    }
}

#[cfg(test)]
pub mod generators {
    use super::*;
    use proptest::prelude::*;

    /// A strategy to generate arbitrary field elements.
    pub fn arb_field() -> impl Strategy<Value = Field> {
        // While not cryptographically secure to reduce such a small element,
        // this is more than enough for our testing purposes.
        any::<u64>().prop_map(|x| Field(x % P))
    }

    /// A strategy to generate non-zero field elements.
    pub fn arb_field_non_zero() -> impl Strategy<Value = Field> {
        arb_field().prop_filter("field elements must be non-zero", |x| x != &Field::zero())
    }

    prop_compose! {
        /// A strategy to generate arbitrary extension field elements.
        pub fn arb_extension()(a in arb_field(), b in arb_field(), c in arb_field()) -> ExtensionField {
            ExtensionField { data: [a, b, c] }
        }
    }

    /// A strategy to generate arbitrary extension field elements.
    pub fn arb_extension_non_zero() -> impl Strategy<Value = ExtensionField> {
        arb_extension().prop_filter("extension field elements must be non-zero", |x| {
            x != &ExtensionField::zero()
        })
    }
}

#[cfg(test)]
pub mod test {
    use super::generators::*;
    use super::*;

    use proptest::prelude::*;

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

    proptest! {
        #[test]
        fn test_extension_addition_commutative(a in arb_extension(), b in arb_extension()) {
            assert_eq!(a + b, b + a);
        }
    }

    proptest! {
        #[test]
        fn test_extension_addition_identity(a in arb_extension()) {
            assert_eq!(a + ExtensionField::zero(), a);
        }
    }

    proptest! {
        #[test]
        fn test_extension_subtraction_with_self_is_zero(a in arb_extension()) {
            assert_eq!(a - a, ExtensionField::zero());
        }
    }

    proptest! {
        #[test]
        fn test_extension_multiplication_commutative(a in arb_extension(), b in arb_extension()) {
            assert_eq!(a * b, b * a);
        }
    }

    proptest! {
        #[test]
        fn test_extension_multiplication_identity(a in arb_extension()) {
            assert_eq!(a * ExtensionField::one(), a);
        }
    }

    proptest! {
        #[test]
        fn test_extension_multiplication_distributive(a in arb_extension(), b in arb_extension(), c in arb_extension()) {
            assert_eq!(a * (b + c), a * b + a * c);
        }
    }

    proptest! {
        #[test]
        fn test_extension_division_by_self_is_one(a in arb_extension_non_zero()) {
            assert_eq!(a / a, ExtensionField::one());
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

    #[test]
    fn test_extension_one_plus_one_is_two() {
        let x = ExtensionField {
            data: [Field::from(1), Field::from(2), Field::from(3)],
        };
        let y = ExtensionField {
            data: [Field::from(2), Field::from(4), Field::from(6)],
        };
        assert_eq!(x + x, y);
    }

    #[test]
    fn test_x_times_x_is_x_squared() {
        let x = ExtensionField {
            data: [Field::zero(), Field::one(), Field::zero()],
        };
        let x2 = ExtensionField {
            data: [Field::zero(), Field::zero(), Field::one()],
        };
        assert_eq!(x * x, x2);
    }
}
