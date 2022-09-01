use auto_ops::impl_op_ex;

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

    #[test]
    fn test_one_plus_one_is_two() {
        assert_eq!(Field::from(1) + Field::from(1), Field::from(2));
    }

    #[test]
    fn test_subtraction_examples() {
        assert_eq!(Field::one() - Field::one(), Field::zero());
        assert_eq!(Field::zero() - Field::one(), Field::from(P - 1));
    }
}
