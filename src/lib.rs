#![allow(private_bounds)]

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
use alloy_primitives::Signed;
#[cfg(all(feature = "malachite", feature = "alloy"))]
use alloy_primitives::Uint;
#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
use bigdecimal::num_bigint::BigUint;
#[cfg(feature = "bigdecimal")]
use bigdecimal::{
    BigDecimal,
    num_bigint::{BigInt, Sign},
};
#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
use malachite::Integer;
#[cfg(feature = "malachite")]
use malachite::{Natural, Rational, num::basic::traits::One};

pub trait Convert {
    fn convert_to<T>(&self) -> T
    where
        Self: Convertable<T>,
    {
        self.convertable_to()
    }
}

impl<A> Convert for A {}

pub trait Convertable<T> {
    fn convertable_to(&self) -> T;
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<BigDecimal> for Uint<BITS, LIMBS> {
    fn convertable_to(&self) -> BigDecimal {
        BigDecimal::from(BigInt::from_bytes_le(Sign::Plus, &self.to_le_bytes_vec()))
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Uint<BITS, LIMBS>> for BigDecimal {
    fn convertable_to(&self) -> Uint<BITS, LIMBS> {
        Uint::<BITS, LIMBS>::from_le_slice(&self.as_bigint_and_scale().0.to_bytes_le().1)
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<BigDecimal> for Signed<BITS, LIMBS> {
    fn convertable_to(&self) -> BigDecimal {
        let sign = if self.is_negative() {
            Sign::Minus
        } else {
            Sign::Plus
        };

        let unsigned_value = self.unsigned_abs();
        BigDecimal::from(BigInt::from_bytes_le(
            sign,
            &unsigned_value.to_le_bytes_vec(),
        ))
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Signed<BITS, LIMBS>> for BigDecimal {
    fn convertable_to(&self) -> Signed<BITS, LIMBS> {
        let (sign, bytes) = self.as_bigint_and_scale().0.to_bytes_le();
        let unsigned = Uint::<BITS, LIMBS>::from_le_slice(&bytes);
        match sign {
            bigdecimal::num_bigint::Sign::Minus => Signed::from_raw(unsigned).wrapping_neg(),
            _ => Signed::from_raw(unsigned),
        }
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Integer> for Signed<BITS, LIMBS> {
    fn convertable_to(&self) -> Integer {
        Integer::from_twos_complement_limbs_asc(self.as_limbs())
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Signed<BITS, LIMBS>> for Integer {
    fn convertable_to(&self) -> Signed<BITS, LIMBS> {
        let limbs = self.to_twos_complement_limbs_asc();
        let limbs_asc_exact: [u64; LIMBS] = if limbs.len() >= LIMBS {
            limbs[0..LIMBS].try_into().unwrap()
        } else {
            let mut l = vec![0u64; LIMBS];
            for (i, limb) in limbs.iter().enumerate() {
                l[i] = *limb;
            }
            if self < &Integer::from(0) {
                for i in limbs.len()..LIMBS {
                    l[i] = u64::MAX;
                }
            }
            l.try_into().unwrap()
        };
        Signed::<BITS, LIMBS>::from_limbs(limbs_asc_exact)
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Uint<BITS, LIMBS>> for BigUint {
    fn convertable_to(&self) -> Uint<BITS, LIMBS> {
        Uint::<BITS, LIMBS>::from_be_slice(&self.to_bytes_be())
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<BigUint> for Uint<BITS, LIMBS> {
    fn convertable_to(&self) -> BigUint {
        BigUint::from_bytes_le(&self.to_le_bytes_vec())
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Signed<BITS, LIMBS>> for BigInt {
    fn convertable_to(&self) -> Signed<BITS, LIMBS> {
        let (sign, bytes) = self.to_bytes_be();
        let unsigned = Uint::<BITS, LIMBS>::from_be_slice(&bytes);
        match sign {
            bigdecimal::num_bigint::Sign::Minus => Signed::from_raw(unsigned).wrapping_neg(),
            _ => Signed::from_raw(unsigned),
        }
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<BigInt> for Signed<BITS, LIMBS> {
    fn convertable_to(&self) -> BigInt {
        let sign = match (self.is_zero(), self.sign()) {
            (true, _) => bigdecimal::num_bigint::Sign::NoSign,
            (_, alloy_primitives::Sign::Negative) => bigdecimal::num_bigint::Sign::Minus,
            (_, alloy_primitives::Sign::Positive) => bigdecimal::num_bigint::Sign::Plus,
        };
        let unsigned_value = self.unsigned_abs();
        BigInt::from_bytes_le(sign, &unsigned_value.to_le_bytes_vec())
    }
}

#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
impl Convertable<Natural> for BigUint {
    fn convertable_to(&self) -> Natural {
        use alloy_primitives::U256;

        let v: U256 = self.convertable_to();
        v.convertable_to()
    }
}

#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
impl Convertable<BigUint> for Natural {
    fn convertable_to(&self) -> BigUint {
        use alloy_primitives::U256;

        let v: U256 = self.convertable_to();
        v.convertable_to()
    }
}

#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
impl Convertable<BigInt> for Natural {
    fn convertable_to(&self) -> BigInt {
        use alloy_primitives::{I256, U256};

        let v: U256 = self.convertable_to();
        I256::from(v).convertable_to()
    }
}

#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
impl Convertable<Integer> for BigInt {
    fn convertable_to(&self) -> Integer {
        use alloy_primitives::I256;

        let v: I256 = self.convertable_to();
        v.convertable_to()
    }
}

#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
impl Convertable<BigInt> for Integer {
    fn convertable_to(&self) -> BigInt {
        use alloy_primitives::I256;

        let v: I256 = self.convertable_to();
        v.convertable_to()
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Natural> for Uint<BITS, LIMBS> {
    fn convertable_to(&self) -> Natural {
        Natural::from_limbs_asc(self.as_limbs())
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Rational> for Uint<BITS, LIMBS> {
    fn convertable_to(&self) -> Rational {
        Rational::from_naturals(self.convertable_to(), Natural::ONE)
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convertable<Uint<BITS, LIMBS>> for Natural {
    fn convertable_to(&self) -> Uint<BITS, LIMBS> {
        Uint::<BITS, LIMBS>::from_limbs_slice(&self.to_limbs_asc())
    }
}

#[cfg(feature = "bigdecimal")]
impl Convertable<BigDecimal> for f64 {
    fn convertable_to(&self) -> BigDecimal {
        BigDecimal::try_from(*self).unwrap()
    }
}

#[cfg(feature = "bigdecimal")]
impl Convertable<u128> for BigDecimal {
    fn convertable_to(&self) -> u128 {
        self.as_bigint_and_scale()
            .0
            .into_owned()
            .try_into()
            .unwrap()
    }
}

#[cfg(feature = "bigdecimal")]
impl Convertable<i128> for BigDecimal {
    fn convertable_to(&self) -> i128 {
        self.as_bigint_and_scale()
            .0
            .into_owned()
            .try_into()
            .unwrap()
    }
}

#[cfg(feature = "bigdecimal")]
impl Convertable<f64> for BigDecimal {
    fn convertable_to(&self) -> f64 {
        self.to_plain_string().parse().unwrap()
    }
}

impl<T, D> Convertable<Option<T>> for Option<D>
where
    D: Convertable<T>,
{
    fn convertable_to(&self) -> Option<T> {
        self.as_ref().map(|v| v.convertable_to())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(all(feature = "bigdecimal", feature = "alloy"))]
    mod uint_bigdecimal_tests {
        use super::*;
        use alloy_primitives::{U128, U256};
        use bigdecimal::Zero;

        #[test]
        fn test_uint256_to_bigdecimal() {
            // Test zero
            let zero = U256::ZERO;
            let bd_zero: BigDecimal = zero.convertable_to();
            assert_eq!(bd_zero, BigDecimal::zero());

            // Test one
            let one = U256::from(1u64);
            let bd_one: BigDecimal = one.convertable_to();
            assert_eq!(bd_one, BigDecimal::from(1u64));

            // Test large number
            let large = U256::from(123456789u64);
            let bd_large: BigDecimal = large.convertable_to();
            assert_eq!(bd_large, BigDecimal::from(123456789u64));

            // Test maximum U256
            let max = U256::MAX;
            let bd_max: BigDecimal = max.convertable_to();
            assert!(bd_max > BigDecimal::zero());
        }

        #[test]
        fn test_bigdecimal_to_uint256() {
            // Test zero
            let bd_zero = BigDecimal::zero();
            let u_zero: U256 = bd_zero.convertable_to();
            assert_eq!(u_zero, U256::ZERO);

            // Test one
            let bd_one = BigDecimal::from(1u64);
            let u_one: U256 = bd_one.convertable_to();
            assert_eq!(u_one, U256::from(1u64));

            // Test large number
            let bd_large = BigDecimal::from(999999999u64);
            let u_large: U256 = bd_large.convertable_to();
            assert_eq!(u_large, U256::from(999999999u64));
        }

        #[test]
        fn test_uint128_bigdecimal_roundtrip() {
            let original = U128::from(42424242u64);
            let bd: BigDecimal = original.convertable_to();
            let back: U128 = bd.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(all(feature = "bigdecimal", feature = "alloy"))]
    mod signed_bigdecimal_tests {
        use super::*;
        use alloy_primitives::{I128, I256};
        use bigdecimal::Zero;

        #[test]
        fn test_signed256_to_bigdecimal() {
            // Test zero
            let zero = I256::ZERO;
            let bd_zero: BigDecimal = zero.convertable_to();
            assert_eq!(bd_zero, BigDecimal::zero());

            // Test positive
            let pos = I256::try_from(12345i64).unwrap();
            let bd_pos: BigDecimal = pos.convertable_to();
            assert_eq!(bd_pos, BigDecimal::from(12345i64));

            // Test negative
            let neg = I256::try_from(-12345i64).unwrap();
            let bd_neg: BigDecimal = neg.convertable_to();
            assert_eq!(bd_neg, BigDecimal::from(-12345i64));

            // Test max and min
            let max = I256::MAX;
            let bd_max: BigDecimal = max.convertable_to();
            assert!(bd_max > BigDecimal::zero());

            let min = I256::MIN;
            let bd_min: BigDecimal = min.convertable_to();
            assert!(bd_min < BigDecimal::zero());
        }

        #[test]
        fn test_bigdecimal_to_signed256() {
            // Test zero
            let bd_zero = BigDecimal::zero();
            let s_zero: I256 = bd_zero.convertable_to();
            assert_eq!(s_zero, I256::ZERO);

            // Test positive
            let bd_pos = BigDecimal::from(999i64);
            let s_pos: I256 = bd_pos.convertable_to();
            assert_eq!(s_pos, I256::try_from(999i64).unwrap());

            // Test negative
            let bd_neg = BigDecimal::from(-999i64);
            let s_neg: I256 = bd_neg.convertable_to();
            assert_eq!(s_neg, I256::try_from(-999i64).unwrap());
        }

        #[test]
        fn test_signed128_bigdecimal_roundtrip() {
            let original = I128::try_from(-123456i64).unwrap();
            let bd: BigDecimal = original.convertable_to();
            let back: I128 = bd.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(all(feature = "malachite", feature = "alloy"))]
    mod signed_integer_tests {
        use super::*;
        use alloy_primitives::{I128, I256};
        use malachite::Integer;

        #[test]
        fn test_signed_to_integer() {
            // Test zero
            let zero = I256::ZERO;
            let int_zero: Integer = zero.convertable_to();
            assert_eq!(int_zero, Integer::from(0));

            // Test positive
            let pos = I256::try_from(789i64).unwrap();
            let int_pos: Integer = pos.convertable_to();
            assert_eq!(int_pos, Integer::from(789));

            // Test negative
            let neg = I256::try_from(-789i64).unwrap();
            let int_neg: Integer = neg.convertable_to();
            assert_eq!(int_neg, Integer::from(-789));
        }

        #[test]
        fn test_integer_to_signed() {
            // Test zero
            let int_zero = Integer::from(0);
            let s_zero: I256 = int_zero.convertable_to();
            assert_eq!(s_zero, I256::ZERO);

            // Test positive
            let int_pos = Integer::from(12345);
            let s_pos: I256 = int_pos.convertable_to();
            assert_eq!(s_pos, I256::try_from(12345i64).unwrap());

            // Test negative
            let int_neg = Integer::from(-12345);
            let s_neg: I256 = int_neg.convertable_to();
            assert_eq!(s_neg, I256::try_from(-12345i64).unwrap());
        }

        #[test]
        fn test_signed128_integer_roundtrip() {
            let original = I128::try_from(9876i64).unwrap();
            let int: Integer = original.convertable_to();
            let back: I128 = int.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(all(feature = "bigdecimal", feature = "alloy"))]
    mod biguint_uint_tests {
        use super::*;
        use alloy_primitives::{U256, U512};
        use bigdecimal::Zero;
        use bigdecimal::num_bigint::BigUint;

        #[test]
        fn test_biguint_to_uint() {
            // Test zero
            let bu_zero = BigUint::zero();
            let u_zero: U256 = bu_zero.convertable_to();
            assert_eq!(u_zero, U256::ZERO);

            // Test one
            let bu_one = BigUint::from(1u32);
            let u_one: U256 = bu_one.convertable_to();
            assert_eq!(u_one, U256::from(1u64));

            // Test large number
            let bu_large = BigUint::from(1000000u64);
            let u_large: U256 = bu_large.convertable_to();
            assert_eq!(u_large, U256::from(1000000u64));
        }

        #[test]
        fn test_uint_to_biguint() {
            // Test zero
            let u_zero = U256::ZERO;
            let bu_zero: BigUint = u_zero.convertable_to();
            assert_eq!(bu_zero, BigUint::zero());

            // Test one
            let u_one = U256::from(1u64);
            let bu_one: BigUint = u_one.convertable_to();
            assert_eq!(bu_one, BigUint::from(1u32));

            // Test large number
            let u_large = U256::from(9999999u64);
            let bu_large: BigUint = u_large.convertable_to();
            assert_eq!(bu_large, BigUint::from(9999999u64));
        }

        #[test]
        fn test_uint512_biguint_roundtrip() {
            let original = U512::from(424242424242u64);
            let bu: BigUint = original.convertable_to();
            let back: U512 = bu.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(all(feature = "bigdecimal", feature = "alloy"))]
    mod bigint_signed_tests {
        use super::*;
        use alloy_primitives::I256;
        use alloy_primitives::aliases::I512;
        use bigdecimal::Zero;
        use bigdecimal::num_bigint::BigInt;

        #[test]
        fn test_bigint_to_signed() {
            // Test zero
            let bi_zero = BigInt::zero();
            let s_zero: I256 = bi_zero.convertable_to();
            assert_eq!(s_zero, I256::ZERO);

            // Test positive
            let bi_pos = BigInt::from(54321i64);
            let s_pos: I256 = bi_pos.convertable_to();
            assert_eq!(s_pos, I256::try_from(54321i64).unwrap());

            // Test negative
            let bi_neg = BigInt::from(-54321i64);
            let s_neg: I256 = bi_neg.convertable_to();
            assert_eq!(s_neg, I256::try_from(-54321i64).unwrap());
        }

        #[test]
        fn test_signed_to_bigint() {
            // Test zero
            let s_zero = I256::ZERO;
            let bi_zero: BigInt = s_zero.convertable_to();
            assert_eq!(bi_zero, BigInt::zero());

            // Test positive
            let s_pos = I256::try_from(11111i64).unwrap();
            let bi_pos: BigInt = s_pos.convertable_to();
            assert_eq!(bi_pos, BigInt::from(11111i64));

            // Test negative
            let s_neg = I256::try_from(-11111i64).unwrap();
            let bi_neg: BigInt = s_neg.convertable_to();
            assert_eq!(bi_neg, BigInt::from(-11111i64));
        }

        #[test]
        fn test_signed512_bigint_roundtrip() {
            let original = I512::try_from(-8765432198765i64).unwrap();
            let bi: BigInt = original.convertable_to();
            let back: I512 = bi.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(all(feature = "bigdecimal", feature = "malachite"))]
    mod natural_biguint_tests {
        use super::*;
        use bigdecimal::Zero;
        use bigdecimal::num_bigint::BigUint;
        use malachite::Natural;

        #[test]
        fn test_natural_to_biguint() {
            // Test zero
            let nat_zero = Natural::from(0u32);
            let bu_zero: BigUint = nat_zero.convertable_to();
            assert_eq!(bu_zero, BigUint::zero());

            // Test one
            let nat_one = Natural::from(1u32);
            let bu_one: BigUint = nat_one.convertable_to();
            assert_eq!(bu_one, BigUint::from(1u32));

            // Test large number
            let nat_large = Natural::from(123456789u64);
            let bu_large: BigUint = nat_large.convertable_to();
            assert_eq!(bu_large, BigUint::from(123456789u64));
        }

        #[test]
        fn test_biguint_to_natural() {
            // Test zero
            let bu_zero = BigUint::zero();
            let nat_zero: Natural = bu_zero.convertable_to();
            assert_eq!(nat_zero, Natural::from(0u32));

            // Test one
            let bu_one = BigUint::from(1u32);
            let nat_one: Natural = bu_one.convertable_to();
            assert_eq!(nat_one, Natural::from(1u32));

            // Test large number
            let bu_large = BigUint::from(987654321u64);
            let nat_large: Natural = bu_large.convertable_to();
            assert_eq!(nat_large, Natural::from(987654321u64));
        }

        #[test]
        fn test_natural_to_bigint() {
            // Test zero
            let nat_zero = Natural::from(0u32);
            let bi_zero: BigInt = nat_zero.convertable_to();
            assert_eq!(bi_zero, BigInt::zero());

            // Test positive number
            let nat_pos = Natural::from(424242u64);
            let bi_pos: BigInt = nat_pos.convertable_to();
            assert_eq!(bi_pos, BigInt::from(424242i64));
        }
    }

    #[cfg(all(feature = "bigdecimal", feature = "malachite"))]
    mod integer_bigint_tests {
        use super::*;
        use bigdecimal::Zero;
        use bigdecimal::num_bigint::BigInt;
        use malachite::Integer;

        #[test]
        fn test_integer_to_bigint() {
            // Test zero
            let int_zero = Integer::from(0);
            let bi_zero: BigInt = int_zero.convertable_to();
            assert_eq!(bi_zero, BigInt::zero());

            // Test positive
            let int_pos = Integer::from(777777);
            let bi_pos: BigInt = int_pos.convertable_to();
            assert_eq!(bi_pos, BigInt::from(777777i64));

            // Test negative
            let int_neg = Integer::from(-777777);
            let bi_neg: BigInt = int_neg.convertable_to();
            assert_eq!(bi_neg, BigInt::from(-777777i64));
        }

        #[test]
        fn test_bigint_to_integer() {
            // Test zero
            let bi_zero = BigInt::zero();
            let int_zero: Integer = bi_zero.convertable_to();
            assert_eq!(int_zero, Integer::from(0));

            // Test positive
            let bi_pos = BigInt::from(555555i64);
            let int_pos: Integer = bi_pos.convertable_to();
            assert_eq!(int_pos, Integer::from(555555));

            // Test negative
            let bi_neg = BigInt::from(-555555i64);
            let int_neg: Integer = bi_neg.convertable_to();
            assert_eq!(int_neg, Integer::from(-555555));
        }

        #[test]
        fn test_integer_bigint_roundtrip() {
            let original = Integer::from(-123123123);
            let bi: BigInt = original.convertable_to();
            let back: Integer = bi.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(all(feature = "malachite", feature = "alloy"))]
    mod uint_natural_rational_tests {
        use super::*;
        use alloy_primitives::{U128, U256};
        use malachite::{Natural, Rational};

        #[test]
        fn test_uint_to_natural() {
            // Test zero
            let u_zero = U256::ZERO;
            let nat_zero: Natural = u_zero.convertable_to();
            assert_eq!(nat_zero, Natural::from(0u32));

            // Test one
            let u_one = U256::from(1u64);
            let nat_one: Natural = u_one.convertable_to();
            assert_eq!(nat_one, Natural::from(1u32));

            // Test large number
            let u_large = U256::from(1234567890u64);
            let nat_large: Natural = u_large.convertable_to();
            assert_eq!(nat_large, Natural::from(1234567890u64));
        }

        #[test]
        fn test_natural_to_uint() {
            // Test zero
            let nat_zero = Natural::from(0u32);
            let u_zero: U256 = nat_zero.convertable_to();
            assert_eq!(u_zero, U256::ZERO);

            // Test one
            let nat_one = Natural::from(1u32);
            let u_one: U256 = nat_one.convertable_to();
            assert_eq!(u_one, U256::from(1u64));

            // Test large number
            let nat_large = Natural::from(9876543210u64);
            let u_large: U256 = nat_large.convertable_to();
            assert_eq!(u_large, U256::from(9876543210u64));
        }

        #[test]
        fn test_uint_to_rational() {
            // Test zero
            let u_zero = U256::ZERO;
            let rat_zero: Rational = u_zero.convertable_to();
            assert_eq!(rat_zero, Rational::from(0u32));

            // Test one
            let u_one = U256::from(1u64);
            let rat_one: Rational = u_one.convertable_to();
            assert_eq!(rat_one, Rational::from(1u32));

            // Test large number (should be integer rational)
            let u_large = U256::from(424242u64);
            let rat_large: Rational = u_large.convertable_to();
            assert_eq!(rat_large, Rational::from(424242u64));
        }

        #[test]
        fn test_uint128_natural_roundtrip() {
            let original = U128::from(888888888u64);
            let nat: Natural = original.convertable_to();
            let back: U128 = nat.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(feature = "bigdecimal")]
    mod f64_bigdecimal_tests {
        use super::*;
        use bigdecimal::BigDecimal;

        #[test]
        fn test_f64_to_bigdecimal() {
            // Test zero
            let f_zero = 0.0f64;
            let bd_zero: BigDecimal = f_zero.convertable_to();
            assert_eq!(bd_zero, BigDecimal::from(0));

            // Test positive integer
            let f_pos = 123.0f64;
            let bd_pos: BigDecimal = f_pos.convertable_to();
            assert_eq!(bd_pos, BigDecimal::from(123));

            // Test negative integer
            let f_neg = -123.0f64;
            let bd_neg: BigDecimal = f_neg.convertable_to();
            assert_eq!(bd_neg, BigDecimal::from(-123));

            // Test decimal - f64 has limited precision, so we need to check approximately
            let f_dec = 123.456f64;
            let bd_dec: BigDecimal = f_dec.convertable_to();
            // Convertable back to f64 and compare with epsilon
            let back_to_f64: f64 = bd_dec.to_string().parse().unwrap();
            assert!((back_to_f64 - f_dec).abs() < 1e-10);
        }

        #[test]
        fn test_bigdecimal_to_f64() {
            // Test zero
            let bd_zero = BigDecimal::from(0);
            let f_zero: f64 = bd_zero.convertable_to();
            assert_eq!(f_zero, 0.0f64);

            // Test positive
            let bd_pos = BigDecimal::from(789);
            let f_pos: f64 = bd_pos.convertable_to();
            assert_eq!(f_pos, 789.0f64);

            // Test negative
            let bd_neg = BigDecimal::from(-789);
            let f_neg: f64 = bd_neg.convertable_to();
            assert_eq!(f_neg, -789.0f64);
        }

        #[test]
        fn test_f64_bigdecimal_roundtrip() {
            // Test simple integer
            let original = 42.0f64;
            let bd: BigDecimal = original.convertable_to();
            let back: f64 = bd.convertable_to();
            assert_eq!(original, back);
        }
    }

    #[cfg(feature = "bigdecimal")]
    mod bigdecimal_integer_tests {
        use super::*;
        use bigdecimal::BigDecimal;

        #[test]
        fn test_bigdecimal_to_u128() {
            // Test zero
            let bd_zero = BigDecimal::from(0);
            let u_zero: u128 = bd_zero.convertable_to();
            assert_eq!(u_zero, 0u128);

            // Test positive
            let bd_pos = BigDecimal::from(12345u64);
            let u_pos: u128 = bd_pos.convertable_to();
            assert_eq!(u_pos, 12345u128);

            // Test large number
            let bd_large = BigDecimal::from(999999999999u64);
            let u_large: u128 = bd_large.convertable_to();
            assert_eq!(u_large, 999999999999u128);
        }

        #[test]
        fn test_bigdecimal_to_i128() {
            // Test zero
            let bd_zero = BigDecimal::from(0);
            let i_zero: i128 = bd_zero.convertable_to();
            assert_eq!(i_zero, 0i128);

            // Test positive
            let bd_pos = BigDecimal::from(54321i64);
            let i_pos: i128 = bd_pos.convertable_to();
            assert_eq!(i_pos, 54321i128);

            // Test negative
            let bd_neg = BigDecimal::from(-54321i64);
            let i_neg: i128 = bd_neg.convertable_to();
            assert_eq!(i_neg, -54321i128);
        }
    }

    #[test]
    fn test_option_convert() {
        // Test Some case with simple conversion
        #[derive(Debug, PartialEq)]
        struct FromType(i32);

        #[derive(Debug, PartialEq)]
        struct ToType(i32);

        impl Convertable<ToType> for FromType {
            fn convertable_to(&self) -> ToType {
                ToType(self.0)
            }
        }

        // Test Some(value)
        let some_from = Some(FromType(42));
        let some_to: Option<ToType> = some_from.convertable_to();
        assert_eq!(some_to, Some(ToType(42)));

        // Test None
        let none_from: Option<FromType> = None;
        let none_to: Option<ToType> = none_from.convertable_to();
        assert_eq!(none_to, None);
    }
}
