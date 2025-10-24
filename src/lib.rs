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
#[cfg(feature = "malachite")]
use malachite::{Natural, Rational, num::basic::traits::One};

pub trait Convert<T> {
    fn convert_to(&self) -> T;
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<BigDecimal> for Uint<BITS, LIMBS> {
    fn convert_to(&self) -> BigDecimal {
        BigDecimal::from(BigInt::from_bytes_le(Sign::Plus, &self.to_le_bytes_vec()))
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<Uint<BITS, LIMBS>> for BigDecimal {
    fn convert_to(&self) -> Uint<BITS, LIMBS> {
        Uint::<BITS, LIMBS>::from_le_slice(&self.as_bigint_and_scale().0.to_bytes_le().1)
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<BigDecimal> for Signed<BITS, LIMBS> {
    fn convert_to(&self) -> BigDecimal {
        let sign = if self.is_negative() {
            Sign::Minus
        } else {
            Sign::Plus
        };

        BigDecimal::from(BigInt::from_bytes_le(sign, &self.to_le_bytes::<32>()))
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<Signed<BITS, LIMBS>> for BigDecimal {
    fn convert_to(&self) -> Signed<BITS, LIMBS> {
        Signed::<BITS, LIMBS>::try_from_le_slice(&self.as_bigint_and_scale().0.to_bytes_le().1)
            .unwrap()
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<Uint<BITS, LIMBS>> for BigUint {
    fn convert_to(&self) -> Uint<BITS, LIMBS> {
        Uint::<BITS, LIMBS>::from_be_slice(&self.to_bytes_be())
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<BigUint> for Uint<BITS, LIMBS> {
    fn convert_to(&self) -> BigUint {
        BigUint::from_bytes_le(&self.to_le_bytes_vec())
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<Signed<BITS, LIMBS>> for BigInt {
    fn convert_to(&self) -> Signed<BITS, LIMBS> {
        Signed::<BITS, LIMBS>::try_from_be_slice(&self.to_bytes_be().1).unwrap()
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<BigInt> for Signed<BITS, LIMBS> {
    fn convert_to(&self) -> BigInt {
        let sign = match (self.is_zero(), self.sign()) {
            (true, _) => bigdecimal::num_bigint::Sign::NoSign,
            (_, alloy_primitives::Sign::Negative) => bigdecimal::num_bigint::Sign::Minus,
            (_, alloy_primitives::Sign::Positive) => bigdecimal::num_bigint::Sign::Plus,
        };
        BigInt::from_bytes_le(sign, &self.to_le_bytes::<32>())
    }
}

#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
impl Convert<Natural> for BigUint {
    fn convert_to(&self) -> Natural {
        let v: alloy_primitives::U256 = self.convert_to();
        v.convert_to()
    }
}

#[cfg(all(feature = "bigdecimal", feature = "malachite"))]
impl Convert<BigUint> for Natural {
    fn convert_to(&self) -> BigUint {
        let v: alloy_primitives::U256 = self.convert_to();
        v.convert_to()
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<Natural> for Uint<BITS, LIMBS> {
    fn convert_to(&self) -> Natural {
        Natural::from_limbs_asc(self.as_limbs())
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<Rational> for Uint<BITS, LIMBS> {
    fn convert_to(&self) -> Rational {
        Rational::from_naturals(self.convert_to(), Natural::ONE)
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl<const BITS: usize, const LIMBS: usize> Convert<Uint<BITS, LIMBS>> for Natural {
    fn convert_to(&self) -> Uint<BITS, LIMBS> {
        Uint::<BITS, LIMBS>::from_limbs_slice(&self.to_limbs_asc())
    }
}

#[cfg(feature = "bigdecimal")]
impl Convert<BigDecimal> for f64 {
    fn convert_to(&self) -> BigDecimal {
        BigDecimal::try_from(*self).unwrap()
    }
}

#[cfg(feature = "bigdecimal")]
impl Convert<u128> for BigDecimal {
    fn convert_to(&self) -> u128 {
        self.as_bigint_and_scale()
            .0
            .into_owned()
            .try_into()
            .unwrap()
    }
}

#[cfg(feature = "bigdecimal")]
impl Convert<i128> for BigDecimal {
    fn convert_to(&self) -> i128 {
        self.as_bigint_and_scale()
            .0
            .into_owned()
            .try_into()
            .unwrap()
    }
}

#[cfg(feature = "bigdecimal")]
impl Convert<f64> for BigDecimal {
    fn convert_to(&self) -> f64 {
        self.to_plain_string().parse().unwrap()
    }
}

impl<T, D> Convert<Option<T>> for Option<D>
where
    D: Convert<T>,
{
    fn convert_to(&self) -> Option<T> {
        self.as_ref().map(|v| v.convert_to())
    }
}
