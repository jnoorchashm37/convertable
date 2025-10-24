#[cfg(feature = "alloy")]
use alloy_primitives::{I256, U160, U256};
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
impl Convert<BigDecimal> for U256 {
    fn convert_to(&self) -> BigDecimal {
        BigDecimal::from(BigInt::from_bytes_le(Sign::Plus, &self.to_le_bytes_vec()))
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl Convert<U256> for BigDecimal {
    fn convert_to(&self) -> U256 {
        U256::from_le_slice(&self.as_bigint_and_scale().0.to_bytes_le().1)
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl Convert<BigDecimal> for U160 {
    fn convert_to(&self) -> BigDecimal {
        BigDecimal::from(BigInt::from_bytes_le(Sign::Plus, &self.to_le_bytes_vec()))
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl Convert<U160> for BigDecimal {
    fn convert_to(&self) -> U160 {
        U160::from_le_slice(&self.as_bigint_and_scale().0.to_bytes_le().1)
    }
}

#[cfg(all(feature = "bigdecimal", feature = "alloy"))]
impl Convert<BigDecimal> for I256 {
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
impl Convert<I256> for BigDecimal {
    fn convert_to(&self) -> I256 {
        I256::try_from_le_slice(&self.as_bigint_and_scale().0.to_bytes_le().1).unwrap()
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl Convert<Natural> for U256 {
    fn convert_to(&self) -> Natural {
        Natural::from_limbs_asc(self.as_limbs())
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl Convert<Rational> for U256 {
    fn convert_to(&self) -> Rational {
        Rational::from_naturals(self.convert_to(), Natural::ONE)
    }
}

#[cfg(all(feature = "malachite", feature = "alloy"))]
impl Convert<U256> for Natural {
    fn convert_to(&self) -> U256 {
        U256::from_limbs_slice(&self.to_limbs_asc())
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
