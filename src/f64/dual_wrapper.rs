use std::{
    f64::consts::PI,
    ops::{Add, Div, Mul, Neg, Sub},
};

use hyperdual::{Dual, Float};
use roots::FloatType;

#[derive(Copy, Clone, Debug, Default, PartialEq, PartialOrd)]
pub(crate) struct D64(pub Dual<f64>);

impl D64 {
    pub fn ln(self) -> Self {
        Self(self.0.ln())
    }
}

impl Add for D64 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0.add(rhs.0))
    }
}

impl Sub for D64 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0.sub(rhs.0))
    }
}

impl Mul for D64 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0.mul(rhs.0))
    }
}

impl Div for D64 {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(self.0.div(rhs.0))
    }
}

impl Add<f64> for D64 {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        Self(self.0.add(rhs))
    }
}

impl Sub<f64> for D64 {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self::Output {
        Self(self.0.sub(rhs))
    }
}

impl Mul<f64> for D64 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0.mul(rhs))
    }
}

impl Div<f64> for D64 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self(self.0.div(rhs))
    }
}

impl Add<D64> for f64 {
    type Output = D64;

    fn add(self, rhs: D64) -> Self::Output {
        D64::from(self).add(rhs)
    }
}

impl Sub<D64> for f64 {
    type Output = D64;

    fn sub(self, rhs: D64) -> Self::Output {
        D64::from(self).sub(rhs)
    }
}

impl Mul<D64> for f64 {
    type Output = D64;

    fn mul(self, rhs: D64) -> Self::Output {
        D64::from(self).mul(rhs)
    }
}

impl Div<D64> for f64 {
    type Output = D64;

    fn div(self, rhs: D64) -> Self::Output {
        D64::from(self).div(rhs)
    }
}

impl Neg for D64 {
    type Output = D64;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl FloatType for D64 {
    fn zero() -> Self {
        Self::from(0.)
    }

    fn one() -> Self {
        Self::from(1.)
    }

    fn one_third() -> Self {
        Self::from(1. / 3.)
    }

    fn pi() -> Self {
        Self::from(PI)
    }

    fn two_third_pi() -> Self {
        Self::from(2. / 3. * PI)
    }

    fn sqrt(self) -> Self {
        Self(self.0.sqrt())
    }

    fn atan(self) -> Self {
        Self(self.0.atan())
    }

    fn acos(self) -> Self {
        Self(self.0.acos())
    }

    fn sin(self) -> Self {
        Self(self.0.sin())
    }

    fn cos(self) -> Self {
        Self(self.0.cos())
    }

    fn abs(self) -> Self {
        Self(self.0.abs())
    }

    fn powf(self, n: Self) -> Self {
        Self(self.0.powf(n.0))
    }
}

impl From<f64> for D64 {
    fn from(value: f64) -> Self {
        Self(value.into())
    }
}

impl From<i16> for D64 {
    fn from(value: i16) -> Self {
        Self(f64::from(value).into())
    }
}
