use std::{
    f64::consts::PI,
    fmt::Display,
    ops::{Add, Div, Mul, RangeInclusive, Sub},
};

use bon::Builder;
use hyperdual::Dual;
use num_complex::{Complex, Complex64};
use ordered_float::OrderedFloat;
use roots::{
    FloatType,
    Roots::{One, Two},
    SimpleConvergency, find_roots_cubic, find_roots_quartic,
};

use crate::{beta::βi, complex::C64};

use super::{F64, beta::MAX_LOOP_ORDER, dual_wrapper::D64};

/// Z boson mass in GeV
pub const MZ: f64 = 91.1879;
/// Value of the strong coupling at the scale of the Z boson mass
pub const ALPHA_S_MZ: f64 = 0.1180;

/// The value of the coupling constant at the given scale in QCD with `nf` active flavours
///
/// The number of flavours has to be in the range `3 <= nf <= 6`.
/// If it is different from 5, the corresponding lighter or
/// heavier quarks will be decoupled using on-shell decoupling
/// at the scale of the quark mass. The highest available
/// perturbative order will be used. Currently, this means
/// five-loop evolution and four-loop decoupling.
///
/// # Panics
///
/// Panics if the number of flavours is outside the range `3 <= nf <= 6`.
pub fn α_s(nf: i32, μ: f64) -> f64 {
    todo!()
}

#[derive(Builder, Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct StrongCoupling {
    #[builder(into, default = OrderedFloat(ALPHA_S_MZ))]
    value: F64,
    #[builder(into, default = OrderedFloat(MZ))]
    μ: F64,
    #[builder(default = 5)]
    nf: i32,
    #[builder(default = MAX_LOOP_ORDER, with = |o: usize| -> Result<_, OutOfRange<usize>> {
        if o <= MAX_LOOP_ORDER {
            Ok(o)
        } else {
            Err(OutOfRange{value: o, range: 0..=MAX_LOOP_ORDER})
        }
    })]
    order: usize,
    #[builder(skip = gen_b_arr(nf, order))]
    b: [F64; MAX_LOOP_ORDER],
    #[builder(skip = get_residues(&b[..order]))]
    residues: [[C64; 2]; MAX_LOOP_ORDER - 1],
}

impl Default for StrongCoupling {
    fn default() -> Self {
        Self::builder().build()
    }
}

fn get_roots(c: &[F64]) -> [C64; MAX_LOOP_ORDER - 1] {
    let mut roots = [C64::from(0.); MAX_LOOP_ORDER - 1];
    match c.len() {
        0 | 1 => {}
        2 => roots[0] = (-c[0] / c[1]).into(),
        3 => {
            let r = roots2([*c[0], *c[1], *c[2]]);
            roots[0] = r[0].into();
            roots[1] = r[1].into();
        }
        4 => {
            let mut c = [c[0], c[1], c[2], c[3]].map(f64::from);
            let One([x0]) = find_roots_cubic(c[3], c[2], c[1], c[0]) else {
                panic!(
                    "Unexpected number of real roots in four-loop β function"
                )
            };
            c[2] += x0 * c[3];
            c[1] += x0 * c[2];
            let r = roots2([c[1], c[2], c[3]]);
            roots[0] = x0.into();
            roots[1] = r[0].into();
            roots[2] = r[1].into();
        }
        5 => {
            let mut c = [c[0], c[1], c[2], c[3], c[4]].map(f64::from);
            let r = find_roots_quartic(c[4], c[3], c[2], c[1], c[0]);
            let Two([x0, x1]) = r else {
                panic!(
                    "Unexpected number of real roots in five-loop β function"
                )
            };
            c[3] += x0 * c[4];
            c[2] += x0 * c[3];
            c[3] += x1 * c[4];
            c[2] += x1 * c[3];
            let r = roots2([c[2], c[3], c[4]]);
            roots[0] = x0.into();
            roots[1] = x1.into();
            roots[2] = r[0].into();
            roots[3] = r[1].into();
        }
        _ => unreachable!(),
    }
    dbg!(roots)
}

fn roots2(c: [f64; 3]) -> [Complex64; 2] {
    let disc = Complex::from(c[1] * c[1] - 4. * c[0] * c[2]).sqrt();
    let den = 1. / (2. * c[2]);
    [(disc - c[1]) * den, (-disc - c[1]) * den]
}

fn get_residues(b: &[F64]) -> [[C64; 2]; MAX_LOOP_ORDER - 1] {
    let roots: [Complex<f64>; _] = get_roots(b).map(|c| c.into());

    let mut res = [[C64::from(0.); 2]; MAX_LOOP_ORDER - 1];
    for (n, αsi) in roots.into_iter().enumerate().take(b.len() - 1) {
        let mut b1_num = Complex::from(0.);
        let mut αsi_pow = Complex::from(1.);
        for bi in &b[1..] {
            b1_num += bi.0 * αsi_pow;
            αsi_pow *= αsi;
        }

        let mut b0_num = Complex::from(0.);
        let mut αsi_pow = Complex::from(1.);
        for bi in &b[2..] {
            b0_num += bi.0 * αsi_pow;
            αsi_pow *= αsi;
        }

        let num = b[1].0 * b1_num - b[0].0 * b0_num;

        let mut den = Complex::from(0.);
        let mut αsi_pow = Complex::from(1.);
        for (n, bi) in b.iter().enumerate().skip(1) {
            den += (n as f64) * bi.0 * αsi_pow;
            αsi_pow *= αsi;
        }

        res[n][0] = αsi.into();
        res[n][1] = (num / den).into();
    }
    res
}

fn gen_b_arr(nf: i32, order: usize) -> [F64; MAX_LOOP_ORDER] {
    let mut b = [OrderedFloat(0.); MAX_LOOP_ORDER];
    let mut pi_pow = 1. / PI;
    for (i, b) in b[..order].iter_mut().enumerate() {
        *b = (βi(i, nf) * pi_pow).into();
        pi_pow /= PI;
    }
    dbg!(b)
}

impl StrongCoupling {
    /// Return the value at the given scale
    pub fn at<F>(&self, μ: F) -> F
    where
        F: Add<f64, Output = F>
            + Div<f64, Output = F>
            + Mul<f64, Output = F>
            + From<f64>
            + Div<Output = F>
            + Mul<Output = F>
            + AlphaSEq
            + FloatType
            + Ln,
    {
        assert_eq!(MAX_LOOP_ORDER, 5);
        let αs0 = f64::from(self.value);
        let μ0 = f64::from(self.μ);
        let log_μ2 = (μ / μ0).ln() * 2.;
        let b = self.b.map(f64::from);
        let residues = self.residues.map(|r| r.map(Complex::from));
        match self.order {
            0 => F::from(αs0),
            1 => F::from(αs0) / (log_μ2 * αs0 * b[0] + 1.),
            o => {
                let mut convergency = SimpleConvergency {
                    eps: F::from(1e-15),
                    max_iter: 30,
                };
                roots::find_root_brent(
                    F::from(1e-100),
                    F::from(1.),
                    |αs| αs.αs_eq(αs0, &b[..o], &residues[..o - 1], log_μ2),
                    &mut convergency,
                )
                .expect("Failed to find strong coupling")
            }
        }
    }

    /// Evolve to the given scale
    pub fn evolve(&mut self, μ: f64) -> f64 {
        let αs = self.at(μ);
        self.value = OrderedFloat(αs);
        self.μ = OrderedFloat(μ);
        αs
    }

    /// Evolve to the given scale and remove the quark with mass m from the theory
    pub const fn decouple_down(&mut self, μ: f64) {
        todo!()
    }

    /// Evolve to the given scale and add the quark with mass m to the theory
    pub const fn decouple_up(&mut self, μ: f64) {
        todo!()
    }

    pub fn nf(&self) -> i32 {
        self.nf
    }

    pub fn order(&self) -> usize {
        self.order
    }
}

pub trait AlphaSEq {
    fn αs_eq(
        self,
        αs0: f64,
        b: &[f64],
        residues: &[[Complex<f64>; 2]],
        log_μ2: Self,
    ) -> Self;

    fn residue_sum(self, residues: &[[Complex<f64>; 2]], b: &[f64]) -> Self;
}

pub trait Ln {
    fn ln(self) -> Self;
}

trait LnTimes: Ln {
    fn ln_times(self, r: Complex<f64>) -> Self;
}

macro_rules! impl_alphas_eq {
    ($t:ty) => {
        impl AlphaSEq for $t {
            fn αs_eq(
                self,
                αs0: f64,
                b: &[f64],
                residues: &[[Complex<f64>; 2]],
                log_μ2: Self,
            ) -> Self {
                let αs = self;
                assert!(b.len() > 1);
                (αs - αs0) * b[0]
                    + αs * αs0
                        * (log_μ2 * b[0] * b[0]
                            + b[1] * (αs0 / αs).ln()
                            + αs.residue_sum(residues, b)
                            - αs0.residue_sum(residues, b))
            }

            fn residue_sum(
                self,
                residues: &[[Complex<f64>; 2]],
                b: &[f64],
            ) -> $t {
                let αs = self;
                assert_eq!(b.len(), residues.len() + 1);
                let res: Complex<$t> = residues
                    .iter()
                    .map(|[αsi, r]| LnTimes::ln_times(αs - *αsi, *r))
                    .reduce(|a, b| Complex {
                        re: a.re + b.re,
                        im: a.im + b.im,
                    })
                    .unwrap();
                res.re
            }
        }
    };
}

impl_alphas_eq!(f64);
impl_alphas_eq!(D64);

impl Ln for Complex<f64> {
    fn ln(self) -> Self {
        <Complex<f64>>::ln(self)
    }
}

impl LnTimes for Complex<f64> {
    fn ln_times(self, r: Complex<f64>) -> Self {
        self.ln() * r
    }
}

impl Ln for f64 {
    fn ln(self) -> Self {
        f64::ln(self)
    }
}

impl Sub<Complex<f64>> for D64 {
    type Output = Complex<D64>;

    fn sub(self, rhs: Complex<f64>) -> Self::Output {
        let Complex { re, im } = rhs;
        let [re, im] = [re, im].map(D64::from);
        Complex { re, im }
    }
}

impl Ln for Complex<D64> {
    fn ln(self) -> Self {
        let Complex { re, im } = self;
        let re = re.0;
        let im = im.0;
        let Complex { re, im } = Complex { re, im }.ln();
        let [re, im] = [re, im].map(D64);
        Complex { re, im }
    }
}

impl Ln for D64 {
    fn ln(self) -> Self {
        D64::ln(self)
    }
}

impl LnTimes for Complex<D64> {
    fn ln_times(self, r: Complex<f64>) -> Self {
        let Complex { re, im } = self.ln();
        let [re, im] = [re, im].map(|x| x.0);
        let ln = Complex { re, im };

        let Complex { re, im } = r;
        let [re, im] = [re, im].map(Dual::from);
        let r = Complex { re, im };

        let Complex { re, im } = ln * r;
        let [re, im] = [re, im].map(D64);
        Complex { re, im }
    }
}

/// An out-of-range value
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct OutOfRange<T> {
    /// Value that is out of range
    pub value: T,
    /// Range of allowed values
    pub range: RangeInclusive<T>,
}

impl<T: Display> Display for OutOfRange<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self { value, range } = self;
        write!(
            f,
            "{value} is not inside range {}..={}",
            range.start(),
            range.end()
        )
    }
}

pub trait PowI {
    fn powi(self, n: i32) -> Self;
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct F<T>(T);

macro_rules! impl_f {
    ($t:ty) => {
        impl Mul for F<$t> {
            type Output = Self;

            fn mul(self, rhs: Self) -> Self::Output {
                Self(self.0.mul(rhs.0))
            }
        }

        impl Mul<f64> for F<$t> {
            type Output = Self;

            fn mul(self, rhs: f64) -> Self::Output {
                Self(self.0.mul(rhs))
            }
        }

        impl Mul<F<$t>> for f64 {
            type Output = F<$t>;

            fn mul(self, rhs: F<$t>) -> Self::Output {
                F(rhs.0.mul(self))
            }
        }

        impl Add for F<$t> {
            type Output = Self;

            fn add(self, rhs: Self) -> Self::Output {
                Self(self.0.add(rhs.0))
            }
        }

        impl Add<f64> for F<$t> {
            type Output = Self;

            fn add(self, rhs: f64) -> Self::Output {
                Self(self.0.add(rhs))
            }
        }

        impl Add<F<$t>> for f64 {
            type Output = F<$t>;

            fn add(self, rhs: F<$t>) -> Self::Output {
                F(rhs.0.add(self))
            }
        }

        impl Sub for F<$t> {
            type Output = Self;

            fn sub(self, rhs: Self) -> Self::Output {
                Self(self.0.sub(rhs.0))
            }
        }

        impl Sub<f64> for F<$t> {
            type Output = Self;

            fn sub(self, rhs: f64) -> Self::Output {
                Self(self.0.sub(rhs))
            }
        }

        impl Sub<F<$t>> for f64 {
            type Output = F<$t>;

            fn sub(self, rhs: F<$t>) -> Self::Output {
                F(-rhs.0.sub(self))
            }
        }

        impl Div for F<$t> {
            type Output = Self;

            fn div(self, rhs: F<$t>) -> Self::Output {
                F(self.0.div(rhs.0))
            }
        }
    };
}

impl_f!(f64);
impl_f!(Dual<f64>);

impl PowI for F<f64> {
    fn powi(self, n: i32) -> Self {
        F(f64::powi(self.0, n))
    }
}

impl PowI for F<Dual<f64>> {
    fn powi(self, n: i32) -> Self {
        F(hyperdual::Float::powi(self.0, n))
    }
}

#[cfg(test)]
mod tests {
    use crate::beta::β;

    use super::*;

    fn d(mut f: impl FnMut(f64) -> f64, x: f64, ε: f64) -> f64 {
        (
            -f(x + 2.*ε)
                + 8.*f(x + ε)
                - 8.*f(x - ε)
                + f(x - 2.*ε)
        )/(12. * ε)
    }

    #[test]
    fn test_rge() {
        const EP_HIGH: f64 = 10.;
        const EP_LOW: f64 = 1e-3;
        for order in 1..=5 {
            eprintln!(
                "------------------------- {order} ------------------------"
            );
            let mut αs0 =
                StrongCoupling::builder().order(order).unwrap().build();
            let nf = αs0.nf;

            let μ_high = 100_000.;
            let αs_high = αs0.at(μ_high);
            dbg!(αs_high);
            let μdμ_αs_high = μ_high
                * hyperdual::differentiate(μ_high, |μ| αs0.at(D64(μ)).0);
            dbg!(μdμ_αs_high);
            //let μdμ_αs_high_alt = μ_high * (αs0.at(μ_high + EP) - αs0.at(μ_high - EP))/(2. * EP);
            let μdμ_αs_high_alt = μ_high * d(|μ| αs0.at(μ), μ_high, EP_HIGH);
            dbg!(μdμ_αs_high_alt);
            dbg!(μdμ_αs_high - μdμ_αs_high_alt);
            let rge_rhs_high = 2. * β(order, nf, αs_high) * αs_high;
            dbg!(rge_rhs_high - μdμ_αs_high_alt);
            // assert!((μdμ_αs_high - rge_rhs_high).abs() < 1e-16);

            dbg!(αs0.at(1.));
            αs0.evolve(μ_high);

            let μ_low = 1.;
            let αs_low = αs0.at(μ_low);
            dbg!(αs_low);
            let μdμ_αs_low =
                μ_low * hyperdual::differentiate(μ_low, |μ| αs0.at(D64(μ)).0);
            dbg!(μdμ_αs_low);
            let μdμ_αs_low_alt = μ_low * d(|μ| αs0.at(μ), μ_low, EP_LOW);
            dbg!(μdμ_αs_low_alt);
            dbg!(μdμ_αs_low - μdμ_αs_low_alt);
            let rge_rhs_low = 2. * β(order, nf, αs_low) * αs_low;
            dbg!(rge_rhs_low);
            dbg!(rge_rhs_low - μdμ_αs_low_alt);
            // assert!((μdμ_αs_low - rge_rhs_low).abs() < 1e-16);
        }
    }
}
