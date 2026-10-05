use std::f64::consts::PI;

use bon::bon;
use num_complex::{Complex, Complex64};
use num_traits::MulAddAssign;
use ordered_float::OrderedFloat;
use roots::SimpleConvergency;

use crate::complex::C64;

use super::{
    F64,
    alpha_s::{ALPHA_S_MZ, MZ, OutOfRange, StrongCoupling},
    beta::{MAX_LOOP_ORDER, β_i},
    gamma_m::γ_m_i,
};

/// QCD parameter evolution
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct QcdParameters {
    α_s: StrongCoupling,
    m: F64,
    m_residues: [C64; MAX_LOOP_ORDER - 1],
}

#[bon]
impl QcdParameters {
    #[builder]
    pub fn new(
        /// Initial value of strong coupling constant
        #[builder(into, default = OrderedFloat(ALPHA_S_MZ))]
        α_s: F64,
        /// Initial mass value in GeV
        #[builder(into)]
        m: F64,
        /// Initial scale in GeV
        #[builder(into, default = OrderedFloat(MZ))]
        μ: F64,
        /// Initial number of active quark flavours.
        #[builder(default = 5)]
        nf: i32,
        /// Loop order used for evolution.
        ///
        /// The argument must not be larger than [MAX_LOOP_ORDER].
        order: usize,
        /// Numerical precision for the evolution.
        #[builder(into, default = OrderedFloat(1e-16))]
        precision: F64,
        /// Maximum number of iterations used in the evolution.
        #[builder(default = 30)]
        max_iter: usize,
    ) -> Result<Self, OutOfRange<usize>> {
        let α_s = StrongCoupling::builder()
            .value(α_s)
            .μ(μ)
            .nf(nf)
            .order(order)?
            .precision(precision)
            .max_iter(max_iter)
            .build();
        let m_residues = get_residues(α_s);
        Ok(Self {
            α_s, m, m_residues
        })
    }
}

/// A set of QCD parameter values, i.e. strong coupling and quark mass
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct ParameterValues {
    α_s: F64,
    m: F64,
}

impl QcdParameters {
    pub fn with_coupling(α_s: StrongCoupling, m: impl Into<F64>) -> Self {
        let m_residues = get_residues(α_s);
        Self {
            α_s,
            m: m.into(),
            m_residues,
        }
    }

    /// Return the parameter values at the given scale in GeV.
    pub fn at(&self, μ: f64) -> ParameterValues {
        assert_eq!(MAX_LOOP_ORDER, 5);
        let order = self.α_s.order();
        let m_0 = f64::from(self.m);
        let μ0 = f64::from(self.α_s.μ_0());
        let α_s = self.α_s.at(μ);
        let α_s_0 = self.α_s.at(μ0);
        let nf = self.α_s.nf();
        let g0 = γ_m_i(0, nf) / β_i(0, nf);
        let m = m_0
            * match order {
                0 => 1.,
                1 => (α_s / α_s_0).powf(g0),
                o => {
                    let residues = self
                        .α_s
                        .residues()
                        .into_iter()
                        .take(o - 1)
                        .map(|r| Complex::from(r[0]))
                        .zip(self.m_residues);
                    let exponent: f64 = residues
                        .map(|(α_s_i, r)| {
                            let r: Complex64 = r.into();
                            let log = ((α_s - α_s_i) / (α_s_0 - α_s_i)).ln();
                            (r * log).re
                        })
                        .sum();
                    (α_s / α_s_0).powf(g0) * exponent.exp()
                }
            };
        ParameterValues {
            α_s: α_s.into(),
            m: m.into(),
        }
    }

    /// Evolve the to the given scale in GeV, setting `μ0 = μ`
    ///
    /// All further calls to [evolve](Self::evolve) and [at](Self::at)
    /// will start from the new scale.
    pub fn evolve(mut self, μ: f64) -> QcdParameters {
        let ParameterValues { α_s, m } = self.at(μ);
        self.m = m.into();
        self.α_s.value = α_s.into();
        self.α_s.μ = μ.into();
        self
    }

    /// Get the scale-invariant mass m(m) in GeV
    pub fn m_of_m(&self) -> f64 {
        let mut convergency = SimpleConvergency {
            eps: self.precision(),
            max_iter: self.max_iter(),
        };
        let m = *self.m;
        roots::find_root_secant(
            m,
            self.m(m),
            |m| self.m(m) - m,
            &mut convergency,
        )
        .expect("Failed to find scale-invariant mass")
    }

    /// Return the mass value at the given scale in GeV.
    pub fn m(&self, μ: f64) -> f64 {
        *self.at(μ).m
    }

    /// Return the coupling value at the given scale in GeV.
    pub fn α_s(&self, μ: f64) -> f64 {
        *self.at(μ).α_s
    }

    /// Get the strong coupling
    pub fn coupling(&self) -> StrongCoupling {
        self.α_s
    }

    /// Get the current number of active quark flavours.
    pub fn nf(&self) -> i32 {
        self.α_s.nf()
    }

    /// Get the loop order used for evolution.
    pub fn order(&self) -> usize {
        self.α_s.order()
    }

    /// Reference scale
    pub fn μ_0(&self) -> f64 {
        self.α_s.μ_0()
    }

    /// Mass at reference scale
    pub fn m_0(&self) -> f64 {
        *self.m
    }

    /// Coupling at reference scale
    pub fn α_s_0(&self) -> f64 {
        self.α_s.at_μ_0()
    }

    /// Numeric precision of evolution
    pub fn precision(&self) -> f64 {
        self.α_s.precision()
    }

    /// Maximum number of iterations for finding the parameter values
    pub fn max_iter(&self) -> usize {
        self.α_s.max_iter()
    }
}

fn get_residues(α_s: StrongCoupling) -> [C64; MAX_LOOP_ORDER - 1] {
    let order = α_s.order();
    let b = α_s.b();
    let g = gen_g_arr(α_s);

    let mut res = [C64::from(0.); MAX_LOOP_ORDER - 1];
    let roots = α_s
        .residues()
        .into_iter()
        .take(order.saturating_sub(1))
        .map(|r| Complex::from(r[0]));
    for (n, α_s_i) in roots.enumerate() {
        let mut b0_sum = Complex::from(0.);
        for g in g[1..order].iter().rev() {
            b0_sum.mul_add_assign(α_s_i, g.into());
        }
        b0_sum *= *b[0];

        let mut g0_sum = Complex::from(0.);
        for b in b[1..order].iter().map(|b| (**b).into()).rev() {
            g0_sum.mul_add_assign(α_s_i, b);
        }
        g0_sum *= g[0];

        let mut den = Complex::from(0.);
        for (n, bi) in b.iter().map(|b| **b).enumerate().skip(1).rev() {
            den.mul_add_assign(α_s_i, ((n as f64) * bi).into());
        }
        res[n] = ((b0_sum - g0_sum) / den).into();
    }
    res
}

fn gen_g_arr(α_s: StrongCoupling) -> [f64; MAX_LOOP_ORDER] {
    let nf = α_s.nf();
    let order = α_s.order();
    let mut g = [0.; MAX_LOOP_ORDER];
    let b0_inv = PI / β_i(0, nf);
    let mut pi_pow = 1. / PI;
    for (i, g) in g[..order].iter_mut().enumerate() {
        *g = γ_m_i(i, nf) * b0_inv * pi_pow;
        pi_pow /= PI;
    }
    g
}
#[cfg(test)]
mod tests {
    use crate::{alpha_s::ALPHA_S_NF4, gamma_m::γ_m};

    use super::*;

    fn d(mut f: impl FnMut(f64) -> f64, x: f64, ε: f64) -> f64 {
        (-f(x + 2. * ε) + 8. * f(x + ε) - 8. * f(x - ε) + f(x - 2. * ε))
            / (12. * ε)
    }

    #[test]
    fn test_rge() {
        const MB: f64 = 4.180;
        const EP_HIGH: f64 = 10.;
        const EP_LOW: f64 = 1e-3;
        for order in 1..=5 {
            let param = QcdParameters::builder()
                .order(order)
                .m(MB)
                .μ(MB)
                .build()
                .unwrap();
            let nf = param.nf();

            let μ_high = 100_000.;
            let μdμ_mb_high = μ_high * d(|μ| param.m(μ), μ_high, EP_HIGH);
            let ParameterValues { α_s, m } = param.at(μ_high);
            let rge_rhs_high = 2. * γ_m(order, nf, dbg!(*α_s)) * dbg!(*m);
            assert!((μdμ_mb_high - rge_rhs_high).abs() < 3e-12);

            let param = param.evolve(μ_high);

            let μ_low = 1.;
            let μdμ_mb_low = μ_low * d(|μ| param.m(μ), μ_low, EP_LOW);
            let ParameterValues { α_s, m } = param.at(μ_low);
            let rge_rhs_low = 2. * γ_m(order, nf, *α_s) * *m;
            assert!((μdμ_mb_low - rge_rhs_low).abs() < 3e-12);
        }
    }

    #[test]
    fn test_si() {
        let alpha_s = ALPHA_S_NF4.evolve(3.);
        let param = QcdParameters::with_coupling(alpha_s, 0.986);
        assert_eq!(param.m_of_m(), 1.2722136772876158);
    }
}
