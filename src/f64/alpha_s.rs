use std::{f64::consts::PI, fmt::Display, ops::RangeInclusive, sync::LazyLock};

use bon::Builder;
use num_complex::{Complex, Complex64};
use ordered_float::OrderedFloat;
use roots::SimpleConvergency;

use crate::{beta::βi, complex::C64, zeta_alpha_s::ζ_α_s};

use super::{
    F64,
    beta::MAX_LOOP_ORDER,
    threshold::{BOTTOM_THRESHOLD, CHARM_THRESHOLD, TOP_THRESHOLD, Threshold},
    zeta_alpha_s::ζ_α_s_inv,
};

/// Z boson mass in GeV
pub const MZ: f64 = 91.1879;
/// Value of the strong coupling at the scale of the Z boson mass
pub const ALPHA_S_MZ: f64 = 0.1180;

static ALPHA_S_NF5: LazyLock<StrongCoupling> =
    LazyLock::new(StrongCoupling::default);

static ALPHA_S_NF4: LazyLock<StrongCoupling> =
    LazyLock::new(|| ALPHA_S_NF5.remove_flavour(BOTTOM_THRESHOLD));
static ALPHA_S_NF3: LazyLock<StrongCoupling> =
    LazyLock::new(|| ALPHA_S_NF4.remove_flavour(CHARM_THRESHOLD));

static ALPHA_S_NF6: LazyLock<StrongCoupling> =
    LazyLock::new(|| ALPHA_S_NF5.add_flavour(TOP_THRESHOLD));

/// The value of the coupling constant at the scale `μ` in QCD with `nf` active flavours
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
    match nf {
        3 => α_s_nf3(μ),
        4 => α_s_nf4(μ),
        5 => α_s_nf5(μ),
        6 => α_s_nf6(μ),
        _ => {
            panic!("Number nf of flavours is {nf}, must be within 3 <= nf <= 6")
        }
    }
}

macro_rules! impl_α_s_nf {
    ($($nf:literal), *) => {
        $(
            pastey::paste!{
                /// The value of the coupling constant at the scale `μ` in QCD with a fixed number of active flavours
                pub fn [<α_s_nf $nf>](μ: f64) -> f64 {
                    [<ALPHA_S_NF $nf>].at(μ)
                }
            }
        )*
    };
}

impl_α_s_nf!(3, 4, 5, 6);

/// The strong coupling
#[derive(Builder, Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct StrongCoupling {
    /// Initial value of the strong coupling.
    #[builder(into, default = OrderedFloat(ALPHA_S_MZ))]
    value: F64,
    /// Scale of the initial strong coupling value in GeV.
    #[builder(into, default = OrderedFloat(MZ))]
    μ: F64,
    /// Initial number of active quark flavours.
    #[builder(default = 5)]
    nf: i32,
    /// Loop order used for coupling evolution.
    ///
    /// The argument must not be larger than [MAX_LOOP_ORDER].
    #[builder(default = MAX_LOOP_ORDER, with = |o: usize| -> Result<_, OutOfRange<usize>> {
        if o <= MAX_LOOP_ORDER {
            Ok(o)
        } else {
            Err(OutOfRange{value: o, range: 0..=MAX_LOOP_ORDER})
        }
    })]
    order: usize,
    /// Numerical precision for the evolution.
    #[builder(into, default = OrderedFloat(1e-16))]
    precision: F64,
    /// Maximum number of iterations used in the evolution.
    #[builder(default = 30)]
    max_iter: usize,
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
            for (lhs, rhs) in roots.iter_mut().zip(r) {
                *lhs = rhs.into()
            }
        }
        4 => {
            let r = roots3([*c[0], *c[1], *c[2], *c[3]]);
            for (lhs, rhs) in roots.iter_mut().zip(r) {
                *lhs = rhs.into()
            }
        }
        5 => {
            let r = roots4([*c[0], *c[1], *c[2], *c[3], *c[4]]);
            for (lhs, rhs) in roots.iter_mut().zip(r) {
                *lhs = rhs.into()
            }
        }
        _ => unreachable!(),
    }
    roots
}

// the `roots` crate only gives us the real-valued roots
// so we have to do re-implement the formulas here

fn roots2(c: [f64; 3]) -> [Complex64; 2] {
    let disc = Complex::from(c[1] * c[1] - 4. * c[0] * c[2]).sqrt();
    let den = 1. / (2. * c[2]);
    [(disc - c[1]) * den, (-disc - c[1]) * den]
}

fn roots3(c: [f64; 4]) -> [Complex64; 3] {
    // from https://en.wikipedia.org/wiki/Cubic_equation
    let [a, b, c, d] = [c[3], c[2], c[1], c[0]];
    let δ0 = b * b - 3. * a * c;
    let δ1 = 2. * b.powi(3) - 9. * a * b * c + 27. * a * a * d;
    if δ0 == 0. && δ1 == 0. {
        let x0 = (b / (-3. * a)).into();
        return [x0; 3];
    }
    let mut c3 = δ1 + Complex::from(δ1 * δ1 - 4. * δ0.powi(3)).sqrt();
    if c3 == Complex64::ZERO {
        c3 = δ1 - Complex::from(δ1 * δ1 - 4. * δ0.powi(3)).sqrt();
        assert_ne!(c3, Complex64::ZERO);
    }
    let mut c = (1. / 2. * c3).cbrt();
    let mut roots = [Complex64::ZERO; 3];
    roots[0] = (b + c + δ0 / c) / (-3. * a);
    let unit_root3 = -1. / 2. * (1. + Complex::from(-3.).sqrt());
    c *= unit_root3;
    roots[1] = (b + c + δ0 / c) / (-3. * a);
    c *= unit_root3;
    roots[2] = (b + c + δ0 / c) / (-3. * a);
    roots
}

fn roots4(c: [f64; 5]) -> [Complex64; 4] {
    // from https://en.wikipedia.org/wiki/Quartic_function
    let [a, b, c, d, e] = [c[4], c[3], c[2], c[1], c[0]];
    let δ0 = c * c - 3. * b * d + 12. * a * e;
    let δ1 =
        2. * c * c * c - 9. * b * c * d + 27. * b * b * e + 27. * a * d * d
            - 72. * a * c * e;
    let p = (a * c - 3. / 8. * b * b) * a.powi(-2);
    let q =
        (1. / 8. * b.powi(3) - 1. / 2. * a * b * c + a * a * d) * a.powi(-3);
    let c3 = δ1 + Complex::from(δ1 * δ1 - 4. * δ0.powi(3)).sqrt();
    let qq = (1. / 2. * c3).cbrt();
    let s = 1. / 2. * (-2. / 3. * p + (qq + δ0 / qq) / (3. * a)).sqrt();
    let sqrt_plus = (-4. * s * s - 2. * p + q / s).sqrt();
    let sqrt_minus = (-4. * s * s - 2. * p - q / s).sqrt();
    let mut roots = [Complex64::ZERO; 4];
    roots[0] = -b / (4. * a) - s + 1. / 2. * sqrt_plus;
    roots[1] = -b / (4. * a) - s - 1. / 2. * sqrt_plus;
    roots[2] = -b / (4. * a) + s + 1. / 2. * sqrt_minus;
    roots[3] = -b / (4. * a) + s - 1. / 2. * sqrt_minus;
    roots
}

fn get_residues(b: &[F64]) -> [[C64; 2]; MAX_LOOP_ORDER - 1] {
    let roots: [Complex<f64>; _] = get_roots(b).map(|c| c.into());

    let mut res = [[C64::from(0.); 2]; MAX_LOOP_ORDER - 1];
    for (n, α_s_i) in roots.into_iter().enumerate().take(b.len() - 1) {
        let mut b1_num = Complex::from(0.);
        let mut α_s_i_pow = Complex::from(1.);
        for bi in &b[1..] {
            b1_num += bi.0 * α_s_i_pow;
            α_s_i_pow *= α_s_i;
        }

        let mut b0_num = Complex::from(0.);
        let mut α_s_i_pow = Complex::from(1.);
        for bi in &b[2..] {
            b0_num += bi.0 * α_s_i_pow;
            α_s_i_pow *= α_s_i;
        }

        let num = b[1].0 * b1_num - b[0].0 * b0_num;

        let mut den = Complex::from(0.);
        let mut α_s_i_pow = Complex::from(1.);
        for (n, bi) in b.iter().enumerate().skip(1) {
            den += (n as f64) * bi.0 * α_s_i_pow;
            α_s_i_pow *= α_s_i;
        }

        res[n][0] = α_s_i.into();
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
    b
}

impl StrongCoupling {
    /// Return the coupling value at the given scale in GeV.
    pub fn at(&self, μ: f64) -> f64 {
        assert_eq!(MAX_LOOP_ORDER, 5);
        let α_s_0 = f64::from(self.value);
        let μ0 = f64::from(self.μ);
        let log_μ2 = (μ / μ0).ln() * 2.;
        let b = self.b.map(f64::from);
        let residues = self.residues.map(|r| r.map(Complex::from));
        match self.order {
            0 => α_s_0,
            1 => α_s_0 / (1. + α_s_0 * b[0] * log_μ2),
            o => {
                let b = &b[..o];
                let residues = &residues[..o - 1];
                let residue_sum_α_s_0 = residue_sum(α_s_0, residues, b);
                let mut convergency = SimpleConvergency {
                    eps: self.precision.into(),
                    max_iter: self.max_iter,
                };
                roots::find_root_secant(
                    1e-100,
                    1.,
                    |α_s| {
                        α_s_eq(
                            α_s,
                            α_s_0,
                            b,
                            residues,
                            log_μ2,
                            residue_sum_α_s_0,
                        )
                    },
                    &mut convergency,
                )
                .expect("Failed to find strong coupling")
            }
        }
    }

    /// Evolve to the given scale in GeV.
    ///
    /// All further calls to [evolve](Self::evolve) and [at](Self::at)
    /// will start from the new scale.
    pub fn evolve(mut self, μ: f64) -> Self {
        let α_s = self.at(μ);
        self.value = OrderedFloat(α_s);
        self.μ = OrderedFloat(μ);
        self
    }

    /// Evolve to the given threshold and remove the quark flavour from the theory.
    pub fn remove_flavour(mut self, to_remove: Threshold) -> Self {
        let Threshold { m, scheme, μ } = to_remove;
        self = self.evolve(*μ);
        let log_μ2_m2 = 2. * (self.μ / m).ln();
        self.value *=
            ζ_α_s(scheme, *self.value, self.nf, log_μ2_m2, self.order);
        self.nf -= 1;
        self.update_b();
        self
    }

    /// Evolve to the given threshold and add the quark flavour from the theory.
    pub fn add_flavour(mut self, to_add: Threshold) -> Self {
        let Threshold { m, scheme, μ } = to_add;
        self = self.evolve(*μ);
        let log_μ2_m2 = 2. * (self.μ / m).ln();
        self.value *=
            ζ_α_s_inv(scheme, *self.value, self.nf, log_μ2_m2, self.order);
        self.nf += 1;
        self.update_b();
        self
    }

    /// Get the current number of active quark flavours.
    pub fn nf(&self) -> i32 {
        self.nf
    }

    /// Get the loop order used for evolution.
    pub fn order(&self) -> usize {
        self.order
    }

    /// Reference scale
    pub fn μ0(&self) -> f64 {
        *self.μ
    }

    /// Value at reference scale
    pub fn at_μ0(&self) -> f64 {
        *self.value
    }

    fn update_b(&mut self) {
        let Self {
            nf,
            order,
            b,
            residues,
            ..
        } = self;
        *b = gen_b_arr(*nf, *order);
        *residues = get_residues(&b[..*order]);
    }
}

fn α_s_eq(
    α_s: f64,
    α_s_0: f64,
    b: &[f64],
    residues: &[[Complex<f64>; 2]],
    log_μ2: f64,
    residue_sum_α_s_0: f64,
) -> f64 {
    assert!(b.len() > 1);
    (α_s - α_s_0) * b[0]
        + α_s
            * α_s_0
            * (log_μ2 * b[0] * b[0]
                + b[1] * (α_s_0 / α_s).ln()
                + residue_sum(α_s, residues, b)
                - residue_sum_α_s_0)
}

fn residue_sum(α_s: f64, residues: &[[Complex<f64>; 2]], b: &[f64]) -> f64 {
    assert_eq!(b.len(), residues.len() + 1);
    let res: Complex<_> = residues
        .iter()
        .map(|[α_s_i, r]| (α_s - *α_s_i).ln() * r)
        .sum();
    res.re
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

#[cfg(test)]
mod tests {
    use crate::{
        MassScheme::{MS, OS, SI},
        beta::β,
    };

    use super::*;

    fn d(mut f: impl FnMut(f64) -> f64, x: f64, ε: f64) -> f64 {
        (-f(x + 2. * ε) + 8. * f(x + ε) - 8. * f(x - ε) + f(x - 2. * ε))
            / (12. * ε)
    }

    #[test]
    fn test_rge() {
        const EP_HIGH: f64 = 10.;
        const EP_LOW: f64 = 1e-3;
        for order in 1..=5 {
            let α_s_0 = StrongCoupling::builder().order(order).unwrap().build();
            let nf = α_s_0.nf;

            let μ_high = 100_000.;
            let α_s_high = α_s_0.at(μ_high);
            let μdμ_α_s_high = μ_high * d(|μ| α_s_0.at(μ), μ_high, EP_HIGH);
            let rge_rhs_high = 2. * β(order, nf, α_s_high) * α_s_high;
            assert!((μdμ_α_s_high - rge_rhs_high).abs() < 2e-13);

            let α_s_0 = α_s_0.evolve(μ_high);

            let μ_low = 1.;
            let α_s_low = α_s_0.at(μ_low);
            let μdμ_α_s_low = μ_low * d(|μ| α_s_0.at(μ), μ_low, EP_LOW);
            let rge_rhs_low = 2. * β(order, nf, α_s_low) * α_s_low;
            assert!((μdμ_α_s_low - rge_rhs_low).abs() < 3e-12);
        }
    }

    // reproduce strong coupling from old C++ QCD library
    const QCD_ALPHA_S_MZ: f64 = 0.1184;
    const QCD_MZ: f64 = 91.18;
    const MB_OS: F64 = OrderedFloat(crate::threshold::MB_OS);

    fn α_s_tst() -> StrongCoupling {
        StrongCoupling::builder()
            .value(QCD_ALPHA_S_MZ)
            .μ(QCD_MZ)
            .build()
    }

    #[test]
    fn test_decouple_os() {
        let mt = OrderedFloat(175.);
        let α_s_nf6 = α_s_tst().add_flavour(Threshold {
            m: mt,
            scheme: OS,
            μ: mt,
        });
        assert_eq!(α_s_nf6.nf, 6);
        assert_eq!(α_s_nf6.at_μ0(), 0.10782666678927451);

        let α_s_nf4 = α_s_tst().remove_flavour(Threshold {
            m: MB_OS,
            scheme: OS,
            μ: MB_OS * 2.,
        });
        assert_eq!(α_s_nf4.nf, 4);
        assert_eq!(α_s_nf4.at(5.), 0.21334916179052535);
    }

    #[test]
    fn test_decouple_1l() {
        let mt = OrderedFloat(175.);
        let α_s_tst = StrongCoupling::builder()
            .value(QCD_ALPHA_S_MZ)
            .μ(QCD_MZ)
            .order(1)
            .unwrap()
            .build();
        let α_s_nf6 = α_s_tst.clone().add_flavour(Threshold {
            m: mt,
            scheme: OS,
            μ: mt,
        });
        assert_eq!(α_s_nf6.nf, 6);
        assert_eq!(α_s_nf6.at_μ0(), 0.10820815172148068);

        let α_s_nf4 = α_s_tst.remove_flavour(Threshold {
            m: MB_OS,
            scheme: OS,
            μ: MB_OS * 2.,
        });
        assert_eq!(α_s_nf4.nf, 4);
        assert_eq!(α_s_nf4.at(5.), 0.20677069643701884);
    }

    #[test]
    fn test_decouple_ms() {
        let mb_10 = 3.610;
        let α_s_nf4 = α_s_tst().remove_flavour(Threshold {
            m: mb_10.into(),
            scheme: MS,
            μ: 10.into(),
        });
        assert_eq!(α_s_nf4.nf, 4);
        assert_eq!(α_s_nf4.at(5.), 0.21339618753543466);
    }

    #[test]
    fn test_decouple_si() {
        let mb_mb = OrderedFloat(4.163);
        let α_s_nf4 = α_s_tst().remove_flavour(Threshold {
            m: mb_mb,
            scheme: SI,
            μ: mb_mb * 2.,
        });
        assert_eq!(α_s_nf4.nf, 4);
        assert_eq!(α_s_nf4.at(5.), 0.21342879414639607);
    }

    #[test]
    fn test_closure() {
        let α_s = StrongCoupling::default();
        let α_s_mz = α_s.evolve(100_000.).evolve(1.).at(MZ);
        assert!((α_s_mz - ALPHA_S_MZ).abs() < 1e-16);

        let α_s = StrongCoupling::default();
        let α_s_mz = α_s
            .add_flavour(TOP_THRESHOLD)
            .remove_flavour(TOP_THRESHOLD)
            .remove_flavour(BOTTOM_THRESHOLD)
            .remove_flavour(CHARM_THRESHOLD)
            .add_flavour(CHARM_THRESHOLD)
            .add_flavour(BOTTOM_THRESHOLD)
            .at(MZ);
        assert!((α_s_mz - ALPHA_S_MZ).abs() < 1e-4);

        let α_s = StrongCoupling::default();
        let α_s_mz = α_s
            .remove_flavour(BOTTOM_THRESHOLD)
            .remove_flavour(CHARM_THRESHOLD)
            .add_flavour(CHARM_THRESHOLD)
            .add_flavour(BOTTOM_THRESHOLD)
            .add_flavour(TOP_THRESHOLD)
            .remove_flavour(TOP_THRESHOLD)
            .at(MZ);
        assert!((α_s_mz - ALPHA_S_MZ).abs() < 1e-4);
    }
}
