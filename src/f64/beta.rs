use std::f64::consts::PI;

use super::colour::{C_A, C_F, D_44_AA, D_44_AF, D_44_FF, N_A, T_F};

pub(crate) const ZETA_3: f64 = 1.202056903159594285399738161511449990765;
pub(crate) const ZETA_4: f64 = PI * PI * PI * PI / 90.;
pub(crate) const ZETA_5: f64 = 1.036927755143369926331365486457034168057;

// beta = - sum beta_ij nf^j as^(i+1)
// as = alpha_s/pi
//taken from Czakon hep-ph/0411261
pub(crate) const BETA_00: f64 = 11. / 12. * C_A;
pub(crate) const BETA_01: f64 = -1. / 3. * T_F;

pub(crate) const BETA_10: f64 = 17. / 24. * C_A * C_A;
pub(crate) const BETA_11: f64 = -(C_F / 4. + 5. / 12. * C_A) * T_F;

pub(crate) const BETA_20: f64 = 2857. / 3456. * C_A * C_A * C_A;
pub(crate) const BETA_21: f64 = -T_F
    * (1415. / 1728. * C_A * C_A + 205. / 576. * C_A * C_F - C_F * C_F / 32.);
pub(crate) const BETA_22: f64 =
    T_F * T_F * (79. / 864. * C_A + 11. / 144. * C_F);

pub(crate) const BETA_30: f64 =
    C_A * C_A * C_A * C_A * (150653. / 124416. - 11. / 576. * ZETA_3)
        + D_44_AA / N_A * (-5. / 144. + 11. / 12. * ZETA_3);

pub(crate) const BETA_31: f64 = T_F
    * (C_A * C_A * C_A * (-39143. / 20736. + 17. / 96. * ZETA_3)
        + C_A * C_A * C_F * (7073. / 62208. - 41. / 144. * ZETA_3)
        + C_A * C_F * C_F * (-1051. / 1728. + 11. / 72. * ZETA_3)
        + C_F * C_F * C_F * 23. / 128.)
    + D_44_AF / N_A * (2. / 9. - 13. / 6. * ZETA_3);

pub(crate) const BETA_32: f64 = T_F
    * T_F
    * (C_A * C_A * (3965. / 10368. + 7. / 72. * ZETA_3)
        + C_A * C_F * (67. / 243. + 7. / 36. * ZETA_3)
        + C_F * C_F * (169. / 864. - 11. / 36. * ZETA_3))
    + D_44_FF / N_A * (-11. / 36. + 2. / 3. * ZETA_3);

pub(crate) const BETA_33: f64 =
    T_F * T_F * T_F * (C_A * 53. / 7776. + C_F * 77. / 3888.);

// taken from Baikov, Chetyrkin, Kühn arXiv:1606.08659
pub(crate) const BETA_40: f64 = (8157455. / 16. + 621885. / 2. * ZETA_3
    - 88209. / 2. * ZETA_4
    - 288090. * ZETA_5)
    / 1024.;

pub(crate) const BETA_41: f64 = (-336460813. / 1944. - 4811164. / 81. * ZETA_3
    + 33935. / 6. * ZETA_4
    + 1358995. / 27. * ZETA_5)
    / 1024.;

pub(crate) const BETA_42: f64 = (25960913. / 1944. + 698531. / 81. * ZETA_3
    - 10526. / 9. * ZETA_4
    - 381760. / 81. * ZETA_5)
    / 1024.;

pub(crate) const BETA_43: f64 = (-630559. / 5832. - 48722. / 243. * ZETA_3
    + 1618. / 27. * ZETA_4
    + 460. / 9. * ZETA_5)
    / 1024.;

pub(crate) const BETA_44: f64 = (1205. / 2916. - 152. / 81. * ZETA_3) / 1024.;

pub const MAX_LOOP_ORDER: usize = 5;

pub(crate) const BETA_COEFFS: [[f64; MAX_LOOP_ORDER]; MAX_LOOP_ORDER] = [
    [BETA_00, BETA_01, 0., 0., 0.],
    [BETA_10, BETA_11, 0., 0., 0.],
    [BETA_20, BETA_21, BETA_22, 0., 0.],
    [BETA_30, BETA_31, BETA_32, BETA_33, 0.],
    [BETA_40, BETA_41, BETA_42, BETA_43, BETA_44],
];

const fn max_nf_pow(order: usize) -> usize {
    if order < 2 { 1 } else { order - 1 }
}

/// Perturbative coefficient of the QCD β function
///
/// The coefficients are defined by β = -Σ βi (αs/π)^(i+1)
pub const fn βi(i: usize, nf: i32) -> f64 {
    let beta_coeff = &BETA_COEFFS[i];
    let mut nf_coeff = max_nf_pow(i);
    let mut res = beta_coeff[nf_coeff];
    while nf_coeff > 0 {
        res *= nf as f64;
        nf_coeff -= 1;
        res += beta_coeff[nf_coeff];
    }
    res
}

/// The QCD β function to the given loop order
pub const fn β(order: usize, nf: i32, αs: f64) -> f64 {
    if order == 0 {
        return 0.;
    }
    let mut i = order - 1;
    let mut res = -βi(i, nf) * (αs / PI);
    while i > 0 {
        i -= 1;
        res -= βi(i, nf);
        res *= αs / PI;
    }
    res
}
