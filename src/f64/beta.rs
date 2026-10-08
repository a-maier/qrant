use std::f64::consts::PI;

use super::{
    colour::{C_A, C_F, D_44_AA, D_44_AF, D_44_FF, N_A, T_F},
    constants::{ZETA_3, ZETA_4, ZETA_5},
};

use const_for::const_for;

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

// taken from arXiv:1701.01404v2
pub(crate) const BETA_40: f64 = (C_A
    * C_A
    * C_A
    * C_A
    * C_A
    * (8296235. / 3888. - 1630. / 81. * ZETA_3 + 121. / 6. * ZETA_4
        - 1045. / 9. * ZETA_5)
    + D_44_AA / N_A
        * C_A
        * (-514. / 3. + 18716. / 3. * ZETA_3
            - 968. * ZETA_4
            - 15400. / 3. * ZETA_5))
    / 1024.;

pub(crate) const BETA_41: f64 = (T_F
    * C_A
    * C_A
    * C_A
    * C_A
    * (-5048959. / 972. + 10505. / 81. * ZETA_3 - 583. / 3. * ZETA_4
        + 1230. * ZETA_5)
    + T_F
        * C_F
        * C_A
        * C_A
        * C_A
        * (8141995. / 1944. + 146. * ZETA_3 + 902. / 3. * ZETA_4
            - 8720. / 3. * ZETA_5)
    + T_F
        * C_F
        * C_F
        * C_A
        * C_A
        * (-548732. / 81. - 50581. / 27. * ZETA_3 - 484. / 3. * ZETA_4
            + 12820. / 3. * ZETA_5)
    + T_F
        * C_F
        * C_F
        * C_F
        * C_A
        * (3717. + 5696. / 3. * ZETA_3 - 7480. / 3. * ZETA_5)
    + T_F * C_F * C_F * C_F * C_F * (-4157. / 6. - 128. * ZETA_3)
    + C_A * D_44_AF / N_A
        * (11312. / 9. - 127736. / 9. * ZETA_3
            + 2288. * ZETA_4
            + 67520. / 9. * ZETA_5)
    + C_F * D_44_AF / N_A
        * (-320. + 1280. / 3. * ZETA_3 + 6400. / 3. * ZETA_5)
    + T_F * D_44_AA / N_A
        * (904. / 9. - 20752. / 9. * ZETA_3
            + 352. * ZETA_4
            + 4000. / 9. * ZETA_5))
    / 1024.;

pub(crate) const BETA_42: f64 = (T_F
    * T_F
    * C_A
    * C_A
    * C_A
    * (843067. / 486. + 18446. / 27. * ZETA_3
        - 104. / 3. * ZETA_4
        - 2200. / 3. * ZETA_5)
    + T_F
        * T_F
        * C_F
        * C_A
        * C_A
        * (5701. / 162. + 26452. / 27. * ZETA_3 - 944. / 3. * ZETA_4
            + 1600. / 3. * ZETA_5)
    + T_F
        * T_F
        * C_F
        * C_F
        * C_A
        * (31583. / 18. - 28628. / 27. * ZETA_3 + 1144. / 3. * ZETA_4
            - 4400. / 3. * ZETA_5)
    + T_F
        * T_F
        * C_F
        * C_F
        * C_F
        * (-5018. / 9. - 2144. / 3. * ZETA_3 + 4640. / 3. * ZETA_5)
    + C_A * D_44_FF / N_A
        * (-7184. / 3. + 40336. / 9. * ZETA_3 - 704. * ZETA_4
            + 2240. / 9. * ZETA_5)
    + C_F * D_44_FF / N_A
        * (4160. / 3. + 5120. / 3. * ZETA_3 - 12800. / 3. * ZETA_5)
    + T_F * D_44_AF / N_A
        * (-3680. / 9. + 40160. / 9. * ZETA_3
            - 832. * ZETA_4
            - 1280. / 9. * ZETA_5))
    / 1024.;

pub(crate) const BETA_43: f64 = (T_F * D_44_FF / N_A
    * (3520. / 9. - 2624. / 3. * ZETA_3 + 256. * ZETA_4 + 1280. / 3. * ZETA_5)
    + T_F
        * T_F
        * T_F
        * C_A
        * C_A
        * (-2077. / 27. - 9736. / 81. * ZETA_3
            + 112. / 3. * ZETA_4
            + 320. / 9. * ZETA_5)
    + T_F
        * T_F
        * T_F
        * C_F
        * C_A
        * (-736. / 81. - 5680. / 27. * ZETA_3 + 224. / 3. * ZETA_4)
    + T_F
        * T_F
        * T_F
        * C_F
        * C_F
        * (-9922. / 81. + 7616. / 27. * ZETA_3 - 352. / 3. * ZETA_4))
    / 1024.;

pub(crate) const BETA_44: f64 = T_F
    * T_F
    * T_F
    * T_F
    * (C_A * (916. / 243. - 640. / 81. * ZETA_3)
        + C_F * (-856. / 243. - 128. / 27. * ZETA_3))
    / 1024.;

/// Loop order to which the QCD β function is implemented.
pub const MAX_LOOP_ORDER: usize = 5;

/// Coefficients of the QCD β function
///
/// `BETA_COEFF[as_pow][nf_pow]` is *minus* the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow`
pub const BETA_COEFF: [[f64; MAX_LOOP_ORDER]; MAX_LOOP_ORDER] = [
    [BETA_00, BETA_01, 0., 0., 0.],
    [BETA_10, BETA_11, 0., 0., 0.],
    [BETA_20, BETA_21, BETA_22, 0., 0.],
    [BETA_30, BETA_31, BETA_32, BETA_33, 0.],
    [BETA_40, BETA_41, BETA_42, BETA_43, BETA_44],
];

const fn max_nf_pow(order: usize) -> usize {
    if order == 0 { 1 } else { order }
}

/// Perturbative coefficient of the QCD β function.
///
/// The coefficients are defined by β = -Σ βi (α_s/π)^(i+1).
pub const fn β_i(i: usize, nf: i32) -> f64 {
    let beta_coeff = &BETA_COEFF[i];
    let max_nf_pow = max_nf_pow(i);
    let mut res = beta_coeff[max_nf_pow];
    const_for! { nf_pow in (0..max_nf_pow).rev() => {
        res = res.mul_add(nf as f64, beta_coeff[nf_pow]);
    }}
    res
}

/// The QCD β function to the given loop order.
pub const fn β(order: usize, nf: i32, α_s: f64) -> f64 {
    if order == 0 {
        return 0.;
    }
    let i = order - 1;
    let α_s_π = α_s / PI;
    let mut res = -β_i(i, nf);
    const_for! { i in (0..i).rev() => {
        res = res.mul_add(α_s_π, -β_i(i, nf));
    }}
    α_s_π * res
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn βi_nf3() {
        const NF: i32 = 3;
        const EP: f64 = 1e-13;
        assert_abs_diff_eq!(β_i(0, NF), 2.25);
        assert_abs_diff_eq!(β_i(1, NF), 4.);
        assert_abs_diff_eq!(
            β_i(2, NF),
            10.05989583333333333333333333333333333333,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            β_i(3, NF),
            47.22803957345199692022899769740179243706,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            β_i(4, NF),
            127.3221746284614255273481830954551696777,
            epsilon = EP
        );
    }

    #[test]
    fn βi_nf5() {
        const NF: i32 = 5;
        const EP: f64 = 1e-12;
        assert_abs_diff_eq!(
            β_i(0, NF),
            1.916666666666666666666666666666666666666666666666666666667,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            β_i(1, NF),
            2.416666666666666666666666666666666666666666666666666666667,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            β_i(2, NF),
            2.826678240740740740740740740740740740740740740740740740741,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            β_i(3, NF),
            18.85217315933944027193522194150208642692,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            β_i(4, NF),
            15.10801977348336322393151931464672088623,
            epsilon = EP
        );
    }
}
