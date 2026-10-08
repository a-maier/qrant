// decoupling coefficients for alpha_s from
// Chetyrkin, Kniehl, Steinhauser hep-ph/9708255
// Chetyrkin, Kuehn, Steinhauser  hep-ph/0004189
// Chetyrkin, Kuehn, Sturm  hep-ph/0512060
// zeta_alpha_s = zeta_g^2
// i.e.
// alpha_s^(nf-1) = zeta_alpha_s * alpha_s^(nf)
// alpha_s^(nf+1) = zeta_alpha_s_inv * alpha_s^(nf)
// lm = log(mu^2/m^2)

use std::f64::consts::{LN_2, PI};

use crate::MassScheme;

use super::constants::{
    A4, A5, LN2_2, LN3_2, LN4_2, LN5_2, PI2, PI4, PI6, ZETA_2, ZETA_3, ZETA_5,
};

use const_for::const_for;
use pastey::paste;

const MAX_LOOP_ORDER: usize = 4;
const MAX_NF_POW: usize = 2;
const MAX_LOG_POW: usize = 4;

pub type ZetaCoeffTable =
    [[[f64; MAX_LOG_POW + 1]; MAX_NF_POW + 1]; MAX_LOOP_ORDER + 1];

const ZERO_ZETA_COEFF_TABLE: ZetaCoeffTable =
    [[[0.; MAX_LOG_POW + 1]; MAX_NF_POW + 1]; MAX_LOOP_ORDER + 1];

const fn max_nf_pow(order: usize) -> usize {
    order.saturating_sub(2)
}

const fn max_log_pow(order: usize, nf_pow: usize) -> usize {
    if nf_pow == 0 { order } else { order - 1 }
}

// coefficients of zeta_alpha_s
// TODO: general colour factors
#[allow(clippy::excessive_precision)]
const T54_3: f64 = -8445.8046390310298;
#[allow(clippy::excessive_precision)]
const T62_2: f64 = -4553.4004372195263;

const DELTA_OS_4_NF0: f64 = -2287952896065241. / 43342154956800.
    + 3657317. / 207360. * ZETA_5
    - 44801124029. / 1567641600. * ZETA_3
    + 151369. / 544320. * ZETA_3 * ZETA_3
    - 121. / 36. * A5
    + 21199949. / 544320. * A4
    - 765521. / 116640. * PI2
    + 1439. / 1296. * PI2 * ZETA_3
    - 389455229. / 1567641600. * PI4
    - 151369. / 30481920. * PI6
    + 983. / 972. * LN_2 * PI2
    - 2057. / 51840. * LN_2 * PI4
    - 18377549. / 13063680. * LN2_2 * PI2
    - 121. / 2592. * LN3_2 * PI2
    + 21199949. / 13063680. * LN4_2
    + 121. / 4320. * LN5_2
    + 82037. / 30965760. * T54_3
    - 151369. / 11612160. * T62_2;

const DELTA_OS_4_NF1: f64 = 318983. / 124416.
    + 115. / 576. * ZETA_5
    + 4778329. / 995328. * ZETA_3
    + 173. / 5184. * A4
    + 583. / 972. * PI2
    - 697709. / 14929920. * PI4
    + 11. / 243. * LN_2 * PI2
    - 1709. / 124416. * LN2_2 * PI2
    + 173. / 124416. * LN4_2;

const DELTA_OS_4_NF2: f64 =
    -140825. / 1492992. - 19. / 1728. * ZETA_3 - 13. / 972. * PI2;

/// Coefficients of ζ_α_s in the on-shell scheme
///
/// `ZETA_ALPHA_OS_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_ALPHA_OS_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[1][0][1] = -1. / 6.;
    res[2][0][0] = -7. / 24.;
    res[2][0][1] = -19. / 24.;
    res[2][0][2] = 1. / 36.;
    res[3][0][0] = -68849. / 124416.
        - 7. / 9. * ZETA_2
        - 2. / 9. * ZETA_2 * LN_2
        - 80507. / 27648. * ZETA_3;
    res[3][0][1] = -4465. / 864.;
    res[3][0][2] = -131. / 576.;
    res[3][0][3] = -1. / 216.;
    res[3][1][0] = 2479. / 31104. + ZETA_2 / 9.;
    res[3][1][1] = 409. / 1728.;
    res[4][0][0] = DELTA_OS_4_NF0;
    res[4][0][1] = -21924511. / 746496.
        - 2703685. / 165888. * ZETA_3
        - 217. / 324. * PI2
        - 31. / 162. * LN_2 * PI2;
    res[4][0][2] = -16921. / 2304.;
    res[4][0][3] = -9013. / 10368.;
    res[4][0][4] = 1. / 1296.;
    res[4][1][0] = DELTA_OS_4_NF1;
    res[4][1][1] = 1117159. / 373248.
        + 132283. / 82944. * ZETA_3
        + 5. / 36. * PI2
        + 1. / 81. * LN_2 * PI2;
    res[4][1][2] = 3577. / 5184.;
    res[4][1][3] = 107. / 1728.;
    res[4][2][0] = DELTA_OS_4_NF2;
    res[4][2][1] = -1679. / 186624. - 1. / 162. * PI2;
    res[4][2][2] = -493. / 20736.;
    res
};

const DELTA_MS_4_NF0: f64 = 224509779902759. / 43342154956800.
    + 4289317. / 207360. * ZETA_5
    - 42159357629. / 1567641600. * ZETA_3
    + 151369. / 544320. * ZETA_3 * ZETA_3
    - 121. / 36. * A5
    + 19667789. / 544320. * A4
    - 452556029. / 1567641600. * PI4
    - 151369. / 30481920. * PI6
    - 2057. / 51840. * LN_2 * PI4
    - 19667789. / 13063680. * LN2_2 * PI2
    - 121. / 2592. * LN3_2 * PI2
    + 19667789. / 13063680. * LN4_2
    + 121. / 4320. * LN5_2
    + 82037. / 30965760. * T54_3
    - 151369. / 11612160. * T62_2;

const DELTA_MS_4_NF1: f64 = -749843. / 373248.
    + 115. / 576. * ZETA_5
    + 3581785. / 995328. * ZETA_3
    + 685. / 5184. * A4
    - 541549. / 14929920. * PI4
    - 685. / 124416. * LN2_2 * PI2
    + 685. / 124416. * LN4_2;

const DELTA_MS_4_NF2: f64 = -271883. / 4478976. + 167. / 5184. * ZETA_3;

/// Coefficients of ζ_α_s in the MSbar scheme
///
/// `ZETA_ALPHA_MS_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_ALPHA_MS_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[1][0][1] = -1. / 6.;
    res[2][0][0] = 11. / 72.;
    res[2][0][1] = -11. / 24.;
    res[2][0][2] = 1. / 36.;
    res[3][0][0] = 575263. / 124416. - 82043. / 27648. * ZETA_3;
    res[3][0][1] = -511. / 288.;
    res[3][0][2] = 23. / 192.;
    res[3][0][3] = -1. / 216.;
    res[3][1][0] = -2633. / 31104.;
    res[3][1][1] = 67. / 576.;
    res[3][1][2] = -1. / 36.;
    res[4][0][0] = DELTA_MS_4_NF0;
    res[4][0][1] = 7639841. / 746496. - 2751301. / 165888. * ZETA_3;
    res[4][0][2] = 5317. / 6912.;
    res[4][0][3] = -1597. / 10368.;
    res[4][0][4] = 1. / 1296.;
    res[4][1][0] = DELTA_MS_4_NF1;
    res[4][1][1] = -137801. / 373248. + 110779. / 82944. * ZETA_3;
    res[4][1][2] = -703. / 5184.;
    res[4][1][3] = -53. / 1728.;
    res[4][2][0] = DELTA_MS_4_NF2;
    res[4][2][1] = 6865. / 186624.;
    res[4][2][2] = -77. / 20736.;
    res[4][2][3] = 1. / 324.;
    res
};

/// Coefficients of ζ_α_s for the scheme-invariant mass m(m)
///
/// `ZETA_ALPHA_SI_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_ALPHA_SI_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[1][0][1] = -1. / 6.;
    res[2][0][0] = 11. / 72.;
    res[2][0][1] = -19. / 24.;
    res[2][0][2] = 1. / 36.;
    res[3][0][0] = 575263. / 124416. - 82043. / 27648. * ZETA_3;
    res[3][0][1] = -131. / 32.;
    res[3][0][2] = -131. / 576.;
    res[3][0][3] = -1. / 216.;
    res[3][1][0] = -2633. / 31104.;
    res[3][1][1] = 281. / 1728.;
    res[4][0][0] = DELTA_MS_4_NF0;
    res[4][0][1] = -2745007. / 746496. - 2751301. / 165888. * ZETA_3;
    res[4][0][2] = -30091. / 6912.;
    res[4][0][3] = -9013. / 10368.;
    res[4][0][4] = 1. / 1296.;
    res[4][1][0] = DELTA_MS_4_NF1;
    res[4][1][1] = 156103. / 373248. + 133819. / 82944. * ZETA_3;
    res[4][1][2] = 59. / 192.;
    res[4][1][3] = 107. / 1728.;
    res[4][2][0] = DELTA_MS_4_NF2;
    res[4][2][1] = 8545. / 186624.;
    res[4][2][2] = -79. / 6912.;
    res
};

/// Coefficients of ζ_α_s^-1 in the on-shell scheme
///
/// `ZETA_ALPHA_OS_INV_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_ALPHA_OS_INV_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[1][0][1] = 1. / 6.;
    res[2][0][0] = 7. / 24.;
    res[2][0][1] = 19. / 24.;
    res[2][0][2] = 1. / 36.;
    res[3][0][0] = 58933. / 124416.
        + 2. / 3. * ZETA_2
        + 2. / 9. * ZETA_2 * LN_2
        + 80507. / 27648. * ZETA_3;
    res[3][0][1] = 8941. / 1728.;
    res[3][0][2] = 511. / 576.;
    res[3][0][3] = 1. / 216.;
    res[3][1][0] = -2479. / 31104. - ZETA_2 / 9.;
    res[3][1][1] = -409. / 1728.;
    res[4][0][0] =
        49. / 192. - DELTA_OS_4_NF2 - DELTA_OS_4_NF1 - DELTA_OS_4_NF0;
    res[4][0][1] = 21084715. / 746496.
        + 2922161. / 165888. * ZETA_3
        + 35. / 54. * PI2
        + 35. / 162. * LN_2 * PI2;
    res[4][0][2] = 47039. / 3456.;
    res[4][0][3] = 14149. / 10368.;
    res[4][0][4] = 1. / 1296.;
    res[4][1][0] = -DELTA_OS_4_NF1 - 2. * DELTA_OS_4_NF2;
    res[4][1][1] = -1140191. / 373248.
        - 132283. / 82944. * ZETA_3
        - 47. / 324. * PI2
        - 1. / 81. * LN_2 * PI2;
    res[4][1][2] = -9115. / 10368.;
    res[4][1][3] = -107. / 1728.;
    res[4][2][0] = -DELTA_OS_4_NF2;
    res[4][2][1] = 1679. / 186624. + 1. / 162. * PI2;
    res[4][2][2] = 493. / 20736.;
    res
};

/// Coefficients of ζ_α_s^-1 in the MSbar scheme
///
/// `ZETA_ALPHA_MS_INV_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_ALPHA_MS_INV_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[1][0][1] = 1. / 6.;
    res[2][0][0] = -11. / 72.;
    res[2][0][1] = 11. / 24.;
    res[2][0][2] = 1. / 36.;
    res[3][0][0] = -564731. / 124416. + 82043. / 27648. * ZETA_3;
    res[3][0][1] = 2645. / 1728.;
    res[3][0][2] = 167. / 576.;
    res[3][0][3] = 1. / 216.;
    res[3][1][0] = 2633. / 31104.;
    res[3][1][1] = -67. / 576.;
    res[3][1][2] = 1. / 36.;
    res[4][0][0] =
        121. / 1728. - DELTA_MS_4_NF2 - DELTA_MS_4_NF1 - DELTA_MS_4_NF0;
    res[4][0][1] = -11093717. / 746496. + 3022001. / 165888. * ZETA_3;
    res[4][0][2] = 1837. / 1152.;
    res[4][0][3] = 2909. / 10368.;
    res[4][0][4] = 1. / 1296.;
    res[4][1][0] = -2. * DELTA_MS_4_NF2 - DELTA_MS_4_NF1;
    res[4][1][1] = 141937. / 373248. - 110779. / 82944. * ZETA_3;
    res[4][1][2] = 277. / 10368.;
    res[4][1][3] = 271. / 5184.;
    res[4][2][0] = -DELTA_MS_4_NF2;
    res[4][2][1] = -6865. / 186624.;
    res[4][2][2] = 77. / 20736.;
    res[4][2][3] = -1. / 324.;
    res
};

/// Coefficients of ζ_α_s^-1 for a scheme-invariant mass m(m)
///
/// `ZETA_ALPHA_SI_INV_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_ALPHA_SI_INV_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[1][0][1] = 1. / 6.;
    res[2][0][0] = -11. / 72.;
    res[2][0][1] = 19. / 24.;
    res[2][0][2] = 1. / 36.;
    res[3][0][0] = -564731. / 124416. + 82043. / 27648. * ZETA_3;
    res[3][0][1] = 2191. / 576.;
    res[3][0][2] = 511. / 576.;
    res[3][0][3] = 1. / 216.;
    res[3][1][0] = 2633. / 31104.;
    res[3][1][1] = -281. / 1728.;
    res[4][0][0] =
        121. / 1728. - DELTA_MS_4_NF2 - DELTA_MS_4_NF1 - DELTA_MS_4_NF0;
    res[4][0][1] = -1531493. / 746496. + 2975921. / 165888. * ZETA_3;
    res[4][0][2] = 33887. / 3456.;
    res[4][0][3] = 14149. / 10368.;
    res[4][0][4] = 1. / 1296.;
    res[4][1][0] = -2. * DELTA_MS_4_NF2 - DELTA_MS_4_NF1;
    res[4][1][1] = -158687. / 373248. - 133819. / 82944. * ZETA_3;
    res[4][1][2] = -515. / 1152.;
    res[4][1][3] = -107. / 1728.;
    res[4][2][0] = -DELTA_MS_4_NF2;
    res[4][2][1] = -8545. / 186624.;
    res[4][2][2] = 79. / 6912.;
    res
};

/// Decoupling constant for strong coupling
///
/// The decoupling constant is defined by α_s^(nf-1) = ζ_α_s * α_s^(nf).
/// That is, it allows to determine the strong coupling in the
/// nf-1 flavour theory from the strong coupling in QCD with
/// nf flavours.
/// `log_μ2_m2 = log(μ^2/m^2)`.
pub const fn ζ_α_s(
    scheme: MassScheme,
    α_s_nf: f64,
    nf: i32,
    log_μ2_m2: f64,
    order: usize,
) -> f64 {
    match scheme {
        MassScheme::OS => ζ_α_s_os(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::MS => ζ_α_s_ms(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::SI => ζ_α_s_si(α_s_nf, nf, log_μ2_m2, order),
    }
}

/// Inverse decoupling constant for strong coupling
///
/// The inverse decoupling constant is defined by
/// α_s^(nf+1) = ζ_α_s^-1 * α_s^(nf).
/// That is, it allows to determine the strong coupling in the
/// nf+1 flavour theory from the strong coupling in QCD with
/// nf flavours.
/// `log_μ2_m2 = log(μ^2/m^2)`.
pub const fn ζ_α_s_inv(
    scheme: MassScheme,
    α_s_nf: f64,
    nf: i32,
    log_μ2_m2: f64,
    order: usize,
) -> f64 {
    match scheme {
        MassScheme::OS => ζ_α_s_os_inv(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::MS => ζ_α_s_ms_inv(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::SI => ζ_α_s_si_inv(α_s_nf, nf, log_μ2_m2, order),
    }
}

macro_rules! impl_ζ_α_s {
    ($scheme:ident) => {
        paste! {
            pub const fn [<ζ_α_s_ $scheme:lower>](
                α_s_nf: f64,
                nf: i32,
                log_μ2_m2: f64,
                order: usize,
            ) -> f64 {
                let α_s_π = α_s_nf / PI;
                let mut res = 0f64;
                const_for! { order in (0..order).rev() => {
                    let zeta_coeff = &[<ZETA_ALPHA_ $scheme _COEFF>][order];
                    let max_nf_pow = max_nf_pow(order);
                    let mut order_res = 0f64;
                    const_for! { nf_pow in (0..max_nf_pow + 1).rev() => {
                        let zeta_coeff = &zeta_coeff[nf_pow];
                        let max_log_pow = max_log_pow(order, nf_pow);
                        let mut nf_res = 0f64;
                        const_for! { log_pow in (0..max_log_pow + 1).rev() => {
                            nf_res = nf_res.mul_add(
                                log_μ2_m2,
                                zeta_coeff[log_pow],
                            );
                        }}
                        order_res = order_res.mul_add(nf as f64, nf_res);
                    }}
                    res = res.mul_add(α_s_π, order_res);
                }}
                res
            }
        }
    };
}

impl_ζ_α_s!(OS);
impl_ζ_α_s!(MS);
impl_ζ_α_s!(SI);
impl_ζ_α_s!(OS_INV);
impl_ζ_α_s!(MS_INV);
impl_ζ_α_s!(SI_INV);
