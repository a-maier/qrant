use std::f64::consts::{LN_2, PI};

use crate::MassScheme;

use super::{
    colour::{C_A, C_F, D_44_FF, N_R, T_F},
    constants::{
        A4, A5, B4, C4, C5, LN2_2, LN3_2, LN4_2, LN5_2, PI2, PI4, PI6, ZETA_2,
        ZETA_3, ZETA_4, ZETA_5, ZETA_6,
    },
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

/// Coefficients of ζ_m in the on-shell scheme
///
/// `ZETA_M_OS_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_M_OS_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;

    res[2][0][0] = 89. / 288. * C_F * T_F;

    res[2][0][1] = -5. / 24. * C_F * T_F;

    res[2][0][2] = 1. / 8. * C_F * T_F;

    res[3][0][0] = -4339. / 7776. * C_F * T_F * T_F
        + 11. / 18. * C_F * T_F * T_F * ZETA_3
        + 16627. / 15552. * C_F * C_A * T_F
        - 2. * C_F * C_A * T_F * A4
        - 629. / 576. * C_F * C_A * T_F * ZETA_3
        - 923. / 576. * C_F * C_F * T_F
        + 4. * C_F * C_F * T_F * A4
        + 57. / 32. * C_F * C_F * T_F * ZETA_3
        + 31. / 1440. * PI4 * C_F * C_A * T_F
        - 11. / 360. * PI4 * C_F * C_F * T_F
        + 1. / 12. * LN2_2 * PI2 * C_F * C_A * T_F
        - 1. / 6. * LN2_2 * PI2 * C_F * C_F * T_F
        - 1. / 12. * LN4_2 * C_F * C_A * T_F
        + 1. / 6. * LN4_2 * C_F * C_F * T_F;

    res[3][0][1] = 35. / 432. * C_F * T_F * T_F + 5. / 64. * C_F * C_A * T_F
        - 3. / 4. * C_F * C_A * T_F * ZETA_3
        - 1. / 64. * C_F * C_F * T_F
        + 3. / 4. * C_F * C_F * T_F * ZETA_3;

    res[3][0][2] = 5. / 72. * C_F * T_F * T_F
        + 29. / 96. * C_F * C_A * T_F
        + 1. / 8. * C_F * C_F * T_F;

    res[3][0][3] = -1. / 36. * C_F * T_F * T_F + 11. / 144. * C_F * C_A * T_F;

    res[3][1][0] =
        1327. / 3888. * C_F * T_F * T_F - 2. / 9. * C_F * T_F * T_F * ZETA_3;

    res[3][1][1] = -53. / 144. * C_F * T_F * T_F;

    res[3][1][3] = -1. / 36. * C_F * T_F * T_F;

    res[4][0][0] = -611. / 384. * D_44_FF / N_R + 40. * D_44_FF / N_R * A4
        - 135. / 32. * D_44_FF / N_R * ZETA_5
        + 973. / 64. * D_44_FF / N_R * ZETA_3
        + 15. / 4. * D_44_FF / N_R * ZETA_3 * ZETA_3
        + 22343981. / 19595520. * C_F * T_F * T_F * T_F
        + 51. / 560. * C_F * T_F * T_F * T_F * ZETA_3
        + 122152133. / 203212800. * C_F * C_A * T_F * T_F
        - 4. * C_F * C_A * T_F * T_F * A5
        + 4979. / 1120. * C_F * C_A * T_F * T_F * A4
        + 641. / 96. * C_F * C_A * T_F * T_F * ZETA_5
        + 33749237. / 11289600. * C_F * C_A * T_F * T_F * ZETA_3
        + 4524863. / 829440. * C_F * C_A * C_A * T_F
        - 173. / 15. * C_F * C_A * C_A * T_F * A5
        - 14539. / 640. * C_F * C_A * C_A * T_F * A4
        + 7551. / 1280. * C_F * C_A * C_A * T_F * ZETA_5
        - 1640279. / 184320. * C_F * C_A * C_A * T_F * ZETA_3
        - 21. / 128. * C_F * C_A * C_A * T_F * ZETA_3 * ZETA_3
        + 321131579. / 67737600. * C_F * C_F * T_F * T_F
        + 8. * C_F * C_F * T_F * T_F * A5
        + 18761. / 1680. * C_F * C_F * T_F * T_F * A4
        - 55. / 4. * C_F * C_F * T_F * T_F * ZETA_5
        + 8847973. / 1881600. * C_F * C_F * T_F * T_F * ZETA_3
        + 632953. / 414720. * C_F * C_F * C_A * T_F
        + 514. / 15. * C_F * C_F * C_A * T_F * A5
        + 11321. / 320. * C_F * C_F * C_A * T_F * A4
        - 77977. / 1920. * C_F * C_F * C_A * T_F * ZETA_5
        + 357209. / 30720. * C_F * C_F * C_A * T_F * ZETA_3
        - 11. / 32. * C_F * C_F * C_A * T_F * ZETA_3 * ZETA_3
        + 9361. / 10240. * C_F * C_F * C_F * T_F
        - 112. / 5. * C_F * C_F * C_F * T_F * A5
        + 3227. / 240. * C_F * C_F * C_F * T_F * A4
        + 10003. / 320. * C_F * C_F * C_F * T_F * ZETA_5
        + 1997. / 480. * C_F * C_F * C_F * T_F * ZETA_3
        + 87. / 32. * C_F * C_F * C_F * T_F * ZETA_3 * ZETA_3
        - 5. / 48. * PI2 * C_F * C_F * T_F * T_F
        + 5. / 144. * PI2 * C_F * C_F * C_A * T_F
        - 25. / 192. * PI2 * C_F * C_F * C_F * T_F
        - 289. / 1152. * PI4 * D_44_FF / N_R
        - 47. / 4320. * PI4 * C_F * T_F * T_F * T_F
        - 8551. / 179200. * PI4 * C_F * C_A * T_F * T_F
        + 759689. / 2764800. * PI4 * C_F * C_A * C_A * T_F
        - 319729. / 2419200. * PI4 * C_F * C_F * T_F * T_F
        - 60439. / 153600. * PI4 * C_F * C_F * C_A * T_F
        - 20897. / 172800. * PI4 * C_F * C_F * C_F * T_F
        - 1. / 252. * PI6 * D_44_FF / N_R
        - 1. / 6048. * PI6 * C_F * C_A * C_A * T_F
        - 85. / 24192. * PI6 * C_F * C_F * C_A * T_F
        + 13. / 12096. * PI6 * C_F * C_F * C_F * T_F
        - 5. / 48. * LN_2 * PI2 * C_F * C_F * C_A * T_F
        + 5. / 24. * LN_2 * PI2 * C_F * C_F * C_F * T_F
        - 49. / 720. * LN_2 * PI4 * C_F * C_A * T_F * T_F
        - 1883. / 10800. * LN_2 * PI4 * C_F * C_A * C_A * T_F
        + 49. / 360. * LN_2 * PI4 * C_F * C_F * T_F * T_F
        + 1321. / 2700. * LN_2 * PI4 * C_F * C_F * C_A * T_F
        - 253. / 900. * LN_2 * PI4 * C_F * C_F * C_F * T_F
        - 5. / 3. * LN2_2 * PI2 * D_44_FF / N_R
        - 4979. / 26880. * LN2_2 * PI2 * C_F * C_A * T_F * T_F
        + 14539. / 15360. * LN2_2 * PI2 * C_F * C_A * C_A * T_F
        - 18761. / 40320. * LN2_2 * PI2 * C_F * C_F * T_F * T_F
        - 11321. / 7680. * LN2_2 * PI2 * C_F * C_F * C_A * T_F
        - 3227. / 5760. * LN2_2 * PI2 * C_F * C_F * C_F * T_F
        - 1. / 18. * LN3_2 * PI2 * C_F * C_A * T_F * T_F
        - 173. / 1080. * LN3_2 * PI2 * C_F * C_A * C_A * T_F
        + 1. / 9. * LN3_2 * PI2 * C_F * C_F * T_F * T_F
        + 257. / 540. * LN3_2 * PI2 * C_F * C_F * C_A * T_F
        - 14. / 45. * LN3_2 * PI2 * C_F * C_F * C_F * T_F
        + 5. / 3. * LN4_2 * D_44_FF / N_R
        + 4979. / 26880. * LN4_2 * C_F * C_A * T_F * T_F
        - 14539. / 15360. * LN4_2 * C_F * C_A * C_A * T_F
        + 18761. / 40320. * LN4_2 * C_F * C_F * T_F * T_F
        + 11321. / 7680. * LN4_2 * C_F * C_F * C_A * T_F
        + 3227. / 5760. * LN4_2 * C_F * C_F * C_F * T_F
        + 1. / 30. * LN5_2 * C_F * C_A * T_F * T_F
        + 173. / 1800. * LN5_2 * C_F * C_A * C_A * T_F
        - 1. / 15. * LN5_2 * C_F * C_F * T_F * T_F
        - 257. / 900. * LN5_2 * C_F * C_F * C_A * T_F
        + 14. / 75. * LN5_2 * C_F * C_F * C_F * T_F;

    res[4][0][1] = 1. / 4. * D_44_FF / N_R
        - 15. / 8. * D_44_FF / N_R * ZETA_3
        - 83. / 2592. * C_F * T_F * T_F * T_F
        + 1. / 18. * C_F * T_F * T_F * T_F * ZETA_3
        - 107467. / 62208. * C_F * C_A * T_F * T_F
        + 5053. / 4608. * C_F * C_A * T_F * T_F * ZETA_3
        + 233903. / 248832. * C_F * C_A * C_A * T_F
        - 11. / 2. * C_F * C_A * C_A * T_F * A4
        + 25. / 16. * C_F * C_A * C_A * T_F * ZETA_5
        - 119723. / 18432. * C_F * C_A * C_A * T_F * ZETA_3
        - 2365. / 6912. * C_F * C_F * T_F * T_F
        + 181. / 256. * C_F * C_F * T_F * T_F * ZETA_3
        - 3025. / 1152. * C_F * C_F * C_A * T_F
        + 11. * C_F * C_F * C_A * T_F * A4
        + 5. / 16. * C_F * C_F * C_A * T_F * ZETA_5
        + 7569. / 1024. * C_F * C_F * C_A * T_F * ZETA_3
        - 17. / 192. * C_F * C_F * C_F * T_F
        - 15. / 8. * C_F * C_F * C_F * T_F * ZETA_5
        + 553. / 256. * C_F * C_F * C_F * T_F * ZETA_3
        + 1. / 8. * PI2 * C_F * C_F * T_F * T_F
        - 1. / 24. * PI2 * C_F * C_F * C_A * T_F
        + 5. / 32. * PI2 * C_F * C_F * C_F * T_F
        + 1. / 240. * PI4 * C_F * C_A * T_F * T_F
        + 407. / 5760. * PI4 * C_F * C_A * C_A * T_F
        - 1. / 240. * PI4 * C_F * C_F * T_F * T_F
        - 55. / 576. * PI4 * C_F * C_F * C_A * T_F
        + 1. / 8. * LN_2 * PI2 * C_F * C_F * C_A * T_F
        - 1. / 4. * LN_2 * PI2 * C_F * C_F * C_F * T_F
        + 11. / 48. * LN2_2 * PI2 * C_F * C_A * C_A * T_F
        - 11. / 24. * LN2_2 * PI2 * C_F * C_F * C_A * T_F
        - 11. / 48. * LN4_2 * C_F * C_A * C_A * T_F
        + 11. / 24. * LN4_2 * C_F * C_F * C_A * T_F;

    res[4][0][2] = -35. / 864. * C_F * T_F * T_F * T_F
        + 3103. / 6912. * C_F * C_A * T_F * T_F
        + 3. / 8. * C_F * C_A * T_F * T_F * ZETA_3
        + 19867. / 13824. * C_F * C_A * C_A * T_F
        - 33. / 32. * C_F * C_A * C_A * T_F * ZETA_3
        + 613. / 2304. * C_F * C_F * T_F * T_F
        - 3. / 8. * C_F * C_F * T_F * T_F * ZETA_3
        + 91. / 192. * C_F * C_F * C_A * T_F
        + 33. / 32. * C_F * C_F * C_A * T_F * ZETA_3
        + 33. / 64. * C_F * C_F * C_F * T_F;

    res[4][0][3] = -5. / 216. * C_F * T_F * T_F * T_F
        - 13. / 96. * C_F * C_A * T_F * T_F
        + 17. / 48. * C_F * C_A * C_A * T_F
        - 17. / 192. * C_F * C_F * T_F * T_F
        + 55. / 384. * C_F * C_F * C_A * T_F;

    res[4][0][4] = 1. / 144. * C_F * T_F * T_F * T_F
        - 11. / 576. * C_F * C_A * T_F * T_F
        + 121. / 2304. * C_F * C_A * C_A * T_F
        + 1. / 128. * C_F * C_F * T_F * T_F;

    res[4][1][0] = -39959. / 20736. * C_F * T_F * T_F * T_F
        + 5. / 144. * C_F * T_F * T_F * T_F * ZETA_3
        - 5095. / 3072. * C_F * C_A * T_F * T_F
        + 4. * C_F * C_A * T_F * T_F * A5
        + 49. / 12. * C_F * C_A * T_F * T_F * A4
        - 253. / 96. * C_F * C_A * T_F * T_F * ZETA_5
        - 65. / 192. * C_F * C_A * T_F * T_F * ZETA_3
        - 35551. / 10368. * C_F * C_F * T_F * T_F
        - 8. * C_F * C_F * T_F * T_F * A5
        - 49. / 6. * C_F * C_F * T_F * T_F * A4
        + 157. / 16. * C_F * C_F * T_F * T_F * ZETA_5
        - 5. / 8. * C_F * C_F * T_F * T_F * ZETA_3
        + 5. / 144. * PI2 * C_F * C_F * T_F * T_F
        + 23. / 1440. * PI4 * C_F * T_F * T_F * T_F
        - 181. / 3840. * PI4 * C_F * C_A * T_F * T_F
        + 1639. / 17280. * PI4 * C_F * C_F * T_F * T_F
        + 49. / 720. * LN_2 * PI4 * C_F * C_A * T_F * T_F
        - 49. / 360. * LN_2 * PI4 * C_F * C_F * T_F * T_F
        - 49. / 288. * LN2_2 * PI2 * C_F * C_A * T_F * T_F
        + 49. / 144. * LN2_2 * PI2 * C_F * C_F * T_F * T_F
        + 1. / 18. * LN3_2 * PI2 * C_F * C_A * T_F * T_F
        - 1. / 9. * LN3_2 * PI2 * C_F * C_F * T_F * T_F
        + 49. / 288. * LN4_2 * C_F * C_A * T_F * T_F
        - 49. / 144. * LN4_2 * C_F * C_F * T_F * T_F
        - 1. / 30. * LN5_2 * C_F * C_A * T_F * T_F
        + 1. / 15. * LN5_2 * C_F * C_F * T_F * T_F;

    res[4][1][1] = 2543. / 3888. * C_F * T_F * T_F * T_F
        - 7. / 9. * C_F * T_F * T_F * T_F * ZETA_3
        - 5155. / 31104. * C_F * C_A * T_F * T_F
        + 2. * C_F * C_A * T_F * T_F * A4
        + 997. / 576. * C_F * C_A * T_F * T_F * ZETA_3
        + 139. / 144. * C_F * C_F * T_F * T_F
        - 4. * C_F * C_F * T_F * T_F * A4
        - 97. / 32. * C_F * C_F * T_F * T_F * ZETA_3
        - 1. / 24. * PI2 * C_F * C_F * T_F * T_F
        - 43. / 1440. * PI4 * C_F * C_A * T_F * T_F
        + 7. / 180. * PI4 * C_F * C_F * T_F * T_F
        - 1. / 12. * LN2_2 * PI2 * C_F * C_A * T_F * T_F
        + 1. / 6. * LN2_2 * PI2 * C_F * C_F * T_F * T_F
        + 1. / 12. * LN4_2 * C_F * C_A * T_F * T_F
        - 1. / 6. * LN4_2 * C_F * C_F * T_F * T_F;

    res[4][1][2] = 35. / 864. * C_F * T_F * T_F * T_F
        - 2581. / 3456. * C_F * C_A * T_F * T_F
        - 101. / 192. * C_F * C_F * T_F * T_F;

    res[4][1][3] =
        -13. / 96. * C_F * C_A * T_F * T_F - 1. / 16. * C_F * C_F * T_F * T_F;

    res[4][1][4] =
        1. / 144. * C_F * T_F * T_F * T_F - 11. / 288. * C_F * C_A * T_F * T_F;

    res[4][2][0] = 17671. / 20736. * C_F * T_F * T_F * T_F
        - 5. / 144. * C_F * T_F * T_F * T_F * ZETA_3
        - 7. / 1440. * PI4 * C_F * T_F * T_F * T_F;

    res[4][2][1] = -3401. / 7776. * C_F * T_F * T_F * T_F
        + 7. / 18. * C_F * T_F * T_F * T_F * ZETA_3;

    res[4][2][2] = 31. / 216. * C_F * T_F * T_F * T_F;

    res[4][2][4] = 1. / 144. * C_F * T_F * T_F * T_F;
    res
};

/// Coefficients of ζ_m in the MSbar scheme
///
/// `ZETA_M_MS_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_M_MS_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;

    res[2][0][0] = 89. / 288. * C_F * T_F;

    res[2][0][1] = -5. / 24. * C_F * T_F;

    res[2][0][2] = 1. / 8. * C_F * T_F;

    res[3][0][0] = -4339. / 7776. * C_F * T_F * T_F
        + 11. / 18. * C_F * T_F * T_F * ZETA_3
        + 16627. / 15552. * C_F * C_A * T_F
        - 2. * C_F * C_A * T_F * A4
        - 629. / 576. * C_F * C_A * T_F * ZETA_3
        - 683. / 576. * C_F * C_F * T_F
        + 4. * C_F * C_F * T_F * A4
        + 57. / 32. * C_F * C_F * T_F * ZETA_3
        + 31. / 1440. * PI4 * C_F * C_A * T_F
        - 11. / 360. * PI4 * C_F * C_F * T_F
        + 1. / 12. * LN2_2 * PI2 * C_F * C_A * T_F
        - 1. / 6. * LN2_2 * PI2 * C_F * C_F * T_F
        - 1. / 12. * LN4_2 * C_F * C_A * T_F
        + 1. / 6. * LN4_2 * C_F * C_F * T_F;

    res[3][0][1] = 35. / 432. * C_F * T_F * T_F + 5. / 64. * C_F * C_A * T_F
        - 3. / 4. * C_F * C_A * T_F * ZETA_3
        - 13. / 64. * C_F * C_F * T_F
        + 3. / 4. * C_F * C_F * T_F * ZETA_3;

    res[3][0][2] = 5. / 72. * C_F * T_F * T_F + 29. / 96. * C_F * C_A * T_F
        - 1. / 4. * C_F * C_F * T_F;

    res[3][0][3] = -1. / 36. * C_F * T_F * T_F + 11. / 144. * C_F * C_A * T_F;

    res[3][1][0] =
        1327. / 3888. * C_F * T_F * T_F - 2. / 9. * C_F * T_F * T_F * ZETA_3;

    res[3][1][1] = -53. / 144. * C_F * T_F * T_F;

    res[3][1][3] = -1. / 36. * C_F * T_F * T_F;

    res[4][0][0] = -611. / 384. * D_44_FF / N_R + 40. * D_44_FF / N_R * A4
        - 135. / 32. * D_44_FF / N_R * ZETA_5
        + 973. / 64. * D_44_FF / N_R * ZETA_3
        + 15. / 4. * D_44_FF / N_R * ZETA_3 * ZETA_3
        + 22343981. / 19595520. * C_F * T_F * T_F * T_F
        + 51. / 560. * C_F * T_F * T_F * T_F * ZETA_3
        + 122152133. / 203212800. * C_F * C_A * T_F * T_F
        - 4. * C_F * C_A * T_F * T_F * A5
        + 4979. / 1120. * C_F * C_A * T_F * T_F * A4
        + 641. / 96. * C_F * C_A * T_F * T_F * ZETA_5
        + 33749237. / 11289600. * C_F * C_A * T_F * T_F * ZETA_3
        + 4524863. / 829440. * C_F * C_A * C_A * T_F
        - 173. / 15. * C_F * C_A * C_A * T_F * A5
        - 14539. / 640. * C_F * C_A * C_A * T_F * A4
        + 7551. / 1280. * C_F * C_A * C_A * T_F * ZETA_5
        - 1640279. / 184320. * C_F * C_A * C_A * T_F * ZETA_3
        - 21. / 128. * C_F * C_A * C_A * T_F * ZETA_3 * ZETA_3
        + 32109731. / 7526400. * C_F * C_F * T_F * T_F
        + 8. * C_F * C_F * T_F * T_F * A5
        + 18761. / 1680. * C_F * C_F * T_F * T_F * A4
        - 55. / 4. * C_F * C_F * T_F * T_F * ZETA_5
        + 8847973. / 1881600. * C_F * C_F * T_F * T_F * ZETA_3
        + 1068103. / 414720. * C_F * C_F * C_A * T_F
        + 514. / 15. * C_F * C_F * C_A * T_F * A5
        + 11321. / 320. * C_F * C_F * C_A * T_F * A4
        - 77977. / 1920. * C_F * C_F * C_A * T_F * ZETA_5
        + 398489. / 30720. * C_F * C_F * C_A * T_F * ZETA_3
        - 11. / 32. * C_F * C_F * C_A * T_F * ZETA_3 * ZETA_3
        + 10301. / 10240. * C_F * C_F * C_F * T_F
        - 112. / 5. * C_F * C_F * C_F * T_F * A5
        + 3227. / 240. * C_F * C_F * C_F * T_F * A4
        + 10003. / 320. * C_F * C_F * C_F * T_F * ZETA_5
        + 1427. / 480. * C_F * C_F * C_F * T_F * ZETA_3
        + 87. / 32. * C_F * C_F * C_F * T_F * ZETA_3 * ZETA_3
        - 289. / 1152. * PI4 * D_44_FF / N_R
        - 47. / 4320. * PI4 * C_F * T_F * T_F * T_F
        - 8551. / 179200. * PI4 * C_F * C_A * T_F * T_F
        + 759689. / 2764800. * PI4 * C_F * C_A * C_A * T_F
        - 319729. / 2419200. * PI4 * C_F * C_F * T_F * T_F
        - 60439. / 153600. * PI4 * C_F * C_F * C_A * T_F
        - 20897. / 172800. * PI4 * C_F * C_F * C_F * T_F
        - 1. / 252. * PI6 * D_44_FF / N_R
        - 1. / 6048. * PI6 * C_F * C_A * C_A * T_F
        - 85. / 24192. * PI6 * C_F * C_F * C_A * T_F
        + 13. / 12096. * PI6 * C_F * C_F * C_F * T_F
        - 49. / 720. * LN_2 * PI4 * C_F * C_A * T_F * T_F
        - 1883. / 10800. * LN_2 * PI4 * C_F * C_A * C_A * T_F
        + 49. / 360. * LN_2 * PI4 * C_F * C_F * T_F * T_F
        + 1321. / 2700. * LN_2 * PI4 * C_F * C_F * C_A * T_F
        - 253. / 900. * LN_2 * PI4 * C_F * C_F * C_F * T_F
        - 5. / 3. * LN2_2 * PI2 * D_44_FF / N_R
        - 4979. / 26880. * LN2_2 * PI2 * C_F * C_A * T_F * T_F
        + 14539. / 15360. * LN2_2 * PI2 * C_F * C_A * C_A * T_F
        - 18761. / 40320. * LN2_2 * PI2 * C_F * C_F * T_F * T_F
        - 11321. / 7680. * LN2_2 * PI2 * C_F * C_F * C_A * T_F
        - 3227. / 5760. * LN2_2 * PI2 * C_F * C_F * C_F * T_F
        - 1. / 18. * LN3_2 * PI2 * C_F * C_A * T_F * T_F
        - 173. / 1080. * LN3_2 * PI2 * C_F * C_A * C_A * T_F
        + 1. / 9. * LN3_2 * PI2 * C_F * C_F * T_F * T_F
        + 257. / 540. * LN3_2 * PI2 * C_F * C_F * C_A * T_F
        - 14. / 45. * LN3_2 * PI2 * C_F * C_F * C_F * T_F
        + 5. / 3. * LN4_2 * D_44_FF / N_R
        + 4979. / 26880. * LN4_2 * C_F * C_A * T_F * T_F
        - 14539. / 15360. * LN4_2 * C_F * C_A * C_A * T_F
        + 18761. / 40320. * LN4_2 * C_F * C_F * T_F * T_F
        + 11321. / 7680. * LN4_2 * C_F * C_F * C_A * T_F
        + 3227. / 5760. * LN4_2 * C_F * C_F * C_F * T_F
        + 1. / 30. * LN5_2 * C_F * C_A * T_F * T_F
        + 173. / 1800. * LN5_2 * C_F * C_A * C_A * T_F
        - 1. / 15. * LN5_2 * C_F * C_F * T_F * T_F
        - 257. / 900. * LN5_2 * C_F * C_F * C_A * T_F
        + 14. / 75. * LN5_2 * C_F * C_F * C_F * T_F;

    res[4][0][1] = 1. / 4. * D_44_FF / N_R
        - 15. / 8. * D_44_FF / N_R * ZETA_3
        - 83. / 2592. * C_F * T_F * T_F * T_F
        + 1. / 18. * C_F * T_F * T_F * T_F * ZETA_3
        - 107467. / 62208. * C_F * C_A * T_F * T_F
        + 5053. / 4608. * C_F * C_A * T_F * T_F * ZETA_3
        + 233903. / 248832. * C_F * C_A * C_A * T_F
        - 11. / 2. * C_F * C_A * C_A * T_F * A4
        + 25. / 16. * C_F * C_A * C_A * T_F * ZETA_5
        - 119723. / 18432. * C_F * C_A * C_A * T_F * ZETA_3
        - 2533. / 6912. * C_F * C_F * T_F * T_F
        + 181. / 256. * C_F * C_F * T_F * T_F * ZETA_3
        - 3529. / 768. * C_F * C_F * C_A * T_F
        + 11. * C_F * C_F * C_A * T_F * A4
        + 5. / 16. * C_F * C_F * C_A * T_F * ZETA_5
        + 8913. / 1024. * C_F * C_F * C_A * T_F * ZETA_3
        + 217. / 768. * C_F * C_F * C_F * T_F
        - 15. / 8. * C_F * C_F * C_F * T_F * ZETA_5
        + 169. / 256. * C_F * C_F * C_F * T_F * ZETA_3
        + 1. / 240. * PI4 * C_F * C_A * T_F * T_F
        + 407. / 5760. * PI4 * C_F * C_A * C_A * T_F
        - 1. / 240. * PI4 * C_F * C_F * T_F * T_F
        - 55. / 576. * PI4 * C_F * C_F * C_A * T_F
        + 11. / 48. * LN2_2 * PI2 * C_F * C_A * C_A * T_F
        - 11. / 24. * LN2_2 * PI2 * C_F * C_F * C_A * T_F
        - 11. / 48. * LN4_2 * C_F * C_A * C_A * T_F
        + 11. / 24. * LN4_2 * C_F * C_F * C_A * T_F;

    res[4][0][2] = -35. / 864. * C_F * T_F * T_F * T_F
        + 3103. / 6912. * C_F * C_A * T_F * T_F
        + 3. / 8. * C_F * C_A * T_F * T_F * ZETA_3
        + 19867. / 13824. * C_F * C_A * C_A * T_F
        - 33. / 32. * C_F * C_A * C_A * T_F * ZETA_3
        + 517. / 2304. * C_F * C_F * T_F * T_F
        - 3. / 8. * C_F * C_F * T_F * T_F * ZETA_3
        - 219. / 128. * C_F * C_F * C_A * T_F
        + 33. / 32. * C_F * C_F * C_A * T_F * ZETA_3
        + 15. / 16. * C_F * C_F * C_F * T_F;

    res[4][0][3] = -5. / 216. * C_F * T_F * T_F * T_F
        - 13. / 96. * C_F * C_A * T_F * T_F
        + 17. / 48. * C_F * C_A * C_A * T_F
        + 7. / 192. * C_F * C_F * T_F * T_F
        - 143. / 384. * C_F * C_F * C_A * T_F;

    res[4][0][4] = 1. / 144. * C_F * T_F * T_F * T_F
        - 11. / 576. * C_F * C_A * T_F * T_F
        + 121. / 2304. * C_F * C_A * C_A * T_F
        + 1. / 128. * C_F * C_F * T_F * T_F;

    res[4][1][0] = -39959. / 20736. * C_F * T_F * T_F * T_F
        + 5. / 144. * C_F * T_F * T_F * T_F * ZETA_3
        - 5095. / 3072. * C_F * C_A * T_F * T_F
        + 4. * C_F * C_A * T_F * T_F * A5
        + 49. / 12. * C_F * C_A * T_F * T_F * A4
        - 253. / 96. * C_F * C_A * T_F * T_F * ZETA_5
        - 65. / 192. * C_F * C_A * T_F * T_F * ZETA_3
        - 15557. / 5184. * C_F * C_F * T_F * T_F
        - 8. * C_F * C_F * T_F * T_F * A5
        - 49. / 6. * C_F * C_F * T_F * T_F * A4
        + 157. / 16. * C_F * C_F * T_F * T_F * ZETA_5
        - 5. / 8. * C_F * C_F * T_F * T_F * ZETA_3
        + 23. / 1440. * PI4 * C_F * T_F * T_F * T_F
        - 181. / 3840. * PI4 * C_F * C_A * T_F * T_F
        + 1639. / 17280. * PI4 * C_F * C_F * T_F * T_F
        + 49. / 720. * LN_2 * PI4 * C_F * C_A * T_F * T_F
        - 49. / 360. * LN_2 * PI4 * C_F * C_F * T_F * T_F
        - 49. / 288. * LN2_2 * PI2 * C_F * C_A * T_F * T_F
        + 49. / 144. * LN2_2 * PI2 * C_F * C_F * T_F * T_F
        + 1. / 18. * LN3_2 * PI2 * C_F * C_A * T_F * T_F
        - 1. / 9. * LN3_2 * PI2 * C_F * C_F * T_F * T_F
        + 49. / 288. * LN4_2 * C_F * C_A * T_F * T_F
        - 49. / 144. * LN4_2 * C_F * C_F * T_F * T_F
        - 1. / 30. * LN5_2 * C_F * C_A * T_F * T_F
        + 1. / 15. * LN5_2 * C_F * C_F * T_F * T_F;

    res[4][1][1] = 2543. / 3888. * C_F * T_F * T_F * T_F
        - 7. / 9. * C_F * T_F * T_F * T_F * ZETA_3
        - 5155. / 31104. * C_F * C_A * T_F * T_F
        + 2. * C_F * C_A * T_F * T_F * A4
        + 997. / 576. * C_F * C_A * T_F * T_F * ZETA_3
        + 319. / 192. * C_F * C_F * T_F * T_F
        - 4. * C_F * C_F * T_F * T_F * A4
        - 97. / 32. * C_F * C_F * T_F * T_F * ZETA_3
        - 43. / 1440. * PI4 * C_F * C_A * T_F * T_F
        + 7. / 180. * PI4 * C_F * C_F * T_F * T_F
        - 1. / 12. * LN2_2 * PI2 * C_F * C_A * T_F * T_F
        + 1. / 6. * LN2_2 * PI2 * C_F * C_F * T_F * T_F
        + 1. / 12. * LN4_2 * C_F * C_A * T_F * T_F
        - 1. / 6. * LN4_2 * C_F * C_F * T_F * T_F;

    res[4][1][2] = 35. / 864. * C_F * T_F * T_F * T_F
        - 2581. / 3456. * C_F * C_A * T_F * T_F
        - 9. / 64. * C_F * C_F * T_F * T_F;

    res[4][1][3] =
        -13. / 96. * C_F * C_A * T_F * T_F + 1. / 8. * C_F * C_F * T_F * T_F;

    res[4][1][4] =
        1. / 144. * C_F * T_F * T_F * T_F - 11. / 288. * C_F * C_A * T_F * T_F;

    res[4][2][0] = 17671. / 20736. * C_F * T_F * T_F * T_F
        - 5. / 144. * C_F * T_F * T_F * T_F * ZETA_3
        - 7. / 1440. * PI4 * C_F * T_F * T_F * T_F;

    res[4][2][1] = -3401. / 7776. * C_F * T_F * T_F * T_F
        + 7. / 18. * C_F * T_F * T_F * T_F * ZETA_3;

    res[4][2][2] = 31. / 216. * C_F * T_F * T_F * T_F;

    res[4][2][4] = 1. / 144. * C_F * T_F * T_F * T_F;
    res
};

/// Coefficients of ζ_m for the scheme-invariant mass m(m)
///
/// `ZETA_M_SI_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_M_SI_COEFF: ZetaCoeffTable = {
    let mut res = ZETA_M_MS_COEFF;

    res[3][0][1] = 35. / 432. * C_F * T_F * T_F + 5. / 64. * C_F * C_A * T_F
        - 3. / 4. * C_F * C_A * T_F * ZETA_3
        - 33. / 64. * C_F * C_F * T_F
        + 3. / 4. * C_F * C_F * T_F * ZETA_3;

    res[3][0][2] = 5. / 72. * C_F * T_F * T_F
        + 29. / 96. * C_F * C_A * T_F
        + 1. / 8. * C_F * C_F * T_F;

    res[4][0][1] = 1. / 4. * D_44_FF / N_R
        - 15. / 8. * D_44_FF / N_R * ZETA_3
        - 83. / 2592. * C_F * T_F * T_F * T_F
        + 1. / 18. * C_F * T_F * T_F * T_F * ZETA_3
        - 107467. / 62208. * C_F * C_A * T_F * T_F
        + 5053. / 4608. * C_F * C_A * T_F * T_F * ZETA_3
        + 233903. / 248832. * C_F * C_A * C_A * T_F
        - 11. / 2. * C_F * C_A * C_A * T_F * A4
        + 25. / 16. * C_F * C_A * C_A * T_F * ZETA_5
        - 119723. / 18432. * C_F * C_A * C_A * T_F * ZETA_3
        - 1693. / 6912. * C_F * C_F * T_F * T_F
        + 181. / 256. * C_F * C_F * T_F * T_F * ZETA_3
        - 11287. / 2304. * C_F * C_F * C_A * T_F
        + 11. * C_F * C_F * C_A * T_F * A4
        + 5. / 16. * C_F * C_F * C_A * T_F * ZETA_5
        + 7761. / 1024. * C_F * C_F * C_A * T_F * ZETA_3
        - 47. / 768. * C_F * C_F * C_F * T_F
        - 15. / 8. * C_F * C_F * C_F * T_F * ZETA_5
        + 457. / 256. * C_F * C_F * C_F * T_F * ZETA_3
        + 1. / 240. * PI4 * C_F * C_A * T_F * T_F
        + 407. / 5760. * PI4 * C_F * C_A * C_A * T_F
        - 1. / 240. * PI4 * C_F * C_F * T_F * T_F
        - 55. / 576. * PI4 * C_F * C_F * C_A * T_F
        + 11. / 48. * LN2_2 * PI2 * C_F * C_A * C_A * T_F
        - 11. / 24. * LN2_2 * PI2 * C_F * C_F * C_A * T_F
        - 11. / 48. * LN4_2 * C_F * C_A * C_A * T_F
        + 11. / 24. * LN4_2 * C_F * C_F * C_A * T_F;

    res[4][0][2] = -35. / 864. * C_F * T_F * T_F * T_F
        + 3103. / 6912. * C_F * C_A * T_F * T_F
        + 3. / 8. * C_F * C_A * T_F * T_F * ZETA_3
        + 19867. / 13824. * C_F * C_A * C_A * T_F
        - 33. / 32. * C_F * C_A * C_A * T_F * ZETA_3
        + 997. / 2304. * C_F * C_F * T_F * T_F
        - 3. / 8. * C_F * C_F * T_F * T_F * ZETA_3
        - 85. / 192. * C_F * C_F * C_A * T_F
        + 33. / 32. * C_F * C_F * C_A * T_F * ZETA_3
        + 33. / 64. * C_F * C_F * C_F * T_F;

    res[4][0][3] = -5. / 216. * C_F * T_F * T_F * T_F
        - 13. / 96. * C_F * C_A * T_F * T_F
        + 17. / 48. * C_F * C_A * C_A * T_F
        - 17. / 192. * C_F * C_F * T_F * T_F
        + 55. / 384. * C_F * C_F * C_A * T_F;

    res[4][1][1] = 2543. / 3888. * C_F * T_F * T_F * T_F
        - 7. / 9. * C_F * T_F * T_F * T_F * ZETA_3
        - 5155. / 31104. * C_F * C_A * T_F * T_F
        + 2. * C_F * C_A * T_F * T_F * A4
        + 997. / 576. * C_F * C_A * T_F * T_F * ZETA_3
        + 689. / 576. * C_F * C_F * T_F * T_F
        - 4. * C_F * C_F * T_F * T_F * A4
        - 97. / 32. * C_F * C_F * T_F * T_F * ZETA_3
        - 43. / 1440. * PI4 * C_F * C_A * T_F * T_F
        + 7. / 180. * PI4 * C_F * C_F * T_F * T_F
        - 1. / 12. * LN2_2 * PI2 * C_F * C_A * T_F * T_F
        + 1. / 6. * LN2_2 * PI2 * C_F * C_F * T_F * T_F
        + 1. / 12. * LN4_2 * C_F * C_A * T_F * T_F
        - 1. / 6. * LN4_2 * C_F * C_F * T_F * T_F;

    res[4][1][2] = 35. / 864. * C_F * T_F * T_F * T_F
        - 2581. / 3456. * C_F * C_A * T_F * T_F
        - 37. / 192. * C_F * C_F * T_F * T_F;

    res[4][1][3] =
        -13. / 96. * C_F * C_A * T_F * T_F - 1. / 16. * C_F * C_F * T_F * T_F;

    res
};

/// Coefficients of ζ_m^-1 in the on-shell scheme
///
/// `ZETA_M_OS_INV_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
// TODO: general colour factors
pub const ZETA_M_OS_INV_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[2][0][0] = -89. / 432.;
    res[2][0][1] = 5. / 36.;
    res[2][0][2] = -1. / 12.;
    res[3][0][0] =
        -1871. / 2916. + B4 / 36. + 407. / 864. * ZETA_3 - 5. / 4. * ZETA_4;
    res[3][0][1] = -299. / 2592. + 5. / 6. * ZETA_3;
    res[3][0][2] = -299. / 432.;
    res[3][0][3] = -35. / 216.;
    res[3][1][0] = -1327. / 11664. + 2. / 27. * ZETA_3;
    res[3][1][1] = 53. / 432.;
    res[3][1][3] = 1. / 108.;
    res[4][0][0] = -122722873709. / 3292047360. + 1924649. / 4354560. * C4
        - 59. / 4860. * C5
        + 5. / 9. * ZETA_2
        + 5. / 27. * LN_2 * ZETA_2
        + 441628331. / 40642560. * ZETA_3
        - 1061. / 576. * ZETA_3 * ZETA_3
        - 16187201. / 580608. * ZETA_4
        + 725. / 108. * LN_2 * ZETA_4
        + 59015. / 1728. * ZETA_5
        + 3935. / 432. * ZETA_6;
    res[4][0][1] = 99319. / 41472. + 11. / 72. * C4
        - 2. / 3. * ZETA_2
        - 2. / 9. * LN_2 * ZETA_2
        + 142463. / 9216. * ZETA_3
        - 481. / 32. * ZETA_4
        - 575. / 72. * ZETA_5;
    res[4][0][2] = -34297. / 3456. + 175. / 48. * ZETA_3;
    res[4][0][3] = -6941. / 2592.;
    res[4][0][4] = -1147. / 3456.;
    res[4][1][0] = 2403419. / 746496. - 49. / 2592. * C4 + C5 / 810.
        - 5. / 54. * ZETA_2
        + 1075. / 1728. * ZETA_3
        + 1225. / 3456. * ZETA_4
        - 49. / 72. * LN_2 * ZETA_4
        - 497. / 288. * ZETA_5;
    res[4][1][1] = -26443. / 93312. - C4 / 108. + ZETA_2 / 9.
        - 599. / 1728. * ZETA_3
        + 163. / 144. * ZETA_4;
    res[4][1][2] = 10237. / 10368.;
    res[4][1][3] = 47. / 288.;
    res[4][1][4] = 17. / 432.;
    res[4][2][0] = -17671. / 124416. + 5. / 864. * ZETA_3 + 7. / 96. * ZETA_4;
    res[4][2][1] = 3401. / 46656. - 7. / 108. * ZETA_3;
    res[4][2][2] = -31. / 1296.;
    res[4][2][4] = -1. / 864.;
    res
};

/// Coefficients of ζ_m^-1 in the MSbar scheme
///
/// `ZETA_M_MS_INV_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_M_MS_INV_COEFF: ZetaCoeffTable = {
    let mut res = ZERO_ZETA_COEFF_TABLE;
    res[0][0][0] = 1.;
    res[2][0][0] = -89. / 432.;
    res[2][0][1] = 5. / 36.;
    res[2][0][2] = -1. / 12.;
    res[3][0][0] =
        -2951. / 2916. + B4 / 36. + 407. / 864. * ZETA_3 - 5. / 4. * ZETA_4;
    res[3][0][1] = 133. / 2592. + 5. / 6. * ZETA_3;
    res[3][0][2] = -155. / 432.;
    res[3][0][3] = -35. / 216.;
    res[3][1][0] = -1327. / 11664. + 2. / 27. * ZETA_3;
    res[3][1][1] = 53. / 432.;
    res[3][1][3] = 1. / 108.;
    res[4][0][0] = -131621265869. / 3292047360. + 1924649. / 4354560. * C4
        - 59. / 4860. * C5
        + 353193131. / 40642560. * ZETA_3
        - 1061. / 576. * ZETA_3 * ZETA_3
        - 16187201. / 580608. * ZETA_4
        + 725. / 108. * LN_2 * ZETA_4
        + 59015. / 1728. * ZETA_5
        + 3935. / 432. * ZETA_6;
    res[4][0][1] = 279367. / 41472. + 11. / 72. * C4 + 42197. / 3072. * ZETA_3
        - 481. / 32. * ZETA_4
        - 575. / 72. * ZETA_5;
    res[4][0][2] = -16193. / 3456. + 175. / 48. * ZETA_3;
    res[4][0][3] = -3161. / 2592.;
    res[4][0][4] = -1147. / 3456.;
    res[4][1][0] = 2261435. / 746496. - 49. / 2592. * C4
        + C5 / 810.
        + 1075. / 1728. * ZETA_3
        + 1225. / 3456. * ZETA_4
        - 49. / 72. * LN_2 * ZETA_4
        - 497. / 288. * ZETA_5;
    res[4][1][1] = -55315. / 93312. - C4 / 108. - 599. / 1728. * ZETA_3
        + 163. / 144. * ZETA_4;
    res[4][1][2] = 8461. / 10368.;
    res[4][1][3] = 23. / 288.;
    res[4][1][4] = 17. / 432.;
    res[4][2][0] = -17671. / 124416. + 5. / 864. * ZETA_3 + 7. / 96. * ZETA_4;
    res[4][2][1] = 3401. / 46656. - 7. / 108. * ZETA_3;
    res[4][2][2] = -31. / 1296.;
    res[4][2][4] = -1. / 864.;
    res
};

/// Coefficients of ζ_m^-1 for a scheme-invariant mass m(m)
///
/// `ZETA_M_SI_INV_COEFF[as_pow][nf_pow][log_pow]` is the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow * log(μ^2/m^2)^log_pow`.
pub const ZETA_M_SI_INV_COEFF: ZetaCoeffTable = {
    let mut res = ZETA_M_MS_INV_COEFF;
    res[3][0][1] = 853. / 2592. + 5. / 6. * ZETA_3;
    res[3][0][2] = -299. / 432.;
    res[4][0][1] = 330503. / 41472. + 11. / 72. * C4 + 47317. / 3072. * ZETA_3
        - 481. / 32. * ZETA_4
        - 575. / 72. * ZETA_5;
    res[4][0][2] = -8531. / 1152. + 175. / 48. * ZETA_3;
    res[4][0][3] = -6941. / 2592.;
    res[4][1][1] = -36019. / 93312. - C4 / 108. - 599. / 1728. * ZETA_3
        + 163. / 144. * ZETA_4;
    res[4][1][2] = 8701. / 10368.;
    res[4][1][3] = 47. / 288.;
    res
};

/// Mass decoupling constant
///
/// The decoupling constant is defined by m^(nf-1) = ζ_m * m^(nf).
/// That is, it allows to determine the mass in the
/// nf-1 flavour theory from the mass in QCD with
/// nf flavours.
/// `log_μ2_m2 = log(μ^2/m^2)`.
pub const fn ζ_m(
    scheme: MassScheme,
    α_s_nf: f64,
    nf: i32,
    log_μ2_m2: f64,
    order: usize,
) -> f64 {
    match scheme {
        MassScheme::OS => ζ_m_os(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::MS => ζ_m_ms(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::SI => ζ_m_si(α_s_nf, nf, log_μ2_m2, order),
    }
}

/// Inverse mass decoupling constant
///
/// The decoupling constant is defined by m^(nf+1) = ζ_m^-1 * m^(nf).
/// That is, it allows to determine the mass in the
/// nf+1 flavour theory from the mass in QCD with
/// nf flavours.
/// `log_μ2_m2 = log(μ^2/m^2)`.
pub const fn ζ_m_inv(
    scheme: MassScheme,
    α_s_nf: f64,
    nf: i32,
    log_μ2_m2: f64,
    order: usize,
) -> f64 {
    match scheme {
        MassScheme::OS => ζ_m_os_inv(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::MS => ζ_m_ms_inv(α_s_nf, nf, log_μ2_m2, order),
        MassScheme::SI => ζ_m_si_inv(α_s_nf, nf, log_μ2_m2, order),
    }
}

/// Alias for [ζ_m]
pub const fn zeta_m(
    scheme: MassScheme,
    α_s_nf: f64,
    nf: i32,
    log_μ2_m2: f64,
    order: usize,
) -> f64 {
    ζ_m(scheme, α_s_nf, nf, log_μ2_m2, order)
}

/// Alias for [ζ_m_inv]
pub const fn zeta_m_inv(
    scheme: MassScheme,
    α_s_nf: f64,
    nf: i32,
    log_μ2_m2: f64,
    order: usize,
) -> f64 {
    ζ_m_inv(scheme, α_s_nf, nf, log_μ2_m2, order)
}

macro_rules! impl_ζ_m {
    ($scheme:ident) => {
        paste! {
            pub const fn [<ζ_m_ $scheme:lower>](
                α_s_nf: f64,
                nf: i32,
                log_μ2_m2: f64,
                order: usize,
            ) -> f64 {
                let α_s_π = α_s_nf / PI;
                let mut res = 0f64;
                const_for! { order in (2..order).rev() => {
                    let zeta_coeff = &[<ZETA_M_ $scheme _COEFF>][order];
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
                1. + α_s_π*α_s_π*res
            }

            pub const fn [<zeta_m_ $scheme:lower>](
                α_s_nf: f64,
                nf: i32,
                log_μ2_m2: f64,
                order: usize,
            ) -> f64 {
                [<ζ_m_ $scheme:lower>](α_s_nf, nf, log_μ2_m2, order)
            }
        }
    };
}

impl_ζ_m!(OS);
impl_ζ_m!(MS);
impl_ζ_m!(SI);
impl_ζ_m!(OS_INV);
impl_ζ_m!(MS_INV);
impl_ζ_m!(SI_INV);

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn test_ζ_m_os() {
        const EP: f64 = 1e-13;
        assert_abs_diff_eq!(ZETA_M_OS_COEFF[0][0][0], 1., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_OS_COEFF[1][0][0], 0., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_OS_COEFF[1][0][1], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[2][0][0],
            0.20601851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[2][0][1],
            -0.1388888888888889,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[2][0][2],
            0.08333333333333333,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[3][0][0],
            1.4525287624862648,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[3][0][1],
            -0.8323468020157112,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[3][0][2],
            0.7384259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[3][0][3],
            0.14351851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[3][1][0],
            0.02472760936815077,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[3][1][1],
            -0.12268518518518519,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_OS_COEFF[3][1][2], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[3][1][3],
            -0.009259259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][0][0],
            1.783112700960487,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][0][1],
            3.325894126200277,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][0][2],
            7.176614089719983,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][0][3],
            2.328317901234568,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][0][4],
            0.30063657407407407,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][1][0],
            -1.6159112400964275,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][1][1],
            -0.6570014862863435,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][1][2],
            -0.9738618827160493,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][1][3],
            -0.16319444444444445,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][1][4],
            -0.037037037037037035,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][2][0],
            0.05615582122474318,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][2][1],
            0.005015862378999758,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][2][2],
            0.023919753086419752,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_OS_COEFF[4][2][3], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_OS_COEFF[4][2][4],
            0.0011574074074074073,
            epsilon = EP
        );
    }

    #[test]
    fn test_ζ_m_ms() {
        const EP: f64 = 1e-13;

        assert_abs_diff_eq!(ZETA_M_MS_COEFF[0][0][0], 1., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_MS_COEFF[1][0][0], 0., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_MS_COEFF[1][0][1], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[2][0][0],
            0.20601851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[2][0][1],
            -0.1388888888888889,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[2][0][2],
            0.08333333333333333,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[3][0][0],
            1.822899132856635,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[3][0][1],
            -0.9990134686823778,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[3][0][2],
            0.4050925925925926,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[3][0][3],
            0.14351851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[3][1][0],
            0.02472760936815077,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[3][1][1],
            -0.12268518518518519,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_MS_COEFF[3][1][2], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[3][1][3],
            -0.009259259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][0][0],
            8.371928509821231,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][0][1],
            -0.9588086344447362,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][0][2],
            1.8317066823125758,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][0][3],
            1.0088734567901234,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][0][4],
            0.30063657407407407,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][1][0],
            -1.5780193327058316,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][1][1],
            -0.1648174541674049,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][1][2],
            -0.8025655864197531,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][1][3],
            -0.0798611111111111,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][1][4],
            -0.037037037037037035,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][2][0],
            0.05615582122474318,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][2][1],
            0.005015862378999758,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][2][2],
            0.023919753086419752,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_MS_COEFF[4][2][3], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_MS_COEFF[4][2][4],
            0.0011574074074074073,
            epsilon = EP
        );
    }

    #[test]
    fn test_ζ_m_si() {
        const EP: f64 = 1e-13;

        assert_abs_diff_eq!(ZETA_M_SI_COEFF[0][0][0], 1., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_SI_COEFF[1][0][0], 0., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_SI_COEFF[1][0][1], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[2][0][0],
            0.20601851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[2][0][1],
            -0.1388888888888889,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[2][0][2],
            0.08333333333333333,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[3][0][0],
            1.822899132856635,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[3][0][1],
            -1.2767912464601556,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[3][0][2],
            0.7384259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[3][0][3],
            0.14351851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[3][1][0],
            0.02472760936815077,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[3][1][1],
            -0.12268518518518519,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_SI_COEFF[3][1][2], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[3][1][3],
            -0.009259259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][0][0],
            8.371928509821231,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][0][1],
            -4.125817053290973,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][0][2],
            4.806243719349613,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][0][3],
            2.328317901234568,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][0][4],
            0.30063657407407407,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][1][0],
            -1.5780193327058316,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][1][1],
            -0.37160757762419505,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][1][2],
            -0.8257137345679012,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][1][3],
            -0.16319444444444445,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][1][4],
            -0.037037037037037035,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][2][0],
            0.05615582122474318,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][2][1],
            0.005015862378999758,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][2][2],
            0.023919753086419752,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_SI_COEFF[4][2][3], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_SI_COEFF[4][2][4],
            0.0011574074074074073,
            epsilon = EP
        );
    }

    #[test]
    fn test_ζ_m_os_inv() {
        const EP: f64 = 1e-13;
        assert_abs_diff_eq!(ZETA_M_OS_INV_COEFF[0][0][0], 1., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_OS_INV_COEFF[1][0][0], 0., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_OS_INV_COEFF[1][0][1], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[2][0][0],
            -0.20601851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[2][0][1],
            0.1388888888888889,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[2][0][2],
            -0.08333333333333333,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[3][0][0],
            -1.4772563718544156,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[3][0][1],
            0.8863591476947236,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[3][0][2],
            -0.6921296296296297,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[3][0][3],
            -0.16203703703703703,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[3][1][0],
            -0.02472760936815077,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[3][1][1],
            0.12268518518518519,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_OS_INV_COEFF[3][1][2], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[3][1][3],
            0.009259259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][0][0],
            -0.30109112125205684,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][0][1],
            -3.7149415236110848,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][0][2],
            -5.541401336860276,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][0][3],
            -2.677854938271605,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][0][4],
            -0.33188657407407407,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][1][0],
            1.5035995976469416,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][1][1],
            0.6346059568442687,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][1][2],
            0.9873649691358025,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][1][3],
            0.16319444444444445,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][1][4],
            0.03935185185185185,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][2][0],
            -0.05615582122474318,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][2][1],
            -0.005015862378999758,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][2][2],
            -0.023919753086419752,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_OS_INV_COEFF[4][2][3], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_OS_INV_COEFF[4][2][4],
            -0.0011574074074074073,
            epsilon = EP
        );
    }

    #[test]
    fn test_ζ_m_ms_inv() {
        const EP: f64 = 1e-13;
        assert_abs_diff_eq!(ZETA_M_MS_INV_COEFF[0][0][0], 1., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_MS_INV_COEFF[1][0][0], 0., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_MS_INV_COEFF[1][0][1], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[2][0][0],
            -0.20601851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[2][0][1],
            0.1388888888888889,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[2][0][2],
            -0.08333333333333333,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[3][0][0],
            -1.8476267422247858,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[3][0][1],
            1.0530258143613902,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[3][0][2],
            -0.3587962962962963,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[3][0][3],
            -0.16203703703703703,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[3][1][0],
            -0.02472760936815077,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[3][1][1],
            0.12268518518518519,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_MS_INV_COEFF[3][1][2], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[3][1][3],
            0.009259259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][0][0],
            -6.744671265486929,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][0][1],
            -0.09371909138130619,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][0][2],
            -0.30297541093434965,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][0][3],
            -1.2195216049382716,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][0][4],
            -0.33188657407407407,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][1][0],
            1.4657076902563455,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][1][1],
            0.14242192472532988,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][1][2],
            0.8160686728395061,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][1][3],
            0.0798611111111111,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][1][4],
            0.03935185185185185,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][2][0],
            -0.05615582122474318,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][2][1],
            -0.005015862378999758,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][2][2],
            -0.023919753086419752,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_MS_INV_COEFF[4][2][3], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_MS_INV_COEFF[4][2][4],
            -0.0011574074074074073,
            epsilon = EP
        );
    }

    #[test]
    fn test_ζ_m_si_inv() {
        const EP: f64 = 1e-13;
        assert_abs_diff_eq!(ZETA_M_SI_INV_COEFF[0][0][0], 1., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_SI_INV_COEFF[1][0][0], 0., epsilon = EP);
        assert_abs_diff_eq!(ZETA_M_SI_INV_COEFF[1][0][1], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[2][0][0],
            -0.20601851851851852,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[2][0][1],
            0.1388888888888889,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[2][0][2],
            -0.08333333333333333,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[3][0][0],
            -1.8476267422247858,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[3][0][1],
            1.3308035921391679,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[3][0][2],
            -0.6921296296296297,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[3][0][3],
            -0.16203703703703703,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[3][1][0],
            -0.02472760936815077,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[3][1][1],
            0.12268518518518519,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_SI_INV_COEFF[3][1][2], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[3][1][3],
            0.009259259259259259,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][0][0],
            -6.744671265486929,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][0][1],
            3.142733771909377,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][0][2],
            -3.022882818341757,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][0][3],
            -2.677854938271605,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][0][4],
            -0.33188657407407407,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][1][0],
            1.4657076902563455,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][1][1],
            0.34921204818212004,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][1][2],
            0.8392168209876543,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][1][3],
            0.16319444444444445,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][1][4],
            0.03935185185185185,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][2][0],
            -0.05615582122474318,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][2][1],
            -0.005015862378999758,
            epsilon = EP
        );
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][2][2],
            -0.023919753086419752,
            epsilon = EP
        );
        assert_abs_diff_eq!(ZETA_M_SI_INV_COEFF[4][2][3], 0., epsilon = EP);
        assert_abs_diff_eq!(
            ZETA_M_SI_INV_COEFF[4][2][4],
            -0.0011574074074074073,
            epsilon = EP
        );
    }
}
