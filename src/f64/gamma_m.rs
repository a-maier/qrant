use std::f64::consts::PI;

use super::{
    colour::{C_A, C_F, D_44_AA, D_44_AF, D_44_FF, N_A, N_R, T_F},
    constants::{ZETA_3, ZETA_4, ZETA_5, ZETA_6, ZETA_7},
};

use const_for::const_for;

const MAX_LOOP_ORDER: usize = 5;
const MAX_NF_POW: usize = 4;

pub type GammaMCoeffTable = [[f64; MAX_NF_POW + 1]; MAX_LOOP_ORDER];

const ZERO_GAMMA_M_COEFF_TABLE: GammaMCoeffTable =
    [[0.; MAX_NF_POW + 1]; MAX_LOOP_ORDER];

/// Highest power of `nf` in γ_m_i
///
/// `i == 0` corresponds to *minus* the one-loop coefficient in γ_m
pub const fn max_nf_pow(i: usize) -> usize {
    i
}

/// Coefficients of γ_m
///
/// `GAMMA_M_COEFF[as_pow][nf_pow]` is *minus* the coefficient of
/// `(α_s/π)^as_pow * nf^nf_pow`.
pub const GAMMA_M_COEFF: GammaMCoeffTable = {
    let mut res = ZERO_GAMMA_M_COEFF_TABLE;
    // taken from arXiv:1709.08541
    res[0][0] = (3. * C_F) / 4.;
    res[1][0] = (3. / 2. * C_F * C_F + 97. / 6. * C_A * C_F) / 16.;
    res[1][1] = (-10. / 3. * C_F * T_F) / 16.;
    res[2][0] = (129. / 2. * C_F * C_F * C_F - 129. / 4. * C_A * C_F * C_F
        + 11413. / 108. * C_A * C_A * C_F)
        / 64.;
    res[2][1] = (-46. * C_F * C_F * T_F + 48. * C_F * C_F * T_F * ZETA_3
        - 556. / 27. * C_A * C_F * T_F
        - 48. * C_A * C_F * T_F * ZETA_3)
        / 64.;
    res[2][2] = (-140. / 27. * C_F * T_F * T_F) / 64.;
    res[3][0] = (-32. * D_44_AF / N_R + 240. * D_44_AF / N_R * ZETA_3
        - 1261. / 8. * C_F * C_F * C_F * C_F
        - 336. * C_F * C_F * C_F * C_F * ZETA_3
        + 15349. / 12. * C_A * C_F * C_F * C_F
        + 316. * C_A * C_F * C_F * C_F * ZETA_3
        - 34045. / 36. * C_A * C_A * C_F * C_F
        + 440. * C_A * C_A * C_F * C_F * ZETA_5
        - 152. * C_A * C_A * C_F * C_F * ZETA_3
        + 70055. / 72. * C_A * C_A * C_A * C_F
        - 440. * C_A * C_A * C_A * C_F * ZETA_5
        + 1418. / 9. * C_A * C_A * C_A * C_F * ZETA_3)
        / 256.;
    res[3][1] = (64. * D_44_FF / N_R
        - 480. * D_44_FF / N_R * ZETA_3
        - 280. / 3. * C_F * C_F * C_F * T_F
        - 480. * C_F * C_F * C_F * T_F * ZETA_5
        + 552. * C_F * C_F * C_F * T_F * ZETA_3
        - 8819. / 27. * C_A * C_F * C_F * T_F
        + 80. * C_A * C_F * C_F * T_F * ZETA_5
        - 264. * C_A * C_F * C_F * T_F * ZETA_4
        + 368. * C_A * C_F * C_F * T_F * ZETA_3
        - 65459. / 162. * C_A * C_A * C_F * T_F
        + 400. * C_A * C_A * C_F * T_F * ZETA_5
        + 264. * C_A * C_A * C_F * T_F * ZETA_4
        - 2684. / 3. * C_A * C_A * C_F * T_F * ZETA_3)
        / 256.;
    res[3][2] = (304. / 27. * C_F * C_F * T_F * T_F
        + 96. * C_F * C_F * T_F * T_F * ZETA_4
        - 160. * C_F * C_F * T_F * T_F * ZETA_3
        + 1342. / 81. * C_A * C_F * T_F * T_F
        - 96. * C_A * C_F * T_F * T_F * ZETA_4
        + 160. * C_A * C_F * T_F * T_F * ZETA_3)
        / 256.;
    res[3][3] = (-664. / 81. * C_F * T_F * T_F * T_F
        + 128. / 9. * C_F * T_F * T_F * T_F * ZETA_3)
        / 256.;
    res[4][0] = (6768. * C_F * D_44_AF / N_R
        - 25440. * C_F * D_44_AF / N_R * ZETA_5
        + 3840. * C_F * D_44_AF / N_R * ZETA_3
        + 358. * C_F * D_44_AA / N_A
        - 1232. * C_F * D_44_AA / N_A * ZETA_7
        - 15640. / 3. * C_F * D_44_AA / N_A * ZETA_5
        + 264. * C_F * D_44_AA / N_A * ZETA_4
        + 4988. / 3. * C_F * D_44_AA / N_A * ZETA_3
        + 3872. * C_F * D_44_AA / N_A * ZETA_3 * ZETA_3
        + 50995. / 8. * C_F * C_F * C_F * C_F * C_F
        + 2080. * C_F * C_F * C_F * C_F * C_F * ZETA_5
        + 848. * C_F * C_F * C_F * C_F * C_F * ZETA_3
        - 76784. / 9. * C_A * D_44_AF / N_R
        + 1232. * C_A * D_44_AF / N_R * ZETA_7
        + 23200. * C_A * D_44_AF / N_R * ZETA_5
        - 1320. * C_A * D_44_AF / N_R * ZETA_4
        - 1208. / 3. * C_A * D_44_AF / N_R * ZETA_3
        - 3872. * C_A * D_44_AF / N_R * ZETA_3 * ZETA_3
        - 2565029. / 144. * C_A * C_F * C_F * C_F * C_F
        + 29600. / 3. * C_A * C_F * C_F * C_F * C_F * ZETA_5
        + 1848. * C_A * C_F * C_F * C_F * C_F * ZETA_4
        - 13060. * C_A * C_F * C_F * C_F * C_F * ZETA_3
        + 2625197. / 72. * C_A * C_A * C_F * C_F * C_F
        - 9240. * C_A * C_A * C_F * C_F * C_F * ZETA_7
        - 99140. / 3. * C_A * C_A * C_F * C_F * C_F * ZETA_5
        - 1738. * C_A * C_A * C_F * C_F * C_F * ZETA_4
        + 24174. * C_A * C_A * C_F * C_F * C_F * ZETA_3
        - 25256617. / 972. * C_A * C_A * C_A * C_F * C_F
        + 6160. * C_A * C_A * C_A * C_F * C_F * ZETA_7
        - 12100. / 3. * C_A * C_A * C_A * C_F * C_F * ZETA_6
        + 374180. / 9. * C_A * C_A * C_A * C_F * C_F * ZETA_5
        + 836. * C_A * C_A * C_A * C_F * C_F * ZETA_4
        - 50642. / 3. * C_A * C_A * C_A * C_F * C_F * ZETA_3
        - 1936. / 3. * C_A * C_A * C_A * C_F * C_F * ZETA_3 * ZETA_3
        + 22663417. / 1944. * C_A * C_A * C_A * C_A * C_F
        + 3080. * C_A * C_A * C_A * C_A * C_F * ZETA_7
        + 12100. / 3. * C_A * C_A * C_A * C_A * C_F * ZETA_6
        - 64405. / 3. * C_A * C_A * C_A * C_A * C_F * ZETA_5
        - 15697. / 18. * C_A * C_A * C_A * C_A * C_F * ZETA_4
        + 418936. / 81. * C_A * C_A * C_A * C_A * C_F * ZETA_3
        + 1936. / 3. * C_A * C_A * C_A * C_A * C_F * ZETA_3 * ZETA_3)
        / 1024.;
    res[4][1] = (3040. / 9. * T_F * D_44_AF / N_R
        + 3520. / 3. * T_F * D_44_AF / N_R * ZETA_5
        - 144. * T_F * D_44_AF / N_R * ZETA_4
        - 14920. / 3. * T_F * D_44_AF / N_R * ZETA_3
        - 1408. * T_F * D_44_AF / N_R * ZETA_3 * ZETA_3
        - 352. * C_F * D_44_FF / N_R
        + 11840. * C_F * D_44_FF / N_R * ZETA_5
        - 4928. * C_F * D_44_FF / N_R * ZETA_3
        - 48797. / 36. * C_F * C_F * C_F * C_F * T_F
        + 6720. * C_F * C_F * C_F * C_F * T_F * ZETA_7
        - 37000. / 3. * C_F * C_F * C_F * C_F * T_F * ZETA_5
        - 672. * C_F * C_F * C_F * C_F * T_F * ZETA_4
        + 6888. * C_F * C_F * C_F * C_F * T_F * ZETA_3
        + 25384. / 9. * C_A * D_44_FF / N_R
        - 35680. / 3. * C_A * D_44_FF / N_R * ZETA_5
        + 2640. * C_A * D_44_FF / N_R * ZETA_4
        - 16720. / 3. * C_A * D_44_FF / N_R * ZETA_3
        + 1408. * C_A * D_44_FF / N_R * ZETA_3 * ZETA_3
        - 406861. / 108. * C_A * C_F * C_F * C_F * T_F
        - 3360. * C_A * C_F * C_F * C_F * T_F * ZETA_7
        + 4400. * C_A * C_F * C_F * C_F * T_F * ZETA_6
        + 5620. * C_A * C_F * C_F * C_F * T_F * ZETA_5
        - 2360. * C_A * C_F * C_F * C_F * T_F * ZETA_4
        - 18013. / 9. * C_A * C_F * C_F * C_F * T_F * ZETA_3
        + 1760. * C_A * C_F * C_F * C_F * T_F * ZETA_3 * ZETA_3
        + 753557. / 1944. * C_A * C_A * C_F * C_F * T_F
        + 560. * C_A * C_A * C_F * C_F * T_F * ZETA_7
        + 2200. / 3. * C_A * C_A * C_F * C_F * T_F * ZETA_6
        - 81464. / 9. * C_A * C_A * C_F * C_F * T_F * ZETA_5
        - 3226. * C_A * C_A * C_F * C_F * T_F * ZETA_4
        + 72194. / 9. * C_A * C_A * C_F * C_F * T_F * ZETA_3
        - 5456. / 3. * C_A * C_A * C_F * C_F * T_F * ZETA_3 * ZETA_3
        - 4994047. / 972. * C_A * C_A * C_A * C_F * T_F
        - 3920. * C_A * C_A * C_A * C_F * T_F * ZETA_7
        - 15400. / 3. * C_A * C_A * C_A * C_F * T_F * ZETA_6
        + 139234. / 9. * C_A * C_A * C_A * C_F * T_F * ZETA_5
        + 54925. / 9. * C_A * C_A * C_A * C_F * T_F * ZETA_4
        - 959759. / 81. * C_A * C_A * C_A * C_F * T_F * ZETA_3
        + 176. / 3. * C_A * C_A * C_A * C_F * T_F * ZETA_3 * ZETA_3)
        / 1024.;
    res[4][2] = (-6896. / 9. * T_F * D_44_FF / N_R
        - 2240. / 3. * T_F * D_44_FF / N_R * ZETA_5
        - 768. * T_F * D_44_FF / N_R * ZETA_4
        + 7312. / 3. * T_F * D_44_FF / N_R * ZETA_3
        + 45253. / 54. * C_F * C_F * C_F * T_F * T_F
        - 1600. * C_F * C_F * C_F * T_F * T_F * ZETA_6
        + 8000. / 3. * C_F * C_F * C_F * T_F * T_F * ZETA_5
        + 1304. * C_F * C_F * C_F * T_F * T_F * ZETA_4
        - 38416. / 9. * C_F * C_F * C_F * T_F * T_F * ZETA_3
        + 896. * C_F * C_F * C_F * T_F * T_F * ZETA_3 * ZETA_3
        + 375373. / 486. * C_A * C_F * C_F * T_F * T_F
        + 800. / 3. * C_A * C_F * C_F * T_F * T_F * ZETA_6
        - 12448. / 9. * C_A * C_F * C_F * T_F * T_F * ZETA_5
        + 1864. * C_A * C_F * C_F * T_F * T_F * ZETA_4
        + 5996. / 9. * C_A * C_F * C_F * T_F * T_F * ZETA_3
        - 6976. / 3. * C_A * C_F * C_F * T_F * T_F * ZETA_3 * ZETA_3
        + 27418. / 243. * C_A * C_A * C_F * T_F * T_F
        + 4000. / 3. * C_A * C_A * C_F * T_F * T_F * ZETA_6
        - 2944. / 3. * C_A * C_A * C_F * T_F * T_F * ZETA_5
        - 9364. / 3. * C_A * C_A * C_F * T_F * T_F * ZETA_4
        + 87722. / 27. * C_A * C_A * C_F * T_F * T_F * ZETA_3
        + 4288. / 3. * C_A * C_A * C_F * T_F * T_F * ZETA_3 * ZETA_3)
        / 1024.;
    res[4][3] = (8966. / 81. * C_F * C_F * T_F * T_F * T_F
        + 512. / 3. * C_F * C_F * T_F * T_F * T_F * ZETA_5
        - 320. * C_F * C_F * T_F * T_F * T_F * ZETA_4
        + 352. / 3. * C_F * C_F * T_F * T_F * T_F * ZETA_3
        + 18667. / 243. * C_A * C_F * T_F * T_F * T_F
        - 4096. / 9. * C_A * C_F * T_F * T_F * T_F * ZETA_5
        + 2176. / 9. * C_A * C_F * T_F * T_F * T_F * ZETA_4
        + 21472. / 81. * C_A * C_F * T_F * T_F * T_F * ZETA_3)
        / 1024.;
    res[4][4] = (-1040. / 81. * C_F * T_F * T_F * T_F * T_F
        + 256. / 9. * C_F * T_F * T_F * T_F * T_F * ZETA_4
        - 1280. / 81. * C_F * T_F * T_F * T_F * T_F * ZETA_3)
        / 1024.;
    res
};

/// Perturbative coefficient of the mass anomalous dimension
///
/// The coefficients are defined by γ_m = -Σ γ_m_i (α_s/π)^(i+1).
pub const fn γ_m_i(i: usize, nf: i32) -> f64 {
    let γ_m_coeff = &GAMMA_M_COEFF[i];
    let max_nf_pow = max_nf_pow(i);
    let mut res = γ_m_coeff[max_nf_pow];
    const_for! { nf_pow in (0..max_nf_pow).rev() => {
        res = res.mul_add(nf as f64, γ_m_coeff[nf_pow]);
    }}
    res
}

/// The mass anomalous dimension to the given loop order.
pub const fn γ_m(order: usize, nf: i32, α_s: f64) -> f64 {
    if order == 0 {
        return 0.;
    }
    let i = order - 1;
    let α_s_π = α_s / PI;
    let mut res = -γ_m_i(i, nf);
    const_for! { i in (0..i).rev() => {
        res = res.mul_add(α_s_π, -γ_m_i(i, nf));
    }}
    α_s_π * res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn γ_m_i_nf3() {
        const NF: i32 = 3;
        const EP: f64 = 1e-12;
        assert!((γ_m_i(0, NF) - 1.).abs() < EP);
        assert!((γ_m_i(1, NF) - 3.791666666666666666666666).abs() < EP);
        assert!((γ_m_i(2, NF) - 12.420204964323236).abs() < EP);
        assert!((γ_m_i(3, NF) - 44.26290846573344).abs() < EP);
        assert!((γ_m_i(4, NF) - 198.90657900228448479).abs() < EP);
    }

    #[test]
    fn γ_m_i_nf5() {
        const NF: i32 = 5;
        const EP: f64 = 1e-12;
        assert!((γ_m_i(0, NF) - 1.).abs() < EP);
        assert!((γ_m_i(1, NF) - 3.513888888888889).abs() < EP);
        assert!((γ_m_i(2, NF) - 7.419863212143666).abs() < EP);
        assert!((γ_m_i(3, NF) - 11.034321242068067).abs() < EP);
        assert!((γ_m_i(4, NF) - 41.820456935958403891).abs() < EP);
    }
}
