use std::f64::consts::{LN_2, PI};

pub(crate) const ZETA_2: f64 = PI * PI / 6.;
pub(crate) const ZETA_3: f64 = 1.202056903159594285399738161511449990765;
pub(crate) const ZETA_4: f64 = PI * PI * PI * PI / 90.;
pub(crate) const ZETA_5: f64 = 1.036927755143369926331365486457034168057;

pub(crate) const LN2_2: f64 = LN_2 * LN_2;
pub(crate) const LN3_2: f64 = LN2_2 * LN_2;
pub(crate) const LN4_2: f64 = LN3_2 * LN_2;
pub(crate) const LN5_2: f64 = LN4_2 * LN_2;
pub(crate) const A4: f64 = 0.5174790616738993863307581618988629456224;
pub(crate) const A5: f64 = 0.5084005792422687074591088492585899413195;

pub(crate) const PI2: f64 = PI * PI;
pub(crate) const PI4: f64 = PI2 * PI2;
pub(crate) const PI6: f64 = PI4 * PI2;
