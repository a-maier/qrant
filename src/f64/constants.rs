use std::f64::consts::{LN_2, PI};

pub(crate) const ZETA_2: f64 = PI * PI / 6.;
pub(crate) const ZETA_3: f64 = 1.202056903159594285399738161511449990765;
pub(crate) const ZETA_4: f64 = PI * PI * PI * PI / 90.;
pub(crate) const ZETA_5: f64 = 1.036927755143369926331365486457034168057;
pub(crate) const ZETA_6: f64 = 1.017343061984449139714517929790920527902;
pub(crate) const ZETA_7: f64 = 1.008349277381922826839797549849796759600;

pub(crate) const LN2_2: f64 = LN_2 * LN_2;
pub(crate) const LN3_2: f64 = LN2_2 * LN_2;
pub(crate) const LN4_2: f64 = LN3_2 * LN_2;
pub(crate) const LN5_2: f64 = LN4_2 * LN_2;
pub(crate) const A4: f64 = 0.5174790616738993863307581618988629456224;
pub(crate) const A5: f64 = 0.5084005792422687074591088492585899413195;

pub(crate) const PI2: f64 = PI * PI;
pub(crate) const PI4: f64 = PI2 * PI2;
pub(crate) const PI6: f64 = PI4 * PI2;

#[allow(clippy::excessive_precision)]
pub(crate) const B4: f64 = -1.762800087073770864061897634679818807215;
// C4 = 24*Li_4(1/2) + log^4(2) - pi^2*log^2(2)
#[allow(clippy::excessive_precision)]
pub(crate) const C4: f64 = 7.908451398072941071188189589256658841231;
// C5 = -360*Li_5(1/2) + 3*log^5(2) - 5*pi^2*log^3(2)
#[allow(clippy::excessive_precision)]
pub(crate) const C5: f64 = -198.9783082886512628182933723717652048993;
