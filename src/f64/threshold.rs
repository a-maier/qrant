use ordered_float::OrderedFloat;

use crate::MassScheme;

use super::F64;

/// Threshold for the decoupling of a quark flavour
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Threshold {
    /// Mass of the quark to be decoupled in GeV
    pub m: F64,
    /// Quark mass scheme
    pub scheme: MassScheme,
    /// Decoupling scale in GeV
    pub μ: F64,
}

/// Default top-quark mass in GeV in the on-shell scheme
///
/// Value taken from [PDG](http://arxiv.org/abs/https://doi.org/10.1142/S0217751X26300115).
pub const MT_OS: f64 = 172.6;

/// Default bottom-quark mass in GeV in the on-shell scheme
pub const MB_OS: f64 = 4.7;

/// Default charm-quark mass in GeV in the on-shell scheme
pub const MC_OS: f64 = 1.5;

/// Default threshold for top-quark decoupling
pub const TOP_THRESHOLD: Threshold = Threshold {
    m: OrderedFloat(MT_OS),
    scheme: MassScheme::OS,
    μ: OrderedFloat(MT_OS),
};

/// Default threshold for bottom-quark decoupling
pub const BOTTOM_THRESHOLD: Threshold = Threshold {
    m: OrderedFloat(MB_OS),
    scheme: MassScheme::OS,
    μ: OrderedFloat(MB_OS),
};

/// Default threshold for charm-quark decoupling
pub const CHARM_THRESHOLD: Threshold = Threshold {
    m: OrderedFloat(MC_OS),
    scheme: MassScheme::OS,
    μ: OrderedFloat(MC_OS),
};
