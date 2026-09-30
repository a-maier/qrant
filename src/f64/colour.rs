/// Dimension of the fundamental representation
pub const N_R: f64 = 3.;

/// Dimension of the adjoint representation
pub const N_A: f64 = 8.;

/// Quadratic Casimir in the fundamental representation
pub const C_F: f64 = 4. / 3.;

/// Quadratic Casimir in the adjoint representation
pub const C_A: f64 = 3.;

/// Trace normalisation
pub const T_F: f64 = 1. / 2.;

/// Product of two symmetrised traces in the adjoint representation
///
/// D_44_AA = (d^abcd_A)^2, where d^abcd_A is the
/// trace over the symmetrised product of four generators in the
/// adjoint representation
pub const D_44_AA: f64 = 135.;

/// Product of two symmetrised traces in the adjoint and fundamental representation
///
/// D_44_AF = (d^abcd_A) (d^abcd_F), where d^abcd_A/d^abcd_F is the
/// trace over the symmetrised product of four generators in the
/// adjoint/fundamental representation
pub const D_44_AF: f64 = 15. / 2.;

/// Product of two symmetrised traces in the fundamental representation
///
/// D_44_FF = (d^abcd_F)^2, where d^abcd_F is the
/// trace over the symmetrised product of four generators in the
/// fundamental representation
pub const D_44_FF: f64 = 5. / 12.;
