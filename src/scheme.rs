/// Supported mass renormalization schemes
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum MassScheme {
    /// On-shell scheme
    OS,
    /// MS-bar scheme
    MS,
    /// Scale-invariant scheme
    ///
    /// The scale-invariant mass is the MS-bar mass
    /// at the scale of the MS-bar mass
    SI,
}

impl MassScheme {
    pub const POLE: Self = Self::OS;
    pub const ON_SHELL: Self = Self::OS;
    pub const MS_BAR: Self = Self::MS;
    pub const SCALE_INVARIANT: Self = Self::SI;
}
