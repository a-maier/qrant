# qrant

A crate for the evolution and decoupling of parameters in
quantum chromodynamics (QCD) in the MS-bar scheme.

In the MS-bar scheme, the strong coupling and the quark masses
depend on the renormalization scale and the number of active
flavours. This crate allows evolving these parameter values
from one scale to another and to decouple quark flavours from
the theory.

By default, the 2026 [PDG](https://pdg.lbl.gov/) parameter
values are used together with five-loop perturbative evolution
and four-loop decoupling.

## Usage

```rust
// Get the value of the strong coupling constant in the four-flavour theory at 3 GeV.
assert_eq!(α_s_nf4(3.), 0.25319569593921754);
// Same in ASCII
assert_eq!(alpha_s_nf4(3.), α_s_nf4(3.));
// Same, for a number of flavours determined at run time
assert_eq!(α_s(4, 3.), α_s_nf4(3.));

// Same, with explicit running and decoupling settings
let mb = 4.7.into();
let bottom_threshold = Threshold{
    m: mb,
    scheme: MassScheme::OS,
    scale: mb,
};
let α_s = StrongCoupling::builder()
    .value(0.1180)
    .μ(91.1879)
    .nf(5)
    .order(5).unwrap()
    .build()
    .remove_flavour(bottom_threshold);
assert_eq!(α_s.at(3.), α_s_nf4(3.));

// Combined strong coupling and mass evolution starting at 3 GeV
let α_s = α_s.evolve(3.);
let mc_3gev = 0.986;
let qcd = QcdParameters::with_coupling(α_s, mc_3gev);
// scale-invariant mass mc(mc)
assert_eq!(qcd.m_of_m(), 1.2722136772876156);
// parameter values at 10 GeV
let ParameterValues{ α_s, m } = qcd.at(10.);
assert_eq!(α_s, 0.17423377517877828);
assert_eq!(m, 0.800826180563377);
```

## Features

`serde`: Enables serialization of `StrongCoupling` and
         `QcdParameters` objects using
         [serde](https://crates.io/crates/serde).

License: GPL-3.0-or-later
