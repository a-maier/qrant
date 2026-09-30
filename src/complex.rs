use std::hash::Hash;

use num_complex::Complex;
use ordered_float::OrderedFloat;

#[derive(Copy, Clone, Debug, Default)]
pub(crate) struct OrderedComplex<T>(pub Complex<OrderedFloat<T>>);

macro_rules! impl_cmp {
    ($t:ty) => {
        impl PartialEq for OrderedComplex<$t> {
            fn eq(&self, other: &Self) -> bool {
                self.0.re == other.0.re && self.0.im == other.0.im
            }
        }

        impl Eq for OrderedComplex<$t> {}

        impl Ord for OrderedComplex<$t> {
            fn cmp(&self, other: &Self) -> std::cmp::Ordering {
                let lhs = self.0;
                let rhs = other.0;
                lhs.to_polar().cmp(&rhs.to_polar())
            }
        }

        impl PartialOrd for OrderedComplex<$t> {
            fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl Hash for OrderedComplex<$t> {
            fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
                self.0.hash(state);
            }
        }

        impl From<$t> for OrderedComplex<$t> {
            fn from(f: $t) -> Self {
                OrderedFloat(f).into()
            }
        }

        impl From<OrderedFloat<$t>> for OrderedComplex<$t> {
            fn from(f: OrderedFloat<$t>) -> Self {
                Self(num_complex::Complex::from(f))
            }
        }

        impl From<Complex<$t>> for OrderedComplex<$t> {
            fn from(f: Complex<$t>) -> Self {
                let Complex { re, im } = f;
                let [re, im] = [re, im].map(|r| OrderedFloat(r));
                Complex { re, im }.into()
            }
        }

        impl From<Complex<OrderedFloat<$t>>> for OrderedComplex<$t> {
            fn from(f: Complex<OrderedFloat<$t>>) -> Self {
                Self(f)
            }
        }

        impl From<OrderedComplex<$t>> for Complex<OrderedFloat<$t>> {
            fn from(f: OrderedComplex<$t>) -> Self {
                f.0
            }
        }

        impl From<OrderedComplex<$t>> for Complex<$t> {
            fn from(f: OrderedComplex<$t>) -> Self {
                Complex {
                    re: f.0.re.into(),
                    im: f.0.im.into(),
                }
            }
        }
    };
}

impl_cmp!(f32);
impl_cmp!(f64);

pub(crate) type C32 = OrderedComplex<f32>;
pub(crate) type C64 = OrderedComplex<f64>;
