//! helper module for working with points in euclidean space

#[repr(transparent)]
#[derive(Clone, Copy)]
pub struct VecN<const L: usize>([f32; L]);

impl<const L: usize> From<[f32; L]> for VecN<L> {
    fn from(value: [f32; L]) -> Self {
        Self(value)
    }
}

impl<const L: usize> From<VecN<L>> for [f32; L] {
    fn from(value: VecN<L>) -> Self {
        value.0
    }
}

impl<const L: usize> std::ops::Add for VecN<L> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(std::array::from_fn(|x| self.0[x] + rhs.0[x]))
    }
}

impl<const L: usize> std::ops::Sub for VecN<L> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self(std::array::from_fn(|x| self.0[x] - rhs.0[x]))
    }
}

impl<const L: usize> std::ops::Mul for VecN<L> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(std::array::from_fn(|x| self.0[x] * rhs.0[x]))
    }
}

impl<const L: usize> std::ops::Div for VecN<L> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self(std::array::from_fn(|x| self.0[x] / rhs.0[x]))
    }
}

impl<const L: usize> std::ops::Mul<f32> for VecN<L> {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        Self(std::array::from_fn(|x| self.0[x] * rhs))
    }
}

impl<const L: usize> std::ops::Div<f32> for VecN<L> {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        Self(std::array::from_fn(|x| self.0[x] / rhs))
    }
}

impl<const L: usize> std::fmt::Debug for VecN<L> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}


impl<const L: usize> VecN<L> {
    pub fn mag(&self) -> f32 {
        self.sqmag().sqrt()
    }

    pub fn sqmag(&self) -> f32 {
        self.0.into_iter().map(|x| x * x).sum()
    }

    /// returns true if `self` and `other` are within `dist` of each other
    pub fn within(self, dist: f32, other: Self) -> bool {
        let sqdist = (self - other).sqmag();
        sqdist < dist * dist
    }


    /// returns true if `self` and `other` are within `dist` of each other, panicking otherwise
    #[track_caller]
    pub fn assert_within(self, dist: f32, other: Self) -> bool {
        #[cold]
        #[inline(never)]
        fn assert_fail<const L: usize>(lhs: VecN<L>, dist: f32, rhs: VecN<L>) -> ! {
            panic!(
            "Assert failed, distance is more than {dist} \n\
            Left hand side: {lhs:?}\n\
            Right hand side: {rhs:?}"
            );
        }
        if !self.within(dist, other) {
            assert_fail(self, dist, other);
        }
        true
    }
}
