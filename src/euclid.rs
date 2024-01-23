#![allow(dead_code)]
use crate::{*, vn::VecN};

pub type Pos = [real_t; 3];
pub type V3 = vn::VecN<3>;

/// merge two Euclidean vertices. This is conservative for now, and I want to look into improving
/// it in the future.
pub fn merge(p1: Pos, r1: real_t, p2: Pos, r2: real_t) -> (Pos, real_t) {
    merge_generic(p1, r1, p2, r2)
}

fn merge_generic<const L: usize>(p1: [f32; L], r1: real_t, p2: [f32; L], r2: real_t) -> ([f32; L], real_t) {
    let p1v: VecN<L> = p1.into();
    let p2v: VecN<L> = p2.into();
    let dist = (p1v - p2v).mag();
    let retp;
    let retr;
    if (r1 - r2).abs() < dist {
        retp = std::array::from_fn(|i| {
            let x = p1[i] + p2[i] + (r1 * (p1[i] - p2[i]) + r2 * (p2[i] - p1[i])) / dist;
            x/2.0
        });
        retr = (dist + r1 + r2) / 2.0;
    } else {
        if r1 > r2 {
            retp = p1;
        } else {
            retp = p2;
        }
        retr = r1.max(r2);
    }
    (retp, retr)
}


#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn merge_far() {
        let p1 = [-1.0, -1.0];
        let r1 = 1.0;
        let p2 = [2.0, 3.0];
        let r2 = 2.0;
        let (p3, r3) = merge_generic(p1, r1, p2, r2);
        let p3v: VecN<2> = p3.into();
        p3v.assert_within(0.001, [0.8, 1.4].into());
        assert_eq!(r3, 4.0);
    }

    #[test]
    fn merge_overlap() {
        let p1 = [-1.0, -1.0];
        let r1 = 6.0;
        let p2 = [2.0, 3.0];
        let r2 = 2.0;
        let (p3, r3) = merge_generic(p1, r1, p2, r2);
        let p3v: VecN<2> = p3.into();
        p3v.assert_within(0.001, [-0.7, -0.6].into());
        assert_eq!(r3, 6.5);
    }

    #[test]
    fn merge_close() {
        let p1 = [-1.0, -1.0];
        let r1 = 7.0;
        let p2 = [2.0, 3.0];
        let r2 = 1.0;
        let (p3, r3) = merge_generic(p1, r1, p2, r2);
        let p3v: VecN<2> = p3.into();
        p3v.assert_within(0.001, p1.into());
        assert_eq!(r3, 7.0);
    }

    #[test]
    fn merge_close_rev() {
        let p1 = [0.0, 0.0];
        let r1 = 4.0;
        let p2 = [0.0, 3.0];
        let r2 = 10.0;
        let (p3, r3) = merge_generic(p1, r1, p2, r2);
        let p3v: VecN<2> = p3.into();
        p3v.assert_within(0.001, p2.into());
        assert_eq!(r3, 10.0);
    }
}
