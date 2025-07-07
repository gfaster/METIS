//! Centralized (standardized) cost estimation and measurement

use std::{fmt, ops::{Add, AddAssign}, time::Duration};

use crate::{strategy::Case, PartScheme};


// Units are "something like" milliseconds
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct RuntimeCost(usize);

impl RuntimeCost {
    #[allow(dead_code)]
    pub const ZERO: Self = RuntimeCost(0);
    #[allow(dead_code)]
    pub const MAX: Self = RuntimeCost(usize::MAX);

    pub fn to_raw(self) -> usize {
        self.0
    }
}

impl fmt::Display for RuntimeCost {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "~{} ms", self.0)
    }
}

// Units are ???? (definitely not stable)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct HumanCost(usize);

impl HumanCost {
    /// returns true if `self` is worse than `improved` by at least `factor`
    #[expect(dead_code)]
    pub fn is_worse_than(self, improved: Self, factor: usize) -> bool {
        self.0.saturating_add(50 + 10 * factor) < improved.0
    }
}

impl fmt::Display for HumanCost {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "~{} human thinking units", self.0)
    }
}


#[derive(Debug, Clone, Copy)]
pub struct CaseCost {
    pub runtime: RuntimeCost,
    pub human: HumanCost,
}

/// an opaque abstraction over constraint distribution. Its job is to attempt to quantitize the
/// complexity of different values for weights
///
/// For example:
/// ```
/// DistCost::from(&[1, 2, 3]) < DistCost::from(&[23, 91, 1]);
/// DistCost::from(&[1, 1, 1]) < DistCost::from(&[2, 2, 2]);
/// DistCost::from(&[1, 1, 1]) < DistCost::from(&[1, 1, 1, 1, 1]);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct DistCost {
    total: u64,
    cnt: u32,
    min: u32,
    max: u32,
    /// trying to describe something along the lines of the average deviation from a running
    /// average. Has zero rigour
    eccentricity: f32,
}

impl DistCost {
    pub fn new(arr: &[usize]) -> Self {
        if arr.is_empty() {
            return Self {
                total: 0,
                cnt: 0,
                min: 0,
                max: 0,
                eccentricity: 0.0,
            }
        }

        let mut total = arr[0] as u64;
        let mut cnt = 1;
        let mut min = arr[0] as u32;
        let mut max = arr[0] as u32;

        // https://en.wikipedia.org/wiki/Exponential_smoothing 
        const ALPHA: f64 = 0.2;
        let mut wavg = arr[0] as f64;

        // this metric is just throwing stuff at the wall and seeing what sticks
        let mut total_eccentricity = 0.0;

        for w in arr.windows(2) {
            total += w[1] as u64;
            cnt += 1;
            min = min.min(w[1] as u32);
            max = max.max(w[1] as u32);
            wavg = ((1.0 - ALPHA) * wavg) + (ALPHA * w[1] as f64);
            total_eccentricity += (wavg - w[1] as f64).abs()
        }

        let eccentricity = (total_eccentricity / cnt as f64).ln_1p() as f32;

        Self {
            total,
            cnt,
            min,
            max,
            eccentricity,
        }
    }

    pub fn scale_values(&self, factor: f64) -> Self {
        assert!(factor > 0.0);

        let Self { total, cnt, min, max, eccentricity } = *self;

        let total = (total as f64 * factor).round() as u64;
        let min = ((min as f64 * factor).floor() as u32).max(1);
        let max = ((max as f64 * factor).ceil() as u32).max(1);

        let total_eccentricity = (eccentricity as f64).exp_m1() * cnt as f64;
        let total_eccentricity = total_eccentricity * factor;
        let eccentricity = (total_eccentricity / cnt as f64).ln_1p() as f32;

        Self { total, cnt, min, max, eccentricity }
    }

    pub fn scale_count(&self, factor: f64) -> Self {
        assert!(factor > 0.0);

        let Self { total, cnt, min, max, eccentricity } = *self;

        let new_cnt = ((cnt as f64) * factor) as u32;
        let total = (total as f64 * factor) as u64;

        Self { total, cnt: new_cnt, min, max, eccentricity }
    }

    /// attempt to condense into a single score
    pub fn score(&self) -> f32 {
        let Self { total, cnt, min, max, eccentricity } = *self;

        if cnt == 0 {
            return 1.0
        }

        // zero range implies this is slightly larger than 0.0
        let range_factor = ((2 + max - min) as f32).ln();

        let len_factor = ((1 + cnt) as f32).ln();

        let avg_factor = (total as f64 / cnt as f64).sqrt() as f32;

        let total = range_factor * avg_factor * (eccentricity + 1.0) * len_factor;

        debug_assert!(total > 0.0, "{self:?} has invalid inv score: {total}\n\t\
            range_factor: {range_factor}\n\t\
            avg_factor: {avg_factor}\n\t\
            eccentricity: {eccentricity}\
            ");

        total
    }
}

impl From<&[usize]> for DistCost {
    fn from(value: &[usize]) -> Self {
        DistCost::new(value)
    }
}

impl From<&Vec<usize>> for DistCost {
    fn from(value: &Vec<usize>) -> Self {
        DistCost::new(value)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CostFactors {
    pub nvtxs: u32,
    pub nedges: u32,
    pub nparts: u32,
    pub contig: bool,
    pub minconn: bool,
    pub niter: u32,
    pub ptype: PartScheme,
    pub ncon: u32,
    pub adjwgt: Option<DistCost>,
    pub vsize: Option<DistCost>,
    pub vwgt: Option<DistCost>,
}


impl CostFactors {
    pub fn from_case(case: &Case) -> Self {
        CostFactors { 
            nvtxs: case.graph.nvtxs() as u32,
            nedges: case.graph.nedges() as u32,
            nparts: case.nparts as u32,
            contig: case.contig,
            niter: case.niter as u32,
            ptype: case.ptype,
            ncon: case.graph.ncon as u32,
            minconn: case.minconn,
            adjwgt: case.graph.adjwgt.as_deref().map(DistCost::from),
            vsize: case.graph.vsize.as_deref().map(DistCost::from),
            vwgt: case.graph.vwgt.as_deref().map(DistCost::from),
        }
    }
}

impl From<&Case> for CostFactors {
    fn from(value: &Case) -> Self {
        Self::from_case(value)
    }
}

impl From<CostFactors> for CaseCost {
    fn from(value: CostFactors) -> Self {
        (&value).into()
    }
}

impl From<&CostFactors> for CaseCost {
    fn from(value: &CostFactors) -> Self {
        CaseCost { 
            runtime: estimate_runtime_cost(value),
            human: estimate_human_cost(value),
        }
    }
}

impl From<&Case> for CaseCost {
    fn from(value: &Case) -> Self {
        CostFactors::from_case(value).into()
    }
}

impl From<Duration> for RuntimeCost {
    fn from(value: Duration) -> Self {
        RuntimeCost(value.as_millis() as usize)
    }
}

impl Add for RuntimeCost {
    type Output = RuntimeCost;

    fn add(self, rhs: RuntimeCost) -> Self::Output {
        RuntimeCost(self.0 + rhs.0)
    }
}

impl Add<Duration> for RuntimeCost {
    type Output = RuntimeCost;

    fn add(self, rhs: Duration) -> Self::Output {
        self + RuntimeCost::from(rhs)
    }
}

impl AddAssign for RuntimeCost {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl AddAssign<Duration> for RuntimeCost {
    fn add_assign(&mut self, rhs: Duration) {
        *self = *self + rhs
    }
}


pub fn estimate_runtime_cost(g: &CostFactors) -> RuntimeCost {
    let nvtxs = g.nvtxs as f64;
    let nedges = g.nedges as f64;

    let avg_degree = nedges / nvtxs;

    // some quick-and-dirty checking suggests that minimizing the average degree is very impactful
    // to runtime (4elt vs web-spam)
    let deg_factor = avg_degree.powi(2);

    // number of vertices itself doesn't seem to have all that much impact, it's primarily on the
    // number of edges. This is a starting point bases on 4elt
    let size_factor = nedges as f64 / 630.0;

    // partition number impact is non-linear, but we'll pretend it isn't
    let nparts_factor = g.nparts as f64;

    // contig is pretty variable - sometimes it's zero cost, but it can get expensive
    let contig_factor = if g.contig && g.ptype.is_kway() { 2.0 } else { 1.0 };

    // estimate based on 4elt and web-spam
    let niter_factor = 1.0 + 0.1 * g.niter as f64;

    // shot-in-the-dark guess
    let vsize_factor = 1.0 + g.vsize.is_some() as u8 as f64;

    // shot-in-the-dark guess
    let minconn_factor = 1.0 + g.minconn as u8 as f64;

    // shot-in-the-dark guess
    let vwgt_factor = g.ncon.max(1).pow(2) as f64;

    // shot-in-the-dark guess
    let adjwgt_factor = if g.adjwgt.is_some() { 1.5 } else { 1.0 };

    // there is zero rhyme nor reason to this calculation
    let vsize_dist_factor = g.vsize.map_or(1.0, |x| x.score());
    let vwgt_dist_factor = g.vwgt.map_or(1.0, |x| x.score());
    let adjwgt_dist_factor = g.adjwgt.map_or(1.0, |x| x.score());
    let dist_factor = (vsize_dist_factor * vwgt_dist_factor * adjwgt_dist_factor + 2.0).ln() as f64;

    let total = dist_factor * adjwgt_factor * minconn_factor * deg_factor * size_factor * nparts_factor * contig_factor * niter_factor * vsize_factor * vwgt_factor;

    assert!(total.is_finite());
    assert!(total > 0.0);

    let total = (total.round() as usize).max(5);

    RuntimeCost(total)
}

pub fn estimate_human_cost(g: &CostFactors) -> HumanCost {
    // generally nice to minimize
    let size_factor = (g.nvtxs * g.nedges) as f64 / 100.0;

    // number of partitions generally isn't a big impediment, but 2 is good
    let nparts_factor = if g.nparts == 2 { 0.5 } else { (g.nparts as f64).ln_1p() * 2.0 };

    // contig doesn't change all too much, but it's nicer to not have it
    let contig_factor = if g.contig && g.ptype.is_kway() { 1.5 } else { 1.0 };

    // it's fine if niter is small, but if it has to be large that kinda sucks
    let niter_factor = if g.niter < 3 { 1.0 } else { 5.0 };

    // not a huge deal, but would rather not have it
    let minconn_factor = 1.0 + g.minconn as u8 as f64;

    // we would very much rather not have multiconstraint
    let ncon_factor = if g.ncon == 0 { 0.75 } else if g.ncon == 1 { 1.0 } else { 5.0 };

    let vsize_factor = g.vsize.map_or(1.0, |x| x.score()) as f64;
    let vwgt_factor = g.vwgt.map_or(1.0, |x| x.score()) as f64;
    let adjwgt_factor = g.adjwgt.map_or(1.0, |x| x.score()) as f64;

    let total = ncon_factor * adjwgt_factor * minconn_factor * size_factor * nparts_factor * contig_factor * niter_factor * vsize_factor * vwgt_factor;

    let total = (total.round() as usize).max(0);

    HumanCost(total)
}

#[cfg(test)]
mod tests {
    use fastrand::Rng;

    use super::*;

    /// dist cost assert less than
    #[track_caller]
    fn dc_alt(l: &[usize], r: &[usize]) {
        let lc = DistCost::new(l);
        let rc = DistCost::new(r);
        let lcs = lc.score();
        let rcs = rc.score();
        // println!("lhs cost: {lc:?} ({lcs})\nrhs cost: {rc:?} ({rcs})\n");
        assert!(lcs < rcs, "lhs cost: {lc:?} ({lcs})\nrhs cost: {rc:?} ({rcs})");
    }
    
    #[test]
    fn dist_cost_sanity_check() {
        dc_alt(&[1, 2, 3], &[23, 91, 1]);
        dc_alt(&[1, 1, 1], &[2, 2, 2]);
        dc_alt(&[1, 1, 1], &[1, 1, 1, 1, 1]);
    }

    #[test]
    fn dist_cost_decrease_entails_decrease() {
        let mut rng = Rng::with_seed(1234);
        let mut vl = |len| {
            let mut rng = rng.fork();
            std::iter::repeat_with(move || rng.usize(..50)).take(len).collect::<Vec<_>>()
        };
        let mut case = |len| {
            let og = vl(len);
            let reduce = og.iter().map(|&x| x.div_ceil(2)).collect::<Vec<_>>();
            dc_alt(&reduce, &og)
        };
        case(3);
        case(10);
        case(20);
        case(20);
        case(20);
        case(20);
    }
}
