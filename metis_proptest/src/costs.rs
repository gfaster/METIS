//! Centralized (standardized) cost estimation and measurement

use std::{fmt, ops::{Add, AddAssign}, time::Duration};

use crate::{strategy::Case, PartScheme};


// Units are "something like" milliseconds
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct RuntimeCost(usize);

impl RuntimeCost {
    pub const ZERO: Self = RuntimeCost(0);

    pub fn from_raw(u: usize) -> Self {
        Self(u)
    }

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

#[derive(Debug, Clone, Copy)]
pub struct CostFactors {
    pub nvtxs: u32,
    pub nedges: u32,
    pub nparts: u32,
    pub contig: bool,
    pub minconn: bool,
    pub niter: u32,
    pub vsize: bool,
    pub ptype: PartScheme,
    pub ncon: u32,
}

impl CostFactors {
    pub fn from_case(case: &Case) -> Self {
        CostFactors { 
            nvtxs: case.graph.nvtxs() as u32,
            nedges: case.graph.nedges() as u32,
            nparts: case.nparts as u32,
            contig: case.contig,
            niter: case.niter as u32,
            vsize: case.graph.vsize.is_some(),
            ptype: case.ptype,
            ncon: case.graph.ncon as u32,
            minconn: case.minconn,
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
    let vsize_factor = 1.0 + g.vsize as u8 as f64;

    // shot-in-the-dark guess
    let minconn_factor = 1.0 + g.minconn as u8 as f64;

    // shot-in-the-dark guess
    let vwgt_factor = g.ncon.max(1).pow(2) as f64;

    let total = minconn_factor * deg_factor * size_factor * nparts_factor * contig_factor * niter_factor * vsize_factor * vwgt_factor;

    debug_assert!(total > 0.0);

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
    let vsize_factor = 1.0 + g.vsize as u8 as f64;

    // not a huge deal, but would rather not have it
    let minconn_factor = 1.0 + g.minconn as u8 as f64;

    // we would very much rather not have multiconstraint
    let vwgt_factor = if g.ncon == 0 { 0.75 } else if g.ncon == 1 { 1.0 } else { 5.0 };

    let total = minconn_factor * size_factor * nparts_factor * contig_factor * niter_factor * vsize_factor * vwgt_factor;

    let total = (total.round() as usize).max(0);

    HumanCost(total)
}
