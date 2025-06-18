use crate::costs::CaseCost;

use super::{Strategy, StrategyIter};

pub struct RemoveVsize;
impl Strategy for RemoveVsize {
    fn dont_backtrack(&self) -> bool { true }

    fn cost(&self, case: &super::Case) -> CaseCost {
        let mut f = case.cost_factors();
        f.vsize = false;
        f.into()
    }

    fn apply(&self, mut case: super::Case) -> Option<impl StrategyIter> {
        if case.graph.vsize.is_none() {
            return None
        }
        case.graph.vsize = None;
        Some(std::iter::once(case))
    }
}

pub struct RemoveAllVWgt;
impl Strategy for RemoveAllVWgt {
    fn dont_backtrack(&self) -> bool { true }

    fn cost(&self, case: &super::Case) -> CaseCost {
        let mut f = case.cost_factors();
        f.ncon = 0;
        f.into()
    }

    fn apply(&self, mut case: super::Case) -> Option<impl StrategyIter> {
        if case.graph.vwgt.is_none() {
            return None
        }
        case.graph.vwgt = None;
        case.graph.ncon = 0;
        Some(std::iter::once(case))
    }
}

pub struct UnsetMinConn;
impl Strategy for UnsetMinConn {
    fn dont_backtrack(&self) -> bool { true }

    fn cost(&self, case: &super::Case) -> CaseCost {
        let mut f = case.cost_factors();
        f.minconn = false;
        f.into()
    }

    fn apply(&self, mut case: super::Case) -> Option<impl StrategyIter> {
        if !case.minconn {
            return None
        }
        case.minconn = false;
        Some(std::iter::once(case))
    }
}

pub struct UnsetContig;
impl Strategy for UnsetContig {
    fn dont_backtrack(&self) -> bool { true }

    fn cost(&self, case: &super::Case) -> CaseCost {
        let mut f = case.cost_factors();
        f.contig = false;
        f.into()
    }

    fn apply(&self, mut case: super::Case) -> Option<impl StrategyIter> {
        if !case.contig {
            return None
        }
        case.contig = false;
        Some(std::iter::once(case))
    }
}
