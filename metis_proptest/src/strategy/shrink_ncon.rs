use std::sync::Arc;

use fastrand::Rng;

use crate::{costs::CaseCost, graph::Graph, utils::RetainIndexed};

use super::{Case, Strategy, StrategyKind};



pub struct ShrinkNcon;

impl Strategy for ShrinkNcon {
    fn is_valid(&self, case: &Case) -> bool {
        case.graph.ncon > 0
    }

    fn kind(&self) -> StrategyKind {
        StrategyKind::Checkpoint
    }

    fn cost(&self, case: &Case) -> CaseCost {
        let mut f = case.cost_factors();
        match f.ncon {
            0 => (),
            1 => {
                f.vwgt = None;
                f.ncon = 0;
            }
            2 => {
                f.vwgt = Some(f.vwgt.expect("ncon ne zero implies vwgt").scale_count(0.5));
                f.ncon = 0;
            }
            3.. => {
                let factor = 2.0 / f.ncon as f64;
                f.vwgt = Some(f.vwgt.expect("ncon ne zero implies vwgt").scale_count(factor));
                f.ncon = 2;
            }
        }
        f.into()
    }

    fn apply(&self, case: Case) -> Option<impl super::StrategyIter> {
        if !self.is_valid(&case) {
            return None
        }

        Some((0..case.graph.ncon).map(move |new_ncon| {
            let mut case = case.clone();
            decrease_ncon(&mut case.graph, new_ncon);
            case
        }))
    }
}

fn decrease_ncon(graph: &mut Graph, new_ncon: usize) {
    debug_assert!(graph.ncon > new_ncon);

    if new_ncon == 0 {
        graph.ncon = 0;
        graph.vwgt = None;
        return
    }

    let vwgt = graph.vwgt.as_mut().map(Arc::make_mut).expect("non-zero ncon should imply vwgt");

    vwgt.retain_indexed(|i, _| i % graph.ncon < new_ncon);

    graph.ncon = new_ncon;
}

pub struct ReduceVwgt {
    pub amt_percent: u8,
}

impl Strategy for ReduceVwgt {
    fn is_valid(&self, case: &Case) -> bool {
        case.graph.vwgt.as_ref().is_some_and(|vwgt| vwgt.iter().any(|&v| 1 < v))
    }

    fn kind(&self) -> StrategyKind {
        StrategyKind::Repeat
    }

    fn cost(&self, case: &Case) -> CaseCost {
        let mut f = case.cost_factors();
        f.vwgt = f.vwgt.map(|x| x.scale_values(0.5));
        f.into()
    }

    fn apply(&self, case: Case) -> Option<impl super::StrategyIter> {
        let vwgt = case.graph.vwgt.clone()?;
        let mut rng = Rng::new();
        let percent = self.amt_percent.min(100);
        Some(std::iter::repeat_with(move || {
            let mut case = case.clone();
            let mut rng = rng.fork();
            let mut vwgt = vwgt.clone();
            for wgt in Arc::make_mut(&mut vwgt) {
                if rng.u8(0..100) < percent {
                    *wgt = wgt.div_ceil(2)
                }
            }
            case.graph.vwgt = Some(vwgt);
            case
        }))
    }
}

pub struct ReduceAdjwgt {
    pub amt_percent: u8,
}

impl Strategy for ReduceAdjwgt {
    fn is_valid(&self, case: &Case) -> bool {
        case.graph.adjwgt.as_ref().is_some_and(|adjwgt| adjwgt.iter().any(|&v| 1 < v))
    }

    fn kind(&self) -> StrategyKind {
        StrategyKind::Repeat
    }

    fn cost(&self, case: &Case) -> CaseCost {
        let mut f = case.cost_factors();
        f.adjwgt = f.adjwgt.map(|x| x.scale_values(0.5));
        f.into()
    }

    fn apply(&self, case: Case) -> Option<impl super::StrategyIter> {
        let mut rng = Rng::new();
        let percent = self.amt_percent.min(100);

        let mut pairs = vec![usize::MAX; case.graph.adjncy.len()];

        {
            let g = &case.graph;
            for (from, tos) in g.edges_iter().enumerate() {
                for (from_off, &to) in tos.iter().enumerate() {
                    debug_assert!(to != from, "graph contains a self-loop");

                    let Some(to_off) = g.xadj[to..(to+1)].iter().position(|&v| v == from) else {
                        continue
                        // eprintln!("{from} -> {to} doesn't have a counterpart ({to} -> {from} is missing)",
                        // from = from + 1, to = to + 1);
                        // g.write(std::io::stderr()).unwrap();
                        // panic!("graph appears to be directed")
                    };

                    let from_idx = g.xadj[from] + from_off;
                    let to_idx = g.xadj[to] + to_off;

                    pairs[from_idx] = to_idx;
                    pairs[to_idx] = from_idx;
                }
            }
        }
        Some(std::iter::repeat_with(move || {
            let mut case = case.clone();
            let mut rng = rng.fork();
            let nedges = case.graph.nedges();
            let adjwgt = Arc::make_mut(case.graph.adjwgt.as_mut().unwrap());
            let pick = (nedges * percent as usize).div_ceil(100);
            let targets = rng.choose_multiple(0..case.graph.adjncy.len(), pick);

            for target in targets {
                adjwgt[target] = adjwgt[target].div_ceil(2);
                if pairs[target] != usize::MAX {
                    adjwgt[pairs[target]] = adjwgt[pairs[target]].div_ceil(2);
                }
            }

            case
        }))
    }
}
