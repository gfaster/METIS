use crate::strategy::DynStrategyIter;
use crate::strategy::DynStrategy;
use crate::strategy::Strategy;
use crate::Case;
use std::collections::BinaryHeap;
use std::sync::Arc;
use std::time::Duration;
use std::time::Instant;
use crate::costs::RuntimeCost;


// have this as a separate struct so I can tweak ordering behavior
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Priority(usize);

impl Priority {
    fn decrease(self) -> Self {
        if self.0 > (usize::MAX >> usize::BITS/2) {
            Self(self.0 / 256)
        } else if self.0 > 100_000_000 {
            Self(self.0 / 8)
        } else if self.0 > 1000 {
            Self(self.0 / 2)
        } else if self.0 > 100 {
            Self(self.0 - 25)
        } else {
            Self(self.0.saturating_sub(10))
        }
    }

    fn is_dead(self) -> bool {
        self.0 == 0
    }

    fn from_est_cost(c: RuntimeCost, existing: Option<Priority>) -> Self {
        let existing = existing.unwrap_or(Priority(100)).0;
        let ret = usize::MAX.saturating_sub(c.to_raw().saturating_mul(100_000).saturating_mul(existing));
        Priority(ret)
    }

    fn zero() -> Self {
        Priority(0)
    }
}

#[derive(Clone)]
struct MinSetPriority<'a> {
    priority: Priority,

    /// index into the chain vector for propagating information to parents
    chain_idx: u32,

    /// index into strats that no further strategies should go before
    backtrack_idx: u32,

    iter: DynStrategyIter<'a>,
}

impl PartialEq for MinSetPriority<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.priority == other.priority
    }
}

impl Eq for MinSetPriority<'_>{ }

impl PartialOrd for MinSetPriority<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.priority.partial_cmp(&other.priority)
    }
}

impl Ord for MinSetPriority<'_> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.priority.cmp(&other.priority)
    }
}

struct InfoChain {
    up: Option<u32>,
    runs: u32,
    crashes: u32,
    total_time: RuntimeCost,
    best: Option<Arc<Case>>,
    best_cost: Option<RuntimeCost>,
}

pub struct MinimizationSet {
    strats: Vec<DynStrategy<'static>>
}

impl MinimizationSet {
    pub fn new() -> Self {
        Self {
            strats: vec![],
        }
    }

    pub fn push(&mut self, i: impl Strategy + 'static) {
        self.strats.push(DynStrategy::new_owned(i))
    }

    pub fn work(&self, case: Case) -> MinimizationSetWorker {
        let acase = Arc::new(case.clone());
        let info = InfoChain {
            up: None,
            runs: 0,
            crashes: 0,
            total_time: RuntimeCost::ZERO,
            best: None,
            best_cost: None,
        };
        let first = MinSetPriority {
            priority: Priority::zero(),
            chain_idx: 0,
            backtrack_idx: 0,
            iter: DynStrategyIter::new(std::iter::once(case)),
        };
        MinimizationSetWorker { 
            set: self,
            info_chain: vec![info],
            queue: BinaryHeap::from([first]),
            best_rt_cost: acase.cost().runtime,
            best_rt: acase,
        }
    }
}

pub struct MinimizationSetWorker<'a> {
    set: &'a MinimizationSet,
    info_chain: Vec<InfoChain>,
    queue: BinaryHeap<MinSetPriority<'a>>,
    best_rt: Arc<Case>,
    best_rt_cost: RuntimeCost,
}

impl MinimizationSetWorker<'_> {
    pub fn run_next(&mut self) -> Option<(Case, bool)> {
        let (mut set, case) = loop {
            let mut it = self.queue.pop()?;

            let Some(case) = it.iter.next() else { continue };
            break (it, case);
        };
        let case_ret = case.clone();
        let case = Arc::new(case);

        eprintln!("Running with priority: {:?} and est cost {}", set.priority, case.cost().runtime);

        let timer = Instant::now();
        let res = case.run_program();
        let cost = RuntimeCost::from(timer.elapsed());
        let crashed = res.is_crash();

        if cost < self.best_rt_cost && crashed {
            self.best_rt_cost = cost;
            self.best_rt = Arc::clone(&case);
        }


        // update info chain
        let mut best_depth = 0;
        let mut total_depth = 0;
        {
            let mut up = Some(set.chain_idx);
            while let Some(i) = up {
                let info = &mut self.info_chain[i as usize];
                info.crashes += crashed as u32;
                info.runs += 1;
                info.total_time += cost;
                if let Some(best_cost) = info.best_cost.as_mut() && cost < *best_cost && crashed {
                    *best_cost = cost;
                    info.best = Some(Arc::clone(&case));
                    best_depth += 1;
                } else if info.best_cost.is_none() && crashed {
                    info.best_cost = Some(cost);
                    info.best = Some(Arc::clone(&case));
                    best_depth += 1;
                }
                total_depth += 1;
                up = info.up;
            }
        }

        if best_depth == 0 {
            // either didn't reproduce crash or wasn't better
            set.priority = set.priority.decrease().decrease();
            if !set.priority.is_dead() {
                self.queue.push(set);
            }
        } else if total_depth <= best_depth + 2 {
            // pretty deep and pretty good, let's explore
            for (i, strat) in self.set.strats[set.backtrack_idx as usize..].iter().enumerate() {
                if !strat.is_valid(&case) { continue }
                let Some(iter) = strat.apply((*case).clone()) else { continue };
                let priority = Priority::from_est_cost(strat.cost(&case).runtime, Some(set.priority));
                let chain_idx = self.info_chain.len() as u32;
                self.info_chain.push(InfoChain { 
                    up: Some(set.chain_idx),
                    runs: 0,
                    crashes: 0,
                    total_time: RuntimeCost::ZERO,
                    best: None,
                    best_cost: None
                });
                let new_set = MinSetPriority {
                    priority,
                    chain_idx,
                    backtrack_idx: if strat.dont_backtrack() { i as u32 } else { set.backtrack_idx },
                    iter,
                };
                self.queue.push(new_set);
            }

            set.priority = if total_depth <= best_depth + 1 {
                set.priority
            } else {
                set.priority.decrease()
            };

            if !set.priority.is_dead() {
                self.queue.push(set);
            }
        } else {
            set.priority = set.priority.decrease();
            if !set.priority.is_dead() {
                self.queue.push(set);
            }
        }

        Some((case_ret, crashed))
    }

    pub fn run_for(&mut self, duration: Duration) -> bool {
        let start = Instant::now();
        let mut improved = false;
        loop {
            let Some((_, better)) = self.run_next() else { break };
            improved |= better;
            if duration < start.elapsed() {
                break
            }
        }
        improved
    }

    pub fn best(&self) -> Arc<Case> {
        self.best_rt.clone()
    }
}
