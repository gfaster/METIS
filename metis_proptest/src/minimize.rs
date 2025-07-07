use crate::costs::CaseCost;
use crate::costs::HumanCost;
use crate::strategy::DynStrategyIter;
use crate::strategy::DynStrategy;
use crate::strategy::Strategy;
use crate::strategy::StrategyKind;
use crate::Case;
use std::collections::BinaryHeap;
use std::fmt;
use std::iter::Peekable;
use std::sync::Arc;
use std::time::Duration;
use std::time::Instant;
use crate::costs::RuntimeCost;


// have this as a separate struct so I can tweak ordering behavior
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Priority(u32);


impl Priority {
    #[must_use]
    fn adjust(self, n: i32) -> Self {
        match n {
            ..0 => self.dec_n(n.abs() as u32),
            0 => self,
            1.. => self.inc_n(n as u32)
        }
    }

    #[must_use]
    fn dec_n(mut self, n: u32) -> Self {
        for _ in 0..n {
            self = self.decrease()
        }
        self
    }

    #[must_use]
    fn decrease(self) -> Self {
        // if self.0 > 10_000_000 {
        //     Self(self.0 / 16)
        // } else if self.0 > 1000 {
        //     Self(self.0 / 4)
        // } else if self.0 > 100 {
        //     Self(self.0 - 25)
        // } else {
        //     Self(self.0.saturating_sub(5))
        // }
        debug_assert!(self <= Self::max());
        let ret = if self.0 * 10 > Self::max().0 {
            self.0 / 2
        } else if self.0 > 11 {
            self.0 * 4 / 5
        } else {
            self.0.saturating_sub(1)
        };
        Self(ret)
    }

    #[must_use]
    fn inc_n(mut self, n: u32) -> Self {
        for _ in 0..n {
            self = self.increase()
        }
        self
    }

    #[must_use]
    fn increase(self) -> Self {
        debug_assert!(self <= Self::max());
        let ret = self.0.max(10).saturating_mul(10);
        Priority(ret).min(Self::max())
    }

    #[must_use]
    fn is_dead(self) -> bool {
        self == Self::zero()
    }

    #[must_use]
    fn from_est_cost(c: RuntimeCost, existing: Option<Priority>) -> Self {
        let existing = existing.unwrap_or_default();
        // let ret = usize::MAX.saturating_sub(c.to_raw().saturating_mul(100_000).saturating_mul(existing));
        // let ret = existing.decrease().decrease().0.max(1).div_ceil(c.to_raw());
        let ret = if c.to_raw() < 500 {
            (existing.0 as usize).saturating_sub(1)
        } else {
            (existing.0 as usize).saturating_sub(c.to_raw())
        };
        Priority(ret.clamp(1, Self::max().0 as usize) as u32)
    }

    #[must_use]
    const fn zero() -> Self {
        Priority(0)
    }

    #[must_use]
    const fn max() -> Self {
        Priority(100_000)
    }
}

impl Default for Priority {
    fn default() -> Self {
        Priority(100)
    }
}

impl fmt::Debug for Priority {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if *self >= Self::max() {
            fmt::Display::fmt(&"MAX", f)
        } else {
            self.0.fmt(f)
        }
    }
}

#[derive(Clone)]
struct MinSetPriority<'a> {
    priority: Priority,

    /// index into the chain vector for propagating information to parents
    chain_idx: u32,

    /// index of this strategy
    strat_idx: u32,

    kind: StrategyKind,

    /// the initial case if `kind` is [`StrategyKind::Checkpoint`]
    checkpoint_first: Option<Arc<Case>>,

    iter: Box<Peekable<DynStrategyIter<'a>>>,
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
    best_rt_cost: RuntimeCost,
    best_hu_cost: HumanCost,
    num_improvements: u32,
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

    pub fn work(&self, case: Case) -> MinimizationSetWorker<'_> {
        let acase = Arc::new(case.clone());
        let CaseCost {
            runtime: best_rt_cost,
            human: best_hu_cost,
        } = acase.cost();

        let info = InfoChain {
            up: None,
            runs: 0,
            crashes: 0,
            total_time: RuntimeCost::ZERO,
            best_rt_cost,
            best_hu_cost,
            num_improvements: 0,
        };
        let first = MinSetPriority {
            priority: Priority::max(),
            chain_idx: 0,
            strat_idx: 0,
            kind: StrategyKind::MustAccept,
            checkpoint_first: None,
            iter: Box::new(DynStrategyIter::new(std::iter::once(case)).peekable()),
        };
        MinimizationSetWorker { 
            set: self,
            info_chain: vec![info],
            queue: BinaryHeap::from([first]),
            best_rt_cost,
            best_rt: Arc::clone(&acase),
            best_hu_cost,
            best_hu: acase,
        }
    }
}

pub struct MinimizationSetWorker<'a> {
    set: &'a MinimizationSet,
    info_chain: Vec<InfoChain>,
    queue: BinaryHeap<MinSetPriority<'a>>,
    best_rt: Arc<Case>,
    best_rt_cost: RuntimeCost,
    best_hu: Arc<Case>,
    best_hu_cost: HumanCost,
}

impl MinimizationSetWorker<'_> {
    // returns true if `strat.apply` returned `Some` (caller needs special handling for if this
    // returns false and kind is Checkpoint)
    fn enqueue(&mut self, case: Arc<Case>, up: u32, strat_idx: u32, kind: StrategyKind, priority: Priority) -> bool {
        assert!((up as usize) < self.info_chain.len());
        if priority.is_dead() {
            return false
        }
        let strat = &self.set.strats[strat_idx as usize];
        let Some(iter) = strat.apply((*case).clone()) else { return false };
        let checkpoint_first = kind.is_checkpoint().then_some(case);
        let chain_idx = self.info_chain.len() as u32;
        self.info_chain.push(InfoChain { 
            up: Some(up),
            runs: 0,
            crashes: 0,
            total_time: RuntimeCost::ZERO,
            best_rt_cost: self.info_chain[up as usize].best_rt_cost,
            best_hu_cost: self.info_chain[up as usize].best_hu_cost,
            num_improvements: 0,
        });
        let new_set = MinSetPriority {
            priority,
            chain_idx,
            strat_idx,
            iter: Box::new(iter.peekable()),
            kind,
            checkpoint_first,
        };
        self.queue.push(new_set);
        true
    }

    fn enqueue_next_valid(&mut self, case: Arc<Case>, up: u32, mut strat_idx: u32, priority: Priority) {
        assert!((up as usize) < self.info_chain.len());
        if priority.is_dead() {
            return
        }
        let (iter, kind) = loop { 
            strat_idx += 1;
            let Some(strat) = self.set.strats.get(strat_idx as usize - 1) else { return };
            if !strat.is_valid(&case) { continue }
            let Some(iter) = strat.apply((*case).clone()) else { continue };
            strat_idx -= 1;
            break (iter, strat.kind())
        };
        let checkpoint_first = kind.is_checkpoint().then_some(case);
        let chain_idx = self.info_chain.len() as u32;
        self.info_chain.push(InfoChain { 
            up: Some(up),
            runs: 0,
            crashes: 0,
            total_time: RuntimeCost::ZERO,
            best_rt_cost: self.info_chain[up as usize].best_rt_cost,
            best_hu_cost: self.info_chain[up as usize].best_hu_cost,
            num_improvements: 0,
        });
        let new_set = MinSetPriority {
            priority,
            chain_idx,
            strat_idx,
            iter: Box::new(iter.peekable()),
            kind,
            checkpoint_first,
        };
        self.queue.push(new_set);
    }

    pub fn run_next(&mut self) -> Option<(Case, bool)> {
        let (mut set, case) = loop {
            let mut it = self.queue.pop()?;

            let Some(case) = it.iter.next() else { continue };
            break (it, case);
        };
        let case_ret = case.clone();
        let case = Arc::new(case);
        let case_est_cost = case.cost();
        let cost_hu = case_est_cost.human;

        eprintln!("Running strat {:2} ({:15?}) with priority: {:?} and est cost {}", set.strat_idx, set.kind, set.priority, case_est_cost.runtime);

        let timer = Instant::now();
        let res = case.run_program();
        let cost_rt = RuntimeCost::from(timer.elapsed());
        let crashed = res.is_crash();
        // if crashed {
        //     eprintln!("\tcrashed!")
        // }
        let improved_rt = cost_rt < self.best_rt_cost && crashed;
        if improved_rt {
            self.best_rt_cost = cost_rt;
            self.best_rt = Arc::clone(&case);
            eprintln!("\tnew best rt")
        }
        let improved_human = cost_hu < self.best_hu_cost && crashed;
        if improved_human {
            self.best_hu_cost = cost_hu;
            self.best_hu = Arc::clone(&case);
            eprintln!("\tnew best human")
        }

        let improvement = improved_rt || improved_human;


        // update info chain
        let mut best_rt_depth = 0;
        let mut best_hu_depth = 0;
        let mut total_depth = 0;
        {
            let mut up = Some(set.chain_idx);
            while let Some(i) = up {
                let info = &mut self.info_chain[i as usize];
                info.crashes += crashed as u32;
                info.runs += 1;
                info.total_time += cost_rt;
                if cost_rt < info.best_rt_cost && crashed {
                    info.best_rt_cost = cost_rt;
                    best_rt_depth += 1;
                }
                if cost_hu < info.best_hu_cost && crashed {
                    info.best_hu_cost = cost_hu;
                    best_hu_depth += 1;
                }
                total_depth += 1;
                up = info.up;
            }
        }
        let best_depth = best_hu_depth.max(best_rt_depth);

        let is_fruitless_search = {
            let info = &self.info_chain[set.chain_idx as usize];
            info.runs > 20 
            && best_depth != total_depth
            && info.num_improvements == 0
        };

        match set.kind {
            StrategyKind::Checkpoint | StrategyKind::MustAccept => {
                if crashed || set.iter.peek().is_none() {
                    // pass on if there's no more or if it was accepted
                    self.queue.clear();

                    if !crashed && set.kind.is_must_accept() {
                        eprintln!("MUST_ACCEPT was not accepted");
                        return None
                    }

                    if crashed {
                        eprintln!("\tcheckpoint pass");
                    } else {
                        eprintln!("\tcheckpoint (fail)");
                    }

                    let next_case = if crashed { case } else { set.checkpoint_first.expect("checkpoint_first is set on checkpoint") };
                    self.enqueue_next_valid(next_case, set.chain_idx, set.strat_idx + 1, Priority::max());
                } else {
                    self.queue.push(set);
                }
            },
            StrategyKind::Ordered | StrategyKind::First | StrategyKind::Repeat if !crashed => {
                set.priority = set.priority.decrease();

                if !set.priority.is_dead() {
                    self.queue.push(set);
                }
            },
            StrategyKind::Ordered | StrategyKind::Repeat if !improved_human => {
                set.priority = set.priority.decrease();

                if !set.priority.is_dead() {
                    self.queue.push(set);
                }
            },
            StrategyKind::First => {
                let child_adj;
                if best_depth == 0 {
                    unreachable!("Strategy with kind First should have non-zero best depth");
                } else if best_depth == total_depth {
                    child_adj = 2;
                } else {
                    child_adj = -1;
                }

                let child_priority = Priority::from_est_cost(cost_rt, Some(set.priority)).adjust(child_adj);
                if !child_priority.is_dead() {
                    self.enqueue_next_valid(case, set.chain_idx, set.strat_idx + 1, child_priority);
                }
            }
            StrategyKind::Ordered => {
                let child_adj;
                let next_adj;
                if best_depth == 0 {
                    child_adj = -2;
                    next_adj = -1;
                } else if best_depth == total_depth {
                    child_adj = 3;
                    next_adj = 3;
                } else {
                    child_adj = -1;
                    next_adj = -1;
                }

                let child_priority = Priority::from_est_cost(cost_rt, Some(set.priority)).adjust(child_adj);
                set.priority = set.priority.adjust(next_adj);

                self.enqueue_next_valid(case, set.chain_idx, set.strat_idx + 1, child_priority);

                if !set.priority.is_dead() {
                    self.queue.push(set);
                }
            }
            StrategyKind::Repeat => {
                let strat = &self.set.strats[set.strat_idx as usize];
                let adj = if improvement { 10 } else { -3 };

                let priority = set.priority.adjust(adj);
                let do_recurse = strat.is_valid(&case) && !is_fruitless_search;
                if do_recurse {
                    self.enqueue(case.clone(), set.chain_idx, set.strat_idx, set.kind, priority);
                }
                self.enqueue_next_valid(case, set.chain_idx, set.strat_idx + 1, priority);

                set.priority = set.priority.adjust(adj - 3);

                if !set.priority.is_dead() {
                    self.queue.push(set);
                }
            },
        }

        Some((case_ret, crashed))
    }

    pub fn run_for(&mut self, duration: Duration) -> Option<bool> {
        let start = Instant::now();
        let initial_best_hu = self.best_hu_cost;
        let initial_best_rt = self.best_rt_cost;
        let mut ran = false;
        loop {
            let Some((_, _)) = self.run_next() else { break };
            ran = true;
            if duration < start.elapsed() {
                break
            }
        }
        ran.then_some(self.best_hu_cost < initial_best_hu || self.best_rt_cost < initial_best_rt)
    }

    pub fn best_rt(&self) -> Arc<Case> {
        self.best_rt.clone()
    }

    pub fn best_hu(&self) -> Arc<Case> {
        self.best_hu.clone()
    }
}
