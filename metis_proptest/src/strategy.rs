use std::{ffi::OsStr, os::unix::process::ExitStatusExt, sync::Arc};

use crate::{costs::{CaseCost, CostFactors}, graph::Graph, Ctype, Iptype, Objtype, PartScheme};

pub mod simple;
pub mod combinator;
pub use combinator::*;
pub mod delete_vtx;
pub mod shrink_ncon;

pub mod strategies {
    #![allow(unused_imports)]

    pub use super::simple::*;
    pub use super::combinator::*;
    pub use super::delete_vtx::*;
    pub use super::shrink_ncon::*;
}

/// program-wide settings
pub struct Settings {
    pub comm: Arc<OsStr>,
    pub nocapture: bool,
}

/// the test case we're trying to reduce
#[derive(Clone)]
pub struct Case {
    pub graph: Graph,
    pub minconn: bool,
    pub contig: bool,
    pub niter: usize,
    pub ncuts: usize,
    pub nparts: usize,
    pub ufactor: u32,
    pub ptype: PartScheme,
    pub iptype: Iptype,
    pub objtype: Objtype,
    pub ctype: Ctype,
    pub seed: i32,
    pub settings: Arc<Settings>,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum StrategyKind {
    /// Provides a major difference in kind that is meaningless to attempt multiple times or to go
    /// back once an accepted case is found. When a checkpoint is hit for the first time, the
    /// first accepted case produced by the strategy iterator (or the original case if none are)
    /// is passed on. Regardless of whether the strategy produces an accepted case, the queue is
    /// cleared and this strategy will not be run again.
    Checkpoint,

    /// Like [`StrategyKind::Checkpoint`], but we abort the whole run if nothing is accepted.
    /// Intended for the initial case.
    MustAccept,

    /// passes each accepted case forward
    Ordered,

    /// like [`StrategyKind::Ordered`], but it only passes the first accepted
    #[expect(dead_code)]
    First,

    /// repeatadly applies until [`Strategy::is_valid`] returns `false`. If the iterator returns
    /// multiple elements, this may result in an exponential explosion of cases.
    Repeat,

}

#[allow(dead_code)]
impl StrategyKind {
    /// Returns `true` if the strategy kind is [`Checkpoint`].
    ///
    /// [`Checkpoint`]: StrategyKind::Checkpoint
    #[must_use]
    pub fn is_checkpoint(&self) -> bool {
        matches!(self, Self::Checkpoint)
    }

    /// Returns `true` if the strategy kind is [`MustAccept`].
    ///
    /// [`MustAccept`]: StrategyKind::MustAccept
    #[must_use]
    pub fn is_must_accept(&self) -> bool {
        matches!(self, Self::MustAccept)
    }

    /// Returns `true` if the strategy kind is [`Repeat`].
    ///
    /// [`Repeat`]: StrategyKind::Repeat
    #[must_use]
    pub fn is_repeat(&self) -> bool {
        matches!(self, Self::Repeat)
    }

    /// Returns `true` if the strategy kind is [`Ordered`].
    ///
    /// [`Ordered`]: StrategyKind::Ordered
    #[must_use]
    pub fn is_ordered(&self) -> bool {
        matches!(self, Self::Ordered)
    }
}

pub trait Strategy {
    /// what kind of strategy is this? See the documentation for [`StrategyKind`] for more details.
    fn kind(&self) -> StrategyKind { StrategyKind::Ordered }

    fn cost(&self, case: &Case) -> CaseCost;

    fn is_valid(&self, case: &Case) -> bool {
        let _ = case;
        true
    }

    /// the intent here is to create an iterator of different applications of the strategy
    fn apply(&self, case: Case) -> Option<impl StrategyIter>;
}

pub trait StrategyIter: Iterator<Item = Case> + Clone
where
{ }

impl<T: Iterator<Item = Case> + Clone> StrategyIter for T
where
{}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ProgramResult {
    NormalExit,
    GracefulFailure,
    Crash,
}

impl ProgramResult {
    /// Returns `true` if the program result is [`GracefulFailure`].
    ///
    /// [`GracefulFailure`]: ProgramResult::GracefulFailure
    #[must_use]
    pub fn is_graceful_failure(&self) -> bool {
        matches!(self, Self::GracefulFailure)
    }

    /// Returns `true` if the program result is [`NormalExit`].
    ///
    /// [`NormalExit`]: ProgramResult::NormalExit
    #[must_use]
    pub fn is_normal_exit(&self) -> bool {
        matches!(self, Self::NormalExit)
    }

    /// Returns `true` if the program result is [`Crash`].
    ///
    /// [`Crash`]: ProgramResult::Crash
    #[must_use]
    pub fn is_crash(&self) -> bool {
        matches!(self, Self::Crash)
    }
}

impl Case {
    pub fn run_program(&self) -> ProgramResult {
        use std::process::{Command, Stdio};

        let mut cmd = Command::new(&self.settings.comm);
        cmd.arg("/proc/self/fd/0");
        cmd.arg(format!("{}", self.nparts));

        cmd.arg(match self.ptype {
            PartScheme::Kway => "-ptype=kway",
            PartScheme::RecursiveBisection => "-ptype=rb",
        });

        cmd.arg(match self.ctype {
            Ctype::Rm => "-ctype=rm",
            Ctype::Shem => "-ctype=shem",
        });

        if self.ptype == PartScheme::Kway {
            cmd.arg(match self.objtype {
                Objtype::Cut => "-objtype=cut",
                Objtype::Volume => "-objtype=vol",
            });
            if self.contig {
                cmd.arg("-contig");
            }
            if self.minconn {
                cmd.arg("-minconn");
            }
        } else {
            cmd.arg(match self.iptype {
                Iptype::Grow => "-iptype=grow",
                Iptype::Random => "-iptype=random",
            });
        }

        cmd.arg(format!("-niter={}", self.niter));
        cmd.arg(format!("-ncuts={}", self.ncuts));
        cmd.arg(format!("-seed={}", self.seed));
        cmd.arg("-nooutput");
        cmd.stdin(Stdio::piped());
        if self.settings.nocapture {
            cmd.stdout(Stdio::inherit());
        } else {
            cmd.stdout(Stdio::null());
        }
        cmd.stderr(Stdio::piped());


        let mut child = cmd.spawn().expect("could not spawn process");

        let w = child.stdin.take().unwrap();
        let _ = self.graph.write(w);

        let output = child.wait_with_output().expect("command wasn't running");

        let status = output.status;
        let stderr = output.stderr;

        if self.settings.nocapture {
            eprintln!("{}", String::from_utf8_lossy(&stderr))
        }

        // stderr is generally super short, so this isn't a big deal
        if output.status.signal().is_some() || stderr.windows(const { b"Assertion".len() }).any(|s| s == b"Assertion") {
            ProgramResult::Crash
        } else if !status.success() {
            // happens if graph is malformed
            ProgramResult::GracefulFailure
        } else {
            ProgramResult::NormalExit
        }
    }

    /// runs the program `trials` times, incrementing the seed each time, returning the number of
    /// crashes
    #[allow(dead_code)]
    pub fn measure_robustness(&self, trials: usize) -> usize {
        let mut cnt = 0;
        for trial in 0..trials {
            let mut case = self.clone();
            case.seed += trial as i32;
            if case.run_program().is_crash() {
                cnt += 1;
            }
        }
        cnt
    }

    /// tries to approximate the cost of partitioning
    pub fn cost(&self) -> CaseCost {
        self.into()
    }

    /// factors that go into computing [`CaseCost`]
    pub fn cost_factors(&self) -> CostFactors {
        self.into()
    }

    pub fn stats(&self) -> String {
        use std::fmt::Write;

        let mut buf = String::with_capacity(200);
        let Case { graph, minconn, contig, niter, ncuts, nparts, ufactor, ptype, iptype, objtype, ctype, seed, settings: _ } = self;
        let nvtxs = graph.nvtxs();
        let nedges = graph.nedges();
        let Graph { xadj: _, adjncy: _, ncon, adjwgt, vwgt: _, vsize } = &graph;
        writeln!(buf, "Graph ({})", self.cost().human).unwrap();
        writeln!(buf, "\tnvtxs: {nvtxs}").unwrap();
        writeln!(buf, "\tnedges: {nedges}").unwrap();
        writeln!(buf, "\tncon: {ncon}").unwrap();
        writeln!(buf, "\tadjwgt: {}", adjwgt.is_some()).unwrap();
        writeln!(buf, "\tvsize: {}", vsize.is_some()).unwrap();
        writeln!(buf, "\tniter: {niter}").unwrap();
        writeln!(buf, "\tncuts: {ncuts}").unwrap();
        writeln!(buf, "\tnparts: {nparts}").unwrap();
        writeln!(buf, "\tufactor: {ufactor}").unwrap();
        writeln!(buf, "\tptype: {ptype:?}").unwrap();
        if self.ptype.is_rb() {
            writeln!(buf, "\tiptype: {iptype:?}").unwrap();
        }
        if self.ptype.is_kway() {
            writeln!(buf, "\tobjtype: {objtype:?}").unwrap();
            writeln!(buf, "\tcontig: {contig}").unwrap();
            writeln!(buf, "\tminconn: {minconn}").unwrap();
        }
        writeln!(buf, "\tctype: {ctype:?}").unwrap();
        writeln!(buf, "\tseed: {seed}").unwrap();
        // writeln!(buf, "\tcommand: {comm}", comm = comm.display()).unwrap();

        buf
    }
}
