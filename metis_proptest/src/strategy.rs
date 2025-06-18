use std::{ffi::OsStr, os::unix::process::ExitStatusExt, sync::Arc};

use crate::{costs::{CaseCost, CostFactors}, graph::Graph, Ctype, Iptype, Objtype, PartScheme};

pub mod simple;
pub mod combinator;
pub use combinator::*;
pub mod delete_vtx;

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

pub trait Strategy {
    /// if [`Strategy::is_valid`] is true for a case, don't go back before this one
    fn dont_backtrack(&self) -> bool { false }

    fn cost(&self, case: &Case) -> CaseCost;

    fn is_valid(&self, case: &Case) -> bool {
        let _= case;
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
        writeln!(buf, "Graph ({})", self.cost().human).unwrap();
        writeln!(buf, "\tnvtxs: {nvtxs}").unwrap();
        writeln!(buf, "\tnedges: {nedges}").unwrap();
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
