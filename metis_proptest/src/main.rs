//! property testing the binaries of the original metis -- this is primarily a test case reducer

use std::{ffi::OsString, path::PathBuf, process::ExitCode, sync::Arc, time::Duration};
use clap::{ValueEnum, Parser};
use fastrand::Rng;
use minimize::MinimizationSet;
use strategy::{Case, Settings, StrategyKind};

mod graph;
mod utils;
mod strategy;
mod costs;
mod minimize;


#[derive(Parser)]
struct Cli {
    /// base graph file
    graph_file: PathBuf,

    /// number of partitions
    nparts: usize,

    /// The METIS standalone program to run
    #[arg(default_value = "gpmetis")]
    comm: OsString,

    /// always run graphchk before running comm
    #[arg(short = 'k', long)]
    always_graphchk: bool,

    #[arg(short, long, default_value = "kway")]
    ptype: PartScheme,

    /// only applies for rb
    #[arg(short, long, default_value = "grow")]
    iptype: Iptype,

    /// only applies for kway
    #[arg(short, long, default_value = "cut")]
    objtype: Objtype,

    /// only applies for kway
    #[arg(short, long)]
    minconn: bool,

    /// only applies for kway
    #[arg(short = 'g', long)]
    contig: bool,

    #[arg(short = 't', long, default_value = "1")]
    ncuts: usize,

    #[arg(short = 'r', long, default_value = "10")]
    niter: usize,

    #[arg(short, long, default_value = "1234")]
    seed: i32,

    #[arg(short, long, default_value = "30")]
    ufactor: u32,

    #[arg(short, long, default_value = "shem")]
    ctype: Ctype,

    #[arg(long)]
    nocapture: bool,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, ValueEnum)]
enum Iptype {
    #[default]
    Grow,
    Random
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, ValueEnum)]
enum PartScheme {
    #[default]
    Kway,
    #[value(alias = "rb")]
    RecursiveBisection
}

impl PartScheme {
    /// Returns `true` if the part scheme is [`Kway`].
    ///
    /// [`Kway`]: PartScheme::Kway
    #[must_use]
    fn is_kway(&self) -> bool {
        matches!(self, Self::Kway)
    }

    /// Returns `true` if the part scheme is [`RecursiveBisection`].
    ///
    /// [`RecursiveBisection`]: PartScheme::RecursiveBisection
    #[must_use]
    fn is_rb(&self) -> bool {
        matches!(self, Self::RecursiveBisection)
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, ValueEnum)]
enum Objtype {
    #[default]
    Cut,
    #[value(alias = "vol")]
    Volume
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, ValueEnum)]
enum Ctype {
    #[default]
    Rm,
    Shem
}

fn build_strats(m: &mut MinimizationSet) {
    use strategy::strategies as s;

    m.push(s::UnsetMinConn);
    m.push(s::UnsetContig);
    m.push(s::RemoveVsize);
    m.push(s::ShrinkNcon);

    m.push(s::DeleteVtxs {
        prob: 2,
        adjust_nparts: true,
        kind: StrategyKind::Checkpoint,
        rng: Rng::new(),
    });
    m.push(s::DeleteVtxs {
        prob: 2,
        adjust_nparts: true,
        kind: StrategyKind::Checkpoint,
        rng: Rng::new(),
    });
    m.push(s::DeleteVtxs {
        prob: 4,
        adjust_nparts: true,
        kind: StrategyKind::Repeat,
        rng: Rng::new(),
    });

    m.push(s::ReduceVwgt {
        amt_percent: 50,
    });
    m.push(s::ReduceAdjwgt {
        amt_percent: 75,
    });
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let Cli { graph_file, nparts, comm, always_graphchk: _, ptype, iptype, objtype, minconn, contig, ncuts, niter, seed, ctype, ufactor, nocapture } = cli;
    let graph = graph::Graph::read_from_path(graph_file).unwrap();
    let initial_case = Case {
        graph,
        minconn,
        contig,
        niter,
        ncuts,
        nparts,
        ufactor,
        ptype,
        iptype,
        objtype,
        ctype,
        seed,
        settings: Arc::new(Settings {
            comm: comm.into(),
            nocapture,
        }),
    };

    let initial_res = initial_case.run_program();
    if initial_res.is_normal_exit() {
        eprintln!("Graph was partioned successfully -- nothing to do");
        return ExitCode::FAILURE
    }
    if initial_res.is_graceful_failure() {
        eprintln!("Graph had graceful failure");
        return ExitCode::FAILURE
    }

    let mut strats = minimize::MinimizationSet::new();
    build_strats(&mut strats);

    let mut worker = strats.work(initial_case);
    let mut dur = Duration::from_secs_f64(0.125);
    let mut iterations_without_improvement = 0;
    for _ in 0..50 {
        eprintln!("working for {:.3} sec...", dur.as_secs_f64());
        let Some(improved) = worker.run_for(dur) else { break };
        if improved {
            iterations_without_improvement = 0;
            eprintln!("Found better graph:");
            eprintln!("{}", worker.best_rt().stats())
        } else {
            dur *= 2;
            dur = dur.min(Duration::from_secs(30));
            iterations_without_improvement += 1;
            if iterations_without_improvement >= 5 {
                break
            }
        };
    }

    eprintln!("best graph found (run time):");
    worker.best_rt().graph.write(std::io::stderr()).unwrap();
    eprintln!("{}", worker.best_rt().stats());

    eprintln!("best graph found (human):");
    worker.best_hu().graph.write(std::io::stdout()).unwrap();
    eprintln!("{}", worker.best_hu().stats());


    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    #[test]
    fn test_cli() {
        Cli::command().debug_assert();
    }
}
