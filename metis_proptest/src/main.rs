//! property testing the binaries of the original metis -- this is primarily a test case reducer

use std::{ffi::OsString, path::PathBuf, process::ExitCode, sync::Arc, time::Duration};
use clap::{ValueEnum, Parser};
use fastrand::Rng;
use strategy::{Case, Settings};

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
    strats.push(strategy::simple::UnsetMinConn);
    strats.push(strategy::simple::UnsetContig);
    strats.push(strategy::simple::RemoveVsize);
    strats.push(strategy::simple::RemoveAllVWgt);
    for i in 0..10_usize {
        strats.push(strategy::delete_vtx::DeleteVtxs {
            prob: 2 + i,
            rng: Rng::new(),
            adjust_nparts: true,
            allow_backtrack: i > 5,
        });
    }

    let mut worker = strats.work(initial_case);
    let mut dur = Duration::from_secs(1);
    for _ in 0..10 {
        eprintln!("working for {:.3} sec...", dur.as_secs_f64());
        if worker.run_for(dur) {
            eprintln!("Found better graph:");
            eprintln!("{}", worker.best().stats())
        } else {
            dur *= 2;
            dur = dur.min(Duration::from_secs(30));
        };
    }

    eprintln!("best graph found:");
    worker.best().graph.write(std::io::stdout()).unwrap();


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
