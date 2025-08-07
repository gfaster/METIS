#![cfg(test)]
#![allow(unused_mut)]

use std::collections::HashMap;
use std::mem::transmute;
use std::path::{Path, PathBuf};
use std::ptr;
use std::sync::{Arc, LazyLock, OnceLock, RwLock};

use crate::bindings::{METIS_NOPTIONS, METIS_OK, idx_t, real_t};
use crate::dyncall::{ab_test, ab_test_eq, ab_test_multi};
use crate::graph_gen::{Csr, GraphBuilder};
use crate::kmetis::METIS_PartGraphKway;
use crate::pmetis::METIS_PartGraphRecursive;
use crate::{Objtype, Optype};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum TestGraph {
    Elt4,
    Youtube,
    Webbase2004,
    WebSpam,
    Luxembourg,
    Orani,
    CitDblp,
}

impl TestGraph {
    #![deny(clippy::wildcard_enum_match_arm)]

    pub const fn is_smallish(self) -> bool {
        match self {
            TestGraph::Elt4 |
            TestGraph::Webbase2004 |
            TestGraph::WebSpam |
            TestGraph::CitDblp => true,
            TestGraph::Orani |
            TestGraph::Youtube |
            TestGraph::Luxembourg => false,
        }
    }

    const fn rev_idx(idx: usize) -> Self {
        Self::try_rev_idx(idx).expect("invalid index")
    }

    const fn try_rev_idx(idx: usize) -> Option<Self> {
        let ret = match idx {
            0 => TestGraph::Elt4,
            1 => TestGraph::Youtube,
            2 => TestGraph::Webbase2004,
            3 => TestGraph::WebSpam,
            4 => TestGraph::Luxembourg,
            5 => TestGraph::Orani,
            6 => TestGraph::CitDblp,
            _ => return None,
        };
        debug_assert!(ret.idx() == idx);
        Some(ret)
    }

    const fn idx(self) -> usize {
        self as usize
    }

    const fn file(self) -> &'static str {
        match self {
            TestGraph::Elt4 => "4elt.graph",
            TestGraph::Youtube => "soc-youtube.graph",
            TestGraph::Webbase2004 => "web-webbase-2001.graph",
            TestGraph::WebSpam => "web-spam.graph",
            TestGraph::Luxembourg => "road-luxembourg-osm.graph",
            TestGraph::Orani => "econ-orani678.graph",
            TestGraph::CitDblp => "cit-DBLP.graph",
        }
    }

    pub fn test_suite() -> impl Iterator<Item = Self> {
        static DO_BIG: LazyLock<bool> =
            LazyLock::new(|| std::env::var_os("DO_BIG") == Some("1".into()));
        let do_big = *DO_BIG;
        Self::ALL
            .into_iter()
            .filter(move |g| g.is_smallish() || do_big)
    }

    pub fn is_contiguous(self) -> bool {
        match self {
            TestGraph::Elt4
            | TestGraph::Youtube
            | TestGraph::Webbase2004
            | TestGraph::WebSpam
            | TestGraph::Luxembourg
            | TestGraph::Orani => true,
            TestGraph::CitDblp => false,
        }
    }

    const COUNT: usize = 7;

    const ALL: [Self; Self::COUNT] = [
        Self::Elt4,
        Self::Youtube,
        Self::Webbase2004,
        Self::WebSpam,
        Self::Luxembourg,
        Self::Orani,
        Self::CitDblp,
    ];
}

/// uses caches to reduce re-parsing
pub(crate) fn read_graph(graph: TestGraph, op: Optype, nparts: usize, ncon: usize) -> GraphBuilder {
    // This is a bit of an awkward synchronization solution, but it's the best I can think of

    fn load_graph_idx(idx: usize) -> Csr {
        let filename = TestGraph::rev_idx(idx).file();
        let path = Path::new("graphs").join(filename);
        let f = std::fs::read(&path)
            .map_err(|e| format!("Could not read graph {path}: {e}", path = path.display()))
            .unwrap();
        GraphBuilder::read_graph_cfg_index(std::io::Cursor::new(&f), Optype::Kmetis, 1, 1, true)
            .unwrap()
            .into_csr()
    }

    const fn graph_of<const N: usize>() -> fn() -> Csr {
        || load_graph_idx(N)
    }

    // probably false sharing, but I don't care enough to fix it
    static GRAPHS: [LazyLock<Csr, fn() -> Csr>; TestGraph::COUNT] = [
        LazyLock::new(graph_of::<0>()),
        LazyLock::new(graph_of::<1>()),
        LazyLock::new(graph_of::<2>()),
        LazyLock::new(graph_of::<3>()),
        LazyLock::new(graph_of::<4>()),
        LazyLock::new(graph_of::<5>()),
        LazyLock::new(graph_of::<6>()),
    ];

    GraphBuilder::from_csr(GRAPHS[graph.idx()].clone(), op, nparts, ncon)
}

/// Runs `f` on each of the small-ish graphs
#[cfg(test)]
pub(crate) fn for_test_suite<F>(op: Optype, nparts: usize, ncon: usize, mut f: F)
where
    F: FnMut(TestGraph, GraphBuilder),
{
    fn inner(op: Optype, nparts: usize, ncon: usize, f: &mut dyn FnMut(TestGraph, GraphBuilder)) {
        for graph in TestGraph::test_suite() {
            fastrand::seed(12513471239123 + graph as u64);
            eprintln!("Testing with {graph:?}");
            let builder = GraphBuilder::test_graph(graph, op, nparts, ncon);
            f(graph, builder);
        }
    }

    inner(op, nparts, ncon, &mut f);
}

/// Runs `f` on each of the small-ish graphs if it passes the filter, and performs `ab_test_single_eq` on partitioning each
#[cfg(test)]
pub(crate) fn ab_test_partition_test_graphs_filter<F>(
    overrides: &str,
    op: Optype,
    nparts: usize,
    ncon: usize,
    mut f: F,
) where
    F: FnMut(TestGraph, GraphBuilder) -> Option<GraphBuilder>,
{
    use crate::dyncall::ab_test_single_eq;

    fn inner(
        overrides: &str,
        op: Optype,
        nparts: usize,
        ncon: usize,
        f: &mut dyn FnMut(TestGraph, GraphBuilder) -> Option<GraphBuilder>,
    ) {
        let mut tested = false;
        for_test_suite(op, nparts, ncon, |tg, graph| {
            if let Some(graph) = f(tg, graph) {
                tested = true;
                ab_test_single_eq(overrides, || graph.clone().call());
            }
        });
        if !tested {
            panic!("no graphs passed the filter")
        }
    }

    inner(overrides, op, nparts, ncon, &mut f);
}

/// Runs `f` on the graph provided, and performs `ab_test_single_eq` on partitioning each
#[cfg(test)]
pub(crate) fn ab_test_partition_test_graph<F>(
    overrides: &str,
    op: Optype,
    nparts: usize,
    ncon: usize,
    tg: TestGraph,
    mut f: F,
) where
    F: FnMut(GraphBuilder) -> GraphBuilder,
{
    use crate::dyncall::ab_test_single_eq;

    fn inner(
        overrides: &str,
        op: Optype,
        nparts: usize,
        ncon: usize,
        tg: TestGraph,
        f: &mut dyn FnMut(GraphBuilder) -> GraphBuilder,
    ) {
        fastrand::seed(85103750);
        eprintln!("Testing with {tg:?}");
        let builder = GraphBuilder::test_graph(tg, op, nparts, ncon);
        let graph = f(builder);
        ab_test_single_eq(overrides, || graph.clone().call());
    }

    inner(overrides, op, nparts, ncon, tg, &mut f);
}

/// Runs `f` on each of the small-ish graphs and performs `ab_test_single_eq` on partitioning each
#[cfg(test)]
pub(crate) fn ab_test_partition_test_graphs<F>(
    overrides: &str,
    op: Optype,
    nparts: usize,
    ncon: usize,
    mut f: F,
) where
    F: FnMut(GraphBuilder) -> GraphBuilder,
{
    ab_test_partition_test_graphs_filter(overrides, op, nparts, ncon, |_, g| Some(f(g)));
}

/// function signature of METIS_PartGraphKway, METIS_PartGraphRecursive
// type PartSig = unsafe extern "C" fn(
//     *mut idx_t,
//     *mut idx_t,
//     *mut idx_t,
//     *mut idx_t,
//     *mut idx_t,
//     *mut idx_t,
//     *mut idx_t,
//     *mut idx_t,
//     *mut real_t,
//     *const real_t,
//     *mut idx_t,
//     *mut idx_t,
//     *mut idx_t,
// ) -> idx_t;

#[test]
fn basic_part_graph_recursive() {
    let mut xadj = &mut [0, 1, 2];
    let mut adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let mut vwgt = ptr::null_mut();
    let mut vsize = ptr::null_mut();
    let mut adjwgt = ptr::null_mut();
    let mut nparts: idx_t = 2;
    let mut tpwgts = ptr::null_mut();
    let mut ubvec = ptr::null_mut();
    let mut objval: idx_t = 0;
    let mut part = [0; 2];
    let mut options = [-1; METIS_NOPTIONS as usize];

    let res = unsafe {
        METIS_PartGraphRecursive(
            &mut nvtxs as *mut _,
            &mut ncon as *mut _,
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            vwgt,
            vsize,
            adjwgt,
            &mut nparts as *mut _,
            tpwgts,
            ubvec,
            options.as_mut_ptr(),
            &mut objval as *mut _,
            part.as_mut_ptr(),
        )
    };

    assert_eq!(res, METIS_OK);
    assert_ne!(part[0], part[1]);
}

#[test]
fn basic_part_graph_kway() {
    // let mut graph = GraphBuilder::new_basic(Optype::Kmetis, 2);
    // graph.with_vtx_degrees(vec![3, 4, 2, 3, 5, 8, 1, 3, 2, 5]);
    // graph.random_edges();
    // assert!(graph.call().is_ok());
}

fn part_graph_and_verify(
    ncon: idx_t,
    options: &mut [i32; METIS_NOPTIONS as usize],
    nparts: idx_t,
    use_vwgt: bool,
    use_adjwgt: bool,
    partfn: Optype,
) {
    let exec = || {
        let mut graph = read_graph(TestGraph::Elt4, partfn, nparts as usize, ncon as usize);

        graph.set_from_options_arr(options);

        if use_vwgt {
            graph.random_vwgt();
        }
        if use_adjwgt {
            graph.random_vwgt();
        }

        let res = graph.call();

        assert!(res.is_ok());

        let (objval, part) = res.unwrap();

        graph.verify_part(objval, &part);
    };
    exec();
    // ab_test_eq("*", exec);
}

macro_rules! partfn_type {
    ($fn:ident) => {
        {
            #[allow(unused_must_use)]
            const _: () = {
                || {
                    #[allow(unused)]
                    use $fn;
                };
            };
            partfn_type!(@inner $fn)
        }
    };
    (@inner METIS_PartGraphKway) => {
        crate::Optype::Kmetis
    };
    (@inner METIS_PartGraphRecursive) => {
        crate::Optype::Pmetis
    };
}

macro_rules! part_test {
    (
    name: $name:ident,
    options: $options:expr,
    nparts: $nparts:literal,
    ncon: $ncon:literal,
    vwgt: $use_vwgt:literal,
    adjwgt: $use_adjwgt:literal,
    partfn: $partfn:ident,
    extra: $extra:ident,
    ) => {
        #[cfg(feature = "extra_tests")]
        #[test]
        fn $name() {
            let mut options = $options;
            part_graph_and_verify(
                $ncon,
                &mut options,
                $nparts,
                $use_vwgt,
                $use_adjwgt,
                partfn_type!($partfn),
            );
        }
    };
    (
    name: $name:ident,
    options: $options:expr,
    nparts: $nparts:literal,
    ncon: $ncon:literal,
    vwgt: $use_vwgt:literal,
    adjwgt: $use_adjwgt:literal,
    partfn: $partfn:ident,
    ) => {
        #[test]
        fn $name() {
            let mut options = $options;
            part_graph_and_verify(
                $ncon,
                &mut options,
                $nparts,
                $use_vwgt,
                $use_adjwgt,
                partfn_type!($partfn),
            );
        }
    };
}

macro_rules! part_test_set {
    (
    set_name: $name:ident,
    options: $options:expr,
    partfn: $partfn:ident,
    ) => {
        mod $name {
            use super::*;

            part_test! {
                name: large_basic,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: false,
                adjwgt: false,
                partfn: $partfn,
            }

            part_test! {
                name: large_1con_vwgt,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: true,
                adjwgt: false,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_2con_vwgt,
                options: $options,
                nparts: 20,
                ncon: 2,
                vwgt: true,
                adjwgt: false,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_1con_vwgt_adjwgt,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: true,
                adjwgt: true,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_2con_vwgt_adjwgt,
                options: $options,
                nparts: 20,
                ncon: 2,
                vwgt: true,
                adjwgt: true,
                partfn: $partfn,
            }

            part_test! {
                name: large_halve,
                options: $options,
                nparts: 2,
                ncon: 1,
                vwgt: false,
                adjwgt: false,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_1con_adjwgt,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: false,
                adjwgt: true,
                partfn: $partfn,
                extra: true,
            }
        }
    };
}

// this could be shrunk but it works
macro_rules! part_test_hyper_set {
    ($partfn:ident as $name:ident => [$({$($first_opt:ident),*}),*]) => {
        mod $name {
            #![allow(non_snake_case, unused_imports)]
            use super::*;
            part_test_hyper_set!(@call_for $partfn : ($(($($first_opt),*)),*) @ ());
        }
    };
    (@call $partfn:ident, $name:ident: ($($option:ident),*)) => {
        // $(part_test_hyper_set!(@dbg $option));*
        part_test_set!(
            set_name: $name,
            options: make_options!($($option)*),
            partfn: $partfn,
        );
    };
    (@call_for $partfn:ident: ($first:tt $(, $rest:tt)*) @ $keep:tt) => {
                // $(part_test_hyper_set!(@dbg $rest))*;
                // part_test_hyper_set!(@dbg $keep);
        part_test_hyper_set!(@call_for_inner $partfn: $first @ ($($rest),*) @ $keep);
    };
    (@call_for_inner $partfn:ident: ($($iter:ident),*) @ () @ $keep:tt) => {
        $(
        // mod $iter {
        //     use super::*;
        // part_test_hyper_set!(@dbg $iter);
            part_test_hyper_set!(@rejoin_call $partfn, $iter: $keep);
        // }
        )*
    };
    (@call_for_inner $partfn:ident: ($($iter:ident),*) @ $rest:tt @ $keep:tt) => {
        $(
            mod $iter {
                use super::*;
        // part_test_hyper_set!(@dbg $keep);
                part_test_hyper_set!(@rejoin $partfn: $rest @ ($iter, $keep));
            }
        )*
    };
    (@rejoin $partfn:ident: $rest:tt @ ($parent:ident, ($($keep:ident),*))) => {
        // part_test_hyper_set!(@dbg $parent);
        // part_test_hyper_set!(@dbg $rest);
        part_test_hyper_set!(@call_for $partfn: $rest @ ($parent $(, $keep)*));
    };
    (@rejoin_call $partfn:ident, $name:ident: ($($keep:ident),*)) => {
        part_test_hyper_set!(@call $partfn, $name: ($name $(, $keep)*));
    };
    (@dbg $var:expr) => {
        part_test_hyper_set!(@dbg_inner $var);
    };
    (@dbg_inner a) => {
    };
}

// TODO: test graphs with odd degree (very high, power law)
// I suspect problems if there is vertex degree > nparts

// Volume + Contig fails even when my code isn't used at all. I should verify this later, but it
// appears to be an issue in the original as well.

// Kmetis only allows Grow and Rb initial partitioning
part_test_hyper_set!(METIS_PartGraphKway as kmetisvol => [{Vol}, {Grow, Rb}, {Rm, Shem}]);
part_test_hyper_set!(METIS_PartGraphKway as kmetiscut => [{Cut}, {Grow, Rb}, {Rm, Shem}, {Contig, None}]);

// Communication volume is illegal on PMETIS routines
// Edge and Node initial partitioning also illegal in PMETIS
part_test_hyper_set!(METIS_PartGraphRecursive as pmetis => [{Cut}, {Grow, Random}, {Rm, Shem}]);

#[test]
fn identical_to_c_kmetis_cut() {
    ab_test_partition_test_graphs("*", Optype::Kmetis, 20, 1, |mut g| {
        g.set_objective(Objtype::Cut);
        g.random_vwgt();
        g.random_adjwgt();
        g
    });
}

#[test]
fn identical_to_c_kmetis_vol() {
    ab_test_partition_test_graphs("*", Optype::Kmetis, 20, 1, |mut g| {
        g.set_objective(Objtype::Vol);
        g.random_vwgt();
        g.random_adjwgt();
        g
    });
}

#[test]
fn identical_to_p_kmetis() {
    ab_test_partition_test_graphs("*", Optype::Pmetis, 20, 1, |mut g| {
        g.random_vwgt();
        g.random_adjwgt();
        g
    });
}

#[test]
#[ignore = "slow"]
fn identical_to_c_kmetis_large() {
    let mut rng = fastrand::Rng::with_seed(3971056);
    for trial in 0..1 {
        let prev_seed = fastrand::get_seed();
        eprintln!("trial: {trial:>3}");

        let thread_seed = rng.fork().get_seed();
        let graph_seed = rng.i32(0..=i32::MAX);
        eprintln!("thread seed: {thread_seed:16x}");
        eprintln!(
            "graph seed:  {graph_seed:16x}",
            graph_seed = graph_seed as usize
        );
        fastrand::seed(thread_seed);

        let mut graph = read_graph(TestGraph::Youtube, Optype::Kmetis, 20, 1);
        graph.random_vwgt(); // takes slightly longer, but exposes some bugs
        graph.set_seed(graph_seed);
        graph.set_objective(Objtype::Cut);

        // graph.enable_dbg(crate::DbgLvl::Info);
        // graph.enable_dbg(crate::DbgLvl::Refine);
        // graph.enable_dbg(crate::DbgLvl::Coarsen);
        // graph.enable_dbg(crate::DbgLvl::Ipart);

        let exec = || {
            let mut graph = graph.clone();

            let res = graph.call();

            let (objval, part) = res.unwrap();

            graph.verify_part(objval, &part);

            (objval, part)
        };
        ab_test_eq("*", exec);
        fastrand::seed(prev_seed);
    }
}

#[test]
fn identical_to_c_kmetis_multiconstraint() {
    ab_test_partition_test_graphs("*", Optype::Kmetis, 20, 2, |mut g| {
        g.random_vwgt();
        g.call().unwrap();
        g
    });
}

#[test]
fn debug_assertions_enabled() {
    if !cfg!(debug_assertions) {
        panic!("debug assertions are disabled")
    }
}

#[test]
fn test_graph_is_contiguous_accurate() {
    for_test_suite(Optype::Kmetis, 1, 1, |tg, g| {
        assert_eq!(tg.is_contiguous(), g.is_contiguous());
    });
}

#[test]
fn test_graph_all_populated() {
    let mut i = 0..;
    let actual: Vec<_> = std::iter::from_fn(|| TestGraph::try_rev_idx(i.next().unwrap()))
        .fuse()
        .collect();
    assert_eq!(actual, TestGraph::ALL);
}
