/*
\file
\brief Functions for computing matchings during graph coarsening

\date Started 7/23/97
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota
\version\verbatim $Id: coarsen.c 20398 2016-11-22 17:17:12Z karypis $ \endverbatim
*/

/*
 * Gavin: This helped me understand what's going on here:
 * https://andreasloukas.blog/2018/05/26/demystifying-graph-coarsening/
 * https://ieeexplore.ieee.org/document/6816682 <= has quite a bit of detail on METIS
 */

// CHECK THE TODO COMMENTS DANGIT
// I know there will be issues, and it's probably going to be there!

use crate::*;

const UNMATCHEDFOR2HOP: f32 = 0.10; /* The fraction of unmatched vertices that triggers 2-hop */
const UNMATCHED: i32 = -1;

/*************************************************************************/
/* This function takes a graph and creates a sequence of coarser graphs.
   It implements the coarsening phase of the multilevel paradigm.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn CoarsenGraph(ctrl: *mut ctrl_t, graph: *mut graph_t) -> *mut graph_t {
    let ctrl = ctrl.as_mut().unwrap();
    // ctrl.dbglvl |= METIS_DBG_COARSEN;
    let mut graph = graph;
    // idx_t i, eqewgts, level=0;

    // ifset!(
    //     ctrl.dbglvl,
    //     METIS_DBG_TIME,
    //     gk_startcputimer(ctrl.CoarsenTmr),
    // );

    /* determine if the weights on the edges are all the same */
    let mut eqewgts = true;
    get_graph_slices!(*graph => adjwgt tvwgt);
    for i in (1)..((*graph).nedges) {
        if adjwgt[0 as usize] != adjwgt[i as usize] {
            eqewgts = false;
            break;
        }
    }

    /* set the maximum allowed coarsest vertex weight */
    for i in (0)..((*graph).ncon) {
        // ctrl.maxvwgt[i as usize] = 1.5 * tvwgt[i as usize] / ctrl.CoarsenTo;
        *ctrl.maxvwgt.add(i as usize) = 3 * tvwgt[i as usize] / (2 * ctrl.CoarsenTo);
    }
    // let mut level = 0;
    loop {
        ifset!(
            ctrl.dbglvl,
            METIS_DBG_COARSEN,
            PrintCGraphStats(ctrl, graph),
        );

        /* allocate memory for cmap, if it has not already been done due to
        multiple cuts */
        if (*graph).cmap.is_null() {
            (*graph).cmap = imalloc(
                (*graph).nvtxs as usize,
                c"CoarsenGraph: graph.cmap".as_ptr(),
            ) as _;
        }

        /* determine which matching scheme you will use */
        match ctrl.ctype {
            METIS_CTYPE_RM => Match_RM(ctrl, graph),
            METIS_CTYPE_SHEM => {
                if eqewgts || (*graph).nedges == 0 {
                    Match_RM(ctrl, graph)
                } else {
                    Match_SHEM(ctrl, graph)
                }
            }
            _ => panic!("Unknown ctype: {}", ctrl.ctype),
        };

        graph::graph_WriteToDisk(ctrl, graph);

        graph = (*graph).coarser;
        eqewgts = false;
        // level += 1;

        debug_assert!(CheckGraph(graph, 0, 1) != 0);

        if !((*graph).nvtxs > ctrl.CoarsenTo
            && ((*graph).nvtxs as real_t) < COARSEN_FRACTION * (*(*graph).finer).nvtxs as real_t
            && (*graph).nedges > (*graph).nvtxs / 2)
        {
            break;
        }
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_COARSEN,
        PrintCGraphStats(ctrl, graph),
    );
    // ifset!(
    //     ctrl.dbglvl,
    //     METIS_DBG_TIME,
    //     gk_stopcputimer(ctrl.CoarsenTmr),
    // );

    return graph;
}

/*************************************************************************/
/* This function takes a graph and creates a sequence of nlevels coarser
   graphs, where_ nlevels is an input parameter.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn CoarsenGraphNlevels(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    nlevels: idx_t,
) -> *mut graph_t {
    let ctrl = ctrl.as_mut().unwrap();
    let mut graph = graph;
    // idx_t i, eqewgts, level;

    // ifset!(
    //     ctrl.dbglvl,
    //     METIS_DBG_TIME,
    //     gk_startcputimer(ctrl.CoarsenTmr),
    // );

    /* determine if the weights on the edges are all the same */
    let mut eqewgts = 1;
    {
        let graph = graph.as_mut().unwrap();
        get_graph_slices!(graph => adjwgt tvwgt);
        for i in (1)..(graph.nedges) {
            if adjwgt[0 as usize] != adjwgt[i as usize] {
                eqewgts = 0;
                break;
            }
        }

        /* set the maximum allowed coarsest vertex weight */
        for i in (0)..(graph.ncon) {
            *ctrl.maxvwgt.add(i as usize) = 3 * tvwgt[i as usize] / (2 * ctrl.CoarsenTo);
        }
    }

    for _level in (0)..(nlevels) {
        ifset!(
            ctrl.dbglvl,
            METIS_DBG_COARSEN,
            PrintCGraphStats(ctrl, graph),
        );

        /* allocate memory for cmap, if it has not already been done due to
        multiple cuts */
        if (*graph).cmap.is_null() {
            (*graph).cmap = imalloc(
                (*graph).nvtxs as usize,
                c"CoarsenGraph: graph.cmap".as_ptr(),
            ) as _;
        }

        /* determine which matching scheme you will use */
        match ctrl.ctype {
            METIS_CTYPE_RM => Match_RM(ctrl, graph),
            METIS_CTYPE_SHEM => {
                if eqewgts != 0 || (*graph).nedges == 0 {
                    Match_RM(ctrl, graph)
                } else {
                    Match_SHEM(ctrl, graph)
                }
            }
            _ => panic!("Unknown ctype: {}", ctrl.ctype),
        };

        graph::graph_WriteToDisk(ctrl, graph);

        graph = (*graph).coarser;
        eqewgts = 0;

        assert!(CheckGraph(graph, 0, 1) != 0);

        if (*graph).nvtxs < ctrl.CoarsenTo
            || ((*graph).nvtxs as real_t) > COARSEN_FRACTION * (*(*graph).finer).nvtxs as real_t
            || (*graph).nedges < (*graph).nvtxs / 2
        {
            break;
        }
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_COARSEN,
        PrintCGraphStats(ctrl, graph),
    );
    // ifset!(
    // ctrl.dbglvl,
    // METIS_DBG_TIME,
    // gk_stopcputimer(ctrl.CoarsenTmr),
    // );

    return graph;
}

/*************************************************************************/
/* This function finds a matching by randomly selecting one of the
   unmatched adjacent vertices.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_RM(ctrl: *mut ctrl_t, graph: *mut graph_t) -> idx_t {
    // println!("Calling Match_RM");
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, pi, ii, j, jj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx,
    // last_unmatched, avgdegree, bnum;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
    // idx_t *match_, *cmap, *degrees, *perm, *tperm;
    // size_t nunmatched=0;

    // WCOREPUSH;

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.MatchTmr));

    // adjwgt declared but not used in original
    get_graph_slices!(graph => adjncy vwgt xadj);
    get_graph_slices_mut!(graph => cmap);
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    // xadj = graph.xadj;
    // vwgt = graph.vwgt;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // cmap = graph.cmap;

    mkslice!(ctrl->maxvwgt, ncon * nvtxs);

    let mut match_: Vec<idx_t> = vec![UNMATCHED as i32; nvtxs];
    let mut perm: Vec<idx_t> = vec![0; nvtxs];
    let mut tperm: Vec<idx_t> = vec![0; nvtxs];
    let mut degrees: Vec<idx_t> = vec![0; nvtxs];

    /* Determine a "random" traversal order that is biased towards
    low-degree vertices */
    irandArrayPermute(nvtxs as idx_t, tperm.as_mut_ptr(), nvtxs as idx_t / 8, 1);

    let avgdegree = (4.0 * (xadj[nvtxs] / nvtxs as idx_t) as real_t) as idx_t;
    for i in (0)..(nvtxs) {
        let bnum = ((1 + xadj[(i + 1) as usize] - xadj[i as usize]) as f64).sqrt() as idx_t;
        degrees[i as usize] = avgdegree.min(bnum);
    }
    bucketsort::BucketSortKeysInc(
        ctrl,
        nvtxs as idx_t,
        avgdegree,
        degrees.as_ptr(),
        tperm.as_ptr(),
        perm.as_mut_ptr(),
    );

    /* Traverse the vertices and compute the matching */
    let mut cnvtxs = 0;
    let mut last_unmatched: usize = 0;
    let mut nunmatched: usize = 0;
    for pi in (0)..(nvtxs) {
        let i = perm[pi as usize] as usize;

        if match_[i as usize] == UNMATCHED {
            /* Unmatched */
            let mut maxidx = i;

            if if ncon == 1 {
                vwgt[i as usize] < maxvwgt[0 as usize]
            } else {
                // ivecle(ncon, vwgt[(i as usize * ncon)..], maxvwgt)
                util::ivecle(&vwgt[cntrng!(i * ncon, ncon)], &maxvwgt[..ncon])
            } {
                /* Deal with island vertices. Find a non-island and match it with.
                The matching ignores ctrl.maxvwgt requirements */
                if xadj[i as usize] == xadj[(i + 1) as usize] {
                    // last_unmatched=gk_max(pi, last_unmatched)+1;
                    last_unmatched = pi.max(last_unmatched) + 1;
                    while last_unmatched < nvtxs {
                        let j = perm[last_unmatched as usize] as usize;
                        if match_[j] == UNMATCHED {
                            maxidx = j;
                            break;
                        }
                        last_unmatched += 1;
                    }
                } else {
                    /* Find a random matching, subject to maxvwgt constraints */
                    if ncon == 1 {
                        /* single constraint version */
                        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                            let k = adjncy[j as usize] as usize;
                            if match_[k] == UNMATCHED
                                && vwgt[i as usize] + vwgt[k] <= maxvwgt[0 as usize]
                            {
                                maxidx = k;
                                break;
                            }
                        }

                        /* If it did not match_, record for a 2-hop matching. */
                        if maxidx == i && 2 * vwgt[i as usize] < maxvwgt[0 as usize] {
                            nunmatched += 1;
                            maxidx = (UNMATCHED as isize) as usize;
                        }
                    } else {
                        /* multi-constraint version */
                        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                            let k = adjncy[j as usize] as usize;
                            if match_[k as usize] == UNMATCHED
                                // && ivecaxpylez(ncon, 1, vwgt[(i as usize * ncon)..], vwgt[(k * ncon)..], maxvwgt))
                                && util::ivecaxpylez(1, &vwgt[cntrng!(i * ncon, ncon)], &vwgt[cntrng!(k * ncon, ncon)], &maxvwgt[..ncon])
                            {
                                maxidx = k;
                                break;
                            }
                        }

                        /* If it did not match_, record for a 2-hop matching. */
                        if maxidx == i
                            // && ivecaxpylez(ncon, 2, vwgt[(i * ncon)..], vwgt[(i * ncon)..], maxvwgt))
                            && util::ivecaxpylez( 2, &vwgt[cntrng!(i * ncon, ncon)], &vwgt[cntrng!(i * ncon, ncon)], &maxvwgt[..ncon])
                        {
                            nunmatched += 1;
                            maxidx = (UNMATCHED as isize) as usize;
                        }
                    }
                }
            }

            if maxidx as idx_t != UNMATCHED {
                cmap[maxidx] = cnvtxs;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
                match_[i as usize] = maxidx as idx_t;
                match_[maxidx as usize] = i as idx_t;
            }
        }
    }

    // println!("nunmatched: {nunmatched}");

    /* see if a 2-hop matching is required/allowed */
    if ctrl.no2hop == 0 && nunmatched as real_t > UNMATCHEDFOR2HOP * nvtxs as real_t {
        // original bound cnvtxs, but that has no effect
        Match_2Hop(
            ctrl,
            graph,
            perm.as_mut_ptr(),
            match_.as_mut_ptr(),
            cnvtxs,
            nunmatched,
        );
    }

    /* match_ the final unmatched vertices with themselves and reorder the vertices
    of the coarse graph for memory-friendly contraction */
    cnvtxs = 0;
    for i in (0)..(nvtxs) {
        if match_[i as usize] == UNMATCHED {
            match_[i as usize] = i as idx_t;
            cmap[i as usize] = cnvtxs;
            cnvtxs += 1;
        } else {
            if i as idx_t <= match_[i as usize] {
                cmap[match_[i as usize] as usize] = cnvtxs;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
            }
        }
    }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.MatchTmr));

    CreateCoarseGraph(ctrl, graph, cnvtxs, match_.as_mut_ptr());

    // WCOREPOP;
    return cnvtxs;
}

/**************************************************************************/
/* This function finds a matching using the HEM heuristic. The vertices
   are visited based on increasing degree to ensure that all vertices are
   given a chance to match_ with something.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_SHEM(ctrl: *mut ctrl_t, graph: *mut graph_t) -> idx_t {
    // println!("Calling Match_SHEM");
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, pi, ii, j, jj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx, maxwgt,
    // last_unmatched, avgdegree, bnum;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
    // idx_t *match_, *cmap, *degrees, *perm, *tperm;
    // size_t nunmatched=0;

    // WCOREPUSH;

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.MatchTmr));

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy adjwgt);
    get_graph_slices_mut!(ctrl, graph => vwgt cmap);
    // xadj = graph.xadj;
    // vwgt = graph.vwgt;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // cmap = graph.cmap;

    mkslice!(ctrl->maxvwgt, ncon * nvtxs);

    let mut match_ = vec![UNMATCHED; nvtxs];
    let mut perm = vec![0; nvtxs];
    let mut tperm = vec![0; nvtxs];
    let mut degrees = vec![0; nvtxs];

    /* Determine a "random" traversal order that is biased towards low-degree vertices */
    irandArrayPermute(nvtxs as idx_t, tperm.as_mut_ptr(), nvtxs as idx_t / 8, 1);

    let avgdegree = (4.0 * (xadj[nvtxs] as real_t / nvtxs as real_t)) as idx_t;
    for i in (0)..(nvtxs) {
        let bnum = ((1 + xadj[(i + 1) as usize] - xadj[i as usize]) as real_t).sqrt() as idx_t;
        degrees[i as usize] = (if bnum > avgdegree { avgdegree } else { bnum }) as idx_t;
    }
    bucketsort::BucketSortKeysInc(
        ctrl,
        nvtxs as idx_t,
        avgdegree,
        degrees.as_ptr(),
        tperm.as_ptr(),
        perm.as_mut_ptr(),
    );

    /* Traverse the vertices and compute the matching */
    let mut cnvtxs = 0;
    let mut last_unmatched = 0;
    let mut nunmatched = 0;
    for pi in (0)..(nvtxs) {
        let i = perm[pi as usize] as usize;

        if match_[i as usize] == UNMATCHED {
            /* Unmatched */
            let mut maxidx = i;
            let mut maxwgt = -1;

            if if ncon == 1 {
                vwgt[i as usize] < maxvwgt[0 as usize]
            } else {
                // ivecle(ncon, &vwgt[i * ncon], maxvwgt)
                util::ivecle(&vwgt[cntrng!(i * ncon, ncon)], &maxvwgt[..ncon])
            } {
                /* Deal with island vertices. Find a non-island and match_ it with.
                The matching ignores ctrl.maxvwgt requirements */
                if xadj[i as usize] == xadj[(i + 1) as usize] {
                    // last_unmatched = gk_max(pi, last_unmatched)+1;
                    last_unmatched = pi.max(last_unmatched) + 1;
                    while last_unmatched < nvtxs {
                        let j = perm[last_unmatched as usize];
                        if match_[j as usize] == UNMATCHED {
                            maxidx = j as usize;
                            break;
                        }
                        last_unmatched += 1;
                    }
                } else {
                    /* Find a heavy-edge matching, subject to maxvwgt constraints */
                    if ncon == 1 {
                        /* single constraint version */
                        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                            let k = adjncy[j as usize] as usize;
                            if match_[k as usize] == UNMATCHED
                                && maxwgt < adjwgt[j as usize]
                                && vwgt[i as usize] + vwgt[k] <= maxvwgt[0 as usize]
                            {
                                maxidx = k;
                                maxwgt = adjwgt[j as usize];
                            }
                        }

                        /* If it did not match_, record for a 2-hop matching. */
                        if maxidx == i && 2 * vwgt[i as usize] < maxvwgt[0 as usize] {
                            nunmatched += 1;
                            maxidx = (UNMATCHED as isize) as usize;
                        }
                    } else {
                        /* multi-constraint version */
                        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                            let k = adjncy[j as usize] as usize;
                            if match_[k] == UNMATCHED
                                // && ivecaxpylez(ncon, 1, vwgt[(i * ncon)..], vwgt[(k * ncon)..], maxvwgt)
                                && util::ivecaxpylez(1, &vwgt[cntrng!(i * ncon, ncon)], &vwgt[cntrng!(k * ncon, ncon)], &maxvwgt[..ncon])
                                && (maxwgt < adjwgt[j as usize]
                                    || (maxwgt == adjwgt[j as usize]
                                        && mcutil::BetterVBalance(
                                            ncon as idx_t,
                                            graph.invtvwgt,
                                            vwgt[(i * ncon)..].as_mut_ptr(),
                                            vwgt[(maxidx * ncon)..].as_mut_ptr(),
                                            vwgt[(k * ncon)..].as_mut_ptr(),
                                        ) != 0))
                            {
                                maxidx = k;
                                maxwgt = adjwgt[j as usize];
                            }
                        }

                        /* If it did not match_, record for a 2-hop matching. */
                        if maxidx == i
                            // && ivecaxpylez(ncon, 2, vwgt[(i * ncon)..], vwgt[(i * ncon)..], maxvwgt))
                            && util::ivecaxpylez(2, &vwgt[cntrng!(i * ncon, ncon)], &vwgt[cntrng!(i * ncon, ncon)], &maxvwgt[..ncon])
                        {
                            nunmatched += 1;
                            maxidx = (UNMATCHED as isize) as usize;
                        }
                    }
                }
            }

            if maxidx as idx_t != UNMATCHED {
                cmap[maxidx as usize] = cnvtxs;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
                match_[i as usize] = maxidx as idx_t;
                match_[maxidx as usize] = i as idx_t;
            }
        }
    }

    //println!("nunmatched: %zu", nunmatched);

    /* see if a 2-hop matching is required/allowed */
    if ctrl.no2hop == 0 && nunmatched > (UNMATCHEDFOR2HOP * nvtxs as real_t) as usize {
        // original bound this to cnvtxs, but that's redundant
        Match_2Hop(
            ctrl,
            graph,
            perm.as_mut_ptr(),
            match_.as_mut_ptr(),
            cnvtxs,
            nunmatched,
        );
    }

    /* match_ the final unmatched vertices with themselves and reorder the vertices
    of the coarse graph for memory-friendly contraction */
    cnvtxs = 0;
    for i in (0)..(nvtxs) {
        if match_[i as usize] == UNMATCHED {
            match_[i as usize] = i as idx_t;
            cmap[i as usize] = cnvtxs;
            cnvtxs += 1;
        } else {
            if i as idx_t <= match_[i as usize] {
                cmap[match_[i as usize] as usize] = cnvtxs;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
            }
        }
    }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.MatchTmr));

    CreateCoarseGraph(ctrl, graph, cnvtxs, match_.as_mut_ptr());

    // WCOREPOP;

    return cnvtxs;
}

/*************************************************************************/
/* This function matches the unmatched vertices using a 2-hop matching
that involves vertices that are two hops away from each other. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_2Hop(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    perm: *mut idx_t,
    match_: *mut idx_t,
    cnvtxs: idx_t,
    nunmatched: usize,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let mut cnvtxs = cnvtxs;
    let mut nunmatched = nunmatched;
    cnvtxs = Match_2HopAny(
        ctrl,
        graph,
        perm,
        match_,
        cnvtxs,
        &mut nunmatched as *mut usize,
        2,
    );
    cnvtxs = Match_2HopAll(
        ctrl,
        graph,
        perm,
        match_,
        cnvtxs,
        &mut nunmatched as *mut usize,
        64,
    );
    if (nunmatched as real_t) > 1.5 * UNMATCHEDFOR2HOP * graph.nvtxs as real_t {
        cnvtxs = Match_2HopAny(
            ctrl,
            graph,
            perm,
            match_,
            cnvtxs,
            &mut nunmatched as *mut usize,
            3,
        );
    }
    if (nunmatched as real_t) > 2.0 * UNMATCHEDFOR2HOP * graph.nvtxs as real_t {
        cnvtxs = Match_2HopAny(
            ctrl,
            graph,
            perm,
            match_,
            cnvtxs,
            &mut nunmatched as *mut usize,
            graph.nvtxs as usize,
        );
    }

    return cnvtxs;
}

/*************************************************************************/
/* This function matches the unmatched vertices whose degree is less than
maxdegree using a 2-hop matching that involves vertices that are two
hops away from each other.
The requirement of the 2-hop matching is a simple non-empty overlap
between the adjancency lists of the vertices. */
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_2HopAny(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    perm: *mut idx_t,
    match_: *mut idx_t,
    cnvtxs: idx_t,
    r_nunmatched: *mut usize,
    maxdegree: usize,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let mut cnvtxs = cnvtxs;
    // idx_t i, pi, ii, j, jj, k, nvtxs;
    // idx_t *xadj, *adjncy, *colptr, *rowind;
    // idx_t *cmap;
    // size_t nunmatched;

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.Aux3Tmr));

    get_graph_slices!(graph => xadj adjncy);
    get_graph_slices_mut!(graph => cmap);
    let nvtxs = graph.nvtxs as usize;
    mkslice_mut!(match_, nvtxs);
    mkslice!(perm, nvtxs);

    let mut nunmatched = *r_nunmatched;

    /*ifset!(ctrl.dbglvl, METIS_DBG_COARSEN, println!("IN: nunmatched: %zu\t", nunmatched)); */

    /* create the inverted index */
    // WCOREPUSH;
    let mut colptr = vec![0; nvtxs + 1];
    for i in (0)..(nvtxs) {
        if match_[i as usize] == UNMATCHED
            && xadj[(i + 1) as usize] - xadj[i as usize] < maxdegree as idx_t
        {
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                colptr[adjncy[j as usize] as usize] += 1;
            }
        }
    }
    util::make_csr(nvtxs, &mut colptr);

    let mut rowind = vec![0; colptr[nvtxs] as usize];
    for pi in (0)..(nvtxs) {
        let i = perm[pi as usize];
        if match_[i as usize] == UNMATCHED
            && xadj[(i + 1) as usize] - xadj[i as usize] < maxdegree as idx_t
        {
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                let p = &mut colptr[adjncy[j as usize] as usize];
                rowind[*p as usize] = i;
                *p += 1
            }
        }
    }
    util::shift_csr(nvtxs, &mut colptr);

    /* compute matchings by going down the inverted index */
    for pi in (0)..(nvtxs) {
        let i = perm[pi as usize];
        if colptr[(i + 1) as usize] - colptr[i as usize] < 2 {
            continue;
        }

        let mut jj = colptr[(i + 1) as usize];
        // for j in (colptr[i as usize])..(jj) {
        let mut j = colptr[i as usize];
        while j < jj {
            if match_[rowind[j as usize] as usize] == UNMATCHED {
                jj -= 1;
                // for (jj--; jj>j; jj--) {
                while j < jj {
                    if match_[rowind[jj as usize] as usize] == UNMATCHED {
                        cmap[rowind[jj as usize] as usize] = cnvtxs;
                        cmap[rowind[j as usize] as usize] = cnvtxs;
                        cnvtxs += 1;
                        match_[rowind[j as usize] as usize] = rowind[jj as usize];
                        match_[rowind[jj as usize] as usize] = rowind[j as usize];
                        nunmatched -= 2;
                        break;
                    }
                    jj -= 1;
                }
            }
            j += 1;
        }
    }
    // WCOREPOP;

    /*ifset!(ctrl.dbglvl, METIS_DBG_COARSEN, println!("OUT: nunmatched: %zu", nunmatched)); */

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.Aux3Tmr));

    *r_nunmatched = nunmatched;
    return cnvtxs;
}

/*************************************************************************/
/* This function matches the unmatched vertices whose degree is less than
   maxdegree using a 2-hop matching that involves vertices that are two
   hops away from each other.
   The requirement of the 2-hop matching is that of identical adjacency
   lists.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_2HopAll(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    perm: *mut idx_t,
    match_: *mut idx_t,
    cnvtxs: idx_t,
    r_nunmatched: *mut usize,
    maxdegree: usize,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let _ctrl = ctrl.as_mut().unwrap();
    // idx_t i, pi, pk, ii, j, jj, k, nvtxs, mask, idegree;
    // idx_t *xadj, *adjncy;
    // idx_t *cmap, *mark;
    // ikv_t *keys;
    // size_t nunmatched, ncand;

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.Aux3Tmr));

    let mut cnvtxs = cnvtxs;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy);
    get_graph_slices_mut!(graph => cmap);
    mkslice!(perm, nvtxs);
    mkslice_mut!(match_, nvtxs);

    let mut nunmatched = *r_nunmatched;
    let mask = idx_t::MAX / maxdegree as idx_t;

    /*ifset!(ctrl.dbglvl, METIS_DBG_COARSEN, println!("IN: nunmatched: %zu\t", nunmatched)); */

    // WCOREPUSH;
    #[derive(Clone, Copy, Default)]
    struct KeyVal {
        key: idx_t,
        val: idx_t,
    }

    /* collapse vertices with identical adjancency lists */
    // keys = ikvwspacemalloc(ctrl, nunmatched);
    let mut keys = vec![KeyVal::default(); nunmatched];
    let mut ncand = 0;
    for pi in (0)..(nvtxs) {
        let i = perm[pi as usize];
        let idegree = xadj[(i + 1) as usize] - xadj[i as usize];
        if match_[i as usize] == UNMATCHED && idegree > 1 && idegree < maxdegree as idx_t {
            let mut k = 0;
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                k += adjncy[j as usize] % mask;
            }
            keys[ncand as usize].val = i;
            keys[ncand as usize].key = (k % mask) * maxdegree as idx_t + idegree;
            ncand += 1;
        }
    }
    keys.sort_by_key(|k| k.key);
    // ikvsorti(ncand, keys);

    let mut mark = vec![0; nvtxs];
    for pi in (0)..(ncand) {
        let i = keys[pi as usize].val;
        if match_[i as usize] != UNMATCHED {
            continue;
        }

        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            mark[adjncy[j as usize] as usize] = i;
        }

        for pk in (pi + 1)..(ncand) {
            let k = keys[pk as usize].val;
            if match_[k as usize] != UNMATCHED {
                continue;
            }

            if keys[pi as usize].key != keys[pk as usize].key {
                break;
            }
            if xadj[(i + 1) as usize] - xadj[i as usize]
                != xadj[(k + 1) as usize] - xadj[k as usize]
            {
                break;
            }

            let mut jj = xadj[k as usize];
            while jj < (xadj[(k + 1) as usize]) {
                if mark[adjncy[jj as usize] as usize] != i {
                    break;
                }
                jj += 1;
            }
            if jj == xadj[(k + 1) as usize] {
                cmap[k as usize] = cnvtxs;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
                match_[i as usize] = k;
                match_[k as usize] = i;
                nunmatched -= 2;
                break;
            }
        }
    }
    // WCOREPOP;

    /*ifset!(ctrl.dbglvl, METIS_DBG_COARSEN, println!("OUT: ncand: %zu, nunmatched: %zu", ncand, nunmatched)); */

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.Aux3Tmr));

    *r_nunmatched = nunmatched;
    return cnvtxs;
}

/*************************************************************************/
/* This function finds a matching by selecting an adjacent vertex based
   on the Jaccard coefficient of the adjaceny lists.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_JC(ctrl: *mut ctrl_t, graph: *mut graph_t) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, pi, ii, iii, j, jj, jjj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx,
    // last_unmatched, avgdegree, bnum;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
    // idx_t *match_, *cmap, *degrees, *perm, *tperm, *vec, *marker;
    // idx_t mytwgt, xtwgt, ctwgt;
    // real_t bscore, score;

    // WCOREPUSH;

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.MatchTmr));

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj vwgt adjncy adjwgt);
    get_graph_slices_mut!(graph => cmap);
    let ncon = graph.ncon as usize;

    // let maxvwgt = ctrl.maxvwgt;
    mkslice!(ctrl->maxvwgt, ncon);

    let mut match_ = vec![UNMATCHED; nvtxs];
    let mut perm = vec![0; nvtxs];
    let mut tperm = vec![0; nvtxs];
    let mut degrees = vec![0; nvtxs];

    irandArrayPermute(nvtxs as idx_t, tperm.as_mut_ptr(), nvtxs as idx_t / 8, 1);

    let avgdegree = (4.0 * (xadj[nvtxs] / nvtxs as idx_t) as real_t) as idx_t;
    for i in (0)..(nvtxs) {
        let bnum = ((1 + xadj[(i + 1) as usize] - xadj[i as usize]) as real_t).sqrt() as idx_t;
        degrees[i as usize] = if bnum > avgdegree { avgdegree } else { bnum };
    }
    bucketsort::BucketSortKeysInc(
        ctrl,
        nvtxs as idx_t,
        avgdegree,
        degrees.as_ptr(),
        tperm.as_ptr(),
        perm.as_mut_ptr(),
    );

    /* point to the wspace vectors that are not needed any more */
    let mut vec = tperm;
    let mut marker = degrees;
    // iset(nvtxs, -1, vec);
    vec.fill(-1);
    // iset(nvtxs, -1, marker);
    marker.fill(-1);

    let mut cnvtxs = 0;
    let last_unmatched = 0;
    for pi in (0)..(nvtxs) {
        let i = perm[pi as usize] as usize;

        if match_[i as usize] == UNMATCHED {
            /* Unmatched */
            let mut maxidx = i as usize;

            if if ncon == 1 {
                vwgt[i as usize] < maxvwgt[0 as usize]
            } else {
                // ivecle(ncon, vwgt[(i * ncon)..], maxvwgt)
                util::ivecle(&vwgt[cntrng!(i * ncon, ncon)], &maxvwgt[..ncon])
            } {
                /* Deal with island vertices. Find a non-island and match_ it with.
                The matching ignores ctrl.maxvwgt requirements */
                if xadj[i as usize] == xadj[(i + 1) as usize] {
                    // last_unmatched = gk_max(pi, last_unmatched)+1;
                    for last_unmatched in (pi.max(last_unmatched) + 1)..(nvtxs) {
                        let j = perm[last_unmatched] as usize;
                        if match_[j] == UNMATCHED {
                            maxidx = j;
                            break;
                        }
                    }
                } else {
                    if ncon == 1 {
                        /* Find a max JC pair, subject to maxvwgt constraints */
                        if xadj[(i + 1) as usize] - xadj[i as usize] < avgdegree {
                            marker[i as usize] = i as idx_t;
                            let mut bscore = 0.0;
                            let mut mytwgt = 0;
                            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                                mytwgt += 1; //adjwgt[j as usize];
                                vec[adjncy[j as usize] as usize] = 1; //adjwgt[j as usize];
                            }

                            /* single constraint pairing */
                            // #ifdef XXX
                            //               for j in (xadj[i as usize])..(xadj[(i+1) as usize]) {
                            //                 ii = adjncy[j as usize];
                            //                 if (marker[ii as usize] == i || match_[ii as usize] != UNMATCHED || vwgt[i as usize]+vwgt[ii as usize] > maxvwgt[0 as usize])
                            //         {
                            //             continue;
                            //         }
                            //
                            //                 ctwgt = xtwgt = 0;
                            //                 for jj in (xadj[ii as usize])..(xadj[(ii+1) as usize]) {
                            //                   xtwgt += adjwgt[jj as usize];
                            //                   if (vec[adjncy[jj as usize] as usize] > 0)
                            //         {
                            //             ctwgt += vec[adjncy[jj as usize] as usize] + adjwgt[jj as usize];
                            //         }
                            //                   else if (adjncy[jj as usize] == i) {
                            //                     ctwgt += adjwgt[jj as usize];
                            //                     xtwgt -= adjwgt[jj as usize];
                            //                   }
                            //                 }
                            //
                            //                 score = 1.0*ctwgt/(mytwgt+xtwgt-ctwgt);
                            //                 if (score > bscore) {
                            //                   bscore = score;
                            //                   maxidx = ii;
                            //                 }
                            //                 marker[ii as usize] = i;
                            //               }
                            // #endif

                            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                                let ii = adjncy[j as usize];
                                for jj in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
                                    let iii = adjncy[jj as usize] as usize;

                                    if marker[iii] == i as idx_t
                                        || match_[iii] != UNMATCHED
                                        || vwgt[i as usize] + vwgt[iii] > maxvwgt[0 as usize]
                                    {
                                        continue;
                                    }

                                    let mut ctwgt = 0;
                                    let mut xtwgt = 0;
                                    for jjj in (xadj[iii])..(xadj[(iii + 1) as usize]) {
                                        xtwgt += 1; //adjwgt[jjj as usize];
                                        if vec[adjncy[jjj as usize] as usize] > 0 {
                                            ctwgt += 2; //vec[adjncy[jjj as usize] as usize] + adjwgt[jjj as usize];
                                        } else if adjncy[jjj as usize] == i as idx_t {
                                            ctwgt += 10 * adjwgt[jjj as usize];
                                        }
                                    }

                                    let score = 1.0 * (ctwgt / (mytwgt + xtwgt)) as real_t;
                                    //println!("{:} {:} {:} %.4f", mytwgt, xtwgt, ctwgt, score);
                                    if score > bscore {
                                        bscore = score;
                                        maxidx = iii;
                                    }
                                    marker[iii] = i as idx_t;
                                }
                            }

                            /* reset vec array */
                            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                                vec[adjncy[j as usize] as usize] = -1;
                            }
                        }
                    } else {
                        /* multi-constraint version */
                        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                            let k = adjncy[j as usize] as usize;
                            if match_[k as usize] == UNMATCHED
                                // && ivecaxpylez(ncon, 1, vwgt[(i * ncon)..], vwgt[(k * ncon)..], maxvwgt))
                                && util::ivecaxpylez(1, &vwgt[cntrng!(i * ncon, ncon)], &vwgt[cntrng!(k * ncon, ncon)], maxvwgt)
                            {
                                maxidx = k;
                                break;
                            }
                        }
                    }
                }
            }

            if maxidx as idx_t != UNMATCHED {
                cmap[maxidx] = cnvtxs;
                cmap[i] = cnvtxs;
                cnvtxs += 1;
                match_[i] = maxidx as idx_t;
                match_[maxidx] = i as idx_t;
            }
        }
    }

    /* match_ the final unmatched vertices with themselves and reorder the vertices
    of the coarse graph for memory-friendly contraction */
    cnvtxs = 0;
    for i in (0)..(nvtxs) {
        if match_[i as usize] == UNMATCHED {
            match_[i as usize] = i as idx_t;
            cmap[i as usize] = cnvtxs;
            cnvtxs += 1;
        } else {
            if i as idx_t <= match_[i as usize] {
                cmap[match_[i as usize] as usize] = cnvtxs;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
            }
        }
    }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.MatchTmr));

    CreateCoarseGraph(ctrl, graph, cnvtxs, match_.as_mut_ptr());

    // WCOREPOP;

    return cnvtxs;
}

/*************************************************************************/
/* This function prints various stats for each graph during coarsening
 */
/*************************************************************************/
#[metis_func]
pub extern "C" fn PrintCGraphStats(ctrl: *mut ctrl_t, graph: *mut graph_t) -> () {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i;
    get_graph_slices!(graph => adjwgt);
    println!(
        "{:10} {:10} {:10} [{:}] [",
        graph.nvtxs,
        graph.nedges,
        adjwgt[..(graph.nedges as usize)].iter().sum::<idx_t>(),
        ctrl.CoarsenTo,
    );

    for i in (0)..(graph.ncon) {
        println!(
            " {:8}:{:8}",
            *ctrl.maxvwgt.add(i as usize),
            *graph.tvwgt.add(i as usize)
        );
    }
    println!(" ]");
}

/*************************************************************************/
/* This function creates the coarser graph. Depending on the size of the
   candidate adjancency lists it either uses a hash table or an array
   to do duplicate detection.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn CreateCoarseGraph(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    cnvtxs: idx_t,
    match_: *mut idx_t,
) -> () {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon,
    // cnedges, v, u, mask;
    // idx_t *xadj, *vwgt, *vsize, *adjncy, *adjwgt;
    // idx_t *cmap, *htable, *dtable;
    // idx_t *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt;
    // graph_t *cgraph;
    // int dovsize, dropedges;
    // idx_t cv, nkeys, droppedewgt;
    // idx_t *keys=NULL, *medianewgts=NULL, *noise=NULL;

    // WCOREPUSH;

    let dovsize = ctrl.objtype == METIS_OBJTYPE_VOL;
    let cnvtxs = cnvtxs as usize;
    let dropedges = ctrl.dropedges;

    let mask = HTLENGTH;

    // ifset!(
    // ctrl.dbglvl,
    // METIS_DBG_TIME,
    // gk_startcputimer(ctrl.ContractTmr),
    // );

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    // xadj = graph.xadj;
    // vwgt = graph.vwgt;
    // vsize = graph.vsize;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // cmap = graph.cmap;

    get_graph_slices!(ctrl, graph => xadj vwgt adjncy adjwgt cmap);
    mkslice_mut!(match_, nvtxs);

    let vsize = if dovsize {
        get_graph_slices!(ctrl, graph => vsize);
        vsize
    } else {
        &[]
    };

    /* Setup structures for dropedges */
    let mut nkeys = 0;
    // these 3 vars are only used when dropedges != 0, but I hope the optimizer will be good
    let mut medianewgts = vec![-1; cnvtxs];
    let mut keys = vec![0; nkeys as usize];
    let mut noise = vec![0; cnvtxs];
    if dropedges != 0 {
        nkeys = 0;
        for v in (0)..(nvtxs) {
            nkeys = nkeys.max(xadj[(v + 1) as usize] - xadj[v as usize]);
        }
        nkeys = 2 * nkeys + 1;

        for v in (0)..(cnvtxs) {
            noise[v as usize] = irandInRange(128);
        }
    }

    /* Initialize the coarser graph */
    debug_assert!(cnvtxs > 0);
    let cgraph = SetupCoarseGraph(graph, cnvtxs as idx_t, dovsize as std::ffi::c_int);
    mkslice_mut!(cxadj: cgraph->xadj, cnvtxs + 1);
    mkslice_mut!(cvwgt: cgraph->vwgt, (*cgraph).ncon as usize * cnvtxs);
    mkslice_mut!(cadjncy: cgraph->adjncy, graph.nedges + 1);
    mkslice_mut!(cadjwgt: cgraph->adjwgt, graph.nedges + 1);
    let cvsize = if dovsize {
        std::slice::from_raw_parts_mut((*cgraph).vsize, cnvtxs)
    } else {
        &mut []
    };
    let mut cadjncy = cadjncy;
    let mut cadjwgt = cadjwgt;
    // cxadj = cgraph.xadj;
    // cvwgt = cgraph.vwgt;
    // cvsize = cgraph.vsize;
    // cadjncy = cgraph.adjncy;
    // cadjwgt = cgraph.adjwgt;

    let mut htable = vec![-1 as idx_t; mask as usize + 1]; /* hash table */
    let mut dtable = vec![-1 as idx_t; cnvtxs]; /* direct table */

    cxadj[0 as usize] = 0;
    let mut cnvtxs: usize = 0;
    let mut cnedges = 0;
    let mut nedges;
    for v in (0)..(nvtxs) {
        let u = match_[v as usize] as usize;
        if u < v {
            continue;
        }

        assert!(cmap[v as usize] == cnvtxs as idx_t);
        assert!(cmap[match_[v as usize] as usize] == cnvtxs as idx_t);

        /* take care of the vertices */
        if ncon == 1 {
            cvwgt[cnvtxs] = vwgt[v as usize];
        } else {
            // icopy(ncon, vwgt[(v * ncon)..], cvwgt[( cnvtxs * ncon)..]);
            cvwgt[cntrng!(cnvtxs * ncon, ncon)].copy_from_slice(&vwgt[cntrng!(v * ncon, ncon)]);
        }

        if dovsize {
            cvsize[cnvtxs] = vsize[v as usize];
        }

        if v != u {
            if ncon == 1 {
                cvwgt[cnvtxs] += vwgt[u as usize];
            } else {
                blas::iaxpy(
                    ncon,
                    1,
                    &vwgt[(u * ncon)..],
                    1,
                    &mut cvwgt[(cnvtxs * ncon)..],
                    1,
                );
            }

            if dovsize {
                cvsize[cnvtxs] += vsize[u as usize];
            }
        }

        /* take care of the edges */
        if (xadj[(v + 1) as usize] - xadj[v as usize] + xadj[(u + 1) as usize] - xadj[u as usize])
            < (mask >> 2)
        {
            /* use mask */
            /* put the ID of the contracted node itself at the start, so that it can be
             * removed easily */
            htable[cnvtxs & mask as usize] = 0;
            cadjncy[0 as usize] = cnvtxs as idx_t;
            nedges = 1;

            {
                let istart = xadj[v as usize] as usize;
                let iend = xadj[(v + 1) as usize] as usize;
                for j in (istart)..(iend) {
                    let k = cmap[adjncy[j as usize] as usize];
                    // for (kk=k&mask; htable[kk as usize]!=-1 && cadjncy[htable[kk as usize] as usize]!=k; kk=((kk+1)&mask));
                    let mut kk = k & mask;
                    while htable[kk as usize] != -1 && cadjncy[htable[kk as usize] as usize] != k {
                        kk = (kk + 1) & mask
                    }
                    let m = htable[kk as usize];
                    if m == -1 {
                        cadjncy[nedges as usize] = k;
                        cadjwgt[nedges as usize] = adjwgt[j as usize];
                        htable[kk as usize] = nedges;
                        nedges += 1;
                    } else {
                        cadjwgt[m as usize] += adjwgt[j as usize];
                    }
                }
            }

            if v != u {
                let istart = xadj[u as usize] as usize;
                let iend = xadj[(u + 1) as usize] as usize;
                for j in (istart)..(iend) {
                    let k = cmap[adjncy[j as usize] as usize];
                    // for (kk=k&mask; htable[kk as usize]!=-1 && cadjncy[htable[kk as usize] as usize]!=k; kk=((kk+1)&mask));
                    let mut kk = k & mask;
                    while htable[kk as usize] != -1 && cadjncy[htable[kk as usize] as usize] != k {
                        kk = (kk + 1) & mask
                    }
                    let m = htable[kk as usize];
                    if m == -1 {
                        cadjncy[nedges as usize] = k;
                        cadjwgt[nedges as usize] = adjwgt[j as usize];
                        htable[kk as usize] = nedges;
                        nedges += 1;
                    } else {
                        cadjwgt[m as usize] += adjwgt[j as usize];
                    }
                }
            }

            /* reset the htable -- reverse order (LIFO) is critical to prevent cadjncy[-1]
             * indexing due to a remove of an earlier entry */
            // for (j=nedges-1; j>=0; j--) {
            for j in (0..nedges).rev() {
                let k = cadjncy[j as usize];
                // for (kk=k&mask; cadjncy[htable[kk as usize] as usize]!=k; kk=((kk+1)&mask));
                let mut kk = k & mask;
                while cadjncy[htable[kk as usize] as usize] != k {
                    kk = (kk + 1) & mask;
                }
                htable[kk as usize] = -1;
                // TODO: Verify this
            }

            /* remove the contracted vertex from the list */
            nedges -= 1;
            cadjncy[0 as usize] = cadjncy[nedges as usize];
            cadjwgt[0 as usize] = cadjwgt[nedges as usize];
        } else {
            nedges = 0;
            let mut istart = xadj[v as usize];
            let mut iend = xadj[(v + 1) as usize];
            for j in (istart)..(iend) {
                let k = cmap[adjncy[j as usize] as usize];
                let m = dtable[k as usize];
                if m == -1 {
                    cadjncy[nedges as usize] = k;
                    cadjwgt[nedges as usize] = adjwgt[j as usize];
                    dtable[k as usize] = nedges;
                    nedges += 1;
                } else {
                    cadjwgt[m as usize] += adjwgt[j as usize];
                }
            }

            if v != u {
                istart = xadj[u as usize];
                iend = xadj[(u + 1) as usize];
                for j in (istart)..(iend) {
                    let k = cmap[adjncy[j as usize] as usize];
                    let m = dtable[k as usize];
                    if m == -1 {
                        cadjncy[nedges as usize] = k;
                        cadjwgt[nedges as usize] = adjwgt[j as usize];
                        dtable[k as usize] = nedges;
                        nedges += 1;
                    } else {
                        cadjwgt[m as usize] += adjwgt[j as usize];
                    }
                }

                /* Remove the contracted self-loop, when present */
                let j = dtable[cnvtxs];
                if j != -1 {
                    assert!(cadjncy[j as usize] == cnvtxs as idx_t);
                    nedges -= 1;
                    cadjncy[j as usize] = cadjncy[nedges as usize];
                    cadjwgt[j as usize] = cadjwgt[nedges as usize];
                    dtable[cnvtxs] = -1;
                }
            }

            /* Zero out the dtable */
            for j in (0)..(nedges) {
                dtable[cadjncy[j as usize] as usize] = -1;
            }
        }

        /* Determine the median weight of the incident edges, which will be used
        to keep an edge (u, v) iff wgt(u, v) >= min(medianewgts[u as usize], medianewgts[v as usize]) */
        if dropedges != 0 {
            assert!(nedges < nkeys, "{:}, {:}", nkeys, nedges);
            medianewgts[cnvtxs] = 8; /* default for island nodes */
            if nedges > 0 {
                for j in (0)..(nedges) {
                    keys[j as usize] = (cadjwgt[j as usize] << 8)
                        + noise[cnvtxs]
                        + noise[cadjncy[j as usize] as usize];
                }
                // isortd(nedges, keys);
                keys.sort_unstable();
                medianewgts[cnvtxs] = keys[(nedges - 1).min(
                    (xadj[(v + 1) as usize] - xadj[v as usize] + xadj[(u + 1) as usize]
                        - xadj[u as usize])
                        >> 1,
                ) as usize];
            }
        }

        cadjncy = &mut cadjncy[nedges as usize..];
        cadjwgt = &mut cadjwgt[nedges as usize..];
        cnedges += nedges;
        cnvtxs += 1;
        cxadj[cnvtxs] = cnedges;
    }

    /* compact the adjacency structure of the coarser graph to keep only +ve edges */
    if dropedges != 0 {
        let mut droppedewgt = 0;

        // cadjncy = cgraph.adjncy;
        // cadjwgt = cgraph.adjwgt;
        mkslice_mut!(cadjncy: cgraph->adjncy, cxadj[cnvtxs]);
        mkslice_mut!(cadjwgt: cgraph->adjwgt, cxadj[cnvtxs]);

        cnedges = 0;
        for u in (0)..(cnvtxs) {
            let istart = cxadj[u as usize];
            let iend = cxadj[(u + 1) as usize];
            for j in (istart)..(iend) {
                let v = cadjncy[j as usize];
                assert!(
                    medianewgts[u as usize] >= 0,
                    "{:} {:}",
                    u,
                    medianewgts[u as usize]
                );
                assert!(
                    medianewgts[v as usize] >= 0,
                    "{:} {:} {:}",
                    v,
                    medianewgts[v as usize],
                    cnvtxs,
                );
                if (cadjwgt[j as usize] << 8) + noise[u as usize] + noise[v as usize]
                    >= medianewgts[u as usize].min(medianewgts[v as usize])
                {
                    cadjncy[cnedges as usize] = cadjncy[j as usize];
                    cadjwgt[cnedges as usize] = cadjwgt[j as usize];
                    cnedges += 1;
                } else {
                    droppedewgt += cadjwgt[j as usize];
                }
            }
            cxadj[u as usize] = cnedges;
        }
        util::shift_csr(cnvtxs, cxadj);

        (*cgraph).droppedewgt = droppedewgt;
    }

    (*cgraph).nedges = cnedges;

    for j in (0)..(ncon) {
        mkslice_mut!(cgraph->tvwgt, ncon);
        mkslice_mut!(cgraph->vwgt, ncon * (*cgraph).nvtxs as usize);
        mkslice_mut!(cgraph->invtvwgt, ncon);
        // tvwgt[j as usize] = isum(cgraph.nvtxs, cgraph.vwgt + j, ncon);
        tvwgt[j] = vwgt[j..].iter().step_by(ncon).sum();
        invtvwgt[j as usize] = 1.0
            / (if tvwgt[j as usize] > 0 {
                tvwgt[j as usize] as real_t
            } else {
                1.0
            });
    }

    ReAdjustMemory(ctrl, graph, cgraph);

    // ifset!(
    // ctrl.dbglvl,
    // METIS_DBG_TIME,
    // gk_stopcputimer(ctrl.ContractTmr),
    // );

    // WCOREPOP;
}

/*************************************************************************/
/* Setup the various arrays for the coarse graph
 */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SetupCoarseGraph(
    graph: *mut graph_t,
    cnvtxs: idx_t,
    dovsize: std::ffi::c_int,
) -> *mut graph_t {
    let graph = graph.as_mut().unwrap();
    let dovsize = dovsize != 0;
    // graph_t *cgraph;

    let cgraph = graph::CreateGraph();
    let cgraph = cgraph.as_mut().unwrap();

    cgraph.nvtxs = cnvtxs;
    cgraph.ncon = graph.ncon;

    cgraph.finer = graph;
    graph.coarser = cgraph;

    /* Allocate memory for the coarser graph.
    NOTE: The +1 in the adjwgt/adjncy is to allow the optimization of self-loop
          detection by adding ahead of time the self-loop. That optimization
          requires a +1 adjncy/adjwgt array for the limit case where_ the
          coarser graph is of the same size of the previous graph. */
    cgraph.xadj = imalloc(cnvtxs as usize + 1, c"SetupCoarseGraph: xadj".as_ptr()) as _;
    cgraph.adjncy = imalloc(
        graph.nedges as usize + 1,
        c"SetupCoarseGraph: adjncy".as_ptr(),
    ) as _;
    cgraph.adjwgt = imalloc(
        graph.nedges as usize + 1,
        c"SetupCoarseGraph: adjwgt".as_ptr(),
    ) as _;
    cgraph.vwgt = imalloc(
        cgraph.ncon as usize * cnvtxs as usize,
        c"SetupCoarseGraph: vwgt".as_ptr(),
    ) as _;
    cgraph.tvwgt = imalloc(cgraph.ncon as usize, c"SetupCoarseGraph: tvwgt".as_ptr()) as _;
    cgraph.invtvwgt = rmalloc(cgraph.ncon as usize, c"SetupCoarseGraph: invtvwgt".as_ptr()) as _;

    if dovsize {
        cgraph.vsize = imalloc(cnvtxs as usize, c"SetupCoarseGraph: vsize".as_ptr()) as _;
    }

    return cgraph;
}

/*************************************************************************/
/* This function re-adjusts the amount of memory that was allocated if
   it will lead to significant savings
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn ReAdjustMemory(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    cgraph: *mut graph_t,
) -> () {
    let graph = graph.as_mut().unwrap();
    let cgraph = cgraph.as_mut().unwrap();
    if cgraph.nedges > 10000 && (cgraph.nedges as real_t) < 0.9 * graph.nedges as real_t {
        cgraph.adjncy = irealloc(
            cgraph.adjncy as _,
            cgraph.nedges as usize,
            c"ReAdjustMemory: adjncy".as_ptr(),
        ) as _;
        cgraph.adjwgt = irealloc(
            cgraph.adjwgt as _,
            cgraph.nedges as usize,
            c"ReAdjustMemory: adjwgt".as_ptr(),
        ) as _;
    }
}
