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

use crate::*;

const UNMATCHEDFOR2HOP: f32 = 0.10; /* The fraction of unmatched vertices that triggers 2-hop */

/*************************************************************************/
/* This function takes a graph and creates a sequence of coarser graphs.
   It implements the coarsening phase of the multilevel paradigm.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn CoarsenGraph(ctrl: *mut ctrl_t, graph: *mut graph_t) -> *mut graph_t {
    // idx_t i, eqewgts, level=0;

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_TIME,
        gk_startcputimer(ctrl.CoarsenTmr),
    );

    /* determine if the weights on the edges are all the same */
    eqewgts = 1;
    for i in (1)..(graph.nedges) {
        if (graph.adjwgt[0] != graph.adjwgt[i]) {
            eqewgts = 0;
            break;
        }
    }

    /* set the maximum allowed coarsest vertex weight */
    for i in (0)..(graph.ncon) {
        ctrl.maxvwgt[i] = 1.5 * graph.tvwgt[i] / ctrl.CoarsenTo;
    }

    loop {
        IFSET(
            ctrl.dbglvl,
            METIS_DBG_COARSEN,
            PrintCGraphStats(ctrl, graph),
        );

        /* allocate memory for cmap, if it has not already been done due to
        multiple cuts */
        if (graph.cmap == NULL) {
            graph.cmap = imalloc(graph.nvtxs, "CoarsenGraph: graph.cmap");
        }

        /* determine which matching scheme you will use */
        match (ctrl.ctype) {
            METIS_CTYPE_RM => Match_RM(ctrl, graph),
            METIS_CTYPE_SHEM => {
                if (eqewgts || graph.nedges == 0) {
                    Match_RM(ctrl, graph);
                } else {
                    Match_SHEM(ctrl, graph);
                }
            }
            _ => gk_errexit(SIGERR, "Unknown ctype: %d\n", ctrl.ctype),
        }

        graph_WriteToDisk(ctrl, graph);

        graph = graph.coarser;
        eqewgts = 0;
        level;
        level += 1;

        ASSERT(CheckGraph(graph, 0, 1));

        if (!(graph.nvtxs > ctrl.CoarsenTo
            && graph.nvtxs < COARSEN_FRACTION * graph.finer.nvtxs
            && graph.nedges > graph.nvtxs / 2))
        {
            break;
        }
    }

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_COARSEN,
        PrintCGraphStats(ctrl, graph),
    );
    IFSET(
        ctrl.dbglvl,
        METIS_DBG_TIME,
        gk_stopcputimer(ctrl.CoarsenTmr),
    );

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
    // idx_t i, eqewgts, level;

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_TIME,
        gk_startcputimer(ctrl.CoarsenTmr),
    );

    /* determine if the weights on the edges are all the same */
    eqewgts = 1;
    for i in (1)..(graph.nedges) {
        if (graph.adjwgt[0] != graph.adjwgt[i]) {
            eqewgts = 0;
            break;
        }
    }

    /* set the maximum allowed coarsest vertex weight */
    for i in (0)..(graph.ncon) {
        ctrl.maxvwgt[i] = 1.5 * graph.tvwgt[i] / ctrl.CoarsenTo;
    }

    for level in (0)..(nlevels) {
        IFSET(
            ctrl.dbglvl,
            METIS_DBG_COARSEN,
            PrintCGraphStats(ctrl, graph),
        );

        /* allocate memory for cmap, if it has not already been done due to
        multiple cuts */
        if (graph.cmap == NULL) {
            graph.cmap = imalloc(graph.nvtxs, "CoarsenGraph: graph.cmap");
        }

        /* determine which matching scheme you will use */
        match (ctrl.ctype) {
            METIS_CTYPE_RM => Match_RM(ctrl, graph),
            METIS_CTYPE_SHEM => {
                if (eqewgts || graph.nedges == 0) {
                    Match_RM(ctrl, graph);
                } else {
                    Match_SHEM(ctrl, graph);
                }
            }
            _ => gk_errexit(SIGERR, "Unknown ctype: %d\n", ctrl.ctype),
        }

        graph_WriteToDisk(ctrl, graph);

        graph = graph.coarser;
        eqewgts = 0;

        ASSERT(CheckGraph(graph, 0, 1));

        if (graph.nvtxs < ctrl.CoarsenTo
            || graph.nvtxs > COARSEN_FRACTION * graph.finer.nvtxs
            || graph.nedges < graph.nvtxs / 2)
        {
            break;
        }
    }

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_COARSEN,
        PrintCGraphStats(ctrl, graph),
    );
    IFSET(
        ctrl.dbglvl,
        METIS_DBG_TIME,
        gk_stopcputimer(ctrl.CoarsenTmr),
    );

    return graph;
}

/*************************************************************************/
/* This function finds a matching by randomly selecting one of the
   unmatched adjacent vertices.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_RM(ctrl: *mut ctrl_t, graph: *mut graph_t) -> idx_t {
    // idx_t i, pi, ii, j, jj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx,
    // last_unmatched, avgdegree, bnum;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
    // idx_t *match__, *cmap, *degrees, *perm, *tperm;
    // size_t nunmatched=0;

    WCOREPUSH;

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.MatchTmr));

    nvtxs = graph.nvtxs;
    ncon = graph.ncon;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;
    cmap = graph.cmap;

    maxvwgt = ctrl.maxvwgt;

    match__ = iset(nvtxs, UNMATCHED, iwspacemalloc(ctrl, nvtxs));
    perm = iwspacemalloc(ctrl, nvtxs);
    tperm = iwspacemalloc(ctrl, nvtxs);
    degrees = iwspacemalloc(ctrl, nvtxs);

    /* Determine a "random" traversal order that is biased towards
    low-degree vertices */
    irandArrayPermute(nvtxs, tperm, nvtxs / 8, 1);

    avgdegree = 4.0 * (xadj[nvtxs] / nvtxs);
    for i in (0)..(nvtxs) {
        bnum = sqrt(1 + xadj[i + 1] - xadj[i]);
        degrees[i] = (if bnum > avgdegree { avgdegree } else { bnum });
    }
    BucketSortKeysInc(ctrl, nvtxs, avgdegree, degrees, tperm, perm);

    /* Traverse the vertices and compute the matching */
    cnvtxs = 0;
    last_unmatched = 0;
    for pi in (0)..(nvtxs) {
        i = perm[pi];

        if (match__[i] == UNMATCHED) {
            /* Unmatched */
            maxidx = i;

            if (if (ncon == 1) {
                vwgt[i] < maxvwgt[0]
            } else {
                ivecle(ncon, vwgt + i * ncon, maxvwgt)
            }) {
                /* Deal with island vertices. Find a non-island and match__ it with.
                The matching ignores ctrl.maxvwgt requirements */
                if (xadj[i] == xadj[i + 1]) {
                    // last_unmatched=gk_max(pi, last_unmatched)+1;
                    for last_unmatched in (gk_max(pi, last_unmatched) + 1)..(nvtxs) {
                        j = perm[last_unmatched];
                        if (match__[j] == UNMATCHED) {
                            maxidx = j;
                            break;
                        }
                    }
                } else {
                    /* Find a random matching, subject to maxvwgt constraints */
                    if (ncon == 1) {
                        /* single constraint version */
                        for j in (xadj[i])..(xadj[i + 1]) {
                            k = adjncy[j];
                            if (match__[k] == UNMATCHED && vwgt[i] + vwgt[k] <= maxvwgt[0]) {
                                maxidx = k;
                                break;
                            }
                        }

                        /* If it did not match__, record for a 2-hop matching. */
                        if (maxidx == i && 2 * vwgt[i] < maxvwgt[0]) {
                            nunmatched;
                            nunmatched += 1;
                            maxidx = UNMATCHED;
                        }
                    } else {
                        /* multi-constraint version */
                        for j in (xadj[i])..(xadj[i + 1]) {
                            k = adjncy[j];
                            if (match__[k] == UNMATCHED
                                && ivecaxpylez(ncon, 1, vwgt + i * ncon, vwgt + k * ncon, maxvwgt))
                            {
                                maxidx = k;
                                break;
                            }
                        }

                        /* If it did not match__, record for a 2-hop matching. */
                        if (maxidx == i
                            && ivecaxpylez(ncon, 2, vwgt + i * ncon, vwgt + i * ncon, maxvwgt))
                        {
                            nunmatched;
                            nunmatched += 1;
                            maxidx = UNMATCHED;
                        }
                    }
                }
            }

            if (maxidx != UNMATCHED) {
                cmap[i] = cmap[maxidx] = cnvtxs;
                cnvtxs += 1;
                match__[i] = maxidx;
                match__[maxidx] = i;
            }
        }
    }

    //printf("nunmatched: %zu\n", nunmatched);

    /* see if a 2-hop matching is required/allowed */
    if (!ctrl.no2hop && nunmatched > UNMATCHEDFOR2HOP * nvtxs) {
        cnvtxs = Match_2Hop(ctrl, graph, perm, match__, cnvtxs, nunmatched);
    }

    /* match__ the final unmatched vertices with themselves and reorder the vertices
    of the coarse graph for memory-friendly contraction */
    cnvtxs = 0;
    for i in (0)..(nvtxs) {
        if (match__[i] == UNMATCHED) {
            match__[i] = i;
            cmap[i] = cnvtxs;
            cnvtxs += 1;
        } else {
            if (i <= match__[i]) {
                cmap[i] = cmap[match__[i]] = cnvtxs;
                cnvtxs += 1;
            }
        }
    }

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.MatchTmr));

    CreateCoarseGraph(ctrl, graph, cnvtxs, match__);

    WCOREPOP;

    return cnvtxs;
}

/**************************************************************************/
/* This function finds a matching using the HEM heuristic. The vertices
   are visited based on increasing degree to ensure that all vertices are
   given a chance to match__ with something.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Match_SHEM(ctrl: *mut ctrl_t, graph: *mut graph_t) -> idx_t {
    // idx_t i, pi, ii, j, jj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx, maxwgt,
    // last_unmatched, avgdegree, bnum;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
    // idx_t *match__, *cmap, *degrees, *perm, *tperm;
    // size_t nunmatched=0;

    WCOREPUSH;

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.MatchTmr));

    nvtxs = graph.nvtxs;
    ncon = graph.ncon;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;
    cmap = graph.cmap;

    maxvwgt = ctrl.maxvwgt;

    match__ = iset(nvtxs, UNMATCHED, iwspacemalloc(ctrl, nvtxs));
    perm = iwspacemalloc(ctrl, nvtxs);
    tperm = iwspacemalloc(ctrl, nvtxs);
    degrees = iwspacemalloc(ctrl, nvtxs);

    /* Determine a "random" traversal order that is biased towards low-degree vertices */
    irandArrayPermute(nvtxs, tperm, nvtxs / 8, 1);

    avgdegree = 4.0 * (xadj[nvtxs] / nvtxs);
    for i in (0)..(nvtxs) {
        bnum = sqrt(1 + xadj[i + 1] - xadj[i]);
        degrees[i] = (if bnum > avgdegree { avgdegree } else { bnum });
    }
    BucketSortKeysInc(ctrl, nvtxs, avgdegree, degrees, tperm, perm);

    /* Traverse the vertices and compute the matching */
    cnvtxs = 0;
    last_unmatched = 0;
    for pi in (0)..(nvtxs) {
        i = perm[pi];

        if (match__[i] == UNMATCHED) {
            /* Unmatched */
            maxidx = i;
            maxwgt = -1;

            if (if ncon == 1 {
                vwgt[i] < maxvwgt[0]
            } else {
                ivecle(ncon, vwgt + i * ncon, maxvwgt)
            }) {
                /* Deal with island vertices. Find a non-island and match__ it with.
                The matching ignores ctrl.maxvwgt requirements */
                if (xadj[i] == xadj[i + 1]) {
                    // last_unmatched = gk_max(pi, last_unmatched)+1;
                    for last_unmatched in (gk_max(pi, last_unmatched) + 1)..(nvtxs) {
                        j = perm[last_unmatched];
                        if (match__[j] == UNMATCHED) {
                            maxidx = j;
                            break;
                        }
                    }
                } else {
                    /* Find a heavy-edge matching, subject to maxvwgt constraints */
                    if (ncon == 1) {
                        /* single constraint version */
                        for j in (xadj[i])..(xadj[i + 1]) {
                            k = adjncy[j];
                            if (match__[k] == UNMATCHED
                                && maxwgt < adjwgt[j]
                                && vwgt[i] + vwgt[k] <= maxvwgt[0])
                            {
                                maxidx = k;
                                maxwgt = adjwgt[j];
                            }
                        }

                        /* If it did not match__, record for a 2-hop matching. */
                        if (maxidx == i && 2 * vwgt[i] < maxvwgt[0]) {
                            nunmatched;
                            nunmatched += 1;
                            maxidx = UNMATCHED;
                        }
                    } else {
                        /* multi-constraint version */
                        for j in (xadj[i])..(xadj[i + 1]) {
                            k = adjncy[j];
                            if (match__[k] == UNMATCHED
                                && ivecaxpylez(ncon, 1, vwgt + i * ncon, vwgt + k * ncon, maxvwgt)
                                && (maxwgt < adjwgt[j]
                                    || (maxwgt == adjwgt[j]
                                        && BetterVBalance(
                                            ncon,
                                            graph.invtvwgt,
                                            vwgt + i * ncon,
                                            vwgt + maxidx * ncon,
                                            vwgt + k * ncon,
                                        ))))
                            {
                                maxidx = k;
                                maxwgt = adjwgt[j];
                            }
                        }

                        /* If it did not match__, record for a 2-hop matching. */
                        if (maxidx == i
                            && ivecaxpylez(ncon, 2, vwgt + i * ncon, vwgt + i * ncon, maxvwgt))
                        {
                            nunmatched;
                            nunmatched += 1;
                            maxidx = UNMATCHED;
                        }
                    }
                }
            }

            if (maxidx != UNMATCHED) {
                cmap[i] = cmap[maxidx] = cnvtxs;
                cnvtxs += 1;
                match__[i] = maxidx;
                match__[maxidx] = i;
            }
        }
    }

    //printf("nunmatched: %zu\n", nunmatched);

    /* see if a 2-hop matching is required/allowed */
    if (!ctrl.no2hop && nunmatched > UNMATCHEDFOR2HOP * nvtxs) {
        cnvtxs = Match_2Hop(ctrl, graph, perm, match__, cnvtxs, nunmatched);
    }

    /* match__ the final unmatched vertices with themselves and reorder the vertices
    of the coarse graph for memory-friendly contraction */
    cnvtxs = 0;
    for i in (0)..(nvtxs) {
        if (match__[i] == UNMATCHED) {
            match__[i] = i;
            cmap[i] = cnvtxs;
            cnvtxs += 1;
        } else {
            if (i <= match__[i]) {
                cmap[i] = cmap[match__[i]] = cnvtxs;
                cnvtxs += 1;
            }
        }
    }

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.MatchTmr));

    CreateCoarseGraph(ctrl, graph, cnvtxs, match__);

    WCOREPOP;

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
    match__: *mut idx_t,
    cnvtxs: idx_t,
    nunmatched: size_t,
) -> idx_t {
    cnvtxs = Match_2HopAny(ctrl, graph, perm, match__, cnvtxs, &nunmatched, 2);
    cnvtxs = Match_2HopAll(ctrl, graph, perm, match__, cnvtxs, &nunmatched, 64);
    if (nunmatched > 1.5 * UNMATCHEDFOR2HOP * graph.nvtxs) {
        cnvtxs = Match_2HopAny(ctrl, graph, perm, match__, cnvtxs, &nunmatched, 3);
    }
    if (nunmatched > 2.0 * UNMATCHEDFOR2HOP * graph.nvtxs) {
        cnvtxs = Match_2HopAny(ctrl, graph, perm, match__, cnvtxs, &nunmatched, graph.nvtxs);
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
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    perm: *mut idx_t,
    match__: *mut idx_t,
    cnvtxs: idx_t,
    r_nunmatched: *mut size_t,
    maxdegree: size_t,
) -> idx_t {
    // idx_t i, pi, ii, j, jj, k, nvtxs;
    // idx_t *xadj, *adjncy, *colptr, *rowind;
    // idx_t *cmap;
    // size_t nunmatched;

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.Aux3Tmr));

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    cmap = graph.cmap;

    nunmatched = *r_nunmatched;

    /*IFSET(ctrl.dbglvl, METIS_DBG_COARSEN, printf("IN: nunmatched: %zu\t", nunmatched)); */

    /* create the inverted index */
    WCOREPUSH;
    colptr = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs + 1));
    for i in (0)..(nvtxs) {
        if (match__[i] == UNMATCHED && xadj[i + 1] - xadj[i] < maxdegree) {
            for j in (xadj[i])..(xadj[i + 1]) {
                colptr[adjncy[j]] += 1;
            }
        }
    }
    MAKECSR(i, nvtxs, colptr);

    rowind = iwspacemalloc(ctrl, colptr[nvtxs]);
    for pi in (0)..(nvtxs) {
        i = perm[pi];
        if (match__[i] == UNMATCHED && xadj[i + 1] - xadj[i] < maxdegree) {
            for j in (xadj[i])..(xadj[i + 1]) {
                rowind[colptr[adjncy[j]]] = i;
                colptr[adjncy[j]] += 1
            }
        }
    }
    SHIFTCSR(i, nvtxs, colptr);

    /* compute matchings by going down the inverted index */
    for pi in (0)..(nvtxs) {
        i = perm[pi];
        if (colptr[i + 1] - colptr[i] < 2) {
            continue;
        }

        jj = colptr[i + 1];
        // for j in (colptr[i])..(jj) {
        j = colptr[i];
        while j < jj {
            if (match__[rowind[j]] == UNMATCHED) {
                jj -= 1;
                // for (jj--; jj>j; jj--) {
                while jj > j {
                    if (match__[rowind[jj]] == UNMATCHED) {
                        cmap[rowind[j]] = cmap[rowind[jj]] = cnvtxs;
                        cnvtxs += 1;
                        match__[rowind[j]] = rowind[jj];
                        match__[rowind[jj]] = rowind[j];
                        nunmatched -= 2;
                        break;
                    }
                    jj -= 1;
                }
                j += 1;
            }
        }
    }
    WCOREPOP;

    /*IFSET(ctrl.dbglvl, METIS_DBG_COARSEN, printf("OUT: nunmatched: %zu\n", nunmatched)); */

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.Aux3Tmr));

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
    match__: *mut idx_t,
    cnvtxs: idx_t,
    r_nunmatched: *mut size_t,
    maxdegree: size_t,
) -> idx_t {
    // idx_t i, pi, pk, ii, j, jj, k, nvtxs, mask, idegree;
    // idx_t *xadj, *adjncy;
    // idx_t *cmap, *mark;
    // ikv_t *keys;
    // size_t nunmatched, ncand;

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.Aux3Tmr));

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    cmap = graph.cmap;

    nunmatched = *r_nunmatched;
    mask = IDX_MAX / maxdegree;

    /*IFSET(ctrl.dbglvl, METIS_DBG_COARSEN, printf("IN: nunmatched: %zu\t", nunmatched)); */

    WCOREPUSH;

    /* collapse vertices with identical adjancency lists */
    keys = ikvwspacemalloc(ctrl, nunmatched);
    ncand = 0;
    for pi in (0)..(nvtxs) {
        i = perm[pi];
        idegree = xadj[i + 1] - xadj[i];
        if (match__[i] == UNMATCHED && idegree > 1 && idegree < maxdegree) {
            k = 0;
            for j in (xadj[i])..(xadj[i + 1]) {
                k += adjncy[j] % mask;
            }
            keys[ncand].val = i;
            keys[ncand].key = (k % mask) * maxdegree + idegree;
            ncand;
            ncand += 1;
        }
    }
    ikvsorti(ncand, keys);

    mark = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
    for pi in (0)..(ncand) {
        i = keys[pi].val;
        if (match__[i] != UNMATCHED) {
            continue;
        }

        for j in (xadj[i])..(xadj[i + 1]) {
            mark[adjncy[j]] = i;
        }

        for pk in (pi + 1)..(ncand) {
            k = keys[pk].val;
            if (match__[k] != UNMATCHED) {
                continue;
            }

            if (keys[pi].key != keys[pk].key) {
                break;
            }
            if (xadj[i + 1] - xadj[i] != xadj[k + 1] - xadj[k]) {
                break;
            }

            for jj in (xadj[k])..(xadj[k + 1]) {
                if (mark[adjncy[jj]] != i) {
                    break;
                }
            }
            if (jj == xadj[k + 1]) {
                cmap[i] = cmap[k] = cnvtxs;
                cnvtxs += 1;
                match__[i] = k;
                match__[k] = i;
                nunmatched -= 2;
                break;
            }
        }
    }
    WCOREPOP;

    /*IFSET(ctrl.dbglvl, METIS_DBG_COARSEN, printf("OUT: ncand: %zu, nunmatched: %zu\n", ncand, nunmatched)); */

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.Aux3Tmr));

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
    // idx_t i, pi, ii, iii, j, jj, jjj, jjinc, k, nvtxs, ncon, cnvtxs, maxidx,
    // last_unmatched, avgdegree, bnum;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *maxvwgt;
    // idx_t *match__, *cmap, *degrees, *perm, *tperm, *vec, *marker;
    // idx_t mytwgt, xtwgt, ctwgt;
    // real_t bscore, score;

    WCOREPUSH;

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.MatchTmr));

    nvtxs = graph.nvtxs;
    ncon = graph.ncon;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;
    cmap = graph.cmap;

    maxvwgt = ctrl.maxvwgt;

    match__ = iset(nvtxs, UNMATCHED, iwspacemalloc(ctrl, nvtxs));
    perm = iwspacemalloc(ctrl, nvtxs);
    tperm = iwspacemalloc(ctrl, nvtxs);
    degrees = iwspacemalloc(ctrl, nvtxs);

    irandArrayPermute(nvtxs, tperm, nvtxs / 8, 1);

    avgdegree = 4.0 * (xadj[nvtxs] / nvtxs);
    for i in (0)..(nvtxs) {
        bnum = sqrt(1 + xadj[i + 1] - xadj[i]);
        degrees[i] = (if bnum > avgdegree { avgdegree } else { bnum });
    }
    BucketSortKeysInc(ctrl, nvtxs, avgdegree, degrees, tperm, perm);

    /* point to the wspace vectors that are not needed any more */
    vec = tperm;
    marker = degrees;
    iset(nvtxs, -1, vec);
    iset(nvtxs, -1, marker);

    cnvtxs = 0;
    last_unmatched = 0;
    for pi in (0)..(nvtxs) {
        i = perm[pi];

        if (match__[i] == UNMATCHED) {
            /* Unmatched */
            maxidx = i;

            if (if ncon == 1 {
                vwgt[i] < maxvwgt[0]
            } else {
                ivecle(ncon, vwgt + i * ncon, maxvwgt)
            }) {
                /* Deal with island vertices. Find a non-island and match__ it with.
                The matching ignores ctrl.maxvwgt requirements */
                if (xadj[i] == xadj[i + 1]) {
                    // last_unmatched = gk_max(pi, last_unmatched)+1;
                    for last_unmatched in (gk_max(pi, last_unmatched) + 1)..(nvtxs) {
                        j = perm[last_unmatched];
                        if (match__[j] == UNMATCHED) {
                            maxidx = j;
                            break;
                        }
                    }
                } else {
                    if (ncon == 1) {
                        /* Find a max JC pair, subject to maxvwgt constraints */
                        if (xadj[i + 1] - xadj[i] < avgdegree) {
                            marker[i] = i;
                            bscore = 0.0;
                            mytwgt = 0;
                            for j in (xadj[i])..(xadj[i + 1]) {
                                mytwgt += 1; //adjwgt[j];
                                vec[adjncy[j]] = 1; //adjwgt[j];
                            }

                            /* single constraint pairing */
                            // #ifdef XXX
                            //               for j in (xadj[i])..(xadj[i+1]) {
                            //                 ii = adjncy[j];
                            //                 if (marker[ii] == i || match__[ii] != UNMATCHED || vwgt[i]+vwgt[ii] > maxvwgt[0])
                            //         {
                            //             continue;
                            //         }
                            //
                            //                 ctwgt = xtwgt = 0;
                            //                 for jj in (xadj[ii])..(xadj[ii+1]) {
                            //                   xtwgt += adjwgt[jj];
                            //                   if (vec[adjncy[jj]] > 0)
                            //         {
                            //             ctwgt += vec[adjncy[jj]] + adjwgt[jj];
                            //         }
                            //                   else if (adjncy[jj] == i) {
                            //                     ctwgt += adjwgt[jj];
                            //                     xtwgt -= adjwgt[jj];
                            //                   }
                            //                 }
                            //
                            //                 score = 1.0*ctwgt/(mytwgt+xtwgt-ctwgt);
                            //                 if (score > bscore) {
                            //                   bscore = score;
                            //                   maxidx = ii;
                            //                 }
                            //                 marker[ii] = i;
                            //               }
                            // #endif

                            for j in (xadj[i])..(xadj[i + 1]) {
                                ii = adjncy[j];
                                for jj in (xadj[ii])..(xadj[ii + 1]) {
                                    iii = adjncy[jj];

                                    if (marker[iii] == i
                                        || match__[iii] != UNMATCHED
                                        || vwgt[i] + vwgt[iii] > maxvwgt[0])
                                    {
                                        continue;
                                    }

                                    ctwgt = xtwgt = 0;
                                    for jjj in (xadj[iii])..(xadj[iii + 1]) {
                                        xtwgt += 1; //adjwgt[jjj];
                                        if (vec[adjncy[jjj]] > 0) {
                                            ctwgt += 2; //vec[adjncy[jjj]] + adjwgt[jjj];
                                        } else if (adjncy[jjj] == i) {
                                            ctwgt += 10 * adjwgt[jjj];
                                        }
                                    }

                                    score = 1.0 * ctwgt / (mytwgt + xtwgt);
                                    //printf("{:} {:} {:} %.4f\n", mytwgt, xtwgt, ctwgt, score);
                                    if (score > bscore) {
                                        bscore = score;
                                        maxidx = iii;
                                    }
                                    marker[iii] = i;
                                }
                            }

                            /* reset vec array */
                            for j in (xadj[i])..(xadj[i + 1]) {
                                vec[adjncy[j]] = -1;
                            }
                        }
                    } else {
                        /* multi-constraint version */
                        for j in (xadj[i])..(xadj[i + 1]) {
                            k = adjncy[j];
                            if (match__[k] == UNMATCHED
                                && ivecaxpylez(ncon, 1, vwgt + i * ncon, vwgt + k * ncon, maxvwgt))
                            {
                                maxidx = k;
                                break;
                            }
                        }
                    }
                }
            }

            if (maxidx != UNMATCHED) {
                cmap[i] = cmap[maxidx] = cnvtxs;
                cnvtxs += 1;
                match__[i] = maxidx;
                match__[maxidx] = i;
            }
        }
    }

    /* match__ the final unmatched vertices with themselves and reorder the vertices
    of the coarse graph for memory-friendly contraction */
    cnvtxs = 0;
    for i in (0)..(nvtxs) {
        if (match__[i] == UNMATCHED) {
            match__[i] = i;
            cmap[i] = cnvtxs;
            cnvtxs += 1;
        } else {
            if (i <= match__[i]) {
                cmap[i] = cmap[match__[i]] = cnvtxs;
                cnvtxs += 1;
            }
        }
    }

    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.MatchTmr));

    CreateCoarseGraph(ctrl, graph, cnvtxs, match__);

    WCOREPOP;

    return cnvtxs;
}

/*************************************************************************/
/* This function prints various stats for each graph during coarsening
 */
/*************************************************************************/
#[metis_func]
pub extern "C" fn PrintCGraphStats(ctrl: *mut ctrl_t, graph: *mut graph_t) -> void {
    // idx_t i;

    printf(
        "{:10} {:10} {:10} [{:}] [",
        graph.nvtxs,
        graph.nedges,
        isum(graph.nedges, graph.adjwgt, 1),
        ctrl.CoarsenTo,
    );

    for i in (0)..(graph.ncon) {
        printf(" {:8}:%8"PRIDX, ctrl.maxvwgt[i], graph.tvwgt[i]);
    }
    printf(" ]\n");
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
    match__: *mut idx_t,
) -> void {
    // idx_t j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon,
    // cnedges, v, u, mask;
    // idx_t *xadj, *vwgt, *vsize, *adjncy, *adjwgt;
    // idx_t *cmap, *htable, *dtable;
    // idx_t *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt;
    // graph_t *cgraph;
    // int dovsize, dropedges;
    // idx_t cv, nkeys, droppedewgt;
    // idx_t *keys=NULL, *medianewgts=NULL, *noise=NULL;

    WCOREPUSH;

    dovsize = (if ctrl.objtype == METIS_OBJTYPE_VOL {
        1
    } else {
        0
    });
    dropedges = ctrl.dropedges;

    mask = HTLENGTH;

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_TIME,
        gk_startcputimer(ctrl.ContractTmr),
    );

    nvtxs = graph.nvtxs;
    ncon = graph.ncon;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    vsize = graph.vsize;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;
    cmap = graph.cmap;

    /* Setup structures for dropedges */
    if (dropedges) {
        nkeys = 0;
        for v in (0)..(nvtxs) {
            nkeys = gk_max(nkeys, xadj[v + 1] - xadj[v]);
        }
        nkeys = 2 * nkeys + 1;

        keys = iwspacemalloc(ctrl, nkeys);
        noise = iwspacemalloc(ctrl, cnvtxs);
        medianewgts = iset(cnvtxs, -1, iwspacemalloc(ctrl, cnvtxs));

        for v in (0)..(cnvtxs) {
            noise[v] = irandInRange(128);
        }
    }

    /* Initialize the coarser graph */
    cgraph = SetupCoarseGraph(graph, cnvtxs, dovsize);
    cxadj = cgraph.xadj;
    cvwgt = cgraph.vwgt;
    cvsize = cgraph.vsize;
    cadjncy = cgraph.adjncy;
    cadjwgt = cgraph.adjwgt;

    htable = iset(mask + 1, -1, iwspacemalloc(ctrl, mask + 1)); /* hash table */
    dtable = iset(cnvtxs, -1, iwspacemalloc(ctrl, cnvtxs)); /* direct table */

    cxadj[0] = cnvtxs = cnedges = 0;
    for v in (0)..(nvtxs) {
        if ((u = match__[v]) < v) {
            continue;
        }

        ASSERT(cmap[v] == cnvtxs);
        ASSERT(cmap[match__[v]] == cnvtxs);

        /* take care of the vertices */
        if (ncon == 1) {
            cvwgt[cnvtxs] = vwgt[v];
        } else {
            icopy(ncon, vwgt + v * ncon, cvwgt + cnvtxs * ncon);
        }

        if (dovsize) {
            cvsize[cnvtxs] = vsize[v];
        }

        if (v != u) {
            if (ncon == 1) {
                cvwgt[cnvtxs] += vwgt[u];
            } else {
                iaxpy(ncon, 1, vwgt + u * ncon, 1, cvwgt + cnvtxs * ncon, 1);
            }

            if (dovsize) {
                cvsize[cnvtxs] += vsize[u];
            }
        }

        /* take care of the edges */
        if ((xadj[v + 1] - xadj[v] + xadj[u + 1] - xadj[u]) < (mask >> 2)) {
            /* use mask */
            /* put the ID of the contracted node itself at the start, so that it can be
             * removed easily */
            htable[cnvtxs & mask] = 0;
            cadjncy[0] = cnvtxs;
            nedges = 1;

            istart = xadj[v];
            iend = xadj[v + 1];
            for j in (istart)..(iend) {
                k = cmap[adjncy[j]];
                // for (kk=k&mask; htable[kk]!=-1 && cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
                kk = k & mask;
                while htable[kk] != -1 && cadjncy[htable[kk]] != k {
                    kk = ((kk + 1) & mask)
                }

                if ((m = htable[kk]) == -1) {
                    cadjncy[nedges] = k;
                    cadjwgt[nedges] = adjwgt[j];
                    htable[kk] = nedges;
                    nedges += 1;
                } else {
                    cadjwgt[m] += adjwgt[j];
                }
            }

            if (v != u) {
                istart = xadj[u];
                iend = xadj[u + 1];
                for j in (istart)..(iend) {
                    k = cmap[adjncy[j]];
                    // for (kk=k&mask; htable[kk]!=-1 && cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
                    kk = k & mask;
                    while htable[kk] != -1 && cadjncy[htable[kk]] != k {
                        kk = ((kk + 1) & mask)
                    }

                    if ((m = htable[kk]) == -1) {
                        cadjncy[nedges] = k;
                        cadjwgt[nedges] = adjwgt[j];
                        htable[kk] = nedges;
                        nedges += 1;
                    } else {
                        cadjwgt[m] += adjwgt[j];
                    }
                }
            }

            /* reset the htable -- reverse order (LIFO) is critical to prevent cadjncy[-1]
             * indexing due to a remove of an earlier entry */
            // for (j=nedges-1; j>=0; j--) {
            for j in 0..=(nedges - 1).rev() {
                k = cadjncy[j];
                // for (kk=k&mask; cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
                kk = k & mask;
                while cadjncy[htable[kk]] != k {
                    kk = ((kk + 1) & mask);
                }
                htable[kk] = -1;
                // TODO: Verify this
            }

            /* remove the contracted vertex from the list */
            nedges -= 1;
            cadjncy[0] = cadjncy[nedges];
            cadjwgt[0] = cadjwgt[nedges];
        } else {
            nedges = 0;
            istart = xadj[v];
            iend = xadj[v + 1];
            for j in (istart)..(iend) {
                k = cmap[adjncy[j]];
                if ((m = dtable[k]) == -1) {
                    cadjncy[nedges] = k;
                    cadjwgt[nedges] = adjwgt[j];
                    dtable[k] = nedges;
                    nedges += 1;
                } else {
                    cadjwgt[m] += adjwgt[j];
                }
            }

            if (v != u) {
                istart = xadj[u];
                iend = xadj[u + 1];
                for j in (istart)..(iend) {
                    k = cmap[adjncy[j]];
                    if ((m = dtable[k]) == -1) {
                        cadjncy[nedges] = k;
                        cadjwgt[nedges] = adjwgt[j];
                        dtable[k] = nedges;
                        nedges += 1;
                    } else {
                        cadjwgt[m] += adjwgt[j];
                    }
                }

                /* Remove the contracted self-loop, when present */
                if ((j = dtable[cnvtxs]) != -1) {
                    ASSERT(cadjncy[j] == cnvtxs);
                    nedges -= 1;
                    cadjncy[j] = cadjncy[nedges];
                    cadjwgt[j] = cadjwgt[nedges];
                    dtable[cnvtxs] = -1;
                }
            }

            /* Zero out the dtable */
            for j in (0)..(nedges) {
                dtable[cadjncy[j]] = -1;
            }
        }

        /* Determine the median weight of the incident edges, which will be used
        to keep an edge (u, v) iff wgt(u, v) >= min(medianewgts[u], medianewgts[v]) */
        if (dropedges) {
            ASSERTP(nedges < nkeys, ("{:}, {:}\n", nkeys, nedges));
            medianewgts[cnvtxs] = 8; /* default for island nodes */
            if (nedges > 0) {
                for j in (0)..(nedges) {
                    keys[j] = (cadjwgt[j] << 8) + noise[cnvtxs] + noise[cadjncy[j]];
                }
                isortd(nedges, keys);
                medianewgts[cnvtxs] = keys[gk_min(
                    nedges - 1,
                    ((xadj[v + 1] - xadj[v] + xadj[u + 1] - xadj[u]) >> 1),
                )];
            }
        }

        cadjncy += nedges;
        cadjwgt += nedges;
        cnedges += nedges;
        cnvtxs += 1;
        cxadj[cnvtxs] = cnedges;
    }

    /* compact the adjacency structure of the coarser graph to keep only +ve edges */
    if (dropedges) {
        droppedewgt = 0;

        cadjncy = cgraph.adjncy;
        cadjwgt = cgraph.adjwgt;

        cnedges = 0;
        for u in (0)..(cnvtxs) {
            istart = cxadj[u];
            iend = cxadj[u + 1];
            for j in (istart)..(iend) {
                v = cadjncy[j];
                ASSERTP(medianewgts[u] >= 0, ("{:} {:}\n", u, medianewgts[u]));
                ASSERTP(
                    medianewgts[v] >= 0,
                    ("{:} {:} {:}\n", v, medianewgts[v], cnvtxs),
                );
                if ((cadjwgt[j] << 8) + noise[u] + noise[v]
                    >= gk_min(medianewgts[u], medianewgts[v]))
                {
                    cadjncy[cnedges] = cadjncy[j];
                    cadjwgt[cnedges] = cadjwgt[j];
                    cnedges += 1;
                } else {
                    droppedewgt += cadjwgt[j];
                }
            }
            cxadj[u] = cnedges;
        }
        SHIFTCSR(j, cnvtxs, cxadj);

        cgraph.droppedewgt = droppedewgt;
    }

    cgraph.nedges = cnedges;

    for j in (0)..(ncon) {
        cgraph.tvwgt[j] = isum(cgraph.nvtxs, cgraph.vwgt + j, ncon);
        cgraph.invtvwgt[j] = 1.0
            / (if cgraph.tvwgt[j] > 0 {
                cgraph.tvwgt[j]
            } else {
                1
            });
    }

    ReAdjustMemory(ctrl, graph, cgraph);

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_TIME,
        gk_stopcputimer(ctrl.ContractTmr),
    );

    WCOREPOP;
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
    // graph_t *cgraph;

    cgraph = CreateGraph();

    cgraph.nvtxs = cnvtxs;
    cgraph.ncon = graph.ncon;

    cgraph.finer = graph;
    graph.coarser = cgraph;

    /* Allocate memory for the coarser graph.
    NOTE: The +1 in the adjwgt/adjncy is to allow the optimization of self-loop
          detection by adding ahead of time the self-loop. That optimization
          requires a +1 adjncy/adjwgt array for the limit case where_ the
          coarser graph is of the same size of the previous graph. */
    cgraph.xadj = imalloc(cnvtxs + 1, "SetupCoarseGraph: xadj");
    cgraph.adjncy = imalloc(graph.nedges + 1, "SetupCoarseGraph: adjncy");
    cgraph.adjwgt = imalloc(graph.nedges + 1, "SetupCoarseGraph: adjwgt");
    cgraph.vwgt = imalloc(cgraph.ncon * cnvtxs, "SetupCoarseGraph: vwgt");
    cgraph.tvwgt = imalloc(cgraph.ncon, "SetupCoarseGraph: tvwgt");
    cgraph.invtvwgt = rmalloc(cgraph.ncon, "SetupCoarseGraph: invtvwgt");

    if (dovsize) {
        cgraph.vsize = imalloc(cnvtxs, "SetupCoarseGraph: vsize");
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
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    cgraph: *mut graph_t,
) -> void {
    if (cgraph.nedges > 10000 && cgraph.nedges < 0.9 * graph.nedges) {
        cgraph.adjncy = irealloc(cgraph.adjncy, cgraph.nedges, "ReAdjustMemory: adjncy");
        cgraph.adjwgt = irealloc(cgraph.adjwgt, cgraph.nedges, "ReAdjustMemory: adjwgt");
    }
}
