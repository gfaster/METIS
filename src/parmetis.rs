/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * parmetis.c
 *
 * This file contains top level routines that are used by ParMETIS
 *
 * Started 10/14/97
 * George
 *
 * $Id: parmetis.c 10481 2011-07-05 18:01:23Z karypis $
 *
 */

use std::ffi::c_int;

use crate::*;

/*************************************************************************/
/* This function is the entry point for the node ND code for ParMETIS.
    The difference between this routine and the standard METIS_NodeND are
    the following

    - It performs at least log2(npes) levels of nested dissection.
    - It stores the size of the log2(npes) top-level separators in the
      sizes array.
*/
/*************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_NodeNDP(
    nvtxs: idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    npes: idx_t,
    options: *const idx_t,
    perm: *mut idx_t,
    iperm: *mut idx_t,
    sizes: *mut idx_t,
) -> c_int {
    let nvtxs = nvtxs as usize;
    mkslice_mut!(perm, nvtxs);
    mkslice_mut!(iperm, nvtxs);
    mkslice_mut!(sizes, 2 * npes - 1);
    // idx_t i, ii, j, l, nnvtxs=0;
    // graph_t *graph;
    // ctrl_t *ctrl;
    // idx_t *cptr, *cind;
    let ctrl = options::SetupCtrl(
        METIS_OP_OMETIS,
        options,
        1,
        3,
        std::ptr::null_mut(),
        std::ptr::null_mut(),
    );
    let Some(ctrl) = ctrl.as_mut() else {
        return METIS_ERROR_INPUT;
    };

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.TotalTmr));
    let mut graph = std::ptr::null_mut();
    let mut nnvtxs = 0;

    /* compress the graph; note that compression only happens if no prunning
    has taken place. */
    let mut cptr = vec![];
    let mut cind = vec![];
    if ctrl.compress != 0 {
        cptr = vec![0; nvtxs + 1];
        cind = vec![0; nvtxs];

        graph = compress::CompressGraph(
            ctrl,
            nvtxs as idx_t,
            xadj,
            adjncy,
            vwgt,
            cptr.as_mut_ptr(),
            cind.as_mut_ptr(),
        );
        if graph.is_null() {
            /* if there was no compression, cleanup the compress flag */
            let _ = std::mem::take(&mut cptr);
            let _ = std::mem::take(&mut cind);
            ctrl.compress = 0;
        } else {
            nnvtxs = (*graph).nvtxs as usize;
        }
    }

    /* if no compression, setup the graph in the normal way. */
    if ctrl.compress == 0 {
        graph = graph::SetupGraph(
            ctrl,
            nvtxs as idx_t,
            1,
            xadj,
            adjncy,
            vwgt,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );
    }
    let graph = graph.as_mut().unwrap();

    /* allocate workspace memory */
    wspace::AllocateWorkSpace(ctrl, graph);

    /* do the nested dissection ordering  */
    // iset(2 * npes - 1, 0, sizes);
    sizes[..2 * npes as usize - 1].fill(0);
    MlevelNestedDissectionP(
        ctrl,
        graph,
        iperm.as_mut_ptr(),
        graph.nvtxs,
        npes,
        0,
        sizes.as_mut_ptr(),
    );

    /* Uncompress the ordering */
    if ctrl.compress != 0 {
        /* construct perm from iperm */
        for i in (0)..(nnvtxs) {
            perm[iperm[i as usize] as usize] = i as idx_t;
        }
        let mut l = 0;
        for ii in 0..nnvtxs {
            let i = perm[ii as usize] as usize;
            for j in (cptr[i])..(cptr[i + 1]) {
                iperm[cind[j as usize] as usize] = l;
                l += 1;
            }
        }

        let _ = std::mem::take(&mut cptr);
        let _ = std::mem::take(&mut cind);
        // gk_free((void **)&cptr, &cind, LTERM);
    }

    for i in (0)..(nvtxs) {
        perm[iperm[i as usize] as usize] = i as idx_t;
    }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

    /* clean up */
    let mut ctrl: *mut _ = ctrl;
    options::FreeCtrl(&mut ctrl);

    return METIS_OK;
}

/*************************************************************************/
/* This function is similar to MlevelNestedDissection with the difference
that it also records separator sizes for the top log2(npes) levels */
/**************************************************************************/
// NOTE(porting): I think this may be removed as this port is not meant to replace ParMETIS (as it
// is not Free Software). That being said, NodeNDP seems like it might be worth keeping?
#[metis_func]
pub extern "C" fn MlevelNestedDissectionP(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    order: *mut idx_t,
    lastvtx: idx_t,
    npes: idx_t,
    cpos: idx_t,
    sizes: *mut idx_t,
) {
    let mut lastvtx = lastvtx;
    // idx_t i, j, nvtxs, nbnd;
    // idx_t *label, *bndind;
    // graph_t *lgraph, *rgraph;
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let nvtxs = graph.nvtxs as usize;

    if nvtxs == 0 {
        let mut graph: *mut _ = graph;
        graph::FreeGraph(&mut graph);
        return;
    }

    ometis::MlevelNodeBisectionMultiple(ctrl, graph);
    get_graph_slices!(ctrl, graph => pwgts);
    mkslice_mut!(sizes, 2 * npes - 1);

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_SEPINFO,
        print!(
            "Nvtxs: {:6}, [{:6} {:6} {:6}]\n",
            graph.nvtxs, pwgts[0], pwgts[1], pwgts[2]
        )
    );

    if cpos < npes - 1 {
        sizes[(2 * npes - 2 - cpos) as usize] = pwgts[2];
        sizes[(2 * npes - 2 - (2 * cpos + 1)) as usize] = pwgts[1];
        sizes[(2 * npes - 2 - (2 * cpos + 2)) as usize] = pwgts[0];
    }

    /* Order the nodes in the separator */
    let nbnd = graph.nbnd as usize;

    get_graph_slices!(graph => bndind label);
    for i in (0)..(nbnd) {
        lastvtx -= 1;
        let ind = bndind[i as usize] as usize;
        let label = label[ind] as usize;
        *order.add(label) = lastvtx;
    }

    let mut lgraph = std::ptr::null_mut();
    let mut rgraph = std::ptr::null_mut();
    ometis::SplitGraphOrder(ctrl, graph, &mut lgraph, &mut rgraph);

    /* Free the memory of the top level graph */
    let mut graph: *mut _ = graph;
    graph::FreeGraph(&mut graph);

    if ((*lgraph).nvtxs > MMDSWITCH || 2 * cpos + 2 < npes - 1) && (*lgraph).nedges > 0 {
        MlevelNestedDissectionP(
            ctrl,
            lgraph,
            order,
            lastvtx - (*rgraph).nvtxs,
            npes,
            2 * cpos + 2,
            sizes.as_mut_ptr(),
        );
    } else {
        ometis::MMDOrder(ctrl, lgraph, order, lastvtx - (*rgraph).nvtxs);
        graph::FreeGraph(&mut lgraph);
    }
    if ((*rgraph).nvtxs > MMDSWITCH || 2 * cpos + 1 < npes - 1) && (*rgraph).nedges > 0 {
        MlevelNestedDissectionP(
            ctrl,
            rgraph,
            order,
            lastvtx,
            npes,
            2 * cpos + 1,
            sizes.as_mut_ptr(),
        );
    } else {
        ometis::MMDOrder(ctrl, rgraph, order, lastvtx);
        graph::FreeGraph(&mut rgraph);
    }
}

/*************************************************************************/
/* This function bisects a graph by computing a vertex separator
*/
/**************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_ComputeVertexSeparator(
    nvtxs: *mut idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    options: *mut idx_t,
    r_sepsize: *mut idx_t,
    part: *mut idx_t,
) -> c_int {
    // idx_t i, j;
    // graph_t *graph;
    // ctrl_t *ctrl;
    let ctrl = options::SetupCtrl(
        METIS_OP_OMETIS,
        options,
        1,
        3,
        std::ptr::null_mut(),
        std::ptr::null_mut(),
    );
    let Some(ctrl) = ctrl.as_mut() else {
        return METIS_ERROR_INPUT;
    };

    util::InitRandom(ctrl.seed);

    let mut graph = graph::SetupGraph(
        ctrl,
        *nvtxs,
        1,
        xadj,
        adjncy,
        vwgt,
        std::ptr::null_mut(),
        std::ptr::null_mut(),
    );

    wspace::AllocateWorkSpace(ctrl, graph);

    /*============================================================
     * Perform the bisection
     *============================================================*/
    ctrl.CoarsenTo = 100;

    ometis::MlevelNodeBisectionMultiple(ctrl, graph);

    *r_sepsize = *(*graph).pwgts.add(2);
    std::ptr::copy_nonoverlapping((*graph).where_, part, *nvtxs as usize);
    // icopy(*nvtxs, graph.where_, part);

    graph::FreeGraph(&mut graph);

    let mut ctrl: *mut _ = ctrl;
    options::FreeCtrl(&mut ctrl);

    return METIS_OK;
}

/*************************************************************************/
/* This function is the entry point of a node-based separator refinement
of the nodes with an hmarker[] of 0. */
/*************************************************************************/
// NOTE(porting): I think this will be removed as this port is not meant to replace ParMETIS (as it
// is not Free Software). It is also not called by any other routine.
#[metis_func(no_pfx)]
pub extern "C" fn METIS_NodeRefine(
    nvtxs: idx_t,
    xadj: *mut idx_t,
    vwgt: *mut idx_t,
    adjncy: *mut idx_t,
    where_: *mut idx_t,
    hmarker: *mut idx_t,
    ubfactor: real_t,
) -> c_int {
    // graph_t *graph;
    // ctrl_t *ctrl;

    /* set up the run time parameters */
    let ctrl = options::SetupCtrl(
        METIS_OP_OMETIS,
        std::ptr::null_mut(),
        1,
        3,
        std::ptr::null_mut(),
        std::ptr::null_mut(),
    );
    let Some(ctrl) = ctrl.as_mut() else {
        return METIS_ERROR_INPUT;
    };

    /* set up the graph */
    let mut graph = graph::SetupGraph(
        ctrl,
        nvtxs,
        1,
        xadj,
        adjncy,
        vwgt,
        std::ptr::null_mut(),
        std::ptr::null_mut(),
    );

    /* allocate workspace memory */
    wspace::AllocateWorkSpace(ctrl, graph);

    /* set up the memory and the input partition */
    srefine::Allocate2WayNodePartitionMemory(ctrl, graph);
    std::ptr::copy_nonoverlapping(where_, (*graph).where_, nvtxs as usize);

    srefine::Compute2WayNodePartitionParams(ctrl, graph);

    FM_2WayNodeRefine1SidedP(ctrl, graph, hmarker, ubfactor, 10);
    /* FM_2WayNodeRefine2SidedP(ctrl, graph, hmarker, ubfactor, 10); */

    std::ptr::copy_nonoverlapping((*graph).where_, where_, nvtxs as usize);
    // icopy(nvtxs, graph.where_, where_);

    graph::FreeGraph(&mut graph);
    let mut ctrl: *mut _ = ctrl;
    options::FreeCtrl(&mut ctrl);

    return METIS_OK;
}

/*************************************************************************/
/* This function performs a node-based 1-sided FM refinement that moves
only nodes whose hmarker[] == -1. It is used by Parmetis. */
/*************************************************************************/
// NOTE(porting): I think this will be removed as this port is not meant to replace ParMETIS (as it
// is not Free Software). It is only called by METIS_NodeRefine
#[metis_func]
pub extern "C" fn FM_2WayNodeRefine1SidedP(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    hmarker: *mut idx_t,
    ubfactor: real_t,
    npasses: idx_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind, nbad, qsize;
    // idx_t *xadj, *vwgt, *adjncy, *where_, *pwgts, *edegrees, *bndind, *bndptr;
    // idx_t *mptr, *mind, *swaps, *inqueue;
    // rpq_t *queue;
    // nrinfo_t *rinfo;
    // idx_t higain, oldgain, mincut, initcut, mincutorder;
    // idx_t pass, from, to, limit;
    // idx_t badmaxpwgt, mindiff, newdiff;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy vwgt);
    get_graph_slices_mut!(ctrl, graph => where_ pwgts nrinfo bndind bndptr);
    debug_assert_eq!(graph.mincut, pwgts[2]);
    let rinfo = nrinfo;
    mkslice_mut!(hmarker, nvtxs);

    let mut queue = pqueue::RPQueue::new(nvtxs);

    let mut inqueue = vec![-1; nvtxs];
    let mut swaps = vec![0; nvtxs];
    let mut mptr = vec![0; nvtxs + 1];
    let mut mind = vec![0; 2 * nvtxs];

    let badmaxpwgt = (ubfactor * ((pwgts[0]).max(pwgts[1])) as real_t) as idx_t;

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        print!(
            "Partitions-N1: [({:6} {:6}) as usize] Nv-Nb[({:6} {:6}) as usize] \
           MaxPwgt[{:6}]. ISep: {:6}\n",
            pwgts[0], pwgts[1], graph.nvtxs, graph.nbnd, badmaxpwgt, graph.mincut
        )
    );

    let mut to = if pwgts[0] < pwgts[1] { 1 } else { 0 };
    for pass in (0)..(npasses) {
        let from = to;
        to = (from + 1) % 2;

        queue.reset();
        // rpqReset(queue);

        let mut mincutorder: idx_t = -1;
        let initcut = graph.mincut;
        let mut mincut = graph.mincut;
        let mut nbnd = graph.nbnd as usize;

        /* use the swaps array in place of the traditional perm array to save memory */
        irandArrayPermute(nbnd as idx_t, swaps.as_mut_ptr(), nbnd as idx_t, 1);
        for ii in (0)..(nbnd) {
            let i = bndind[swaps[ii as usize] as usize];
            assert!(where_[i as usize] == 2);
            if hmarker[i as usize] == -1 || hmarker[i as usize] == to {
                // rpqInsert(queue, i, vwgt[i as usize]-rinfo[i as usize].edegrees[from as usize]);
                queue.insert(
                    i,
                    (vwgt[i as usize] - rinfo[i as usize].edegrees[from as usize]) as real_t,
                );
                inqueue[i as usize] = pass;
            }
        }
        let qsize = queue.length();

        debug_assert!(debug::CheckNodeBnd(graph, nbnd as idx_t) != 0);
        debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

        let limit = nbnd;

        /******************************************************
         * Get into the FM loop
         *******************************************************/
        mptr[0] = 0;
        let mut nmind = 0;
        let mut nbad = 0;
        let mut mindiff = (pwgts[0] - pwgts[1]).abs();
        let mut nswaps = -1;
        // for nswaps in (0)..(nvtxs)
        while {
            nswaps += 1;
            nswaps
        } < nvtxs as idx_t
        {
            let Some(higain) = queue.pop() else { break };

            assert!(bndptr[higain as usize] != -1);

            /* The following check is to ensure we break out if there is a possibility
            of over-running the mind array.  */
            if nmind + xadj[(higain + 1) as usize] - xadj[higain as usize] >= 2 * nvtxs as idx_t - 1
            {
                break;
            }

            inqueue[higain as usize] = -1;

            if pwgts[to as usize] + vwgt[higain as usize] > badmaxpwgt {
                /* Skip this vertex */
                if nbad > limit {
                    break;
                } else {
                    nbad += 1;
                    nswaps -= 1;
                    continue;
                }
            }

            pwgts[2] -= vwgt[higain as usize] - rinfo[higain as usize].edegrees[from as usize];

            let newdiff = (pwgts[to as usize] + vwgt[higain as usize]
                - (pwgts[from as usize] - rinfo[higain as usize].edegrees[from as usize]))
                .abs();
            if pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff) {
                mincut = pwgts[2];
                mincutorder = nswaps as idx_t;
                mindiff = newdiff;
                nbad = 0;
            } else {
                if nbad > limit {
                    nbad += 1;
                    pwgts[2] +=
                        vwgt[higain as usize] - rinfo[higain as usize].edegrees[from as usize];
                    break; /* No further improvement, break out */
                }
            }

            BNDDelete!(nbnd, bndind, bndptr, higain as usize);
            pwgts[to as usize] += vwgt[higain as usize];
            where_[higain as usize] = to;
            swaps[nswaps as usize] = higain;

            /**********************************************************
             * Update the degrees of the affected nodes
             ***********************************************************/
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;
                if where_[k] == 2 {
                    /* For the in-separator vertices modify their edegree[to as usize] */
                    rinfo[k].edegrees[to as usize] += vwgt[higain as usize];
                } else if where_[k] == from {
                    /* This vertex is pulled into the separator */
                    debug_assert!(bndptr[k] == -1, "{:} {:} {:}\n", k, bndptr[k], where_[k]);
                    BNDInsert!(nbnd, bndind, bndptr, k);

                    mind[(nmind) as usize] = k; /* Keep track for rollback */
                    nmind += 1;
                    where_[k] = 2;
                    pwgts[from as usize] -= vwgt[k];

                    rinfo[k].edegrees = [0, 0];
                    for jj in (xadj[k])..(xadj[(k + 1) as usize]) {
                        let kk = adjncy[jj as usize] as usize;
                        if where_[kk] != 2 {
                            rinfo[k].edegrees[where_[kk] as usize] += vwgt[kk];
                        } else {
                            let oldgain = vwgt[kk] - rinfo[kk].edegrees[from as usize];
                            rinfo[kk].edegrees[from as usize] -= vwgt[k];

                            /* Update the gain of this node if it was not skipped */
                            if inqueue[kk] == pass {
                                queue.update(kk as idx_t, (oldgain + vwgt[k]) as real_t);
                            }
                        }
                    }

                    /* Insert the new vertex into the priority queue. Safe due to one-sided moves */
                    if hmarker[k] == -1 || hmarker[k] == to {
                        queue.insert(
                            k as idx_t,
                            (vwgt[k] - rinfo[k].edegrees[from as usize]) as real_t,
                        );
                        inqueue[k] = pass;
                    }
                }
            }
            mptr[(nswaps + 1) as usize] = nmind;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "Moved {:6} to {:3}, Gain: {:5} [{:5}] \t[({:5} {:5} {:5}) as usize] [({:3} {:2}) as usize]\n",
                    higain,
                    to,
                    (vwgt[higain as usize] - rinfo[higain as usize].edegrees[from as usize]),
                    vwgt[higain as usize],
                    pwgts[0],
                    pwgts[1],
                    pwgts[2],
                    nswaps,
                    limit
                )
            );
        }

        // silence unused modification
        let _ = nbad;

        /****************************************************************
         * Roll back computation
         *****************************************************************/
        for nswaps in ((mincutorder + 1)..nswaps).rev() {
            // for (nswaps--; nswaps>mincutorder; nswaps--)
            let higain = swaps[nswaps as usize] as usize;

            debug_assert!(debug::CheckNodePartitionParams(graph) != 0);
            debug_assert!(where_[higain] == to);

            inc_dec!(pwgts[2], pwgts[to as usize], vwgt[higain]);
            where_[higain] = 2;
            BNDInsert!(nbnd, bndind, bndptr, higain);

            rinfo[higain].edegrees = [0, 0];
            for j in (xadj[higain])..(xadj[higain + 1]) {
                let k = adjncy[j as usize] as usize;
                if where_[k] == 2 {
                    rinfo[k].edegrees[to as usize] -= vwgt[higain];
                } else {
                    rinfo[higain].edegrees[where_[k] as usize] += vwgt[k];
                }
            }

            /* Push nodes out of the separator */
            for j in (mptr[nswaps as usize])..(mptr[(nswaps + 1) as usize]) {
                let k = mind[j as usize] as usize;
                assert!(where_[k] == 2);
                where_[k] = from;
                inc_dec!(pwgts[from as usize], pwgts[2], vwgt[k]);
                BNDDelete!(nbnd, bndind, bndptr, k);
                for jj in (xadj[k])..(xadj[(k + 1) as usize]) {
                    let kk = adjncy[jj as usize] as usize;
                    if where_[kk] == 2 {
                        rinfo[kk].edegrees[from as usize] += vwgt[k];
                    }
                }
            }
        }

        debug_assert!(mincut == pwgts[2]);

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_REFINE,
            print!(
                "\tMinimum sep: {:6} at {:5}, PWGTS: [({:6} {:6}) as usize], NBND: {:6}, QSIZE: {:6}\n",
                mincut, mincutorder, pwgts[0], pwgts[1], nbnd, qsize
            )
        );

        graph.mincut = mincut;
        graph.nbnd = nbnd as idx_t;

        if pass % 2 == 1 && (mincutorder == -1 || mincut >= initcut) {
            break;
        }
    }

    // rpqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************/
/* This function performs a node-based (two-sided) FM refinement that
moves only nodes whose hmarker[] == -1. It is used by Parmetis. */
/*************************************************************************/
// NOTE(porting): I think this will be removed as this port is not meant to replace ParMETIS (as it
// is not Free Software). It is only called by METIS_NodeRefine
#[metis_func]
pub extern "C" fn FM_2WayNodeRefine2SidedP(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    hmarker: *mut idx_t,
    ubfactor: real_t,
    npasses: idx_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind;
    // idx_t *xadj, *vwgt, *adjncy, *where_, *pwgts, *edegrees, *bndind, *bndptr;
    // idx_t *mptr, *mind, *moved, *swaps;
    // rpq_t *queues[2];
    // nrinfo_t *rinfo;
    // idx_t higain, oldgain, mincut, initcut, mincutorder;
    // idx_t pass, to, other, limit;
    // idx_t badmaxpwgt, mindiff, newdiff;
    // idx_t u[2], g[2];

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy vwgt);
    get_graph_slices_mut!(ctrl, graph => nrinfo pwgts where_ bndind bndptr);
    let rinfo = nrinfo;
    mkslice!(hmarker, nvtxs);

    let mut queues = [pqueue::RPQueue::new(nvtxs), pqueue::RPQueue::new(nvtxs)];

    let mut moved = vec![0; nvtxs as usize];
    let mut swaps = vec![0; nvtxs as usize];
    let mut mptr = vec![0; nvtxs + 1 as usize];
    let mut mind = vec![0; 2 * nvtxs as usize];

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        print!(
            "Partitions: [({:6} {:6}) as usize] Nv-Nb[({:6} {:6}) as usize]. ISep: {:6}\n",
            pwgts[0], pwgts[1], graph.nvtxs, graph.nbnd, graph.mincut
        )
    );

    let badmaxpwgt = (ubfactor * pwgts[0].max(pwgts[1]) as real_t) as idx_t;

    for pass in (0)..(npasses) {
        moved.fill(-1);
        queues[0].reset();
        queues[1].reset();

        let mut mincutorder: idx_t = -1;
        let mut mincut = graph.mincut;
        let initcut = graph.mincut;
        let mut nbnd = graph.nbnd as usize;

        /* use the swaps array in place of the traditional perm array to save memory */
        irandArrayPermute(nbnd as idx_t, swaps.as_mut_ptr(), nbnd as idx_t, 1);
        for ii in (0)..(nbnd) {
            let i = bndind[swaps[ii] as usize] as usize;
            assert!(where_[i] == 2);
            if hmarker[i] == -1 {
                queues[0].insert(i as idx_t, (vwgt[i] - rinfo[i].edegrees[1]) as real_t);
                queues[1].insert(i as idx_t, (vwgt[i] - rinfo[i].edegrees[0]) as real_t);
                moved[i] = -5;
            } else if hmarker[i] != 2 {
                queues[hmarker[i] as usize].insert(
                    i as idx_t,
                    (vwgt[i] - rinfo[i].edegrees[(hmarker[i] as usize + 1) % 2]) as f32,
                );
                moved[i] = -(10 + hmarker[i]);
            }
        }

        debug_assert!(debug::CheckNodeBnd(graph, nbnd as idx_t) != 0);
        debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

        let limit = nbnd;

        /******************************************************
         * Get into the FM loop
         *******************************************************/
        mptr[0] = 0;
        let mut nmind = 0;
        let mut mindiff = (pwgts[0] - pwgts[1]).abs();
        // was never read
        // let mut to = if pwgts[0] < pwgts[1] { 0 } else { 1 };
        let mut to;
        let mut u_ = [-1, -1];
        let mut nswaps_ = usize::MAX;
        for nswaps in (0)..(nvtxs) {
            nswaps_ = nswaps;
            let mut g = [0, 0];
            match (queues[0].peek_val(), queues[1].peek_val()) {
                (Some(u_0), Some(u_1)) => {
                    let u = [u_0 as usize, u_1 as usize];
                    g[0] = vwgt[u[0]] - rinfo[u[0]].edegrees[1];
                    g[1] = vwgt[u[1]] - rinfo[u[1]].edegrees[0];

                    to = if g[0] > g[1] {
                        0
                    } else if g[0] < g[1] {
                        1
                    } else {
                        pass % 2
                    };

                    if pwgts[to as usize] + vwgt[u[to as usize] as usize] > badmaxpwgt {
                        to = (to + 1) % 2;
                    }
                }
                (None, None) => break,
                (Some(u_0), _) if pwgts[0] + vwgt[u_0 as usize] <= badmaxpwgt => to = 0,
                (_, Some(u_1)) if pwgts[1] + vwgt[u_1 as usize] <= badmaxpwgt => to = 1,
                _ => break,
            };
            u_[0] = queues[0].peek_val().unwrap_or(-1);
            u_[1] = queues[1].peek_val().unwrap_or(-1);

            let other = (to + 1) % 2;

            let higain = queues[to as usize].pop().unwrap();

            /* Delete its matching entry in the other queue */
            if moved[higain as usize] == -5 {
                queues[other as usize].delete(higain);
            }

            assert!(bndptr[higain as usize] != -1);

            /* The following check is to ensure we break out if there is a possibility
            of over-running the mind array.  */
            if nmind + xadj[(higain + 1) as usize] - xadj[higain as usize] >= 2 * nvtxs as idx_t - 1
            {
                break;
            }

            pwgts[2] -= vwgt[higain as usize] - rinfo[higain as usize].edegrees[other as usize];

            let newdiff = (pwgts[to as usize] + vwgt[higain as usize]
                - (pwgts[other as usize] - rinfo[higain as usize].edegrees[other as usize]))
                .abs();
            if pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff) {
                mincut = pwgts[2];
                mincutorder = nswaps as idx_t;
                mindiff = newdiff;
            } else {
                // TODO: it seems like this is not quite correct on first run when mincutorder = -1
                if nswaps as idx_t - mincutorder > limit as idx_t {
                    pwgts[2] +=
                        vwgt[higain as usize] - rinfo[higain as usize].edegrees[other as usize];
                    break; /* No further improvement, break out */
                }
            }

            BNDDelete!(nbnd, bndind, bndptr, higain as usize);
            pwgts[to as usize] += vwgt[higain as usize];
            where_[higain as usize] = to;
            moved[higain as usize] = nswaps as idx_t;
            swaps[nswaps as usize] = higain;

            /**********************************************************
             * Update the degrees of the affected nodes
             ***********************************************************/
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;
                if where_[k] == 2 {
                    /* For the in-separator vertices modify their edegree[to as usize] */
                    let oldgain = vwgt[k] - rinfo[k].edegrees[to as usize];
                    rinfo[k].edegrees[to as usize] += vwgt[higain as usize];
                    if moved[k] == -5 || moved[k] == -(10 + other) {
                        queues[other as usize]
                            .update(k as idx_t, (oldgain - vwgt[higain as usize]) as real_t);
                    }
                } else if where_[k] == other {
                    /* This vertex is pulled into the separator */
                    debug_assert!(bndptr[k] == -1, "{:} {:} {:}\n", k, bndptr[k], where_[k]);
                    BNDInsert!(nbnd, bndind, bndptr, k);

                    mind[(nmind) as usize] = k; /* Keep track for rollback */
                    nmind += 1;
                    where_[k] = 2;
                    pwgts[other as usize] -= vwgt[k];

                    let mut edegrees = rinfo[k].edegrees;
                    edegrees[0] = 0;
                    edegrees[1] = 0;
                    for jj in (xadj[k])..(xadj[(k + 1) as usize]) {
                        let kk = adjncy[jj as usize] as usize;
                        if where_[kk] != 2 {
                            edegrees[where_[kk] as usize] += vwgt[kk];
                        } else {
                            let oldgain = vwgt[kk] - rinfo[kk].edegrees[other as usize];
                            rinfo[kk].edegrees[other as usize] -= vwgt[k];
                            if moved[kk] == -5 || moved[kk] == -(10 + to) {
                                queues[to as usize]
                                    .update(kk as idx_t, (oldgain + vwgt[k]) as real_t);
                            }
                        }
                    }

                    /* Insert the new vertex into the priority queue (if it has not been moved). */
                    if moved[k] == -1 && (hmarker[k] == -1 || hmarker[k] == to) {
                        queues[to as usize]
                            .insert(k as idx_t, (vwgt[k] - edegrees[other as usize]) as real_t);
                        moved[k] = -(10 + to);
                    }
                    // #ifdef FULLMOVES  /* this does not work as well as the above partial one */
                    //           if (moved[k] == -1) {
                    //             if (hmarker[k] == -1) {
                    //               rpqInsert(queues[0], k, vwgt[k]-edegrees[1]);
                    //               rpqInsert(queues[1], k, vwgt[k]-edegrees[0]);
                    //               moved[k] = -5;
                    //             }
                    //             else if (hmarker[k] != 2) {
                    //               rpqInsert(queues[hmarker[k] as usize], k, vwgt[k]-edegrees[(hmarker[k]+1)%2]);
                    //               moved[k] = -(10+hmarker[k]);
                    //             }
                    //           }
                    // #endif
                }
            }
            mptr[(nswaps + 1) as usize] = nmind;

            // FIXME: This doesn't seem quite right
            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "Moved {:6} to {:3}, Gain: {:5} [{:5}] \
                    [{:4} {:4}] \t[{:5} {:5} {:5}]\n",
                    higain,
                    to,
                    g[to as usize],
                    g[other as usize],
                    vwgt[u_[to as usize] as usize],
                    vwgt[u_[other as usize] as usize],
                    pwgts[0],
                    pwgts[1],
                    pwgts[2]
                )
            );
        }

        /****************************************************************
         * Roll back computation
         *****************************************************************/
        // for (nswaps--; nswaps>mincutorder; nswaps--)
        let nswaps = if nswaps_ == usize::MAX {
            nvtxs as idx_t
        } else {
            nswaps_ as idx_t
        };
        for nswaps in ((mincutorder + 1)..nswaps).rev() {
            let higain = swaps[nswaps as usize] as usize;

            debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

            let to = where_[higain as usize] as usize;
            let other = (to + 1) % 2;
            inc_dec!(pwgts[2], pwgts[to as usize], vwgt[higain as usize]);
            where_[higain as usize] = 2;
            BNDInsert!(nbnd, bndind, bndptr, higain);

            let mut edegrees = rinfo[higain as usize].edegrees;
            edegrees[0] = 0;
            edegrees[1] = 0;
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;
                if where_[k] == 2 {
                    rinfo[k].edegrees[to as usize] -= vwgt[higain as usize];
                } else {
                    edegrees[where_[k] as usize] += vwgt[k];
                }
            }

            /* Push nodes out of the separator */
            for j in (mptr[nswaps as usize])..(mptr[(nswaps + 1) as usize]) {
                let k = mind[j as usize] as usize;
                debug_assert!(where_[k] == 2);
                where_[k] = other as idx_t;
                inc_dec!(pwgts[other as usize], pwgts[2], vwgt[k]);
                BNDDelete!(nbnd, bndind, bndptr, k);
                for jj in (xadj[k])..(xadj[(k + 1) as usize]) {
                    let kk = adjncy[jj as usize] as usize;
                    if where_[kk] == 2 {
                        rinfo[kk].edegrees[other as usize] += vwgt[k];
                    }
                }
            }
        }

        assert!(mincut == pwgts[2]);

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_REFINE,
            print!(
                "\tMinimum sep: {:6} at {:5}, PWGTS: [({:6} {:6}) as usize], NBND: {:6}\n",
                mincut, mincutorder, pwgts[0], pwgts[1], nbnd
            )
        );

        graph.mincut = mincut;
        graph.nbnd = nbnd as idx_t;

        if mincutorder == -1 || mincut >= initcut {
            break;
        }
    }

    // rpqDestroy(queues[0]);
    // rpqDestroy(queues[1]);

    // WCOREPOP;
}

/*************************************************************************/
/* This function computes a cache-friendly permutation of each partition.
    The resulting permutation is retuned in old2new, which is a vector of
    size nvtxs such for vertex i, old2new[i as usize] is its new vertex number.
*/
/**************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_CacheFriendlyReordering(
    nvtxs: idx_t,
    xadj: *const idx_t,
    adjncy: *const idx_t,
    part: *mut idx_t,
    old2new: *mut idx_t,
) -> c_int {
    // idx_t i, j, k, first, last, lastlevel, maxdegree, nparts;
    // idx_t *cot, *pos, *pwgts;
    // ikv_t *levels;

    mkslice!(xadj, nvtxs + 1);
    mkslice!(adjncy, xadj[nvtxs as usize]);
    mkslice_mut!(part, nvtxs);
    mkslice_mut!(old2new, nvtxs);

    util::InitRandom(123);

    /* This array ([C as usize]losed[O as usize]pen[T as usize]odo => cot) serves three purposes.
    Positions from [(0...first) is the current iperm[) as usize] vector of the explored vertices;
    Positions from [first...last) is the OPEN list (i.e., visited vertices);
    Positions from [last...nvtxs) is the todo list. */
    // cot = iincset(nvtxs, 0, imalloc(nvtxs, "METIS_CacheFriendlyReordering: cor"));
    let mut cot: Vec<idx_t> = (0..nvtxs).collect();

    /* This array will function like pos + touched of the CC method */
    // pos = iincset(nvtxs, 0, imalloc(nvtxs, "METIS_CacheFriendlyReordering: pos"));
    let mut pos: Vec<idx_t> = (0..nvtxs).collect();

    /* pick a random starting vertex */
    {
        let i = irandInRange(nvtxs);
        pos[0] = i;
        cot[0] = i;
        pos[i as usize] = 0;
        cot[i as usize] = 0;
    }

    let nvtxs = nvtxs as usize;

    /* compute a BFS ordering */
    let mut first: usize = 0;
    let mut last: usize = 0;
    let mut lastlevel: idx_t = 0;
    let mut maxdegree: idx_t = 0;
    while first < nvtxs {
        if first == last {
            /* Find another starting vertex */
            let k = cot[last] as usize;
            assert!(pos[k] >= 0);
            lastlevel -= 1;
            pos[k] = lastlevel; /* mark node as being visited by assigning its current level (-ve) */
            last += 1;
        }

        let i = cot[(first) as usize] as usize;
        first += 1;
        maxdegree = if maxdegree < xadj[(i + 1) as usize] - xadj[i as usize] {
            xadj[(i + 1) as usize] - xadj[i as usize]
        } else {
            maxdegree
        };
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            let k = adjncy[j as usize] as usize;
            /* if a node has been already been visited, its pos[] will be -1 */
            if pos[k] >= 0 {
                /* pos[k] is the location within cot of where_ k resides (it is in the 'todo' part);
                put in that location cot[last as usize] that we are about to overwrite
                and update pos[cot[last as usize] as usize] to reflect that. */
                cot[pos[k] as usize] = cot[last as usize];
                pos[cot[last as usize] as usize] = pos[k];

                cot[(last) as usize] = k as idx_t; /* put node at the end of the "queue" */
                last += 1;
                pos[k] = pos[i as usize] - 1; /* mark node as being visited by assigning to next level */
                lastlevel = pos[k]; /* for correctly advancing the levels in case of disconnected graphs */
            }
        }
    }
    //  print!("lastlevel: %d\n", (int)-lastlevel);

    /* sort based on decreasing level and decreasing degree (RCM) */
    // levels = ikvmalloc(nvtxs, "METIS_CacheFriendlyReordering: levels");
    // maxdegree += 1;
    // for i in (0)..(nvtxs) {
    //   levels[i as usize].val = i;
    //   levels[i as usize].key = -pos[i as usize]*maxdegree + xadj[(i+1) as usize]-xadj[i as usize];
    // }
    // ikvsortd(nvtxs, levels);

    maxdegree += 1;
    let mut levels: Vec<_> = (0..nvtxs)
        .map(|i| {
            (
                -pos[i as usize] * maxdegree + xadj[(i + 1) as usize] - xadj[i as usize],
                i as idx_t,
            )
        })
        .collect();
    levels.sort_unstable_by_key(|&(key, _val)| std::cmp::Reverse(key));

    /* figure out the partitions */
    // nparts = imax(nvtxs, part, 1)+1;
    // pwgts  = ismalloc(nparts+1, 0, "METIS_CacheFriendlyReordering: pwgts");
    let nparts = part.iter().copied().max().unwrap_or_default() as usize + 1;
    let mut pwgts = vec![0; nparts+1];

    for i in (0)..(nvtxs) {
        pwgts[part[i as usize] as usize] += 1;
    }
    util::make_csr(nparts, &mut pwgts);

    for i in (0)..(nvtxs) {
        let pwgt = &mut pwgts[part[levels[i].1 as usize] as usize];
        old2new[levels[i as usize].1 as usize] = *pwgt;
        *pwgt += 1;
    }

    // #ifdef XXX
    //   for i in (0)..(nvtxs) {
    //     for j in xadj[i as usize]..xadj[(i+1) as usize] {
    //       print!("COO: %d %d\n", (int)old2new[i as usize], (int)old2new[adjncy[j as usize] as usize]);
    //     }
    //   }
    // #endif

    // gk_free((void **)&cot, &pos, &levels, &pwgts, LTERM);

    return METIS_OK;
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use crate::{dyncall::ab_test_single_eq, tests::{ab_test_partition_test_graphs, for_test_suite}};
    use super::*;


    #[test]
    fn ab_METIS_CacheFriendlyReordering() {
        for_test_suite(Optype::Kmetis, 8, 1, |_tg, mut g| {
            g.random_adjwgt();
            let (_objval, mut part) = g.call().unwrap();
            let csr = g.into_csr();
            let nvtxs = csr.nvtxs();
            let (xadj, adjncy) = csr.as_parts();

            ab_test_single_eq("METIS_CacheFriendlyReordering:rs", || {
                let mut old2new = vec![idx_t::MAX; nvtxs];
                unsafe {
                    METIS_CacheFriendlyReordering(nvtxs as idx_t, xadj.as_ptr(), adjncy.as_ptr(), part.as_mut_ptr(), old2new.as_mut_ptr());
                }
                old2new
            });
        });
    }

    #[test]
    fn ab_METIS_NodeNDP() {
        // TODO: I realized I didn't do tests with varying nseps in ometis
        for_test_suite(Optype::Ometis, 3, 1, |_tg, mut g| {
            g.random_vwgt();
            g.set_nseps(3);
            ab_test_single_eq("METIS_NodeNDP:rs", || {
                g.call_ndp(8)
            });
        });

        for_test_suite(Optype::Ometis, 3, 1, |_tg, mut g| {
            g.random_vwgt();
            ab_test_single_eq("METIS_NodeNDP:rs", || {
                g.call_ndp(16)
            });
        });
    }

    #[test]
    fn ab_METIS_NodeNDP_nocompress() {
        for_test_suite(Optype::Ometis, 3, 1, |_tg, mut g| {
            g.random_vwgt();
            g.set_compress(false);
            g.set_nseps(3);
            ab_test_single_eq("METIS_NodeNDP:rs", || {
                g.call_ndp(8)
            });
        });
    }

    #[test]
    fn ab_METIS_ComputeVertexSeparator() {
        ab_test_partition_test_graphs("METIS_ComputeVertexSeparator:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g.compute_vertex_separator(true);
            g
        });
    }
}
