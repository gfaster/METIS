/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * sfm.c
 *
 * This file contains code that implements an FM-based separator refinement
 *
 * Started 8/1/97
 * George
 *
 * $Id: sfm.c 10874 2011-10-17 23:13:00Z karypis $
 *
 */

use crate::*;

/*************************************************************************/
/* This function performs a node-based FM refinement */
/**************************************************************************/
#[metis_func]
pub extern "C" fn FM_2WayNodeRefine2Sided(ctrl: *mut ctrl_t, graph: *mut graph_t, niter: idx_t) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    // idx_t i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind;
    // idx_t *xadj, *vwgt, *adjncy, *where_, *pwgts, *edegrees, *bndind, *bndptr;
    // idx_t *mptr, *mind, *moved, *swaps;
    // rpq_t *queues[2];
    // nrinfo_t *rinfo;
    // idx_t higain, oldgain, mincut, initcut, mincutorder;
    // idx_t pass, to, other, limit;
    // idx_t badmaxpwgt, mindiff, newdiff;
    // idx_t u[2], g[2];
    // real_t mult;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy vwgt);

    get_graph_slices_mut!(ctrl, graph => bndind bndptr where_ pwgts nrinfo);
    let rinfo = nrinfo;

    let mut queues = [pqueue::RPQueue::new(nvtxs), pqueue::RPQueue::new(nvtxs)];

    let mut moved: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut swaps: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut mptr: Vec<idx_t> = vec![0; nvtxs + 1 as usize];
    let mut mind: Vec<idx_t> = vec![0; 2 * nvtxs as usize];

    let mult = 0.5 * *ctrl.ubfactors;
    let badmaxpwgt = (mult * (pwgts[0] + pwgts[1] + pwgts[2]) as real_t) as idx_t;

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        print!(
            "Partitions-N2: [{:6} {:6}] Nv-Nb[{:6} {:6}]. ISep: {:6}\n",
            pwgts[0], pwgts[1], graph.nvtxs, graph.nbnd, graph.mincut
        )
    );

    for pass in (0)..(niter) {
        moved.fill(-1);
        queues[0].reset();
        queues[1].reset();

        // this is currently idx_t, but I'm not sure if it makes sense for this to be negative,
        // there is some potentially funky interections that don't expect that
        let mut mincutorder: idx_t = -1;
        let initcut = graph.mincut;
        let mut mincut = graph.mincut;
        let mut nbnd = graph.nbnd as usize;

        /* use the swaps array in place of the traditional perm array to save memory */
        // remember: 1 in the last argument means auto fill
        irandArrayPermute(nbnd as idx_t, swaps.as_mut_ptr(), nbnd as idx_t, 1);
        for ii in (0)..(nbnd) {
            let i = bndind[swaps[ii as usize] as usize];
            debug_assert!(where_[i as usize] == 2);
            queues[0].insert(
                i,
                (vwgt[i as usize] - rinfo[i as usize].edegrees[1]) as real_t,
            );
            queues[1].insert(
                i,
                (vwgt[i as usize] - rinfo[i as usize].edegrees[0]) as real_t,
            );
        }

        debug_assert!(debug::CheckNodeBnd(graph, nbnd as idx_t) != 0);
        debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

        let limit = if ctrl.compress != 0 {
            (5 * nbnd).min(400)
        } else {
            (2 * nbnd).min(300)
        };

        /******************************************************
         * Get into the FM loop
         *******************************************************/
        mptr[0] = 0;
        let mut nmind = 0;
        let mut mindiff = (pwgts[0] - pwgts[1]).abs();
        // these are declared here so they can be printed below
        let mut nswaps: usize = 0;
        // for nswaps in (0)..(nvtxs)
        while nswaps < nvtxs {
            let mut to;
            let mut g: [idx_t; 2] = [-1, -1];
            let u = [
                queues[0].peek_val().unwrap_or(-1),
                queues[1].peek_val().unwrap_or(-1),
            ];
            match (queues[0].peek_val(), queues[1].peek_val()) {
                (Some(u_0), Some(u_1)) => {
                    let u = [u_0 as usize, u_1 as usize];
                    g[0] = vwgt[u[0] as usize] - rinfo[u[0]].edegrees[1];
                    g[1] = vwgt[u[1] as usize] - rinfo[u[1]].edegrees[0];

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
            }

            let other = (to + 1) % 2;

            let higain = (queues[to as usize]).pop().unwrap() as usize;
            if moved[higain as usize] == -1 {
                /* Delete if it was in the separator originally */
                queues[other as usize].delete(higain as idx_t);
            }

            debug_assert!(bndptr[higain as usize] != -1);

            /* The following check is to ensure we break out if there is a possibility
            of over-running the mind array.  */
            if nmind + xadj[(higain + 1) as usize] as usize - xadj[higain as usize] as usize
                >= 2 * nvtxs - 1
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
                // FIXME: I think there's something funky here with initial value of mincutorder
                if nswaps as idx_t - mincutorder > 2 * limit as idx_t
                    || (nswaps as idx_t - mincutorder > limit as idx_t
                        && (pwgts[2] as real_t) > 1.10 * mincut as real_t)
                {
                    pwgts[2] +=
                        vwgt[higain as usize] - rinfo[higain as usize].edegrees[other as usize];
                    break; /* No further improvement, break out */
                }
            }

            BNDDelete!(nbnd, bndind, bndptr, higain);
            pwgts[to as usize] += vwgt[higain as usize];
            where_[higain as usize] = to;
            moved[higain as usize] = nswaps as idx_t;
            swaps[nswaps as usize] = higain as idx_t;

            /**********************************************************
             * Update the degrees of the affected nodes
             ***********************************************************/
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;
                if where_[k] == 2 {
                    /* For the in-separator vertices modify their edegree[to as usize] */
                    let oldgain = vwgt[k] - rinfo[k].edegrees[to as usize];
                    rinfo[k].edegrees[to as usize] += vwgt[higain as usize];
                    if moved[k] == -1 || moved[k] == -(2 + other) {
                        queues[other as usize]
                            .update(k as idx_t, (oldgain - vwgt[higain as usize]) as real_t);
                    }
                } else if where_[k] == other {
                    /* This vertex is pulled into the separator */
                    debug_assert!(bndptr[k] == -1, "{:} {:} {:}", k, bndptr[k], where_[k]);
                    BNDInsert!(nbnd, bndind, bndptr, k);

                    mind[(nmind) as usize] = k as idx_t; /* Keep track for rollback */
                    nmind += 1;
                    where_[k] = 2;
                    pwgts[other as usize] -= vwgt[k];

                    // edegrees = rinfo[k].edegrees;
                    rinfo[k].edegrees = [0, 0];
                    for jj in (xadj[k])..(xadj[(k + 1) as usize]) {
                        let kk = adjncy[jj as usize] as usize;
                        if where_[kk] != 2 {
                            rinfo[k].edegrees[where_[kk] as usize] += vwgt[kk];
                        } else {
                            let oldgain = vwgt[kk] - rinfo[kk].edegrees[other as usize];
                            rinfo[kk].edegrees[other as usize] -= vwgt[k];
                            if moved[kk] == -1 || moved[kk] == -(2 + to) {
                                queues[to as usize]
                                    .update(kk as idx_t, (oldgain + vwgt[k]) as real_t);
                            }
                        }
                    }

                    /* Insert the new vertex into the priority queue. Only one side! */
                    if moved[k] == -1 {
                        queues[to as usize].insert(
                            k as idx_t,
                            (vwgt[k] - rinfo[k].edegrees[other as usize]) as real_t,
                        );
                        moved[k] = -(2 + to);
                    }
                }
            }
            mptr[(nswaps + 1) as usize] = nmind as idx_t;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "Moved {:6} to {:3}, Gain: {:5} [{:5}] [{:4} {:4}] \t[{:5} {:5} {:5}]\n",
                    higain,
                    to,
                    g[to as usize],
                    g[other as usize],
                    vwgt[u[to as usize] as usize],
                    vwgt[u[other as usize] as usize],
                    pwgts[0],
                    pwgts[1],
                    pwgts[2]
                )
            );

            nswaps += 1;
        }

        /****************************************************************
         * Roll back computation
         *****************************************************************/
        // for (nswaps--; nswaps>mincutorder; nswaps--)
        for nswaps in ((mincutorder + 1) as usize..nswaps).rev() {
            let higain = swaps[nswaps as usize] as usize;

            debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

            let to = where_[higain as usize];
            let other = (to + 1) % 2;
            inc_dec!(pwgts[2], pwgts[to as usize], vwgt[higain as usize]);
            where_[higain as usize] = 2;
            BNDInsert!(nbnd, bndind, bndptr, higain);

            // edegrees = rinfo[higain as usize].edegrees;
            rinfo[higain as usize].edegrees = [0, 0];
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;
                if where_[k] == 2 {
                    rinfo[k].edegrees[to as usize] -= vwgt[higain as usize];
                } else {
                    rinfo[higain].edegrees[where_[k] as usize] += vwgt[k];
                }
            }

            /* Push nodes out of the separator */
            for j in (mptr[nswaps as usize])..(mptr[(nswaps + 1) as usize]) {
                let k = mind[j as usize] as usize;
                debug_assert!(where_[k] == 2);
                where_[k] = other;
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

        debug_assert!(mincut == pwgts[2]);

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_REFINE,
            print!(
                "\tMinimum sep: {:6} at {:5}, PWGTS: [{:6} {:6}], NBND: {:6}\n",
                mincut, mincutorder, pwgts[0], pwgts[1], nbnd
            )
        );

        graph.mincut = mincut;
        graph.nbnd = nbnd as idx_t;

        if mincutorder == -1 || mincut >= initcut {
            break;
        }
    }

    // WCOREPOP;
}

/*************************************************************************/
/* This function performs a node-based FM refinement.
    Each refinement iteration is split into two sub-iterations.
    In each sub-iteration only moves to one of the left/right partitions
    is allowed; hence, it is one-sided.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn FM_2WayNodeRefine1Sided(ctrl: *mut ctrl_t, graph: *mut graph_t, niter: idx_t) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    // idx_t i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind, iend;
    // idx_t *xadj, *vwgt, *adjncy, *where_, *pwgts, *edegrees, *bndind, *bndptr;
    // idx_t *mptr, *mind, *swaps;
    // rpq_t *queue;
    // nrinfo_t *rinfo;
    // idx_t higain, mincut, initcut, mincutorder;
    // idx_t pass, to, other, limit;
    // idx_t badmaxpwgt, mindiff, newdiff;
    // real_t mult;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy vwgt);

    get_graph_slices_mut!(ctrl, graph => bndind bndptr where_ pwgts nrinfo);
    let rinfo = nrinfo;

    let mut queue = pqueue::RPQueue::new(nvtxs);

    let mut swaps: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut mptr: Vec<idx_t> = vec![0; nvtxs + 1 as usize];
    let mut mind: Vec<idx_t> = vec![0; 2 * nvtxs as usize];

    let mult = 0.5 * *ctrl.ubfactors;
    let badmaxpwgt = (mult * (pwgts[0] + pwgts[1] + pwgts[2]) as real_t) as idx_t;

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        print!(
            "Partitions-N1: [{:6} {:6}] Nv-Nb[{:6} {:6}]. ISep: {:6}\n",
            pwgts[0], pwgts[1], graph.nvtxs, graph.nbnd, graph.mincut
        )
    );

    let mut to = if pwgts[0] < pwgts[1] { 1 } else { 0 };
    for pass in (0)..(2 * niter) {
        /* the 2*niter is for the two sides */
        let other = to;
        to = (to + 1) % 2;

        queue.reset();

        let mut mincutorder: idx_t = -1;
        let initcut = graph.mincut;
        let mut mincut = graph.mincut;
        let mut nbnd = graph.nbnd as usize;

        /* use the swaps array in place of the traditional perm array to save memory */
        // remember 1 in flag populats slice
        irandArrayPermute(nbnd as idx_t, swaps.as_mut_ptr(), nbnd as idx_t, 1);
        for ii in (0)..(nbnd) {
            let i = bndind[swaps[ii as usize] as usize] as usize;
            debug_assert!(where_[i as usize] == 2);
            queue.insert(
                i as idx_t,
                (vwgt[i as usize] - rinfo[i as usize].edegrees[other as usize]) as real_t,
            );
        }

        debug_assert!(debug::CheckNodeBnd(graph, nbnd as idx_t) != 0);
        debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

        let limit = if ctrl.compress != 0 {
            (5 * nbnd).min(500)
        } else {
            (3 * nbnd).min(300)
        };

        /******************************************************
         * Get into the FM loop
         *******************************************************/
        // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.Aux3Tmr));
        mptr[0] = 0;
        let mut nmind: usize = 0;
        let mut mindiff = (pwgts[0] - pwgts[1]).abs();
        let mut nswaps = 0;
        // for nswaps in (0)..(nvtxs)
        while nswaps < nvtxs {
            let Some(higain) = queue.pop() else { break };
            let higain = higain as usize;

            debug_assert!(bndptr[higain as usize] != -1);

            /* The following check is to ensure we break out if there is a possibility
            of over-running the mind array.  */
            if nmind + xadj[higain + 1] as usize - xadj[higain as usize] as usize >= 2 * nvtxs - 1 {
                break;
            }

            if pwgts[to as usize] + vwgt[higain as usize] > badmaxpwgt {
                /* No point going any further. Balance will be bad */
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
            } else if nswaps as idx_t - mincutorder > 3 * limit as idx_t
                || (nswaps as idx_t - mincutorder > limit as idx_t
                    && (pwgts[2] as real_t) > 1.10 * mincut as real_t)
            {
                pwgts[2] +=
                    vwgt[higain as usize] - rinfo[higain as usize].edegrees[other as usize];
                break; /* No further improvement, break out */
            }

            BNDDelete!(nbnd, bndind, bndptr, higain);
            pwgts[to as usize] += vwgt[higain as usize];
            where_[higain as usize] = to;
            swaps[nswaps as usize] = higain as idx_t;

            /**********************************************************
             * Update the degrees of the affected nodes
             ***********************************************************/
            // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.Aux1Tmr));
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;

                if where_[k] == 2 {
                    /* For the in-separator vertices modify their edegree[to as usize] */
                    rinfo[k].edegrees[to as usize] += vwgt[higain as usize];
                } else if where_[k] == other {
                    /* This vertex is pulled into the separator */
                    debug_assert!(bndptr[k] == -1, "{:} {:} {:}", k, bndptr[k], where_[k]);
                    BNDInsert!(nbnd, bndind, bndptr, k);

                    mind[nmind as usize] = k as idx_t; /* Keep track for rollback */
                    nmind += 1;
                    where_[k] = 2;
                    pwgts[other as usize] -= vwgt[k];

                    // edegrees = rinfo[k].edegrees;
                    rinfo[k].edegrees = [0, 0];
                    for jj in xadj[k]..xadj[(k + 1) as usize] {
                        let kk = adjncy[jj as usize] as usize;
                        if where_[kk] != 2 {
                            rinfo[k].edegrees[where_[kk] as usize] += vwgt[kk];
                        } else {
                            rinfo[kk].edegrees[other as usize] -= vwgt[k];

                            /* Since the moves are one-sided this vertex has not been moved yet */
                            queue.update(
                                kk as idx_t,
                                (vwgt[kk] - rinfo[kk].edegrees[other as usize]) as real_t,
                            );
                        }
                    }

                    /* Insert the new vertex into the priority queue. Safe due to one-sided moves */
                    queue.insert(
                        k as idx_t,
                        (vwgt[k] - rinfo[k].edegrees[other as usize]) as real_t,
                    );
                }
            }
            mptr[(nswaps + 1) as usize] = nmind as idx_t;
            // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.Aux1Tmr));

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "Moved {:6} to {:3}, Gain: {:5} [{:5}] \t[{:5} {:5} {:5}] [{:3} {:2}]\n",
                    higain,
                    to,
                    (vwgt[higain as usize] - rinfo[higain as usize].edegrees[other as usize]),
                    vwgt[higain as usize],
                    pwgts[0],
                    pwgts[1],
                    pwgts[2],
                    nswaps,
                    limit
                )
            );

            nswaps += 1;
        }
        // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.Aux3Tmr));

        /****************************************************************
         * Roll back computation
         *****************************************************************/
        // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.Aux2Tmr));
        // for (nswaps--; nswaps>mincutorder; nswaps--)
        for nswaps in ((mincutorder + 1) as usize..nswaps).rev() {
            let higain = swaps[nswaps as usize] as usize;

            debug_assert!(debug::CheckNodePartitionParams(graph) != 0);
            debug_assert!(where_[higain as usize] == to);

            inc_dec!(pwgts[2], pwgts[to as usize], vwgt[higain as usize]);
            where_[higain as usize] = 2;
            BNDInsert!(nbnd, bndind, bndptr, higain);

            // edegrees = rinfo[higain as usize].edegrees;
            rinfo[higain as usize].edegrees = [0, 0];
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;
                if where_[k] == 2 {
                    rinfo[k].edegrees[to as usize] -= vwgt[higain as usize];
                } else {
                    rinfo[higain].edegrees[where_[k] as usize] += vwgt[k];
                }
            }

            /* Push nodes out of the separator */
            for j in (mptr[nswaps as usize])..(mptr[(nswaps + 1) as usize]) {
                let k = mind[j as usize] as usize;
                debug_assert!(where_[k] == 2);
                where_[k] = other;
                inc_dec!(pwgts[other as usize], pwgts[2], vwgt[k]);
                BNDDelete!(nbnd, bndind, bndptr, k);
                for jj in xadj[k]..xadj[(k + 1) as usize] {
                    let kk = adjncy[jj as usize] as usize;
                    if where_[kk] == 2 {
                        rinfo[kk].edegrees[other as usize] += vwgt[k];
                    }
                }
            }
        }
        // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.Aux2Tmr));

        debug_assert!(mincut == pwgts[2]);

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_REFINE,
            print!(
                "\tMinimum sep: {:6} at {:5}, PWGTS: [{:6} {:6}], NBND: {:6}\n",
                mincut, mincutorder, pwgts[0], pwgts[1], nbnd
            )
        );

        graph.mincut = mincut;
        graph.nbnd = nbnd as idx_t;

        if pass % 2 == 1 && (mincutorder == -1 || mincut >= initcut) {
            break;
        }
    }

    // WCOREPOP;
}

/*************************************************************************/
/* This function balances the left/right partitions of a separator
tri-section */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FM_2WayNodeBalance(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, gain;
    // idx_t badmaxpwgt, higain, oldgain, pass, to, other;
    // idx_t *xadj, *vwgt, *adjncy, *where_, *pwgts, *edegrees, *bndind, *bndptr;
    // idx_t *perm, *moved;
    // rpq_t *queue;
    // nrinfo_t *rinfo;
    // real_t mult;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy vwgt);

    get_graph_slices_mut!(ctrl, graph => bndind bndptr where_ pwgts nrinfo);
    let rinfo = nrinfo;

    let mult = 0.5 * *ctrl.ubfactors;

    let badmaxpwgt = (mult * (pwgts[0] + pwgts[1]) as real_t) as idx_t;
    if (pwgts[0]).max(pwgts[1]) < badmaxpwgt {
        return;
    }
    if (pwgts[0] - pwgts[1]).abs() < 3 * *graph.tvwgt / nvtxs as idx_t {
        return;
    }

    // WCOREPUSH;

    let to = if pwgts[0] < pwgts[1] { 0 } else { 1 };
    let other = (to + 1) % 2;

    let mut queue = pqueue::RPQueue::new(nvtxs);

    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut moved: Vec<idx_t> = vec![-1; nvtxs as usize];

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        print!(
            "Partitions: [{:6} {:6}] Nv-Nb[{:6} {:6}]. ISep: {:6} [B]\n",
            pwgts[0], pwgts[1], graph.nvtxs, graph.nbnd, graph.mincut
        )
    );

    let mut nbnd = graph.nbnd as usize;
    // putting 1 in for the flag assigns incrementally
    irandArrayPermute(nbnd as idx_t, perm.as_mut_ptr(), nbnd as idx_t, 1);
    for ii in (0)..(nbnd) {
        let i = bndind[perm[ii as usize] as usize] as usize;
        debug_assert!(where_[i as usize] == 2);
        queue.insert(
            i as idx_t,
            (vwgt[i as usize] - rinfo[i as usize].edegrees[other as usize]) as real_t,
        );
    }

    debug_assert!(debug::CheckNodeBnd(graph, nbnd as idx_t) != 0);
    debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

    /******************************************************
     * Get into the FM loop
     *******************************************************/

    // for nswaps in (0)..(nvtxs)
    let mut nswaps: usize = 0;
    while nswaps < nvtxs {
        let Some(higain) = queue.pop() else { break };
        let higain = higain as usize;

        moved[higain as usize] = 1;

        let gain = vwgt[higain as usize] - rinfo[higain as usize].edegrees[other as usize];
        let badmaxpwgt = (mult * (pwgts[0] + pwgts[1]) as real_t) as idx_t;

        /* break if other is now underwight */
        if pwgts[to as usize] > pwgts[other as usize] {
            break;
        }

        /* break if balance is achieved and no +ve or zero gain */
        if gain < 0 && pwgts[other as usize] < badmaxpwgt {
            break;
        }

        /* skip this vertex if it will violate balance on the other side */
        if pwgts[to as usize] + vwgt[higain as usize] > badmaxpwgt {
            // I don't think this does anything except change the debug output
            nswaps += 1;
            continue;
        }

        debug_assert!(bndptr[higain as usize] != -1);

        pwgts[2] -= gain;

        BNDDelete!(nbnd, bndind, bndptr, higain);
        pwgts[to as usize] += vwgt[higain as usize];
        where_[higain as usize] = to;

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_MOVEINFO,
            print!(
                "Moved {:6} to {:3}, Gain: {:3}, \t[{:5} {:5} {:5}]\n",
                higain,
                to,
                vwgt[higain as usize] - rinfo[higain as usize].edegrees[other as usize],
                pwgts[0],
                pwgts[1],
                pwgts[2]
            )
        );

        /**********************************************************
         * Update the degrees of the affected nodes
         ***********************************************************/
        for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
            let k = adjncy[j as usize] as usize;
            if where_[k] == 2 {
                /* For the in-separator vertices modify their edegree[to as usize] */
                rinfo[k].edegrees[to as usize] += vwgt[higain as usize];
            } else if where_[k] == other {
                /* This vertex is pulled into the separator */
                debug_assert!(bndptr[k] == -1, "{:} {:} {:}", k, bndptr[k], where_[k]);
                BNDInsert!(nbnd, bndind, bndptr, k);

                where_[k] = 2;
                pwgts[other as usize] -= vwgt[k];

                rinfo[k].edegrees = [0, 0];
                for jj in (xadj[k])..(xadj[(k + 1) as usize]) {
                    let kk = adjncy[jj as usize] as usize;
                    if where_[kk] != 2 {
                        rinfo[k].edegrees[where_[kk] as usize] += vwgt[kk];
                    } else {
                        debug_assert!(bndptr[kk] != -1);
                        let oldgain = vwgt[kk] - rinfo[kk].edegrees[other as usize];
                        rinfo[kk].edegrees[other as usize] -= vwgt[k];

                        if moved[kk] == -1 {
                            queue.update(kk as idx_t, (oldgain + vwgt[k]) as real_t);
                        }
                    }
                }

                /* Insert the new vertex into the priority queue */
                queue.insert(
                    k as idx_t,
                    (vwgt[k] - rinfo[k].edegrees[other as usize]) as real_t,
                );
            }
        }

        nswaps += 1;
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        print!(
            "\tBalanced sep: {:6} at {:4}, PWGTS: [{:6} {:6}], NBND: {:6}\n",
            pwgts[2], nswaps, pwgts[0], pwgts[1], nbnd
        )
    );

    graph.mincut = pwgts[2];
    graph.nbnd = nbnd as idx_t;

    // WCOREPOP;
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use super::*;
    use crate::tests::ab_test_partition_test_graphs;

    #[test]
    fn ab_FM_2WayNodeRefine2Sided() {
        ab_test_partition_test_graphs(
            "FM_2WayNodeRefine2Sided:rs",
            Optype::Ometis,
            3,
            1,
            |mut g| {
                g.random_vwgt();
                g.set_rtype(graph_gen::OmetisRtype::Sep2Sided);
                g
            },
        );
    }

    #[test]
    fn ab_FM_2WayNodeRefine1Sided() {
        ab_test_partition_test_graphs(
            "FM_2WayNodeRefine1Sided:rs",
            Optype::Ometis,
            3,
            1,
            |mut g| {
                g.random_vwgt();
                g.set_rtype(graph_gen::OmetisRtype::Sep1Sided);
                g
            },
        );
    }

    #[test]
    fn ab_FM_2WayNodeBalance() {
        ab_test_partition_test_graphs("FM_2WayNodeBalance:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g.set_rtype(graph_gen::OmetisRtype::Sep1Sided);
            g
        });
    }
}
