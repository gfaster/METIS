/*
\file
\brief Functions for the edge-based FM refinement

\date Started 7/23/97
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota
\version\verbatim $Id: fm.c 10187 2011-06-13 13:46:57Z karypis $ \endverbatim
*/

use crate::pqueue::{IPQueue, RPQueue};
use crate::*;

/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
#[metis_func]
pub extern "C" fn FM_2WayRefine(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niter: idx_t,
) {
    if (*graph).ncon == 1 {
        FM_2WayCutRefine(ctrl, graph, ntpwgts, niter);
    } else {
        FM_Mc2WayCutRefine(ctrl, graph, ntpwgts, niter);
    }
}

/*************************************************************************/
/* This function performs a cut-focused FM refinement */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FM_2WayCutRefine(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niter: idx_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, pass, me, limit, tmp;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_, *id, *ed, *bndptr, *bndind, *pwgts;
    // idx_t *moved, *swaps, *perm;
    // rpq_t *queues[2 as usize];
    // idx_t higain, mincut, mindiff, origdiff, initcut, newcut, mincutorder, avgvwgt;
    // idx_t tpwgts[2 as usize];

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj vwgt adjncy tvwgt adjwgt);
    get_graph_slices_mut!(ctrl, graph => where_ id ed pwgts bndptr bndind);
    mkslice!(ntpwgts, 2 * graph.ncon);

    let mut moved: Vec<idx_t> = vec![-1; nvtxs as usize as usize];
    let mut swaps: Vec<idx_t> = vec![0; nvtxs as usize as usize];
    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize as usize];

    let tpwgts = {
        let tpwgt0 = (tvwgt[0 as usize] as real_t * ntpwgts[0 as usize]) as idx_t;
        [tpwgt0, tvwgt[0 as usize] - tpwgt0]
    };

    // let limit = gk_min(gk_max(0.01 * nvtxs, 15), 100);
    let limit = (0.01 * nvtxs as f64).clamp(15.0, 100.0) as idx_t;
    let avgvwgt = ((pwgts[0 as usize] + pwgts[1 as usize]) / 20)
        .min(2 * (pwgts[0 as usize] + pwgts[1 as usize]) / nvtxs as idx_t);

    // queues[0 as usize] = rpqCreate(nvtxs);
    // queues[1 as usize] = rpqCreate(nvtxs);
    let mut queues: [_; 2] = std::array::from_fn(|_| RPQueue::new(nvtxs));

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        Print2WayRefineStats(ctrl, graph, ntpwgts.as_ptr(), 0.0, -2),
    );

    let origdiff = (tpwgts[0 as usize] - pwgts[0 as usize]).abs();
    // iset(nvtxs, -1, moved);
    for _pass in (0)..(niter) {
        /* Do a number of passes */
        // rpqReset(queues[0 as usize]);
        // rpqReset(queues[1 as usize]);
        queues[0].reset();
        queues[1].reset();

        let mut mincutorder: idx_t = -1;
        let mut newcut = graph.mincut;
        let initcut = graph.mincut;
        let mut mincut = graph.mincut;
        // newcut = mincut = initcut = graph.mincut;
        let mut mindiff = (tpwgts[0 as usize] - pwgts[0 as usize]).abs();

        debug_assert_eq!(debug::ComputeCut(graph, where_.as_ptr()), graph.mincut);
        debug_assert!(debug::CheckBnd(graph) != 0);

        /* Insert boundary nodes in the priority queues */
        let mut nbnd = graph.nbnd as usize;
        irandArrayPermute(nbnd as idx_t, perm.as_mut_ptr(), nbnd as idx_t, 1);
        for ii in (0)..(nbnd) {
            let i = perm[ii as usize] as usize;
            assert!(ed[bndind[i as usize] as usize] > 0 || id[bndind[i as usize] as usize] == 0);
            assert!(bndptr[bndind[i as usize] as usize] != -1);
            // rpqInsert(
            //     queues[where_[bndind[i as usize] as usize]],
            //     bndind[i as usize],
            //     ed[bndind[i as usize] as usize] - id[bndind[i as usize]],
            // );
            queues[where_[bndind[i as usize] as usize] as usize].insert(
                bndind[i],
                (ed[bndind[i as usize] as usize] - id[bndind[i as usize] as usize]) as real_t,
            );
        }

        let mut nswaps: idx_t = 0;
        // for nswaps in (0)..(nvtxs) {
        while nswaps < nvtxs as idx_t {
            let from = if tpwgts[0 as usize] - pwgts[0 as usize]
                < tpwgts[1 as usize] - pwgts[1 as usize]
            {
                0
            } else {
                1
            };
            let to = (from + 1) % 2;

            // if ((higain = rpqGetTop(queues[from as usize])) == -1) {
            //     break;
            // }
            let Some(higain) = queues[from].pop() else {
                break;
            };
            let higain = higain as usize;
            assert!(bndptr[higain as usize] != -1);

            newcut -= ed[higain as usize] - id[higain as usize];
            inc_dec!(
                pwgts[to as usize],
                pwgts[from as usize],
                vwgt[higain as usize]
            );

            if (newcut < mincut
                && (tpwgts[0 as usize] - pwgts[0 as usize]).abs() <= origdiff + avgvwgt)
                || (newcut == mincut && (tpwgts[0 as usize] - pwgts[0 as usize]).abs() < mindiff)
            {
                mincut = newcut;
                mindiff = (tpwgts[0 as usize] - pwgts[0 as usize]).abs();
                mincutorder = nswaps as idx_t;
            } else if nswaps - mincutorder > limit {
                /* We hit the limit, undo last move */
                // dead assignment
                // newcut += ed[higain as usize] - id[higain as usize];
                inc_dec!(
                    pwgts[from as usize],
                    pwgts[to as usize],
                    vwgt[higain as usize]
                );
                break;
            }

            where_[higain as usize] = to as idx_t;
            moved[higain as usize] = nswaps;
            swaps[nswaps as usize] = higain as idx_t;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                println!(
                    "Moved {:6} from {:}. [({:3} {:3})] {:5} [({:4} {:4})]",
                    higain,
                    from,
                    ed[higain as usize] - id[higain as usize],
                    vwgt[higain as usize],
                    newcut,
                    pwgts[0 as usize],
                    pwgts[1 as usize],
                ),
            );

            /**************************************************************
             * Update the id[i as usize]/ed[i as usize] values of the affected nodes
             ***************************************************************/
            std::mem::swap(&mut id[higain as usize], &mut ed[higain as usize]);
            if ed[higain as usize] == 0 && xadj[higain as usize] < xadj[(higain + 1) as usize] {
                BNDDelete!(nbnd, bndind, bndptr, higain);
            }

            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;

                let kwgt = if to as idx_t == where_[k as usize] {
                    adjwgt[j as usize]
                } else {
                    -adjwgt[j as usize]
                };
                inc_dec!(id[k as usize], ed[k as usize], kwgt);

                /* Update its boundary information and queue position */
                if bndptr[k as usize] != -1 {
                    /* If k was a boundary vertex */
                    if ed[k as usize] == 0 {
                        /* Not a boundary vertex any more */
                        BNDDelete!(nbnd, bndind, bndptr, k);
                        if moved[k as usize] == -1
                        /* Remove it if in the queues */
                        {
                            // rpqDelete(queues[where_[k as usize] as usize], k);
                            (queues[where_[k as usize] as usize].delete(k as idx_t));
                        }
                    } else {
                        /* If it has not been moved, update its position in the queue */
                        if moved[k as usize] == -1 {
                            // rpqUpdate(queues[where_[k as usize] as usize], k, ed[k as usize] - id[k as usize]);
                            queues[where_[k as usize] as usize]
                                .update(k as idx_t, (ed[k as usize] - id[k as usize]) as real_t);
                        }
                    }
                } else {
                    if ed[k as usize] > 0 {
                        /* It will now become a boundary vertex */
                        BNDInsert!(nbnd, bndind, bndptr, k);
                        if moved[k as usize] == -1 {
                            queues[where_[k as usize] as usize]
                                .insert(k as idx_t, (ed[k as usize] - id[k as usize]) as real_t);
                        }
                    }
                }
            }
            nswaps += 1;
        }

        /****************************************************************
         * Roll back computations
         *****************************************************************/
        for i in (0)..(nswaps) {
            moved[swaps[i as usize] as usize] = -1; /* reset moved array */
        }
        nswaps -= 1;
        // for (nswaps--; nswaps>mincutorder; nswaps--) {
        while nswaps as idx_t > mincutorder {
            let higain = swaps[nswaps as usize] as usize;
            let to = (where_[higain as usize] + 1) % 2;
            where_[higain as usize] = to;
            std::mem::swap(&mut id[higain as usize], &mut ed[higain as usize]);
            if ed[higain as usize] == 0
                && bndptr[higain as usize] != -1
                && xadj[higain as usize] < xadj[(higain + 1) as usize]
            {
                BNDDelete!(nbnd, bndind, bndptr, higain);
            } else if ed[higain as usize] > 0 && bndptr[higain as usize] == -1 {
                BNDInsert!(nbnd, bndind, bndptr, higain);
            }

            inc_dec!(
                pwgts[to as usize],
                pwgts[((to + 1) % 2) as usize],
                vwgt[higain as usize]
            );
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize] as usize;

                let kwgt = if to == where_[k as usize] {
                    adjwgt[j as usize]
                } else {
                    -adjwgt[j as usize]
                };
                inc_dec!(id[k as usize], ed[k as usize], kwgt);

                if bndptr[k as usize] != -1 && ed[k as usize] == 0 {
                    BNDDelete!(nbnd, bndind, bndptr, k);
                }
                if bndptr[k as usize] == -1 && ed[k as usize] > 0 {
                    BNDInsert!(nbnd, bndind, bndptr, k);
                }
            }
            nswaps -= 1;
        }

        graph.mincut = mincut;
        graph.nbnd = nbnd as idx_t;

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_REFINE,
            Print2WayRefineStats(ctrl, graph, ntpwgts.as_ptr(), 0.0, mincutorder),
        );

        if mincutorder <= 0 || mincut == initcut {
            break;
        }
    }

    // rpqDestroy(queues[0 as usize]);
    // rpqDestroy(queues[1 as usize]);

    // WCOREPOP;
}

/*************************************************************************/
/* This function performs a cut-focused multi-constraint FM refinement */
/*************************************************************************/
#[metis_func]
pub extern "C" fn FM_Mc2WayCutRefine(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niter: idx_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, l, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, pass,
    // me, limit, tmp, cnum;
    // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *pwgts, *where_, *id, *ed,
    // *bndptr, *bndind;
    // idx_t *moved, *swaps, *perm, *qnum;
    // idx_t higain, mincut, initcut, newcut, mincutorder;
    // real_t *invtvwgt, *ubfactors, *minbalv, *newbalv;
    // real_t origbal, minbal, newbal, rgain, ffactor;
    // rpq_t **queues;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    // xadj = graph.xadj;
    // vwgt = graph.vwgt;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // invtvwgt = graph.invtvwgt;
    // where_ = graph.where_;
    // id = graph.id;
    // ed = graph.ed;
    // pwgts = graph.pwgts;
    // bndptr = graph.bndptr;
    // bndind = graph.bndind;
    get_graph_slices!(ctrl, graph => xadj vwgt adjncy adjwgt invtvwgt);
    get_graph_slices_mut!(ctrl, graph => where_ id ed pwgts bndptr bndind);

    let mut moved: Vec<idx_t> = vec![-1; nvtxs as usize as usize];
    let mut swaps: Vec<idx_t> = vec![0; nvtxs as usize as usize];
    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize as usize];
    let mut qnum: Vec<idx_t> = vec![0; nvtxs as usize as usize];
    let mut ubfactors = vec![0.0; ncon as usize as usize];
    let mut newbalv = vec![0.0; ncon as usize as usize];
    let mut minbalv = vec![0.0; ncon as usize as usize];
    mkslice!(ctrl->pijbm, ctrl.nparts as usize * ncon);

    // let limit = gk_min(gk_max(0.01 * nvtxs, 25), 150);
    // let limit = (nvtxs / 100).clamp(25, 150);
    let limit = (nvtxs as f64 * 0.01).clamp(25.0, 150.0) as idx_t;

    /* Determine a fudge factor to allow the refinement routines to get out
    of tight balancing constraints. */
    let ffactor = (0.5 / 20.0_f64.max(nvtxs as f64)) as real_t;

    /* Initialize the queues */
    // queues = (rpq_t **)wspacemalloc(ctrl, 2*ncon*sizeof(rpq_t *));
    for i in (0)..(nvtxs) {
        qnum[i as usize] = util::iargmax_nrm(ncon, vwgt[(i * ncon)..].as_ptr(), invtvwgt.as_ptr());
    }
    let mut queues = Vec::from_iter(std::iter::repeat_with(|| RPQueue::new(nvtxs)).take(2 * ncon));

    /* Determine the unbalance tolerance for each constraint. The tolerance is
    equal to the maximum of the original load imbalance and the user-supplied
    allowed tolerance. The rationale behind this approach is to allow the
    refinement routine to improve the cut, without having to worry about fixing
    load imbalance problems. The load imbalance is addressed by the balancing
    routines. */
    let origbal = mcutil::ComputeLoadImbalanceDiffVec(
        graph,
        2,
        ctrl.pijbm,
        ctrl.ubfactors,
        ubfactors.as_mut_ptr(),
    );
    for i in (0)..(ncon) {
        ubfactors[i as usize] = if ubfactors[i as usize] > 0.0 {
            *ctrl.ubfactors.add(i as usize) + ubfactors[i as usize]
        } else {
            *ctrl.ubfactors.add(i as usize)
        };
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        Print2WayRefineStats(ctrl, graph, ntpwgts, origbal, -2),
    );

    // iset(nvtxs, -1, moved);
    for _pass in (0)..(niter) {
        /* Do a number of passes */
        for i in (0)..(2 * ncon) {
            queues[i as usize].reset();
        }

        let mut mincutorder: idx_t = -1;
        let mut newcut = graph.mincut;
        let mut mincut = graph.mincut;
        let initcut = graph.mincut;
        // newcut = mincut = initcut = graph.mincut;

        let mut minbal = mcutil::ComputeLoadImbalanceDiffVec(
            graph,
            2,
            ctrl.pijbm,
            ubfactors.as_ptr(),
            minbalv.as_mut_ptr(),
        );

        assert!(debug::ComputeCut(graph, where_.as_ptr()) == graph.mincut);
        assert!(debug::CheckBnd(graph) != 0);

        /* Insert boundary nodes in the priority queues */
        let mut nbnd = graph.nbnd as usize;
        irandArrayPermute(nbnd as idx_t, perm.as_mut_ptr(), nbnd as idx_t / 5, 1);
        for ii in (0)..(nbnd) {
            let i = bndind[perm[ii as usize] as usize] as usize;
            assert!(ed[i as usize] > 0 || id[i as usize] == 0);
            assert!(bndptr[i as usize] != -1);
            //rgain = 1.0*(ed[i as usize]-id[i as usize])/sqrt(vwgt[(i*ncon+qnum[i) as usize] as usize]+1);
            //rgain = (if ed[i as usize]-id[i as usize] > 0  {  1.0*(ed[i as usize]-id[i as usize])/sqrt(vwgt[(i*ncon+qnum[i) as usize] as usize]+1)  } else {  ed[i as usize]-id[i as usize] });
            let rgain: real_t = (ed[i as usize] - id[i as usize]) as real_t;
            (queues[2 * qnum[i as usize] as usize + where_[i as usize] as usize])
                .insert(i as idx_t, rgain);
        }

        // for nswaps in (0)..(nvtxs) {
        let mut nswaps: idx_t = 0;
        while nswaps < nvtxs as idx_t {
            // SelectQueue(graph, ctrl.pijbm, ubfactors, queues, &from, &cnum);
            let (from, cnum) = pqueue::select_queue(graph, pijbm, &ubfactors, &queues);

            if from == -1 {
                break;
            }
            let from = from as usize;
            let cnum = cnum as usize;
            let to = (from as usize + 1) % 2;
            let Some(higain) = queues[(2 * cnum + from) as usize].pop() else {
                break;
            };
            let higain = higain as usize;
            assert!(bndptr[higain as usize] != -1);

            newcut -= ed[higain as usize] - id[higain as usize];

            // blas::iaxpy(ncon, 1, vwgt + higain * ncon, 1, pwgts + to * ncon, 1);
            blas::iaxpy(
                ncon,
                1,
                &vwgt[cntrng!(higain * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(to * ncon, ncon)],
                1,
            );
            // iaxpy(ncon, -1, vwgt + higain * ncon, 1, pwgts + from * ncon, 1);
            blas::iaxpy(
                ncon,
                -1,
                &vwgt[cntrng!(higain * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(from * ncon, ncon)],
                1,
            );
            let newbal = mcutil::ComputeLoadImbalanceDiffVec(
                graph,
                2,
                ctrl.pijbm,
                ubfactors.as_ptr(),
                newbalv.as_mut_ptr(),
            );

            if (newcut < mincut && newbal <= ffactor)
                || (newcut == mincut
                    && (newbal < minbal
                        || (newbal == minbal
                            && mcutil::BetterBalance2Way(
                                ncon as idx_t,
                                minbalv.as_ptr(),
                                newbalv.as_ptr(),
                            ) != 0)))
            {
                mincut = newcut;
                minbal = newbal;
                mincutorder = nswaps as idx_t;
                // rcopy(ncon, newbalv, minbalv);
                minbalv.copy_from_slice(&newbalv);
            } else if nswaps as idx_t - mincutorder > limit as idx_t {
                /* We hit the limit, undo last move */
                // dead assignment
                // newcut += ed[higain as usize] - id[higain as usize];

                // iaxpy(ncon, 1, vwgt + higain * ncon, 1, pwgts + from * ncon, 1);
                blas::iaxpy(
                    ncon,
                    1,
                    &vwgt[cntrng!(higain * ncon, ncon)],
                    1,
                    &mut pwgts[cntrng!(from * ncon, ncon)],
                    1,
                );
                // iaxpy(ncon, -1, vwgt + higain * ncon, 1, pwgts + to * ncon, 1);
                blas::iaxpy(
                    ncon,
                    -1,
                    &vwgt[cntrng!(higain * ncon, ncon)],
                    1,
                    &mut pwgts[cntrng!(to * ncon, ncon)],
                    1,
                );
                break;
            }

            where_[higain as usize] = to as idx_t;
            moved[higain as usize] = nswaps as idx_t;
            swaps[nswaps as usize] = higain as idx_t;

            if (ctrl.dbglvl & METIS_DBG_MOVEINFO) != 0 {
                println!(
                    "Moved{:6} from {:}({:}) Gain:{:5}, Cut:{:5}, NPwgts:",
                    higain,
                    from,
                    cnum,
                    ed[higain as usize] - id[higain as usize],
                    newcut,
                );
                for l in (0)..(ncon) {
                    println!(
                        "({:.3} {:.3})",
                        pwgts[l as usize] as real_t * invtvwgt[l as usize],
                        pwgts[(ncon + l) as usize] as real_t * invtvwgt[l as usize],
                    );
                }
                println!(
                    " {:+.3} LB: {:.3}({:+.3})",
                    minbal,
                    mcutil::ComputeLoadImbalance(graph, 2, ctrl.pijbm),
                    newbal,
                );
            }

            /**************************************************************
             * Update the id[i as usize]/ed[i as usize] values of the affected nodes
             ***************************************************************/
            std::mem::swap(&mut id[higain as usize], &mut ed[higain as usize]);
            if ed[higain as usize] == 0 && xadj[higain as usize] < xadj[(higain + 1) as usize] {
                BNDDelete!(nbnd, bndind, bndptr, higain);
            }

            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize];

                let kwgt = if to == where_[k as usize] as usize {
                    adjwgt[j as usize]
                } else {
                    -adjwgt[j as usize]
                };
                inc_dec!(id[k as usize], ed[k as usize], kwgt);

                /* Update its boundary information and queue position */
                if bndptr[k as usize] != -1 {
                    /* If k was a boundary vertex */
                    if ed[k as usize] == 0 {
                        /* Not a boundary vertex any more */
                        BNDDelete!(nbnd, bndind, bndptr, k as usize);
                        if moved[k as usize] == -1
                        /* Remove it if in the queues */
                        {
                            queues[2 * qnum[k as usize] as usize + where_[k as usize] as usize]
                                .delete(k);
                        }
                    } else {
                        /* If it has not been moved, update its position in the queue */
                        if moved[k as usize] == -1 {
                            //rgain = 1.0*(ed[k as usize]-id[k as usize])/sqrt(vwgt[(k*ncon+qnum[k) as usize] as usize]+1);
                            //rgain = (ed[k as usize]-id[k as usize] > 0 ?
                            //              1.0*(ed[k as usize]-id[k as usize])/sqrt(vwgt[(k*ncon+qnum[k) as usize] as usize]+1) : ed[k as usize]-id[k as usize]);
                            let rgain = ed[k as usize] - id[k as usize];
                            queues[2 * qnum[k as usize] as usize + where_[k as usize] as usize]
                                .update(k, rgain as real_t);
                        }
                    }
                } else {
                    if ed[k as usize] > 0 {
                        /* It will now become a boundary vertex */
                        BNDInsert!(nbnd, bndind, bndptr, k as usize);
                        if moved[k as usize] == -1 {
                            //rgain = 1.0*(ed[k as usize]-id[k as usize])/sqrt(vwgt[(k*ncon+qnum[k) as usize] as usize]+1);
                            //rgain = (ed[k as usize]-id[k as usize] > 0 ?
                            //              1.0*(ed[k as usize]-id[k as usize])/sqrt(vwgt[(k*ncon+qnum[k) as usize] as usize]+1) : ed[k as usize]-id[k as usize]);
                            let rgain = ed[k as usize] - id[k as usize];
                            queues[2 * qnum[k as usize] as usize + where_[k as usize] as usize]
                                .insert(k, rgain as real_t);
                        }
                    }
                }
            }
            nswaps += 1;
        }

        /****************************************************************
         * Roll back computations
         *****************************************************************/
        for i in (0)..(nswaps) {
            moved[swaps[i as usize] as usize] = -1; /* reset moved array */
        }
        nswaps -= 1;
        // for (nswaps--; nswaps>mincutorder; nswaps--) {
        while nswaps > mincutorder {
            let higain = swaps[nswaps as usize] as usize;

            // to = where_[higain as usize] = (where_[higain as usize] + 1) % 2;
            let to = (where_[higain as usize] as usize + 1) % 2;
            where_[higain as usize] = to as idx_t;
            std::mem::swap(&mut id[higain as usize], &mut ed[higain as usize]);
            if ed[higain as usize] == 0
                && bndptr[higain as usize] != -1
                && xadj[higain as usize] < xadj[(higain + 1) as usize]
            {
                BNDDelete!(nbnd, bndind, bndptr, higain);
            } else if ed[higain as usize] > 0 && bndptr[higain as usize] == -1 {
                BNDInsert!(nbnd, bndind, bndptr, higain);
            }

            // iaxpy(ncon, 1, vwgt + higain * ncon, 1, pwgts + to * ncon, 1);
            blas::iaxpy(
                ncon,
                1,
                &vwgt[cntrng!(higain * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(to * ncon, ncon)],
                1,
            );
            // iaxpy( ncon, -1, vwgt + higain * ncon, 1, pwgts + ((to + 1) % 2) * ncon, 1,);
            let from = (to + 1) % 2;
            blas::iaxpy(
                ncon,
                -1,
                &vwgt[cntrng!(higain * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(from * ncon, ncon)],
                1,
            );
            for j in (xadj[higain as usize])..(xadj[(higain + 1) as usize]) {
                let k = adjncy[j as usize];

                let kwgt = if to == where_[k as usize] as usize {
                    adjwgt[j as usize]
                } else {
                    -adjwgt[j as usize]
                };
                inc_dec!(id[k as usize], ed[k as usize], kwgt);

                if bndptr[k as usize] != -1 && ed[k as usize] == 0 {
                    BNDDelete!(nbnd, bndind, bndptr, k as usize);
                }
                if bndptr[k as usize] == -1 && ed[k as usize] > 0 {
                    BNDInsert!(nbnd, bndind, bndptr, k as usize);
                }
            }
            nswaps -= 1;
        }

        graph.mincut = mincut;
        graph.nbnd = nbnd as idx_t;

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_REFINE,
            Print2WayRefineStats(ctrl, graph, ntpwgts, minbal, mincutorder),
        );

        if mincutorder <= 0 || mincut == initcut {
            break;
        }
    }

    // WCOREPOP;
}

/*************************************************************************/
/* Prints statistics about the refinement */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Print2WayRefineStats(
    ctrl: *const ctrl_t,
    graph: *const graph_t,
    ntpwgts: *const real_t,
    deltabal: real_t,
    mincutorder: idx_t,
) {
    // int i;

    let graph = graph.as_ref().unwrap();
    let ctrl = ctrl.as_ref().unwrap();
    let ncon = graph.ncon as usize;
    get_graph_slices!(ctrl, graph => invtvwgt pwgts);
    mkslice!(ntpwgts, 2 * ncon);
    if mincutorder == -2 {
        println!("Parts: ");
        println!(
            "Nv-Nb[({:5} {:5})] ICut: {:6}",
            graph.nvtxs, graph.nbnd, graph.mincut,
        );
        print!(" [");
        for i in (0)..(ncon) {
            println!(
                "({:.3} {:.3} T:{:.3} {:.3})",
                pwgts[i as usize] as real_t * invtvwgt[i as usize],
                pwgts[(ncon + i) as usize] as real_t * invtvwgt[i as usize],
                ntpwgts[i as usize],
                ntpwgts[(ncon + i) as usize],
            );
        }
        print!(
            "] LB: {:.3}({:+.3})",
            mcutil::ComputeLoadImbalance(graph, 2, ctrl.pijbm),
            deltabal,
        );
    } else {
        println!(
            "\tMincut: {:6} at {:5} NBND {:6} NPwgts: [",
            graph.mincut, mincutorder, graph.nbnd,
        );
        for i in (0)..(graph.ncon) {
            println!(
                "({:.3} {:.3})",
                pwgts[i as usize] as real_t * invtvwgt[i as usize],
                pwgts[(graph.ncon + i) as usize] as real_t * invtvwgt[i as usize],
            );
        }
        println!(
            "] LB: {:.3}({:+.3})",
            mcutil::ComputeLoadImbalance(graph, 2, ctrl.pijbm),
            deltabal,
        );
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use super::*;
    use crate::dyncall::ab_test_single_eq;
    use crate::graph_gen::GraphBuilder;
    use crate::tests::ab_test_partition_test_graphs;

    #[test]
    fn ab_FM_Mc2WayCutRefine() {
        ab_test_partition_test_graphs("FM_Mc2WayCutRefine:rs", Optype::Kmetis, 20, 2, |mut g| {
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    fn ab_FM_2WayCutRefine() {
        ab_test_partition_test_graphs("FM_2WayCutRefine:rs", Optype::Kmetis, 20, 1, |mut g| {
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }
}
