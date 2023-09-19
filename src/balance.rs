/*
\file
\brief Functions for the edge-based balancing

\date Started 7/23/97
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota
\version\verbatim $Id: balance.c 10187 2011-06-13 13:46:57Z karypis $ \endverbatim
*/

use crate::*;

/*************************************************************************
* This function is the entry poidx_t of the bisection balancing algorithms.
**************************************************************************/
#[metis_func]
pub extern "C" fn Balance2Way(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    if (ComputeLoadImbalanceDiff(graph, 2, ctrl.pijbm, ctrl.ubfactors) <= 0.0) {
        return;
    }
    if (graph.ncon == 1) {
        /* return right away if the balance is OK */
        {
            get_graph_slices!(ctrl, graph => tvwgt pwgts);
            if ((*ntpwgts * tvwgt[0] as f32 - pwgts[0] as f32).abs()
                < 3 as f32 * tvwgt[0] as f32 / graph.nvtxs as f32)
            {
                return;
            }
        }

        if (graph.nbnd > 0) {
            Bnd2WayBalance(ctrl, graph, ntpwgts);
        } else {
            General2WayBalance(ctrl, graph, ntpwgts);
        }
    } else {
        McGeneral2WayBalance(ctrl, graph, ntpwgts);
    }
}

/*************************************************************************
* This function balances two partitions by moving boundary nodes
* from the domain that is overweight to the one that is underweight.
**************************************************************************/
#[metis_func]
pub extern "C" fn Bnd2WayBalance(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, pass, me, tmp;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_, *id, *ed, *bndptr, *bndind, *pwgts;
    // idx_t *moved, *perm;
    // rpq_t *queue;
    // idx_t higain, mincut, mindiff;
    // idx_t tpwgts[2];

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj vwgt adjncy adjwgt tvwgt);
    get_graph_slices_mut!(ctrl, graph => where_ id ed pwgts bndind bndptr);

    let mut moved: Vec<idx_t> = vec![0; nvtxs];
    let mut perm: Vec<idx_t> = vec![0; nvtxs];

    mkslice!(ntpwgts, graph.ncon); // maybe * 2

    /* Determine from which domain you will be moving data */
    let mut tpwgts = [0; 2];
    tpwgts[0] = (tvwgt[0] as f32 * ntpwgts[0]) as idx_t;
    let tpwgts = [tpwgts[0], tvwgt[0] - tpwgts[0]];
    let mindiff = (tpwgts[0] - pwgts[0]).abs();
    let from = (if pwgts[0] < tpwgts[0] { 1 } else { 0 }) as usize;
    let to = (from + 1) % 2;

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        println!(
            "Partitions: [{:6} {:6}] T[{:6} {:6}], Nv-Nb[{:6} {:6}]. ICut: {:6} [B]\n",
            pwgts[0], pwgts[1], tpwgts[0], tpwgts[1], graph.nvtxs, graph.nbnd, graph.mincut,
        )
    );

    // queue = rpqCreate(nvtxs);
    // TODO: do we want this an rpq or ipq?
    let mut queue = pqueue::RPQueue::new(nvtxs);

    // iset(nvtxs, -1, moved);
    moved.fill(-1);

    assert!(ComputeCut(graph, where_.as_ptr()) == graph.mincut);
    assert!(debug::CheckBnd(graph) != 0);

    /* Insert the boundary nodes of the proper partition whose size is OK in the priority queue */
    let mut nbnd = graph.nbnd as usize;
    irandArrayPermute(nbnd as idx_t, perm.as_mut_ptr(), nbnd as idx_t / 5, 1);
    for ii in (0)..(nbnd) {
        let i = perm[ii as usize] as usize;
        let ind = bndind[i] as usize;
        assert!(ed[ind] > 0 || id[ind] == 0);
        assert!(bndptr[ind] != -1);
        if where_[ind] == from as idx_t && vwgt[ind] <= mindiff {
            queue.insert(ind as idx_t, ed[ind] as f32 - id[ind] as f32);
        }
    }

    let mut mincut = graph.mincut;
    for nswaps in (0)..(nvtxs) {
        // let higain = rpqGetTop(queue)
        let higain = queue.get_top();
        if (higain == None) {
            break;
        }
        let higain = higain.unwrap() as usize;

        assert!(bndptr[higain] != -1);

        if (pwgts[to] + vwgt[higain] > tpwgts[to]) {
            break;
        }

        mincut -= (ed[higain] - id[higain]);
        inc_dec!(pwgts[to], pwgts[from], vwgt[higain]);

        where_[higain] = to as idx_t;
        moved[higain] = nswaps as idx_t;

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_MOVEINFO,
            println!(
                "Moved {:6} from {:}. [{:3} {:3}] {:5} [{:4} {:4}]\n",
                higain,
                from,
                ed[higain] - id[higain],
                vwgt[higain],
                mincut,
                pwgts[0],
                pwgts[1],
            )
        );

        /**************************************************************
         * Update the id[i]/ed[i] values of the affected nodes
         ***************************************************************/
        std::mem::swap(&mut id[higain], &mut ed[higain]);
        if (ed[higain] == 0 && xadj[higain] < xadj[higain + 1]) {
            BNDDelete!(nbnd, bndind, bndptr, higain);
        }

        for j in (xadj[higain])..(xadj[higain + 1]) {
            let j = j as usize;
            let k = adjncy[j] as usize;
            let kwgt = (if to as idx_t == where_[k] {
                adjwgt[j]
            } else {
                -adjwgt[j]
            });
            inc_dec!(id[k], ed[k], kwgt);

            /* Update its boundary information and queue position */
            if (bndptr[k] != -1) {
                /* If k was a boundary vertex */
                if (ed[k] == 0) {
                    /* Not a boundary vertex any more */
                    let mut nbnd = nbnd as usize;
                    BNDDelete!(nbnd, bndind, bndptr, k);
                    if (moved[k] == -1 && where_[k] == from as idx_t && vwgt[k] <= mindiff)
                    /* Remove it if in the queues */
                    {
                        // rpqDelete(queue, k);
                        queue.delete(k as idx_t);
                    }
                } else {
                    /* If it has not been moved, update its position in the queue */
                    if (moved[k] == -1 && where_[k] == from as idx_t && vwgt[k] <= mindiff) {
                        // rpqUpdate(queue, k, ed[k] - id[k]);
                        queue.update(k as idx_t, ed[k] as real_t - id[k] as real_t);
                    }
                }
            } else {
                if (ed[k] > 0) {
                    /* It will now become a boundary vertex */
                    BNDInsert!(nbnd, bndind, bndptr, k);
                    if (moved[k] == -1 && where_[k] == from as idx_t && vwgt[k] <= mindiff) {
                        // rpqInsert(queue, k, ed[k] - id[k]);
                        queue.insert(k as idx_t, ed[k] as real_t - id[k] as real_t);
                    }
                }
            }
        }
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        println!(
            "\tMinimum cut: {:6}, PWGTS: [{:6} {:6}], NBND: {:6}\n",
            mincut, pwgts[0], pwgts[1], nbnd,
        )
    );

    graph.mincut = mincut;
    graph.nbnd = nbnd as idx_t;

    // rpqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************
* This function balances two partitions by moving the highest gain
* (including negative gain) vertices to the other domain.
* It is used only when the unbalance is due to non contiguous
* subdomains. That is, the are no boundary vertices.
* It moves vertices from the domain that is overweight to the one that
* is underweight.
**************************************************************************/
#[metis_func]
pub extern "C" fn General2WayBalance(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, kwgt, nvtxs, nbnd, nswaps, from, to, pass, me, tmp;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_, *id, *ed, *bndptr, *bndind, *pwgts;
    // idx_t *moved, *perm;
    // rpq_t *queue;
    // idx_t higain, mincut, mindiff;
    // idx_t tpwgts[2];

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj vwgt adjncy adjwgt);
    get_graph_slices_mut!(ctrl, graph => where_ ed id pwgts bndind bndptr);

    let mut moved: Vec<idx_t> = vec![0; nvtxs];
    let mut perm: Vec<idx_t> = vec![0; nvtxs];

    /* Determine from which domain you will be moving data */
    let tpwgts: [idx_t; 2] = {
        let mut tpwgts = [0; 2];
        get_graph_slices!(graph => tvwgt);
        mkslice!(ntpwgts, 1);
        tpwgts[0] = (tvwgt[0] as real_t * ntpwgts[0]) as idx_t;
        tpwgts[1] = tvwgt[0] - tpwgts[0];
        tpwgts
    };
    let mindiff = (tpwgts[0] - pwgts[0]).abs();
    let from = (if pwgts[0] < tpwgts[0] { 1 } else { 0 }) as usize;
    let to = (from + 1) % 2;

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        println!(
            "Partitions: [{:6} {:6}] T[{:6} {:6}], Nv-Nb[{:6} {:6}]. ICut: {:6} [B]\n",
            pwgts[0], pwgts[1], tpwgts[0], tpwgts[1], graph.nvtxs, graph.nbnd, graph.mincut,
        )
    );

    // queue = rpqCreate(nvtxs);
    let mut queue = pqueue::RPQueue::new(nvtxs);

    // iset(nvtxs, -1, moved);
    moved.fill(-1);

    assert!(ComputeCut(graph, where_.as_ptr()) == graph.mincut);
    assert!(debug::CheckBnd(graph) != 0);

    /* Insert the nodes of the proper partition whose size is OK in the priority queue */
    irandArrayPermute(nvtxs as idx_t, perm.as_mut_ptr(), nvtxs as idx_t / 5, 1);
    for ii in 0..nvtxs {
        let i = perm[ii] as usize;
        if where_[i] == from as idx_t && vwgt[i] <= mindiff {
            // rpqInsert(queue, i, ed[i] - id[i]);
            queue.insert(i as idx_t, ed[i] as real_t - id[i] as real_t);
        }
    }

    let mut mincut = graph.mincut;
    let mut nbnd = graph.nbnd;
    for nswaps in (0)..(nvtxs) {
        let higain = queue.get_top();
        if higain == None || higain == Some(-1) {
            break;
        }
        let higain = higain.unwrap() as usize;

        if pwgts[to] + vwgt[higain] > tpwgts[to as usize] {
            break;
        }

        mincut -= ed[higain] - id[higain];
        inc_dec!(pwgts[to], pwgts[from], vwgt[higain]);

        where_[higain] = to as idx_t;
        moved[higain] = nswaps as idx_t;

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_MOVEINFO,
            println!(
                "Moved {:6} from {:}. [{:3} {:3}] {:5} [{:4} {:4}]\n",
                higain,
                from,
                ed[higain] - id[higain],
                vwgt[higain],
                mincut,
                pwgts[0],
                pwgts[1],
            )
        );

        /**************************************************************
         * Update the id[i]/ed[i] values of the affected nodes
         ***************************************************************/
        std::mem::swap(&mut id[higain], &mut ed[higain]);
        if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain + 1]) {
            BNDDelete!(nbnd, bndind, bndptr, higain);
        }

        if (ed[higain] > 0 && bndptr[higain] == -1) {
            BNDInsert!(nbnd, bndind, bndptr, higain);
        }

        for j in (xadj[higain])..(xadj[higain + 1]) {
            let j = j as usize;
            let k = adjncy[j] as usize;

            let kwgt = (if to == where_[k] as usize {
                adjwgt[j]
            } else {
                -adjwgt[j]
            });
            inc_dec!(id[k], ed[k], kwgt);

            /* Update the queue position */
            if (moved[k] == -1 && where_[k] == from as idx_t && vwgt[k] <= mindiff) {
                // rpqUpdate(queue, k, ed[k] - id[k]);
                queue.update(k as idx_t, (ed[k] - id[k]) as f32)
            }

            /* Update its boundary information */
            if (ed[k] == 0 && bndptr[k] != -1) {
                BNDDelete!(nbnd, bndind, bndptr, k);
            } else if (ed[k] > 0 && bndptr[k] == -1) {
                BNDInsert!(nbnd, bndind, bndptr, k);
            }
        }
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        println!(
            "\tMinimum cut: {:6}, PWGTS: [{:6} {:6}], NBND: {:6}\n",
            mincut, pwgts[0], pwgts[1], nbnd,
        )
    );

    graph.mincut = mincut;
    graph.nbnd = nbnd;

    // rpqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
#[metis_func]
pub extern "C" fn McGeneral2WayBalance(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, l, kwgt, nvtxs, ncon, nbnd, nswaps, from, to, pass,
    //       me, limit, tmp, cnum;
    // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *where_, *pwgts, *id, *ed, *bndptr, *bndind;
    // idx_t *moved, *swaps, *perm, *qnum, *qsizes;
    // idx_t higain, mincut, newcut, mincutorder;
    // real_t *invtvwgt, *minbalv, *newbalv, minbal, newbal;
    // rpq_t **queues;

    // WCOREPUSH;

    get_graph_slices!(ctrl, graph => xadj vwgt adjncy adjwgt invtvwgt);
    get_graph_slices_mut!(ctrl, graph => pwgts where_ id ed bndind bndptr);
    let ncon = graph.ncon as usize;
    let nvtxs = graph.nvtxs as usize;

    let mut moved: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut swaps = vec![0; nvtxs as usize];
    let mut perm = vec![0; nvtxs as usize];
    let mut qnum = vec![0; nvtxs as usize];
    let mut newbalv = vec![0.0; ncon as usize];
    let mut minbalv: Vec<real_t> = vec![0.0; ncon as usize];
    let mut qsizes: Vec<idx_t> = vec![0; 2 * ncon as usize];

    mkslice!(ntpwgts, ncon * 2);

    // let limit = gk_min(gk_max(0.01 * nvtxs, 15), 100);
    let limit = ((0.01 * nvtxs as real_t) as idx_t).clamp(15, 100);

    /* Initialize the queues */
    // let mut queues = vec![0; 2*ncon as usize];
    let mut queues = Vec::with_capacity(2 * ncon as usize);
    for i in (0)..(2 * ncon as usize) {
        // queues[i] = rpqCreate(nvtxs);
        queues[i] = pqueue::RPQueue::new(nvtxs);
        qsizes[i] = 0;
    }

    for i in (0)..(nvtxs) {
        qnum[i] = util::iargmax_nrm(ncon as usize, &vwgt[i * ncon as usize], invtvwgt.as_ptr());
        qsizes[2 * qnum[i] as usize + where_[i] as usize] += 1;
    }

    /* for the empty queues, move into them vertices from other queues */
    for from in 0..2 {
        for j in 0..(ncon as usize) {
            if qsizes[2 * j + from] == 0 {
                for i in (0)..(nvtxs) {
                    if (where_[i] != from as idx_t) {
                        continue;
                    }

                    let k = util::iargmax2_nrm(
                        ncon as usize,
                        &vwgt[i * ncon as usize],
                        invtvwgt.as_ptr(),
                    ) as usize;
                    if (k == j
                        && qsizes[2 * qnum[i] as usize + from] > qsizes[2 * j + from]
                        && vwgt[i * ncon as usize + qnum[i] as usize] as real_t
                            * (invtvwgt[qnum[i] as usize] as real_t)
                            < 1.3 * vwgt[i * ncon as usize + j] as real_t * invtvwgt[j] as real_t)
                    {
                        qsizes[2 * qnum[i] as usize + from] -= 1;
                        qsizes[2 * j + from] += 1;
                        qnum[i] = j as idx_t;
                    }
                }
            }
        }
    }

    let mut minbal =
        ComputeLoadImbalanceDiffVec(graph, 2, ctrl.pijbm, ctrl.ubfactors, minbalv.as_mut_ptr());
    assert!(minbal > 0.0);

    let mut newcut = graph.mincut;
    let mut mincut = graph.mincut;
    let mut mincutorder: idx_t = -1;

    if (ctrl.dbglvl & METIS_DBG_REFINE != 0) {
        println!("Parts: [");
        for l in (0)..(ncon) {
            println!(
                "({:6} {:6} {:.3} {:.3}) ",
                pwgts[l as usize],
                pwgts[(ncon + l) as usize],
                ntpwgts[l as usize],
                ntpwgts[(ncon + l) as usize],
            );
        }

        println!(
            "] Nv-Nb[{:5}, {:5}]. ICut: {:6}, LB: {:+.3} [B]\n",
            graph.nvtxs, graph.nbnd, graph.mincut, minbal,
        );
    }

    // iset(nvtxs, -1, moved);
    moved.fill(-1);

    assert!(ComputeCut(graph, where_.as_ptr()) == graph.mincut);
    assert!(debug::CheckBnd(graph) != 0);

    /* Insert all nodes in the priority queues */
    let mut nbnd = graph.nbnd;
    irandArrayPermute(nvtxs as idx_t, perm.as_mut_ptr(), nvtxs as idx_t / 10, 1);
    for ii in (0)..(nvtxs) {
        let i = perm[ii] as usize;
        // rpqInsert(queues[2 * qnum[i] + where_[i]], i, ed[i] - id[i]);
        queues[(2 * qnum[i] + where_[i]) as usize]
            .insert(i as idx_t, ed[i] as real_t - id[i] as real_t);
    }

    let mut nswaps = 0;
    for nswaps_ in (0)..(nvtxs) {
        nswaps = nswaps_;
        if (minbal <= 0.0) {
            break;
        }

        let mut from: idx_t = 0;
        let mut cnum: idx_t = 0;
        pqueue::select_queue(
            graph,
            std::slice::from_raw_parts(ctrl.pijbm, (*ctrl).nparts as usize * ncon),
            std::slice::from_raw_parts(ctrl.ubfactors, ncon),
            &mut queues,
            &mut from,
            &mut cnum,
        );
        let to = (from + 1) % 2;

        // if (from == -1 || (higain = rpqGetTop(queues[(2 * cnum + from) as usize])) == -1) {
        let higain = queues[(2 * cnum + from) as usize].get_top();
        if (from == -1 || higain == None) {
            break;
        }
        let higain = higain.unwrap() as usize;
        newcut -= (ed[higain] - id[higain]);

        blas::iaxpy(
            ncon,
            1,
            &vwgt[(higain * ncon)..],
            1,
            &mut pwgts[(to as usize * ncon)..],
            1,
        );
        blas::iaxpy(
            ncon,
            -1,
            &vwgt[(higain * ncon)..],
            1,
            &mut pwgts[(from as usize * ncon)..],
            1,
        );
        let newbal = ComputeLoadImbalanceDiffVec(graph, 2, ctrl.pijbm, ctrl.ubfactors, newbalv.as_mut_ptr());

        if (newbal < minbal
            || (newbal == minbal
                && (newcut < mincut
                    || (newcut == mincut
                        && BetterBalance2Way(
                            ncon as idx_t,
                            minbalv.as_mut_ptr(),
                            newbalv.as_mut_ptr(),
                        ) != 0))))
        {
            mincut = newcut;
            minbal = newbal;
            mincutorder = nswaps as idx_t;
            // rcopy(ncon, newbalv, minbalv);
            minbalv.copy_from_slice(&newbalv);
        } else if (nswaps as idx_t - mincutorder > limit) {
            /* We hit the limit, undo last move */
            newcut += (ed[higain] - id[higain]);
            blas::iaxpy(
                ncon,
                1,
                &vwgt[(higain * ncon)..],
                1,
                &mut pwgts[(from as usize * ncon)..],
                1,
            );
            blas::iaxpy(
                ncon,
                -1,
                &vwgt[(higain * ncon)..],
                1,
                &mut pwgts[(to as usize * ncon)..],
                1,
            );
            break;
        }

        where_[higain] = to;
        moved[higain] = nswaps as idx_t;
        swaps[nswaps] = higain;

        if (ctrl.dbglvl & METIS_DBG_MOVEINFO != 0) {
            println!(
                "Moved {:6} from {:}({:}). Gain: {:5}, Cut: {:5}, NPwgts: ",
                higain,
                from,
                cnum,
                ed[higain] - id[higain],
                newcut,
            );
            for l in (0)..(ncon) {
                println!("({:6}, {:6}) ", pwgts[l], pwgts[ncon + l]);
            }

            println!(", {:+.3} LB: {:+.3}\n", minbal, newbal);
        }

        /**************************************************************
         * Update the id[i]/ed[i] values of the affected nodes
         ***************************************************************/
        std::mem::swap(&mut id[higain], &mut ed[higain]);
        if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain + 1]) {
            BNDDelete!(nbnd, bndind, bndptr, higain);
        }

        if (ed[higain] > 0 && bndptr[higain] == -1) {
            BNDInsert!(nbnd, bndind, bndptr, higain);
        }

        for j in xadj[higain]..xadj[higain + 1] {
            let j = j as usize;
            let k = adjncy[j] as usize;

            let kwgt = (if to == where_[k] {
                adjwgt[j]
            } else {
                -adjwgt[j]
            });
            inc_dec!(id[k], ed[k], kwgt);

            /* Update the queue position */
            if (moved[k] == -1) {
                queues[2 * qnum[k] as usize + where_[k] as usize]
                    .update(k as idx_t, (ed[k] - id[k]) as real_t);
            }

            /* Update its boundary information */
            if (ed[k] == 0 && bndptr[k] != -1) {
                BNDDelete!(nbnd, bndind, bndptr, k);
            } else if (ed[k] > 0 && bndptr[k] == -1) {
                BNDInsert!(nbnd, bndind, bndptr, k);
            }
        }
    }

    /****************************************************************
     * Roll back computations
     *****************************************************************/
    // for (nswaps--; nswaps>mincutorder; nswaps--) {
    nswaps -= 1;
    while nswaps as idx_t > mincutorder {
        nswaps -= 1;
        let higain = swaps[nswaps];

        // let to = where_[higain] = (where_[higain] + 1) % 2;
        where_[higain] = (where_[higain] + 1) % 2;
        let to = where_[higain];
        std::mem::swap(&mut id[higain], &mut ed[higain]);
        if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain + 1]) {
            BNDDelete!(nbnd, bndind, bndptr, higain);
        } else if (ed[higain] > 0 && bndptr[higain] == -1) {
            BNDInsert!(nbnd, bndind, bndptr, higain);
        }

        blas::iaxpy(
            ncon,
            1,
            &vwgt[(higain * ncon)..],
            1,
            &mut pwgts[(to as usize * ncon)..],
            1,
        );
        blas::iaxpy(
            ncon,
            -1,
            &vwgt[(higain * ncon)..],
            1,
            &mut pwgts[(((to as usize + 1) % 2) * ncon)..],
            1,
        );
        for j in (xadj[higain])..(xadj[higain + 1]) {
            let j = j as usize;
            let k = adjncy[j] as usize;

            let kwgt = (if to == where_[k] {
                adjwgt[j]
            } else {
                -adjwgt[j]
            });
            inc_dec!(id[k], ed[k], kwgt);

            if (bndptr[k] != -1 && ed[k] == 0) {
                BNDDelete!(nbnd, bndind, bndptr, k);
            }

            if (bndptr[k] == -1 && ed[k] > 0) {
                BNDInsert!(nbnd, bndind, bndptr, k);
            }
        }
    }

    if (ctrl.dbglvl & METIS_DBG_REFINE != 0) {
        println!(
            "\tMincut: {:6} at {:5}, NBND: {:6}, NPwgts: [",
            mincut, mincutorder, nbnd,
        );
        for l in (0)..(ncon) {
            println!("({:6}, {:6}) ", pwgts[l], pwgts[ncon + l]);
        }

        println!("], LB: {:.3}\n", ComputeLoadImbalance(graph, 2, ctrl.pijbm));
    }

    graph.mincut = mincut;
    graph.nbnd = nbnd;

    // for i in (0)..(2 * ncon) {
    //     rpqDestroy(queues[i]);
    // }

    // WCOREPOP;
}
