/*
\file
\brief Routines for k-way refinement

\date Started 7/28/97
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: kwayfm.c 17513 2014-08-05 16:20:50Z dominique $
*/

use crate::*;

/*************************************************************************/
/* Top-level routine for k-way partitioning refinement. This routine just
calls the appropriate refinement routine based on the objectives and
constraints. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Greedy_KWayOptimize(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    niter: idx_t,
    ffactor: real_t,
    omode: idx_t,
) {
    match (ctrl.objtype) {
        METIS_OBJTYPE_CUT => {
            if (graph.ncon == 1) {
                Greedy_KWayCutOptimize(ctrl, graph, niter, ffactor, omode);
            } else {
                Greedy_McKWayCutOptimize(ctrl, graph, niter, ffactor, omode);
            }
        }

        METIS_OBJTYPE_VOL => {
            if (graph.ncon == 1) {
                Greedy_KWayVolOptimize(ctrl, graph, niter, ffactor, omode);
            } else {
                Greedy_McKWayVolOptimize(ctrl, graph, niter, ffactor, omode);
            }
        }

        _ => gk_errexit(SIGERR, "Unknown objtype of %d\n", ctrl.objtype),
    }
}

/*************************************************************************/
/* K-way partitioning optimization in which the vertices are visited in
    decreasing ed/sqrt(nnbrs)-id order. Note this is just an
    approximation, as the ed is often split across different subdomains
    and the sqrt(nnbrs) is just a crude approximation.

  \param graph is the graph that is being refined.
  \param niter is the number of refinement iterations.
  \param ffactor is the \em fudge-factor for allowing positive gain moves
         to violate the max-pwgt constraint.
  \param omode is the type of optimization that will performed among
         OMODE_REFINE and OMODE_BALANCE

*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Greedy_KWayCutOptimize(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    niter: idx_t,
    ffactor: real_t,
    omode: idx_t,
) {
    /* Common variables to all types of kway-refinement/balancing routines */
    // idx_t i, ii, iii, j, k, l, pass, nvtxs, nparts, gain;
    // idx_t from, me, to, oldcut, vwgt;
    // idx_t *xadj, *adjncy, *adjwgt;
    // idx_t *where_, *pwgts, *perm, *bndptr, *bndind, *minpwgts, *maxpwgts;
    // idx_t nmoved, nupd, *vstatus, *updptr, *updind;
    // idx_t maxndoms, *safetos=std::ptr::null_mut(), *nads=std::ptr::null_mut(), *doms=std::ptr::null_mut(), **adids=std::ptr::null_mut(), **adwgts=std::ptr::null_mut();
    // idx_t *bfslvl=std::ptr::null_mut(), *bfsind=std::ptr::null_mut(), *bfsmrk=std::ptr::null_mut();
    // idx_t bndtype = (omode == OMODE_REFINE ? BNDTYPE_REFINE : BNDTYPE_BALANCE);
    // real_t *tpwgts, ubfactor;

    /* Edgecut-specific/different variables */
    // idx_t nbnd, oldnnbrs;
    // rpq_t *queue;
    // real_t rgain;
    // ckrinfo_t *myrinfo;
    // cnbr_t *mynbrs;

    ffactor = 0.0;
    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;

    bndind = graph.bndind;
    bndptr = graph.bndptr;

    where_ = graph.where_;
    pwgts = graph.pwgts;

    nparts = ctrl.nparts;
    tpwgts = ctrl.tpwgts;

    /* Setup the weight intervals of the various subdomains */
    minpwgts = vec![0; nparts as usize];
    maxpwgts = vec![0; nparts as usize];

    if (omode == OMODE_BALANCE) {
        ubfactor = ctrl.ubfactors[0 as usize];
    } else {
        ubfactor = gk_max(
            ctrl.ubfactors[0 as usize],
            ComputeLoadImbalance(graph, nparts, ctrl.pijbm),
        );
    }

    for i in (0)..(nparts) {
        maxpwgts[i as usize] = tpwgts[i as usize] * graph.tvwgt[0 as usize] * ubfactor;
        minpwgts[i as usize] = tpwgts[i as usize] * graph.tvwgt[0 as usize] * (1.0 / ubfactor);
    }

    perm = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    safetos = vec![2; nparts as usize];

    if (ctrl.minconn) {
        ComputeSubDomainGraph(ctrl, graph);

        nads = ctrl.nads;
        adids = ctrl.adids;
        adwgts = ctrl.adwgts;
        doms = iset(nparts, 0, ctrl.pvec1);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    vstatus = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    updptr = vec![-1; nvtxs as usize];
    updind = vec![0; nvtxs as usize];

    if (ctrl.contig) {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    }

    if (ctrl.dbglvl & METIS_DBG_REFINE) {
        print!("{}: [({:6} {:6}) as usize]-[({:6} {:6}) as usize], Bal: {:5.3}," 
            " Nv-Nb[({:6} {:6}) as usize], Cut: %6"PRIDX,
            (omode == OMODE_REFINE ? "GRC" : "GBC"),
            pwgts[(iargmin(nparts, pwgts,1)) as usize], imax(nparts, pwgts,1), minpwgts[0 as usize], maxpwgts[0 as usize], 
            ComputeLoadImbalance(graph, nparts, ctrl.pijbm), 
            graph.nvtxs, graph.nbnd, graph.mincut);
        if (ctrl.minconn) {
            print!(
                ", Doms: [({:3} {:4}) as usize]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    queue = rpqCreate(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        assert!(ComputeCut(graph, where_) == graph.mincut);
        if (omode == OMODE_REFINE) {
            assert!(CheckBnd2(graph));
        }

        if (omode == OMODE_BALANCE) {
            /* Check to see if things are out of balance, given the tolerance */
            for i in (0)..(nparts) {
                if (pwgts[i as usize] > maxpwgts[i as usize]
                    || pwgts[i as usize] < minpwgts[i as usize])
                {
                    break;
                }
            }
            if (i == nparts)
            /* Things are balanced. Return right away */
            {
                break;
            }
        }

        oldcut = graph.mincut;
        nbnd = graph.nbnd;
        nupd = 0;

        if (ctrl.minconn) {
            maxndoms = imax(nparts, nads, 1);
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(nbnd, perm, nbnd / 4, 1);
        for ii in (0)..(nbnd) {
            i = bndind[perm[ii as usize] as usize];
            rgain = (if graph.ckrinfo[i as usize].nnbrs > 0 {
                1.0 * graph.ckrinfo[i as usize].ed / sqrt(graph.ckrinfo[i as usize].nnbrs)
            } else {
                0.0
            }) - graph.ckrinfo[i as usize].id;
            rpqInsert(queue, i, rgain);
            vstatus[i as usize] = VPQSTATUS_PRESENT;
            ListInsert(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        nmoved = 0;
        for iii in (0).. {
            if ((i = rpqGetTop(queue)) == -1) {
                break;
            }
            vstatus[i as usize] = VPQSTATUS_EXTRACTED;

            myrinfo = graph.ckrinfo + i;
            mynbrs = ctrl.cnbrpool + myrinfo.inbr;

            from = where_[i as usize];
            vwgt = graph.vwgt[i as usize];

            // #ifdef XXX
            //       /* Prevent moves that make 'from' domain underbalanced */
            //       if (omode == OMODE_REFINE) {
            //         if (myrinfo.id > 0 && pwgts[from as usize]-vwgt < minpwgts[from as usize])  {
            //           continue;
            // }
            //       }
            //       else { /* OMODE_BALANCE */
            //         if (pwgts[from as usize]-vwgt < minpwgts[from as usize])  {
            //           continue;
            // }
            //       }
            // #endif

            if (ctrl.contig && IsArticulationNode(i, xadj, adjncy, where_, bfslvl, bfsind, bfsmrk))
            {
                continue;
            }

            if (ctrl.minconn) {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            /* Find the most promising subdomain to move to */
            if (omode == OMODE_REFINE) {
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    if (((mynbrs[k as usize].ed > myrinfo.id)
                        && ((pwgts[from as usize] - vwgt >= minpwgts[from as usize])
                            || (tpwgts[from as usize] * pwgts[to as usize]
                                < tpwgts[to as usize] * (pwgts[from as usize] - vwgt)))
                        && ((pwgts[to as usize] + vwgt <= maxpwgts[to as usize])
                            || (tpwgts[from as usize] * pwgts[to as usize]
                                < tpwgts[to as usize] * (pwgts[from as usize] - vwgt))))
                        || ((mynbrs[k as usize].ed == myrinfo.id)
                            && (tpwgts[from as usize] * pwgts[to as usize]
                                < tpwgts[to as usize] * (pwgts[from as usize] - vwgt))))
                    {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    if (((mynbrs[j as usize].ed > mynbrs[k as usize].ed)
                        && ((pwgts[from as usize] - vwgt >= minpwgts[from as usize])
                            || (tpwgts[from as usize] * pwgts[to as usize]
                                < tpwgts[to as usize] * (pwgts[from as usize] - vwgt)))
                        && ((pwgts[to as usize] + vwgt <= maxpwgts[to as usize])
                            || (tpwgts[from as usize] * pwgts[to as usize]
                                < tpwgts[to as usize] * (pwgts[from as usize] - vwgt))))
                        || ((mynbrs[j as usize].ed == mynbrs[k as usize].ed)
                            && (tpwgts[mynbrs[k as usize].pid] * pwgts[to as usize]
                                < tpwgts[to as usize] * pwgts[mynbrs[k as usize].pid])))
                    {
                        k = j;
                    }
                }

                to = mynbrs[k as usize].pid;

                gain = mynbrs[k as usize].ed - myrinfo.id;
                /*
                if (!(gain > 0
                      || (gain == 0
                          && (pwgts[from as usize] >= maxpwgts[from as usize]
                              || tpwgts[to as usize]*pwgts[from as usize] > tpwgts[from as usize]*(pwgts[to as usize]+vwgt)
                              || (iii%2 == 0 && safetos[to as usize] == 2)
                             )
                         )
                     )
                   )
                  continue;
                */
            } else {
                /* OMODE_BALANCE */
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    /* the correctness of the following test follows from the correctness
                    of the similar test in the subsequent loop */
                    if (from >= nparts
                        || tpwgts[from as usize] * pwgts[to as usize]
                            < tpwgts[to as usize] * (pwgts[from as usize] - vwgt))
                    {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    if (tpwgts[mynbrs[k as usize].pid] * pwgts[to as usize]
                        < tpwgts[to as usize] * pwgts[mynbrs[k as usize].pid])
                    {
                        k = j;
                    }
                }

                to = mynbrs[k as usize].pid;

                //if (pwgts[from as usize] < maxpwgts[from as usize] && pwgts[to as usize] > minpwgts[to as usize] &&
                //    mynbrs[k as usize].ed-myrinfo.id < 0)
                //  continue;
            }

            /*=====================================================================
             * If we got here, we can now move the vertex from 'from' to 'to'
             *======================================================================*/
            graph.mincut -= mynbrs[k as usize].ed - myrinfo.id;
            nmoved += 1;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "\t\tMoving {:6} from {:3}/{:} to {:3}/{:} [({:6} {:6}) as usize]. Gain: {:4}. Cut: {:6}\n",
                    i,
                    from,
                    safetos[from as usize],
                    to,
                    safetos[to as usize],
                    pwgts[from as usize],
                    pwgts[to as usize],
                    mynbrs[k as usize].ed - myrinfo.id,
                    graph.mincut
                )
            );

            /* Update the subdomain connectivity information */
            if (ctrl.minconn) {
                /* take care of i's move itself */
                UpdateEdgeSubDomainGraph(
                    ctrl,
                    from,
                    to,
                    myrinfo.id - mynbrs[k as usize].ed,
                    &maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    me = where_[adjncy[j as usize] as usize];
                    if (me != from && me != to) {
                        UpdateEdgeSubDomainGraph(ctrl, from, me, -adjwgt[j as usize], &maxndoms);
                        UpdateEdgeSubDomainGraph(ctrl, to, me, adjwgt[j as usize], &maxndoms);
                    }
                }
            }

            /* Update ID/ED and BND related information for the moved vertex */
            INC_DEC(pwgts[to as usize], pwgts[from as usize], vwgt);
            UpdateMovedVertexInfoAndBND(
                i, from, k, to, myrinfo, mynbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                ii = adjncy[j as usize];
                me = where_[ii as usize];
                myrinfo = graph.ckrinfo + ii;

                oldnnbrs = myrinfo.nnbrs;

                UpdateAdjacentVertexInfoAndBND(
                    ctrl,
                    ii,
                    xadj[(ii + 1) as usize] - xadj[ii as usize],
                    me,
                    from,
                    to,
                    myrinfo,
                    adjwgt[j as usize],
                    nbnd,
                    bndptr,
                    bndind,
                    bndtype,
                );

                UpdateQueueInfo(
                    queue, vstatus, ii, me, from, to, myrinfo, oldnnbrs, nupd, updptr, updind,
                    bndtype,
                );

                assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }
        }

        graph.nbnd = nbnd;

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            assert!(updptr[updind[i as usize] as usize] != -1);
            assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if (ctrl.dbglvl & METIS_DBG_REFINE) {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}."
              " Nmoves: {:5}, Cut: {:6}, Vol: %6"PRIDX,
              pwgts[(iargmin(nparts, pwgts,1)) as usize], imax(nparts, pwgts,1),
              ComputeLoadImbalance(graph, nparts, ctrl.pijbm), 
              graph.nbnd, nmoved, graph.mincut, ComputeVolume(graph, where_));
            if (ctrl.minconn) {
                print!(
                    ", Doms: [({:3} {:4}) as usize]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if (nmoved == 0 || (omode == OMODE_REFINE && graph.mincut == oldcut)) {
            break;
        }
    }

    rpqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************/
/* K-way refinement that minimizes the communication volume. This is a
    greedy routine and the vertices are visited in decreasing gv order.

  \param graph is the graph that is being refined.
  \param niter is the number of refinement iterations.
  \param ffactor is the \em fudge-factor for allowing positive gain moves
         to violate the max-pwgt constraint.

*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Greedy_KWayVolOptimize(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    niter: idx_t,
    ffactor: real_t,
    omode: idx_t,
) {
    /* Common variables to all types of kway-refinement/balancing routines */
    // idx_t i, ii, iii, j, k, l, pass, nvtxs, nparts, gain;
    // idx_t from, me, to, oldcut, vwgt;
    // idx_t *xadj, *adjncy;
    // idx_t *where_, *pwgts, *perm, *bndptr, *bndind, *minpwgts, *maxpwgts;
    // idx_t nmoved, nupd, *vstatus, *updptr, *updind;
    // idx_t maxndoms, *safetos=std::ptr::null_mut(), *nads=std::ptr::null_mut(), *doms=std::ptr::null_mut(), **adids=std::ptr::null_mut(), **adwgts=std::ptr::null_mut();
    // idx_t *bfslvl=std::ptr::null_mut(), *bfsind=std::ptr::null_mut(), *bfsmrk=std::ptr::null_mut();
    // idx_t bndtype = (omode == OMODE_REFINE ? BNDTYPE_REFINE : BNDTYPE_BALANCE);
    // real_t *tpwgts;

    /* Volume-specific/different variables */
    // ipq_t *queue;
    // idx_t oldvol, xgain;
    // idx_t *vmarker, *pmarker, *modind;
    // vkrinfo_t *myrinfo;
    // vnbr_t *mynbrs;

    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    bndptr = graph.bndptr;
    bndind = graph.bndind;
    where_ = graph.where_;
    pwgts = graph.pwgts;

    nparts = ctrl.nparts;
    tpwgts = ctrl.tpwgts;

    /* Setup the weight intervals of the various subdomains */
    minpwgts = vec![0; nparts as usize];
    maxpwgts = vec![0; nparts as usize];

    for i in (0)..(nparts) {
        maxpwgts[i as usize] =
            ctrl.tpwgts[i as usize] * graph.tvwgt[0 as usize] * ctrl.ubfactors[0 as usize];
        minpwgts[i as usize] =
            ctrl.tpwgts[i as usize] * graph.tvwgt[0 as usize] * (1.0 / ctrl.ubfactors[0 as usize]);
    }

    perm = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    safetos = vec![2; nparts as usize];

    if (ctrl.minconn) {
        ComputeSubDomainGraph(ctrl, graph);

        nads = ctrl.nads;
        adids = ctrl.adids;
        adwgts = ctrl.adwgts;
        doms = iset(nparts, 0, ctrl.pvec1);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    vstatus = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    updptr = vec![-1; nvtxs as usize];
    updind = vec![0; nvtxs as usize];

    if (ctrl.contig) {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    }

    /* Vol-refinement specific working arrays */
    modind = vec![0; nvtxs as usize];
    vmarker = vec![0; nvtxs as usize];
    pmarker = vec![-1; nparts as usize];

    if (ctrl.dbglvl & METIS_DBG_REFINE) {
        print!("{}: [({:6} {:6}) as usize]-[({:6} {:6}) as usize], Bal: %5.3"PRREAL
         ", Nv-Nb[({:6} {:6}) as usize], Cut: {:5}, Vol: %5"PRIDX,
         (omode == OMODE_REFINE ? "GRV" : "GBV"),
         pwgts[(iargmin(nparts, pwgts,1)) as usize], imax(nparts, pwgts,1), minpwgts[0 as usize], maxpwgts[0 as usize], 
         ComputeLoadImbalance(graph, nparts, ctrl.pijbm), 
         graph.nvtxs, graph.nbnd, graph.mincut, graph.minvol);
        if (ctrl.minconn) {
            print!(
                ", Doms: [({:3} {:4}) as usize]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    queue = ipqCreate(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        assert!(ComputeVolume(graph, where_) == graph.minvol);

        if (omode == OMODE_BALANCE) {
            /* Check to see if things are out of balance, given the tolerance */
            for i in (0)..(nparts) {
                if (pwgts[i as usize] > maxpwgts[i as usize]) {
                    break;
                }
            }
            if (i == nparts)
            /* Things are balanced. Return right away */
            {
                break;
            }
        }

        oldcut = graph.mincut;
        oldvol = graph.minvol;
        nupd = 0;

        if (ctrl.minconn) {
            maxndoms = imax(nparts, nads, 1);
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(graph.nbnd, perm, graph.nbnd / 4, 1);
        for ii in (0)..(graph.nbnd) {
            i = bndind[perm[ii as usize] as usize];
            ipqInsert(queue, i, graph.vkrinfo[i as usize].gv);
            vstatus[i as usize] = VPQSTATUS_PRESENT;
            ListInsert(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        nmoved = 0;
        for iii in (0).. {
            if ((i = ipqGetTop(queue)) == -1) {
                break;
            }
            vstatus[i as usize] = VPQSTATUS_EXTRACTED;

            myrinfo = graph.vkrinfo + i;
            mynbrs = ctrl.vnbrpool + myrinfo.inbr;

            from = where_[i as usize];
            vwgt = graph.vwgt[i as usize];

            /* Prevent moves that make 'from' domain underbalanced */
            if (omode == OMODE_REFINE) {
                if (myrinfo.nid > 0 && pwgts[from as usize] - vwgt < minpwgts[from as usize]) {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                if (pwgts[from as usize] - vwgt < minpwgts[from as usize]) {
                    continue;
                }
            }

            if (ctrl.contig && IsArticulationNode(i, xadj, adjncy, where_, bfslvl, bfsind, bfsmrk))
            {
                continue;
            }

            if (ctrl.minconn) {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            xgain = (if myrinfo.nid == 0 && myrinfo.ned > 0 {
                graph.vsize[i as usize]
            } else {
                0
            });

            /* Find the most promising subdomain to move to */
            if (omode == OMODE_REFINE) {
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    gain = mynbrs[k as usize].gv + xgain;
                    if (gain >= 0
                        && pwgts[to as usize] + vwgt <= maxpwgts[to as usize] + ffactor * gain)
                    {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    gain = mynbrs[j as usize].gv + xgain;
                    if ((mynbrs[j as usize].gv > mynbrs[k as usize].gv
                        && pwgts[to as usize] + vwgt <= maxpwgts[to as usize] + ffactor * gain)
                        || (mynbrs[j as usize].gv == mynbrs[k as usize].gv
                            && mynbrs[j as usize].ned > mynbrs[k as usize].ned
                            && pwgts[to as usize] + vwgt <= maxpwgts[to as usize])
                        || (mynbrs[j as usize].gv == mynbrs[k as usize].gv
                            && mynbrs[j as usize].ned == mynbrs[k as usize].ned
                            && tpwgts[mynbrs[k as usize].pid] * pwgts[to as usize]
                                < tpwgts[to as usize] * pwgts[mynbrs[k as usize].pid]))
                    {
                        k = j;
                    }
                }
                to = mynbrs[k as usize].pid;

                assert!(xgain + mynbrs[k as usize].gv >= 0);

                j = 0;
                if (xgain + mynbrs[k as usize].gv > 0 || mynbrs[k as usize].ned - myrinfo.nid > 0) {
                    j = 1;
                } else if (mynbrs[k as usize].ned - myrinfo.nid == 0) {
                    if ((iii % 2 == 0 && safetos[to as usize] == 2)
                        || pwgts[from as usize] >= maxpwgts[from as usize]
                        || tpwgts[from as usize] * (pwgts[to as usize] + vwgt)
                            < tpwgts[to as usize] * pwgts[from as usize])
                    {
                        j = 1;
                    }
                }
                if (j == 0) {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    if (pwgts[to as usize] + vwgt <= maxpwgts[to as usize]
                        || tpwgts[from as usize] * (pwgts[to as usize] + vwgt)
                            <= tpwgts[to as usize] * pwgts[from as usize])
                    {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    if (tpwgts[mynbrs[k as usize].pid] * pwgts[to as usize]
                        < tpwgts[to as usize] * pwgts[mynbrs[k as usize].pid])
                    {
                        k = j;
                    }
                }
                to = mynbrs[k as usize].pid;

                if (pwgts[from as usize] < maxpwgts[from as usize]
                    && pwgts[to as usize] > minpwgts[to as usize]
                    && (xgain + mynbrs[k as usize].gv < 0
                        || (xgain + mynbrs[k as usize].gv == 0
                            && mynbrs[k as usize].ned - myrinfo.nid < 0)))
                {
                    continue;
                }
            }

            /*=====================================================================
             * If we got here, we can now move the vertex from 'from' to 'to'
             *======================================================================*/
            INC_DEC(pwgts[to as usize], pwgts[from as usize], vwgt);
            graph.mincut -= mynbrs[k as usize].ned - myrinfo.nid;
            graph.minvol -= (xgain + mynbrs[k as usize].gv);
            where_[i as usize] = to;
            nmoved += 1;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!("\t\tMoving {:6} from {:3} to {:3}. "
                 "Gain: [({:4} {:4}) as usize]. Cut: {:6}, Vol: {:6}\n", 
              i, from, to, xgain+mynbrs[k as usize].gv, mynbrs[k as usize].ned-myrinfo.nid, 
              graph.mincut, graph.minvol)
            );

            /* Update the subdomain connectivity information */
            if (ctrl.minconn) {
                /* take care of i's move itself */
                UpdateEdgeSubDomainGraph(
                    ctrl,
                    from,
                    to,
                    myrinfo.nid - mynbrs[k as usize].ned,
                    &maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    me = where_[adjncy[j as usize] as usize];
                    if (me != from && me != to) {
                        UpdateEdgeSubDomainGraph(ctrl, from, me, -1, &maxndoms);
                        UpdateEdgeSubDomainGraph(ctrl, to, me, 1, &maxndoms);
                    }
                }
            }

            /* Update the id/ed/gains/bnd/queue of potentially affected nodes */
            KWayVolUpdate(
                ctrl, graph, i, from, to, queue, vstatus, &nupd, updptr, updind, bndtype, vmarker,
                pmarker, modind,
            );

            /*CheckKWayVolPartitionParams(ctrl, graph); */
        }

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            assert!(updptr[updind[i as usize] as usize] != -1);
            assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if (ctrl.dbglvl & METIS_DBG_REFINE) {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}."
              " Nmoves: {:5}, Cut: {:6}, Vol: %6"PRIDX,
              pwgts[(iargmin(nparts, pwgts,1)) as usize], imax(nparts, pwgts,1),
              ComputeLoadImbalance(graph, nparts, ctrl.pijbm), 
              graph.nbnd, nmoved, graph.mincut, graph.minvol);
            if (ctrl.minconn) {
                print!(
                    ", Doms: [({:3} {:4}) as usize]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if (nmoved == 0
            || (omode == OMODE_REFINE && graph.minvol == oldvol && graph.mincut == oldcut))
        {
            break;
        }
    }

    ipqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************/
/* K-way partitioning optimization in which the vertices are visited in
    decreasing ed/sqrt(nnbrs)-id order. Note this is just an
    approximation, as the ed is often split across different subdomains
    and the sqrt(nnbrs) is just a crude approximation.

  \param graph is the graph that is being refined.
  \param niter is the number of refinement iterations.
  \param ffactor is the \em fudge-factor for allowing positive gain moves
         to violate the max-pwgt constraint.
  \param omode is the type of optimization that will performed among
         OMODE_REFINE and OMODE_BALANCE


*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Greedy_McKWayCutOptimize(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    niter: idx_t,
    ffactor: real_t,
    omode: idx_t,
) {
    /* Common variables to all types of kway-refinement/balancing routines */
    // idx_t i, ii, iii, j, k, l, pass, nvtxs, ncon, nparts, gain;
    // idx_t from, me, to, cto, oldcut;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt;
    // idx_t *where_, *pwgts, *perm, *bndptr, *bndind, *minpwgts, *maxpwgts;
    // idx_t nmoved, nupd, *vstatus, *updptr, *updind;
    // idx_t maxndoms, *safetos=std::ptr::null_mut(), *nads=std::ptr::null_mut(), *doms=std::ptr::null_mut(), **adids=std::ptr::null_mut(), **adwgts=std::ptr::null_mut();
    // idx_t *bfslvl=std::ptr::null_mut(), *bfsind=std::ptr::null_mut(), *bfsmrk=std::ptr::null_mut();
    // idx_t bndtype = (omode == OMODE_REFINE ? BNDTYPE_REFINE : BNDTYPE_BALANCE);
    // real_t *ubfactors, *pijbm;
    // real_t origbal;

    /* Edgecut-specific/different variables */
    // idx_t nbnd, oldnnbrs;
    // rpq_t *queue;
    // real_t rgain;
    // ckrinfo_t *myrinfo;
    // cnbr_t *mynbrs;

    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;

    bndind = graph.bndind;
    bndptr = graph.bndptr;

    where_ = graph.where_;
    pwgts = graph.pwgts;

    nparts = ctrl.nparts;
    pijbm = ctrl.pijbm;

    /* Determine the ubfactors. The method used is different based on omode.
    When OMODE_BALANCE, the ubfactors are those supplied by the user.
    When OMODE_REFINE, the ubfactors are the max of the current partition
    and the user-specified ones. */
    ubfactors = vec![0.0; ncon as usize];
    ComputeLoadImbalanceVec(graph, nparts, pijbm, ubfactors);
    origbal = rvecmaxdiff(ncon, ubfactors, ctrl.ubfactors);
    if (omode == OMODE_BALANCE) {
        rcopy(ncon, ctrl.ubfactors, ubfactors);
    } else {
        for i in (0)..(ncon) {
            ubfactors[i as usize] = (if ubfactors[i as usize] > ctrl.ubfactors[i as usize] {
                ubfactors[i as usize]
            } else {
                ctrl.ubfactors[i as usize]
            });
        }
    }

    /* Setup the weight intervals of the various subdomains */
    minpwgts = vec![0; nparts * ncon as usize];
    maxpwgts = vec![0; nparts * ncon as usize];

    for i in (0)..(nparts) {
        for j in (0)..(ncon) {
            maxpwgts[(i * ncon + j) as usize] = ctrl.tpwgts[(i * ncon + j) as usize]
                * graph.tvwgt[j as usize]
                * ubfactors[j as usize];
            /*minpwgts[(i*ncon+j) as usize]  = ctrl.tpwgts[(i*ncon+j) as usize]*graph.tvwgt[j as usize]*(.9/ubfactors[j as usize]);*/
            minpwgts[(i * ncon + j) as usize] =
                ctrl.tpwgts[(i * ncon + j) as usize] * graph.tvwgt[j as usize] * 0.2;
        }
    }

    perm = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    safetos = vec![2; nparts as usize];

    if (ctrl.minconn) {
        ComputeSubDomainGraph(ctrl, graph);

        nads = ctrl.nads;
        adids = ctrl.adids;
        adwgts = ctrl.adwgts;
        doms = iset(nparts, 0, ctrl.pvec1);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    vstatus = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    updptr = vec![-1; nvtxs as usize];
    updind = vec![0; nvtxs as usize];

    if (ctrl.contig) {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    }

    if (ctrl.dbglvl & METIS_DBG_REFINE) {
        print!("{}: [({:6} {:6} {:6}) as usize], Bal: {:5.3}({:.3})," 
            " Nv-Nb[({:6} {:6}) as usize], Cut: {:6}, ({:})",
            (omode == OMODE_REFINE ? "GRC" : "GBC"),
            imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), imax(nparts*ncon, maxpwgts,1),
            ComputeLoadImbalance(graph, nparts, pijbm), origbal,
            graph.nvtxs, graph.nbnd, graph.mincut, niter);
        if (ctrl.minconn) {
            print!(
                ", Doms: [({:3} {:4}) as usize]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    queue = rpqCreate(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        assert!(ComputeCut(graph, where_) == graph.mincut);
        if (omode == OMODE_REFINE) {
            assert!(CheckBnd2(graph));
        }

        /* In balancing mode, exit as soon as balance is reached */
        if (omode == OMODE_BALANCE && IsBalanced(ctrl, graph, 0)) {
            break;
        }

        oldcut = graph.mincut;
        nbnd = graph.nbnd;
        nupd = 0;

        if (ctrl.minconn) {
            maxndoms = imax(nparts, nads, 1);
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(nbnd, perm, nbnd / 4, 1);
        for ii in (0)..(nbnd) {
            i = bndind[perm[ii as usize] as usize];
            rgain = (if graph.ckrinfo[i as usize].nnbrs > 0 {
                1.0 * graph.ckrinfo[i as usize].ed / sqrt(graph.ckrinfo[i as usize].nnbrs)
            } else {
                0.0
            }) - graph.ckrinfo[i as usize].id;
            rpqInsert(queue, i, rgain);
            vstatus[i as usize] = VPQSTATUS_PRESENT;
            ListInsert(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        nmoved = 0;
        for iii in (0).. {
            if ((i = rpqGetTop(queue)) == -1) {
                break;
            }
            vstatus[i as usize] = VPQSTATUS_EXTRACTED;

            myrinfo = graph.ckrinfo + i;
            mynbrs = ctrl.cnbrpool + myrinfo.inbr;

            from = where_[i as usize];

            /* Prevent moves that make 'from' domain underbalanced */
            if (omode == OMODE_REFINE) {
                if (myrinfo.id > 0
                    && !ivecaxpygez(
                        ncon,
                        -1,
                        vwgt + i * ncon,
                        pwgts + from * ncon,
                        minpwgts + from * ncon,
                    ))
                {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                if (!ivecaxpygez(
                    ncon,
                    -1,
                    vwgt + i * ncon,
                    pwgts + from * ncon,
                    minpwgts + from * ncon,
                )) {
                    continue;
                }
            }

            if (ctrl.contig && IsArticulationNode(i, xadj, adjncy, where_, bfslvl, bfsind, bfsmrk))
            {
                continue;
            }

            if (ctrl.minconn) {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            /* Find the most promising subdomain to move to */
            if (omode == OMODE_REFINE) {
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    gain = mynbrs[k as usize].ed - myrinfo.id;
                    if (gain >= 0
                        && ivecaxpylez(
                            ncon,
                            1,
                            vwgt + i * ncon,
                            pwgts + to * ncon,
                            maxpwgts + to * ncon,
                        ))
                    {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                cto = to;
                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    if ((mynbrs[j as usize].ed > mynbrs[k as usize].ed
                        && ivecaxpylez(
                            ncon,
                            1,
                            vwgt + i * ncon,
                            pwgts + to * ncon,
                            maxpwgts + to * ncon,
                        ))
                        || (mynbrs[j as usize].ed == mynbrs[k as usize].ed
                            && BetterBalanceKWay(
                                ncon,
                                vwgt + i * ncon,
                                ubfactors,
                                1,
                                pwgts + cto * ncon,
                                pijbm + cto * ncon,
                                1,
                                pwgts + to * ncon,
                                pijbm + to * ncon,
                            )))
                    {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                gain = mynbrs[k as usize].ed - myrinfo.id;
                if (!(gain > 0
                    || (gain == 0
                        && (BetterBalanceKWay(
                            ncon,
                            vwgt + i * ncon,
                            ubfactors,
                            -1,
                            pwgts + from * ncon,
                            pijbm + from * ncon,
                            1,
                            pwgts + to * ncon,
                            pijbm + to * ncon,
                        ) || (iii % 2 == 0 && safetos[to as usize] == 2)))))
                {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    if (ivecaxpylez(
                        ncon,
                        1,
                        vwgt + i * ncon,
                        pwgts + to * ncon,
                        maxpwgts + to * ncon,
                    ) || BetterBalanceKWay(
                        ncon,
                        vwgt + i * ncon,
                        ubfactors,
                        -1,
                        pwgts + from * ncon,
                        pijbm + from * ncon,
                        1,
                        pwgts + to * ncon,
                        pijbm + to * ncon,
                    )) {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                cto = to;
                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    if (BetterBalanceKWay(
                        ncon,
                        vwgt + i * ncon,
                        ubfactors,
                        1,
                        pwgts + cto * ncon,
                        pijbm + cto * ncon,
                        1,
                        pwgts + to * ncon,
                        pijbm + to * ncon,
                    )) {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                if (mynbrs[k as usize].ed - myrinfo.id < 0
                    && !BetterBalanceKWay(
                        ncon,
                        vwgt + i * ncon,
                        ubfactors,
                        -1,
                        pwgts + from * ncon,
                        pijbm + from * ncon,
                        1,
                        pwgts + to * ncon,
                        pijbm + to * ncon,
                    ))
                {
                    continue;
                }
            }

            /*=====================================================================
             * If we got here, we can now move the vertex from 'from' to 'to'
             *======================================================================*/
            graph.mincut -= mynbrs[k as usize].ed - myrinfo.id;
            nmoved += 1;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "\t\tMoving {:6} to {:3}. Gain: {:4}. Cut: {:6}\n",
                    i,
                    to,
                    mynbrs[k as usize].ed - myrinfo.id,
                    graph.mincut
                )
            );

            /* Update the subdomain connectivity information */
            if (ctrl.minconn) {
                /* take care of i's move itself */
                UpdateEdgeSubDomainGraph(
                    ctrl,
                    from,
                    to,
                    myrinfo.id - mynbrs[k as usize].ed,
                    &maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    me = where_[adjncy[j as usize] as usize];
                    if (me != from && me != to) {
                        UpdateEdgeSubDomainGraph(ctrl, from, me, -adjwgt[j as usize], &maxndoms);
                        UpdateEdgeSubDomainGraph(ctrl, to, me, adjwgt[j as usize], &maxndoms);
                    }
                }
            }

            /* Update ID/ED and BND related information for the moved vertex */
            iaxpy(ncon, 1, vwgt + i * ncon, 1, pwgts + to * ncon, 1);
            iaxpy(ncon, -1, vwgt + i * ncon, 1, pwgts + from * ncon, 1);
            UpdateMovedVertexInfoAndBND(
                i, from, k, to, myrinfo, mynbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                ii = adjncy[j as usize];
                me = where_[ii as usize];
                myrinfo = graph.ckrinfo + ii;

                oldnnbrs = myrinfo.nnbrs;

                UpdateAdjacentVertexInfoAndBND(
                    ctrl,
                    ii,
                    xadj[(ii + 1) as usize] - xadj[ii as usize],
                    me,
                    from,
                    to,
                    myrinfo,
                    adjwgt[j as usize],
                    nbnd,
                    bndptr,
                    bndind,
                    bndtype,
                );

                UpdateQueueInfo(
                    queue, vstatus, ii, me, from, to, myrinfo, oldnnbrs, nupd, updptr, updind,
                    bndtype,
                );

                assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }
        }

        graph.nbnd = nbnd;

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            assert!(updptr[updind[i as usize] as usize] != -1);
            assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if (ctrl.dbglvl & METIS_DBG_REFINE) {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}."
              " Nmoves: {:5}, Cut: {:6}, Vol: %6"PRIDX,
              imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), 
              ComputeLoadImbalance(graph, nparts, pijbm), 
              graph.nbnd, nmoved, graph.mincut, ComputeVolume(graph, where_));
            if (ctrl.minconn) {
                print!(
                    ", Doms: [({:3} {:4}) as usize]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if (nmoved == 0 || (omode == OMODE_REFINE && graph.mincut == oldcut)) {
            break;
        }
    }

    rpqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************/
/* K-way refinement that minimizes the communication volume. This is a
    greedy routine and the vertices are visited in decreasing gv order.

  \param graph is the graph that is being refined.
  \param niter is the number of refinement iterations.
  \param ffactor is the \em fudge-factor for allowing positive gain moves
         to violate the max-pwgt constraint.

*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Greedy_McKWayVolOptimize(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    niter: idx_t,
    ffactor: real_t,
    omode: idx_t,
) {
    /* Common variables to all types of kway-refinement/balancing routines */
    // idx_t i, ii, iii, j, k, l, pass, nvtxs, ncon, nparts, gain;
    // idx_t from, me, to, cto, oldcut;
    // idx_t *xadj, *vwgt, *adjncy;
    // idx_t *where_, *pwgts, *perm, *bndptr, *bndind, *minpwgts, *maxpwgts;
    // idx_t nmoved, nupd, *vstatus, *updptr, *updind;
    // idx_t maxndoms, *safetos=std::ptr::null_mut(), *nads=std::ptr::null_mut(), *doms=std::ptr::null_mut(), **adids=std::ptr::null_mut(), **adwgts=std::ptr::null_mut();
    // idx_t *bfslvl=std::ptr::null_mut(), *bfsind=std::ptr::null_mut(), *bfsmrk=std::ptr::null_mut();
    // idx_t bndtype = (omode == OMODE_REFINE ? BNDTYPE_REFINE : BNDTYPE_BALANCE);
    // real_t *ubfactors, *pijbm;
    // real_t origbal;

    /* Volume-specific/different variables */
    // ipq_t *queue;
    // idx_t oldvol, xgain;
    // idx_t *vmarker, *pmarker, *modind;
    // vkrinfo_t *myrinfo;
    // vnbr_t *mynbrs;

    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    bndptr = graph.bndptr;
    bndind = graph.bndind;
    where_ = graph.where_;
    pwgts = graph.pwgts;

    nparts = ctrl.nparts;
    pijbm = ctrl.pijbm;

    /* Determine the ubfactors. The method used is different based on omode.
    When OMODE_BALANCE, the ubfactors are those supplied by the user.
    When OMODE_REFINE, the ubfactors are the max of the current partition
    and the user-specified ones. */
    ubfactors = vec![0.0; ncon as usize];
    ComputeLoadImbalanceVec(graph, nparts, pijbm, ubfactors);
    origbal = rvecmaxdiff(ncon, ubfactors, ctrl.ubfactors);
    if (omode == OMODE_BALANCE) {
        rcopy(ncon, ctrl.ubfactors, ubfactors);
    } else {
        for i in (0)..(ncon) {
            ubfactors[i as usize] = (if ubfactors[i as usize] > ctrl.ubfactors[i as usize] {
                ubfactors[i as usize]
            } else {
                ctrl.ubfactors[i as usize]
            });
        }
    }

    /* Setup the weight intervals of the various subdomains */
    minpwgts = vec![0; nparts * ncon as usize];
    maxpwgts = vec![0; nparts * ncon as usize];

    for i in (0)..(nparts) {
        for j in (0)..(ncon) {
            maxpwgts[(i * ncon + j) as usize] = ctrl.tpwgts[(i * ncon + j) as usize]
                * graph.tvwgt[j as usize]
                * ubfactors[j as usize];
            /*minpwgts[(i*ncon+j) as usize]  = ctrl.tpwgts[(i*ncon+j) as usize]*graph.tvwgt[j as usize]*(.9/ubfactors[j as usize]); */
            minpwgts[(i * ncon + j) as usize] =
                ctrl.tpwgts[(i * ncon + j) as usize] * graph.tvwgt[j as usize] * 0.2;
        }
    }

    perm = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    safetos = vec![2; nparts as usize];

    if (ctrl.minconn) {
        ComputeSubDomainGraph(ctrl, graph);

        nads = ctrl.nads;
        adids = ctrl.adids;
        adwgts = ctrl.adwgts;
        doms = iset(nparts, 0, ctrl.pvec1);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    vstatus = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    updptr = vec![-1; nvtxs as usize];
    updind = vec![0; nvtxs as usize];

    if (ctrl.contig) {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    }

    /* Vol-refinement specific working arrays */
    modind = vec![0; nvtxs as usize];
    vmarker = vec![0; nvtxs as usize];
    pmarker = vec![-1; nparts as usize];

    if (ctrl.dbglvl & METIS_DBG_REFINE) {
        print!("{}: [({:6} {:6} {:6}) as usize], Bal: {:5.3}({:.3}),"
         ", Nv-Nb[({:6} {:6}) as usize], Cut: {:5}, Vol: {:5}, ({:})",
         (omode == OMODE_REFINE ? "GRV" : "GBV"),
         imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), imax(nparts*ncon, maxpwgts,1),
         ComputeLoadImbalance(graph, nparts, pijbm), origbal,
         graph.nvtxs, graph.nbnd, graph.mincut, graph.minvol, niter);
        if (ctrl.minconn) {
            print!(
                ", Doms: [({:3} {:4}) as usize]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    queue = ipqCreate(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        assert!(ComputeVolume(graph, where_) == graph.minvol);

        /* In balancing mode, exit as soon as balance is reached */
        if (omode == OMODE_BALANCE && IsBalanced(ctrl, graph, 0)) {
            break;
        }

        oldcut = graph.mincut;
        oldvol = graph.minvol;
        nupd = 0;

        if (ctrl.minconn) {
            maxndoms = imax(nparts, nads, 1);
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(graph.nbnd, perm, graph.nbnd / 4, 1);
        for ii in (0)..(graph.nbnd) {
            i = bndind[perm[ii as usize] as usize];
            ipqInsert(queue, i, graph.vkrinfo[i as usize].gv);
            vstatus[i as usize] = VPQSTATUS_PRESENT;
            ListInsert(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        nmoved = 0;
        for iii in (0).. {
            if ((i = ipqGetTop(queue)) == -1) {
                break;
            }
            vstatus[i as usize] = VPQSTATUS_EXTRACTED;

            myrinfo = graph.vkrinfo + i;
            mynbrs = ctrl.vnbrpool + myrinfo.inbr;

            from = where_[i as usize];

            /* Prevent moves that make 'from' domain underbalanced */
            if (omode == OMODE_REFINE) {
                if (myrinfo.nid > 0
                    && !ivecaxpygez(
                        ncon,
                        -1,
                        vwgt + i * ncon,
                        pwgts + from * ncon,
                        minpwgts + from * ncon,
                    ))
                {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                if (!ivecaxpygez(
                    ncon,
                    -1,
                    vwgt + i * ncon,
                    pwgts + from * ncon,
                    minpwgts + from * ncon,
                )) {
                    continue;
                }
            }

            if (ctrl.contig && IsArticulationNode(i, xadj, adjncy, where_, bfslvl, bfsind, bfsmrk))
            {
                continue;
            }

            if (ctrl.minconn) {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            xgain = (if myrinfo.nid == 0 && myrinfo.ned > 0 {
                graph.vsize[i as usize]
            } else {
                0
            });

            /* Find the most promising subdomain to move to */
            if (omode == OMODE_REFINE) {
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    gain = mynbrs[k as usize].gv + xgain;
                    if (gain >= 0
                        && ivecaxpylez(
                            ncon,
                            1,
                            vwgt + i * ncon,
                            pwgts + to * ncon,
                            maxpwgts + to * ncon,
                        ))
                    {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                cto = to;
                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    gain = mynbrs[j as usize].gv + xgain;
                    if ((mynbrs[j as usize].gv > mynbrs[k as usize].gv
                        && ivecaxpylez(
                            ncon,
                            1,
                            vwgt + i * ncon,
                            pwgts + to * ncon,
                            maxpwgts + to * ncon,
                        ))
                        || (mynbrs[j as usize].gv == mynbrs[k as usize].gv
                            && mynbrs[j as usize].ned > mynbrs[k as usize].ned
                            && ivecaxpylez(
                                ncon,
                                1,
                                vwgt + i * ncon,
                                pwgts + to * ncon,
                                maxpwgts + to * ncon,
                            ))
                        || (mynbrs[j as usize].gv == mynbrs[k as usize].gv
                            && mynbrs[j as usize].ned == mynbrs[k as usize].ned
                            && BetterBalanceKWay(
                                ncon,
                                vwgt + i * ncon,
                                ubfactors,
                                1,
                                pwgts + cto * ncon,
                                pijbm + cto * ncon,
                                1,
                                pwgts + to * ncon,
                                pijbm + to * ncon,
                            )))
                    {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                j = 0;
                if (xgain + mynbrs[k as usize].gv > 0 || mynbrs[k as usize].ned - myrinfo.nid > 0) {
                    j = 1;
                } else if (mynbrs[k as usize].ned - myrinfo.nid == 0) {
                    if ((iii % 2 == 0 && safetos[to as usize] == 2)
                        || BetterBalanceKWay(
                            ncon,
                            vwgt + i * ncon,
                            ubfactors,
                            -1,
                            pwgts + from * ncon,
                            pijbm + from * ncon,
                            1,
                            pwgts + to * ncon,
                            pijbm + to * ncon,
                        ))
                    {
                        j = 1;
                    }
                }
                if (j == 0) {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                for k in ((myrinfo.nnbrs - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[k as usize].pid]) {
                        continue;
                    }
                    if (ivecaxpylez(
                        ncon,
                        1,
                        vwgt + i * ncon,
                        pwgts + to * ncon,
                        maxpwgts + to * ncon,
                    ) || BetterBalanceKWay(
                        ncon,
                        vwgt + i * ncon,
                        ubfactors,
                        -1,
                        pwgts + from * ncon,
                        pijbm + from * ncon,
                        1,
                        pwgts + to * ncon,
                        pijbm + to * ncon,
                    )) {
                        break;
                    }
                }
                if (k < 0) {
                    continue;
                } /* break out if you did not find a candidate */

                cto = to;
                for j in ((k - 1)..=(0)).rev() {
                    if (!safetos[to = mynbrs[j as usize].pid]) {
                        continue;
                    }
                    if (BetterBalanceKWay(
                        ncon,
                        vwgt + i * ncon,
                        ubfactors,
                        1,
                        pwgts + cto * ncon,
                        pijbm + cto * ncon,
                        1,
                        pwgts + to * ncon,
                        pijbm + to * ncon,
                    )) {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                if ((xgain + mynbrs[k as usize].gv < 0
                    || (xgain + mynbrs[k as usize].gv == 0
                        && mynbrs[k as usize].ned - myrinfo.nid < 0))
                    && !BetterBalanceKWay(
                        ncon,
                        vwgt + i * ncon,
                        ubfactors,
                        -1,
                        pwgts + from * ncon,
                        pijbm + from * ncon,
                        1,
                        pwgts + to * ncon,
                        pijbm + to * ncon,
                    ))
                {
                    continue;
                }
            }

            /*=====================================================================
             * If we got here, we can now move the vertex from 'from' to 'to'
             *======================================================================*/
            graph.mincut -= mynbrs[k as usize].ned - myrinfo.nid;
            graph.minvol -= (xgain + mynbrs[k as usize].gv);
            where_[i as usize] = to;
            nmoved += 1;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!("\t\tMoving {:6} from {:3} to {:3}. "
                 "Gain: [({:4} {:4}) as usize]. Cut: {:6}, Vol: {:6}\n", 
              i, from, to, xgain+mynbrs[k as usize].gv, mynbrs[k as usize].ned-myrinfo.nid, 
              graph.mincut, graph.minvol)
            );

            /* Update the subdomain connectivity information */
            if (ctrl.minconn) {
                /* take care of i's move itself */
                UpdateEdgeSubDomainGraph(
                    ctrl,
                    from,
                    to,
                    myrinfo.nid - mynbrs[k as usize].ned,
                    &maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    me = where_[adjncy[j as usize] as usize];
                    if (me != from && me != to) {
                        UpdateEdgeSubDomainGraph(ctrl, from, me, -1, &maxndoms);
                        UpdateEdgeSubDomainGraph(ctrl, to, me, 1, &maxndoms);
                    }
                }
            }

            /* Update pwgts */
            iaxpy(ncon, 1, vwgt + i * ncon, 1, pwgts + to * ncon, 1);
            iaxpy(ncon, -1, vwgt + i * ncon, 1, pwgts + from * ncon, 1);

            /* Update the id/ed/gains/bnd/queue of potentially affected nodes */
            KWayVolUpdate(
                ctrl, graph, i, from, to, queue, vstatus, &nupd, updptr, updind, bndtype, vmarker,
                pmarker, modind,
            );

            /*CheckKWayVolPartitionParams(ctrl, graph); */
        }

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            assert!(updptr[updind[i as usize] as usize] != -1);
            assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if (ctrl.dbglvl & METIS_DBG_REFINE) {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}."
              " Nmoves: {:5}, Cut: {:6}, Vol: %6"PRIDX,
              imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), 
              ComputeLoadImbalance(graph, nparts, pijbm), 
              graph.nbnd, nmoved, graph.mincut, graph.minvol);
            if (ctrl.minconn) {
                print!(
                    ", Doms: [({:3} {:4}) as usize]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if (nmoved == 0
            || (omode == OMODE_REFINE && graph.minvol == oldvol && graph.mincut == oldcut))
        {
            break;
        }
    }

    ipqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************/
/* This function performs an approximate articulation vertex test.
It assumes that the bfslvl, bfsind, and bfsmrk arrays are initialized
appropriately. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn IsArticulationNode(
    i: idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    where_: *mut idx_t,
    bfslvl: *mut idx_t,
    bfsind: *mut idx_t,
    bfsmrk: *mut idx_t,
) -> idx_t {
    // idx_t ii, j, k=0, head, tail, nhits, tnhits, from, BFSDEPTH=5;

    from = where_[i as usize];

    /* Determine if the vertex is safe to move from a contiguity standpoint */
    tnhits = 0;
    for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
        if (where_[adjncy[j as usize] as usize] == from) {
            assert!(bfsmrk[adjncy[j as usize] as usize] == 0);
            assert!(bfslvl[adjncy[j as usize] as usize] == 0);
            bfsmrk[k = adjncy[j as usize] as usize] = 1;
            tnhits += 1;
        }
    }

    /* Easy cases */
    if (tnhits == 0) {
        return 0;
    }
    if (tnhits == 1) {
        bfsmrk[k as usize] = 0;
        return 0;
    }

    assert!(bfslvl[i as usize] == 0);
    bfslvl[i as usize] = 1;

    bfsind[0 as usize] = k; /* That was the last one from the previous loop */
    bfslvl[k as usize] = 1;
    bfsmrk[k as usize] = 0;
    head = 0;
    tail = 1;

    /* Do a limited BFS traversal to see if you can get to all the other nodes */
    nhits = 1;
    while head < tail {
        ii = bfsind[(head) as usize];
        head += 1;
        for j in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
            if (where_[k = adjncy[j as usize] as usize] == from) {
                if (bfsmrk[k as usize]) {
                    bfsmrk[k as usize] = 0;
                    nhits += 1;
                    if (nhits == tnhits) {
                        break;
                    }
                }
                if (bfslvl[k as usize] == 0 && bfslvl[ii as usize] < BFSDEPTH) {
                    bfsind[(tail) as usize] = k;
                    tail += 1;
                    bfslvl[k as usize] = bfslvl[ii as usize] + 1;
                }
            }
        }
        if (nhits == tnhits) {
            break;
        }
    }

    /* Reset the various BFS related arrays */
    bfslvl[i as usize] = 0;
    for j in (0)..(tail) {
        bfslvl[bfsind[j as usize] as usize] = 0;
    }

    /* Reset the bfsmrk array for the next vertex when has not already being cleared */
    if (nhits < tnhits) {
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            if (where_[adjncy[j as usize] as usize] == from) {
                bfsmrk[adjncy[j as usize] as usize] = 0;
            }
        }
    }

    return (nhits != tnhits);
}

/*************************************************************************/
/*
 This function updates the edge and volume gains due to a vertex movement.
 v from 'from' to 'to'.

 \param ctrl is the control structure.
 \param graph is the graph being partitioned.
 \param v is the vertex that is moving.
 \param from is the original partition of v.
 \param to is the new partition of v.
 \param queue is the priority queue. If the queue is std::ptr::null_mut(), no priority-queue
        related updates are performed.
 \param vstatus is an array that marks the status of the vertex in terms
        of the priority queue. If queue is std::ptr::null_mut(), this parameter is ignored.
 \param r_nqupd is the number of vertices that have been inserted/removed
        from the queue. If queue is std::ptr::null_mut(), this parameter is ignored.
 \param updptr stores the index of each vertex in updind. If queue is std::ptr::null_mut(),
        this parameter is ignored.
 \param updind is the list of vertices that have been inserted/removed from
        the queue. If queue is std::ptr::null_mut(), this parameter is ignored.
 \param vmarker is of size nvtxs and is used internally as a temporary array.
        On entry and return all of its entries are 0.
 \param pmarker is of size nparts and is used internally as a temporary marking
        array. On entry and return all of its entries are -1.
 \param modind is an array of size nvtxs and is used to keep track of the
        list of vertices whose gains need to be updated.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn KWayVolUpdate(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    v: idx_t,
    from: idx_t,
    to: idx_t,
    queue: *mut ipq_t,
    vstatus: *mut idx_t,
    r_nupd: *mut idx_t,
    updptr: *mut idx_t,
    updind: *mut idx_t,
    bndtype: idx_t,
    vmarker: *mut idx_t,
    pmarker: *mut idx_t,
    modind: *mut idx_t,
) {
    // idx_t i, ii, iii, j, jj, k, kk, l, u, nmod, other, me, myidx;
    // idx_t *xadj, *vsize, *adjncy, *where_;
    // vkrinfo_t *myrinfo, *orinfo;
    // vnbr_t *mynbrs, *onbrs;

    xadj = graph.xadj;
    adjncy = graph.adjncy;
    vsize = graph.vsize;
    where_ = graph.where_;

    myrinfo = graph.vkrinfo + v;
    mynbrs = ctrl.vnbrpool + myrinfo.inbr;

    /*======================================================================
     * Remove the contributions on the gain made by 'v'.
     *=====================================================================*/
    for k in (0)..(myrinfo.nnbrs) {
        pmarker[mynbrs[k as usize].pid] = k;
    }
    pmarker[from as usize] = k;

    myidx = pmarker[to as usize]; /* Keep track of the index in mynbrs of the 'to' domain */

    for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
        ii = adjncy[j as usize];
        other = where_[ii as usize];
        orinfo = graph.vkrinfo + ii;
        onbrs = ctrl.vnbrpool + orinfo.inbr;

        if (other == from) {
            for k in (0)..(orinfo.nnbrs) {
                if (pmarker[onbrs[k as usize].pid] == -1) {
                    onbrs[k as usize].gv += vsize[v as usize];
                }
            }
        } else {
            assert!(pmarker[other as usize] != -1);

            if (mynbrs[pmarker[other as usize] as usize].ned > 1) {
                for k in (0)..(orinfo.nnbrs) {
                    if (pmarker[onbrs[k as usize].pid] == -1) {
                        onbrs[k as usize].gv += vsize[v as usize];
                    }
                }
            } else {
                /* There is only one connection */
                for k in (0)..(orinfo.nnbrs) {
                    if (pmarker[onbrs[k as usize].pid] != -1) {
                        onbrs[k as usize].gv -= vsize[v as usize];
                    }
                }
            }
        }
    }

    for k in (0)..(myrinfo.nnbrs) {
        pmarker[mynbrs[k as usize].pid] = -1;
    }
    pmarker[from as usize] = -1;

    /*======================================================================
     * Update the id/ed of vertex 'v'
     *=====================================================================*/
    if (myidx == -1) {
        myidx = myrinfo.nnbrs += 1;
        assert!(myidx < xadj[(v + 1) as usize] - xadj[v as usize]);
        mynbrs[myidx as usize].ned = 0;
    }
    myrinfo.ned += myrinfo.nid - mynbrs[myidx as usize].ned;
    SWAP(myrinfo.nid, mynbrs[myidx as usize].ned, j);
    if (mynbrs[myidx as usize].ned == 0) {
        myrinfo -= 1;
        mynbrs[myidx as usize] = mynbrs[(myrinfo.nnbrs) as usize];
    } else {
        mynbrs[myidx as usize].pid = from;
    }

    /*======================================================================
     * Update the degrees of adjacent vertices and their volume gains
     *=====================================================================*/
    vmarker[v as usize] = 1;
    modind[0 as usize] = v;
    nmod = 1;
    for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
        ii = adjncy[j as usize];
        me = where_[ii as usize];

        if (!vmarker[ii as usize]) {
            /* The marking is done for boundary and max gv calculations */
            vmarker[ii as usize] = 2;
            modind[(nmod) as usize] = ii;
            nmod += 1;
        }

        myrinfo = graph.vkrinfo + ii;
        if (myrinfo.inbr == -1) {
            myrinfo.inbr = vnbrpoolGetNext(ctrl, xadj[(ii + 1) as usize] - xadj[ii as usize]);
        }
        mynbrs = ctrl.vnbrpool + myrinfo.inbr;

        if (me == from) {
            INC_DEC(myrinfo.ned, myrinfo.nid, 1);
        } else if (me == to) {
            INC_DEC(myrinfo.nid, myrinfo.ned, 1);
        }

        /* Remove the edgeweight from the 'pid == from' entry of the vertex */
        if (me != from) {
            for k in (0)..(myrinfo.nnbrs) {
                if (mynbrs[k as usize].pid == from) {
                    if (mynbrs[k as usize].ned == 1) {
                        myrinfo -= 1;
                        mynbrs[k as usize] = mynbrs[(myrinfo.nnbrs) as usize];
                        vmarker[ii as usize] = 1; /* You do a complete .gv calculation */

                        /* All vertices adjacent to 'ii' need to be updated */
                        for jj in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
                            u = adjncy[jj as usize];
                            other = where_[u as usize];
                            orinfo = graph.vkrinfo + u;
                            onbrs = ctrl.vnbrpool + orinfo.inbr;

                            for kk in (0)..(orinfo.nnbrs) {
                                if (onbrs[kk as usize].pid == from) {
                                    onbrs[kk as usize].gv -= vsize[ii as usize];
                                    if (!vmarker[u as usize]) {
                                        /* Need to update boundary etc */
                                        vmarker[u as usize] = 2;
                                        modind[(nmod) as usize] = u;
                                        nmod += 1;
                                    }
                                    break;
                                }
                            }
                        }
                    } else {
                        mynbrs[k as usize].ned;
                        ned -= 1;

                        /* Update the gv due to single 'ii' connection to 'from' */
                        if (mynbrs[k as usize].ned == 1) {
                            /* find the vertex 'u' that 'ii' was connected into 'from' */
                            for jj in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
                                u = adjncy[jj as usize];
                                other = where_[u as usize];

                                if (other == from) {
                                    orinfo = graph.vkrinfo + u;
                                    onbrs = ctrl.vnbrpool + orinfo.inbr;

                                    /* The following is correct because domains in common
                                    between ii and u will lead to a reduction over the
                                    previous gain, where_as domains only in u but not in
                                    ii, will lead to no change as opposed to the earlier
                                    increase */
                                    for kk in (0)..(orinfo.nnbrs) {
                                        onbrs[kk as usize].gv += vsize[ii as usize];
                                    }

                                    if (!vmarker[u as usize]) {
                                        /* Need to update boundary etc */
                                        vmarker[u as usize] = 2;
                                        modind[(nmod) as usize] = u;
                                        nmod += 1;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }

        /* Add the edgeweight to the 'pid == to' entry of the vertex */
        if (me != to) {
            for k in (0)..(myrinfo.nnbrs) {
                if (mynbrs[k as usize].pid == to) {
                    mynbrs[k as usize].ned += 1;

                    /* Update the gv due to non-single 'ii' connection to 'to' */
                    if (mynbrs[k as usize].ned == 2) {
                        /* find the vertex 'u' that 'ii' was connected into 'to' */
                        for jj in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
                            u = adjncy[jj as usize];
                            other = where_[u as usize];

                            if (u != v && other == to) {
                                orinfo = graph.vkrinfo + u;
                                onbrs = ctrl.vnbrpool + orinfo.inbr;
                                for kk in (0)..(orinfo.nnbrs) {
                                    onbrs[kk as usize].gv -= vsize[ii as usize];
                                }

                                if (!vmarker[u as usize]) {
                                    /* Need to update boundary etc */
                                    vmarker[u as usize] = 2;
                                    modind[(nmod) as usize] = u;
                                    nmod += 1;
                                }
                                break;
                            }
                        }
                    }
                    break;
                }
            }

            if (k == myrinfo.nnbrs) {
                mynbrs[myrinfo.nnbrs].pid = to;
                mynbrs[(myrinfo.nnbrs) as usize].ned = 1;
                nnbrs += 1;
                vmarker[ii as usize] = 1; /* You do a complete .gv calculation */

                /* All vertices adjacent to 'ii' need to be updated */
                for jj in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
                    u = adjncy[jj as usize];
                    other = where_[u as usize];
                    orinfo = graph.vkrinfo + u;
                    onbrs = ctrl.vnbrpool + orinfo.inbr;

                    for kk in (0)..(orinfo.nnbrs) {
                        if (onbrs[kk as usize].pid == to) {
                            onbrs[kk as usize].gv += vsize[ii as usize];
                            if (!vmarker[u as usize]) {
                                /* Need to update boundary etc */
                                vmarker[u as usize] = 2;
                                modind[(nmod) as usize] = u;
                                nmod += 1;
                            }
                            break;
                        }
                    }
                }
            }
        }

        assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
    }

    /*======================================================================
     * Add the contributions on the volume gain due to 'v'
     *=====================================================================*/
    myrinfo = graph.vkrinfo + v;
    mynbrs = ctrl.vnbrpool + myrinfo.inbr;
    for k in (0)..(myrinfo.nnbrs) {
        pmarker[mynbrs[k as usize].pid] = k;
    }
    pmarker[to as usize] = k;

    for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
        ii = adjncy[j as usize];
        other = where_[ii as usize];
        orinfo = graph.vkrinfo + ii;
        onbrs = ctrl.vnbrpool + orinfo.inbr;

        if (other == to) {
            for k in (0)..(orinfo.nnbrs) {
                if (pmarker[onbrs[k as usize].pid] == -1) {
                    onbrs[k as usize].gv -= vsize[v as usize];
                }
            }
        } else {
            assert!(pmarker[other as usize] != -1);

            if (mynbrs[pmarker[other as usize] as usize].ned > 1) {
                for k in (0)..(orinfo.nnbrs) {
                    if (pmarker[onbrs[k as usize].pid] == -1) {
                        onbrs[k as usize].gv -= vsize[v as usize];
                    }
                }
            } else {
                /* There is only one connection */
                for k in (0)..(orinfo.nnbrs) {
                    if (pmarker[onbrs[k as usize].pid] != -1) {
                        onbrs[k as usize].gv += vsize[v as usize];
                    }
                }
            }
        }
    }
    for k in (0)..(myrinfo.nnbrs) {
        pmarker[mynbrs[k as usize].pid] = -1;
    }
    pmarker[to as usize] = -1;

    /*======================================================================
     * Recompute the volume information of the 'hard' nodes, and update the
     * max volume gain for all the modified vertices and the priority queue
     *=====================================================================*/
    for iii in (0)..(nmod) {
        i = modind[iii as usize];
        me = where_[i as usize];

        myrinfo = graph.vkrinfo + i;
        mynbrs = ctrl.vnbrpool + myrinfo.inbr;

        if (vmarker[i as usize] == 1) {
            /* Only complete gain updates go through */
            for k in (0)..(myrinfo.nnbrs) {
                mynbrs[k as usize].gv = 0;
            }

            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                ii = adjncy[j as usize];
                other = where_[ii as usize];
                orinfo = graph.vkrinfo + ii;
                onbrs = ctrl.vnbrpool + orinfo.inbr;

                for kk in (0)..(orinfo.nnbrs) {
                    pmarker[onbrs[kk as usize].pid] = kk;
                }
                pmarker[other as usize] = 1;

                if (me == other) {
                    /* Find which domains 'i' is connected and 'ii' is not and update their gain */
                    for k in (0)..(myrinfo.nnbrs) {
                        if (pmarker[mynbrs[k as usize].pid] == -1) {
                            mynbrs[k as usize].gv -= vsize[ii as usize];
                        }
                    }
                } else {
                    assert!(pmarker[me as usize] != -1);

                    /* I'm the only connection of 'ii' in 'me' */
                    if (onbrs[pmarker[me as usize] as usize].ned == 1) {
                        /* Increase the gains for all the common domains between 'i' and 'ii' */
                        for k in (0)..(myrinfo.nnbrs) {
                            if (pmarker[mynbrs[k as usize].pid] != -1) {
                                mynbrs[k as usize].gv += vsize[ii as usize];
                            }
                        }
                    } else {
                        /* Find which domains 'i' is connected and 'ii' is not and update their gain */
                        for k in (0)..(myrinfo.nnbrs) {
                            if (pmarker[mynbrs[k as usize].pid] == -1) {
                                mynbrs[k as usize].gv -= vsize[ii as usize];
                            }
                        }
                    }
                }

                for kk in (0)..(orinfo.nnbrs) {
                    pmarker[onbrs[kk as usize].pid] = -1;
                }
                pmarker[other as usize] = -1;
            }
        }

        /* Compute the overall gv for that node */
        myrinfo.gv = IDX_MIN;
        for k in (0)..(myrinfo.nnbrs) {
            if (mynbrs[k as usize].gv > myrinfo.gv) {
                myrinfo.gv = mynbrs[k as usize].gv;
            }
        }

        /* Add the xtra gain due to id == 0 */
        if (myrinfo.ned > 0 && myrinfo.nid == 0) {
            myrinfo.gv += vsize[i as usize];
        }

        /*======================================================================
         * Maintain a consistent boundary
         *=====================================================================*/
        if (bndtype == BNDTYPE_REFINE) {
            if (myrinfo.gv >= 0 && graph.bndptr[i as usize] == -1) {
                BNDInsert(graph.nbnd, graph.bndind, graph.bndptr, i);
            }

            if (myrinfo.gv < 0 && graph.bndptr[i as usize] != -1) {
                BNDDelete(graph.nbnd, graph.bndind, graph.bndptr, i);
            }
        } else {
            if (myrinfo.ned > 0 && graph.bndptr[i as usize] == -1) {
                BNDInsert(graph.nbnd, graph.bndind, graph.bndptr, i);
            }

            if (myrinfo.ned == 0 && graph.bndptr[i as usize] != -1) {
                BNDDelete(graph.nbnd, graph.bndind, graph.bndptr, i);
            }
        }

        /*======================================================================
         * Update the priority queue appropriately (if allowed)
         *=====================================================================*/
        if (queue != std::ptr::null_mut()) {
            if (vstatus[i as usize] != VPQSTATUS_EXTRACTED) {
                if (graph.bndptr[i as usize] != -1) {
                    /* In-boundary vertex */
                    if (vstatus[i as usize] == VPQSTATUS_PRESENT) {
                        ipqUpdate(queue, i, myrinfo.gv);
                    } else {
                        ipqInsert(queue, i, myrinfo.gv);
                        vstatus[i as usize] = VPQSTATUS_PRESENT;
                        ListInsert(*r_nupd, updind, updptr, i);
                    }
                } else {
                    /* Off-boundary vertex */
                    if (vstatus[i as usize] == VPQSTATUS_PRESENT) {
                        ipqDelete(queue, i);
                        vstatus[i as usize] = VPQSTATUS_NOTPRESENT;
                        ListDelete(*r_nupd, updind, updptr, i);
                    }
                }
            }
        }

        vmarker[i as usize] = 0;
    }
}

/*************************************************************************/
/* K-way partitioning optimization in which the vertices are visited in
    decreasing ed/sqrt(nnbrs)-id order. Note this is just an
    approximation, as the ed is often split across different subdomains
    and the sqrt(nnbrs) is just a crude approximation.

  \param graph is the graph that is being refined.
  \param niter is the number of refinement iterations.
  \param ffactor is the \em fudge-factor for allowing positive gain moves
         to violate the max-pwgt constraint.
  \param omode is the type of optimization that will performed among
         OMODE_REFINE and OMODE_BALANCE


*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Greedy_KWayEdgeStats(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    /* Common variables to all types of kway-refinement/balancing routines */
    // idx_t i, ii, iii, j, k, l, nvtxs, nparts, gain, u, v, uw, vw;
    // idx_t *xadj, *adjncy, *adjwgt, *vwgt;
    // idx_t *where_, *pwgts, *bndptr, *bndind, *minpwgts, *maxpwgts;
    // idx_t nbnd;
    // ckrinfo_t *urinfo, *vrinfo;
    // cnbr_t *unbrs, *vnbrs;
    // real_t *tpwgts, ubfactor;

    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    vwgt = graph.vwgt;
    adjwgt = graph.adjwgt;

    bndind = graph.bndind;
    bndptr = graph.bndptr;

    where_ = graph.where_;
    pwgts = graph.pwgts;

    nparts = ctrl.nparts;
    tpwgts = ctrl.tpwgts;

    /* Setup the weight intervals of the various subdomains */
    minpwgts = vec![0; nparts as usize];
    maxpwgts = vec![0; nparts as usize];

    ubfactor = ctrl.ubfactors[0 as usize];
    for i in (0)..(nparts) {
        maxpwgts[i as usize] = tpwgts[i as usize] * graph.tvwgt[0 as usize] * ubfactor;
        minpwgts[i as usize] = tpwgts[i as usize] * graph.tvwgt[0 as usize] * (0.95 / ubfactor);
    }

    /* go and determine the positive gain valid swaps */
    nbnd = graph.nbnd;

    for ii in (0)..(nbnd) {
        u = bndind[ii as usize];
        uw = where_[u as usize];

        urinfo = graph.ckrinfo + u;
        unbrs = ctrl.cnbrpool + urinfo.inbr;

        for j in (xadj[u as usize])..(xadj[(u + 1) as usize]) {
            v = adjncy[j as usize];
            vw = where_[v as usize];

            vrinfo = graph.ckrinfo + v;
            vnbrs = ctrl.cnbrpool + vrinfo.inbr;

            if (uw == vw) {
                continue;
            }
            if (pwgts[uw as usize] - vwgt[u as usize] + vwgt[v as usize] > maxpwgts[uw as usize]
                || pwgts[vw as usize] - vwgt[v as usize] + vwgt[u as usize] > maxpwgts[vw as usize])
            {
                continue;
            }

            for k in ((urinfo.nnbrs - 1)..=(0)).rev() {
                if (unbrs[k as usize].pid == vw) {
                    break;
                }
            }
            if (k < 0) {
                print!("Something went wrong!\n");
            }
            gain = unbrs[k as usize].ed - urinfo.id;

            for k in ((vrinfo.nnbrs - 1)..=(0)).rev() {
                if (vnbrs[k as usize].pid == uw) {
                    break;
                }
            }
            if (k < 0) {
                print!("Something went wrong!\n");
            }
            gain += vnbrs[k as usize].ed - vrinfo.id;

            gain -= 2 * adjwgt[j as usize];

            if (gain > 0) {
                print!(
                    "  Gain: {:} for moving ({:}, {:}) between ({:}, {:})\n",
                    gain, u, v, uw, vw
                );
            }
        }
    }

    // WCOREPOP;
}

/*************************************************************************/
/* K-way partitioning optimization in which the vertices are visited in
    random order and the best edge is selected to swap its incident vertices

  \param graph is the graph that is being refined.
  \param niter is the number of refinement iterations.

*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn Greedy_KWayEdgeCutOptimize(ctrl: *mut ctrl_t, graph: *mut graph_t, niter: idx_t) {
    /* Common variables to all types of kway-refinement/balancing routines */
    // idx_t ii, j, k, pass, nvtxs, nparts, u, v, uw, vw, gain, bestgain, jbest;
    // idx_t from, me, to, oldcut, nmoved;
    // idx_t *xadj, *adjncy, *adjwgt, *vwgt;
    // idx_t *where_, *pwgts, *perm, *bndptr, *bndind, *minpwgts, *maxpwgts;
    // idx_t bndtype = BNDTYPE_REFINE;
    // real_t *tpwgts, ubfactor;

    /* Edgecut-specific/different variables */
    // idx_t nbnd, oldnnbrs;
    // ckrinfo_t *myrinfo, *urinfo, *vrinfo;
    // cnbr_t *unbrs, *vnbrs;

    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;
    vwgt = graph.vwgt;

    bndind = graph.bndind;
    bndptr = graph.bndptr;

    where_ = graph.where_;
    pwgts = graph.pwgts;

    nparts = ctrl.nparts;
    tpwgts = ctrl.tpwgts;

    /* Setup the weight intervals of the various subdomains */
    minpwgts = vec![0; nparts as usize];
    maxpwgts = vec![0; nparts as usize];

    ubfactor = gk_max(
        ctrl.ubfactors[0 as usize],
        ComputeLoadImbalance(graph, nparts, ctrl.pijbm),
    );
    for k in (0)..(nparts) {
        maxpwgts[k as usize] = tpwgts[k as usize] * graph.tvwgt[0 as usize] * ubfactor;
        minpwgts[k as usize] = tpwgts[k as usize] * graph.tvwgt[0 as usize] * (1.0 / ubfactor);
    }

    perm = vec![0; nvtxs as usize];

    if (ctrl.dbglvl & METIS_DBG_REFINE) {
        print!("GRE: [({:6} {:6}) as usize]-[({:6} {:6}) as usize], Bal: {:5.3}," 
            " Nv-Nb[({:6} {:6}) as usize], Cut: {:6}\n",
            pwgts[(iargmin(nparts, pwgts,1)) as usize], imax(nparts, pwgts,1), minpwgts[0 as usize], maxpwgts[0 as usize], 
            ComputeLoadImbalance(graph, nparts, ctrl.pijbm), 
            graph.nvtxs, graph.nbnd, graph.mincut);
    }

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        GKassert!(ComputeCut(graph, where_) == graph.mincut);

        oldcut = graph.mincut;
        nbnd = graph.nbnd;
        nmoved = 0;

        /* Insert the boundary vertices in the priority queue */
        /* Visit the vertices in random order and see if you can swap them */
        irandArrayPermute(nvtxs, perm, nbnd, 1);
        for ii in (0)..(nvtxs) {
            if (bndptr[u = perm[ii as usize] as usize] == -1) {
                continue;
            }

            uw = where_[u as usize];

            urinfo = graph.ckrinfo + u;
            unbrs = ctrl.cnbrpool + urinfo.inbr;

            bestgain = 0;
            jbest = -1;
            for j in (xadj[u as usize])..(xadj[(u + 1) as usize]) {
                v = adjncy[j as usize];
                vw = where_[v as usize];

                if (uw == vw) {
                    continue;
                }
                if (pwgts[uw as usize] - vwgt[u as usize] + vwgt[v as usize]
                    > maxpwgts[uw as usize]
                    || pwgts[vw as usize] - vwgt[v as usize] + vwgt[u as usize]
                        > maxpwgts[vw as usize])
                {
                    continue;
                }
                if (pwgts[uw as usize] - vwgt[u as usize] + vwgt[v as usize]
                    < minpwgts[uw as usize]
                    || pwgts[vw as usize] - vwgt[v as usize] + vwgt[u as usize]
                        < minpwgts[vw as usize])
                {
                    continue;
                }

                vrinfo = graph.ckrinfo + v;
                vnbrs = ctrl.cnbrpool + vrinfo.inbr;

                gain = -2 * adjwgt[j as usize];

                for k in ((urinfo.nnbrs - 1)..=(0)).rev() {
                    if (unbrs[k as usize].pid == vw) {
                        break;
                    }
                }
                GKassert!(k >= 0);
                gain += unbrs[k as usize].ed - urinfo.id;

                for k in ((vrinfo.nnbrs - 1)..=(0)).rev() {
                    if (vnbrs[k as usize].pid == uw) {
                        break;
                    }
                }
                GKassert!(k >= 0);
                gain += vnbrs[k as usize].ed - vrinfo.id;

                if (gain > bestgain && vnbrs[k as usize].ed > adjwgt[j as usize]) {
                    bestgain = gain;
                    jbest = j;
                }
            }

            if (jbest == -1) {
                continue;
            } /* no valid positive swap */

            /*=====================================================================
             * If we got here, we can now swap the vertices
             *======================================================================*/
            v = adjncy[jbest as usize];
            vw = where_[v as usize];

            vrinfo = graph.ckrinfo + v;
            vnbrs = ctrl.cnbrpool + vrinfo.inbr;

            /* move u to v's partition */
            for k in ((urinfo.nnbrs - 1)..=(0)).rev() {
                if (unbrs[k as usize].pid == vw) {
                    break;
                }
            }
            GKassert!(k >= 0);

            from = uw;
            to = vw;

            graph.mincut -= unbrs[k as usize].ed - urinfo.id;
            nmoved += 1;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "\t\tMoving {:6} from {:3} to {:3} [({:6} {:6}) as usize]. Gain: {:4}. Cut: {:6}\n",
                    u,
                    from,
                    to,
                    pwgts[from as usize],
                    pwgts[to as usize],
                    unbrs[k as usize].ed - urinfo.id,
                    graph.mincut
                )
            );

            /* Update ID/ED and BND related information for the moved vertex */
            INC_DEC(pwgts[to as usize], pwgts[from as usize], vwgt[u as usize]);
            UpdateMovedVertexInfoAndBND(
                u, from, k, to, urinfo, unbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[u as usize])..(xadj[(u + 1) as usize]) {
                ii = adjncy[j as usize];
                me = where_[ii as usize];
                myrinfo = graph.ckrinfo + ii;

                oldnnbrs = myrinfo.nnbrs;

                UpdateAdjacentVertexInfoAndBND(
                    ctrl,
                    ii,
                    xadj[(ii + 1) as usize] - xadj[ii as usize],
                    me,
                    from,
                    to,
                    myrinfo,
                    adjwgt[j as usize],
                    nbnd,
                    bndptr,
                    bndind,
                    bndtype,
                );

                assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }

            /* move v to u's partition */
            for k in ((vrinfo.nnbrs - 1)..=(0)).rev() {
                if (vnbrs[k as usize].pid == uw) {
                    break;
                }
            }
            GKassert!(k >= 0);
            // #ifdef XXX
            //       if (k < 0) { /* that was removed, go and re-insert it */
            //         k = vrinfo.nnbrs += 1;
            //         vnbrs[k as usize].pid = uw;
            //         vnbrs[k as usize].ed  = 0;
            //       }
            // #endif

            from = vw;
            to = uw;

            graph.mincut -= vnbrs[k as usize].ed - vrinfo.id;
            nmoved += 1;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!(
                    "\t\tMoving {:6} from {:3} to {:3} [({:6} {:6}) as usize]. Gain: {:4}. Cut: {:6}\n",
                    v,
                    from,
                    to,
                    pwgts[from as usize],
                    pwgts[to as usize],
                    vnbrs[k as usize].ed - vrinfo.id,
                    graph.mincut
                )
            );

            /* Update ID/ED and BND related information for the moved vertex */
            INC_DEC(pwgts[to as usize], pwgts[from as usize], vwgt[v as usize]);
            UpdateMovedVertexInfoAndBND(
                v, from, k, to, vrinfo, vnbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
                ii = adjncy[j as usize];
                me = where_[ii as usize];
                myrinfo = graph.ckrinfo + ii;

                oldnnbrs = myrinfo.nnbrs;

                UpdateAdjacentVertexInfoAndBND(
                    ctrl,
                    ii,
                    xadj[(ii + 1) as usize] - xadj[ii as usize],
                    me,
                    from,
                    to,
                    myrinfo,
                    adjwgt[j as usize],
                    nbnd,
                    bndptr,
                    bndind,
                    bndtype,
                );

                assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }
        }

        graph.nbnd = nbnd;

        if (ctrl.dbglvl & METIS_DBG_REFINE) {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}."
              " Nmoves: {:5}, Cut: {:6}, Vol: {:6}\n",
              pwgts[(iargmin(nparts, pwgts,1)) as usize], imax(nparts, pwgts,1),
              ComputeLoadImbalance(graph, nparts, ctrl.pijbm), 
              graph.nbnd, nmoved, graph.mincut, ComputeVolume(graph, where_));
        }

        if (nmoved == 0 || graph.mincut == oldcut) {
            break;
        }
    }

    // WCOREPOP;
}
