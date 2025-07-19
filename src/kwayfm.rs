/*
\file
\brief Routines for k-way refinement

\date Started 7/28/97
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: kwayfm.c 17513 2014-08-05 16:20:50Z dominique $
*/

use crate::*;

/// The vertex is in the queue
const VPQSTATUS_PRESENT: idx_t = 1;
/// The vertex has been extracted from the queue
const VPQSTATUS_EXTRACTED: idx_t = 2;
/// The vertex is not present in the queue and has not been extracted before
const VPQSTATUS_NOTPRESENT: idx_t = 3;

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
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    match ctrl.objtype {
        METIS_OBJTYPE_CUT => {
            if graph.ncon == 1 {
                Greedy_KWayCutOptimize(ctrl, graph, niter, ffactor, omode);
            } else {
                Greedy_McKWayCutOptimize(ctrl, graph, niter, ffactor, omode);
            }
        }

        METIS_OBJTYPE_VOL => {
            if graph.ncon == 1 {
                Greedy_KWayVolOptimize(ctrl, graph, niter, ffactor, omode);
            } else {
                Greedy_McKWayVolOptimize(ctrl, graph, niter, ffactor, omode);
            }
        }

        _ => panic!("Unknown objtype of {}", ctrl.objtype),
    }
}

fn imax(n: usize, slice: &[idx_t], step: usize) -> idx_t {
    debug_assert_eq!(n, slice.len());
    slice.iter().copied().step_by(step).max().unwrap_or_default()
}

fn imin(n: usize, slice: &[idx_t], step: usize) -> idx_t {
    debug_assert_eq!(n, slice.len());
    slice.iter().copied().step_by(step).min().unwrap_or_default()
}

fn isum(n: usize, slice: &[idx_t], step: usize) -> idx_t {
    debug_assert_eq!(n, slice.len());
    slice.iter().copied().step_by(step).sum()
}

/// little helper to replace all the times we get a `ckrinfo_t` + `&[cnbr_t]`
#[inline]
unsafe fn crinfos<'a>(graph_ckrinfo: *const ckrinfo_t, ctrl_cnbrpool: *const cnbr_t, idx: usize) -> (&'a ckrinfo_t, &'a [cnbr_t]) {
    let info = &*graph_ckrinfo.add(idx);
    debug_assert!(info.inbr != -1);
    let slice = std::slice::from_raw_parts(ctrl_cnbrpool.add(info.inbr as usize), info.nnbrs as usize);
    (info, slice)
}

/// little helper to replace all the times we get a `vkrinfo_t` `&[vnbr_t]`
#[inline]
unsafe fn vrinfos<'a>(graph_vkrinfo: *const vkrinfo_t, ctrl_vnbrpool: *const vnbr_t, idx: usize) -> (&'a vkrinfo_t, &'a [vnbr_t]) {
    let info = &*graph_vkrinfo.add(idx);
    debug_assert!(info.inbr != -1);
    let slice = std::slice::from_raw_parts(ctrl_vnbrpool.add(info.inbr as usize), info.nnbrs as usize);
    (info, slice)
}

/// little helper to replace all the times we get a `vkrinfo_t` `&mut [vnbr_t]`. Note that if we
/// want to modify `vkrinfo_t.nnbrs`, then we have to reconstruct in order to access that new
/// element (or to not oob the old in case of removal)
#[inline]
unsafe fn vrinfos_mut<'a>(graph_vkrinfo: *mut vkrinfo_t, ctrl_vnbrpool: *mut vnbr_t, idx: usize) -> (&'a mut vkrinfo_t, &'a mut [vnbr_t]) {
    let info = &mut *graph_vkrinfo.add(idx);
    debug_assert!(info.inbr != -1);
    let slice = std::slice::from_raw_parts_mut(ctrl_vnbrpool.add(info.inbr as usize), info.nnbrs as usize);
    (info, slice)
}

/// function to replace the identically named macro
#[allow(non_snake_case)]
unsafe fn UpdateAdjacentVertexInfoAndBND (ctrl: &mut ctrl_t, vid: usize, adjlen: idx_t, me: usize, from: usize,
    to: usize, myrinfo: &mut ckrinfo_t, ewgt: idx_t, nbnd: &mut usize, bndptr: &mut [idx_t], bndind: &mut [idx_t], bndtype: i32)  {
        // idx_t k;
        // cnbr_t *mynbrs;

        if myrinfo.inbr == -1 {
            myrinfo.inbr = cnbrpoolGetNext(ctrl, adjlen as idx_t);
            myrinfo.nnbrs = 0;
        }
        debug_assert!(debug::CheckRInfo(ctrl, myrinfo) != 0);

        let mynbrs = std::slice::from_raw_parts_mut(
            ctrl.cnbrpool.add(myrinfo.inbr as usize),
            myrinfo.nnbrs as usize + 1,
        );

        /* Update global ID/ED and boundary */
        if me as idx_t == from as idx_t {
            inc_dec!(myrinfo.ed, myrinfo.id, ewgt);
            if bndtype == BNDTYPE_REFINE {
                if myrinfo.ed - myrinfo.id >= 0 && bndptr[vid] == -1 {
                    BNDInsert!(*nbnd, bndind, bndptr, vid);
                }
            } else {
                if myrinfo.ed > 0 && bndptr[vid] == -1 {
                    BNDInsert!(*nbnd, bndind, bndptr, vid);
                }
            }
        } else if me == to {
            inc_dec!(myrinfo.id, myrinfo.ed, ewgt);
            if bndtype == BNDTYPE_REFINE {
                if myrinfo.ed - myrinfo.id < 0 && bndptr[vid] != -1 {
                    BNDDelete!(*nbnd, bndind, bndptr, vid);
                }
            } else {
                if myrinfo.ed <= 0 && bndptr[vid] != -1 {
                    BNDDelete!(*nbnd, bndind, bndptr, vid);
                }
            }
        }

        /* Remove contribution from the .ed of 'from' */
        if me != from {
            for k in 0..myrinfo.nnbrs as usize {
                if mynbrs[k].pid == from as idx_t {
                    if mynbrs[k].ed == (ewgt) {
                        myrinfo.nnbrs -= 1;
                        mynbrs[k] = mynbrs[myrinfo.nnbrs as usize];
                    } else {
                        mynbrs[k].ed -= ewgt;
                    }
                    break;
                }
            }
        }

        /* Add contribution to the .ed of 'to' */
        if me != to {
            let mut k: usize = 0;
            while k < myrinfo.nnbrs as usize {
                if mynbrs[k].pid == to as idx_t {
                    mynbrs[k].ed += ewgt;
                    break;
                }
                k += 1;
            }
            if k == myrinfo.nnbrs as usize {
                mynbrs[k].pid = to as idx_t;
                mynbrs[k].ed = ewgt as idx_t;
                myrinfo.nnbrs += 1;
            }
        }

        debug_assert!(debug::CheckRInfo(ctrl, myrinfo) != 0);
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
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
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

    let bndtype = if omode == OMODE_BALANCE {
        BNDTYPE_REFINE
    } else {
        BNDTYPE_BALANCE
    };

    ffactor = 0.0;
    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    let nparts = ctrl.nparts as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy adjwgt bndind bndptr where_ pwgts ckrinfo vwgt);
    mkslice!(ctrl->tpwgts, nparts);

    /* Setup the weight intervals of the various subdomains */
    let mut minpwgts: Vec<idx_t> = vec![0; nparts as usize];
    let mut maxpwgts: Vec<idx_t> = vec![0; nparts as usize];

    let ubfactor = if omode == OMODE_BALANCE {
        *ctrl.ubfactors
    } else {
        (*ctrl.ubfactors).max(mcutil::ComputeLoadImbalance(graph, nparts as idx_t, ctrl.pijbm))
    };

    for i in (0)..(nparts) {
        maxpwgts[i] = (tpwgts[i] * *graph.tvwgt as f32 * ubfactor) as idx_t;
        minpwgts[i] = (tpwgts[i] * *graph.tvwgt as f32 * (1.0 / ubfactor)) as idx_t;
    }

    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    let mut safetos: Vec<idx_t> = vec![2; nparts as usize];

    let mut nads = &mut [][..];
    let mut adids = &mut [][..];
    let mut adwgts = &mut [][..];
    let mut doms = &mut [][..];
    if ctrl.minconn != 0 {
        minconn::ComputeSubDomainGraph(ctrl, graph);

        nads = std::slice::from_raw_parts_mut(ctrl.nads, nparts);
        adids = std::slice::from_raw_parts_mut(ctrl.adids, nparts);
        adwgts = std::slice::from_raw_parts_mut(ctrl.adwgts, nparts);
        doms = std::slice::from_raw_parts_mut(ctrl.pvec1, nparts);
        doms.fill(0);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    let mut vstatus: Vec<idx_t> = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    let mut updptr: Vec<idx_t> = vec![-1; nvtxs as usize];
    let mut updind: Vec<idx_t> = vec![0; nvtxs as usize];

    let mut bfslvl: Vec<idx_t> = vec![];
    let mut bfsind: Vec<idx_t> = vec![];
    let mut bfsmrk: Vec<idx_t> = vec![];
    if ctrl.contig != 0 {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    }

    if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
        print!("{}: [({:6} {:6}) as usize]-[({:6} {:6}) as usize], Bal: {:5.3}, \
            Nv-Nb[({:6} {:6}) as usize], Cut: {:6}",
            (if omode == OMODE_REFINE { "GRC" } else { "GBC" }),
            pwgts[util::iargmin(pwgts, 1) as usize], imax(nparts, pwgts,1), minpwgts[0], maxpwgts[0], 
            mcutil::ComputeLoadImbalance(graph, nparts as idx_t, ctrl.pijbm), 
            graph.nvtxs, graph.nbnd, graph.mincut);
        if ctrl.minconn != 0 {
            print!(
                ", Doms: [({:3} {:4}) as usize]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    let mut queue = pqueue::RPQueue::new(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        debug_assert!(debug::ComputeCut(graph, where_.as_ptr()) == graph.mincut);
        if omode == OMODE_REFINE {
            debug_assert!(debug::CheckBnd2(graph) != 0);
        }

        if omode == OMODE_BALANCE {
            // /* Check to see if things are out of balance, given the tolerance */
            // for i in (0)..(nparts) {
            //     if (pwgts[i as usize] > maxpwgts[i as usize]
            //         || pwgts[i as usize] < minpwgts[i as usize])
            //     {
            //         break;
            //     }
            // }
            // if (i == nparts)
            // /* Things are balanced. Return right away */
            // {
            //     break;
            // }

            // Check to see if things are out of balance, given the tolerance, and return right
            // away if things are balanced
            if !(0..nparts).any(|i| pwgts[i] > maxpwgts[i]
                || pwgts[i] < minpwgts[i]) {break}
        }

        let oldcut = graph.mincut;
        let nbnd = graph.nbnd as usize;
        let nupd = 0;

        let maxndoms;
        if ctrl.minconn != 0 {
            maxndoms = imax(nparts, nads, 1);
        } else {
            maxndoms = -1;
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(nbnd as idx_t, perm.as_mut_ptr(), nbnd as idx_t / 4, 1);
        for ii in (0)..(nbnd) {
            let i = bndind[perm[ii as usize] as usize] as usize;
            let rgain: real_t = ((if ckrinfo[i].nnbrs > 0 {
                1.0 * ckrinfo[i].ed as f64 / (ckrinfo[i].nnbrs as f64).sqrt()
            } else {
                0.0
            }) - ckrinfo[i as usize].id as f64) as real_t;
            queue.insert(i as idx_t, rgain);
            vstatus[i as usize] = VPQSTATUS_PRESENT;
            ListInsert!(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        let mut nmoved = 0;
        for iii in (0).. {
            // if ((i = rpqGetTop(queue)) == -1) {
            //     break;
            // }
            let Some(i) = queue.pop() else { break };
            let i = i as usize;

            vstatus[i] = VPQSTATUS_EXTRACTED;

            let (myrinfo, mynbrs) = crinfos(graph.ckrinfo, ctrl.cnbrpool, i);

            let from = where_[i as usize] as usize;
            let vwgt = vwgt[i as usize];

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

            if ctrl.contig != 0
            && IsArticulationNode(nvtxs as idx_t, i as idx_t, xadj.as_mut_ptr(), adjncy.as_mut_ptr(), where_.as_mut_ptr(), bfslvl.as_mut_ptr(), bfsind.as_mut_ptr(), bfsmrk.as_mut_ptr()) != 0
            {
                continue;
            }

            if ctrl.minconn != 0 {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            /* Find the most promising subdomain to move to */
            let gain;
            let mut to: usize = usize::MAX;
            let k = if omode == OMODE_REFINE {
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     if (!safetos[to = mynbrs[k as usize].pid]) {
                //         continue;
                //     }
                //     if (((mynbrs[k as usize].ed > myrinfo.id)
                //         && ((pwgts[from as usize] - vwgt >= minpwgts[from as usize])
                //             || (tpwgts[from as usize] * pwgts[to as usize]
                //                 < tpwgts[to as usize] * (pwgts[from as usize] - vwgt)))
                //         && ((pwgts[to as usize] + vwgt <= maxpwgts[to as usize])
                //             || (tpwgts[from as usize] * pwgts[to as usize]
                //                 < tpwgts[to as usize] * (pwgts[from as usize] - vwgt))))
                //         || ((mynbrs[k as usize].ed == myrinfo.id)
                //             && (tpwgts[from as usize] * pwgts[to as usize]
                //                 < tpwgts[to as usize] * (pwgts[from as usize] - vwgt))))
                //     {
                //         break;
                //     }
                // }
                // if (k < 0) {
                //     continue;
                // }

                let Some(mut k) = (0..myrinfo.nnbrs as usize)
                    .filter(|&k| {
                        let to = mynbrs[k].pid as usize;
                        safetos[to] != 0
                    })
                    .rfind(|&k| {
                        let to = mynbrs[k].pid as usize;
                        let from = from as usize;
                        ((mynbrs[k].ed > myrinfo.id)
                        && ((pwgts[from] - vwgt >= minpwgts[from])
                        || (tpwgts[from] * (pwgts[to] as real_t) < tpwgts[to] * (pwgts[from] - vwgt) as real_t))
                        && ((pwgts[to] + vwgt <= maxpwgts[to])
                        || (tpwgts[from] * (pwgts[to] as real_t) < tpwgts[to] * (pwgts[from] - vwgt) as real_t)))
                        || ((mynbrs[k].ed == myrinfo.id)
                        && (tpwgts[from] * (pwgts[to] as real_t) < tpwgts[to] * (pwgts[from] - vwgt) as real_t))
                    })
                else {
                    // break out if you did not find a candidate
                    continue;
                };

                for j in (0..k).rev() {
                    to = mynbrs[j as usize].pid as usize;
                    if safetos[to] == 0 {
                        continue;
                    }
                    if ((mynbrs[j].ed > mynbrs[k].ed)
                    && ((pwgts[from] - vwgt >= minpwgts[from])
                    || ((tpwgts[from] * pwgts[to] as real_t)
                    < tpwgts[to] * (pwgts[from] - vwgt) as real_t))
                    && ((pwgts[to] + vwgt <= maxpwgts[to])
                    || (tpwgts[from] * (pwgts[to] as real_t)
                    < tpwgts[to] * (pwgts[from] - vwgt) as real_t)))
                    || ((mynbrs[j].ed == mynbrs[k].ed)
                    && (tpwgts[mynbrs[k].pid as usize] * (pwgts[to] as real_t)
                    < tpwgts[to] * pwgts[mynbrs[k].pid as usize] as real_t))
                    {
                        k = j;
                    }
                }

                to = mynbrs[k as usize].pid as usize;

                gain = mynbrs[k as usize].ed - myrinfo.id;

                k
            } else {
                /* OMODE_BALANCE */
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     if (!safetos[to = mynbrs[k as usize].pid]) {
                //         continue;
                //     }
                //     /* the correctness of the following test follows from the correctness
                //     of the similar test in the subsequent loop */
                //     if (from >= nparts
                //         || tpwgts[from as usize] * pwgts[to as usize]
                //             < tpwgts[to as usize] * (pwgts[from as usize] - vwgt))
                //     {
                //         break;
                //     }
                // }
                // if (k < 0) {
                //     continue;
                // } /* break out if you did not find a candidate */

                let Some(mut k) = (0..myrinfo.nnbrs as usize)
                    .filter(|&k| {
                        let to = mynbrs[k].pid as usize;
                        safetos[to] != 0
                    })
                    .rfind(|&k| {
                        // the correctness of the following test follows from the correctness of
                        // the similar test in the subsequent loop
                        let to = mynbrs[k].pid as usize;
                        from >= nparts
                            || tpwgts[from] * (pwgts[to] as real_t)
                                < tpwgts[to] * (pwgts[from] - vwgt) as real_t
                    }) else {
                    // break out if you did not find a candidate
                    continue
                };

                for j in (0..k).rev() {
                    let to = mynbrs[j as usize].pid as usize;
                    if safetos[to] == 0 {
                        continue;
                    }
                    if tpwgts[mynbrs[k].pid as usize] * (pwgts[to] as real_t)
                    < tpwgts[to] * pwgts[mynbrs[k].pid as usize] as real_t
                    {
                        k = j;
                    }
                }

                to = mynbrs[k as usize].pid as usize;

                //if (pwgts[from as usize] < maxpwgts[from as usize] && pwgts[to as usize] > minpwgts[to as usize] &&
                //    mynbrs[k as usize].ed-myrinfo.id < 0)
                //  continue;

                k
            };

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
            if ctrl.minconn != 0 {
                /* take care of i's move itself */
                minconn::UpdateEdgeSubDomainGraph(
                    ctrl,
                    from as idx_t,
                    to as idx_t,
                    myrinfo.id - mynbrs[k as usize].ed,
                    &mut maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    let me = where_[adjncy[j as usize] as usize] as usize;
                    if me != from && me != to {
                        minconn::UpdateEdgeSubDomainGraph(ctrl, from as idx_t, me as idx_t, -adjwgt[j as usize], &mut maxndoms);
                        minconn::UpdateEdgeSubDomainGraph(ctrl, to as idx_t, me as idx_t, adjwgt[j as usize], &mut maxndoms);
                    }
                }
            }

            /* Update ID/ED and BND related information for the moved vertex */
            inc_dec!(pwgts[to as usize], pwgts[from as usize], vwgt);
            UpdateMovedVertexInfoAndBND!(
                i, from, k, to, myrinfo, mynbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                let ii = adjncy[j as usize] as usize;
                let me = where_[ii as usize] as usize;
                let myrinfo = &mut *graph.ckrinfo.add(ii);

                let oldnnbrs = myrinfo.nnbrs;

                UpdateAdjacentVertexInfoAndBND(
                    ctrl,
                    ii,
                    (xadj[ii + 1] - xadj[ii]) as usize,
                    me,
                    from,
                    to,
                    myrinfo,
                    adjwgt[j as usize],
                    &mut nbnd,
                    &mut bndptr,
                    &mut bndind,
                    bndtype,
                );

                UpdateQueueInfo(
                    queue, vstatus, ii, me, from, to, myrinfo, oldnnbrs, nupd, updptr, updind,
                    bndtype,
                );

                debug_assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }
        }

        graph.nbnd = nbnd as idx_t;

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            debug_assert!(updptr[updind[i as usize] as usize] != -1);
            debug_assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
            print!("\t[{:6} {:6}], Bal: {:5.3}, Nb: {:6}. \
              Nmoves: {:5}, Cut: {:6}, Vol: {:6}",
              pwgts[util::iargmin(pwgts, 1)], imax(nparts, pwgts,1),
              mcutil::ComputeLoadImbalance(graph, nparts as idx_t, ctrl.pijbm), 
              graph.nbnd, nmoved, graph.mincut, debug::ComputeVolume(graph, where_.as_ptr()));
            if ctrl.minconn != 0 {
                print!(
                    ", Doms: [{:3} {:4}]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if nmoved == 0 || (omode == OMODE_REFINE && graph.mincut == oldcut) {
            break;
        }
    }

    // rpqDestroy(queue);

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
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
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

    let bndtype = if omode == OMODE_BALANCE {
        BNDTYPE_REFINE
    } else {
        BNDTYPE_BALANCE
    };

    // WCOREPUSH;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy bndptr bndind where_ pwgts vwgt tvwgt vsize);

    let nparts = ctrl.nparts as usize;
    mkslice!(ctrl->tpwgts, nparts);

    /* Setup the weight intervals of the various subdomains */
    let mut minpwgts: Vec<idx_t> = vec![0; nparts as usize];
    let mut maxpwgts: Vec<idx_t> = vec![0; nparts as usize];

    mkslice!(ctrl->tpwgts, nparts);
    for i in (0)..(nparts) {
        maxpwgts[i as usize] =
            (tpwgts[i as usize] * tvwgt[0] as real_t * *ctrl.ubfactors) as idx_t;
        minpwgts[i as usize] =
            (tpwgts[i as usize] * tvwgt[0] as real_t * (1.0 / *ctrl.ubfactors)) as idx_t;
    }

    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    let mut safetos: Vec<idx_t> = vec![2; nparts as usize];

    let mut nads = &mut [][..];
    let mut adids = &mut [][..];
    let mut adwgts = &mut [][..];
    let mut doms = &mut [][..];
    if ctrl.minconn != 0 {
        minconn::ComputeSubDomainGraph(ctrl, graph);

        nads = std::slice::from_raw_parts_mut(ctrl.nads, nparts);
        adids = std::slice::from_raw_parts_mut(ctrl.adids, nparts);
        adwgts = std::slice::from_raw_parts_mut(ctrl.adwgts, nparts);
        doms = std::slice::from_raw_parts_mut(ctrl.pvec1, nparts);
        doms.fill(0);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    let mut vstatus: Vec<idx_t> = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    let mut updptr: Vec<idx_t> = vec![-1; nvtxs as usize];
    let mut updind: Vec<idx_t> = vec![0; nvtxs as usize];

    let mut bfslvl: Vec<idx_t>;
    let mut bfsind: Vec<idx_t>;
    let mut bfsmrk: Vec<idx_t>;
    if ctrl.contig != 0 {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    } else {
        bfslvl = vec![];
        bfsind = vec![];
        bfsmrk = vec![];
    }

    /* Vol-refinement specific working arrays */
    let mut modind: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut vmarker: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut pmarker: Vec<idx_t> = vec![-1; nparts as usize];

    if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
        print!("{}: [{:6} {:6}]-[{:6} {:6}], Bal: {:5.3}\
            , Nv-Nb[{:6} {:6}], Cut: {:5}, Vol: {:5}",
            (if omode == OMODE_REFINE { "GRV" } else { "GBV" }),
            pwgts[(util::iargmin(&pwgts, 1)) as usize], imax(nparts, pwgts,1), minpwgts[0], maxpwgts[0], 
            mcutil::ComputeLoadImbalance(graph, nparts as idx_t, ctrl.pijbm), 
            graph.nvtxs, graph.nbnd, graph.mincut, graph.minvol);
        if ctrl.minconn != 0 {
            print!(
                ", Doms: [({:3} {:4}) as usize]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    // queue = ipqCreate(nvtxs);
    let mut queue = pqueue::IPQueue::new(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        debug_assert!(debug::ComputeVolume(graph, where_.as_ptr()) == graph.minvol);

        if omode == OMODE_BALANCE {
            // /* Check to see if things are out of balance, given the tolerance */
            // for i in (0)..(nparts) {
            //     if (pwgts[i as usize] > maxpwgts[i as usize]) {
            //         break;
            //     }
            // }
            // if (i == nparts)
            // /* Things are balanced. Return right away */
            // {
            //     break;
            // }

            // Check to see if things are out of balance, given the tolerance
            if pwgts.iter().enumerate().all(|(i, &p)| p <= maxpwgts[i]) {
                // Things are balanced. Return right away
                break
            }
        }

        let oldcut = graph.mincut;
        let oldvol = graph.minvol;
        let mut nupd = 0;

        let mut maxndoms;
        if ctrl.minconn != 0 {
            maxndoms = imax(nparts, nads, 1);
        } else {
            maxndoms = -1;
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(graph.nbnd, perm.as_mut_ptr(), graph.nbnd / 4, 1);
        for ii in (0)..(graph.nbnd) {
            let i = bndind[perm[ii as usize] as usize] as usize;
            queue.insert(i as idx_t, (*graph.vkrinfo.add(i)).gv);
            vstatus[i] = VPQSTATUS_PRESENT;
            ListInsert!(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        let mut nmoved = 0;
        let mut to: usize;
        for iii in (0).. {
            let Some(i) = queue.pop() else { break };
            let i = i as usize;

            vstatus[i] = VPQSTATUS_EXTRACTED;

            let (myrinfo, mynbrs) = vrinfos(graph.vkrinfo, ctrl.vnbrpool, i);

            let from = where_[i as usize] as usize;
            let vwgt = vwgt[i as usize];

            /* Prevent moves that make 'from' domain underbalanced */
            if omode == OMODE_REFINE {
                if myrinfo.nid > 0 && pwgts[from as usize] - vwgt < minpwgts[from as usize] {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                if pwgts[from as usize] - vwgt < minpwgts[from as usize] {
                    continue;
                }
            }

            if ctrl.contig != 0 && IsArticulationNode(nvtxs as idx_t, i as idx_t, xadj.as_ptr(), adjncy.as_ptr(), where_.as_ptr(), bfslvl.as_mut_ptr(), bfsind.as_mut_ptr(), bfsmrk.as_mut_ptr()) != 0
            {
                continue;
            }

            if ctrl.minconn != 0 {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            let xgain = if myrinfo.nid == 0 && myrinfo.ned > 0 {
                vsize[i as usize]
            } else {
                0
            };

            /* Find the most promising subdomain to move to */
            let mut to;
            let k = if omode == OMODE_REFINE {
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     to = mynbrs[k as usize].pid as usize;
                //     if (!safetos[to]) {
                //         continue;
                //     }
                //     let gain = mynbrs[k as usize].gv + xgain;
                //     if (gain >= 0
                //         && pwgts[to as usize] + vwgt <= maxpwgts[to as usize] + ffactor * gain)
                //     {
                //         break;
                //     }
                // }

                // if (k < 0) {
                //     continue;
                // } /* break out if you did not find a candidate */

                let Some(mut k) = (0..myrinfo.nnbrs as usize).filter(|&k| {
                    let to = mynbrs[k].pid as usize;
                    safetos[to] != 0
                }).rfind(|&k| {
                        to = mynbrs[k].pid as usize;
                        let gain = mynbrs[k].gv + xgain;
                        
                        gain >= 0 
                            && (pwgts[to as usize] + vwgt) as real_t <= maxpwgts[to as usize] as real_t + ffactor * gain as real_t
                    }) else { continue };

                for j in (0..k).rev() {
                    let to = mynbrs[j as usize].pid as usize;
                    if safetos[to] == 0 {
                        continue;
                    }
                    let gain = mynbrs[j].gv + xgain;
                    if (mynbrs[j].gv > mynbrs[k].gv
                    && (pwgts[to] + vwgt) as real_t <= maxpwgts[to] as real_t + ffactor * gain as real_t)
                    || (mynbrs[j].gv == mynbrs[k].gv
                    && mynbrs[j].ned > mynbrs[k].ned
                    && pwgts[to] + vwgt <= maxpwgts[to])
                    || (mynbrs[j].gv == mynbrs[k].gv
                    && mynbrs[j].ned == mynbrs[k].ned
                    && tpwgts[mynbrs[k].pid as usize] * (pwgts[to] as real_t)
                    < tpwgts[to] * pwgts[mynbrs[k].pid as usize] as real_t)
                    {
                        k = j;
                    }
                }
                to = mynbrs[k as usize].pid as usize;

                debug_assert!(xgain + mynbrs[k as usize].gv >= 0);

                let mut j = 0;
                if xgain + mynbrs[k as usize].gv > 0 || mynbrs[k as usize].ned - myrinfo.nid > 0 {
                    j = 1;
                } else if mynbrs[k as usize].ned - myrinfo.nid == 0 {
                    if (iii % 2 == 0 && safetos[to as usize] == 2)
                    || pwgts[from as usize] >= maxpwgts[from as usize]
                    || tpwgts[from as usize] * ((pwgts[to as usize] + vwgt) as real_t)
                    < tpwgts[to as usize] * pwgts[from as usize] as real_t
                    {
                        j = 1;
                    }
                }
                if j == 0 {
                    continue;
                }

                k
            } else {
                /* OMODE_BALANCE */
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     if !safetos[to = mynbrs[k as usize].pid] {
                //         continue;
                //     }
                //     if pwgts[to as usize] + vwgt <= maxpwgts[to as usize]
                //     || tpwgts[from as usize] * (pwgts[to as usize] + vwgt)
                //     <= tpwgts[to as usize] * pwgts[from as usize]
                //     {
                //         break;
                //     }
                // }
                // if k < 0 {
                //     continue;
                // } /* break out if you did not find a candidate */

                let Some(mut k) = (0..myrinfo.nnbrs as usize).filter(|&k| {
                    safetos[mynbrs[k].pid as usize] != 0
                }).rfind(|&k| {
                        let to = mynbrs[k].pid as usize;
                        pwgts[to as usize] + vwgt <= maxpwgts[to as usize]
                        || tpwgts[from as usize] * (pwgts[to as usize] + vwgt) as real_t
                        <= tpwgts[to as usize] * pwgts[from as usize] as real_t
                    }) else { continue };

                for j in (0..k).rev() {
                    let to = mynbrs[j as usize].pid as usize;
                    if safetos[to] == 0 {
                        continue;
                    }
                    if tpwgts[mynbrs[k].pid as usize] * (pwgts[to] as real_t)
                    < tpwgts[to] * pwgts[mynbrs[k].pid as usize] as real_t
                    {
                        k = j;
                    }
                }
                to = mynbrs[k as usize].pid as usize;

                if pwgts[from as usize] < maxpwgts[from as usize]
                && pwgts[to as usize] > minpwgts[to as usize]
                && (xgain + mynbrs[k as usize].gv < 0
                || (xgain + mynbrs[k as usize].gv == 0
                && mynbrs[k as usize].ned - myrinfo.nid < 0))
                {
                    continue;
                }

                k
            };

            /*=====================================================================
             * If we got here, we can now move the vertex from 'from' to 'to'
             *======================================================================*/
            inc_dec!(pwgts[to as usize], pwgts[from as usize], vwgt);
            graph.mincut -= mynbrs[k as usize].ned - myrinfo.nid;
            graph.minvol -= xgain + mynbrs[k as usize].gv;
            where_[i as usize] = to as idx_t;
            nmoved += 1;

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_MOVEINFO,
                print!("\t\tMoving {:6} from {:3} to {:3}. \
                Gain: [{:4} {:4}]. Cut: {:6}, Vol: {:6}\n", 
                    i, from, to, xgain+mynbrs[k as usize].gv, mynbrs[k as usize].ned-myrinfo.nid, 
              graph.mincut, graph.minvol)
            );

            /* Update the subdomain connectivity information */
            if ctrl.minconn != 0 {
                /* take care of i's move itself */
                minconn::UpdateEdgeSubDomainGraph(
                    ctrl,
                    from as idx_t,
                    to as idx_t,
                    myrinfo.nid - mynbrs[k as usize].ned,
                    &mut maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    let me = where_[adjncy[j as usize] as usize] as usize;
                    if me != from && me != to {
                        minconn::UpdateEdgeSubDomainGraph(ctrl, from as idx_t, me as idx_t, -1, &mut maxndoms);
                        minconn::UpdateEdgeSubDomainGraph(ctrl, to as idx_t, me as idx_t, 1, &mut maxndoms);
                    }
                }
            }

            /* Update the id/ed/gains/bnd/queue of potentially affected nodes */
            KWayVolUpdate(
                ctrl, graph, i, from, to, &mut queue, &mut vstatus, &mut nupd, &mut updptr, &mut updind, bndtype, &mut vmarker,
                &mut pmarker, &mut modind,
            );

            /*CheckKWayVolPartitionParams(ctrl, graph); */
        }

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            debug_assert!(updptr[updind[i as usize] as usize] != -1);
            debug_assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
            print!("\t[{:6} {:6}], Bal: {:5.3}, Nb: {:6}. \
              Nmoves: {:5}, Cut: {:6}, Vol: {:6}",
              pwgts[(iargmin(nparts, pwgts,1)) as usize], imax(nparts, pwgts,1),
              mcutil::ComputeLoadImbalance(graph, nparts, ctrl.pijbm), 
              graph.nbnd, nmoved, graph.mincut, graph.minvol);
            if ctrl.minconn != 0 {
                print!(
                    ", Doms: [{:3} {:4}]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if nmoved == 0
        || (omode == OMODE_REFINE && graph.minvol == oldvol && graph.mincut == oldcut)
        {
            break;
        }
    }


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
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
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

    let bndtype = if omode == OMODE_BALANCE {
        BNDTYPE_REFINE
    } else {
        BNDTYPE_BALANCE
    };


    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    let nparts = ctrl.nparts as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy vwgt adjwgt bndind bndptr where_ pwgts ckrinfo);

    mkslice_mut!(ctrl->pijbm, nparts * ncon);

    /* Determine the ubfactors. The method used is different based on omode.
    When OMODE_BALANCE, the ubfactors are those supplied by the user.
    When OMODE_REFINE, the ubfactors are the max of the current partition
    and the user-specified ones. */
    let mut ubfactors: Vec<real_t> = vec![0.0; ncon as usize];
    mcutil::ComputeLoadImbalanceVec(graph, nparts as idx_t, pijbm.as_ptr(), ubfactors.as_mut_ptr());
    let origbal = mcutil::rvecmaxdiff(ncon as idx_t, ubfactors.as_mut_ptr(), ctrl.ubfactors);

    {
        let ctrl_ubfactors = std::slice::from_raw_parts(ctrl.ubfactors, ncon);
        if omode == OMODE_BALANCE {
            ubfactors.copy_from_slice(ctrl_ubfactors);
        } else {
            for i in (0)..(ncon) {
                ubfactors[i as usize] = if ubfactors[i as usize] > ctrl_ubfactors[i as usize] {
                    ubfactors[i as usize]
                } else {
                        ctrl_ubfactors[i as usize]
                    };
            }
        }
    }

    /* Setup the weight intervals of the various subdomains */
    let mut minpwgts: Vec<idx_t> = vec![0; nparts * ncon as usize];
    let mut maxpwgts: Vec<idx_t> = vec![0; nparts * ncon as usize];

    for i in (0)..(nparts) {
        mkslice!(ctrl->tpwgts, nparts * ncon);
        mkslice!(graph->tvwgt, ncon);
        for j in (0)..(ncon) {
            maxpwgts[(i * ncon + j) as usize] = (tpwgts[(i * ncon + j) as usize]
                * tvwgt[j as usize] as real_t
                * ubfactors[j as usize]) as idx_t;
            /*minpwgts[(i*ncon+j) as usize]  = ctrl.tpwgts[(i*ncon+j) as usize]*graph.tvwgt[j as usize]*(.9/ubfactors[j as usize]);*/
            minpwgts[(i * ncon + j) as usize] =
                (tpwgts[(i * ncon + j) as usize] * tvwgt[j as usize] as real_t * 0.2) as idx_t;
        }
    }

    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    let mut safetos: Vec<idx_t> = vec![2; nparts as usize];

    let mut nads = &mut [][..];
    let mut adids = &mut [][..];
    let mut adwgts = &mut [][..];
    let mut doms = &mut [][..];
    if ctrl.minconn != 0 {
        minconn::ComputeSubDomainGraph(ctrl, graph);

        nads = std::slice::from_raw_parts_mut(ctrl.nads, nparts);
        adids = std::slice::from_raw_parts_mut(ctrl.adids, nparts);
        adwgts = std::slice::from_raw_parts_mut(ctrl.adwgts, nparts);
        doms = std::slice::from_raw_parts_mut(ctrl.pvec1, nparts);
        doms.fill(0);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    let mut vstatus: Vec<idx_t> = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    let mut updptr: Vec<idx_t> = vec![-1; nvtxs as usize];
    let mut updind: Vec<idx_t> = vec![0; nvtxs as usize];

    let mut bfslvl: Vec<idx_t>;
    let mut bfsind: Vec<idx_t>;
    let mut bfsmrk: Vec<idx_t>;
    if ctrl.contig != 0 {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    } else {
        bfslvl = vec![];
        bfsind = vec![];
        bfsmrk = vec![];
    }

    if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
        print!("{}: [{:6} {:6} {:6}], Bal: {:5.3}({:.3}), \
            Nv-Nb[{:6} {:6}], Cut: {:6}, ({:})",
            (if omode == OMODE_REFINE { "GRC" } else { "GBC" }),
            imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), imax(nparts*ncon, &maxpwgts, 1),
            mcutil::ComputeLoadImbalance(graph, nparts as idx_t, pijbm), origbal,
            graph.nvtxs, graph.nbnd, graph.mincut, niter);
        if ctrl.minconn != 0 {
            print!(
                ", Doms: [({:3} {:4}) as usize]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    let mut queue = pqueue::RPQueue::new(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        debug_assert!(debug::ComputeCut(graph, where_.as_ptr()) == graph.mincut);
        if omode == OMODE_REFINE {
            debug_assert!(debug::CheckBnd2(graph) != 0);
        }

        /* In balancing mode, exit as soon as balance is reached */
        if omode == OMODE_BALANCE && kwayrefine::IsBalanced(ctrl, graph, 0.0) != 0 {
            break;
        }

        let oldcut = graph.mincut;
        let nbnd = graph.nbnd as usize;
        let mut nupd = 0;

        let mut maxndoms;
        if ctrl.minconn != 0 {
            maxndoms = imax(nparts, nads, 1);
        } else {
            maxndoms = -1;
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(nbnd as idx_t, perm.as_mut_ptr(), nbnd as idx_t / 4, 1);
        for ii in (0)..(nbnd) {
            let i = bndind[perm[ii] as usize] as usize;
            let rgain = (if ckrinfo[i].nnbrs > 0 {
                1.0 * ckrinfo[i].ed as f64 / (ckrinfo[i].nnbrs as f64).sqrt()
            } else {
                0.0
            }) - ckrinfo[i].id as f64;
            queue.insert(i as idx_t, rgain as real_t);
            vstatus[i] = VPQSTATUS_PRESENT;
            ListInsert!(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        let mut nmoved = 0;
        for iii in (0).. {
            let Some(i) = queue.pop() else { break };
            let i = i as usize;
            vstatus[i as usize] = VPQSTATUS_EXTRACTED;

            let (myrinfo, mynbrs) = crinfos(graph.ckrinfo, ctrl.cnbrpool, i);

            let from = where_[i as usize] as usize;

            /* Prevent moves that make 'from' domain underbalanced */
            if omode == OMODE_REFINE {
                if myrinfo.id > 0
                && !util::ivecaxpygez(
                    -1,
                    &vwgt[cntrng!(i * ncon, ncon)],
                    &pwgts[cntrng!(from * ncon, ncon)],
                    &minpwgts[cntrng!(from * ncon, ncon)],
                )
                {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                if !util::ivecaxpygez(
                    -1,
                    &vwgt[cntrng!(i * ncon, ncon)],
                    &pwgts[cntrng!(from * ncon, ncon)],
                    &minpwgts[cntrng!(from * ncon, ncon)],
                ) {
                    continue;
                }
            }

            if ctrl.contig != 0 && IsArticulationNode(nvtxs as idx_t, i as idx_t, xadj.as_ptr(), adjncy.as_ptr(), where_.as_ptr(), bfslvl.as_mut_ptr(), bfsind.as_mut_ptr(), bfsmrk.as_mut_ptr()) != 0
            {
                continue;
            }

            if ctrl.minconn != 0 {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            /* Find the most promising subdomain to move to */
            let mut to = -1;
            let k = if (omode == OMODE_REFINE) {
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     if (!safetos[to = mynbrs[k as usize].pid]) {
                //         continue;
                //     }
                //     gain = mynbrs[k as usize].ed - myrinfo.id;
                //     if (gain >= 0
                //         && ivecaxpylez(
                //             ncon,
                //             1,
                //             vwgt + i * ncon,
                //             pwgts + to * ncon,
                //             maxpwgts + to * ncon,
                //         ))
                //     {
                //         break;
                //     }
                // }
                // if (k < 0) {
                //     continue;
                // } /* break out if you did not find a candidate */

                let Some(mut k) = (0..myrinfo.nnbrs as usize).filter(|k| {
                    to = mynbrs[k as usize].pid as usize;
                    safetos[to] != 0
                }).rfind(|&k| {
                    to = mynbrs[k as usize].pid as usize;
                    let gain = mynbrs[k as usize].ed - myrinfo.id;
                    gain >= 0
                    && util::ivecaxpylez(
                        1,
                        &vwgt[cntrng!(i * ncon, ncon)],
                        &pwgts[cntrng!(to * ncon, ncon)],
                        &maxpwgts[cntrng!(to * ncon, ncon)],
                    )
                    }) else {
                    // break out if you did not find a candidate
                    continue
                };

                let mut cto = to;
                for j in (0..k).rev() {
                    let to = mynbrs[j as usize].pid as usize;
                    if safetos[to] == 0 {
                        continue;
                    }
                    if (mynbrs[j].ed > mynbrs[k].ed
                    && util::ivecaxpylez(
                        1,
                        &vwgt[cntrng!(i * ncon, ncon)],
                        &pwgts[cntrng!(to * ncon, ncon)],
                        &maxpwgts[cntrng!(to * ncon, ncon)],
                    ))
                    || (mynbrs[j as usize].ed == mynbrs[k as usize].ed
                    && mcutil::BetterBalanceKWay(
                        ncon as idx_t,
                        vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                        ubfactors.as_ptr(),
                        1,
                        pwgts[cntrng!(cto * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(cto * ncon, ncon)].as_ptr(),
                        1,
                        pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                    ) != 0)
                    {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                let gain = mynbrs[k as usize].ed - myrinfo.id;
                if !(gain > 0
                || (gain == 0
                && (mcutil::BetterBalanceKWay(
                    ncon,
                    vwgt + i * ncon,
                    ubfactors,
                    -1,
                    pwgts + from * ncon,
                    pijbm + from * ncon,
                    1,
                    pwgts + to * ncon,
                    pijbm + to * ncon,
                ) != 0 || (iii % 2 == 0 && safetos[to as usize] == 2))))
                {
                    continue;
                }

                k
            } else {
                /* OMODE_BALANCE */
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     if (!safetos[to = mynbrs[k as usize].pid]) {
                //         continue;
                //     }
                //     if (ivecaxpylez(
                //         ncon,
                //         1,
                //         vwgt + i * ncon,
                //         pwgts + to * ncon,
                //         maxpwgts + to * ncon,
                //     ) || BetterBalanceKWay(
                //         ncon,
                //         vwgt + i * ncon,
                //         ubfactors,
                //         -1,
                //         pwgts + from * ncon,
                //         pijbm + from * ncon,
                //         1,
                //         pwgts + to * ncon,
                //         pijbm + to * ncon,
                //     )) {
                //         break;
                //     }
                // }
                // if (k < 0) {
                //     continue;
                // } /* break out if you did not find a candidate */

                let Some(k) = (0..myrinfo.nnbrs as usize).filter(|k| {
                    to = mynbrs[k].pid as usize;
                    safetos[to] != 0
                }).rfind(|k| {
                    to = mynbrs[k].pid as usize;
                        util::ivecaxpylez(
                            1,
                            &vwgt[cntrng!(i * ncon, ncon)],
                            &pwgts[cntrng!(to * ncon, ncon)],
                            &maxpwgts[cntrng!(to * ncon, ncon)],
                        ) || mcutil::BetterBalanceKWay(
                            ncon as idx_t,
                            vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                            ubfactors,
                            -1,
                            pwgts[cntrng!(from * ncon, ncon)].as_ptr(),
                            pijbm[cntrng!(from * ncon, ncon)].as_ptr(),
                            1,
                            pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                            pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                        ) != 0
                    }
                    ) else {
                    // break out if you did not find a candidate
                    continue
                };

                let mut cto = to;
                for j in (0..k).rev() {
                    let to = mynbrs[j as usize].pid;
                    if safetos[to] == 0 {
                        continue;
                    }
                    if mcutil::BetterBalanceKWay(
                        ncon as idx_t,
                        vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                        ubfactors,
                        1,
                        pwgts[cntrng!(cto * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(cto * ncon, ncon)].as_ptr(),
                        1,
                        pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                    ) != 0 {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                if mynbrs[k as usize].ed - myrinfo.id < 0
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
                )
                {
                    continue;
                }

                k
            };

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
            if ctrl.minconn {
                /* take care of i's move itself */
                minconn::UpdateEdgeSubDomainGraph(
                    ctrl,
                    from,
                    to,
                    myrinfo.id - mynbrs[k as usize].ed,
                    &mut maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    let me = where_[adjncy[j as usize] as usize] as usize;
                    if me != from && me != to {
                        minconn::UpdateEdgeSubDomainGraph(ctrl, from as idx_t, me as idx_t, -adjwgt[j as usize], &mut maxndoms);
                        minconn::UpdateEdgeSubDomainGraph(ctrl, to as idx_t, me as idx_t, adjwgt[j as usize], &mut maxndoms);
                    }
                }
            }

            /* Update ID/ED and BND related information for the moved vertex */
            blas::iaxpy(ncon, 1, vwgt + i * ncon, 1, pwgts + to * ncon, 1);
            blas::iaxpy(ncon, -1, vwgt + i * ncon, 1, pwgts + from * ncon, 1);
            UpdateMovedVertexInfoAndBND!(
                i, from, k, to, myrinfo, mynbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                let ii = adjncy[j as usize];
                let me = where_[ii as usize] as usize;
                myrinfo = graph.ckrinfo + ii;

                let oldnnbrs = myrinfo.nnbrs;

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

                debug_assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }
        }

        graph.nbnd = nbnd;

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            debug_assert!(updptr[updind[i as usize] as usize] != -1);
            debug_assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}. \
              Nmoves: {:5}, Cut: {:6}, Vol: {:6}",
              imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), 
              ComputeLoadImbalance(graph, nparts, pijbm), 
              graph.nbnd, nmoved, graph.mincut, ComputeVolume(graph, where_));
            if ctrl.minconn {
                print!(
                    ", Doms: [({:3} {:4}) as usize]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if nmoved == 0 || (omode == OMODE_REFINE && graph.mincut == oldcut) {
            break;
        }
    }

    // rpqDestroy(queue);

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
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
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

    let bndtype = if omode == OMODE_BALANCE {
        BNDTYPE_REFINE
    } else {
        BNDTYPE_BALANCE
    };


    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    get_graph_slices!(ctrl, graph => xadj vwgt adjncy bndptr bndind where_ pwgts vsize);

    let nparts = ctrl.nparts as usize;

    mkslice!(ctrl->pijbm, nparts * ncon);

    /* Determine the ubfactors. The method used is different based on omode.
    When OMODE_BALANCE, the ubfactors are those supplied by the user.
    When OMODE_REFINE, the ubfactors are the max of the current partition
    and the user-specified ones. */
    let mut ubfactors: Vec<real_t> = vec![0.0; ncon as usize];
    mcutil::ComputeLoadImbalanceVec(graph, nparts as idx_t, pijbm.as_ptr(), ubfactors.as_mut_ptr());
    let origbal = mcutil::rvecmaxdiff(ncon as idx_t, ubfactors.as_mut_ptr(), ctrl.ubfactors);
    {
        mkslice!(ctrl_ubfactors: ctrl->ubfactors, ncon);
        if omode == OMODE_BALANCE {
            ubfactors.copy_from_slice(ctrl_ubfactors);
        } else {
            for i in (0)..(ncon) {
                ubfactors[i as usize] = if ubfactors[i as usize] > ctrl_ubfactors[i as usize] {
                    ubfactors[i as usize]
                } else {
                        ctrl_ubfactors[i as usize]
                    };
            }
        }
    }

    /* Setup the weight intervals of the various subdomains */
    let mut minpwgts: Vec<idx_t> = vec![0; nparts * ncon as usize];
    let mut maxpwgts: Vec<idx_t> = vec![0; nparts * ncon as usize];

    for i in (0)..(nparts) {
        mkslice!(graph->tvwgt, ncon);
        mkslice!(ctrl->tpwgts, ncon * nparts);
        for j in (0)..(ncon) {
            maxpwgts[(i * ncon + j) as usize] = tpwgts[(i * ncon + j) as usize]
                * tvwgt[j as usize]
                * ubfactors[j as usize];
            /*minpwgts[(i*ncon+j) as usize]  = ctrl.tpwgts[(i*ncon+j) as usize]*graph.tvwgt[j as usize]*(.9/ubfactors[j as usize]); */
            minpwgts[(i * ncon + j) as usize] =
                tpwgts[(i * ncon + j) as usize] * tvwgt[j as usize] * 0.2;
        }
    }

    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize];

    /* This stores the valid target subdomains. It is used when ctrl.minconn to
    control the subdomains to which moves are allowed to be made.
    When ctrl.minconn is false, the default values of 2 allow all moves to
    go through and it does not interfere with the zero-gain move selection. */
    let mut safetos: Vec<idx_t> = vec![2; nparts as usize];

    let mut nads = &mut [][..];
    let mut adids = &mut [][..];
    let mut adwgts = &mut [][..];
    let mut doms = &mut [][..];
    if ctrl.minconn != 0 {
        minconn::ComputeSubDomainGraph(ctrl, graph);

        nads = std::slice::from_raw_parts_mut(ctrl.nads, nparts);
        adids = std::slice::from_raw_parts_mut(ctrl.adids, nparts);
        adwgts = std::slice::from_raw_parts_mut(ctrl.adwgts, nparts);
        doms = std::slice::from_raw_parts_mut(ctrl.pvec1, nparts);
        doms.fill(0);
    }

    /* Setup updptr, updind like boundary info to keep track of the vertices whose
    vstatus's need to be reset at the end of the inner iteration */
    let mut vstatus: Vec<idx_t> = vec![VPQSTATUS_NOTPRESENT; nvtxs as usize];
    let mut updptr: Vec<idx_t> = vec![-1; nvtxs as usize];
    let mut updind: Vec<idx_t> = vec![0; nvtxs as usize];

    let mut bfslvl: Vec<idx_t>;
    let mut bfsind: Vec<idx_t>;
    let mut bfsmrk: Vec<idx_t>;
    if ctrl.contig != 0 {
        /* The arrays that will be used for limited check of articulation points */
        bfslvl = vec![0; nvtxs as usize];
        bfsind = vec![0; nvtxs as usize];
        bfsmrk = vec![0; nvtxs as usize];
    } else {
        bfslvl = vec![];
        bfsind = vec![];
        bfsmrk = vec![];
    }

    /* Vol-refinement specific working arrays */
    let mut modind: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut vmarker: Vec<idx_t> = vec![0; nvtxs as usize];
    let mut pmarker: Vec<idx_t> = vec![-1; nparts as usize];

    if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
        print!("{}: [{:6} {:6} {:6}], Bal: {:5.3}({:.3}), \
            Nv-Nb[{:6} {:6}], Cut: {:5}, Vol: {:5}, ({:})",
            (if omode == OMODE_REFINE { "GRV" } else { "GBV" }),
            imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), imax(nparts*ncon, &maxpwgts, 1),
            mcutil::ComputeLoadImbalance(graph, nparts as idx_t, pijbm.as_ptr()), origbal,
            graph.nvtxs, graph.nbnd, graph.mincut, graph.minvol, niter);
        if ctrl.minconn != 0 {
            print!(
                ", Doms: [{:3} {:4}]",
                imax(nparts, nads, 1),
                isum(nparts, nads, 1)
            );
        }
        print!("\n");
    }

    let mut queue = pqueue::IPQueue::new(nvtxs);

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        debug_assert!(debug::ComputeVolume(graph, where_.as_ptr()) == graph.minvol);

        /* In balancing mode, exit as soon as balance is reached */
        if omode == OMODE_BALANCE && kwayrefine::IsBalanced(ctrl, graph, 0.0) != 0 {
            break;
        }

        let oldcut = graph.mincut;
        let oldvol = graph.minvol;
        let mut nupd = 0;

        let maxndoms;
        if ctrl.minconn != 0 {
            maxndoms = imax(nparts, nads, 1);
        } else {
            maxndoms = idx_t::MAX;
        }

        /* Insert the boundary vertices in the priority queue */
        irandArrayPermute(graph.nbnd, perm.as_mut_ptr(), graph.nbnd / 4, 1);
        for ii in (0)..(graph.nbnd) {
            let i = bndind[perm[ii as usize] as usize] as usize;
            get_graph_slices!(graph => vkrinfo);
            queue.insert(i as idx_t, vkrinfo[i].gv);
            vstatus[i as usize] = VPQSTATUS_PRESENT;
            ListInsert!(nupd, updind, updptr, i);
        }

        /* Start extracting vertices from the queue and try to move them */
        let mut nmoved = 0;
        for iii in (0).. {
            let Some(i) = queue.pop() else { break };
            let i = i as usize;

            vstatus[i] = VPQSTATUS_EXTRACTED;

            let myrinfo = &*graph.vkrinfo.add(i);
            let mynbrs = std::slice::from_raw_parts(ctrl.vnbrpool.add(myrinfo.inbr as usize), myrinfo.nnbrs as usize);

            let from = where_[i as usize] as usize;

            /* Prevent moves that make 'from' domain underbalanced */
            if omode == OMODE_REFINE {
                if myrinfo.nid > 0
                && !util::ivecaxpygez(
                    -1,
                    &vwgt[cntrng!(i * ncon, ncon)],
                    &pwgts[cntrng!(from * ncon, ncon)],
                    &minpwgts[cntrng!(from * ncon, ncon)],
                )
                {
                    continue;
                }
            } else {
                /* OMODE_BALANCE */
                if !util::ivecaxpygez(
                    -1,
                    &vwgt[cntrng!(i * ncon, ncon)],
                    &pwgts[cntrng!(from * ncon, ncon)],
                    &minpwgts[cntrng!(from * ncon, ncon)],
                ) {
                    continue;
                }
            }

            if ctrl.contig != 0 && IsArticulationNode(nvtxs as idx_t, i as idx_t, xadj.as_ptr(), adjncy.as_ptr(), where_.as_ptr(), bfslvl.as_mut_ptr(), bfsind.as_mut_ptr(), bfsmrk.as_mut_ptr()) != 0
            {
                continue;
            }

            if ctrl.minconn {
                SelectSafeTargetSubdomains(myrinfo, mynbrs, nads, adids, maxndoms, safetos, doms);
            }

            let xgain = if myrinfo.nid == 0 && myrinfo.ned > 0 {
                graph.vsize[i as usize]
            } else {
                0
            };

            /* Find the most promising subdomain to move to */
            let mut to;
            let k = if omode == OMODE_REFINE {
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     if !safetos[to = mynbrs[k as usize].pid] {
                //         continue;
                //     }
                //     gain = mynbrs[k as usize].gv + xgain;
                //     if gain >= 0
                //     && ivecaxpylez(
                //         ncon,
                //         1,
                //         vwgt + i * ncon,
                //         pwgts + to * ncon,
                //         maxpwgts + to * ncon,
                //     )
                //     {
                //         break;
                //     }
                // }
                // if k < 0 {
                //     continue;
                // } /* break out if you did not find a candidate */

                let Some(mut k) = (0..myrinfo.nnbrs as usize)
                    .filter(|&k| safetos[mynbrs[k].pid as usize] != 0)
                    .rfind(|&k| {
                        to = mynbrs[k].pid as usize;
                        let gain = mynbrs[k as usize].gv + xgain;
                        gain >= 0
                        && util::ivecaxpylez(
                            1,
                            &vwgt[cntrng!(i * ncon, ncon)],
                            &pwgts[cntrng!(to * ncon, ncon)],
                            &maxpwgts[cntrng!(to * ncon, ncon)],
                        )
                    }) else { break };

                let mut cto = to;
                for j in (0..k).rev() {
                    let to = mynbrs[j as usize].pid as usize;
                    if safetos[to] == 0 {
                        continue;
                    }
                    let gain = mynbrs[j as usize].gv + xgain;
                    if (mynbrs[j as usize].gv > mynbrs[k as usize].gv
                    && util::ivecaxpylez(
                        1,
                        &vwgt[cntrng!(i * ncon, ncon)],
                        &pwgts[cntrng!(to * ncon, ncon)],
                        &maxpwgts[cntrng!(to * ncon, ncon)],
                    ))
                    || (mynbrs[j as usize].gv == mynbrs[k as usize].gv
                    && mynbrs[j as usize].ned > mynbrs[k as usize].ned
                    && util::ivecaxpylez(
                        1,
                        &vwgt[cntrng!(i * ncon, ncon)],
                        &pwgts[cntrng!(to * ncon, ncon)],
                        &maxpwgts[cntrng!(to * ncon, ncon)],
                    ))
                    || (mynbrs[j as usize].gv == mynbrs[k as usize].gv
                    && mynbrs[j as usize].ned == mynbrs[k as usize].ned
                    && mcutil::BetterBalanceKWay(
                        ncon as idx_t,
                        vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                        ubfactors,
                        1,
                        pwgts[cntrng!(cto * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(cto * ncon, ncon)].as_ptr(),
                        1,
                        pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                    ) != 0)
                    {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                let mut j = 0;
                if xgain + mynbrs[k as usize].gv > 0 || mynbrs[k as usize].ned - myrinfo.nid > 0 {
                    j = 1;
                } else if (mynbrs[k as usize].ned - myrinfo.nid == 0) {
                    if (iii % 2 == 0 && safetos[to as usize] == 2)
                    || mcutil::BetterBalanceKWay(
                        ncon as idx_t,
                        vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                        ubfactors,
                        -1,
                        pwgts[cntrng!(from * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(from * ncon, ncon)].as_ptr(),
                        1,
                        pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                    ) != 0
                    {
                        j = 1;
                    }
                }
                if j == 0 {
                    continue;
                }

                k
            } else {
                /* OMODE_BALANCE */
                // for k in (0..=(myrinfo.nnbrs - 1)).rev() {
                //     if !safetos[to = mynbrs[k as usize].pid] {
                //         continue;
                //     }
                //     if ivecaxpylez(
                //         ncon,
                //         1,
                //         vwgt + i * ncon,
                //         pwgts + to * ncon,
                //         maxpwgts + to * ncon,
                //     ) || BetterBalanceKWay(
                //         ncon,
                //         vwgt + i * ncon,
                //         ubfactors,
                //         -1,
                //         pwgts + from * ncon,
                //         pijbm + from * ncon,
                //         1,
                //         pwgts + to * ncon,
                //         pijbm + to * ncon,
                //     ) {
                //         break;
                //     }
                // }
                // if k < 0 {
                //     continue;
                // } /* break out if you did not find a candidate */

                let Some(mut k) = (0..myrinfo.nnbrs as usize)
                    .filter(|&k| safetos[mynbrs[k].pid as usize] != 0)
                    .rfind(|&k| {
                        util::ivecaxpylez(
                        1,
                        &vwgt[cntrng!(i * ncon, ncon)],
                        &pwgts[cntrng!(to * ncon, ncon)],
                        &maxpwgts[cntrng!(to * ncon, ncon)],
                    ) || mcutil::BetterBalanceKWay(
                        ncon as idx_t,
                        vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                        ubfactors.as_ptr(),
                        -1,
                        pwgts[cntrng!(from * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(from * ncon, ncon)].as_ptr(),
                        1,
                        pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                    ) != 0
                    }) else { continue };

                let mut cto = to;
                for j in (0..=(k - 1)).rev() {
                    let to = mynbrs[j as usize].pid as usize;
                    if safetos[to] == 0 {
                        continue;
                    }
                    if mcutil::BetterBalanceKWay(
                        ncon as idx_t,
                        vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                        ubfactors.as_ptr(),
                        1,
                        pwgts[cntrng!(cto * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(cto * ncon, ncon)].as_ptr(),
                        1,
                        pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                        pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                    ) != 0 {
                        k = j;
                        cto = to;
                    }
                }
                to = cto;

                if (xgain + mynbrs[k as usize].gv < 0
                || (xgain + mynbrs[k as usize].gv == 0
                && mynbrs[k as usize].ned - myrinfo.nid < 0))
                && mcutil::BetterBalanceKWay(
                    ncon as idx_t,
                    vwgt[cntrng!(i * ncon, ncon)].as_ptr(),
                    ubfactors.as_ptr(),
                    -1,
                    pwgts[cntrng!(from * ncon, ncon)].as_ptr(),
                    pijbm[cntrng!(from * ncon, ncon)].as_ptr(),
                    1,
                    pwgts[cntrng!(to * ncon, ncon)].as_ptr(),
                    pijbm[cntrng!(to * ncon, ncon)].as_ptr(),
                ) == 0
                {
                    continue;
                }

                k
            };

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
                print!("\t\tMoving {:6} from {:3} to {:3}. \
                Gain: [({:4} {:4}) as usize]. Cut: {:6}, Vol: {:6}\n", 
              i, from, to, xgain+mynbrs[k as usize].gv, mynbrs[k as usize].ned-myrinfo.nid, 
              graph.mincut, graph.minvol)
            );

            /* Update the subdomain connectivity information */
            if ctrl.minconn != 0 {
                /* take care of i's move itself */
                minconn::UpdateEdgeSubDomainGraph(
                    ctrl,
                    from,
                    to,
                    myrinfo.nid - mynbrs[k as usize].ned,
                    &maxndoms,
                );

                /* take care of the adjacent vertices */
                for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                    let me = where_[adjncy[j as usize] as usize] as usize;
                    if me != from && me != to {
                        minconn::UpdateEdgeSubDomainGraph(ctrl, from as idx_t, me as idx_t, -1, &mut maxndoms);
                        minconn::UpdateEdgeSubDomainGraph(ctrl, to as idx_t, me as idx_t, 1, &mut maxndoms);
                    }
                }
            }

            /* Update pwgts */
            blas::iaxpy(ncon, 1, &vwgt[cntrng!(i * ncon, ncon)], 1, &mut pwgts[cntrng!(to * ncon, ncon)], 1);
            blas::iaxpy(ncon, -1, &vwgt[cntrng!(i * ncon, ncon)], 1, &mut pwgts[cntrng!(from * ncon, ncon)], 1);

            /* Update the id/ed/gains/bnd/queue of potentially affected nodes */
            KWayVolUpdate(
                ctrl, graph, i, from, to, queue, vstatus, &nupd, updptr, updind, bndtype, vmarker,
                pmarker, modind,
            );

            /*CheckKWayVolPartitionParams(ctrl, graph); */
        }

        /* Reset the vstatus and associated data structures */
        for i in (0)..(nupd) {
            debug_assert!(updptr[updind[i as usize] as usize] != -1);
            debug_assert!(vstatus[updind[i as usize] as usize] != VPQSTATUS_NOTPRESENT);
            vstatus[updind[i as usize] as usize] = VPQSTATUS_NOTPRESENT;
            updptr[updind[i as usize] as usize] = -1;
        }

        if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}. \
              Nmoves: {:5}, Cut: {:6}, Vol: {:6}",
              imin(nparts*ncon, pwgts,1), imax(nparts*ncon, pwgts,1), 
              mcutil::ComputeLoadImbalance(graph, nparts as idx_t, pijbm.as_ptr()), 
              graph.nbnd, nmoved, graph.mincut, graph.minvol);
            if ctrl.minconn != 0 {
                print!(
                    ", Doms: [({:3} {:4}) as usize]",
                    imax(nparts, nads, 1),
                    isum(nparts, nads, 1)
                );
            }
            print!("\n");
        }

        if nmoved == 0
        || (omode == OMODE_REFINE && graph.minvol == oldvol && graph.mincut == oldcut)
        {
            break;
        }
    }

    // ipqDestroy(queue);

    // WCOREPOP;
}

/*************************************************************************/
/* This function performs an approximate articulation vertex test.
It assumes that the bfslvl, bfsind, and bfsmrk arrays are initialized
appropriately. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn IsArticulationNode(
    nvtxs: idx_t,
    i: idx_t,
    xadj: *const idx_t,
    adjncy: *const idx_t,
    where_: *const idx_t,
    bfslvl: *mut idx_t,
    bfsind: *mut idx_t,
    bfsmrk: *mut idx_t,
) -> idx_t {
    // idx_t ii, j, k=0, head, tail, nhits, tnhits, from, BFSDEPTH=5;
    const BFSDEPTH: idx_t = 5;

    let nvtxs = nvtxs as usize;
    let i = i as usize;
    mkslice!(xadj, nvtxs + 1);
    mkslice!(adjncy, xadj[nvtxs]);
    mkslice!(where_, nvtxs);
    mkslice_mut!(bfsmrk, nvtxs);
    mkslice_mut!(bfslvl, nvtxs);
    mkslice_mut!(bfsind, nvtxs);

    let from = where_[i as usize] as usize;

    /* Determine if the vertex is safe to move from a contiguity standpoint */
    let mut tnhits = 0;
    let mut k = 0;
    for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
        if where_[adjncy[j as usize] as usize] == from as idx_t {
            debug_assert!(bfsmrk[adjncy[j as usize] as usize] == 0);
            debug_assert!(bfslvl[adjncy[j as usize] as usize] == 0);
            k = adjncy[j as usize] as usize;
            bfsmrk[k] = 1;
            tnhits += 1;
        }
    }

    /* Easy cases */
    if tnhits == 0 {
        return 0;
    }
    if tnhits == 1 {
        bfsmrk[k as usize] = 0;
        return 0;
    }

    debug_assert!(bfslvl[i as usize] == 0);
    bfslvl[i as usize] = 1;

    bfsind[0] = k; /* That was the last one from the previous loop */
    bfslvl[k as usize] = 1;
    bfsmrk[k as usize] = 0;
    let mut head = 1;
    let mut tail = 1;

    /* Do a limited BFS traversal to see if you can get to all the other nodes */
    let mut nhits = 1;
    while head < tail {
        let ii = bfsind[head] as usize;
        head += 1;
        for j in (xadj[ii])..(xadj[ii + 1]) {
            k = adjncy[j as usize] as usize;
            if where_[k] == from as idx_t {
                if bfsmrk[k] != 0 {
                    bfsmrk[k] = 0;
                    nhits += 1;
                    if nhits == tnhits {
                        break;
                    }
                }
                if bfslvl[k] == 0 && bfslvl[ii] < BFSDEPTH {
                    bfsind[tail] = k as idx_t;
                    tail += 1;
                    bfslvl[k] = bfslvl[ii] + 1;
                }
            }
        }
        if nhits == tnhits {
            break;
        }
    }

    /* Reset the various BFS related arrays */
    bfslvl[i as usize] = 0;
    for j in (0)..(tail) {
        bfslvl[bfsind[j as usize] as usize] = 0;
    }

    /* Reset the bfsmrk array for the next vertex when has not already being cleared */
    if nhits < tnhits {
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            if where_[adjncy[j as usize] as usize] == from as idx_t {
                bfsmrk[adjncy[j as usize] as usize] = 0;
            }
        }
    }

    return (nhits != tnhits) as idx_t;
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
// NOTE(porting): this can't be a metis_func since I don't want to have to make queue abi
// compatible
#[allow(non_snake_case)]
pub unsafe fn KWayVolUpdate(
    ctrl: &mut ctrl_t,
    graph: &mut graph_t,
    v: usize,
    from: usize,
    to: usize,
    queue: Option<&mut pqueue::IPQueue>,
    vstatus: Option<&mut [idx_t]>,
    r_nupd: Option<&mut idx_t>,
    updptr: Option<&mut [idx_t]>,
    updind: Option<&mut [idx_t]>,
    bndtype: idx_t,
    vmarker: &mut [idx_t],
    pmarker: &mut [idx_t],
    modind: &mut [idx_t],
) {
    let mut queue = queue;
    let mut vstatus = vstatus;
    let mut r_nupd = r_nupd;
    let mut updptr = updptr;
    let mut updind = updind;
    // idx_t i, ii, iii, j, jj, k, kk, l, u, nmod, other, me, myidx;
    // idx_t *xadj, *vsize, *adjncy, *where_;
    // vkrinfo_t *myrinfo, *orinfo;
    // vnbr_t *mynbrs, *onbrs;

    get_graph_slices!(graph => xadj adjncy vsize where_);

    let (myrinfo, mynbrs) = vrinfos(graph.vkrinfo, ctrl.vnbrpool, v);

    /*======================================================================
     * Remove the contributions on the gain made by 'v'.
     *=====================================================================*/
    for k in (0)..(myrinfo.nnbrs as usize) {
        pmarker[mynbrs[k as usize].pid as usize] = k as idx_t;
    }
    pmarker[from as usize] = myrinfo.nnbrs;

    // has to be signed
    let myidx: idx_t = pmarker[to as usize]; /* Keep track of the index in mynbrs of the 'to' domain */

    for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
        let ii = adjncy[j as usize] as usize;
        let other = where_[ii as usize] as usize;
        let (orinfo, onbrs) = vrinfos(graph.vkrinfo, ctrl.vnbrpool, ii);

        if other == from {
            for k in (0)..(orinfo.nnbrs) {
                if pmarker[onbrs[k as usize].pid as usize] == -1 {
                    onbrs[k as usize].gv += vsize[v as usize];
                }
            }
        } else {
            debug_assert!(pmarker[other as usize] != -1);

            if mynbrs[pmarker[other as usize] as usize].ned > 1 {
                for k in (0)..(orinfo.nnbrs as usize) {
                    if pmarker[onbrs[k].pid as usize] == -1 {
                        onbrs[k].gv += vsize[v as usize];
                    }
                }
            } else {
                /* There is only one connection */
                for k in (0)..(orinfo.nnbrs as usize) {
                    if pmarker[onbrs[k].pid as usize] != -1 {
                        onbrs[k].gv -= vsize[v as usize];
                    }
                }
            }
        }
    }

    for k in (0)..(myrinfo.nnbrs as usize) {
        pmarker[mynbrs[k].pid as usize] = -1;
    }
    pmarker[from as usize] = -1;

    /*======================================================================
     * Update the id/ed of vertex 'v'
     *=====================================================================*/
    if myidx == -1 {
        myidx = myrinfo.nnbrs;
        myrinfo.nnbrs += 1;
        debug_assert!(myidx < xadj[(v + 1) as usize] - xadj[v as usize]);
        mynbrs[myidx as usize].ned = 0;
    }
    myrinfo.ned += myrinfo.nid - mynbrs[myidx as usize].ned;
    std::mem::swap(&mut myrinfo.nid, &mut mynbrs[myidx as usize].ned);
    if mynbrs[myidx as usize].ned == 0 {
        myrinfo.nnbrs -= 1;
        mynbrs[myidx as usize] = mynbrs[myrinfo.nnbrs as usize];
    } else {
        mynbrs[myidx as usize].pid = from as idx_t;
    }

    /*======================================================================
     * Update the degrees of adjacent vertices and their volume gains
     *=====================================================================*/
    vmarker[v] = 1;
    modind[0] = v as idx_t;
    let mut nmod: usize = 1;
    for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
        let ii = adjncy[j as usize] as usize;
        let me = where_[ii as usize] as usize;

        if vmarker[ii] == 0 {
            /* The marking is done for boundary and max gv calculations */
            vmarker[ii] = 2;
            modind[(nmod) as usize] = ii as idx_t;
            nmod += 1;
        }

        let myrinfo = graph.vkrinfo.add(ii);
        if (*myrinfo).inbr == -1 {
            (*myrinfo).inbr = vnbrpoolGetNext(ctrl, xadj[(ii + 1) as usize] - xadj[ii]);
        }
        let (myrinfo, mynbrs) = vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, ii);

        if me == from {
            inc_dec!(myrinfo.ned, myrinfo.nid, 1);
        } else if me == to {
            inc_dec!(myrinfo.nid, myrinfo.ned, 1);
        }

        /* Remove the edgeweight from the 'pid == from' entry of the vertex */
        if me != from {
            for k in (0)..(myrinfo.nnbrs) {
                if mynbrs[k as usize].pid == from as idx_t {
                    if mynbrs[k as usize].ned == 1 {
                        myrinfo.nnbrs -= 1;
                        mynbrs[k as usize] = mynbrs[myrinfo.nnbrs as usize];
                        vmarker[ii] = 1; /* You do a complete .gv calculation */

                        /* All vertices adjacent to 'ii' need to be updated */
                        for jj in (xadj[ii])..(xadj[(ii + 1) as usize]) {
                            let u = adjncy[jj as usize] as usize;
                            let other = where_[u as usize] as usize;
                            let (orinfo, onbrs) = vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, u);

                            for kk in (0)..(orinfo.nnbrs) {
                                if onbrs[kk as usize].pid == from as idx_t {
                                    onbrs[kk as usize].gv -= vsize[ii];
                                    if vmarker[u as usize] == 0 {
                                        /* Need to update boundary etc */
                                        vmarker[u as usize] = 2;
                                        modind[(nmod) as usize] = u as idx_t;
                                        nmod += 1;
                                    }
                                    break;
                                }
                            }
                        }
                    } else {
                        mynbrs[k as usize].ned -= 1;

                        /* Update the gv due to single 'ii' connection to 'from' */
                        if mynbrs[k as usize].ned == 1 {
                            /* find the vertex 'u' that 'ii' was connected into 'from' */
                            for jj in (xadj[ii])..(xadj[(ii + 1) as usize]) {
                                let u = adjncy[jj as usize] as usize;
                                let other = where_[u as usize] as usize;

                                if other == from {
                                    let (orinfo, onbrs) = vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, u);

                                    /* The following is correct because domains in common
                                    between ii and u will lead to a reduction over the
                                    previous gain, where_as domains only in u but not in
                                    ii, will lead to no change as opposed to the earlier
                                    increase */
                                    for nbr in onbrs {
                                        nbr.gv += vsize[ii];
                                    }

                                    if vmarker[u as usize] == 0 {
                                        /* Need to update boundary etc */
                                        vmarker[u as usize] = 2;
                                        modind[nmod as usize] = u as idx_t;
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
        if me != to {
            let mut found = false;
            for k in (0)..(myrinfo.nnbrs) {
                if mynbrs[k as usize].pid == to as idx_t {
                    mynbrs[k as usize].ned += 1;

                    /* Update the gv due to non-single 'ii' connection to 'to' */
                    if mynbrs[k as usize].ned == 2 {
                        /* find the vertex 'u' that 'ii' was connected into 'to' */
                        for jj in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
                            let u = adjncy[jj as usize] as usize;
                            let other = where_[u as usize] as usize;

                            if u != v && other == to {
                                let (orinfo, onbrs) = vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, u);
                                for kk in (0)..(orinfo.nnbrs) {
                                    onbrs[kk as usize].gv -= vsize[ii as usize];
                                }

                                if vmarker[u] == 0 {
                                    /* Need to update boundary etc */
                                    vmarker[u] = 2;
                                    modind[nmod as usize] = u as idx_t;
                                    nmod += 1;
                                }
                                found = true;
                                break;
                            }
                        }
                    }
                    found = true;
                    break;
                }
            }

            if !found {
                mynbrs[myrinfo.nnbrs as usize].pid = to as idx_t;
                mynbrs[(myrinfo.nnbrs) as usize].ned = 1;
                myrinfo.nnbrs += 1;
                vmarker[ii as usize] = 1; /* You do a complete .gv calculation */

                /* All vertices adjacent to 'ii' need to be updated */
                for jj in (xadj[ii as usize])..(xadj[(ii + 1) as usize]) {
                    let u = adjncy[jj as usize] as usize;
                    let other = where_[u];
                    let (orinfo, onbrs) = vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, u);

                    for kk in (0)..(orinfo.nnbrs) {
                        if onbrs[kk as usize].pid == to as idx_t {
                            onbrs[kk as usize].gv += vsize[ii as usize];
                            if vmarker[u as usize] == 0 {
                                /* Need to update boundary etc */
                                vmarker[u as usize] = 2;
                                modind[nmod as usize] = u as idx_t;
                                nmod += 1;
                            }
                            break;
                        }
                    }
                }
            }
        }

        debug_assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
    }

    /*======================================================================
     * Add the contributions on the volume gain due to 'v'
     *=====================================================================*/
    let (myrinfo, mynbrs) = vrinfos(graph.vkrinfo, ctrl.vnbrpool, v);
    for k in (0)..(myrinfo.nnbrs) {
        pmarker[mynbrs[k as usize].pid as usize] = k;
    }
    pmarker[to as usize] = myrinfo.nnbrs;

    for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
        let ii = adjncy[j as usize] as usize;
        let other = where_[ii] as usize;
        let (orinfo, onbrs) = vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, ii);

        if other == to {
            for k in (0)..(orinfo.nnbrs) {
                if pmarker[onbrs[k as usize].pid as usize] == -1 {
                    onbrs[k as usize].gv -= vsize[v as usize];
                }
            }
        } else {
            debug_assert!(pmarker[other as usize] != -1);

            if mynbrs[pmarker[other as usize] as usize].ned > 1 {
                for k in (0)..(orinfo.nnbrs) {
                    if pmarker[onbrs[k as usize].pid as usize] == -1 {
                        onbrs[k as usize].gv -= vsize[v as usize];
                    }
                }
            } else {
                /* There is only one connection */
                for k in (0)..(orinfo.nnbrs) {
                    if pmarker[onbrs[k as usize].pid as usize] != -1 {
                        onbrs[k as usize].gv += vsize[v as usize];
                    }
                }
            }
        }
    }
    for k in (0)..(myrinfo.nnbrs) {
        pmarker[mynbrs[k as usize].pid as usize] = -1;
    }
    pmarker[to as usize] = -1;

    /*======================================================================
     * Recompute the volume information of the 'hard' nodes, and update the
     * max volume gain for all the modified vertices and the priority queue
     *=====================================================================*/
    for iii in (0)..(nmod) {
        let i = modind[iii] as usize;
        let me = where_[i as usize] as usize;

        let (myrinfo, mynbrs) = vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, i);

        if vmarker[i as usize] == 1 {
            /* Only complete gain updates go through */
            for k in (0)..(myrinfo.nnbrs) {
                mynbrs[k as usize].gv = 0;
            }

            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                let ii = adjncy[j as usize] as usize;
                let other = where_[ii as usize] as usize;
                let (orinfo, onbrs) = vrinfos(graph.vkrinfo, ctrl.vnbrpool, ii);

                for kk in (0)..(orinfo.nnbrs) {
                    pmarker[onbrs[kk as usize].pid as usize] = kk;
                }
                pmarker[other as usize] = 1;

                if me == other {
                    /* Find which domains 'i' is connected and 'ii' is not and update their gain */
                    for k in (0)..(myrinfo.nnbrs) {
                        if pmarker[mynbrs[k as usize].pid as usize] == -1 {
                            mynbrs[k as usize].gv -= vsize[ii as usize];
                        }
                    }
                } else {
                    debug_assert!(pmarker[me as usize] != -1);

                    /* I'm the only connection of 'ii' in 'me' */
                    if onbrs[pmarker[me as usize] as usize].ned == 1 {
                        /* Increase the gains for all the common domains between 'i' and 'ii' */
                        for k in (0)..(myrinfo.nnbrs) {
                            if pmarker[mynbrs[k as usize].pid as usize] != -1 {
                                mynbrs[k as usize].gv += vsize[ii as usize];
                            }
                        }
                    } else {
                        /* Find which domains 'i' is connected and 'ii' is not and update their gain */
                        for k in (0)..(myrinfo.nnbrs) {
                            if pmarker[mynbrs[k as usize].pid as usize] == -1 {
                                mynbrs[k as usize].gv -= vsize[ii as usize];
                            }
                        }
                    }
                }

                for kk in (0)..(orinfo.nnbrs) {
                    pmarker[onbrs[kk as usize].pid as usize] = -1;
                }
                pmarker[other as usize] = -1;
            }
        }

        /* Compute the overall gv for that node */
        myrinfo.gv = idx_t::MIN;
        for k in (0)..(myrinfo.nnbrs) {
            if mynbrs[k as usize].gv > myrinfo.gv {
                myrinfo.gv = mynbrs[k as usize].gv;
            }
        }

        /* Add the xtra gain due to id == 0 */
        if myrinfo.ned > 0 && myrinfo.nid == 0 {
            myrinfo.gv += vsize[i as usize];
        }

        /*======================================================================
         * Maintain a consistent boundary
         *=====================================================================*/
        {
            get_graph_slices_mut!(graph => bndind bndptr);
            if bndtype == BNDTYPE_REFINE {
                if myrinfo.gv >= 0 && bndptr[i as usize] == -1 {
                    BNDInsert!(graph.nbnd, bndind, bndptr, i);
                }

                if myrinfo.gv < 0 && bndptr[i as usize] != -1 {
                    BNDDelete!(graph.nbnd, bndind, bndptr, i);
                }
            } else {
                if myrinfo.ned > 0 && bndptr[i as usize] == -1 {
                    BNDInsert!(graph.nbnd, bndind, bndptr, i);
                }

                if myrinfo.ned == 0 && bndptr[i as usize] != -1 {
                    BNDDelete!(graph.nbnd, bndind, bndptr, i);
                }
            }
        }

        /*======================================================================
         * Update the priority queue appropriately (if allowed)
         *=====================================================================*/
        if let Some(queue) = queue.as_mut() {
            get_graph_slices!(graph => bndptr);
            let vstatus = &mut **vstatus.as_mut().unwrap();
            let updind = &mut **updind.as_mut().unwrap();
            let updptr = &mut **updptr.as_mut().unwrap();
            let r_nupd = &mut **r_nupd.as_mut().unwrap();
            if vstatus[i as usize] != VPQSTATUS_EXTRACTED {
                if bndptr[i as usize] != -1 {
                    /* In-boundary vertex */
                    if vstatus[i as usize] == VPQSTATUS_PRESENT {
                        queue.update(i as idx_t, myrinfo.gv);
                    } else {
                        queue.insert(i as idx_t, myrinfo.gv);
                        vstatus[i as usize] = VPQSTATUS_PRESENT;
                        ListInsert!(*r_nupd, updind, updptr, i);
                    }
                } else {
                    /* Off-boundary vertex */
                    if vstatus[i as usize] == VPQSTATUS_PRESENT {
                        queue.delete(i as idx_t);
                        vstatus[i as usize] = VPQSTATUS_NOTPRESENT;
                        ListDelete!(*r_nupd, updind, updptr, i);
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
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
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
    get_graph_slices!(ctrl, graph => xadj adjncy vwgt adjwgt bndind where_ pwgts);

    let nparts = ctrl.nparts as usize;
    mkslice!(ctrl->tpwgts, nparts);

    /* Setup the weight intervals of the various subdomains */
    let mut minpwgts: Vec<idx_t> = vec![0; nparts as usize];
    let mut maxpwgts: Vec<idx_t> = vec![0; nparts as usize];

    let ubfactor = *ctrl.ubfactors;
    for i in (0)..(nparts) {
        maxpwgts[i as usize] = (tpwgts[i as usize] * *graph.tvwgt as f32 * ubfactor) as idx_t;
        minpwgts[i as usize] = (tpwgts[i as usize] * *graph.tvwgt as f32 * (0.95 / ubfactor)) as idx_t;
    }

    /* go and determine the positive gain valid swaps */
    let nbnd = graph.nbnd as usize;

    for ii in (0)..(nbnd) {
        let u = bndind[ii as usize] as usize;
        let uw = where_[u as usize] as usize;

        let (urinfo, unbrs) = crinfos(graph.ckrinfo, ctrl.cnbrpool, u);

        for j in (xadj[u as usize])..(xadj[(u + 1) as usize]) {
            let v = adjncy[j as usize] as usize;
            let vw = where_[v as usize] as usize;

            let (vrinfo, vnbrs) = crinfos(graph.ckrinfo, ctrl.cnbrpool, v);

            if uw == vw {
                continue;
            }

            if pwgts[uw as usize] - vwgt[u as usize] + vwgt[v as usize] > maxpwgts[uw as usize]
            || pwgts[vw as usize] - vwgt[v as usize] + vwgt[u as usize] > maxpwgts[vw as usize]
            {
                continue;
            }

            // TODO: the original was liable to index out of bounds here (both of thse loops that
            // was replaced with rfind) if we didn't find an approprite neighbour. I should file a
            // PR

            let Some(unbr) = unbrs.iter().rfind(|unbr| unbr.pid == vw as idx_t) else {
                println!("Something went wrong!");
                continue
            };
            let mut gain = unbr.ed - urinfo.id;

            let Some(vnbr) = vnbrs.iter().rfind(|vnbr| vnbr.pid == uw as idx_t) else {
                println!("Something went wrong!");
                continue
            };
            gain += vnbr.ed - vrinfo.id;

            gain -= 2 * adjwgt[j as usize];

            if gain > 0 {
                println!(
                    "  Gain: {:} for moving ({:}, {:}) between ({:}, {:})",
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
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
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

    let bndtype = BNDTYPE_REFINE;

    /* Link the graph fields */
    let nvtxs = graph.nvtxs as usize;
    let nparts = ctrl.nparts as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy adjwgt vwgt bndind bndptr where_ pwgts tvwgt);
    mkslice!(ctrl->tpwgts, nparts);

    /* Setup the weight intervals of the various subdomains */
    let mut minpwgts: Vec<idx_t> = vec![0; nparts as usize];
    let mut maxpwgts: Vec<idx_t> = vec![0; nparts as usize];

    let ubfactor = (*ctrl.ubfactors).max(mcutil::ComputeLoadImbalance(graph, nparts as idx_t, ctrl.pijbm));
    for k in (0)..(nparts) {
        maxpwgts[k as usize] = (tpwgts[k as usize] * tvwgt[0] as real_t * ubfactor) as idx_t;
        minpwgts[k as usize] = (tpwgts[k as usize] * tvwgt[0] as real_t * (1.0 / ubfactor)) as idx_t;
    }

    let mut perm: Vec<idx_t> = vec![0; nvtxs as usize];

    if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
        print!("GRE: [({:6} {:6}) as usize]-[({:6} {:6}) as usize], Bal: {:5.3}, \
            Nv-Nb[({:6} {:6}) as usize], Cut: {:6}\n",
            pwgts[(util::iargmin(pwgts, 1)) as usize], imax(nparts, pwgts,1), minpwgts[0], maxpwgts[0], 
            mcutil::ComputeLoadImbalance(graph, nparts as idx_t, ctrl.pijbm), 
            graph.nvtxs, graph.nbnd, graph.mincut);
    }

    /*=====================================================================
     * The top-level refinement loop
     *======================================================================*/
    for pass in (0)..(niter) {
        // this is a GKASSERT (always runs) in the original
        assert!(debug::ComputeCut(graph, where_.as_ptr()) == graph.mincut);

        let oldcut = graph.mincut;
        let nbnd = graph.nbnd as usize;
        let mut nmoved = 0;

        /* Insert the boundary vertices in the priority queue */
        /* Visit the vertices in random order and see if you can swap them */
        irandArrayPermute(nvtxs as idx_t, perm.as_mut_ptr(), nbnd as idx_t, 1);
        for ii in (0)..(nvtxs) {
            let u = perm[ii as usize] as usize;
            if bndptr[u] == -1 {
                continue;
            }

            let uw = where_[u as usize] as usize;

            let (urinfo, unbrs) = crinfos(graph.ckrinfo, ctrl.cnbrpool, u);

            let mut bestgain = 0;
            let mut jbest = None;
            for j in (xadj[u as usize])..(xadj[(u + 1) as usize]) {
                let v = adjncy[j as usize] as usize;
                let vw = where_[v as usize] as usize ;

                if uw == vw {
                    continue;
                }
                if pwgts[uw as usize] - vwgt[u as usize] + vwgt[v as usize]
                > maxpwgts[uw as usize]
                || pwgts[vw as usize] - vwgt[v as usize] + vwgt[u as usize]
                > maxpwgts[vw as usize]
                {
                    continue;
                }
                if pwgts[uw as usize] - vwgt[u as usize] + vwgt[v as usize]
                < minpwgts[uw as usize]
                || pwgts[vw as usize] - vwgt[v as usize] + vwgt[u as usize]
                < minpwgts[vw as usize]
                {
                    continue;
                }

                let (vrinfo, vnbrs) = crinfos(graph.ckrinfo, ctrl.cnbrpool, v);

                let mut gain = -2 * adjwgt[j as usize];


                let Some(unbr) = unbrs.iter().rfind(|unbr| unbr.pid == vw as idx_t) else {
                    // original debug_assert message is "k >= 0"
                    panic!("no neighbor pid matches vw")
                };
                gain += unbr.ed - urinfo.id;

                let Some(vnbr) = vnbrs.iter().rfind(|vnbr| vnbr.pid == uw as idx_t) else {
                    // original debug_assert message is "k >= 0"
                    panic!("no neighbor pid matches uw")
                };
                gain += vnbr.ed - vrinfo.id;

                if gain > bestgain && vnbr.ed > adjwgt[j as usize] {
                    bestgain = gain;
                    jbest = Some(j);
                }
            }

            let Some(jbest) = jbest else {
                // no valid positive swap
                continue;
            };

            /*=====================================================================
             * If we got here, we can now swap the vertices
             *======================================================================*/
            let v = adjncy[jbest as usize] as usize;
            let vw = where_[v as usize] as usize;

            let (vrinfo, vnbrs) = crinfos(graph.ckrinfo, ctrl.cnbrpool, v);

            /* move u to v's partition */
            let Some(k) = unbrs.iter().rposition(|nbr| nbr.pid == vw as idx_t) else {
                // original debug_assert message: "k >= 0"
                unreachable!("should have already found unbr with unbr.pid == vw");
            };

            let from = uw;
            let to = vw;

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
            inc_dec!(pwgts[to as usize], pwgts[from as usize], vwgt[u as usize]);
            UpdateMovedVertexInfoAndBND!(
                u, from, k, to, urinfo, unbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[u as usize])..(xadj[(u + 1) as usize]) {
                let ii = adjncy[j as usize] as usize;
                let me = where_[ii as usize] as usize;
                let myrinfo = &mut *graph.ckrinfo.add(ii);

                let oldnnbrs = myrinfo.nnbrs;

                UpdateAdjacentVertexInfoAndBND(
                    ctrl,
                    ii,
                    xadj[(ii + 1) as usize] - xadj[ii as usize],
                    me,
                    from,
                    to,
                    myrinfo,
                    adjwgt[j as usize],
                    &mut nbnd,
                    bndptr,
                    bndind,
                    bndtype,
                );

                debug_assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }

            /* move v to u's partition */
            // for k in (0..=(vrinfo.nnbrs - 1)).rev() {
            //     if (vnbrs[k as usize].pid == uw) {
            //         break;
            //     }
            // }
            // debug_assert!(k >= 0);
            let k = (0..vrinfo.nnbrs as usize)
                .rfind(|&k| vnbrs[k].pid == uw as idx_t)
                .expect("k >= 0");

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
            inc_dec!(pwgts[to as usize], pwgts[from as usize], vwgt[v as usize]);
            UpdateMovedVertexInfoAndBND!(
                v, from, k, to, vrinfo, vnbrs, where_, nbnd, bndptr, bndind, bndtype,
            );

            /* Update the degrees of adjacent vertices */
            for j in (xadj[v as usize])..(xadj[(v + 1) as usize]) {
                let ii = adjncy[j as usize] as usize;
                let me = where_[ii as usize];
                let myrinfo = &mut *graph.ckrinfo.add(ii);

                let oldnnbrs = myrinfo.nnbrs;

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

                debug_assert!(myrinfo.nnbrs <= xadj[(ii + 1) as usize] - xadj[ii as usize]);
            }
        }

        graph.nbnd = nbnd as idx_t;

        if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
            print!("\t[({:6} {:6}) as usize], Bal: {:5.3}, Nb: {:6}. \
              Nmoves: {:5}, Cut: {:6}, Vol: {:6}\n",
              pwgts[(util::iargmin(pwgts, 1)) as usize], imax(nparts, pwgts,1),
              mcutil::ComputeLoadImbalance(graph, nparts as idx_t, ctrl.pijbm), 
              graph.nbnd, nmoved, graph.mincut, debug::ComputeVolume(graph, where_.as_ptr()));
        }

        if nmoved == 0 || graph.mincut == oldcut {
            break;
        }
    }

    // WCOREPOP;
}
