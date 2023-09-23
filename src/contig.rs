/*
\file
\brief Functions that deal with eliminating disconnected partitions

\date Started 7/15/98
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: contig.c 10513 2011-07-07 22:06:03Z karypis $
*/

use crate::*;

/*************************************************************************/
/* This function finds the connected components induced by the
    partitioning vector.

    \param graph is the graph structure
    \param where_ is the partitioning vector. If this is NULL, then the
           entire graph is treated to belong into a single partition.
    \param cptr is the ptr structure of the CSR representation of the
           components. The length of this vector must be graph.nvtxs+1.
    \param cind is the indices structure of the CSR representation of
           the components. The length of this vector must be graph.nvtxs.

    \returns the number of components that it found.

    \note The cptr and cind parameters can be NULL, in which case only the
          number of connected components is returned.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn FindPartitionInducedComponents(
    graph: *mut graph_t,
    where_: *mut idx_t,
    cptr: *mut idx_t,
    cind: *mut idx_t,
) -> idx_t {
    // idx_t i, ii, j, jj, k, me=0, nvtxs, first, last, nleft, ncmps;
    // idx_t *xadj, *adjncy;
    // idx_t *touched, *perm, *todo;
    // idx_t mustfree_ccsr=0, mustfree_where_=0;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    adjncy = graph.adjncy;

    /* Deal with NULL supplied cptr/cind vectors */
    if (cptr == NULL) {
        cptr = imalloc(nvtxs + 1, "FindPartitionInducedComponents: cptr");
        cind = imalloc(nvtxs, "FindPartitionInducedComponents: cind");
        mustfree_ccsr = 1;
    }

    /* Deal with NULL supplied where_ vector */
    if (where_ == NULL) {
        where_ = ismalloc(nvtxs, 0, "FindPartitionInducedComponents: where_");
        mustfree_where_ = 1;
    }

    /* Allocate memory required for the BFS traversal */
    perm = iincset(
        nvtxs,
        0,
        imalloc(nvtxs, "FindPartitionInducedComponents: perm"),
    );
    todo = iincset(
        nvtxs,
        0,
        imalloc(nvtxs, "FindPartitionInducedComponents: todo"),
    );
    touched = ismalloc(nvtxs, 0, "FindPartitionInducedComponents: touched");

    /* Find the connected componends induced by the partition */
    ncmps = -1;
    first = last = 0;
    nleft = nvtxs;
    while (nleft > 0) {
        if (first == last) {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps] = first;
            ASSERT(touched[todo[0]] == 0);
            i = todo[0];
            cind[last] = i;
            last += 1;
            touched[i] = 1;
            me = where_[i];
        }

        i = cind[first];
        first += 1;
        k = perm[i];
        nleft -= 1;
        j = todo[k] = todo[nleft];
        perm[j] = k;

        for j in (xadj[i])..(xadj[i + 1]) {
            k = adjncy[j];
            if (where_[k] == me && !touched[k]) {
                cind[last] = k;
                last += 1;
                touched[k] = 1;
            }
        }
    }
    ncmps += 1;
    cptr[ncmps] = first;

    if (mustfree_ccsr) {
        panic!("gk_free((void **)&cptr, &cind, LTERM)");
    }
    if (mustfree_where_) {
        panic!("gk_free((void **)&where_, LTERM)");
    }

    panic!("gk_free((void **)&perm, &todo, &touched, LTERM)");

    return ncmps;
}

/*************************************************************************/
/* This function computes a permutation of the vertices based on a
    breadth-first-traversal. It can be used for re-ordering the graph
    to reduce its bandwidth for better cache locality.

    \param ctrl is the control structure
    \param graph is the graph structure
    \param perm is the array that upon completion, perm[i] will store
           the ID of the vertex that corresponds to the ith vertex in the
           re-ordered graph.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn ComputeBFSOrdering(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    bfsperm: *mut idx_t,
) -> void {
    // idx_t i, j, k, nvtxs, first, last;
    // idx_t *xadj, *adjncy, *perm;

    WCOREPUSH;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    adjncy = graph.adjncy;

    /* Allocate memory required for the BFS traversal */
    perm = iincset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));

    iincset(nvtxs, 0, bfsperm); /* this array will also store the vertices
                                still to be processed */

    /* Find the connected componends induced by the partition */
    first = last = 0;
    while (first < nvtxs) {
        if (first == last) {
            /* Find another starting vertex */
            k = bfsperm[last];
            ASSERT(perm[k] != -1);
            perm[k] = -1; /* mark node as being visited */
            last;
            last += 1;
        }

        i = bfsperm[first];
        first += 1;
        for j in (xadj[i])..(xadj[i + 1]) {
            k = adjncy[j];
            /* if a node has been already been visited, its perm[] will be -1 */
            if (perm[k] != -1) {
                /* perm[k] is the location within bfsperm of where_ k resides;
                put in that location bfsperm[last] that we are about to
                overwrite and update perm[bfsperm[last]] to reflect that. */
                bfsperm[perm[k]] = bfsperm[last];
                perm[bfsperm[last]] = perm[k];

                bfsperm[last] = k; /* put node at the end of the "queue" */
                last += 1;
                perm[k] = -1; /* mark node as being visited */
            }
        }
    }

    WCOREPOP;
}

/*************************************************************************/
/* This function checks whether a graph is contiguous or not.
 */
/**************************************************************************/
#[metis_func]
pub extern "C" fn IsConnected(graph: *mut graph_t, report: idx_t) -> idx_t {
    // idx_t ncmps;

    ncmps = FindPartitionInducedComponents(graph, NULL, NULL, NULL);

    if (ncmps != 1 && report) {
        printf(
            "The graph is not connected. It has {:} connected components.\n",
            ncmps,
        );
    }

    return (ncmps == 1);
}

/*************************************************************************/
/* This function checks whether or not partition pid is contiguous
 */
/*************************************************************************/
#[metis_func]
pub extern "C" fn IsConnectedSubdomain(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    pid: idx_t,
    report: idx_t,
) -> idx_t {
    // idx_t i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
    // idx_t *xadj, *adjncy, *where_, *touched, *queue;
    // idx_t *cptr;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    where_ = graph.where_;

    touched = ismalloc(nvtxs, 0, "IsConnected: touched");
    queue = imalloc(nvtxs, "IsConnected: queue");
    cptr = imalloc(nvtxs + 1, "IsConnected: cptr");

    nleft = 0;
    for i in (0)..(nvtxs) {
        if (where_[i] == pid) {
            nleft;
            nleft += 1;
        }
    }

    for i in (0)..(nvtxs) {
        if (where_[i] == pid) {
            break;
        }
    }

    touched[i] = 1;
    queue[0] = i;
    first = 0;
    last = 1;

    cptr[0] = 0; /* This actually points to queue */
    ncmps = 0;
    while (first != nleft) {
        if (first == last) {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps] = first;
            for i in (0)..(nvtxs) {
                if (where_[i] == pid && !touched[i]) {
                    break;
                }
            }
            queue[last] = i;
            last += 1;
            touched[i] = 1;
        }

        i = queue[first];
        first += 1;
        for j in (xadj[i])..(xadj[i + 1]) {
            k = adjncy[j];
            if (where_[k] == pid && !touched[k]) {
                queue[last] = k;
                last += 1;
                touched[k] = 1;
            }
        }
    }
    ncmps += 1;
    cptr[ncmps] = first;

    if (ncmps > 1 && report) {
        printf(
            "The graph has {:} connected components in partition {:}:\t",
            ncmps,
            pid,
        );
        for i in (0)..(ncmps) {
            wgt = 0;
            for j in (cptr[i])..(cptr[i + 1]) {
                wgt += graph.vwgt[queue[j]];
            }
            printf("[{:5} {:5}] ", cptr[i + 1] - cptr[i], wgt);
            /*
                  if (cptr[i+1]-cptr[i] == 1)
            {
            printf("[{:} {:}] ", queue[cptr[i]], xadj[queue[cptr[i]]+1]-xadj[queue[cptr[i]]]);
            }
                  */
        }
        printf("\n");
    }

    panic!("gk_free((void **)&touched, &queue, &cptr, LTERM)");

    // return (ncmps == 1 ? 1 : 0);
    if ncmps == 1 {
        1
    } else {
        0
    }
}

/*************************************************************************/
/* This function identifies the number of connected components in a graph
    that result after removing the vertices that belong to the vertex
    separator (i.e., graph.where_[i] == 2).
    The connected component memberships are returned in the CSR-style
    pair of arrays cptr, cind.
*/
/**************************************************************************/
#[metis_func]
pub extern "C" fn FindSepInducedComponents(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    cptr: *mut idx_t,
    cind: *mut idx_t,
) -> idx_t {
    // idx_t i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
    // idx_t *xadj, *adjncy, *where_, *touched, *queue;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    where_ = graph.where_;

    touched = ismalloc(nvtxs, 0, "IsConnected: queue");

    for i in (0)..(graph.nbnd) {
        touched[graph.bndind[i]] = 1;
    }

    queue = cind;

    nleft = 0;
    for i in (0)..(nvtxs) {
        if (where_[i] != 2) {
            nleft;
            nleft += 1;
        }
    }

    for i in (0)..(nvtxs) {
        if (where_[i] != 2) {
            break;
        }
    }

    touched[i] = 1;
    queue[0] = i;
    first = 0;
    last = 1;
    cptr[0] = 0; /* This actually points to queue */
    ncmps = 0;

    while (first != nleft) {
        if (first == last) {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps] = first;
            for i in (0)..(nvtxs) {
                if (!touched[i]) {
                    break;
                }
            }
            queue[last] = i;
            last += 1;
            touched[i] = 1;
        }

        i = queue[first];
        first += 1;
        for j in (xadj[i])..(xadj[i + 1]) {
            k = adjncy[j];
            if (!touched[k]) {
                queue[last] = k;
                last += 1;
                touched[k] = 1;
            }
        }
    }
    ncmps += 1;
    cptr[ncmps] = first;

    panic!("gk_free((void **)&touched, LTERM)");

    return ncmps;
}

/*************************************************************************/
/* This function finds all the connected components induced by the
partitioning vector in graph.where_ and tries to push them around to
remove some of them. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn EliminateComponents(ctrl: *mut ctrl_t, graph: *mut graph_t) -> void {
    // idx_t i, ii, j, jj, k, me, nparts, nvtxs, ncon, ncmps, other,
    //       ncand, target;
    // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *where_, *pwgts;
    // idx_t *cptr, *cind, *cpvec, *pcptr, *pcind, *cwhere_;
    // idx_t cid, bestcid, *cwgt, *bestcwgt;
    // idx_t ntodo, oldntodo, *todo;
    // rkv_t *cand;
    // real_t *tpwgts;
    // idx_t *vmarker=NULL, *pmarker=NULL, *modind=NULL;  /* volume specific work arrays */
    WCOREPUSH;

    nvtxs = graph.nvtxs;
    ncon = graph.ncon;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    vwgt = graph.vwgt;
    // adjwgt = (ctrl.objtype == METIS_OBJTYPE_VOL ? NULL : graph.adjwgt);
    adjwgt = if ctrl.objtype == METIS_OBJTYPE_VOL {
        std::ptr::null()
    } else {
        graph.adjwgt
    };

    where_ = graph.where_;
    pwgts = graph.pwgts;

    nparts = ctrl.nparts;
    tpwgts = ctrl.tpwgts;

    cptr = iwspacemalloc(ctrl, nvtxs + 1);
    cind = iwspacemalloc(ctrl, nvtxs);

    ncmps = FindPartitionInducedComponents(graph, where_, cptr, cind);

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_CONTIGINFO,
        printf(
            "I found {:} components, for this {:}-way partition\n",
            ncmps,
            nparts,
        ),
    );

    /* There are more components than partitions */
    if (ncmps > nparts) {
        cwgt = iwspacemalloc(ctrl, ncon);
        bestcwgt = iwspacemalloc(ctrl, ncon);
        cpvec = iwspacemalloc(ctrl, nparts);
        pcptr = iset(nparts + 1, 0, iwspacemalloc(ctrl, nparts + 1));
        pcind = iwspacemalloc(ctrl, ncmps);
        cwhere_ = iset(nvtxs, -1, iwspacemalloc(ctrl, nvtxs));
        todo = iwspacemalloc(ctrl, ncmps);
        // cand     = (rkv_t *)wspacemalloc(ctrl, nparts*sizeof(rkv_t));
        cand = wspacemalloc(ctrl, nparts * sizeof(rkv_t));

        if (ctrl.objtype == METIS_OBJTYPE_VOL) {
            /* Vol-refinement specific working arrays */
            modind = iwspacemalloc(ctrl, nvtxs);
            vmarker = iset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
            pmarker = iset(nparts, -1, iwspacemalloc(ctrl, nparts));
        }

        /* Get a CSR representation of the components-2-partitions mapping */
        for i in (0)..(ncmps) {
            pcptr[where_[cind[cptr[i]]]] += 1;
        }
        MAKECSR(i, nparts, pcptr);
        for i in (0)..(ncmps) {
            let p = &mut pcptr[where_[cind[cptr[i]]]];
            pcind[*p] = i;
            *p += 1;
        }
        SHIFTCSR(i, nparts, pcptr);

        /* Assign the heaviest component of each partition to its original partition */
        ntodo = 0;
        for i in (0)..(nparts) {
            if (pcptr[i + 1] - pcptr[i] == 1) {
                bestcid = pcind[pcptr[i]];
            } else {
                bestcid = -1;
                for j in (pcptr[i])..(pcptr[i + 1]) {
                    cid = pcind[j];
                    iset(ncon, 0, cwgt);
                    for ii in (cptr[cid])..(cptr[cid + 1]) {
                        iaxpy(ncon, 1, vwgt + cind[ii] * ncon, 1, cwgt, 1);
                    }
                    if (bestcid == -1 || isum(ncon, bestcwgt, 1) < isum(ncon, cwgt, 1)) {
                        bestcid = cid;
                        icopy(ncon, cwgt, bestcwgt);
                    }
                }
                /* Keep track of those that need to be dealt with */
                for j in (pcptr[i])..(pcptr[i + 1]) {
                    if (pcind[j] != bestcid) {
                        todo[ntodo] = pcind[j];
                        ntodo += 1;
                    }
                }
            }

            for j in (cptr[bestcid])..(cptr[bestcid + 1]) {
                ASSERT(where_[cind[j]] == i);
                cwhere_[cind[j]] = i;
            }
        }

        while (ntodo > 0) {
            oldntodo = ntodo;
            for i in (0)..(ntodo) {
                cid = todo[i];
                me = where_[cind[cptr[cid]]]; /* Get the domain of this component */

                /* Determine the weight of the block to be moved */
                iset(ncon, 0, cwgt);
                for j in (cptr[cid])..(cptr[cid + 1]) {
                    iaxpy(ncon, 1, vwgt + cind[j] * ncon, 1, cwgt, 1);
                }

                IFSET(
                    ctrl.dbglvl,
                    METIS_DBG_CONTIGINFO,
                    printf(
                        "Trying to move {:} [{:}] from {:}\n",
                        cid,
                        isum(ncon, cwgt, 1),
                        me,
                    ),
                );

                /* Determine the connectivity */
                iset(nparts, 0, cpvec);
                for j in (cptr[cid])..(cptr[cid + 1]) {
                    ii = cind[j];
                    for jj in (xadj[ii])..(xadj[ii + 1]) {
                        if (cwhere_[adjncy[jj]] != -1) {
                            cpvec[cwhere_[adjncy[jj]]] +=
                                (if adjwgt != 0 { adjwgt[jj] } else { 1 });
                        }
                    }
                }

                /* Put the neighbors into a cand[] array for sorting */
                ncand = 0;
                for j in (0)..(nparts) {
                    if (cpvec[j] > 0) {
                        cand[ncand].key = cpvec[j];
                        cand[ncand].val = j;
                        ncand += 1;
                    }
                }
                if (ncand == 0) {
                    continue;
                }

                rkvsortd(ncand, cand);

                /* Limit the moves to only the top candidates, which are defined as
                those with connectivity at least 50% of the best.
                This applies only when ncon=1, as for multi-constraint, balancing
                will be hard. */
                if (ncon == 1) {
                    for j in (1)..(ncand) {
                        if (cand[j].key < 0.5 * cand[0].key) {
                            break;
                        }
                    }
                    ncand = j;
                }

                /* Now among those, select the one with the best balance */
                target = cand[0].val;
                for j in (1)..(ncand) {
                    if (BetterBalanceKWay(
                        ncon,
                        cwgt,
                        ctrl.ubfactors,
                        1,
                        pwgts + target * ncon,
                        ctrl.pijbm + target * ncon,
                        1,
                        pwgts + cand[j].val * ncon,
                        ctrl.pijbm + cand[j].val * ncon,
                    )) {
                        target = cand[j].val;
                    }
                }

                IFSET(
                    ctrl.dbglvl,
                    METIS_DBG_CONTIGINFO,
                    printf(
                        "\tMoving it to {:} [{:}] [{:}]\n",
                        target,
                        cpvec[target],
                        ncand,
                    ),
                );

                /* Note that as a result of a previous movement, a connected component may
                now will like to stay to its original partition */
                if (target != me) {
                    match (ctrl.objtype) {
                        METIS_OBJTYPE_CUT => {
                            MoveGroupContigForCut(ctrl, graph, target, cid, cptr, cind);
                        }

                        METIS_OBJTYPE_VOL => {
                            MoveGroupContigForVol(
                                ctrl, graph, target, cid, cptr, cind, vmarker, pmarker, modind,
                            );
                        }

                        _ => gk_errexit(SIGERR, "Unknown objtype %d\n", ctrl.objtype),
                    }
                }

                /* Update the cwhere_ vector */
                for j in (cptr[cid])..(cptr[cid + 1]) {
                    cwhere_[cind[j]] = target;
                }

                ntodo -= 1;
                todo[i] = todo[ntodo];
            }
            if (oldntodo == ntodo) {
                IFSET(
                    ctrl.dbglvl,
                    METIS_DBG_CONTIGINFO,
                    printf("Stopped at ntodo: {:}\n", ntodo),
                );
                break;
            }
        }

        for i in (0)..(nvtxs) {
            ASSERT(where_[i] == cwhere_[i]);
        }
    }

    WCOREPOP;
}

/*************************************************************************/
/* This function moves a collection of vertices and updates their rinfo
 */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MoveGroupContigForCut(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    to: idx_t,
    gid: idx_t,
    ptr: *mut idx_t,
    ind: *mut idx_t,
) -> void {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, iii, j, jj, k, l, nvtxs, nbnd, from, me;
    // idx_t *xadj, *adjncy, *adjwgt, *where_, *bndptr, *bndind;
    // ckrinfo_t *myrinfo;
    // cnbr_t *mynbrs;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;

    where_ = graph.where_;
    bndptr = graph.bndptr;
    bndind = graph.bndind;

    nbnd = graph.nbnd;

    for iii in (ptr[gid])..(ptr[gid + 1]) {
        i = ind[iii];
        from = where_[i];

        myrinfo = graph.ckrinfo + i;
        if (myrinfo.inbr == -1) {
            myrinfo.inbr = cnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
            myrinfo.nnbrs = 0;
        }
        mynbrs = ctrl.cnbrpool + myrinfo.inbr;

        /* find the location of 'to' in myrinfo or create it if it is not there */
        for k in (0)..(myrinfo.nnbrs) {
            if (mynbrs[k].pid == to) {
                break;
            }
        }
        if (k == myrinfo.nnbrs) {
            mynbrs[k].pid = to;
            mynbrs[k].ed = 0;
            myrinfo.nnbrs;
            nnbrs += 1;
        }

        graph.mincut -= mynbrs[k].ed - myrinfo.id;

        /* Update ID/ED and BND related information for the moved vertex */
        iaxpy( graph.ncon, 1, graph.vwgt + i * graph.ncon, 1, graph.pwgts + to * graph.ncon, 1,);
        iaxpy( graph.ncon, -1, graph.vwgt + i * graph.ncon, 1, graph.pwgts + from * graph.ncon, 1,);
        UpdateMovedVertexInfoAndBND(
            i,
            from,
            k,
            to,
            myrinfo,
            mynbrs,
            where_,
            nbnd,
            bndptr,
            bndind,
            BNDTYPE_REFINE,
        );

        /* Update the degrees of adjacent vertices */
        for j in (xadj[i])..(xadj[i + 1]) {
            ii = adjncy[j];
            me = where_[ii];
            myrinfo = graph.ckrinfo + ii;

            UpdateAdjacentVertexInfoAndBND(
                ctrl,
                ii,
                xadj[ii + 1] - xadj[ii],
                me,
                from,
                to,
                myrinfo,
                adjwgt[j],
                nbnd,
                bndptr,
                bndind,
                BNDTYPE_REFINE,
            );
        }

        ASSERT(CheckRInfo(ctrl, graph.ckrinfo + i));
    }

    graph.nbnd = nbnd;
}

/*************************************************************************/
/* This function moves a collection of vertices and updates their rinfo
 */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MoveGroupContigForVol(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    to: idx_t,
    gid: idx_t,
    ptr: *mut idx_t,
    ind: *mut idx_t,
    vmarker: *mut idx_t,
    pmarker: *mut idx_t,
    modind: *mut idx_t,
) -> void {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, iii, j, jj, k, l, nvtxs, from, me, other, xgain;
    // idx_t *xadj, *vsize, *adjncy, *where_;
    // vkrinfo_t *myrinfo, *orinfo;
    // vnbr_t *mynbrs, *onbrs;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    vsize = graph.vsize;
    adjncy = graph.adjncy;
    where_ = graph.where_;

    for iii in (ptr[gid])..(ptr[gid + 1]) {
        i = ind[iii];
        from = where_[i];

        myrinfo = graph.vkrinfo + i;
        if (myrinfo.inbr == -1) {
            myrinfo.inbr = vnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
            myrinfo.nnbrs = 0;
        }
        mynbrs = ctrl.vnbrpool + myrinfo.inbr;

        xgain = (if myrinfo.nid == 0 && myrinfo.ned > 0 {
            vsize[i]
        } else {
            0
        });

        /* find the location of 'to' in myrinfo or create it if it is not there */
        for k in (0)..(myrinfo.nnbrs) {
            if (mynbrs[k].pid == to) {
                break;
            }
        }
        if (k == myrinfo.nnbrs) {
            if (myrinfo.nid > 0) {
                xgain -= vsize[i];
            }

            /* determine the volume gain resulting from that move */
            for j in (xadj[i])..(xadj[i + 1]) {
                ii = adjncy[j];
                other = where_[ii];
                orinfo = graph.vkrinfo + ii;
                onbrs = ctrl.vnbrpool + orinfo.inbr;
                ASSERT(other != to);

                if (from == other) {
                    /* Same subdomain vertex: Decrease the gain if 'to' is a new neighbor. */
                    for l in (0)..(orinfo.nnbrs) {
                        if (onbrs[l].pid == to) {
                            break;
                        }
                    }
                    if (l == orinfo.nnbrs) {
                        xgain -= vsize[ii];
                    }
                } else {
                    /* Remote vertex: increase if 'to' is a new subdomain */
                    for l in (0)..(orinfo.nnbrs) {
                        if (onbrs[l].pid == to) {
                            break;
                        }
                    }
                    if (l == orinfo.nnbrs) {
                        xgain -= vsize[ii];
                    }

                    /* Remote vertex: decrease if i is the only connection to 'from' */
                    for l in (0)..(orinfo.nnbrs) {
                        if (onbrs[l].pid == from && onbrs[l].ned == 1) {
                            xgain += vsize[ii];
                            break;
                        }
                    }
                }
            }
            graph.minvol -= xgain;
            graph.mincut -= -myrinfo.nid;
        } else {
            graph.minvol -= (xgain + mynbrs[k].gv);
            graph.mincut -= mynbrs[k].ned - myrinfo.nid;
        }

        /* Update where_ and pwgts */
        where_[i] = to;
        iaxpy(
            graph.ncon,
            1,
            graph.vwgt + i * graph.ncon,
            1,
            graph.pwgts + to * graph.ncon,
            1,
        );
        iaxpy(
            graph.ncon,
            -1,
            graph.vwgt + i * graph.ncon,
            1,
            graph.pwgts + from * graph.ncon,
            1,
        );

        /* Update the id/ed/gains/bnd of potentially affected nodes */
        KWayVolUpdate(
            ctrl,
            graph,
            i,
            from,
            to,
            NULL,
            NULL,
            NULL,
            NULL,
            NULL,
            BNDTYPE_REFINE,
            vmarker,
            pmarker,
            modind,
        );

        /*CheckKWayVolPartitionParams(ctrl, graph);*/
    }

    ASSERT(ComputeCut(graph, where_) == graph.mincut);
    ASSERTP(
        ComputeVolume(graph, where_) == graph.minvol,
        ("{:} {:}\n", ComputeVolume(graph, where_), graph.minvol),
    );
}
