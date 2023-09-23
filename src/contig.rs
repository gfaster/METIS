/*
\file
\brief Functions that deal with eliminating disconnected partitions

\date Started 7/15/98
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: contig.c 10513 2011-07-07 22:06:03Z karypis $
*/

use core::{slice, ptr};

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
    where_: *const idx_t,
    cptr: *mut idx_t,
    cind: *mut idx_t,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    // idx_t i, ii, j, jj, k, me=0, nvtxs, first, last, nleft, ncmps;
    // idx_t *xadj, *adjncy;
    // idx_t *touched, *perm, *todo;
    // idx_t mustfree_ccsr=0, mustfree_where_=0;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy);
    // xadj = graph.xadj;
    // adjncy = graph.adjncy;

    /* Deal with NULL supplied cptr/cind vectors */
    // if (cptr == NULL) {
    //     cptr = imalloc(nvtxs + 1, "FindPartitionInducedComponents: cptr");
    //     cind = imalloc(nvtxs, "FindPartitionInducedComponents: cind");
    //     mustfree_ccsr = 1;
    // }
    let mut cptr_v: Option<Vec<idx_t>> = None;
    let cptr = if cptr.is_null() {
        cptr_v = Some(vec![0; nvtxs + 1]);
        &mut cptr_v.unwrap()[..]
    } else {
        slice::from_raw_parts_mut(cptr, nvtxs + 1)
    };
    let mut cind_v: Option<Vec<idx_t>> = None;
    let cind = if cind.is_null() {
        cind_v = Some(vec![0; nvtxs]);
        &mut cind_v.unwrap()[..]
    } else {
        slice::from_raw_parts_mut(cind, nvtxs)
    };

    /* Deal with NULL supplied where_ vector */
    // if (where_ == NULL) {
    //     where_ = ismalloc(nvtxs, 0, "FindPartitionInducedComponents: where_");
    //     mustfree_where_ = 1;
    // }
    let mut where_v: Option<Vec<idx_t>> = None;
    let where_ = if where_.is_null() {
        where_v = Some(vec![0; nvtxs]);
        &where_v.unwrap()[..]
    } else {
        slice::from_raw_parts(where_, nvtxs)
    };

    /* Allocate memory required for the BFS traversal */
    // perm = iincet( nvtxs, 0, imalloc(nvtxs, "FindPartitionInducedComponents: perm"),);
    let perm = Vec::from_iter(0..(nvtxs as usize));
    // todo = iincset( nvtxs, 0, imalloc(nvtxs, "FindPartitionInducedComponents: todo"),);
    let todo = Vec::from_iter(0..(nvtxs as usize));
    // touched = ismalloc(nvtxs, 0, "FindPartitionInducedComponents: touched");
    let touched = vec![0; nvtxs];

    /* Find the connected componends induced by the partition */
    let mut ncmps: idx_t = -1;
    let mut first= 0;
    let mut last= 0;
    let mut me = 0;
    let mut nleft = nvtxs;
    while (nleft > 0) {
        if (first == last) {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps as usize] = first;
            ASSERT(touched[todo[0]] == 0);
            let i = todo[0] as usize;
            cind[last as usize] = i as idx_t;
            last += 1;
            touched[i] = 1;
            me = where_[i];
        }

        let i = cind[first as usize] as usize;
        first += 1;
        let k = perm[i];
        nleft -= 1;
        let j = todo[nleft];
        todo[k as usize] = j;
        perm[j] = k;

        for j in (xadj[i])..(xadj[i + 1]) {
            k = adjncy[j as usize] as usize;
            if (where_[k] == me && touched[k] == 0) {
                cind[last as usize] = k as idx_t;
                last += 1;
                touched[k] = 1;
            }
        }
    }
    ncmps += 1;
    cptr[ncmps as usize] = first;

    // if (mustfree_ccsr) {
    //     panic!("gk_free((void **)&cptr, &cind, LTERM)");
    // }
    // if (mustfree_where_) {
    //     panic!("gk_free((void **)&where_, LTERM)");
    // }
    //
    // panic!("gk_free((void **)&perm, &todo, &touched, LTERM)");

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
) -> () {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, k, nvtxs, first, last;
    // idx_t *xadj, *adjncy, *perm;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy);
    // xadj = graph.xadj;
    // adjncy = graph.adjncy;

    /* Allocate memory required for the BFS traversal */
    // perm = iincset(nvtxs, 0, vec![0; nvtxs) as usize];
    let mut perm = Vec::from_iter(0..(nvtxs as idx_t));


    mkslice_mut!(bfsperm, nvtxs);
    // iincset(nvtxs, 0, bfsperm); /* this array will also store the vertices
    //                             still to be processed */
    bfsperm.iter_mut().enumerate().map(|(i, v)| *v = i as idx_t);


    /* Find the connected componends induced by the partition */
    let mut first = 0;
    let mut last = 0;
    while (first < nvtxs) {
        if (first == last) {
            /* Find another starting vertex */
            let k = bfsperm[last] as usize;
            ASSERT(perm[k] != -1);
            perm[k] = -1; /* mark node as being visited */
            last += 1;
        }

        let i = bfsperm[first] as usize;
        first += 1;
        for j in (xadj[i])..(xadj[i + 1]) {
            let k = adjncy[j as usize] as usize;
            /* if a node has been already been visited, its perm[] will be -1 */
            if (perm[k] != -1) {
                /* perm[k] is the location within bfsperm of where_ k resides;
                put in that location bfsperm[last] that we are about to
                overwrite and update perm[bfsperm[last]] to reflect that. */
                bfsperm[perm[k] as usize] = bfsperm[last];
                perm[bfsperm[last] as usize] = perm[k];

                bfsperm[last] = k as idx_t; /* put node at the end of the "queue" */
                last += 1;
                perm[k] = -1; /* mark node as being visited */
            }
        }
    }

    // WCOREPOP;
}

/*************************************************************************/
/* This function checks whether a graph is contiguous or not.
 */
/**************************************************************************/
#[metis_func]
pub extern "C" fn IsConnected(graph: *mut graph_t, report: idx_t) -> idx_t {
    // idx_t ncmps;

    let ncmps = FindPartitionInducedComponents(graph, ptr::null_mut(), ptr::null_mut(), ptr::null_mut());

    if (ncmps != 1 && report != 0) {
        printf(
            "The graph is not connected. It has {:} connected components.\n",
            ncmps,
        );
    }

    return (ncmps == 1) as idx_t;
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
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
    // idx_t *xadj, *adjncy, *where_, *touched, *queue;
    // idx_t *cptr;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy where_);

    // touched = ismalloc(nvtxs, 0, "IsConnected: touched");
    // queue = imalloc(nvtxs, "IsConnected: queue");
    // cptr = imalloc(nvtxs + 1, "IsConnected: cptr");
    let mut touched = vec![0; nvtxs];
    let mut queue = vec![0; nvtxs];
    let mut cptr = vec![0; nvtxs + 1];

    let mut nleft = 0;
    for i in (0)..(nvtxs) {
        if (where_[i] == pid) {
            nleft += 1;
        }
    }

    let mut i = 0;
    for ii in (0)..(nvtxs) {
        i = ii;
        if (where_[i] == pid) {
            break;
        }
    }

    touched[i as usize] = 1;
    queue[0] = i;
    let mut first = 0;
    let mut last = 1;

    cptr[0] = 0; /* This actually points to queue */
    let mut ncmps = 0;
    while (first != nleft) {
        if (first == last) {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps] = first;
            for i in (0)..(nvtxs) {
                if (where_[i] == pid && touched[i] == 0) {
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
            let k = adjncy[j as usize] as usize;
            if (where_[k] == pid && touched[k] == 0) {
                queue[last] = k;
                last += 1;
                touched[k] = 1;
            }
        }
    }
    ncmps += 1;
    cptr[ncmps] = first;

    if (ncmps > 1 && report != 0) {
        printf(
            "The graph has {:} connected components in partition {:}:\t",
            ncmps,
            pid,
        );
        get_graph_slices!(graph => vwgt);
        for i in (0)..(ncmps) {
            let mut wgt = 0;
            for j in (cptr[i])..(cptr[i + 1]) {
                wgt += vwgt[queue[j]];
            }
            print!("[{:5} {:5}] ", cptr[i + 1] - cptr[i], wgt);
            /*
                  if (cptr[i+1]-cptr[i] == 1)
            {
            printf("[{:} {:}] ", queue[cptr[i]], xadj[queue[cptr[i]]+1]-xadj[queue[cptr[i]]]);
            }
                  */
        }
        printf("\n");
    }

    // panic!("gk_free((void **)&touched, &queue, &cptr, LTERM)");

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
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, k, nvtxs, first, last, nleft, ncmps, wgt;
    // idx_t *xadj, *adjncy, *where_, *touched, *queue;

    let nvtxs = graph.nvtxs as usize;
    // xadj = graph.xadj;
    // adjncy = graph.adjncy;
    // where_ = graph.where_;
    get_graph_slices!(graph => xadj adjncy where_);

    // touched = ismalloc(nvtxs, 0, "IsConnected: queue");
    let mut touched = vec![0; nvtxs];

    {
        get_graph_slices!(graph => bndind);
        for i in (0)..(graph.nbnd) {
            touched[bndind[i as usize] as usize] = 1;
        }
    }

    // queue = cind;
    mkslice!(queue: cind, nvtxs);
    mkslice!(cptr, nvtxs);

    let mut nleft = 0;
    for i in (0)..(nvtxs) {
        if (where_[i] != 2) {
            nleft;
            nleft += 1;
        }
    }

    {
    let mut i = 0;
    for ii in (0)..(nvtxs) {
        i = ii;
        if (where_[i] != 2) {
            break;
        }
    }
    queue[0] = i as i32;
    touched[i] = 1;
}

    let mut first = 0;
    let mut last = 1;
    cptr[0] = 0; /* This actually points to queue */
    let mut ncmps = 0;

    while (first != nleft) {
        if (first == last) {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps] = first;
            let mut i = 0;
            for ii in (0)..(nvtxs) {
                i = ii;
                if (touched[i] == 0) {
                    break;
                }
            }
            queue[last as usize] = i as idx_t;
            last += 1;
            touched[i] = 1;
        }

        let i = queue[first as usize] as usize;
        first += 1;
        for j in (xadj[i])..(xadj[i + 1]) {
            let k = adjncy[j as usize] as usize;
            if (touched[k] == 0) {
                queue[last as usize] = k as idx_t;
                last += 1;
                touched[k] = 1;
            }
        }
    }
    ncmps += 1;
    cptr[ncmps] = first;

    // panic!("gk_free((void **)&touched, LTERM)");

    return ncmps as idx_t;
}

/*************************************************************************/
/* This function finds all the connected components induced by the
partitioning vector in graph.where_ and tries to push them around to
remove some of them. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn EliminateComponents(ctrl: *mut ctrl_t, graph: *mut graph_t) -> () {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, jj, k, me, nparts, nvtxs, ncon, ncmps, other,
    //       ncand, target;
    // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *where_, *pwgts;
    // idx_t *cptr, *cind, *cpvec, *pcptr, *pcind, *cwhere_;
    // idx_t cid, bestcid, *cwgt, *bestcwgt;
    // idx_t ntodo, oldntodo, *todo;
    // rkv_t *cand;
    // real_t *tpwgts;
    // idx_t *vmarker=NULL, *pmarker=NULL, *modind=NULL;  /* volume specific work arrays */
    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    // xadj = graph.xadj;
    // adjncy = graph.adjncy;
    // vwgt = graph.vwgt;
    get_graph_slices!(ctrl, graph=> xadj adjncy vwgt where_ pwgts);

    // adjwgt = (ctrl.objtype == METIS_OBJTYPE_VOL ? NULL : graph.adjwgt);
    let adjwgt = if ctrl.objtype == METIS_OBJTYPE_VOL {
        std::ptr::null()
    } else {
        graph.adjwgt
    };
    #[derive(Clone, Copy)]
    struct KeyVal {
        key: real_t,
        val: idx_t,
    }


    // where_ = graph.where_;
    // pwgts = graph.pwgts;

    let nparts = ctrl.nparts as usize;
    let tpwgts = ctrl.tpwgts;

    // cptr = vec![0; nvtxs + 1 as usize];
    // cind = vec![0; nvtxs as usize];
    let mut cptr = vec![0; nvtxs + 1];
    let mut cind = vec![0; nvtxs];

    let ncmps = FindPartitionInducedComponents(graph, where_.as_ptr(), cptr.as_mut_ptr(), cind.as_mut_ptr()) as usize;

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
        let mut cwgt = vec![0; ncon];
        let mut bestcwgt = vec![0; ncon];
        let mut cpvec = vec![0; nparts];
        let mut pcptr = vec![0; (nparts + 1)];
        let mut pcind: Vec<idx_t> = vec![0; ncmps];
        let mut cwhere_: Vec<idx_t> = vec![-1; nvtxs];
        let mut todo = vec![0; ncmps];
        // cand     = (rkv_t *)wspacemalloc(ctrl, nparts*sizeof(rkv_t));
        let mut cand = vec![KeyVal { key: 0.0, val: 0}; nparts];

        let mut modind = Vec::new();
        let mut vmarker = Vec::new();
        let mut pmarker = Vec::new();
        if (ctrl.objtype == METIS_OBJTYPE_VOL) {
            /* Vol-refinement specific working arrays */
            modind = vec![0; nvtxs];
            vmarker = vec![0; nvtxs];
            pmarker = vec![-1; nparts];
        }

        /* Get a CSR representation of the components-2-partitions mapping */
        for i in (0)..(ncmps) {
            pcptr[where_[cind[cptr[i] as usize] as usize] as usize] += 1;
        }
        MAKECSR(i, nparts, pcptr);
        for i in (0)..(ncmps) {
            let p = &mut pcptr[where_[cind[cptr[i] as usize] as usize] as usize];
            pcind[*p] = i as idx_t;
            *p += 1;
        }
        SHIFTCSR(i, nparts, pcptr);

        /* Assign the heaviest component of each partition to its original partition */
        let mut ntodo = 0;
        for i in (0)..(nparts) {
            let bestcid: idx_t;
            if (pcptr[i + 1] - pcptr[i] == 1) {
                bestcid = pcind[pcptr[i] as usize];
            } else {
                bestcid = -1;
                for j in (pcptr[i])..(pcptr[i + 1]) {
                    let cid = pcind[j as usize];
                    iset(ncon, 0, cwgt);
                    for ii in (cptr[cid as usize])..(cptr[cid as usize + 1]) {
                        blas::iaxpy(ncon, 1, &vwgt[(cind[ii as usize] as usize * ncon)..], 1, &mut cwgt[..], 1);
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

            for j in (cptr[bestcid as usize])..(cptr[bestcid as usize + 1]) {
                ASSERT(where_[cind[j as usize] as usize] == i as idx_t);
                cwhere_[cind[j as usize] as usize] = i as idx_t;
            }
        }

        while (ntodo > 0) {
            let oldntodo = ntodo;
            for i in (0)..(ntodo) {
                let cid = todo[i];
                let me = where_[cind[cptr[cid as usize] as usize] as usize]; /* Get the domain of this component */

                /* Determine the weight of the block to be moved */
                iset(ncon, 0, cwgt);
                for j in (cptr[cid as usize])..(cptr[cid as usize + 1]) {
                    blas::iaxpy(ncon, 1, &vwgt[(cind[j as usize] as usize * ncon)..], 1, &mut cwgt[..], 1);
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
                for j in (cptr[cid as usize])..(cptr[cid as usize + 1]) {
                    let ii = cind[j as usize] as usize;
                    for jj in (xadj[ii])..(xadj[ii + 1]) {
                        let jj = jj as usize;
                        if (cwhere_[adjncy[jj] as usize] != -1) {
                            cpvec[cwhere_[adjncy[jj] as usize] as usize] +=
                        {
                                // mkslice!(adjwgt, graph.nedges);
                                (if !adjwgt.is_null() { *adjwgt.add(jj) } else { 1 })
                            }
                        }
                    }
                }

                /* Put the neighbors into a cand[] array for sorting */
                let mut ncand = 0;
                for j in (0)..(nparts) {
                    if (cpvec[j] > 0) {
                        cand[ncand].key = cpvec[j] as real_t;
                        cand[ncand].val = j as idx_t;
                        ncand += 1;
                    }
                }
                if (ncand == 0) {
                    continue;
                }

                // rkvsortd(ncand, cand);
                cand.sort_unstable_by(|l, r| l.key.partial_cmp(&r.key).unwrap());

                /* Limit the moves to only the top candidates, which are defined as
                those with connectivity at least 50% of the best.
                This applies only when ncon=1, as for multi-constraint, balancing
                will be hard. */
                if (ncon == 1) {
                    let mut j = 0;
                    for jj in (1)..(ncand) {
                        j = jj;
                        if (cand[j].key < 0.5 * cand[0].key) {
                            break;
                        }
                    }
                    ncand = j;
                }

                /* Now among those, select the one with the best balance */
                let target = cand[0].val as usize;
                for j in (1)..(ncand) {
                    // ctrl.pijbm is split, but shouldn't actually alias. I'll have to fix this for
                    // real later
                    if (BetterBalanceKWay(
                        ncon,
                        cwgt,
                        ctrl.ubfactors,
                        1,
                        &pwgts[target * ncon..],
                        ctrl.pijbm.add(target * ncon),
                        1,
                        &pwgts[ cand[j].val as usize * ncon..],
                        ctrl.pijbm.add( cand[j].val as usize * ncon),
                    ) != 0) {
                        target = cand[j].val as usize;
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
                if (target != me as usize) {
                    match (ctrl.objtype) {
                        METIS_OBJTYPE_CUT => {
                            MoveGroupContigForCut(ctrl, graph, target, cid, cptr, cind);
                        }

                        METIS_OBJTYPE_VOL => {
                            MoveGroupContigForVol(
                                ctrl, graph, target, cid, cptr, cind, vmarker, pmarker, modind,
                            );
                        }

                        _ => panic!("Unknown objtype {}", ctrl.objtype),
                    }
                }

                /* Update the cwhere_ vector */
                for j in (cptr[cid as usize])..(cptr[cid as usize + 1]) {
                    cwhere_[cind[j as usize] as usize] = target as idx_t;
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

    // WCOREPOP;
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
){
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, iii, j, jj, k, l, nvtxs, nbnd, from, me;
    // idx_t *xadj, *adjncy, *adjwgt, *where_, *bndptr, *bndind;
    // ckrinfo_t *myrinfo;
    // cnbr_t *mynbrs;

    let nvtxs = graph.nvtxs as usize;
    // xadj = graph.xadj;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    get_graph_slices!(graph => xadj adjncy adjwgt);

    // where_ = graph.where_;
    // bndptr = graph.bndptr;
    // bndind = graph.bndind;
    get_graph_slices_mut!(graph => where_ bndptr bndind);

    mkslice!(ptr, nvtxs);
    mkslice!(ind, nvtxs);
    let gid = gid as usize;

    let nbnd = graph.nbnd as usize;

    for iii in (ptr[gid])..(ptr[gid + 1]) {
        let i = ind[iii as usize] as usize;
        let from = where_[i];

        let myrinfo = graph.ckrinfo.add(i as usize).as_mut().unwrap();
        if (myrinfo.inbr == -1) {
            myrinfo.inbr = cnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
            myrinfo.nnbrs = 0;
        }
        let mynbrs = std::slice::from_raw_parts_mut(ctrl.cnbrpool.add( myrinfo.inbr as usize), myrinfo.nnbrs as usize + 1);

        /* find the location of 'to' in myrinfo or create it if it is not there */
        let mut k = 0;
        for kk in (0)..(myrinfo.nnbrs) {
            k = kk as usize;
            if (mynbrs[k].pid == to) {
                break;
            }
        }
        if (k == myrinfo.nnbrs as usize) {
            mynbrs[k].pid = to;
            mynbrs[k].ed = 0;
            myrinfo.nnbrs+=1;
        }

        graph.mincut -= mynbrs[k].ed - myrinfo.id;

        /* Update ID/ED and BND related information for the moved vertex */
        {
            get_graph_slices!(graph => vwgt);
            get_graph_slices_mut!(ctrl, graph => pwgts);
            let ncon = graph.ncon;
            blas::iaxpy( ncon as usize, 1, &vwgt [cntrng!(i * ncon, ncon)], 1, &mut pwgts[cntrng!(to * ncon, ncon)], 1,);
            blas::iaxpy( ncon as usize, -1, &vwgt [cntrng!(i * ncon, ncon)], 1,&mut pwgts[cntrng!(to * ncon, ncon)], 1,);
        }
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
        for j in (xadj[i as usize])..(xadj[i as usize + 1]) {
            let j = j as usize;
            let ii = adjncy[j] as usize;
            let me = where_[ii];
            let myrinfo = graph.ckrinfo.add(ii);

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

        ASSERT(CheckRInfo(ctrl, graph.ckrinfo.add(i as usize)));
    }

    graph.nbnd = nbnd as idx_t;
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

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj vsize adjncy where_);
    // xadj = graph.xadj;
    // vsize = graph.vsize;
    // adjncy = graph.adjncy;
    // where_ = graph.where_;

    mkslice!(ptr, nvtxs);
    mkslice!(ind, nvtxs);
    let gid = gid as usize;

    for iii in (ptr[gid])..(ptr[gid + 1]) {
        let i = ind[iii as usize] as usize;
        let from = where_[i];

        let myrinfo = graph.vkrinfo.add(i).as_mut().unwrap();
        if (myrinfo.inbr == -1) {
            myrinfo.inbr = vnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
            myrinfo.nnbrs = 0;
        }
        let mynbrs = slice::from_raw_parts_mut(ctrl.vnbrpool.add( myrinfo.inbr as usize), myrinfo.nnbrs as usize + 1);

        let xgain = (if myrinfo.nid == 0 && myrinfo.ned > 0 {
            vsize[i]
        } else {
            0
        });

        /* find the location of 'to' in myrinfo or create it if it is not there */
        let mut k = 0;
        for kk in (0)..(myrinfo.nnbrs as usize) {
            k = kk;
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
                let ii = adjncy[j as usize] as usize;
                let other = where_[ii];
                let orinfo = *graph.vkrinfo.add(ii);
                let onbrs = slice::from_raw_parts_mut(ctrl.vnbrpool.add(orinfo.inbr as usize), orinfo.nnbrs as usize + 1);
                ASSERT(other != to);

                if (from == other) {
                    /* Same subdomain vertex: Decrease the gain if 'to' is a new neighbor. */

                    let mut l = 0;
                    for ll in (0)..(orinfo.nnbrs) {
                        l = ll;
                        if (onbrs[l as usize].pid == to) {
                            break;
                        }
                    }
                    if (l == orinfo.nnbrs) {
                        xgain -= vsize[ii];
                    }
                } else {
                    /* Remote vertex: increase if 'to' is a new subdomain */
                    {
                    let mut l = 0;
                    for ll in (0)..(orinfo.nnbrs) {
                        l = ll;
                        if (onbrs[l as usize].pid == to) {
                            break;
                        }
                    }
                    if (l == orinfo.nnbrs) {
                        xgain -= vsize[ii];
                    }
                    }

                    /* Remote vertex: decrease if i is the only connection to 'from' */
                   for l in (0)..(orinfo.nnbrs as usize) {
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
        // iaxpy( graph.ncon, 1, graph.vwgt + i * graph.ncon, 1, graph.pwgts + to * graph.ncon, 1,);
        // iaxpy( graph.ncon, -1, graph.vwgt + i * graph.ncon, 1, graph.pwgts + from * graph.ncon, 1,);
        {
            get_graph_slices!(graph => vwgt);
            get_graph_slices_mut!(ctrl, graph => pwgts);
            let ncon = graph.ncon;
            blas::iaxpy( ncon as usize, 1, &vwgt [cntrng!(i as idx_t * ncon, ncon)], 1, &mut pwgts[cntrng!(to * ncon, ncon)], 1,);
            blas::iaxpy( ncon as usize, -1, &vwgt [cntrng!(i as idx_t * ncon, ncon)], 1,&mut pwgts[cntrng!(to * ncon, ncon)], 1,);
        }

        /* Update the id/ed/gains/bnd of potentially affected nodes */
        KWayVolUpdate(
            ctrl,
            graph,
            i,
            from,
            to,
            ptr::null(),
            ptr::null(),
            ptr::null(),
            ptr::null(),
            ptr::null(),
            BNDTYPE_REFINE,
            vmarker,
            pmarker,
            modind,
        );

        /*CheckKWayVolPartitionParams(ctrl, graph);*/
    }

    ASSERT(ComputeCut(graph, where_.as_ptr()) == graph.mincut);
    ASSERTP(
        ComputeVolume(graph, where_.as_ptr()) == graph.minvol,
        ("{:} {:}\n", ComputeVolume(graph, where_.as_ptr()), graph.minvol),
    );
}
