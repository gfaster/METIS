/*
\file
\brief Functions that deal with eliminating disconnected partitions

\date Started 7/15/98
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: contig.c 10513 2011-07-07 22:06:03Z karypis $
*/

use core::{ptr, slice};
use std::cmp::Reverse;

use crate::*;
use kwayfm::{crinfos, crinfos_mut, vrinfos, vrinfos_mut};

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

    slice_default!(mut cptr, [0; nvtxs + 1]);
    slice_default!(mut cind, [0; nvtxs]);
    slice_default!(where_, [0; nvtxs]);

    /* Allocate memory required for the BFS traversal */
    // perm = iincet( nvtxs, 0, imalloc(nvtxs, "FindPartitionInducedComponents: perm"),);
    let mut perm = Vec::from_iter(0..(nvtxs as idx_t));
    // todo = iincset( nvtxs, 0, imalloc(nvtxs, "FindPartitionInducedComponents: todo"),);
    let mut todo = Vec::from_iter(0..(nvtxs as idx_t));
    // touched = ismalloc(nvtxs, 0, "FindPartitionInducedComponents: touched");
    let mut touched = vec![false; nvtxs];

    /* Find the connected componends induced by the partition */
    let mut ncmps: idx_t = -1;
    let mut first = 0;
    let mut last = 0;
    let mut me = 0;
    let mut nleft = nvtxs;
    while nleft > 0 {
        if first == last {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps as usize] = first;
            debug_assert!(!touched[todo[0] as usize]);
            let i = todo[0] as usize;
            cind[last as usize] = i as idx_t;
            last += 1;
            touched[i] = true;
            me = where_[i];
        }

        let i = cind[first as usize] as usize;
        first += 1;
        let k = perm[i];
        nleft -= 1;
        let j = todo[nleft];
        todo[k as usize] = j;
        perm[j as usize] = k;

        for j in (xadj[i])..(xadj[i + 1]) {
            let k = adjncy[j as usize] as usize;
            if where_[k] == me && !touched[k] {
                cind[last as usize] = k as idx_t;
                last += 1;
                touched[k] = true;
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

/// Computes a permutation of the vertices based on a
/// breadth-first-traversal. It can be used for re-ordering the graph
/// to reduce its bandwidth for better cache locality. It is unused by METIS.
///
/// - `ctrl` is the control structure
/// - `graph` is the graph structure
/// - `perm` is the array that upon completion, perm[i] will store the ID of the vertex that
///    corresponds to the ith vertex in the re-ordered graph.
#[metis_func]
pub extern "C" fn ComputeBFSOrdering(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    bfsperm: *mut idx_t,
) -> () {
    let graph = graph.as_mut().unwrap();
    // let ctrl = ctrl.as_mut().unwrap();
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
    bfsperm
        .iter_mut()
        .enumerate()
        .map(|(i, v)| *v = i as idx_t)
        .last();

    /* Find the connected componends induced by the partition */
    let mut first = 0;
    let mut last = 0;
    while first < nvtxs {
        if first == last {
            /* Find another starting vertex */
            let k = bfsperm[last] as usize;
            assert!(perm[k] != -1);
            perm[k] = -1; /* mark node as being visited */
            last += 1;
        }

        let i = bfsperm[first] as usize;
        first += 1;
        for j in (xadj[i])..(xadj[i + 1]) {
            let k = adjncy[j as usize] as usize;
            /* if a node has been already been visited, its perm[] will be -1 */
            if perm[k] != -1 {
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

/// This function checks whether a graph is contiguous or not.
#[metis_func]
pub extern "C" fn IsConnected(graph: *mut graph_t, report: idx_t) -> idx_t {
    // idx_t ncmps;

    let ncmps =
        FindPartitionInducedComponents(graph, ptr::null_mut(), ptr::null_mut(), ptr::null_mut());

    if ncmps != 1 && report != 0 {
        println!(
            "The graph is not connected. It has {:} connected components.\n",
            ncmps,
        );
    }

    return (ncmps == 1) as idx_t;
}

/// This function checks whether or not partition pid is contiguous
///
/// Unused in METIS
#[metis_func]
pub extern "C" fn IsConnectedSubdomain(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    pid: idx_t,
    report: idx_t,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    // let ctrl = ctrl.as_mut().unwrap();
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
        if where_[i] == pid {
            nleft += 1;
        }
    }

    {
        let mut i = 0;
        while i < nvtxs {
            if where_[i] == pid {
                break;
            }
            i += 1;
        }
        touched[i as usize] = 1;
        queue[0] = i;
    }

    let mut first = 0;
    let mut last = 1;

    cptr[0] = 0; /* This actually points to queue */
    let mut ncmps = 0;
    while first != nleft {
        if first == last {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps] = first;
            let mut i = 0;
            while i < nvtxs {
                if where_[i] == pid && touched[i] == 0 {
                    break;
                }
                i += 1;
            }
            queue[last] = i;
            last += 1;
            touched[i] = 1;
        }

        let i = queue[first];
        first += 1;
        for j in (xadj[i])..(xadj[i + 1]) {
            let k = adjncy[j as usize] as usize;
            if where_[k] == pid && touched[k] == 0 {
                queue[last] = k;
                last += 1;
                touched[k] = 1;
            }
        }
    }
    ncmps += 1;
    cptr[ncmps] = first;

    if ncmps > 1 && report != 0 {
        println!(
            "The graph has {:} connected components in partition {:}:\t",
            ncmps, pid,
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
            println!("[{:} {:}] ", queue[cptr[i]], xadj[queue[cptr[i]]+1]-xadj[queue[cptr[i]]]);
            }
                  */
        }
        println!("\n");
    }

    // return (ncmps == 1 ? 1 : 0);
    if ncmps == 1 { 1 } else { 0 }
}

/// This function identifies the number of connected components in a graph
/// that result after removing the vertices that belong to the vertex
/// separator (i.e., graph.where_[i] == 2).
///
/// The connected component memberships are returned in the CSR-style
/// pair of arrays cptr, cind.
#[metis_func]
pub extern "C" fn FindSepInducedComponents(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    cptr: *mut idx_t,
    cind: *mut idx_t,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
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
    mkslice_mut!(queue: cind, nvtxs);
    mkslice_mut!(cptr, nvtxs);

    let mut nleft = 0;
    for i in (0)..(nvtxs) {
        if where_[i] != 2 {
            nleft += 1;
        }
    }

    {
        let mut i = 0;
        for ii in (0)..(nvtxs) {
            i = ii;
            if where_[i] != 2 {
                break;
            }
        }
        queue[0] = i as idx_t;
        touched[i] = 1;
    }

    let mut first = 0;
    let mut last = 1;
    cptr[0] = 0; /* This actually points to queue */
    let mut ncmps = 0;

    while first != nleft {
        if first == last {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps] = first;
            let mut i = 0;
            for ii in (0)..(nvtxs) {
                i = ii;
                if touched[i] == 0 {
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
            if touched[k] == 0 {
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

/// This function finds all the connected components induced by the partitioning vector in
/// `graph.where_` and tries to push them around to remove some of them.
#[metis_func]
pub extern "C" fn EliminateComponents(ctrl: *mut ctrl_t, graph: *mut graph_t) -> () {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // ctrl.dbglvl |= METIS_DBG_CONTIGINFO;
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
    get_graph_slices!(ctrl, graph => xadj adjncy vwgt where_ pwgts);

    // adjwgt = (ctrl.objtype == METIS_OBJTYPE_VOL ? NULL : graph.adjwgt);
    let adjwgt = if ctrl.objtype == METIS_OBJTYPE_VOL {
        None
    } else {
        Some(slice::from_raw_parts(graph.adjwgt, graph.nedges as usize))
    };

    #[derive(Clone, Copy)]
    struct KeyVal {
        key: real_t,
        val: idx_t,
    }

    // where_ = graph.where_;
    // pwgts = graph.pwgts;

    let nparts = ctrl.nparts as usize;
    // unused in original
    // let tpwgts = ctrl.tpwgts;

    let mut cptr = vec![0; nvtxs + 1];
    let mut cind = vec![0; nvtxs];

    let ncmps = FindPartitionInducedComponents(
        graph,
        where_.as_ptr(),
        cptr.as_mut_ptr(),
        cind.as_mut_ptr(),
    ) as usize;

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_CONTIGINFO,
        println!(
            "I found {:} components, for this {:}-way partition\n",
            ncmps, nparts,
        ),
    );

    /* There are more components than partitions */
    if ncmps > nparts {
        let mut cwgt = vec![0; ncon].into_boxed_slice();
        // let mut bestcwgt = vec![0; ncon];
        let mut cpvec = vec![0; nparts].into_boxed_slice();
        let mut pcptr = vec![0; nparts + 1].into_boxed_slice();
        let mut pcind: Box<[idx_t]> = vec![0; ncmps].into_boxed_slice();
        let mut cwhere: Box<[idx_t]> = vec![-1; nvtxs].into_boxed_slice();
        let mut todo = vec![0; ncmps];
        // cand     = (rkv_t *)wspacemalloc(ctrl, nparts*sizeof(rkv_t));
        let mut cand = vec![KeyVal { key: 0.0, val: 0 }; nparts].into_boxed_slice();

        let mut modind = Vec::new();
        let mut vmarker = Vec::new();
        let mut pmarker = Vec::new();
        if ctrl.objtype == METIS_OBJTYPE_VOL {
            /* Vol-refinement specific working arrays */
            modind = vec![0; nvtxs];
            vmarker = vec![0; nvtxs];
            pmarker = vec![-1; nparts];
        }

        /* Get a CSR representation of the components-to-partitions mapping */
        for i in (0)..(ncmps) {
            pcptr[where_[cind[cptr[i] as usize] as usize] as usize] += 1;
        }
        util::make_csr(nparts, &mut pcptr);
        for i in (0)..(ncmps) {
            let p = &mut pcptr[where_[cind[cptr[i] as usize] as usize] as usize];
            pcind[*p as usize] = i as idx_t;
            *p += 1;
        }
        // nparts = pcptr.len() - 1
        util::shift_csr(nparts, &mut pcptr);

        /* Assign the heaviest component of each partition to its original partition */
        let mut ntodo = 0;
        for i in (0)..(nparts) {
            let mut bestcid: idx_t;
            let mut bestcwgt = 0;
            if pcptr[i + 1] - pcptr[i] == 1 {
                bestcid = pcind[pcptr[i] as usize];
            } else {
                bestcid = -1;
                for j in (pcptr[i])..(pcptr[i + 1]) {
                    let cid = pcind[j as usize];
                    // iset(ncon, 0, cwgt);
                    cwgt.fill(0);
                    for ii in cptr[cid as usize]..cptr[cid as usize + 1] {
                        blas::iaxpy(
                            ncon,
                            1,
                            &vwgt[(cind[ii as usize] as usize * ncon)..],
                            1,
                            &mut cwgt[..],
                            1,
                        );
                    }
                    if bestcid == -1 || bestcwgt < cwgt.iter().sum() {
                        bestcid = cid;
                        // icopy(ncon, cwgt, bestcwgt);
                        // bestcwgt.copy_from_slice(&cwgt);
                        bestcwgt = cwgt.iter().sum();
                    }
                }
                /* Keep track of those that need to be dealt with */
                for j in (pcptr[i])..(pcptr[i + 1]) {
                    if pcind[j as usize] != bestcid {
                        todo[ntodo] = pcind[j as usize];
                        ntodo += 1;
                    }
                }
            }

            for j in (cptr[bestcid as usize])..(cptr[bestcid as usize + 1]) {
                debug_assert_eq!(where_[cind[j as usize] as usize], i as idx_t);
                cwhere[cind[j as usize] as usize] = i as idx_t;
            }
        }

        while ntodo > 0 {
            let oldntodo = ntodo;
            let mut i = 0;
            // for i in (0)..(ntodo) {}
            while i < ntodo {
                let cid = todo[i];
                let me = where_[cind[cptr[cid as usize] as usize] as usize]; /* Get the domain of this component */

                /* Determine the weight of the block to be moved */
                // iset(ncon, 0, cwgt);
                cwgt.fill(0);
                for j in (cptr[cid as usize])..(cptr[cid as usize + 1]) {
                    blas::iaxpy(
                        ncon,
                        1,
                        &vwgt[(cind[j as usize] as usize * ncon)..],
                        1,
                        &mut cwgt[..],
                        1,
                    );
                }

                ifset!(
                    ctrl.dbglvl,
                    METIS_DBG_CONTIGINFO,
                    println!(
                        "Trying to move {:} [{:}] from {:}\n",
                        cid,
                        cwgt.iter().sum::<idx_t>(),
                        me,
                    ),
                );

                /* Determine the connectivity */
                // iset(nparts, 0, cpvec);
                cpvec.fill(0);
                for j in (cptr[cid as usize])..(cptr[cid as usize + 1]) {
                    let ii = cind[j as usize] as usize;
                    for jj in (xadj[ii])..(xadj[ii + 1]) {
                        let jj = jj as usize;
                        if cwhere[adjncy[jj] as usize] != -1 {
                            cpvec[cwhere[adjncy[jj] as usize] as usize] += {
                                // mkslice!(adjwgt, graph.nedges);
                                adjwgt.map_or(1, |adjwgt| adjwgt[jj])
                            }
                        }
                    }
                }

                /* Put the neighbors into a cand[] array for sorting */
                let mut ncand = 0;
                for j in (0)..(nparts) {
                    if cpvec[j] > 0 {
                        cand[ncand].key = cpvec[j] as real_t;
                        cand[ncand].val = j as idx_t;
                        ncand += 1;
                    }
                }
                if ncand == 0 {
                    i += 1;
                    continue;
                }

                // rkvsortd(ncand, cand);
                // cand.sort_unstable_by(|l, r| l.key.partial_cmp(&r.key).unwrap());
                cand[..ncand].sort_unstable_by(|l, r| l.key.total_cmp(&r.key).reverse());

                /* Limit the moves to only the top candidates, which are defined as
                those with connectivity at least 50% of the best.
                This applies only when ncon=1, as for multi-constraint, balancing
                will be hard. */
                if ncon == 1 {
                    ncand = cand[..ncand]
                        .iter()
                        .position(|x| x.key < 0.5 * cand[0].key)
                        .unwrap_or(ncand);
                }

                /* Now among those, select the one with the best balance */
                let mut target = cand[0].val as usize;
                for j in (1)..(ncand) {
                    // ctrl.pijbm is split, but shouldn't actually alias. I'll have to fix this for
                    // real later
                    if mcutil::BetterBalanceKWay(
                        ncon as idx_t,
                        cwgt.as_mut_ptr(),
                        ctrl.ubfactors,
                        1,
                        pwgts[(target * ncon)..].as_ptr(),
                        ctrl.pijbm.add(target * ncon),
                        1,
                        pwgts[(cand[j].val as usize * ncon)..].as_ptr(),
                        ctrl.pijbm.add(cand[j].val as usize * ncon),
                    ) != 0
                    {
                        target = cand[j].val as usize;
                    }
                }

                ifset!(
                    ctrl.dbglvl,
                    METIS_DBG_CONTIGINFO,
                    println!(
                        "\tMoving it to {:} [{:}] [{:}]\n",
                        target, cpvec[target], ncand,
                    ),
                );

                /* Note that as a result of a previous movement, a connected component may
                now will like to stay to its original partition */
                if target != me as usize {
                    match ctrl.objtype {
                        METIS_OBJTYPE_CUT => {
                            MoveGroupContigForCut(
                                ctrl,
                                graph,
                                target as idx_t,
                                cid,
                                cptr.as_ptr(),
                                cind.as_ptr(),
                            );
                        }

                        METIS_OBJTYPE_VOL => {
                            MoveGroupContigForVol(
                                ctrl,
                                graph,
                                target as idx_t,
                                cid,
                                cptr.as_ptr(),
                                cind.as_ptr(),
                                vmarker.as_mut_ptr(),
                                pmarker.as_mut_ptr(),
                                modind.as_mut_ptr(),
                            );
                        }

                        _ => panic!("Unknown objtype {}", ctrl.objtype),
                    }
                }

                /* Update the cwhere_ vector */
                for j in (cptr[cid as usize])..(cptr[cid as usize + 1]) {
                    cwhere[cind[j as usize] as usize] = target as idx_t;
                }

                // TODO: todo should be added and removed from instead of doing this funny
                // indexing, but I have to be careful of the exact pattern here. Notice that we
                // actually skip some entries in todo here, and leave them for another iteration of
                // the more outer loop
                ntodo -= 1;
                todo[i] = todo[ntodo];
                i += 1;
            }
            if oldntodo == ntodo {
                ifset!(
                    ctrl.dbglvl,
                    METIS_DBG_CONTIGINFO,
                    println!("Stopped at ntodo: {:}", ntodo),
                );
                break;
            }
        }

        for i in (0)..(nvtxs) {
            debug_assert_eq!(where_[i], cwhere[i],);
        }
    }

    debug_assert!(FindPartitionInducedComponents(graph, where_.as_ptr(), std::ptr::null_mut(), std::ptr::null_mut()) <= nparts as idx_t);

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
    ptr: *const idx_t,
    ind: *const idx_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, iii, j, jj, k, l, nvtxs, nbnd, from, me;
    // idx_t *xadj, *adjncy, *adjwgt, *where_, *bndptr, *bndind;
    // ckrinfo_t *myrinfo;
    // cnbr_t *mynbrs;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy adjwgt);
    get_graph_slices_mut!(graph => where_ bndptr bndind);

    mkslice!(ptr, nvtxs);
    mkslice!(ind, nvtxs);
    let gid = gid as usize;

    let mut nbnd = graph.nbnd as usize;

    for iii in (ptr[gid])..(ptr[gid + 1]) {
        let i = ind[iii as usize] as usize;
        let from = where_[i];

        let myrinfo = &mut *graph.ckrinfo.add(i as usize);
        if myrinfo.inbr == -1 {
            myrinfo.inbr = cnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
            myrinfo.nnbrs = 0;
        }
        let (myrinfo, mynbrs) = crinfos_mut(graph.ckrinfo, ctrl.cnbrpool, i);

        /* find the location of 'to' in myrinfo or create it if it is not there */
        let mut k = 0;
        while k < myrinfo.nnbrs as usize {
            if mynbrs[k].pid == to {
                break;
            }
            k += 1;
        }
        if k == myrinfo.nnbrs as usize {
            mynbrs[k].pid = to;
            mynbrs[k].ed = 0;
            myrinfo.nnbrs += 1;
        }
        let (myrinfo, mynbrs) = crinfos_mut(graph.ckrinfo, ctrl.cnbrpool, i);

        graph.mincut -= mynbrs[k].ed - myrinfo.id;

        /* Update ID/ED and BND related information for the moved vertex */
        {
            get_graph_slices!(graph => vwgt);
            get_graph_slices_mut!(ctrl, graph => pwgts);
            let ncon = graph.ncon as usize;
            blas::iaxpy(
                ncon as usize,
                1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(to as usize * ncon, ncon)],
                1,
            );
            blas::iaxpy(
                ncon as usize,
                -1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(from as usize * ncon, ncon)],
                1,
            );
        }
        kwayfm::UpdateMovedVertexInfoAndBND(
            i,
            from as usize,
            k,
            to as usize,
            myrinfo,
            mynbrs,
            where_,
            &mut nbnd,
            bndptr,
            bndind,
            BNDTYPE_REFINE
        );

        /* Update the degrees of adjacent vertices */
        for j in (xadj[i as usize])..(xadj[i as usize + 1]) {
            let j = j as usize;
            let ii = adjncy[j] as usize;
            let me = where_[ii];
            let (myrinfo, _mynbrs) = crinfos_mut(graph.ckrinfo, ctrl.cnbrpool, ii);

            kwayfm::UpdateAdjacentVertexInfoAndBND(
                ctrl,
                ii,
                xadj[ii + 1] - xadj[ii],
                me as usize,
                from as usize,
                to as usize,
                myrinfo,
                adjwgt[j],
                &mut nbnd,
                bndptr,
                bndind,
                BNDTYPE_REFINE
            );
        }

        debug_assert!(debug::CheckRInfo(ctrl, graph.ckrinfo.add(i as usize)) != 0);
    }

    graph.nbnd = nbnd as idx_t;
    debug_assert!(debug::CheckBnd2(graph) != 0);
    debug_assert_eq!(debug::ComputeCut(graph, where_.as_ptr()), graph.mincut);
    // assert!(debug::CheckBnd(graph) != 0);
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
    ptr: *const idx_t,
    ind: *const idx_t,
    vmarker: *mut idx_t,
    pmarker: *mut idx_t,
    modind: *mut idx_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj vsize adjncy);
    get_graph_slices_mut!(graph => where_);

    mkslice!(ptr, nvtxs);
    mkslice!(ind, nvtxs);
    let gid = gid as usize;

    for iii in (ptr[gid])..(ptr[gid + 1]) {
        let i = ind[iii as usize] as usize;
        let from = where_[i];

        {
            let myrinfo = &mut *graph.vkrinfo.add(i);
            if myrinfo.inbr == -1 {
                myrinfo.inbr = vnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
                myrinfo.nnbrs = 0;
            }
        }
        let (myrinfo, mynbrs) = vrinfos(graph.vkrinfo, ctrl.vnbrpool, i);

        let mut xgain = if myrinfo.nid == 0 && myrinfo.ned > 0 {
            vsize[i]
        } else {
            0
        };

        /* find the location of 'to' in myrinfo or create it if it is not there */
        if let Some(k) = mynbrs.iter().position(|rinfo| rinfo.pid == to) {
            graph.minvol -= xgain + mynbrs[k].gv;
            graph.mincut -= mynbrs[k].ned - myrinfo.nid;
        } else {
            if myrinfo.nid > 0 {
                xgain -= vsize[i];
            }

            /* determine the volume gain resulting from that move */
            for j in (xadj[i])..(xadj[i + 1]) {
                let ii = adjncy[j as usize] as usize;
                let other = where_[ii];
                let (_orinfo, onbrs) = vrinfos(graph.vkrinfo, ctrl.vnbrpool, ii);
                assert!(other != to);

                if from == other {
                    /* Same subdomain vertex: Decrease the gain if 'to' is a new neighbor. */
                    if !onbrs.iter().any(|rinfo| rinfo.pid == to) {
                        xgain -= vsize[ii];
                    }
                } else {
                    /* Remote vertex: increase if 'to' is a new subdomain */
                    if !onbrs.iter().any(|rinfo| rinfo.pid == to) {
                        xgain -= vsize[ii]; // GAVIN: Should this be +?
                    }

                    /* Remote vertex: decrease if i is the only connection to 'from' */
                    if onbrs
                        .iter()
                        .any(|rinfo| rinfo.pid == from && rinfo.ned == 1)
                    {
                        xgain += vsize[ii];
                    }
                }
            }
            graph.minvol -= xgain;
            graph.mincut -= -myrinfo.nid;
        }

        /* Update where_ and pwgts */
        where_[i] = to;
        // iaxpy( graph.ncon, 1, graph.vwgt + i * graph.ncon, 1, graph.pwgts + to * graph.ncon, 1,);
        // iaxpy( graph.ncon, -1, graph.vwgt + i * graph.ncon, 1, graph.pwgts + from * graph.ncon, 1,);
        {
            get_graph_slices!(graph => vwgt);
            get_graph_slices_mut!(ctrl, graph => pwgts);
            let ncon = graph.ncon as usize;
            blas::iaxpy(
                ncon as usize,
                1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(to as usize * ncon, ncon)],
                1,
            );
            blas::iaxpy(
                ncon as usize,
                -1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(from as usize * ncon, ncon)],
                1,
            );
        }

        /* Update the id/ed/gains/bnd of potentially affected nodes */
        mkslice_mut!(vmarker, nvtxs);
        mkslice_mut!(pmarker, ctrl.nparts);
        mkslice_mut!(modind, nvtxs);
        kwayfm::KWayVolUpdate(
            ctrl,
            graph,
            i,
            from as usize,
            to as usize,
            None,
            None,
            None,
            None,
            None,
            BNDTYPE_REFINE,
            vmarker,
            pmarker,
            modind,
        );

        /*CheckKWayVolPartitionParams(ctrl, graph);*/
    }

    debug_assert_eq!(
        debug::ComputeCutUnweighted(graph, where_.as_ptr()),
        graph.mincut
    );
    debug_assert_eq!(debug::ComputeVolume(graph, where_.as_ptr()), graph.minvol);
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use crate::tests::{ab_test_partition_test_graphs, ab_test_partition_test_graphs_filter};

    use super::*;
    use dyncall::ab_test_single_eq;
    use graph_gen::GraphBuilder;

    #[test]
    fn ab_FindPartitionInducedComponents() {
        ab_test_single_eq("FindPartitionInducedComponents:rs", || {
            fastrand::seed(123);
            // note that setting ncon to 2 causes an abort in C-only
            let mut g = GraphBuilder::new(Optype::Kmetis, 16, 1);
            g.set_seed(123);
            g.edge_list(
                std::iter::repeat_with(|| (fastrand::i32(0..=50), fastrand::i32(0..=50))).take(230),
            );
            g.set_contig(true);
            g.random_adjwgt();
            g.call().unwrap()
        });
    }

    #[test]
    fn ab_FindSepInducedComponents() {
        ab_test_partition_test_graphs_filter(
            "FindSepInducedComponents:rs",
            Optype::Ometis,
            3,
            1,
            |tg, mut g| {
                tg.is_contiguous().then(|| {
                    g.random_vwgt();
                    g.set_ccorder(true);
                    g
                })
            },
        );
    }

    #[test]
    fn ab_EliminateComponents_cut() {
        ab_test_single_eq("EliminateComponents:rs", || {
            fastrand::seed(123);
            let mut g = GraphBuilder::new(Optype::Kmetis, 16, 1);
            g.set_seed(123);
            g.edge_list(
                std::iter::repeat_with(|| (fastrand::i32(0..=50), fastrand::i32(0..=50))).take(230),
            );
            g.set_contig(true);
            g.set_objective(Objtype::Cut);
            g.random_adjwgt();
            g.call().unwrap()
        });
    }

    #[test]
    fn vol_c_contig_regression() {
        // This failed in earlier versions for two reasons:
        // 1) If refinement completely eliminated a partition, the check of
        //    FindPartitionInducedSubdomains == nparts would fail. This was fixed by changing the
        //    '==' to '<='. Perhaps it would be good to assert that refinement doesn't eliminate
        //    partitions in general, but it definitely doesn't make sense to only die when contig
        //    is set
        // 2) Volume refinement info and mindegree uses unweigted internal and external degrees,
        //    but certain places in the code assumed that it was weighted (most notably in
        //    ProjectKWayPartition). I changed this to always use unweighted, but this comes with
        //    an expected perf regression.
        let _k = dyncall::set_local_overrides("*");
        fastrand::seed(1230);
        let mut g = GraphBuilder::new(Optype::Kmetis, 16, 1);
        g.set_seed(1230);
        g.edge_list(
            std::iter::repeat_with(|| (fastrand::i32(0..=500), fastrand::i32(0..=500))).take(2300),
        );
        g.set_objective(Objtype::Vol);
        g.set_contig(true);
        g.random_vsize();
        g.random_vwgt();
        g.call().unwrap();
    }

    #[test]
    fn ab_EliminateComponents_vol() {
        ab_test_single_eq("EliminateComponents:rs", || {
            fastrand::seed(131);
            let mut g = GraphBuilder::new(Optype::Kmetis, 16, 1);
            g.set_seed(125);
            g.edge_list(
                std::iter::repeat_with(|| (fastrand::i32(0..=50), fastrand::i32(0..=50))).take(230),
            );
            assert!(g.is_contiguous());
            g.set_contig(true);
            g.set_objective(Objtype::Vol);
            g.call().unwrap()
        });
    }

    #[test]
    #[ignore = "need better graph"]
    fn ab_MoveGroupContigForCut() {
        ab_test_partition_test_graphs_filter(
            "MoveGroupContigForCut:rs",
            Optype::Kmetis,
            5,
            1,
            |tg, mut g| {
                tg.is_contiguous().then(move || {
                    g.set_objective(Objtype::Cut);
                    g.set_contig(true);
                    g.random_adjwgt();
                    g
                })
            },
        );
        // ab_test_single_eq("MoveGroupContigForCut:rs", || {
        //     fastrand::seed(123);
        //     let mut g = GraphBuilder::new(Optype::Kmetis, 16, 1);
        //     g.set_seed(123);
        //     g.edge_list(
        //         std::iter::repeat_with(|| (fastrand::i32(0..=50), fastrand::i32(0..=50))).take(230),
        //     );
        //     g.set_contig(true);
        //     g.random_adjwgt();
        //     g.call().unwrap()
        // });
    }

    // #[ignore = "broken on original"]
    #[test]
    fn ab_MoveGroupContigForVol() {
        ab_test_single_eq("MoveGroupContigForVol:rs", || {
            fastrand::seed(131);
            let mut g = GraphBuilder::new(Optype::Kmetis, 16, 1);
            g.set_seed(125);
            g.edge_list(
                std::iter::repeat_with(|| (fastrand::i32(0..=50), fastrand::i32(0..=50))).take(230),
            );
            assert!(g.is_contiguous());
            g.set_contig(true);
            // g.enable_dbg(DbgLvl::Refine);
            g.set_objective(Objtype::Vol);
            // g.write_graph({
            //     std::io::BufWriter::new(std::fs::OpenOptions::new().create(true).write(true).open("bug.graph").unwrap())
            // }).unwrap();
            // g.random_vsize();
            // g.random_vwgt();
            g.call().unwrap()
        });
    }
}
