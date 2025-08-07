/*
\file
\brief Driving routines for multilevel k-way refinement

\date   Started 7/28/1997
\author George
\author  Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: kwayrefine.c 20398 2016-11-22 17:17:12Z karypis $
*/

use crate::*;

use core::ptr;
use std::slice;

/// This function is the entry point of cut-based refinement
#[metis_func]
pub extern "C" fn RefineKWay(ctrl: *mut ctrl_t, orggraph: *mut graph_t, graph: *mut graph_t) {
    let mut nlevels;
    let contig = (*ctrl).contig;
    let mut ptr: *mut graph_t;
    let mut graph = graph;

    /* Determine how many levels are there */
    ptr = graph;
    nlevels = 0;
    while ptr != orggraph {
        nlevels += 1;
        ptr = (*ptr).finer;
    }

    /* Compute the parameters of the coarsest graph */
    ComputeKWayPartitionParams(ctrl, graph);

    /* Try to minimize the sub-domain connectivity */
    if (*ctrl).minconn != 0 {
        minconn::EliminateSubDomainEdges(ctrl, graph);
    }

    /* Deal with contiguity constraints at the beginning */
    if contig != 0
        && contig::FindPartitionInducedComponents(
            graph,
            (*graph).where_,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        ) > (*ctrl).nparts
    {
        contig::EliminateComponents(ctrl, graph);

        ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
        kwayfm::Greedy_KWayOptimize(ctrl, graph, 5, 0.0, OMODE_BALANCE);

        ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
        kwayfm::Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 0.0, OMODE_REFINE);

        (*ctrl).contig = 0;
    }

    /* Refine each successively finer graph */
    for i in 0.. {
        if (*ctrl).minconn != 0 && i == nlevels / 2 {
            minconn::EliminateSubDomainEdges(ctrl, graph);
        }

        // IFSET((*ctrl).dbglvl, METIS_DBG_TIME, gk_startcputimer((*ctrl).RefTmr));

        if 2 * i >= nlevels && IsBalanced(ctrl, graph, 0.02) == 0 {
            ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
            kwayfm::Greedy_KWayOptimize(ctrl, graph, 1, 0.0, OMODE_BALANCE);
            ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
        }

        kwayfm::Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 5.0, OMODE_REFINE);

        // IFSET((*ctrl).dbglvl, METIS_DBG_TIME, gk_stopcputimer((*ctrl).RefTmr));

        /* Deal with contiguity constraints in the middle */
        if contig != 0
            && i == nlevels / 2
            && contig::FindPartitionInducedComponents(
                graph,
                (*graph).where_,
                ptr::null_mut(),
                ptr::null_mut(),
            ) > (*ctrl).nparts
        {
            contig::EliminateComponents(ctrl, graph);

            if IsBalanced(ctrl, graph, 0.02) == 0 {
                (*ctrl).contig = 1;
                ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
                kwayfm::Greedy_KWayOptimize(ctrl, graph, 5, 0.0, OMODE_BALANCE);

                ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
                kwayfm::Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 0.0, OMODE_REFINE);
                (*ctrl).contig = 0;
            }
        }

        if graph == orggraph {
            break;
        }

        graph = (*graph).finer;

        graph::graph_ReadFromDisk(ctrl, graph);

        // IFSET(
        //     (*ctrl).dbglvl,
        //     METIS_DBG_TIME,
        //     gk_startcputimer((*ctrl).ProjectTmr),
        // );
        assert!(!(*graph).vwgt.is_null());

        ProjectKWayPartition(ctrl, graph);
        // IFSET(
        //     (*ctrl).dbglvl,
        //     METIS_DBG_TIME,
        //     gk_stopcputimer((*ctrl).ProjectTmr),
        // );
    }

    /* Deal with contiguity requirement at the end */
    (*ctrl).contig = contig;
    if contig != 0
        && contig::FindPartitionInducedComponents(
            graph,
            (*graph).where_,
            ptr::null_mut(),
            ptr::null_mut(),
        ) > (*ctrl).nparts
    {
        contig::EliminateComponents(ctrl, graph);
    }

    if IsBalanced(ctrl, graph, 0.0) == 0 {
        ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
        kwayfm::Greedy_KWayOptimize(ctrl, graph, 10, 0.0, OMODE_BALANCE);

        ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
        kwayfm::Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 0.0, OMODE_REFINE);
    }

    if (*ctrl).contig != 0 {
        debug_assert!(
            contig::FindPartitionInducedComponents(
                graph,
                (*graph).where_,
                ptr::null_mut(),
                ptr::null_mut()
            ) <= (*ctrl).nparts,
        );
    }

    // IFSET(
    //     (*ctrl).dbglvl,
    //     METIS_DBG_TIME,
    //     gk_stopcputimer((*ctrl).UncoarsenTmr),
    // );
}

/// This function allocates memory for the k-way cut-based refinement
#[metis_func]
pub fn AllocateKWayPartitionMemory(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    graph.pwgts = imalloc(
        (ctrl.nparts * graph.ncon) as usize,
        c"AllocateKWayPartitionMemory: pwgts".as_ptr(),
    ) as _;
    graph.where_ = imalloc(
        graph.nvtxs as usize,
        c"AllocateKWayPartitionMemory: where".as_ptr(),
    ) as _;
    graph.bndptr = imalloc(
        graph.nvtxs as usize,
        c"AllocateKWayPartitionMemory: bndptr".as_ptr(),
    ) as _;
    graph.bndind = imalloc(
        graph.nvtxs as usize,
        c"AllocateKWayPartitionMemory: bndind".as_ptr(),
    ) as _;

    match ctrl.objtype {
        METIS_OBJTYPE_CUT => {
            graph.ckrinfo = gk_malloc(
                graph.nvtxs as usize * std::mem::size_of::<ckrinfo_t>(),
                c"AllocateKWayPartitionMemory: ckrinfo".as_ptr(),
            ) as *mut _;
        }

        METIS_OBJTYPE_VOL => {
            graph.vkrinfo = gk_malloc(
                graph.nvtxs as usize * std::mem::size_of::<vkrinfo_t>(),
                c"AllocateKWayVolPartitionMemory: vkrinfo".as_ptr(),
            ) as *mut _;

            /* This is to let the cut-based -minconn and -contig large-scale graph
            changes to go through */
            // wtf this is terrible??? they aren't even the same size??
            // FIXME: figure out wtf is happening here
            graph.ckrinfo = graph.vkrinfo as *mut _;
        }

        _ => {
            panic!("Unknown objtype of {}", ctrl.objtype);
        }
    };
}

/// This function computes the initial id/ed  for cut-based partitioning
///
/// I suspect bug in this function (volume)
#[metis_func]
pub extern "C" fn ComputeKWayPartitionParams(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    let mut nbnd: idx_t;
    let mut mincut: idx_t;

    let nparts: idx_t = ctrl.nparts;

    let nvtxs: idx_t = graph.nvtxs;
    let ncon: idx_t = graph.ncon;

    get_graph_slices!(ctrl, graph => xadj adjncy vwgt where_ adjwgt);
    get_graph_slices_mut!(ctrl, graph => pwgts bndptr bndind ckrinfo);
    pwgts.fill(0);
    bndptr.fill(-1);

    nbnd = 0;
    mincut = 0;

    /* Compute pwgts (partition weights) */
    if ncon == 1 {
        for i in 0..(nvtxs as usize) {
            assert!(where_[i] >= 0 && where_[i] < nparts);
            pwgts[where_[i] as usize] += vwgt[i];
        }
    } else {
        for i in 0..(nvtxs as usize) {
            let me: idx_t = where_[i];
            for j in 0..ncon {
                pwgts[(me * ncon + j) as usize] += vwgt[(i as idx_t * ncon + j) as usize];
            }
        }
    }

    /* Compute the required info for refinement */
    match ctrl.objtype {
        METIS_OBJTYPE_CUT => {
            // cnbrpool operations may realloc, so we need to be careful
            ckrinfo.fill_with(std::default::Default::default);

            // memset(rsgraph.ckrinfo, 0, sizeof(ckrinfo_t) * nvtxs);
            // ckrinfo_t * myrinfo;
            // let myrinfo: &mut ckrinfo_t;

            // cnbr_t * mynbrs;
            // let mynbrs: &mut [cnbr_t];

            cnbrpoolReset(ctrl);

            for i in 0..(nvtxs as usize) {
                let me = where_[i];
                // myrinfo = rsgraph.ckrinfo + i;
                let myrinfo = &mut ckrinfo[i];

                for j in (xadj[i] as usize)..(xadj[i + 1] as usize) {
                    if me == where_[adjncy[j] as usize] {
                        myrinfo.id += adjwgt[j];
                    } else {
                        myrinfo.ed += adjwgt[j];
                    }
                }

                /* Time to compute the particular external degrees */
                if myrinfo.ed > 0 {
                    mincut += myrinfo.ed;

                    myrinfo.inbr = cnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);

                    for j in xadj[i]..xadj[i + 1] {
                        let j = j as usize;
                        let other = where_[adjncy[j] as usize];
                        let mynbrs = slice::from_raw_parts_mut(
                            ctrl.cnbrpool.add(myrinfo.inbr as usize),
                            myrinfo.nnbrs as usize + 1,
                        );
                        if me != other {
                            let mut k = 0;
                            while k < myrinfo.nnbrs as usize {
                                if mynbrs[k].pid == other {
                                    mynbrs[k].ed += adjwgt[j];
                                    break;
                                }
                                k += 1;
                            }
                            if k == myrinfo.nnbrs as usize {
                                mynbrs[k].pid = other;
                                mynbrs[k].ed = adjwgt[j];
                                myrinfo.nnbrs += 1;
                            }
                        }
                    }

                    debug_assert!(myrinfo.nnbrs <= xadj[i + 1] - xadj[i]);

                    /* Only ed-id>=0 nodes are considered to be in the boundary */
                    if myrinfo.ed - myrinfo.id >= 0 {
                        BNDInsert!(nbnd, bndind, bndptr, i);
                    }
                } else {
                    myrinfo.inbr = -1;
                }
            }

            graph.mincut = mincut / 2;
            graph.nbnd = nbnd;

            debug_assert!(debug::CheckBnd2(graph) != 0);
        }

        METIS_OBJTYPE_VOL => {
            {
                // memset((*graph).vkrinfo, 0, sizeof(vkrinfo_t) * nvtxs);
                // cnbrpool operations may realloc, so we can't leave this be
                get_graph_slices_mut!(graph => vkrinfo);
                vkrinfo.fill_with(std::default::Default::default);
            }
            // memset((*graph).ckrinfo, 0, sizeof(ckrinfo_t) * nvtxs);
            // vkrinfo_t * myrinfo;
            // vnbr_t *mynbrs;

            vnbrpoolReset(ctrl);

            /* Compute now the id/ed degrees */
            for i in 0..(nvtxs as usize) {
                let me = where_[i];
                let myrinfo;
                {
                    get_graph_slices_mut!(graph => vkrinfo);
                    myrinfo = &mut vkrinfo[i];
                }

                for j in xadj[i]..xadj[i + 1] {
                    if me == where_[adjncy[j as usize] as usize] {
                        myrinfo.nid += 1;
                    } else {
                        myrinfo.ned += 1;
                    }
                }
                debug_assert_eq!(myrinfo.nnbrs, 0);
                /* Time to compute the particular external degrees */
                if myrinfo.ned > 0 {
                    mincut += myrinfo.ned;

                    myrinfo.inbr = vnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
                    let mut mynbrs = slice::from_raw_parts_mut(
                        ctrl.vnbrpool.add(myrinfo.inbr as usize),
                        myrinfo.nnbrs as usize + 1,
                    );

                    for j in xadj[i]..xadj[i + 1] {
                        let j = j as usize;
                        let other = where_[adjncy[j] as usize];
                        if me != other {
                            let mut k = 0;
                            while k < (myrinfo.nnbrs as usize) {
                                if mynbrs[k].pid == other {
                                    mynbrs[k].ned += 1;
                                    break;
                                }
                                k += 1;
                            }
                            // add as neighbor if not encountered
                            if k == myrinfo.nnbrs as usize {
                                mynbrs[k].gv = 0;
                                mynbrs[k].pid = other;
                                mynbrs[k].ned = 1;
                                myrinfo.nnbrs += 1;
                                mynbrs = slice::from_raw_parts_mut(
                                    ctrl.vnbrpool.add(myrinfo.inbr as usize),
                                    myrinfo.nnbrs as usize + 1,
                                );
                            }
                        }
                    }
                    assert!(myrinfo.nnbrs <= xadj[i + 1] - xadj[i]);
                } else {
                    myrinfo.inbr = -1;
                }
            }
            graph.mincut = mincut / 2;

            ComputeKWayVolGains(ctrl, graph);
            debug_assert_eq!(graph.minvol, debug::ComputeVolume(graph, graph.where_));
        }
        _ => panic!("Unknown objtype of {}", ctrl.objtype),
    }
}

/// This function projects a partition, and at the same time computes the parameters for refinement.
#[metis_func]
pub extern "C" fn ProjectKWayPartition(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    let mut nbnd: idx_t;

    let dropedges = ctrl.dropedges;

    let nparts: idx_t = ctrl.nparts;

    let cgraph: *mut graph_t = graph.coarser;
    // TODO: I think this maybe should be cgraph.nvtxs
    mkslice_mut!(cwhere: cgraph->where_, (*cgraph).nvtxs);

    if ctrl.objtype == METIS_OBJTYPE_CUT {
        assert!(debug::CheckBnd2(cgraph) != 0);
    } else {
        assert!((*cgraph).minvol == debug::ComputeVolume(cgraph, (*cgraph).where_));
    }

    /* free the coarse graph's structure (reduce maxmem) */
    // disabling this since it violates the normally strict wspace rules - so it's hard to deal
    // with over ffi
    graph::FreeSData(cgraph);

    let nvtxs: idx_t = graph.nvtxs;

    AllocateKWayPartitionMemory(ctrl, graph);

    get_graph_slices!(ctrl, *graph => xadj adjncy adjwgt);
    get_graph_slices_mut!(ctrl, *graph => cmap bndind bndptr where_);
    bndptr.fill(-1);

    let mut htable: Vec<idx_t> = vec![-1; nparts as usize];

    /* Compute the required info for refinement */
    match ctrl.objtype {
        METIS_OBJTYPE_CUT => {
            /* go through and project partition and compute id/ed for the nodes */
            for i in 0..(nvtxs as usize) {
                let k = cmap[i] as usize;
                where_[i] = cwhere[k];
                cmap[i] = if dropedges != 0 {
                    1
                } else {
                    (*(*cgraph).ckrinfo.add(k)).ed
                }; /* For optimization */
            }

            {
                // memset((*graph).ckrinfo, 0, sizeof(ckrinfo_t) * nvtxs);
                mkslice_mut!(graph->ckrinfo, nvtxs);
                ckrinfo.fill_with(std::default::Default::default);
            }

            cnbrpoolReset(ctrl);

            nbnd = 0;
            for i in 0..(nvtxs as usize) {
                let istart = xadj[i] as usize;
                let iend = xadj[i + 1] as usize;

                let myrinfo = graph.ckrinfo.add(i).as_mut().unwrap();

                if cmap[i] == 0 {
                    /* Interior node. Note that cmap[i] = crinfo[cmap[i]].ed */
                    let mut tid = 0;
                    for j in istart..iend {
                        tid += adjwgt[j];
                    }

                    myrinfo.id = tid;
                    myrinfo.inbr = -1;
                } else {
                    /* Potentially an interface node */
                    myrinfo.inbr = cnbrpoolGetNext(ctrl, iend as idx_t - istart as idx_t);
                    let mut mynbrs = slice::from_raw_parts_mut(
                        ctrl.cnbrpool.add(myrinfo.inbr as usize),
                        myrinfo.nnbrs as usize + 1,
                    );

                    let me = where_[i] as usize;
                    let mut tid = 0;
                    let mut ted = 0;
                    let mut k;
                    for j in istart..iend {
                        let other = where_[adjncy[j] as usize] as usize;
                        if me == other {
                            tid += adjwgt[j];
                        } else {
                            ted += adjwgt[j];
                            k = htable[other];
                            if k == -1 {
                                htable[other] = myrinfo.nnbrs;
                                mynbrs[myrinfo.nnbrs as usize].pid = other as idx_t;
                                mynbrs[myrinfo.nnbrs as usize].ed = adjwgt[j];
                                myrinfo.nnbrs += 1;
                                mynbrs = slice::from_raw_parts_mut(
                                    mynbrs.as_mut_ptr(),
                                    myrinfo.nnbrs as usize + 1,
                                );
                            } else {
                                mynbrs[k as usize].ed += adjwgt[j];
                            }
                        }
                    }
                    myrinfo.id = tid;
                    myrinfo.ed = ted;

                    /* Remove space for edegrees if it was interior */
                    if ted == 0 {
                        ctrl.nbrpoolcpos -= (nparts as usize).min(iend - istart);
                        myrinfo.inbr = -1;
                    } else {
                        if ted - tid >= 0 {
                            BNDInsert!(nbnd, bndind, bndptr, i);
                        }

                        for j in 0..myrinfo.nnbrs {
                            htable[mynbrs[j as usize].pid as usize] = -1;
                        }
                    }
                }
            }

            graph.nbnd = nbnd;
            debug_assert!(debug::CheckBnd2(graph) != 0);
            graph.mincut = if dropedges != 0 {
                debug::ComputeCut(graph, where_.as_mut_ptr())
            } else {
                (*cgraph).mincut
            };
        }

        METIS_OBJTYPE_VOL => {
            /* go through and project partition and compute id/ed for the nodes */
            let mut k;
            for i in 0..(nvtxs as usize) {
                k = cmap[i];
                where_[i] = cwhere[k as usize];
                cmap[i] = if dropedges != 0 {
                    1
                } else {
                    (*(*cgraph).vkrinfo.add(k as usize)).ned
                }; /* For optimization */
            }

            {
                get_graph_slices_mut!(graph => vkrinfo);
                vkrinfo.fill_with(std::default::Default::default); // zeroed
            }
            vnbrpoolReset(ctrl);

            for i in 0..(nvtxs as usize) {
                let istart = xadj[i] as usize;
                let iend = xadj[i + 1] as usize;
                let myrinfo = graph.vkrinfo.add(i).as_mut().unwrap();

                if cmap[i] == 0 {
                    /* Note that cmap[i] = crinfo[cmap[i]].ed */
                    myrinfo.nid = (iend - istart) as idx_t;
                    myrinfo.inbr = -1;
                } else {
                    /* Potentially an interface node */
                    myrinfo.inbr = vnbrpoolGetNext(ctrl, (iend - istart) as idx_t);
                    // mynbrs = (*ctrl).vnbrpool + myrinfo.inbr;
                    debug_assert!(myrinfo.inbr >= 0);

                    let me = where_[i];
                    let mut tid = 0;
                    let mut ted = 0;
                    for j in istart..iend {
                        let other = where_[adjncy[j] as usize];
                        let mynbrs = slice::from_raw_parts_mut(
                            ctrl.vnbrpool.add(myrinfo.inbr as usize),
                            // (nvtxs - myrinfo.inbr) as usize,
                            myrinfo.nnbrs as usize + 1,
                        );
                        if me == other {
                            tid += 1;
                        } else {
                            ted += 1;
                            k = htable[other as usize];
                            if k == -1 {
                                htable[other as usize] = myrinfo.nnbrs;
                                mynbrs[myrinfo.nnbrs as usize].gv = 0;
                                mynbrs[myrinfo.nnbrs as usize].pid = other;
                                mynbrs[myrinfo.nnbrs as usize].ned = 1;
                                myrinfo.nnbrs += 1;
                            } else {
                                mynbrs[k as usize].ned += 1;
                            }
                        }
                    }
                    myrinfo.nid = tid;
                    myrinfo.ned = ted;

                    /* Remove space for edegrees if it was interior */
                    if ted == 0 {
                        ctrl.nbrpoolcpos -= (nparts as usize).min(iend - istart);
                        myrinfo.inbr = -1;
                    } else {
                        let mynbrs = slice::from_raw_parts_mut(
                            ctrl.vnbrpool.add(myrinfo.inbr as usize),
                            // (nvtxs - myrinfo.inbr) as usize,
                            myrinfo.nnbrs as usize + 1,
                        );
                        for j in 0..myrinfo.nnbrs {
                            htable[mynbrs[j as usize].pid as usize] = -1;
                        }
                    }
                }
            }

            ComputeKWayVolGains(ctrl, graph);

            debug_assert_eq!(graph.minvol, debug::ComputeVolume(graph, graph.where_));
            /* Can't use dropedges since volume rinfo does not keep track of     */
            /* weight -- cgraph's mincut will generally not be valid for graph.  */
            /* Doing it this way is a performance pessimization so perhaps it    */
            /* would be better to match cut refinement's behavior in using rinfo */
            // TODO: can I use weigted cut on cgraph?
            // PERF: this was pessimized from the original (but this is more correct I think)
            // graph.mincut = if dropedges != 0 { 
            //     debug::ComputeCutUnweighted(graph, where_.as_mut_ptr()) 
            // } else { (*cgraph).mincut };
            graph.mincut = debug::ComputeCutUnweighted(graph, where_.as_mut_ptr());
        }

        _ => panic!("Unknown objtype of {}", ctrl.objtype),
    }

    {
        mkslice_mut!(dst: graph->pwgts, nparts * graph.ncon);
        mkslice_mut!(src: cgraph->pwgts, nparts * graph.ncon);
        dst.copy_from_slice(src);
        // remember: icopy(n, src, dest) -> reverse of memcpy
        // icopy(nparts * (*graph).ncon, (*cgraph).pwgts, (*graph).pwgts);
    }

    graph::FreeGraph((&mut graph.coarser) as *mut *mut graph_t);
}

/// This function computes the boundary definition for balancing.
#[metis_func]
pub extern "C" fn ComputeKWayBoundary(ctrl: *mut ctrl_t, graph: *mut graph_t, bndtype: idx_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let mut nbnd;

    mkslice_mut!(graph->bndptr, graph.nvtxs);
    bndptr.fill(-1);

    nbnd = 0;
    match ctrl.objtype {
        METIS_OBJTYPE_CUT => {
            /* Compute the boundary */
            get_graph_slices_mut!(ctrl, graph => bndind bndptr);
            get_graph_slices!(ctrl, graph => ckrinfo);
            let nvtxs = graph.nvtxs as usize;
            if bndtype == BNDTYPE_REFINE {
                for i in 0..nvtxs {
                    if ckrinfo[i].ed > 0 && ckrinfo[i].ed - ckrinfo[i].id >= 0 {
                        BNDInsert!(nbnd, bndind, bndptr, i);
                    }
                }
            } else {
                /* BNDTYPE_BALANCE */
                for i in 0..nvtxs {
                    if ckrinfo[i].ed > 0 {
                        BNDInsert!(nbnd, bndind, bndptr, i);
                    }
                }
            }
        }

        METIS_OBJTYPE_VOL => {
            /* Compute the boundary */
            get_graph_slices_mut!(ctrl, graph => bndind bndptr);
            get_graph_slices!(ctrl, graph => vkrinfo);
            let nvtxs = graph.nvtxs as usize;
            if bndtype == BNDTYPE_REFINE {
                for i in 0..nvtxs {
                    if vkrinfo[i].gv >= 0 {
                        BNDInsert!(nbnd, bndind, bndptr, i);
                    }
                }
            } else {
                /* BNDTYPE_BALANCE */
                for i in 0..nvtxs {
                    if vkrinfo[i].ned > 0 {
                        BNDInsert!(nbnd, bndind, bndptr, i);
                    }
                }
            }
        }
        _ => panic!("Unknown objtype of {}\n", ctrl.objtype),
    }

    graph.nbnd = nbnd;
}

/// This function computes the initial gains in the communication volume
#[metis_func]
pub extern "C" fn ComputeKWayVolGains(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let ctrl: &mut ctrl_t = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    let nparts: idx_t = ctrl.nparts;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj vsize adjncy);
    get_graph_slices_mut!(ctrl, graph => where_ bndind bndptr vkrinfo);

    // bndptr = iset(nvtxs, -1, graph.bndptr);
    bndptr.fill(-1);

    // ophtable = iset(nparts, -1, iwspacemalloc(ctrl, nparts));
    // other partition h_____ table? Referred to as "marker vector"
    let mut ophtable: Vec<idx_t> = vec![-1; nparts as usize];

    /* Compute the volume gains */
    graph.minvol = 0;
    graph.nbnd = 0;
    for i in 0..nvtxs {
        // basically myrinfo
        vkrinfo[i].gv = idx_t::MIN;

        // basically myrinfo
        if vkrinfo[i].nnbrs > 0 {
            let myrinfo = &vkrinfo[i];
            let me = where_[i];
            debug_assert!(myrinfo.inbr != -1);
            // mynbrs = (*ctrl).vnbrpool + myrinfo.inbr;
            let mynbrs = slice::from_raw_parts_mut(
                ctrl.vnbrpool.add(myrinfo.inbr as usize),
                // changing this to debug segfaults - majority of tests fail in this file
                // ctrl.nbrpoolsize - myrinfo.inbr as usize,
                myrinfo.nnbrs as usize + 1,
            );

            graph.minvol += myrinfo.nnbrs * vsize[i];

            for j in xadj[i]..xadj[i + 1] {
                let j = j as usize;
                let ii = adjncy[j] as usize;
                let other = where_[ii];
                let orinfo = vkrinfo[ii].clone();

                for k in 0..orinfo.nnbrs {
                    // FIX: move construction into here to prevent invalid pointer arithmetic
                    debug_assert!(orinfo.inbr != -1);
                    let onbrs = slice::from_raw_parts_mut(
                        ctrl.vnbrpool.add(orinfo.inbr as usize),
                        orinfo.nnbrs as usize,
                    );
                    ophtable[onbrs[k as usize].pid as usize] = k;
                }
                ophtable[other as usize] = 1; /* this is to simplify coding */

                if me == other {
                    /* Find which domains 'i' is connected to but 'ii' is not
                    and update their gain */
                    for k in 0..myrinfo.nnbrs {
                        if ophtable[mynbrs[k as usize].pid as usize] == -1 {
                            mynbrs[k as usize].gv -= vsize[ii];
                        }
                    }
                } else {
                    debug_assert!(orinfo.inbr != -1);
                    let onbrs = slice::from_raw_parts_mut(
                        ctrl.vnbrpool.add(orinfo.inbr as usize),
                        orinfo.nnbrs as usize,
                    );
                    let me = me as usize;
                    assert!(
                        ophtable[me] != -1,
                        concat!(
                            "vtx {i} is adj to vtx {ii} but {ii}'s neighboring subdomains' ",
                            "partitions do not include {i}'s (part = {me}).\n",
                            "vtx {ii}'s neighbors' partitions are {nbr:?}",
                        ),
                        i = i,
                        ii = ii,
                        me = me,
                        nbr = (onbrs.iter().map(|d| d.pid).collect::<Vec<_>>())
                    );

                    if (onbrs[ophtable[me] as usize]).ned == 1 {
                        /* I'm the only connection of 'ii' in 'me' */
                        /* Increase the gains for all the common domains between 'i' and 'ii' */
                        for k in 0..myrinfo.nnbrs {
                            let k = k as usize;
                            if ophtable[mynbrs[k].pid as usize] != -1 {
                                mynbrs[k].gv += vsize[ii];
                            }
                        }
                    } else {
                        /* Find which domains 'i' is connected to and 'ii' is not
                        and update their gain */
                        for k in 0..myrinfo.nnbrs {
                            let k = k as usize;
                            if ophtable[mynbrs[k].pid as usize] == -1 {
                                mynbrs[k].gv -= vsize[ii];
                            }
                        }
                    }
                }

                /* Reset the marker vector */
                for k in 0..orinfo.nnbrs {
                    let onbrs = slice::from_raw_parts_mut(
                        ctrl.vnbrpool.add(orinfo.inbr as usize),
                        orinfo.nnbrs as usize,
                    );
                    ophtable[onbrs[k as usize].pid as usize] = -1;
                }
                ophtable[other as usize] = -1;
            }
            let myrinfo = &mut vkrinfo[i];

            /* Compute the max vgain */
            for k in 0..myrinfo.nnbrs {
                let k = k as usize;
                if mynbrs[k].gv > myrinfo.gv {
                    myrinfo.gv = mynbrs[k].gv;
                }
            }

            /* Add the extra gain due to id == 0 */
            if myrinfo.ned > 0 && myrinfo.nid == 0 {
                myrinfo.gv += vsize[i];
            }
        }

        // myrinfo
        if vkrinfo[i].gv >= 0 {
            BNDInsert!(graph.nbnd, bndind, bndptr, i);
        }
    }
}

/// This function checks if the partition weights are within the balance constraints
#[metis_func]
pub extern "C" fn IsBalanced(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ffactor: real_t,
) -> std::ffi::c_int {
    (mcutil::ComputeLoadImbalanceDiff(graph, (*ctrl).nparts, (*ctrl).pijbm, (*ctrl).ubfactors)
        <= ffactor) as _
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use std::convert::identity;

    use super::*;
    use dyncall::ab_test_single_eq;
    use crate::tests::{ab_test_partition_test_graphs, ab_test_partition_test_graphs_filter};
    use graph_gen::GraphBuilder;

    #[test]
    fn kwayrefine_contig_mc_abort_regression() {
        // See: contig::tests::vol_c_contig_regression
        fastrand::seed(123);
        // note that setting ncon to 2 causes an abort even in C-only
        let mut g = GraphBuilder::new(Optype::Kmetis, 16, 2);
        g.set_seed(123);
        g.edge_list(
            std::iter::repeat_with(|| (fastrand::i32(0..=50), fastrand::i32(0..=50))).take(2300),
        );
        g.set_contig(true);
        g.random_adjwgt();
        g.call().unwrap();
    }

    #[test]
    fn ab_RefineKWay() {
        ab_test_partition_test_graphs("RefineKWay:rs", Optype::Kmetis, 20, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_RefineKWay_contig() {
        ab_test_partition_test_graphs_filter("RefineKWay:rs", Optype::Kmetis, 20, 1, |tg, mut g| {
            tg.is_contiguous().then(|| {
                g.set_contig(true);
                g.random_vwgt();
                g
            })
        });
    }

    #[test]
    fn ab_AllocateKWayPartitionMemory() {
        // this is a test since there's weird things happening in this function right now
        let func = "AllocateKWayPartitionMemory:rs";
        ab_test_partition_test_graphs(func, Optype::Kmetis, 2, 1, identity);
        ab_test_partition_test_graphs(func, Optype::Kmetis, 2, 1, |mut g| {
            g.set_objective(Objtype::Vol);
            g
        });
    }

    #[test]
    fn ab_ComputeKWayPartitionParams() {
        ab_test_partition_test_graphs("ComputeKWayPartitionParams:rs", Optype::Kmetis, 20, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_ProjectKWayPartition_cut() {
        ab_test_partition_test_graphs("ProjectKWayPartition:rs", Optype::Kmetis, 20, 1, |mut g| {
            g.set_objective(Objtype::Cut);
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_ProjectKWayPartition_vol() {
        ab_test_partition_test_graphs("ProjectKWayPartition:rs", Optype::Kmetis, 20, 1, |mut g| {
            g.set_objective(Objtype::Vol);
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_ComputeKWayBoundary_vol() {
        ab_test_partition_test_graphs("ComputeKWayBoundary:rs", Optype::Kmetis, 20, 1, |mut g| {
            // g.set_contig(true);
            g.set_objective(Objtype::Vol);
            g
        });
    }

    #[test]
    fn ab_ComputeKWayBoundary_cut() {
        ab_test_partition_test_graphs("ComputeKWayBoundary:rs", Optype::Kmetis, 40, 1, |mut g| {
            g.set_ufactor(5);
            // g.set_contig(true);
            g.set_objective(Objtype::Cut);
            g
        });
    }

    #[test]
    fn ab_ComputeKWayVolGains() {
        ab_test_partition_test_graphs("ComputeKWayVolGains:rs", Optype::Kmetis, 5, 1, |mut g| {
            g.set_objective(Objtype::Vol);
            g.random_vwgt();
            g
        });
    }
}
