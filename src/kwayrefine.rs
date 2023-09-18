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
pub extern "C" fn RefineKWay(ctrl: *mut ctrl_t, orggraph: *mut graph_t, graph: *mut graph_t) -> () {
    eprintln!("Called RefineKWay");
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
        EliminateSubDomainEdges(ctrl, graph);
    }

    /* Deal with contiguity constraints at the beginning */
    if contig != 0
        && FindPartitionInducedComponents(
            graph,
            (*graph).where_,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        ) > (*ctrl).nparts
    {
        EliminateComponents(ctrl, graph);

        ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
        Greedy_KWayOptimize(ctrl, graph, 5, 0.0, OMODE_BALANCE);

        ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
        Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 0.0, OMODE_REFINE);

        (*ctrl).contig = 0;
    }

    /* Refine each successively finer graph */
    for i in 0.. {
        if (*ctrl).minconn != 0 && i == nlevels / 2 {
            EliminateSubDomainEdges(ctrl, graph);
        }

        // IFSET((*ctrl).dbglvl, METIS_DBG_TIME, gk_startcputimer((*ctrl).RefTmr));

        if 2 * i >= nlevels && IsBalanced(ctrl, graph, 0.02) == 0 {
            ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
            Greedy_KWayOptimize(ctrl, graph, 1, 0.0, OMODE_BALANCE);
            ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
        }

        Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 5.0, OMODE_REFINE);

        // IFSET((*ctrl).dbglvl, METIS_DBG_TIME, gk_stopcputimer((*ctrl).RefTmr));

        /* Deal with contiguity constraints in the middle */
        if contig != 0 && i == nlevels / 2 {
            if FindPartitionInducedComponents(
                graph,
                (*graph).where_,
                ptr::null_mut(),
                ptr::null_mut(),
            ) > (*ctrl).nparts
            {
                EliminateComponents(ctrl, graph);

                if IsBalanced(ctrl, graph, 0.02) == 0 {
                    (*ctrl).contig = 1;
                    ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
                    Greedy_KWayOptimize(ctrl, graph, 5, 0.0, OMODE_BALANCE);

                    ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
                    Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 0.0, OMODE_REFINE);
                    (*ctrl).contig = 0;
                }
            }
        }

        if graph == orggraph {
            break;
        }

        graph = (*graph).finer;

        graph_ReadFromDisk(ctrl, graph);

        // IFSET(
        //     (*ctrl).dbglvl,
        //     METIS_DBG_TIME,
        //     gk_startcputimer((*ctrl).ProjectTmr),
        // );
        assert!((*graph).vwgt != ptr::null_mut());

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
        && FindPartitionInducedComponents(graph, (*graph).where_, ptr::null_mut(), ptr::null_mut())
            > (*ctrl).nparts
    {
        EliminateComponents(ctrl, graph);
    }

    if IsBalanced(ctrl, graph, 0.0) == 0 {
        ComputeKWayBoundary(ctrl, graph, BNDTYPE_BALANCE);
        Greedy_KWayOptimize(ctrl, graph, 10, 0.0, OMODE_BALANCE);

        ComputeKWayBoundary(ctrl, graph, BNDTYPE_REFINE);
        Greedy_KWayOptimize(ctrl, graph, (*ctrl).niter, 0.0, OMODE_REFINE);
    }

    if (*ctrl).contig != 0 {
        assert!(
            FindPartitionInducedComponents(
                graph,
                (*graph).where_,
                ptr::null_mut(),
                ptr::null_mut()
            ) == (*ctrl).nparts,
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
pub fn AllocateKWayPartitionMemory(ctrl: *mut ctrl_t, graph: *mut graph_t) -> () {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    graph.pwgts = imalloc(
        (ctrl.nparts * graph.ncon) as usize,
        "AllocateKWayPartitionMemory: pwgts\0".as_ptr()
    ) as _;
    graph.where_ = imalloc(
        graph.nvtxs as usize,
        "AllocateKWayPartitionMemory: where\0".as_ptr(),
    ) as _;
    graph.bndptr = imalloc(
        graph.nvtxs as usize,
        "AllocateKWayPartitionMemory: bndptr\0".as_ptr(),
    ) as _;
    graph.bndind = imalloc(
        graph.nvtxs as usize,
        "AllocateKWayPartitionMemory: bndind\0".as_ptr(),
    ) as _;

    match ctrl.objtype {
        METIS_OBJTYPE_CUT => {
            graph.ckrinfo = gk_malloc(
                graph.nvtxs as usize * std::mem::size_of::<ckrinfo_t>(),
                "AllocateKWayPartitionMemory: ckrinfo\0".as_ptr(),
            ) as *mut _;
        }

        METIS_OBJTYPE_VOL => {
            graph.vkrinfo = gk_malloc(
                graph.nvtxs as usize * std::mem::size_of::<vkrinfo_t>(),
                "AllocateKWayVolPartitionMemory: vkrinfo\0".as_ptr(),
            ) as *mut _;

            /* This is to let the cut-based -minconn and -contig large-scale graph
            changes to go through */
            // wtf this is terrible??? they aren't even the same size??
            // TODO: figure out tf is happening here
            graph.ckrinfo = graph.vkrinfo as *mut _;
        }

        _ => {
            panic!("Unknown objtype of {}", ctrl.objtype);
        }
    };
}

/// This function computes the initial id/ed  for cut-based partitioning
#[metis_func]
pub extern "C" fn ComputeKWayPartitionParams(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    eprintln!("Called ComputeKWayPartitionParams");
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    let nvtxs: idx_t;
    let ncon: idx_t;
    let nparts: idx_t;
    let mut nbnd: idx_t;
    let mut mincut: idx_t;

    nparts = ctrl.nparts;

    nvtxs = graph.nvtxs;
    ncon = graph.ncon;

    get_graph_slices!(ctrl, graph => xadj adjncy vwgt where_ adjwgt);
    get_graph_slices_mut!(ctrl, graph => pwgts bndptr bndind ckrinfo);
    pwgts.fill(0);
    bndptr.fill(-1);

    nbnd = 0;
    mincut = 0;

    /* Compute pwgts (Gavin: parition weights?) */
    if ncon == 1 {
        for i in 0..(nvtxs as usize) {
            assert!(where_[i] >= 0 && where_[i] < nparts);
            pwgts[where_[i] as usize] += vwgt[i];
        }
    } else {
        for i in 0..(nvtxs as usize) {
            let me = where_[i];
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
                            (*ctrl).cnbrpool.add(myrinfo.inbr as usize),
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

                    assert!(myrinfo.nnbrs <= xadj[i + 1] - xadj[i]);

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

            assert!(debug::CheckBnd2(graph) != 0);
        }

        METIS_OBJTYPE_VOL => {
            {
                // memset((*graph).vkrinfo, 0, sizeof(vkrinfo_t) * nvtxs);
                // cnbrpool operations may realloc, so we can't leave this be
                let amyrinfo: &mut [vkrinfo_t] =
                    slice::from_raw_parts_mut(graph.vkrinfo, nvtxs as usize);
                amyrinfo.fill_with(std::default::Default::default);
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

                /* Time to compute the particular external degrees */
                if myrinfo.ned > 0 {
                    mincut += myrinfo.ned;

                    myrinfo.inbr = vnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
                    let mynbrs = slice::from_raw_parts_mut(
                        ctrl.vnbrpool.add(myrinfo.inbr as usize),
                        myrinfo.nnbrs as usize,
                    );

                    for j in xadj[0]..xadj[i + 1] {
                        let j = j as usize;
                        let other = where_[adjncy[j] as usize];
                        if me != other {
                            let mut k = 0;
                            for kk in 0..(myrinfo.nnbrs as usize) {
                                k = kk; // used after loop
                                if mynbrs[k].pid == other {
                                    mynbrs[k].ned += 1;
                                    break;
                                }
                            }
                            if k == myrinfo.nnbrs as usize {
                                mynbrs[k].gv = 0;
                                mynbrs[k].pid = other;
                                mynbrs[k].ned = 1;
                                myrinfo.nnbrs += 1;
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
            assert!(graph.minvol == ComputeVolume(graph, graph.where_));
        }
        _ => panic!("Unknown objtype of {}", ctrl.objtype),
    }
}

/// This function projects a partition, and at the same time computes the parameters for refinement.
#[metis_func]
pub extern "C" fn ProjectKWayPartition(ctrl: *mut ctrl_t, graph: *mut graph_t) -> () {
    eprintln!("Called ProjectKWayPartition");
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let nvtxs: idx_t;
    let mut nbnd: idx_t;
    let nparts: idx_t;

    let cgraph: *mut graph_t;
    let dropedges;

    dropedges = (*ctrl).dropedges;

    nparts = (*ctrl).nparts;

    cgraph = (*graph).coarser;
    // TODO: I think this may be should be cgraph.nvtxs
    mkslice_mut!(cwhere: cgraph->where_, (*graph).nvtxs as usize);

    if (*ctrl).objtype == METIS_OBJTYPE_CUT {
        assert!(debug::CheckBnd2(cgraph) != 0);
    } else {
        assert!((*cgraph).minvol == ComputeVolume(cgraph, (*cgraph).where_));
    }

    /* free the coarse graph's structure (reduce maxmem) */
    FreeSData(cgraph);

    nvtxs = (*graph).nvtxs;

    AllocateKWayPartitionMemory(ctrl, graph);

    get_graph_slices!(ctrl, *graph => xadj adjncy adjwgt);
    get_graph_slices_mut!(ctrl, *graph => cmap bndind bndptr where_);
    bndptr.fill(-1);

    let mut htable: Vec<idx_t> = vec![-1; nparts as usize];

    /* Compute the required info for refinement */
    match (*ctrl).objtype {
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

                let myrinfo = (*graph).ckrinfo.add(i).as_mut().unwrap();

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

            (*graph).nbnd = nbnd;
            assert!(debug::CheckBnd2(graph) != 0);
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
                let myrinfo = (*graph).vkrinfo.add(i).as_mut().unwrap();

                if cmap[i] == 0 {
                    /* Note that cmap[i] = crinfo[cmap[i]].ed */
                    myrinfo.nid = (iend - istart) as idx_t;
                    myrinfo.inbr = -1;
                } else {
                    /* Potentially an interface node */
                    myrinfo.inbr = vnbrpoolGetNext(ctrl, (iend - istart) as idx_t);
                    // mynbrs = (*ctrl).vnbrpool + myrinfo.inbr;
                    let mynbrs = slice::from_raw_parts_mut(
                        (*ctrl).vnbrpool.add(myrinfo.inbr as usize),
                        nvtxs as usize - myrinfo.inbr as usize,
                    );

                    let me = where_[i];
                    let mut tid = 0;
                    let mut ted = 0;
                    for j in istart..iend {
                        let other = where_[adjncy[j] as usize];
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
                        (*ctrl).nbrpoolcpos -= (nparts as usize).min(iend - istart);
                        myrinfo.inbr = -1;
                    } else {
                        for j in 0..myrinfo.nnbrs {
                            htable[mynbrs[j as usize].pid as usize] = -1;
                        }
                    }
                }
            }

            ComputeKWayVolGains(ctrl, graph);

            assert!((*graph).minvol == ComputeVolume(graph, (*graph).where_));
        }

        _ => panic!("Unknown objtype of {}", (*ctrl).objtype),
    }

    (*graph).mincut = if dropedges != 0 {
        ComputeCut(graph, where_.as_mut_ptr())
    } else {
        (*cgraph).mincut
    };
    {
        mkslice_mut!(dst: graph->pwgts, nparts * (*graph).ncon);
        mkslice_mut!(src: cgraph->pwgts, nparts * (*graph).ncon);
        dst.copy_from_slice(src);
        // remember: icopy(n, src, dest) -> reverse of memcpy
        // icopy(nparts * (*graph).ncon, (*cgraph).pwgts, (*graph).pwgts);
    }

    FreeGraph((&mut graph.coarser) as *mut *mut graph_t);
}

/// This function computes the boundary definition for balancing.
#[metis_func]
pub extern "C" fn ComputeKWayBoundary(ctrl: *mut ctrl_t, graph: *mut graph_t, bndtype: idx_t) {
    eprintln!("Called ComputeKWayBoundary");
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
    eprintln!("Called ComputeKWayVolGains");
    let ctrl: &mut ctrl_t = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    let nparts: idx_t;

    nparts = (*ctrl).nparts;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj vsize adjncy);
    get_graph_slices_mut!(ctrl, graph => where_ bndind bndptr vkrinfo);

    // bndptr = iset(nvtxs, -1, graph.bndptr);
    bndptr.fill(-1);

    // ophtable = iset(nparts, -1, iwspacemalloc(ctrl, nparts));
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
            // mynbrs = (*ctrl).vnbrpool + myrinfo.inbr;
            let mynbrs = slice::from_raw_parts_mut(
                (*ctrl).vnbrpool.add(myrinfo.inbr as usize),
                ctrl.nbrpoolsize - myrinfo.inbr as usize,
            );

            graph.minvol += myrinfo.nnbrs * vsize[i];

            for j in xadj[i]..xadj[i + 1] {
                let j = j as usize;
                let ii = adjncy[j] as usize;
                let other = where_[ii];
                let orinfo = vkrinfo[ii].clone();
                let onbrs = slice::from_raw_parts_mut(
                    ctrl.vnbrpool.add(orinfo.inbr as usize),
                    orinfo.nnbrs as usize,
                );

                for k in 0..orinfo.nnbrs {
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
                    let me = me as usize;
                    assert!(ophtable[me] != -1);

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
    return (ComputeLoadImbalanceDiff(graph, (*ctrl).nparts, (*ctrl).pijbm, (*ctrl).ubfactors)
        <= ffactor) as _;
}
