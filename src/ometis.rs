/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ometis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id: ometis.c 10513 2011-07-07 22:06:03Z karypis $
 *
 */

use std::ffi::c_int;

use crate::*;

/*************************************************************************/
/* This function is the entry point for the multilevel nested dissection
    ordering code. At each bisection, a node-separator is computed using
    a node-based refinement approach.

    \param nvtxs is the number of vertices in the graph.
    \param xadj is of length nvtxs+1 marking the start of the adjancy
           list of each vertex in adjncy.
    \param adjncy stores the adjacency lists of the vertices. The adjnacy
           list of a vertex should not contain the vertex itself.
    \param vwgt is an array of size nvtxs storing the weight of each
           vertex. If vwgt is std::ptr::null_mut(), then the vertices are considered
           to have unit weight.
    \param numflag is either 0 or 1 indicating that the numbering of
           the vertices starts from 0 or 1, respectively.
    \param options is an array of size METIS_NOPTIONS used to pass
           various options impacting the of the algorithm. A std::ptr::null_mut()
           value indicates use of default options.
    \param perm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A'[i] = A[perm[i]].
    \param iperm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A[i] = A'[iperm[i]].
*/
/*************************************************************************/
#[metis_func(no_pfx)]
pub extern "C" fn METIS_NodeND(
    nvtxs: *const idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    options: *const idx_t,
    perm: *mut idx_t,
    iperm: *mut idx_t,
) -> c_int {
    // int sigrval=0, renumber=0;
    // idx_t i, ii, j, l, nnvtxs=0;
    // graph_t *graph=std::ptr::null_mut();
    // ctrl_t *ctrl;
    // idx_t *cptr, *cind, *piperm;
    // int numflag = 0;

    /* set up malloc cleaning code and signal catchers */
    if gk_malloc_init() == 0 {
        eprintln!("metis memory error");
        return METIS_ERROR_MEMORY;
    }

    let nvtxs = *nvtxs as usize;
    let mut nnvtxs: usize = 0;

    mkslice_mut!(perm, nvtxs);
    mkslice_mut!(iperm, nvtxs);

    //   gk_sigtrap();
    //
    //   if ((sigrval = gk_sigcatch()) != 0)  {
    //     goto SIGTHROW;
    // }

    /* set up the run time parameters */
    let ctrl = options::SetupCtrl(
        METIS_OP_OMETIS,
        options,
        1,
        3,
        std::ptr::null_mut(),
        std::ptr::null_mut(),
    );

    let Some(mut ctrl) = ctrl.as_mut() else {
        eprintln!("failed setting up ctrl");
        // gk_siguntrap();
        return METIS_ERROR_INPUT;
    };

    /* if required, change the numbering to 0 */
    // if (ctrl.numflag == 1) {
    //     Change2CNumbering(*nvtxs, xadj, adjncy);
    //     renumber = 1;
    // }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

    let mut graph = std::ptr::null_mut();
    let mut cptr: Vec<idx_t> = vec![];
    let mut cind: Vec<idx_t> = vec![];
    let mut piperm: Vec<idx_t> = vec![];
    /* prune the dense columns */
    if ctrl.pfactor > 0.0 {
        piperm = vec![0; nvtxs];

        graph = compress::PruneGraph(
            ctrl,
            nvtxs as idx_t,
            xadj,
            adjncy,
            vwgt,
            piperm.as_mut_ptr(),
            ctrl.pfactor,
        );
        if graph.is_null() {
            /* if there was no prunning, cleanup the pfactor */
            std::mem::take(&mut piperm);
            ctrl.pfactor = 0.0;
        } else {
            nnvtxs = (*graph).nvtxs as usize;
            ctrl.compress = 0; /* disable compression if prunning took place */
        }
    }

    if ctrl.compress != 0 {
        /* compress the graph; note that compression only happens if not prunning
        has taken place. */
        cptr = vec![0; nvtxs + 1];
        cind = vec![0; nvtxs];

        graph = compress::CompressGraph(
            ctrl,
            nvtxs as idx_t,
            xadj,
            adjncy,
            vwgt,
            cptr.as_mut_ptr(),
            cind.as_mut_ptr(),
        );
        if graph.is_null() {
            /* if there was no compression, cleanup the compress flag */
            std::mem::take(&mut cptr);
            std::mem::take(&mut cind);
            ctrl.compress = 0;
        } else {
            nnvtxs = (*graph).nvtxs as usize;
            ctrl.cfactor = 1.0 * (nvtxs as real_t) / (nnvtxs as real_t);
            if ctrl.cfactor > 1.5 && ctrl.nseps == 1 {
                ctrl.nseps = 2;
            }
            //ctrl.nseps = (idx_t)(ctrl.cfactor*ctrl.nseps);
        }
    }

    if ctrl.pfactor == 0.0 && ctrl.compress == 0 {
        /* if no prunning and no compression, setup the graph in the normal way. */
        graph = graph::SetupGraph(
            ctrl,
            nvtxs as idx_t,
            1,
            xadj,
            adjncy,
            vwgt,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );
    }

    let graph = graph.as_mut().unwrap();

    debug_assert!(checkgraph::CheckGraph(graph, ctrl.numflag, 1) != 0);

    /* allocate workspace memory */
    wspace::AllocateWorkSpace(ctrl, graph);

    /* do the nested dissection ordering  */
    if ctrl.ccorder != 0 {
        MlevelNestedDissectionCC(ctrl, graph, iperm.as_mut_ptr(), graph.nvtxs);
    } else {
        MlevelNestedDissection(ctrl, graph, iperm.as_mut_ptr(), graph.nvtxs);
    }

    if ctrl.pfactor > 0.0 {
        /* Order any prunned vertices */
        /* Use perm as an auxiliary array */
        let perm = &mut perm[..nnvtxs];
        perm.copy_from_slice(&iperm[..nnvtxs]);

        for i in (0)..(nnvtxs) {
            iperm[piperm[i as usize] as usize] = perm[i as usize];
        }
        for i in (nnvtxs)..(nvtxs) {
            iperm[piperm[i as usize] as usize] = i as idx_t;
        }

        std::mem::take(&mut piperm);
    } else if ctrl.compress != 0 {
        /* Uncompress the ordering */
        /* construct perm from iperm */
        for i in (0)..(nnvtxs) {
            perm[iperm[i as usize] as usize] = i as idx_t;
        }
        let mut l = 0;
        for ii in 0..nnvtxs {
            let i = perm[ii as usize];
            for j in (cptr[i as usize])..(cptr[(i + 1) as usize]) {
                iperm[cind[j as usize] as usize] = l;
                l += 1;
            }
        }

        std::mem::take(&mut cptr);
        std::mem::take(&mut cind);
    }

    for i in (0)..(nvtxs) {
        perm[iperm[i as usize] as usize] = i as idx_t;
    }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

    /* clean up */
    options::FreeCtrl((&raw mut ctrl).cast());

    // SIGTHROW:
    /* if required, change the numbering back to 1 */
    // I'm not doing fortran renaming
    //   if (renumber) {
    //     Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
    // }

    // gk_siguntrap();
    gk_malloc_cleanup(0);

    return util::metis_rcode(0);
}

/*************************************************************************/
/* This is the driver for the recursive tri-section of a graph into the
   left, separator, and right partitions. The graphs correspond to the
   left and right parts are further tri-sected in a recursive fashion.
   The nodes in the separator are ordered at the end of the left & right
   nodes.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNestedDissection(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    order: *mut idx_t,
    lastvtx: idx_t,
) {
    let mut lastvtx = lastvtx;
    // idx_t i, j, nvtxs, nbnd;
    // idx_t *label, *bndind;
    // graph_t *lgraph, *rgraph;

    let mut graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    MlevelNodeBisectionMultiple(ctrl, graph);

    // ifset!(
    //     ctrl.dbglvl,
    //     METIS_DBG_SEPINFO,
    //     print!(
    //         "Nvtxs: {:6}, [({:6} {:6} {:6}) as usize]\n",
    //         graph.nvtxs, graph.pwgts[0], graph.pwgts[1], graph.pwgts[2]
    //     )
    // );

    /* Order the nodes in the separator */
    let nbnd = graph.nbnd as usize;
    get_graph_slices!(graph => bndind label);
    let label: &[idx_t] = label; // this is a nop but rust_analyzer is confused
    for &ind in &bndind[..nbnd] {
        lastvtx -= 1;
        // I'm not able to make this a slice since this is called recursively and the size of order
        // is determined by the call
        *order.add(label[ind as usize] as usize) = lastvtx;
    }

    let mut lgraph = std::ptr::null_mut();
    let mut rgraph = std::ptr::null_mut();
    SplitGraphOrder(ctrl, graph, &mut lgraph, &mut rgraph);
    let mut lgraph = lgraph.as_mut().unwrap();
    let mut rgraph = rgraph.as_mut().unwrap();

    /* Free the memory of the top level graph */
    graph::FreeGraph((&raw mut graph).cast());

    /* Recurse on lgraph first, as its lastvtx depends on rgraph.nvtxs, which
    will not be defined upon return from MlevelNestedDissection. */
    if lgraph.nvtxs > MMDSWITCH && lgraph.nedges > 0 {
        MlevelNestedDissection(ctrl, lgraph, order, lastvtx - rgraph.nvtxs);
    } else {
        MMDOrder(ctrl, lgraph, order, lastvtx - rgraph.nvtxs);
        graph::FreeGraph((&raw mut lgraph).cast());
    }
    if rgraph.nvtxs > MMDSWITCH && rgraph.nedges > 0 {
        MlevelNestedDissection(ctrl, rgraph, order, lastvtx);
    } else {
        MMDOrder(ctrl, rgraph, order, lastvtx);
        graph::FreeGraph((&raw mut rgraph).cast());
    }
}

/*************************************************************************/
/* This routine is similar to its non 'CC' counterpart. The difference is
    that after each tri-section, the connected components of the original
    graph that result after removing the separator vertises are ordered
    independently (i.e., this may lead to more than just the left and
    the right subgraphs).
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNestedDissectionCC(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    order: *mut idx_t,
    lastvtx: idx_t,
) {
    let mut lastvtx = lastvtx;
    // idx_t i, j, nvtxs, nbnd, ncmps, rnvtxs, snvtxs;
    // idx_t *label, *bndind;
    // idx_t *cptr, *cind;
    // graph_t **sgraphs;

    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    let nvtxs = graph.nvtxs as usize;

    MlevelNodeBisectionMultiple(ctrl, graph);

    // ifset!(
    //     ctrl.dbglvl,
    //     METIS_DBG_SEPINFO,
    //     print!(
    //         "Nvtxs: {:6}, [({:6} {:6} {:6}) as usize]\n",
    //         graph.nvtxs, graph.pwgts[0], graph.pwgts[1], graph.pwgts[2]
    //     )
    // );

    /* Order the nodes in the separator */
    let nbnd = graph.nbnd as usize;
    get_graph_slices!(graph => bndind label);
    // bndind = graph.bndind;
    // label = graph.label;
    for i in (0)..(nbnd) {
        lastvtx -= 1;
        // I'm not able to make this a slice since this is called recursively and the size of order
        // is determined by the call
        *order.add(label[bndind[i as usize] as usize] as usize) = lastvtx;
    }

    // WCOREPUSH;
    let mut cptr = vec![0; nvtxs + 1 as usize];
    let mut cind = vec![0; nvtxs as usize];
    let ncmps = contig::FindSepInducedComponents(ctrl, graph, cptr.as_mut_ptr(), cind.as_mut_ptr())
        as usize;

    // if ctrl.dbglvl & METIS_DBG_INFO != 0 {
    //     if ncmps > 2 {
    //         print!("  Bisection resulted in {:} connected components\n", ncmps);
    //     }
    // }

    let sgraphs = SplitGraphOrderCC(
        ctrl,
        graph,
        ncmps as idx_t,
        cptr.as_mut_ptr(),
        cind.as_mut_ptr(),
    );
    mkslice_mut!(sgraphs, ncmps);

    // WCOREPOP;

    /* Free the memory of the top level graph */
    let mut graph: *mut _ = graph;
    graph::FreeGraph((&raw mut graph).cast());

    /* Go and process the subgraphs */
    let mut rnvtxs = 0;
    for i in 0..ncmps {
        /* Save the number of vertices in sgraphs[i] because sgraphs[i] is freed
        inside MlevelNestedDissectionCC, and as such it will be undefined. */
        // Porting: this is why not making &mut
        let snvtxs = (*sgraphs[i as usize]).nvtxs;

        if snvtxs > MMDSWITCH && (*sgraphs[i as usize]).nedges > 0 {
            MlevelNestedDissectionCC(
                ctrl,
                sgraphs[i as usize],
                order,
                lastvtx - rnvtxs,
            );
        } else {
            MMDOrder(
                ctrl,
                sgraphs[i as usize],
                order,
                lastvtx - rnvtxs,
            );
            graph::FreeGraph(&mut sgraphs[i as usize]);
        }
        rnvtxs += snvtxs;
    }

    gk::free_slice_guard(sgraphs);
}

/*************************************************************************/
/* This function performs multilevel node bisection (i.e., tri-section).
It performs multiple bisections and selects the best. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNodeBisectionMultiple(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    // idx_t i, mincut;
    // idx_t *bestwhere_;
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    /* if the graph is small, just find a single vertex separator */
    if ctrl.nseps == 1 || graph.nvtxs < (if ctrl.compress != 0 { 1000 } else { 2000 }) {
        MlevelNodeBisectionL2(ctrl, graph, LARGENIPARTS);
        return;
    }

    // WCOREPUSH;

    let mut bestwhere_ = vec![0; graph.nvtxs as usize];

    let mut mincut = *graph.tvwgt;
    for i in (0)..(ctrl.nseps) {
        MlevelNodeBisectionL2(ctrl, graph, LARGENIPARTS);
        get_graph_slices!(graph => where_);

        if i == 0 || graph.mincut < mincut {
            mincut = graph.mincut;
            if i < ctrl.nseps - 1 {
                // icopy(graph.nvtxs, graph.where_, bestwhere_);
                bestwhere_.copy_from_slice(where_);
            }
        }

        if mincut == 0 {
            break;
        }

        if i < ctrl.nseps - 1 {
            graph::FreeRData(graph);
        }
    }

    if mincut != graph.mincut {
        get_graph_slices_mut!(graph => where_);
        where_.copy_from_slice(&bestwhere_);
        // icopy(graph.nvtxs, bestwhere_, graph.where_);
        srefine::Compute2WayNodePartitionParams(ctrl, graph);
    }

    // WCOREPOP;
}

pub const MLEVEL_NODE_BISECTION_L2_CUTOFF: usize = 5000;

/*************************************************************************/
/* This function performs multilevel node bisection (i.e., tri-section).
It performs multiple bisections and selects the best. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNodeBisectionL2(ctrl: *mut ctrl_t, graph: *mut graph_t, niparts: idx_t) {
    // idx_t i, mincut, nruns=5;
    // graph_t *cgraph;
    // idx_t *bestwhere_;
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    /* if the graph is small, just find a single vertex separator */
    if graph.nvtxs < MLEVEL_NODE_BISECTION_L2_CUTOFF as idx_t {
        MlevelNodeBisectionL1(ctrl, graph, niparts);
        return;
    }

    // WCOREPUSH;

    // ctrl.CoarsenTo = gk_max(100, graph.nvtxs / 30);
    ctrl.CoarsenTo = (graph.nvtxs / 30).max(100);

    let cgraph = coarsen::CoarsenGraphNlevels(ctrl, graph, 4);
    let cgraph = cgraph.as_mut().unwrap();

    let mut bestwhere_ = vec![0; cgraph.nvtxs as usize];

    let mut mincut = *graph.tvwgt;
    let nruns = 5;
    for i in (0)..(nruns) {
        MlevelNodeBisectionL1(ctrl, cgraph, (0.7 * niparts as real_t) as idx_t);

        if i == 0 || cgraph.mincut < mincut {
            mincut = cgraph.mincut;
            if i < nruns - 1 {
                get_graph_slices!(cgraph => where_);
                bestwhere_.copy_from_slice(where_);
                // icopy(cgraph.nvtxs, cgraph.where_, bestwhere_);
            }
        }

        if mincut == 0 {
            break;
        }

        if i < nruns - 1 {
            graph::FreeRData(cgraph);
        }
    }

    if mincut != cgraph.mincut {
        // icopy(cgraph.nvtxs, bestwhere_, cgraph.where_);
        get_graph_slices_mut!(cgraph => where_);
        where_.copy_from_slice(&bestwhere_);
    }

    // WCOREPOP;

    srefine::Refine2WayNode(ctrl, graph, cgraph);
}

/*************************************************************************/
/* The top-level routine of the actual multilevel node bisection */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNodeBisectionL1(ctrl: *mut ctrl_t, graph: *mut graph_t, niparts: idx_t) {
    // graph_t *cgraph;
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    ctrl.CoarsenTo = (graph.nvtxs / 8).clamp(40, 100);

    let cgraph = coarsen::CoarsenGraph(ctrl, graph);
    let cgraph = cgraph.as_mut().unwrap();

    let niparts = (if cgraph.nvtxs <= ctrl.CoarsenTo {
        niparts / 2
    } else {
        niparts
    })
    .max(1);
    /*niparts = (cgraph.nvtxs <= ctrl.CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);*/
    initpart::InitSeparator(ctrl, cgraph, niparts);

    srefine::Refine2WayNode(ctrl, graph, cgraph);
}

/*************************************************************************/
/* This function takes a graph and a tri-section (left, right, separator)
    and splits it into two graphs.

    This function relies on the fact that adjwgt is all equal to 1.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn SplitGraphOrder(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    r_lgraph: *mut *mut graph_t,
    r_rgraph: *mut *mut graph_t,
) {
    // idx_t i, ii, j, k, l, istart, iend, mypart, nvtxs, snvtxs[3], snedges[3];
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label, *where_, *bndptr, *bndind;
    // idx_t *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *slabel[2];
    // idx_t *rename;
    // idx_t *auxadjncy;
    // graph_t *lgraph, *rgraph;

    // WCOREPUSH;

    let graph = graph.as_mut().unwrap();

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.SplitTmr));

    let nvtxs = graph.nvtxs as usize;
    // xadj = graph.xadj;
    // vwgt = graph.vwgt;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // label = graph.label;
    // where_ = graph.where_;
    // bndptr = graph.bndptr;
    // bndind = graph.bndind;
    get_graph_slices!(graph => xadj vwgt adjncy label where_ bndind);
    get_graph_slices_mut!(graph => bndptr);

    let mut rename: Vec<idx_t> = vec![0; nvtxs as usize];

    let mut snvtxs = [0_usize; 3];
    let mut snedges = [0_usize; 3];
    for i in (0)..(nvtxs) {
        let k = where_[i as usize] as usize;
        rename[i as usize] = snvtxs[k] as idx_t;
        snvtxs[k] += 1;
        snedges[k] += xadj[(i + 1) as usize] as usize - xadj[i as usize] as usize;
    }

    let lgraph = graph::SetupSplitGraph(graph, snvtxs[0] as idx_t, snedges[0] as idx_t);
    let lgraph = lgraph.as_mut().unwrap();

    let rgraph = graph::SetupSplitGraph(graph, snvtxs[1] as idx_t, snedges[1] as idx_t);
    let rgraph = rgraph.as_mut().unwrap();
    let sxadj = [
        std::slice::from_raw_parts_mut(lgraph.xadj, snvtxs[0] + 1),
        std::slice::from_raw_parts_mut(rgraph.xadj, snvtxs[1] + 1),
    ];
    let sadjncy = [
        std::slice::from_raw_parts_mut(lgraph.adjncy, snedges[0]),
        std::slice::from_raw_parts_mut(rgraph.adjncy, snedges[1]),
    ];
    let sadjwgt = [
        std::slice::from_raw_parts_mut(lgraph.adjwgt, snedges[0]),
        std::slice::from_raw_parts_mut(rgraph.adjwgt, snedges[1]),
    ];
    let ncon = graph.ncon as usize;
    let svwgt = [
        std::slice::from_raw_parts_mut(lgraph.vwgt, snvtxs[0] * ncon),
        std::slice::from_raw_parts_mut(rgraph.vwgt, snvtxs[1] * ncon),
    ];
    let slabel = [
        std::slice::from_raw_parts_mut(lgraph.label, snvtxs[0]),
        std::slice::from_raw_parts_mut(rgraph.label, snvtxs[1]),
    ];

    /* Go and use bndptr to also mark the boundary nodes in the two partitions */
    for ii in (0)..(graph.nbnd) {
        let i = bndind[ii as usize] as usize;
        for j in (xadj[i])..(xadj[i + 1]) {
            bndptr[adjncy[j as usize] as usize] = 1;
        }
    }

    snvtxs[0] = 0;
    snvtxs[1] = 0;
    snedges[0] = 0;
    snedges[1] = 0;
    sxadj[0][0] = 0;
    sxadj[1][0] = 0;
    for i in (0)..(nvtxs) {
        let mypart = where_[i as usize] as usize;
        if mypart == 2 {
            continue;
        }

        let istart = xadj[i as usize] as usize;
        let iend = xadj[(i + 1) as usize] as usize;
        if bndptr[i as usize] == -1 {
            /* This is an interior vertex */
            // let auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
            // for j in istart..iend {
            //     auxadjncy[j as usize] = adjncy[j as usize];
            // }
            let len = iend - istart;
            sadjncy[mypart][cntrng!(snedges[mypart], len)].copy_from_slice(&adjncy[istart..iend]);
            snedges[mypart] += iend - istart;
        } else {
            let auxadjncy = &mut *sadjncy[mypart];
            let mut l = snedges[mypart];
            for j in (istart)..(iend) {
                let k = adjncy[j as usize];
                if where_[k as usize] == mypart as idx_t {
                    auxadjncy[(l) as usize] = k;
                    l += 1;
                }
            }
            snedges[mypart] = l;
        }

        svwgt[mypart][snvtxs[mypart] as usize] = vwgt[i as usize];
        slabel[mypart][snvtxs[mypart] as usize] = label[i as usize];
        snvtxs[mypart] += 1;
        sxadj[mypart][snvtxs[mypart] as usize] = snedges[mypart as usize] as idx_t;
    }

    for mypart in (0)..(2) {
        let iend = snedges[mypart as usize];
        sadjwgt[mypart][..iend].fill(1);
        // iset(iend, 1, sadjwgt[mypart as usize]);

        for i in (0)..(iend) {
            sadjncy[mypart][i] = rename[sadjncy[mypart][i] as usize];
        }
    }

    lgraph.nvtxs = snvtxs[0] as idx_t;
    lgraph.nedges = snedges[0] as idx_t;
    rgraph.nvtxs = snvtxs[1] as idx_t;
    rgraph.nedges = snedges[1] as idx_t;

    graph::SetupGraph_tvwgt(lgraph);
    graph::SetupGraph_tvwgt(rgraph);

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.SplitTmr));

    *r_lgraph = lgraph;
    *r_rgraph = rgraph;

    // WCOREPOP;
}

/*************************************************************************/
/* This function takes a graph and generates a set of graphs, each of
    which is a connected component in the original graph.

    This function relies on the fact that adjwgt is all equal to 1.

    \param ctrl stores run state info.
    \param graph is the graph to be split.
    \param ncmps is the number of connected components.
    \param cptr is an array of size ncmps+1 that marks the start and end
           locations of the vertices in cind that make up the respective
           components (i.e., cptr, cind is in CSR format).
    \param cind is an array of size equal to the number of vertices in
           the original graph and stores the vertices that belong to each
           connected component.

    \returns an array of subgraphs corresponding to the extracted subgraphs.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn SplitGraphOrderCC(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ncmps: idx_t,
    cptr: *mut idx_t,
    cind: *mut idx_t,
) -> *mut *mut graph_t {
    // idx_t i, ii, iii, j, k, l, istart, iend, mypart, nvtxs, snvtxs, snedges;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label, *where_, *bndptr, *bndind;
    // idx_t *sxadj, *svwgt, *sadjncy, *sadjwgt, *slabel;
    // idx_t *rename;
    // idx_t *auxadjncy;
    // graph_t **sgraphs;

    // WCOREPUSH;
    let graph = graph.as_mut().unwrap();

    let ncmps = ncmps as usize;

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.SplitTmr));

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj vwgt adjncy label where_ bndind);
    get_graph_slices_mut!(graph => bndptr);
    let ncmps = ncmps as usize;

    mkslice_mut!(cptr, nvtxs + 1);
    mkslice_mut!(cind, nvtxs);

    /* Go and use bndptr to also mark the boundary nodes in the two partitions */
    for ii in (0)..(graph.nbnd) {
        let i = bndind[ii as usize] as usize;
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            bndptr[adjncy[j as usize] as usize] = 1;
        }
    }

    let mut rename: Vec<idx_t> = vec![0; nvtxs as usize];

    let sgraphs_ret: *mut *mut graph_t = gk_malloc(
        size_of::<*mut graph_t>() * ncmps,
        c"SplitGraphOrderCC: sgraphs".as_ptr(),
    )
    .cast();
    mkslice_mut!(sgraphs: sgraphs_ret, ncmps);

    /* Go and split the graph a component at a time */
    for iii in (0)..(ncmps) {
        irandArrayPermute(
            cptr[iii + 1] - cptr[iii],
            cind[cptr[iii as usize] as usize..].as_mut_ptr(),
            cptr[iii + 1] - cptr[iii],
            0,
        );
        let mut snvtxs: usize = 0;
        let mut snedges: idx_t = 0;
        for j in (cptr[iii as usize])..(cptr[(iii + 1) as usize]) {
            let i = cind[j as usize];
            rename[i as usize] = snvtxs as idx_t;
            snvtxs += 1;
            snedges += xadj[(i + 1) as usize] - xadj[i as usize];
        }

        sgraphs[iii as usize] = graph::SetupSplitGraph(graph, snvtxs as idx_t, snedges as idx_t);
        let sgraphs_iii = sgraphs[iii].as_mut().unwrap();

        // let sxadj = sgraphs_iii.xadj;
        // let svwgt = sgraphs_iii.vwgt;
        // let sadjncy = sgraphs_iii.adjncy;
        // let sadjwgt = sgraphs_iii.adjwgt;
        // let slabel = sgraphs_iii.label;
        let sxadj = get_graph_slice_mut!(sgraphs_iii => xadj);
        let svwgt = get_graph_slice_mut!(sgraphs_iii => vwgt);
        // haven't populated xadj yet, so we need to do this manually
        let sadjncy = std::slice::from_raw_parts_mut(sgraphs_iii.adjncy, snedges as usize);
        let sadjwgt = std::slice::from_raw_parts_mut(sgraphs_iii.adjwgt, snedges as usize);
        let slabel = get_graph_slice_mut!(sgraphs_iii => label);

        let mut snvtxs = 0;
        let mut snedges = 0;
        sxadj[0] = 0;
        for ii in (cptr[iii as usize])..(cptr[(iii + 1) as usize]) {
            let ii = ii as usize;
            let i = cind[ii as usize] as usize;

            let istart = xadj[i] as usize;
            let iend = xadj[i + 1] as usize;
            if bndptr[i as usize] == -1 {
                /* This is an interior vertex */
                // auxadjncy = sadjncy + snedges - istart;
                // for j in istart..iend {
                //     auxadjncy[j as usize] = adjncy[j as usize];
                // }
                let len = iend - istart;
                sadjncy[cntrng!(snedges, len)].copy_from_slice(&adjncy[istart..iend]);
                snedges += len;
            } else {
                let mut l = snedges;
                for j in (istart)..(iend) {
                    let k = adjncy[j as usize];
                    if where_[k as usize] != 2 {
                        sadjncy[(l) as usize] = k;
                        l += 1;
                    }
                }
                snedges = l;
            }
            svwgt[snvtxs] = vwgt[i];
            slabel[snvtxs] = label[i];
            snvtxs += 1;
            sxadj[snvtxs] = snedges as idx_t;
        }

        // iset(snedges, 1, sadjwgt);
        sadjwgt[..snedges].fill(1);

        // for i in (0)..(snedges) {
        //     sadjncy[i as usize] = rename[sadjncy[i as usize] as usize];
        // }
        for adj in &mut sadjncy[..snedges] {
            *adj = rename[*adj as usize];
        }

        sgraphs_iii.nvtxs = snvtxs as idx_t;
        sgraphs_iii.nedges = snedges as idx_t;

        graph::SetupGraph_tvwgt(sgraphs_iii);
    }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.SplitTmr));

    // WCOREPOP;

    return sgraphs_ret;
}

/*************************************************************************/
/* This function uses MMD to order the graph. The vertices are numbered
from lastvtx downwards. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MMDOrder(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    order: *mut idx_t,
    lastvtx: idx_t,
) {
    // idx_t i, j, k, nvtxs, nofsub, firstvtx;
    // idx_t *xadj, *adjncy, *label;
    // idx_t *perm, *iperm, *head, *qsize, *list, *marker;

    // WCOREPUSH;
    let graph = graph.as_mut().unwrap();

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices_mut!(graph => xadj adjncy);
    get_graph_slices!(graph => label);

    let mut perm = vec![0; nvtxs + 5];
    let mut iperm = vec![0; nvtxs + 5];
    let mut head = vec![0; nvtxs + 5];
    let mut qsize = vec![0; nvtxs + 5];
    let mut list = vec![0; nvtxs + 5];
    let mut marker = vec![0; nvtxs + 5];

    let mut nofsub = 0;
    mmd::genmmd(
        nvtxs as idx_t,
        xadj.as_mut_ptr(),
        adjncy.as_mut_ptr(),
        iperm.as_mut_ptr(),
        perm.as_mut_ptr(),
        1,
        head.as_mut_ptr(),
        qsize.as_mut_ptr(),
        list.as_mut_ptr(),
        marker.as_mut_ptr(),
        idx_t::MAX,
        &mut nofsub,
    );

    // label = graph.label;
    let firstvtx = lastvtx - nvtxs as idx_t;
    for i in (0)..(nvtxs) {
        // can't use a slice since size is from the original graph
        // I am a bit suspicious of this though
        *order.add(label[i] as usize) = firstvtx + iperm[i] - 1;
    }

    // WCOREPOP;
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]

    use super::*;
    use crate::{dyncall::{ab_test_eq, ab_test_single_eq}, graph_gen::{Csr, GraphBuilder}, tests::{ab_test_partition_test_graphs, TestGraph}};

    #[test]
    fn ab_basic_METIS_NodeND() {
        ab_test_partition_test_graphs("METIS_NodeND:rs", Optype::Ometis, 3, 1, std::convert::identity);
    }

    #[test]
    fn ab_METIS_NodeND() {
        ab_test_partition_test_graphs("METIS_NodeND:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_METIS_NodeND_nocompress() {
        ab_test_partition_test_graphs("METIS_NodeND:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g.set_compress(false);
            g
        });
    }

    #[test]
    fn ab_METIS_NodeND_prune() {
        ab_test_partition_test_graphs("METIS_NodeND:rs", Optype::Ometis, 3, 1, |mut g| {
            g.set_pfactor(1);
            g.random_vwgt();
            g
        });
        ab_test_partition_test_graphs("METIS_NodeND:rs", Optype::Ometis, 3, 1, |mut g| {
            g.set_pfactor(5);
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_MlevelNestedDissect_norm() {
        ab_test_partition_test_graphs("MlevelNestedDissection:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_MlevelNestedDissectCC() {
        ab_test_partition_test_graphs("MlevelNestedDissectionCC:rs", Optype::Ometis, 3, 1, |mut g| {
            g.set_ccorder(true);
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_MlevelNodeBisectionMultiple() {
        ab_test_partition_test_graphs("MlevelNodeBisectionMultiple:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_MlevelNodeBisectionL2() {
        ab_test_partition_test_graphs("MlevelNodeBisectionL2:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_MlevelNodeBisectionL1() {
        ab_test_partition_test_graphs("MlevelNodeBisectionL1:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_SplitGraphOrder() {
        ab_test_partition_test_graphs("SplitGraphOrder:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_SplitGraphOrderCC() {
        ab_test_partition_test_graphs("SplitGraphOrderCC:rs", Optype::Ometis, 3, 1, |mut g| {
            g.set_ccorder(true);
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_MMDOrder() {
        ab_test_partition_test_graphs("MMDOrder:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn metis_issue_26() {
        // https://github.com/KarypisLab/METIS/issues/26
        // these are directly the values in the issue, but they contain self-loops
        // I'm not sure if the self-loops caused the issue reported, but I can't seem to reproduce
        // it even with ASAN
        #[rustfmt::skip]
        let base_xadj = [
            0, 8, 13, 21, 32, 37, 48, 55, 61, 62, 63, 67, 71, 83, 84, 91, 101, 102, 108, 113, 125,
            130, 135, 136, 141, 153, 158, 166, 170, 171, 178
        ];

        #[rustfmt::skip]
        let base_adjncy = [
            0, 3, 5, 10, 17, 26, 27, 29, 1, 4, 18, 20, 25, 2, 6, 7, 12, 14, 15, 19, 24, 0, 3, 5,
            11, 12, 15, 17, 19, 24, 26, 29, 1, 4, 18, 20, 25, 0, 3, 5, 11, 12, 15, 17, 19, 24, 26,
            29, 2, 6, 12, 14, 15, 19, 24, 2, 7, 12, 15, 19, 24, 8, 9, 0, 10, 26, 27, 3, 5, 11, 29,
            2, 3, 5, 6, 7, 12, 14, 15, 19, 21, 23, 24, 13, 2, 6, 12, 14, 15, 19, 24, 2, 3, 5, 6, 7,
            12, 14, 15, 19, 24, 16, 0, 3, 5, 17, 26, 29, 1, 4, 18, 20, 25, 2, 3, 5, 6, 7, 12, 14,
            15, 19, 21, 23, 24, 1, 4, 18, 20, 25, 12, 19, 21, 23, 24, 22, 12, 19, 21, 23, 24, 2, 3,
            5, 6, 7, 12, 14, 15, 19, 21, 23, 24, 1, 4, 18, 20, 25, 0, 3, 5, 10, 17, 26, 27, 29, 0,
            10, 26, 27, 28, 0, 3, 5, 11, 17, 26, 29
        ];

        let mut xadj = Vec::with_capacity(base_xadj.len());
        let mut adjncy = Vec::with_capacity(base_adjncy.len());
        for (i, w) in base_xadj.windows(2).enumerate() {
            xadj.push(adjncy.len() as idx_t);
            adjncy.extend(base_adjncy[w[0] as usize .. w[1] as usize].iter().copied().filter(|&e| e as usize != i));
        }
        xadj.push(adjncy.len() as idx_t);

        let graph = GraphBuilder::from_csr(Csr::from_parts(xadj, adjncy), Optype::Ometis, 3, 1);
        // graph.enable_dbg(DbgLvl::Info);
        ab_test_single_eq("CompressGraph:rs", || {
            graph.clone().call().unwrap()
        });
    }
}
