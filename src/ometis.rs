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
           the original and permuted matrices, then A'[i as usize] = A[perm[i as usize] as usize].
    \param iperm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A[i as usize] = A'[iperm[i as usize] as usize].
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn METIS_NodeND(
    nvtxs: *mut idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    options: *mut idx_t,
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
    if (!gk_malloc_init()) {
        return METIS_ERROR_MEMORY;
    }

    //   gk_sigtrap();
    //
    //   if ((sigrval = gk_sigcatch()) != 0)  {
    //     goto SIGTHROW;
    // }

    /* set up the run time parameters */
    ctrl = SetupCtrl(
        METIS_OP_OMETIS,
        options,
        1,
        3,
        std::ptr::null_mut(),
        std::ptr::null_mut(),
    );
    if (!ctrl) {
        gk_siguntrap();
        return METIS_ERROR_INPUT;
    }

    /* if required, change the numbering to 0 */
    // if (ctrl.numflag == 1) {
    //     Change2CNumbering(*nvtxs, xadj, adjncy);
    //     renumber = 1;
    // }

    ifset!(ctrl.dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
    ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

    /* prune the dense columns */
    if (ctrl.pfactor > 0.0) {
        piperm = imalloc(*nvtxs, "OMETIS: piperm");

        graph = PruneGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, piperm, ctrl.pfactor);
        if (graph == std::ptr::null_mut()) {
            /* if there was no prunning, cleanup the pfactor */
            gk_free_one(&mut piperm, LTERM);
            ctrl.pfactor = 0.0;
        } else {
            nnvtxs = graph.nvtxs;
            ctrl.compress = 0; /* disable compression if prunning took place */
        }
    }

    /* compress the graph; note that compression only happens if not prunning
    has taken place. */
    if (ctrl.compress) {
        cptr = imalloc(*nvtxs + 1, "OMETIS: cptr");
        cind = imalloc(*nvtxs, "OMETIS: cind");

        graph = CompressGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, cptr, cind);
        if (graph == std::ptr::null_mut()) {
            /* if there was no compression, cleanup the compress flag */
            gk_free_one(&mut cptr);
            gk_free_one(&mut cind);
            ctrl.compress = 0;
        } else {
            nnvtxs = graph.nvtxs;
            ctrl.cfactor = 1.0 * (*nvtxs) / nnvtxs;
            if (ctrl.cfactor > 1.5 && ctrl.nseps == 1) {
                ctrl.nseps = 2;
            }
            //ctrl.nseps = (idx_t)(ctrl.cfactor*ctrl.nseps);
        }
    }

    /* if no prunning and no compression, setup the graph in the normal way. */
    if (ctrl.pfactor == 0.0 && ctrl.compress == 0) {
        graph = SetupGraph(
            ctrl,
            *nvtxs,
            1,
            xadj,
            adjncy,
            vwgt,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );
    }

    assert!(CheckGraph(graph, ctrl.numflag, 1));

    /* allocate workspace memory */
    AllocateWorkSpace(ctrl, graph);

    /* do the nested dissection ordering  */
    if (ctrl.ccorder) {
        MlevelNestedDissectionCC(ctrl, graph, iperm, graph.nvtxs);
    } else {
        MlevelNestedDissection(ctrl, graph, iperm, graph.nvtxs);
    }

    if (ctrl.pfactor > 0.0) {
        /* Order any prunned vertices */
        icopy(nnvtxs, iperm, perm); /* Use perm as an auxiliary array */
        for i in (0)..(nnvtxs) {
            iperm[piperm[i as usize] as usize] = perm[i as usize];
        }
        for i in (nnvtxs)..(*nvtxs) {
            iperm[piperm[i as usize] as usize] = i;
        }

        gk_free_one(&mut piperm);
    } else if (ctrl.compress) {
        /* Uncompress the ordering */
        /* construct perm from iperm */
        for i in (0)..(nnvtxs) {
            perm[iperm[i as usize] as usize] = i;
        }
        let mut l = 0;
        for ii in 0..nnvtxs {
            i = perm[ii as usize];
            for j in (cptr[i as usize])..(cptr[(i + 1) as usize]) {
                iperm[cind[j as usize] as usize] = l += 1;
            }
        }

        gk_free_one(&mut cptr);
        gk_free_one(&mut cind);
    }

    for i in (0)..(*nvtxs) {
        perm[iperm[i as usize] as usize] = i;
    }

    ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
    ifset!(ctrl.dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

    /* clean up */
    FreeCtrl(&ctrl);

    // SIGTHROW:
    /* if required, change the numbering back to 1 */
    // I'm not doing fortran renaming
    //   if (renumber) {
    //     Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
    // }

    // gk_siguntrap();
    gk_malloc_cleanup(0);

    return util::metis_rcode(sigrval);
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
    // idx_t i, j, nvtxs, nbnd;
    // idx_t *label, *bndind;
    // graph_t *lgraph, *rgraph;

    let nvtxs = graph.nvtxs as usize;

    MlevelNodeBisectionMultiple(ctrl, graph);

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_SEPINFO,
        print!(
            "Nvtxs: {:6}, [({:6} {:6} {:6}) as usize]\n",
            graph.nvtxs, graph.pwgts[0 as usize], graph.pwgts[1 as usize], graph.pwgts[2 as usize]
        )
    );

    /* Order the nodes in the separator */
    nbnd = graph.nbnd;
    bndind = graph.bndind;
    label = graph.label;
    for i in (0)..(nbnd) {
        lastvtx -= 1;
        order[label[bndind[i as usize] as usize]] = lastvtx;
    }

    SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

    /* Free the memory of the top level graph */
    FreeGraph(&graph);

    /* Recurse on lgraph first, as its lastvtx depends on rgraph.nvtxs, which
    will not be defined upon return from MlevelNestedDissection. */
    if (lgraph.nvtxs > MMDSWITCH && lgraph.nedges > 0) {
        MlevelNestedDissection(ctrl, lgraph, order, lastvtx - rgraph.nvtxs);
    } else {
        MMDOrder(ctrl, lgraph, order, lastvtx - rgraph.nvtxs);
        FreeGraph(&lgraph);
    }
    if (rgraph.nvtxs > MMDSWITCH && rgraph.nedges > 0) {
        MlevelNestedDissection(ctrl, rgraph, order, lastvtx);
    } else {
        MMDOrder(ctrl, rgraph, order, lastvtx);
        FreeGraph(&rgraph);
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
    // idx_t i, j, nvtxs, nbnd, ncmps, rnvtxs, snvtxs;
    // idx_t *label, *bndind;
    // idx_t *cptr, *cind;
    // graph_t **sgraphs;

    let nvtxs = graph.nvtxs as usize;

    MlevelNodeBisectionMultiple(ctrl, graph);

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_SEPINFO,
        print!(
            "Nvtxs: {:6}, [({:6} {:6} {:6}) as usize]\n",
            graph.nvtxs, graph.pwgts[0 as usize], graph.pwgts[1 as usize], graph.pwgts[2 as usize]
        )
    );

    /* Order the nodes in the separator */
    nbnd = graph.nbnd;
    bndind = graph.bndind;
    label = graph.label;
    for i in (0)..(nbnd) {
        lastvtx -= 1;
        order[label[bndind[i as usize] as usize]] = lastvtx;
    }

    // WCOREPUSH;
    cptr = vec![0; nvtxs + 1 as usize];
    cind = vec![0; nvtxs as usize];
    ncmps = FindSepInducedComponents(ctrl, graph, cptr, cind);

    if (ctrl.dbglvl & METIS_DBG_INFO) {
        if (ncmps > 2) {
            print!("  Bisection resulted in {:} connected components\n", ncmps);
        }
    }

    sgraphs = SplitGraphOrderCC(ctrl, graph, ncmps, cptr, cind);

    // WCOREPOP;

    /* Free the memory of the top level graph */
    FreeGraph(&graph);

    /* Go and process the subgraphs */
    rnvtxs = 0;
    for i in 0..ncmps {
        /* Save the number of vertices in sgraphs[i as usize] because sgraphs[i as usize] is freed
        inside MlevelNestedDissectionCC, and as such it will be undefined. */
        snvtxs = sgraphs[i as usize].nvtxs;

        if (sgraphs[i as usize].nvtxs > MMDSWITCH && sgraphs[i as usize].nedges > 0) {
            MlevelNestedDissectionCC(ctrl, sgraphs[i as usize], order, lastvtx - rnvtxs);
        } else {
            MMDOrder(ctrl, sgraphs[i as usize], order, lastvtx - rnvtxs);
            FreeGraph(&sgraphs[i as usize]);
        }
        rnvtxs += snvtxs;
    }

    gk_free_one(&mut sgraphs, LTERM);
}

/*************************************************************************/
/* This function performs multilevel node bisection (i.e., tri-section).
It performs multiple bisections and selects the best. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNodeBisectionMultiple(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    // idx_t i, mincut;
    // idx_t *bestwhere_;

    /* if the graph is small, just find a single vertex separator */
    if (ctrl.nseps == 1 || graph.nvtxs < (if ctrl.compress { 1000 } else { 2000 })) {
        MlevelNodeBisectionL2(ctrl, graph, LARGENIPARTS);
        return;
    }

    // WCOREPUSH;

    bestwhere_ = vec![0; graph.nvtxs as usize];

    mincut = graph.tvwgt[0 as usize];
    for i in (0)..(ctrl.nseps) {
        MlevelNodeBisectionL2(ctrl, graph, LARGENIPARTS);

        if (i == 0 || graph.mincut < mincut) {
            mincut = graph.mincut;
            if (i < ctrl.nseps - 1) {
                icopy(graph.nvtxs, graph.where_, bestwhere_);
            }
        }

        if (mincut == 0) {
            break;
        }

        if (i < ctrl.nseps - 1) {
            FreeRData(graph);
        }
    }

    if (mincut != graph.mincut) {
        icopy(graph.nvtxs, bestwhere_, graph.where_);
        Compute2WayNodePartitionParams(ctrl, graph);
    }

    // WCOREPOP;
}

/*************************************************************************/
/* This function performs multilevel node bisection (i.e., tri-section).
It performs multiple bisections and selects the best. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNodeBisectionL2(ctrl: *mut ctrl_t, graph: *mut graph_t, niparts: idx_t) {
    // idx_t i, mincut, nruns=5;
    // graph_t *cgraph;
    // idx_t *bestwhere_;

    /* if the graph is small, just find a single vertex separator */
    if (graph.nvtxs < 5000) {
        MlevelNodeBisectionL1(ctrl, graph, niparts);
        return;
    }

    // WCOREPUSH;

    ctrl.CoarsenTo = gk_max(100, graph.nvtxs / 30);

    cgraph = CoarsenGraphNlevels(ctrl, graph, 4);

    bestwhere_ = vec![0; cgraph.nvtxs as usize];

    mincut = graph.tvwgt[0 as usize];
    for i in (0)..(nruns) {
        MlevelNodeBisectionL1(ctrl, cgraph, 0.7 * niparts);

        if (i == 0 || cgraph.mincut < mincut) {
            mincut = cgraph.mincut;
            if (i < nruns - 1) {
                icopy(cgraph.nvtxs, cgraph.where_, bestwhere_);
            }
        }

        if (mincut == 0) {
            break;
        }

        if (i < nruns - 1) {
            FreeRData(cgraph);
        }
    }

    if (mincut != cgraph.mincut) {
        icopy(cgraph.nvtxs, bestwhere_, cgraph.where_);
    }

    // WCOREPOP;

    Refine2WayNode(ctrl, graph, cgraph);
}

/*************************************************************************/
/* The top-level routine of the actual multilevel node bisection */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelNodeBisectionL1(ctrl: *mut ctrl_t, graph: *mut graph_t, niparts: idx_t) {
    // graph_t *cgraph;

    ctrl.CoarsenTo = graph.nvtxs / 8;
    if (ctrl.CoarsenTo > 100) {
        ctrl.CoarsenTo = 100;
    } else if (ctrl.CoarsenTo < 40) {
        ctrl.CoarsenTo = 40;
    }

    cgraph = CoarsenGraph(ctrl, graph);

    niparts = gk_max(
        1,
        (if cgraph.nvtxs <= ctrl.CoarsenTo {
            niparts / 2
        } else {
            niparts
        }),
    );
    /*niparts = (cgraph.nvtxs <= ctrl.CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);*/
    InitSeparator(ctrl, cgraph, niparts);

    Refine2WayNode(ctrl, graph, cgraph);
}

/*************************************************************************/
/* This function takes a graph and a tri-section (left, right, separator)
    and splits it into two graphs.

    This function relies on the fact that adjwgt is all equal to 1.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn SplitGraphOrder(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    r_lgraph: *mut *mut graph_t,
    r_rgraph: *mut *mut graph_t,
) {
    // idx_t i, ii, j, k, l, istart, iend, mypart, nvtxs, snvtxs[3 as usize], snedges[3 as usize];
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label, *where_, *bndptr, *bndind;
    // idx_t *sxadj[2 as usize], *svwgt[2 as usize], *sadjncy[2 as usize], *sadjwgt[2 as usize], *slabel[2 as usize];
    // idx_t *rename;
    // idx_t *auxadjncy;
    // graph_t *lgraph, *rgraph;

    // WCOREPUSH;

    ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.SplitTmr));

    let nvtxs = graph.nvtxs as usize;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;
    label = graph.label;
    where_ = graph.where_;
    bndptr = graph.bndptr;
    bndind = graph.bndind;
    assert!(bndptr != std::ptr::null_mut());

    rename = vec![0; nvtxs as usize];

    snvtxs[0 as usize] = 0;
    snvtxs[1 as usize] = 0;
    snvtxs[2 as usize] = 0;
    snedges[0 as usize] = 0;
    snedges[1 as usize] = 0;
    snedges[2 as usize] = 0;
    for i in (0)..(nvtxs) {
        k = where_[i as usize];
        rename[i as usize] = snvtxs[k as usize] += 1;
        snedges[k as usize] += xadj[(i + 1) as usize] - xadj[i as usize];
    }

    lgraph = SetupSplitGraph(graph, snvtxs[0 as usize], snedges[0 as usize]);
    sxadj[0 as usize] = lgraph.xadj;
    svwgt[0 as usize] = lgraph.vwgt;
    sadjncy[0 as usize] = lgraph.adjncy;
    sadjwgt[0 as usize] = lgraph.adjwgt;
    slabel[0 as usize] = lgraph.label;

    rgraph = SetupSplitGraph(graph, snvtxs[1 as usize], snedges[1 as usize]);
    sxadj[1 as usize] = rgraph.xadj;
    svwgt[1 as usize] = rgraph.vwgt;
    sadjncy[1 as usize] = rgraph.adjncy;
    sadjwgt[1 as usize] = rgraph.adjwgt;
    slabel[1 as usize] = rgraph.label;

    /* Go and use bndptr to also mark the boundary nodes in the two partitions */
    for ii in (0)..(graph.nbnd) {
        i = bndind[ii as usize];
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            bndptr[adjncy[j as usize] as usize] = 1;
        }
    }

    snvtxs[0 as usize] = 0;
    snvtxs[1 as usize] = 0;
    snedges[0 as usize] = 0;
    snedges[1 as usize] = 0;
    sxadj[0 as usize][0 as usize] = 0;
    sxadj[1 as usize][0 as usize] = 0;
    for i in (0)..(nvtxs) {
        mypart = where_[i as usize];
        if (mypart == 2) {
            continue;
        }

        istart = xadj[i as usize];
        iend = xadj[(i + 1) as usize];
        if (bndptr[i as usize] == -1) {
            /* This is an interior vertex */
            auxadjncy = sadjncy[mypart as usize] + snedges[mypart as usize] - istart;
            for j in istart..iend {
                auxadjncy[j as usize] = adjncy[j as usize];
            }
            snedges[mypart as usize] += iend - istart;
        } else {
            auxadjncy = sadjncy[mypart as usize];
            l = snedges[mypart as usize];
            for j in (istart)..(iend) {
                k = adjncy[j as usize];
                if (where_[k as usize] == mypart) {
                    auxadjncy[(l) as usize] = k;
                    l += 1;
                }
            }
            snedges[mypart as usize] = l;
        }

        svwgt[mypart as usize][snvtxs[mypart as usize] as usize] = vwgt[i as usize];
        slabel[mypart as usize][snvtxs[mypart as usize] as usize] = label[i as usize];
        snvtxs += 1;
        sxadj[mypart as usize][snvtxs[mypart as usize] as usize] = snedges[mypart as usize];
    }

    for mypart in (0)..(2) {
        iend = snedges[mypart as usize];
        iset(iend, 1, sadjwgt[mypart as usize]);

        auxadjncy = sadjncy[mypart as usize];
        for i in (0)..(iend) {
            auxadjncy[i as usize] = rename[auxadjncy[i as usize] as usize];
        }
    }

    lgraph.nvtxs = snvtxs[0 as usize];
    lgraph.nedges = snedges[0 as usize];
    rgraph.nvtxs = snvtxs[1 as usize];
    rgraph.nedges = snedges[1 as usize];

    SetupGraph_tvwgt(lgraph);
    SetupGraph_tvwgt(rgraph);

    ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.SplitTmr));

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
    ctrl: *mut ctrl_t,
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

    ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.SplitTmr));

    let nvtxs = graph.nvtxs as usize;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;
    label = graph.label;
    where_ = graph.where_;
    bndptr = graph.bndptr;
    bndind = graph.bndind;
    assert!(bndptr != std::ptr::null_mut());

    /* Go and use bndptr to also mark the boundary nodes in the two partitions */
    for ii in (0)..(graph.nbnd) {
        i = bndind[ii as usize];
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            bndptr[adjncy[j as usize] as usize] = 1;
        }
    }

    let mut rename = vec![0; nvtxs as usize];

    let sgraphs: *mut *mut graph_t = gk_malloc(
        size_of::<*mut graph_t>() * ncmps,
        "SplitGraphOrderCC: sgraphs",
    );

    /* Go and split the graph a component at a time */
    for iii in (0)..(ncmps) {
        irandArrayPermute(
            cptr[(iii + 1) as usize] - cptr[iii as usize],
            cind + cptr[iii as usize],
            cptr[(iii + 1) as usize] - cptr[iii as usize],
            0,
        );
        snvtxs = snedges = 0;
        for j in (cptr[iii as usize])..(cptr[(iii + 1) as usize]) {
            i = cind[j as usize];
            rename[i as usize] = snvtxs += 1;
            snedges += xadj[(i + 1) as usize] - xadj[i as usize];
        }

        sgraphs[iii as usize] = SetupSplitGraph(graph, snvtxs, snedges);

        sxadj = sgraphs[iii as usize].xadj;
        svwgt = sgraphs[iii as usize].vwgt;
        sadjncy = sgraphs[iii as usize].adjncy;
        sadjwgt = sgraphs[iii as usize].adjwgt;
        slabel = sgraphs[iii as usize].label;

        snvtxs = snedges = sxadj[0 as usize] = 0;
        for ii in (cptr[iii as usize])..(cptr[(iii + 1) as usize]) {
            i = cind[ii as usize];

            istart = xadj[i as usize];
            iend = xadj[(i + 1) as usize];
            if (bndptr[i as usize] == -1) {
                /* This is an interior vertex */
                auxadjncy = sadjncy + snedges - istart;
                for j in istart..iend {
                    auxadjncy[j as usize] = adjncy[j as usize];
                }
                snedges += iend - istart;
            } else {
                l = snedges;
                for j in (istart)..(iend) {
                    k = adjncy[j as usize];
                    if (where_[k as usize] != 2) {
                        sadjncy[(l) as usize] = k;
                        l += 1;
                    }
                }
                snedges = l;
            }

            svwgt[snvtxs as usize] = vwgt[i as usize];
            slabel[snvtxs as usize] = label[i as usize];
            snvtxs += 1;
            sxadj[(snvtxs) as usize] = snedges;
        }

        iset(snedges, 1, sadjwgt);
        for i in (0)..(snedges) {
            sadjncy[i as usize] = rename[sadjncy[i as usize] as usize];
        }

        sgraphs[iii as usize].nvtxs = snvtxs;
        sgraphs[iii as usize].nedges = snedges;

        SetupGraph_tvwgt(sgraphs[iii as usize]);
    }

    ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.SplitTmr));

    // WCOREPOP;

    return sgraphs;
}

/*************************************************************************/
/* This function uses MMD to order the graph. The vertices are numbered
from lastvtx downwards. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MMDOrder(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    order: *mut idx_t,
    lastvtx: idx_t,
) {
    // idx_t i, j, k, nvtxs, nofsub, firstvtx;
    // idx_t *xadj, *adjncy, *label;
    // idx_t *perm, *iperm, *head, *qsize, *list, *marker;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    xadj = graph.xadj;
    adjncy = graph.adjncy;

    /* Relabel the vertices so that it starts from 1 */
    k = xadj[nvtxs as usize];
    for i in (0)..(k) {
        adjncy[i as usize] += 1;
    }
    for i in (0)..(nvtxs + 1) {
        xadj[i as usize] += 1;
    }

    perm = vec![0; nvtxs + 5 as usize];
    iperm = vec![0; nvtxs + 5 as usize];
    head = vec![0; nvtxs + 5 as usize];
    qsize = vec![0; nvtxs + 5 as usize];
    list = vec![0; nvtxs + 5 as usize];
    marker = vec![0; nvtxs + 5 as usize];

    genmmd(
        nvtxs, xadj, adjncy, iperm, perm, 1, head, qsize, list, marker, IDX_MAX, &nofsub,
    );

    label = graph.label;
    firstvtx = lastvtx - nvtxs;
    for i in (0)..(nvtxs) {
        order[label[i as usize] as usize] = firstvtx + iperm[i as usize] - 1;
    }

    /* Relabel the vertices so that it starts from 0 */
    for i in (0)..(nvtxs + 1) {
        xadj[i as usize] -= 1;
    }
    k = xadj[nvtxs as usize];
    for i in (0)..(k) {
        adjncy[i as usize] -= 1;
    }

    // WCOREPOP;
}
