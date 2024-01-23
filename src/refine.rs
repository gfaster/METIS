/*
\file
\brief This file contains the driving routines for multilevel refinement

\date   Started 7/24/1997
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version\verbatim $Id: refine.c 14362 2013-05-21 21:35:23Z karypis $ \endverbatim
*/

use crate::*;

/*************************************************************************/
/* This function is the entry point of refinement */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Refine2Way(
    ctrl: *mut ctrl_t,
    orggraph: *mut graph_t,
    graph: *mut graph_t,
    tpwgts: *mut real_t,
) {
    let ctrl = ctrl.as_mut().unwrap();
    let mut graph = graph;

    // ifset!(
    //     ctrl.dbglvl,
    //     METIS_DBG_TIME,
    //     gk_startcputimer(ctrl.UncoarsenTmr)
    // );

    /* Compute the parameters of the coarsest graph */
    Compute2WayPartitionParams(ctrl, graph);

    loop {
        assert!(debug::CheckBnd(graph) != 0);

        // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.RefTmr));

        balance::Balance2Way(ctrl, graph, tpwgts);

        FM_2WayRefine(ctrl, graph, tpwgts, ctrl.niter);

        // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.RefTmr));

        if graph == orggraph {
            break;
        }

        graph = (*graph).finer;
        graph::graph_ReadFromDisk(ctrl, graph);

        // ifset!(
        //     ctrl.dbglvl,
        //     METIS_DBG_TIME,
        //     gk_startcputimer(ctrl.ProjectTmr)
        // );
        Project2WayPartition(ctrl, graph);
        // ifset!(
        //     ctrl.dbglvl,
        //     METIS_DBG_TIME,
        //     gk_stopcputimer(ctrl.ProjectTmr)
        // );
    }

    // ifset!(
    //     ctrl.dbglvl,
    //     METIS_DBG_TIME,
    //     gk_stopcputimer(ctrl.UncoarsenTmr)
    // );
}

/*************************************************************************/
/* This function allocates memory for 2-way edge refinement */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Allocate2WayPartitionMemory(_ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    // idx_t nvtxs, ncon;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    graph.pwgts = imalloc(2 * ncon, "Allocate2WayPartitionMemory: pwgts\0".as_ptr()) as _;
    graph.where_ = imalloc(nvtxs, "Allocate2WayPartitionMemory: where_\0".as_ptr()) as _;
    graph.bndptr = imalloc(nvtxs, "Allocate2WayPartitionMemory: bndptr\0".as_ptr()) as _;
    graph.bndind = imalloc(nvtxs, "Allocate2WayPartitionMemory: bndind\0".as_ptr()) as _;
    graph.id = imalloc(nvtxs, "Allocate2WayPartitionMemory: id\0".as_ptr()) as _;
    graph.ed = imalloc(nvtxs, "Allocate2WayPartitionMemory: ed\0".as_ptr()) as _;
}

/*************************************************************************/
/* This function computes the initial id/ed */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Compute2WayPartitionParams(_ctrl: *mut ctrl_t, graph: *mut graph_t) {
    // idx_t i, j, nvtxs, ncon, nbnd, mincut, istart, iend, tid, ted, me;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *pwgts;
    // idx_t *where_, *bndptr, *bndind, *id, *ed;

    let graph = graph.as_mut().unwrap();
    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    // xadj = graph.xadj;
    // vwgt = graph.vwgt;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    //
    // where_ = graph.where_;
    // id = graph.id;
    // ed = graph.ed;
    //
    //
    // pwgts = iset(2 * ncon, 0, graph.pwgts);
    // bndptr = iset(nvtxs, -1, graph.bndptr);
    // bndind = graph.bndind;
    get_graph_slices!(graph => xadj vwgt adjncy adjwgt tvwgt);
    get_graph_slices_mut!(graph => where_ id ed bndptr bndind);

    // use mkslice for pwgts since this is strictly bisectioning, and the allocation allocates 2 *
    // ncon (see above)
    mkslice_mut!(graph->pwgts, 2 * ncon);

    bndptr.fill(-1);
    pwgts.fill(0);

    /* Compute pwgts */
    if ncon == 1 {
        for i in (0)..(nvtxs) {
            assert!(where_[i as usize] >= 0 && where_[i as usize] <= 1);
            pwgts[where_[i as usize] as usize] += vwgt[i as usize];
        }
        assert!(pwgts[0 as usize] + pwgts[1 as usize] == tvwgt[0 as usize]);
    } else {
        for i in (0)..(nvtxs) {
            let me = where_[i as usize] as usize;
            for j in (0)..(ncon) {
                pwgts[(me * ncon + j) as usize] += vwgt[(i * ncon + j) as usize];
            }
        }
    }

    /* Compute the required info for refinement  */
    let mut nbnd = 0;
    let mut mincut = 0;
    for i in (0)..(nvtxs) {
        let istart = xadj[i as usize] as usize;
        let iend = xadj[(i + 1) as usize] as usize;

        let me = where_[i as usize] as usize;
        let mut tid = 0;
        let mut ted = 0;

        for j in (istart)..(iend) {
            if me == where_[adjncy[j as usize] as usize] as usize {
                tid += adjwgt[j as usize] as usize;
            } else {
                ted += adjwgt[j as usize] as usize;
            }
        }
        id[i as usize] = tid as idx_t;
        ed[i as usize] = ted as idx_t;

        if ted > 0 || istart == iend {
            BNDInsert!(nbnd, bndind, bndptr, i);
            mincut += ted;
        }
    }

    graph.mincut = mincut as idx_t / 2;
    graph.nbnd = nbnd;
}

/*************************************************************************/
/* Projects a partition and computes the refinement params. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Project2WayPartition(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    // idx_t i, j, istart, iend, nvtxs, nbnd, me, tid, ted;
    // idx_t *xadj, *adjncy, *adjwgt;
    // idx_t *cmap, *where_, *bndptr, *bndind;
    // idx_t *cwhere_, *cbndptr;
    // idx_t *id, *ed;
    // graph_t *cgraph;
    // int dropedges;

    Allocate2WayPartitionMemory(ctrl, graph);
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    let dropedges = ctrl.dropedges;

    assert!(!graph.coarser.is_null());
    let cgraph = &*graph.coarser;
    mkslice!(cwhere_: cgraph->where_, cgraph.nvtxs);
    mkslice!(cbndptr: cgraph->bndptr, cgraph.nvtxs);

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj adjncy adjwgt);
    get_graph_slices_mut!(ctrl, graph => cmap id ed bndptr bndind where_ pwgts);

    // bndptr = iset(nvtxs, -1, graph.bndptr);
    // bndind = graph.bndind;
    bndptr.fill(-1);

    /* Project the partition and record which of these nodes came from the
    coarser boundary */
    for i in (0)..(nvtxs) {
        let j = cmap[i as usize];
        where_[i as usize] = cwhere_[j as usize];
        cmap[i as usize] = if dropedges != 0 {
            0
        } else {
            cbndptr[j as usize]
        };
    }

    /* Compute the refinement information of the nodes */
    let mut nbnd = 0;
    for i in (0)..(nvtxs) {
        let istart = xadj[i as usize] as usize;
        let iend = xadj[(i + 1) as usize] as usize;

        let mut tid = 0;
        let mut ted = 0;
        if cmap[i as usize] == -1 {
            /* Interior node. Note that cmap[i as usize] = cbndptr[cmap[i as usize] as usize] */
            for j in (istart)..(iend) {
                tid += adjwgt[j as usize];
            }
        } else {
            /* Potentially an interface node */
            let me = where_[i as usize];
            for j in (istart)..(iend) {
                if me == where_[adjncy[j as usize] as usize] {
                    tid += adjwgt[j as usize];
                } else {
                    ted += adjwgt[j as usize];
                }
            }
        }
        id[i as usize] = tid;
        ed[i as usize] = ted;

        if ted > 0 || istart == iend {
            BNDInsert!(nbnd, bndind, bndptr, i);
        }
    }
    graph.mincut = if dropedges != 0 {
        ComputeCut(graph, where_.as_ptr())
    } else {
        cgraph.mincut
    };
    graph.nbnd = nbnd;

    /* copy pwgts */
    // icopy(2 * graph.ncon, cgraph.pwgts, graph.pwgts);
    mkslice!(cpwgts: cgraph->pwgts, 2 * graph.ncon);
    pwgts[..(2 * graph.ncon as usize)].copy_from_slice(cpwgts);

    graph::FreeGraph(&mut graph.coarser);
    graph.coarser = std::ptr::null_mut();
}
