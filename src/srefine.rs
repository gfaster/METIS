/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * srefine.c
 *
 * This file contains code for the separator refinement algorithms
 *
 * Started 8/1/97
 * George
 *
 * $Id: srefine.c 14362 2013-05-21 21:35:23Z karypis $
 *
 */

use crate::*;

/*************************************************************************/
/* This function is the entry point of the separator refinement.
It does not perform any refinement on graph, but it starts by first
projecting it to the next level finer graph and proceeds from there. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Refine2WayNode(ctrl: *mut ctrl_t, orggraph: *mut graph_t, graph: *mut graph_t) {
    let ctrl = ctrl.as_mut().unwrap();

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.UncoarsenTmr));

    if graph == orggraph {
        Compute2WayNodePartitionParams(ctrl, graph);
    } else {
        let mut first = true;
        let mut graph = graph;
        while first || graph != orggraph {
            first = false;
            graph = (*graph).finer;

            graph::graph_ReadFromDisk(ctrl, graph);

            // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.ProjectTmr));
            Project2WayNodePartition(ctrl, graph);
            // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.ProjectTmr));

            // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.RefTmr));
            sfm::FM_2WayNodeBalance(ctrl, graph);

            debug_assert!(debug::CheckNodePartitionParams(graph) != 0);

            match ctrl.rtype {
                METIS_RTYPE_SEP2SIDED => sfm::FM_2WayNodeRefine2Sided(ctrl, graph, ctrl.niter),
                METIS_RTYPE_SEP1SIDED => sfm::FM_2WayNodeRefine1Sided(ctrl, graph, ctrl.niter),
                _ => panic!("Unknown rtype of {}\n", ctrl.rtype),
            };
            // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.RefTmr));
        }
    }

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.UncoarsenTmr));
}

/*************************************************************************/
/* This function allocates memory for 2-way node-based refinement */
/**************************************************************************/
#[metis_func]
pub extern "C" fn Allocate2WayNodePartitionMemory(_ctrl: *mut ctrl_t, graph: *mut graph_t) {
    // idx_t nvtxs;
    let graph = graph.as_mut().unwrap();
    let nvtxs = graph.nvtxs as usize;
    graph.pwgts = gk::imalloc(3, c"Allocate2WayNodePartitionMemory: pwgts").as_buf_ptr();
    graph.where_ = gk::imalloc(nvtxs, c"Allocate2WayNodePartitionMemory: where_").as_buf_ptr();
    graph.bndptr = gk::imalloc(nvtxs, c"Allocate2WayNodePartitionMemory: bndptr").as_buf_ptr();
    graph.bndind = gk::imalloc(nvtxs, c"Allocate2WayNodePartitionMemory: bndind").as_buf_ptr();
    graph.nrinfo = gk_malloc(
        nvtxs * size_of::<nrinfo_t>(),
        c"Allocate2WayNodePartitionMemory: nrinfo".as_ptr(),
    )
    .cast();
}

/*************************************************************************/
/* This function computes the edegrees[] to the left & right sides */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Compute2WayNodePartitionParams(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, nvtxs, nbnd;
    // idx_t *xadj, *adjncy, *vwgt;
    // idx_t *where_, *pwgts, *bndind, *bndptr, *edegrees;
    // nrinfo_t *rinfo;
    // idx_t me, other;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj vwgt adjncy);

    get_graph_slices_mut!(ctrl, graph => where_ nrinfo pwgts bndind bndptr);
    let rinfo = nrinfo;
    debug_assert_eq!(pwgts.len(), 3);
    pwgts.fill(0);
    bndptr.fill(-1);

    /*------------------------------------------------------------
    / Compute now the separator external degrees
    /------------------------------------------------------------*/
    let mut nbnd: usize = 0;
    for i in (0)..(nvtxs) {
        let me = where_[i as usize];
        debug_assert!(me >= 0 && me <= 2);

        pwgts[me as usize] += vwgt[i as usize];

        if me == 2 {
            /* If it is on the separator do some computations */
            BNDInsert!(nbnd, bndind, bndptr, i);

            rinfo[i].edegrees = [0, 0];

            for j in (xadj[i])..(xadj[i + 1]) {
                let other = where_[adjncy[j as usize] as usize] as usize;
                if other != 2 {
                    rinfo[i].edegrees[other as usize] += vwgt[adjncy[j as usize] as usize];
                }
            }
        }
    }

    debug_assert!(debug::CheckNodeBnd(graph, nbnd as idx_t) != 0);

    graph.mincut = pwgts[2 as usize];
    graph.nbnd = nbnd as idx_t;
}

/*************************************************************************/
/* This function projects the node-based bisection */
/*************************************************************************/
#[metis_func]
pub extern "C" fn Project2WayNodePartition(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    // idx_t i, j, nvtxs;
    // idx_t *cmap, *where_, *cwhere_;
    // graph_t *cgraph;
    let graph = graph.as_mut().unwrap();
    let cgraph = graph.coarser;
    // let cwhere_ = (*cgraph).where_;
    let cwhere_ = get_graph_slice!(*cgraph => where_);

    let nvtxs = graph.nvtxs as usize;

    Allocate2WayNodePartitionMemory(ctrl, graph);
    get_graph_slices_mut!(graph => where_);
    get_graph_slices!(graph => cmap);

    /* Project the partition */
    for i in (0)..(nvtxs) {
        where_[i as usize] = cwhere_[cmap[i as usize] as usize];
        debug_assert!(
            where_[i as usize] >= 0 && where_[i as usize] <= 2,
            "{:} {:} {:} {:}",
            i,
            cmap[i as usize],
            where_[i as usize],
            cwhere_[cmap[i as usize] as usize]
        );
    }

    graph::FreeGraph(&mut graph.coarser);
    graph.coarser = std::ptr::null_mut();

    Compute2WayNodePartitionParams(ctrl, graph);
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use super::*;
    use crate::tests::ab_test_partition_test_graphs;

    #[test]
    fn ab_Refine2WayNode() {
        ab_test_partition_test_graphs("Refine2WayNode:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_Compute2WayNodePartitionParams() {
        ab_test_partition_test_graphs(
            "Compute2WayNodePartitionParams:rs",
            Optype::Ometis,
            3,
            1,
            |mut g| {
                g.random_vwgt();
                g
            },
        );
    }

    #[test]
    fn ab_Project2WayNodePartition() {
        ab_test_partition_test_graphs(
            "Project2WayNodePartition:rs",
            Optype::Ometis,
            3,
            1,
            |mut g| {
                g.random_vwgt();
                g
            },
        );
    }
}
