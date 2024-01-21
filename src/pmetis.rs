/**
\file
\brief This file contains the top level routines for the multilevel recursive bisection
       algorithm PMETIS.

\date   Started 7/24/1997
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version\verbatim $Id: pmetis.c 10513 2011-07-07 22:06:03Z karypis $ \endverbatim
*/
use crate::{debug::CheckBnd, *};

/*************************************************************************/
/* \ingroup api
    \brief Recursive partitioning routine.

    This function computes a partitioning of a graph based on multilevel
    recursive bisection. It can be used to partition a graph into \e k
    parts. The objective of the partitioning is to minimize the edgecut
    subject to one or more balancing constraints.

    \param[in] nvtxs is the number of vertices in the graph.

    \param[in] ncon is the number of balancing constraints. For the standard
           partitioning problem in which each vertex is either unweighted
           or has a single weight, ncon should be 1.

    \param[in] xadj is an array of size nvtxs+1 used to specify the starting
           positions of the adjacency structure of the vertices in the
           adjncy array.

    \param[in] adjncy is an array of size to the sum of the degrees of the
           graph that stores for each vertex the set of vertices that
           is adjacent to.

    \param[in] vwgt is an array of size nvtxs*ncon that stores the weights
           of the vertices for each constraint. The ncon weights for the
           ith vertex are stored in the ncon consecutive locations starting
           at vwgt[i*ncon]. When ncon==1, a NULL value can be passed indicating
           that all the vertices in the graph have the same weight.

    \param[in] adjwgt is an array of size equal to adjncy, specifying the weight
           for each edge (i.e., adjwgt[j] corresponds to the weight of the
           edge stored in adjncy[j]).
           A NULL value can be passed indicating that all the edges in the
           graph have the same weight.

    \param[in] nparts is the number of desired partitions.

    \param[in] tpwgts is an array of size nparts*ncon that specifies the
           desired weight for each part and constraint. The \e{target partition
           weight} for the ith part and jth constraint is specified
           at tpwgts[i*ncon+j] (the numbering of i and j starts from 0).
           For each constraint, the sum of the tpwgts[] entries must be
           1.0 (i.e., \f$ \sum_i tpwgts[i*ncon+j] = 1.0 \f$).
           A NULL value can be passed indicating that the graph should
           be equally divided among the parts.

    \param[in] ubvec is an array of size ncon that specifies the allowed
           load imbalance tolerance for each constraint.
           For the ith part and jth constraint the allowed weight is the
           ubvec[j]*tpwgts[i*ncon+j] fraction of the jth's constraint total
           weight. The load imbalances must be greater than 1.0.
           A NULL value can be passed indicating that the load imbalance
           tolerance for each constraint should be 1.001 (for ncon==1)
           or 1.01 (for ncon>1).

    \params[in] options is the array for passing additional parameters
           in order to customize the behaviour of the partitioning
           algorithm.

    \params[out] edgecut stores the cut of the partitioning.

    \params[out] part is an array of size nvtxs used to store the
           computed partitioning. The partition number for the ith
           vertex is stored in part[i]. Based on the numflag parameter,
           the numbering of the parts starts from either 0 or 1.


    \returns
      \retval METIS_OK  indicates that the function returned normally.
      \retval METIS_ERROR_INPUT indicates an input error.
      \retval METIS_ERROR_MEMORY indicates that it could not allocate
              the required memory.

*/
/*************************************************************************/
#[allow(non_snake_case)]
pub unsafe extern "C" fn METIS_PartGraphRecursive(
    nvtxs: *mut idx_t,
    ncon: *mut idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    vsize: *mut idx_t,
    adjwgt: *mut idx_t,
    nparts: *mut idx_t,
    tpwgts: *mut real_t,
    ubvec: *const real_t,
    options: *mut idx_t,
    objval: *mut idx_t,
    part: *mut idx_t,
) -> std::ffi::c_int {
    // int sigrval=0, renumber=0;
    // graph_t *graph;
    // ctrl_t *ctrl;

    /* set up malloc cleaning code and signal catchers */
    if gk_malloc_init() == 0 {
        eprintln!("memory init error");
        return METIS_ERROR_MEMORY;
    }

    /* set up the run parameters */
    let ctrl = SetupCtrl(METIS_OP_PMETIS, options, *ncon, *nparts, tpwgts, ubvec);
    if ctrl.is_null() {
        eprintln!("input error");
        return METIS_ERROR_INPUT;
    }
    let ctrl = ctrl.as_mut().unwrap();

    /* if required, change the numbering to 0 */
    // no.
    // if (ctrl.numflag == 1) {
    //   Change2CNumbering(*nvtxs, xadj, adjncy);
    //   renumber = 1;
    // }

    /* set up the graph */
    let graph = graph::SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);
    let graph = graph.as_mut().unwrap();

    debug_assert!(debug::check_adj(graph));
    /* allocate workspace memory */
    AllocateWorkSpace(ctrl, graph);

    /* start the partitioning */
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

    // iset(*nvtxs, 0, part);
    mkslice_mut!(part, *nparts);
    part.fill(0);
    *objval = if *nparts == 1 {
        0
    } else {
        MlevelRecursiveBisection(ctrl, graph, *nparts, part.as_mut_ptr(), ctrl.tpwgts, 0)
    };

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

    /* clean up */
    FreeCtrl(&mut (ctrl as *mut ctrl_t) as *mut *mut ctrl_t);

    gk_malloc_cleanup(0);

    return util::metis_rcode(0);
}

/*************************************************************************/
/* This function is the top-level driver of the recursive bisection
routine. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MlevelRecursiveBisection(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    nparts: idx_t,
    part: *mut idx_t,
    tpwgts: *mut real_t,
    fpart: idx_t,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, nvtxs, ncon, objval;
    // idx_t *label, *where_;
    // graph_t *lgraph, *rgraph;
    // real_t wsum, *tpwgts2;

    let nparts = nparts as usize;
    let nvtxs = graph.nvtxs as usize;
    if nvtxs == 0 {
        print!(concat!(
            "\t***Cannot bisect a graph with 0 vertices!\n",
            "\t***You are trying to partition a graph into too many parts!\n"
        ));
        return 0;
    }

    let ncon = graph.ncon as usize;

    /* determine the weights of the two partitions as a function of the weight of the
    target partition weights */
    // WCOREPUSH;
    let mut tpwgts2 = vec![0.0; 2 * ncon as usize];
    mkslice_mut!(tpwgts, ncon * ctrl.nparts as usize);
    for i in (0)..(ncon) {
        // tpwgts2[i] = rsum((nparts >> 1), tpwgts + i, ncon);
        tpwgts2[i] = tpwgts[i..]
            .iter()
            .step_by(ncon)
            .take((nparts as usize) >> 1)
            .sum();
        tpwgts2[ncon + i] = 1.0 - tpwgts2[i];
    }

    /* perform the bisection */
    let mut objval = MultilevelBisect(ctrl, graph, tpwgts2.as_mut_ptr());

    // WCOREPOP;

    // not a slice because idt we know the len
    // mkslice_mut!(part, nvtxs);
    assert!(!part.is_null());

    // label = graph.label;
    // where_ = graph.where_;
    get_graph_slices_mut!(graph => label where_);
    for i in (0)..(nvtxs) {
        // part[label[i] as usize] = where_[i] + fpart;
        *part.add(label[i] as usize) = where_[i] + fpart;
    }

    let mut lgraph: *mut graph_t = std::ptr::null_mut();
    let mut rgraph: *mut graph_t = std::ptr::null_mut();

    if nparts > 2 {
        SplitGraphPart(
            ctrl,
            graph,
            &mut lgraph as *mut *mut graph_t,
            &mut rgraph as *mut *mut graph_t,
        );
    }

    /* Free the memory of the top level graph */
    graph::FreeGraph(&mut (graph as *mut graph_t) as *mut *mut graph_t);

    /* Scale the fractions in the tpwgts according to the true weight */
    for i in (0)..(ncon) {
        // wsum = rsum((nparts >> 1), tpwgts + i, ncon);
        let wsum: real_t = tpwgts[i..]
            .iter()
            .step_by(ncon)
            .take((nparts as usize) >> 1)
            .sum();
        // rscale((nparts >> 1), 1.0 / wsum, tpwgts + i, ncon);
        tpwgts[i..]
            .iter_mut()
            .step_by(ncon)
            .take((nparts as usize) >> 1)
            .map(|w| *w *= 1.0 / wsum)
            .last();
        // rscale(
        //     nparts - (nparts >> 1),
        //     1.0 / (1.0 - wsum),
        //     tpwgts + (nparts >> 1) * ncon + i,
        //     ncon,
        // );
        tpwgts[((nparts >> 1) * ncon + i)..]
            .iter_mut()
            .step_by(ncon)
            .take(nparts - (nparts >> 1))
            .map(|w| *w *= 1.0 / (1.0 - wsum))
            .last();
    }

    /* Do the recursive call */
    if nparts > 3 {
        objval += MlevelRecursiveBisection(
            ctrl,
            lgraph,
            (nparts >> 1) as idx_t,
            part,
            tpwgts.as_mut_ptr(),
            fpart,
        );
        objval += MlevelRecursiveBisection(
            ctrl,
            rgraph,
            (nparts - (nparts >> 1)) as idx_t,
            part,
            tpwgts[((nparts >> 1) * ncon)..].as_mut_ptr(),
            (fpart as usize + (nparts >> 1)) as idx_t,
        );
    } else if nparts == 3 {
        graph::FreeGraph(&mut lgraph as *mut *mut graph_t);
        objval += MlevelRecursiveBisection(
            ctrl,
            rgraph,
            (nparts - (nparts >> 1)) as idx_t,
            part,
            tpwgts[(nparts >> 1) * ncon..].as_mut_ptr(),
            (fpart as usize + (nparts >> 1)) as idx_t,
        );
    }

    return objval;
}

/*************************************************************************/
/* This function performs a multilevel bisection */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MultilevelBisect(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    tpwgts: *mut real_t,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, niparts, bestobj=0, curobj=0, *bestwhere_=NULL;
    // graph_t *cgraph;
    // real_t bestbal=0.0, curbal=0.0;

    Setup2WayBalMultipliers(ctrl, graph, tpwgts);

    // WCOREPUSH;

    let mut bestobj = 0;

    let mut bestwhere_: Vec<_>;
    if ctrl.ncuts > 1 {
        bestwhere_ = vec![0; graph.nvtxs as usize];
    } else {
        bestwhere_ = vec![];
    }

    let mut bestbal = 0.0;
    let mut curobj = 0;
    for i in (0)..(ctrl.ncuts) {
        let cgraph = coarsen::CoarsenGraph(ctrl, graph);
        let cgraph = cgraph.as_mut().unwrap();

        let niparts = if cgraph.nvtxs <= ctrl.CoarsenTo {
            SMALLNIPARTS
        } else {
            LARGENIPARTS
        };
        initpart::Init2WayPartition(ctrl, cgraph, tpwgts, niparts);

        refine::Refine2Way(ctrl, graph, cgraph, tpwgts);

        curobj = graph.mincut;
        let curbal = ComputeLoadImbalanceDiff(graph, 2, ctrl.pijbm, ctrl.ubfactors);

        if i == 0
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal)
        {
            bestobj = curobj;
            bestbal = curbal;
            if i < ctrl.ncuts - 1 {
                // icopy(graph.nvtxs, graph.where_, bestwhere_);
                get_graph_slices!(graph => where_);
                bestwhere_.copy_from_slice(where_);
            }
        }

        if bestobj == 0 {
            break;
        }

        if i < ctrl.ncuts - 1 {
            graph::FreeRData(graph);
        }
    }

    if bestobj != curobj {
        // icopy(graph.nvtxs, bestwhere_, graph.where_);
        get_graph_slices_mut!(graph => where_);
        where_.copy_from_slice(&bestwhere_);

        refine::Compute2WayPartitionParams(ctrl, graph);
    }

    // WCOREPOP;

    return bestobj;
}

/*************************************************************************/
/* This function splits a graph into two based on its bisection */
/*************************************************************************/
#[metis_func]
pub extern "C" fn SplitGraphPart(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    r_lgraph: *mut *mut graph_t,
    r_rgraph: *mut *mut graph_t,
) -> () {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, k, l, istart, iend, mypart, nvtxs, ncon, snvtxs[2], snedges[2];
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label, *where_, *bndptr;
    // idx_t *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *slabel[2];
    // idx_t *rename;
    // idx_t *auxadjncy, *auxadjwgt;
    // graph_t *lgraph, *rgraph;

    // WCOREPUSH;

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.SplitTmr));

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon;

    assert!(!graph.bndptr.is_null());
    debug_assert!(debug::check_adj(graph));

    get_graph_slices!(graph => xadj vwgt adjncy adjwgt where_ bndptr label);

    debug_assert!(debug::CheckBnd(graph) != 0);

    let mut rename = vec![0; nvtxs as usize];

    let mut snvtxs = [0, 0];
    let mut snedges = [0, 0];
    for i in (0)..(nvtxs) {
        let mypart = where_[i] as usize;
        rename[i] = snvtxs[mypart];
        snvtxs[mypart] += 1;
        snedges[mypart] += xadj[i + 1] - xadj[i];
    }

    let lgraph = graph::SetupSplitGraph(graph, snvtxs[0], snedges[0]);
    // sxadj[0] = lgraph.xadj;
    // svwgt[0] = lgraph.vwgt;
    // sadjncy[0] = lgraph.adjncy;
    // sadjwgt[0] = lgraph.adjwgt;
    // slabel[0] = lgraph.label;

    let rgraph = graph::SetupSplitGraph(graph, snvtxs[1], snedges[1]);
    // sxadj[1] = rgraph.xadj;
    // svwgt[1] = rgraph.vwgt;
    // sadjncy[1] = rgraph.adjncy;
    // sadjwgt[1] = rgraph.adjwgt;
    // slabel[1] = rgraph.label;
    let sxadj = {
        [
            std::slice::from_raw_parts_mut((*lgraph).xadj, snvtxs[0] as usize + 1),
            std::slice::from_raw_parts_mut((*rgraph).xadj, snvtxs[1] as usize + 1),
        ]
    };
    let svwgt = {
        [
            std::slice::from_raw_parts_mut((*lgraph).vwgt, (snvtxs[0] * ctrl.ncon) as usize),
            std::slice::from_raw_parts_mut((*rgraph).vwgt, (snvtxs[1] * ctrl.ncon) as usize),
        ]
    };
    let sadjncy = {
        [
            std::slice::from_raw_parts_mut((*lgraph).adjncy, snedges[0] as usize),
            std::slice::from_raw_parts_mut((*rgraph).adjncy, snedges[1] as usize),
        ]
    };
    let sadjwgt = {
        [
            std::slice::from_raw_parts_mut((*lgraph).adjwgt, snedges[0] as usize),
            std::slice::from_raw_parts_mut((*rgraph).adjwgt, snedges[1] as usize),
        ]
    };
    let slabel = {
        [
            std::slice::from_raw_parts_mut((*lgraph).label, snvtxs[0] as usize),
            std::slice::from_raw_parts_mut((*rgraph).label, snvtxs[1] as usize),
        ]
    };

    let mut snvtxs = [0, 0];
    let mut snedges = [0, 0];
    sxadj[0][0] = 0;
    sxadj[1][0] = 0;
    for i in (0)..(nvtxs) {
        let mypart = where_[i] as usize;

        let istart = xadj[i];
        let iend = xadj[i + 1];
        if bndptr[i] == -1 {
            /* This is an interior vertex */

            // let auxadjncy = &mut sadjncy[mypart][(snedges[mypart] - istart) as usize..];
            // let auxadjwgt = &mut sadjwgt[mypart][(snedges[mypart] - istart) as usize..];
            // for j in istart..iend {
            //     let j = j as usize;
            //     auxadjncy[j] = adjncy[j];
            //     auxadjwgt[j] = adjwgt[j];
            // }

            // TODO: look at disasm for this
            let len = iend - istart;
            // assert!(len >= 0);
            // let auxadjncy = &mut sadjncy[mypart][cntrng!(snedges[mypart], len)];
            // let auxadjwgt = &mut sadjwgt[mypart][cntrng!(snedges[mypart], len)];
            // for j in 0..len {
            //     let j = j as usize;
            //     auxadjncy[j] = adjncy[j + istart as usize];
            //     auxadjwgt[j] = adjwgt[j + istart as usize];
            // }
            for j in istart..iend {
                debug_assert_eq!(where_[adjncy[j as usize] as usize], mypart as idx_t);
            }
            sadjncy[mypart][cntrng!(snedges[mypart], len)]
                .copy_from_slice(&adjncy[cntrng!(istart, len)]);
            sadjwgt[mypart][cntrng!(snedges[mypart], len)]
                .copy_from_slice(&adjwgt[cntrng!(istart, len)]);
            snedges[mypart] += len;
        } else {
            // auxadjncy = sadjncy[mypart];
            // auxadjwgt = sadjwgt[mypart];
            let mut l = snedges[mypart] as usize;
            for j in (istart)..(iend) {
                let k = adjncy[j as usize];
                if where_[k as usize] == mypart as idx_t {
                    sadjncy[mypart][l] = k;
                    sadjwgt[mypart][l] = adjwgt[j as usize];
                    l += 1
                }
            }
            snedges[mypart] = l as idx_t;
        }

        /* copy vertex weights */
        for k in (0)..(ncon) {
            svwgt[mypart][(snvtxs[mypart] * ncon + k) as usize] =
                vwgt[(i as idx_t * ncon + k) as usize];
        }

        slabel[mypart][snvtxs[mypart] as usize] = label[i];
        snvtxs[mypart] += 1;
        sxadj[mypart][snvtxs[mypart] as usize] = snedges[mypart];
    }

    for mypart in (0)..(2) {
        let iend = sxadj[mypart][snvtxs[mypart] as usize];
        // let auxadjncy = &mut sadjncy[mypart];
        for i in (0)..(iend as usize) {
            sadjncy[mypart][i] = rename[sadjncy[mypart][i] as usize];
        }
    }

    (*lgraph).nedges = snedges[0];
    (*rgraph).nedges = snedges[1];

    graph::SetupGraph_tvwgt(lgraph);
    graph::SetupGraph_tvwgt(rgraph);

    debug_assert!(debug::check_adj(lgraph.as_ref().unwrap()));
    debug_assert!(debug::check_adj(rgraph.as_ref().unwrap()));

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.SplitTmr));

    *r_lgraph = lgraph;
    *r_rgraph = rgraph;

    // WCOREPOP;
}
