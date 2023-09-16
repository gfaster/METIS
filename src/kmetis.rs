/*!
\file
\brief The top-level routines for  multilevel k-way partitioning that minimizes
       the edge cut.

\date   Started 7/28/1997
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota
\version\verbatim $Id: kmetis.c 20398 2016-11-22 17:17:12Z karypis $ \endverbatim
*/

use std::{os::raw::c_void, ptr};

use crate::*;

macro_rules! goto_sigthrow {
    () => {
        gk_siguntrap();
        gk_malloc_cleanup(0);

        return metis_rcode(sigrval);
    };
}

/*************************************************************************/
/* This function is the entry point for MCKMETIS */
/*************************************************************************/
#[metis_func]
pub fn METIS_PartGraphKway(
    nvtxs: *mut idx_t,
    ncon: *mut idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    vsize: *mut idx_t,
    adjwgt: *mut idx_t,
    nparts: *mut idx_t,
    tpwgts: *mut real_t,
    ubvec: *mut real_t,
    options: *mut idx_t,
    objval: *mut idx_t,
    part: *mut idx_t,
) -> std::ffi::c_int {
    let sigrval: int = 0;
    let renumber: int = 0;
    graph_t * graph;
    ctrl_t * ctrl;

    /* set up malloc cleaning code and signal catchers */
    if (!gk_malloc_init()) {
        return METIS_ERROR_MEMORY;
    }
    gk_sigtrap();

    if ((sigrval = gk_sigcatch()) != 0) {
        goto_sigthrow!();
    }
    /* set up the run parameters */
    ctrl = SetupCtrl(METIS_OP_KMETIS, options, *ncon, *nparts, tpwgts, ubvec);
    if (!ctrl) {
        gk_siguntrap();
        return METIS_ERROR_INPUT;
    }

    /* if required, change the numbering to 0 */
    if (ctrl.numflag == 1) {
        Change2CNumbering(*nvtxs, xadj, adjncy);
        renumber = 1;
    }

    /* set up the graph */
    graph = SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);

    /* set up multipliers for making balance computations easier */
    SetupKWayBalMultipliers(ctrl, graph);

    /* set various run parameters that depend on the graph */
    ctrl.CoarsenTo = gk_max((*nvtxs) / (40 * gk_log2(*nparts)), 30 * (*nparts));
    ctrl.nIparts = (if ctrl.nIparts != -1 {
        ctrl.nIparts
    } else {
        (if ctrl.CoarsenTo == 30 * (*nparts) {
            4
        } else {
            5
        })
    });

    /* take care contiguity requests for disconnected graphs */
    if (ctrl.contig && !IsConnected(graph, 0)) {
        panic!(
            "METIS Error: A contiguous partition is requested for a non-contiguous input graph.\n"
        );
    }
    /* allocate workspace memory */
    AllocateWorkSpace(ctrl, graph);

    /* start the partitioning */
    IFSET(ctrl.dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

    iset(*nvtxs, 0, part);
    if (ctrl.dbglvl & 512) {
        *objval = (if *nparts == 1 {
            0
        } else {
            BlockKWayPartitioning(ctrl, graph, part)
        });
    } else {
        *objval = (if *nparts == 1 {
            0
        } else {
            MlevelKWayPartitioning(ctrl, graph, part)
        });
    }
    IFSET(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
    IFSET(ctrl.dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

    /* clean up */
    FreeCtrl(&mut ctrl);
}

/*************************************************************************/
/* This function computes a k-way partitioning of a graph that minimizes
    the specified objective function.

    \param ctrl is the control structure
    \param graph is the graph to be partitioned
    \param part is the vector that on return will store the partitioning

    \returns the objective value of the partitioning. The partitioning
             itself is stored in the part vector.
*/
/*************************************************************************/
#[metis_func]
pub fn MlevelKWayPartitioning(ctrl: *mut ctrl_t, graph: *mut graph_t, part: *mut idx_t) -> idx_t {
    let i: idx_t;
    let j: idx_t;
    let objval: idx_t = 0;
    let curobj: idx_t = 0;
    let bestobj: idx_t = 0;
    let curbal: real_t = 0.0;
    let bestbal: real_t = 0.0;
    let cgraph: *mut graph_t;
    let status: int;

    for i in 0..ctrl.ncuts {
        cgraph = CoarsenGraph(ctrl, graph);

        IFSET(
            ctrl.dbglvl,
            METIS_DBG_TIME,
            gk_startcputimer(ctrl.InitPartTmr),
        );
        AllocateKWayPartitionMemory(ctrl, cgraph);

        /* Release the work space */
        FreeWorkSpace(ctrl);

        /* Compute the initial partitioning */
        InitKWayPartitioning(ctrl, cgraph);

        /* Re-allocate the work space */
        AllocateWorkSpace(ctrl, graph);
        AllocateRefinementWorkSpace(ctrl, graph.nedges, 2 * cgraph.nedges);

        IFSET(
            ctrl.dbglvl,
            METIS_DBG_TIME,
            gk_stopcputimer(ctrl.InitPartTmr),
        );
        IFSET(
            ctrl.dbglvl,
            METIS_DBG_IPART,
            print!(
                "Initial {:}-way partitioning cut: {:}\n",
                ctrl.nparts, objval
            ),
        );

        RefineKWay(ctrl, graph, cgraph);

        match (ctrl.objtype) {
            METIS_OBJTYPE_CUT => {
                curobj = graph.mincut;
            }
            METIS_OBJTYPE_VOL => {
                curobj = graph.minvol;
            }
            _ => panic!("Unknown objtype: {}\n", ctrl.objtype),
        }

        curbal = ComputeLoadImbalanceDiff(graph, ctrl.nparts, ctrl.pijbm, ctrl.ubfactors);

        if (i == 0
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal))
        {
            icopy(graph.nvtxs, graph.where_, part);
            bestobj = curobj;
            bestbal = curbal;
        }

        FreeRData(graph);

        if (bestobj == 0) {
            break;
        }
    }

    FreeGraph(&mut graph);

    return bestobj;
}

/*************************************************************************/
/* This function computes the initial k-way partitioning using PMETIS
*/
/*************************************************************************/
#[metis_func]
pub fn InitKWayPartitioning(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let i: idx_t;
    let ntrials: idx_t;
    let options: [idx_t; METIS_NOPTIONS];
    let curobj: idx_t = 0;
    let bestobj: idx_t = 0;
    let bestwhere_: *mut idx_t = ptr::null_mut();
    let ubvec: *mut real_t = ptr::null_mut();
    let status: std::ffi::c_int;

    METIS_SetDefaultOptions(options);
    //options[METIS_OPTION_NITER]     = 10;
    options[METIS_OPTION_NITER] = ctrl.niter;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NO2HOP] = ctrl.no2hop;
    options[METIS_OPTION_ONDISK] = ctrl.ondisk;
    options[METIS_OPTION_DROPEDGES] = ctrl.dropedges;
    //options[METIS_OPTION_DBGLVL]    = ctrl.dbglvl;

    ubvec = rmalloc(graph.ncon, "InitKWayPartitioning: ubvec");
    for i in 0..graph.ncon {
        ubvec[i] = ctrl.ubfactors[i].pow(1.0 / log(ctrl.nparts)) as real_t;
    }

    match (ctrl.objtype) {
        METIS_OBJTYPE_CUT | METIS_OBJTYPE_VOL => {
            options[METIS_OPTION_NCUTS] = ctrl.nIparts;
            status = METIS_PartGraphRecursive(
                &mut graph.nvtxs,
                &mut graph.ncon,
                graph.xadj,
                graph.adjncy,
                graph.vwgt,
                graph.vsize,
                graph.adjwgt,
                &mut ctrl.nparts,
                ctrl.tpwgts,
                ubvec,
                options,
                &mut curobj,
                graph.where_,
            );

            if (status != METIS_OK) {
                panic!("Failed during initial partitioning\n");
            }
        }

        // #ifdef XXX /* This does not seem to help */
        //     METIS_OBJTYPE_VOL =>
        //       bestwhere_ = imalloc(graph.nvtxs, "InitKWayPartitioning: bestwhere_");
        //       options[METIS_OPTION_NCUTS] = 2;
        //
        //       ntrials = (ctrl.nIparts+1)/2;
        //       for i in 0..ntrials {
        //         status = METIS_PartGraphRecursive(&mut graph.nvtxs, &graph.ncon,
        //                      graph.xadj, graph.adjncy, graph.vwgt, graph.vsize,
        //                      graph.adjwgt, &mut ctrl.nparts, ctrl.tpwgts, ubvec,
        //                      options, &mut curobj, graph.where_);
        //         if (status != METIS_OK)
        //           panic!("Failed during initial partitioning\n");
        //
        //         curobj = ComputeVolume(graph, graph.where_);
        //
        //         if (i == 0 || bestobj > curobj) {
        //           bestobj = curobj;
        //           if (i < ntrials-1)
        //             icopy(graph.nvtxs, graph.where_, bestwhere_);
        //         }
        //
        //         if (bestobj == 0)
        //           break;
        //       }
        //       if (bestobj != curobj)
        //         icopy(graph.nvtxs, bestwhere_, graph.where_);
        //
        //       break;
        // #endif
        _ => panic!("Unknown objtype: {}\n", ctrl.objtype),
    }

    gk_free(&mut ubvec as *mut *mut c_void, &mut bestwhere_, LTERM);
}

/*************************************************************************/
/* This function computes a k-way partitioning of a graph that minimizes
    the specified objective function.

    \param ctrl is the control structure
    \param graph is the graph to be partitioned
    \param part is the vector that on return will store the partitioning

    \returns the objective value of the partitioning. The partitioning
             itself is stored in the part vector.
*/
/*************************************************************************/
#[metis_func]
pub fn BlockKWayPartitioning(ctrl: *mut ctrl_t, graph: *mut graph_t, part: *mut idx_t) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let i: idx_t;
    let ii: idx_t;
    let j: idx_t;
    let nvtxs: idx_t;
    let objval: idx_t = 0;
    let vwgt: *mut idx_t;
    let nparts: idx_t;
    let mynparts: idx_t;
    let fpwgts: *mut idx_t;
    let cpwgts: *mut idx_t;
    let fpart: *mut idx_t;
    let perm: *mut idx_t;
    ipq_t * queue;

    WCOREPUSH;

    nvtxs = graph.nvtxs;
    vwgt = graph.vwgt;

    nparts = ctrl.nparts;

    mynparts = gk_min(100 * nparts, sqrt(nvtxs));

    for i in 0..nvtxs {
        part[i] = i % nparts;
    }
    irandArrayPermute(nvtxs, part, 4 * nvtxs, 0);
    print!("Random cut: {}\n", (int)ComputeCut(graph, part));

    /* create the initial multi-section */
    mynparts = GrowMultisection(ctrl, graph, mynparts, part);

    /* balance using label-propagation and refine using a randomized greedy strategy */
    BalanceAndRefineLP(ctrl, graph, mynparts, part);

    /* determine the size of the fine partitions */
    fpwgts = iset(mynparts, 0, iwspacemalloc(ctrl, mynparts));
    for i in 0..nvtxs {
        fpwgts[part[i]] += vwgt[i];
    }

    /* create and initialize the queue that will determine
    where_ to put the next one */
    cpwgts = iset(nparts, 0, iwspacemalloc(ctrl, nparts));
    queue = ipqCreate(nparts);
    for i in 0..nparts {
        ipqInsert(queue, i, 0);
    }

    /* assign the fine partitions into the coarse partitions */
    fpart = iwspacemalloc(ctrl, mynparts);
    perm = iwspacemalloc(ctrl, mynparts);
    irandArrayPermute(mynparts, perm, mynparts, 1);
    for ii in 0..mynparts {
        i = perm[ii];
        j = ipqSeeTopVal(queue);
        fpart[i] = j;
        cpwgts[j] += fpwgts[i];
        ipqUpdate(queue, j, -cpwgts[j]);
    }
    ipqDestroy(queue);

    for i in 0..nparts {
        print!("cpwgts[{}] = {}\n", i, cpwgts[i]);
    }
    for i in 0..nvtxs {
        part[i] = fpart[part[i]];
    }
    WCOREPOP;

    return ComputeCut(graph, part);
}

/*************************************************************************/
/* This function takes a graph and produces a bisection by using a region
    growing algorithm. The resulting bisection is refined using FM.
    The resulting partition is returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn GrowMultisection(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    nparts: idx_t,
    where_: *mut idx_t,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let i: idx_t;
    let j: idx_t;
    let k: idx_t;
    let l: idx_t;
    let nvtxs: idx_t;
    let nleft: idx_t;
    let first: idx_t;
    let last: idx_t;
    let xadj: *mut idx_t;
    let vwgt: *mut idx_t;
    let adjncy: *mut idx_t;
    let queue: *mut idx_t;
    let tvwgt: idx_t;
    let maxpwgt: idx_t;
    let pwgts: *mut idx_t;

    WCOREPUSH;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    vwgt = graph.xadj;
    adjncy = graph.adjncy;

    queue = iwspacemalloc(ctrl, nvtxs);

    /* Select the seeds for the nparts-way BFS */
    nleft = 0;
    for i in 0..nvtxs {
        if (xadj[i + 1] - xadj[i] > 1)
        /* a seed's degree should be > 1 */
        {
            where_[nleft] = i;
            nleft += 1;
        }
    }
    nparts = gk_min(nparts, nleft);
    for i in 0..nparts {
        j = irandInRange(nleft);
        queue[i] = where_[j];
        where_[j] = --nleft;
    }

    pwgts = iset(nparts, 0, iwspacemalloc(ctrl, nparts));
    tvwgt = isum(nvtxs, vwgt, 1);
    maxpwgt = (1.5 * tvwgt) / nparts;

    iset(nvtxs, -1, where_);
    for i in 0..nparts {
        where_[queue[i]] = i;
        pwgts[i] = vwgt[queue[i]];
    }

    first = 0;
    last = nparts;
    nleft = nvtxs - nparts;

    /* Start the BFS from queue to get a partition */
    while (first < last) {
        i = queue[first];
        first += 1;
        l = where_[i];
        if (pwgts[l] > maxpwgt) {
            continue;
        }
        for j in xadj[i]..xadj[i + 1] {
            k = adjncy[j];
            if (where_[k] == -1) {
                if (pwgts[l] + vwgt[k] > maxpwgt) {
                    break;
                }
                pwgts[l] += vwgt[k];
                where_[k] = l;
                queue[last] = k;
                last += 1;
                nleft -= 1;
            }
        }
    }

    /* Assign the unassigned vertices randomly to the nparts partitions */
    if (nleft > 0) {
        for i in 0..nvtxs {
            if (where_[i] == -1) {
                where_[i] = irandInRange(nparts);
            }
        }
    }

    WCOREPOP;

    return nparts;
}

/*************************************************************************/
/* This function balances the partitioning using label propagation.
*/
/*************************************************************************/
#[metis_func]
pub fn BalanceAndRefineLP(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    nparts: idx_t,
    where_: *mut idx_t,
) -> () {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    let ii: idx_t;
    let i: idx_t;
    let j: idx_t;
    let k: idx_t;
    let u: idx_t;
    let v: idx_t;
    let nvtxs: idx_t;
    let iter: idx_t;
    let xadj: *mut idx_t;
    let vwgt: *mut idx_t;
    let adjncy: *mut idx_t;
    let adjwgt: *mut idx_t;
    let tvwgt: idx_t;
    let pwgts: *mut idx_t;
    let maxpwgt: idx_t;
    let minpwgt: idx_t;
    let perm: *mut idx_t;
    let from: idx_t;
    let to: idx_t;
    let nmoves: idx_t;
    let nnbrs: idx_t;
    let nbrids: *mut idx_t;
    let nbrwgts: *mut idx_t;
    let nbrmrks: *mut idx_t;
    let ubfactor: real_t;

    WCOREPUSH;

    nvtxs = graph.nvtxs;
    xadj = graph.xadj;
    vwgt = graph.vwgt;
    adjncy = graph.adjncy;
    adjwgt = graph.adjwgt;

    pwgts = iset(nparts, 0, iwspacemalloc(ctrl, nparts));

    ubfactor = I2RUBFACTOR(ctrl.ufactor);
    tvwgt = isum(nvtxs, vwgt, 1);
    maxpwgt = (ubfactor * tvwgt) / nparts;
    minpwgt = (1.0 * tvwgt) / (ubfactor * nparts);

    for i in 0..nvtxs {
        pwgts[where_[i]] += vwgt[i];
    }
    /* for randomly visiting the vertices */
    perm = iincset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));

    /* for keeping track of adjacent partitions */
    nbrids = iwspacemalloc(ctrl, nparts);
    nbrwgts = iset(nparts, 0, iwspacemalloc(ctrl, nparts));
    nbrmrks = iset(nparts, -1, iwspacemalloc(ctrl, nparts));

    /* perform a fixed number of balancing LP iterations */
    if (ctrl.dbglvl & METIS_DBG_REFINE) {
        print!(
            "BLP: nparts: {:}, min-max: [{:}, {:}], bal: {:7.4}, cut: {:9}\n",
            nparts,
            minpwgt,
            maxpwgt,
            1.0 * imax(nparts, pwgts, 1) * nparts / tvwgt,
            ComputeCut(graph, where_)
        );
    }
    for iter in 0..ctrl.niter {
        if (imax(nparts, pwgts, 1) * nparts < ubfactor * tvwgt) {
            break;
        }
        irandArrayPermute(nvtxs, perm, nvtxs / 8, 1);
        nmoves = 0;

        for ii in 0..nvtxs {
            u = perm[ii];

            from = where_[u];
            if (pwgts[from] - vwgt[u] < minpwgt) {
                continue;
            }
            nnbrs = 0;
            for j in xadj[u]..xadj[u + 1] {
                v = adjncy[j];
                to = where_[v];

                if (pwgts[to] + vwgt[u] > maxpwgt) {
                    continue; /* skip if 'to' is overweight */
                }
                if ((k = nbrmrks[to]) == -1) {
                    k = nnbrs;
                    nbrmrks[to] = nnbrs;
                    nnbrs += 1;
                    nbrids[k] = to;
                }
                nbrwgts[k] += xadj[v + 1] - xadj[v];
            }
            if (nnbrs == 0) {
                continue;
            }
            to = nbrids[iargmax(nnbrs, nbrwgts, 1)];
            if (from != to) {
                where_[u] = to;
                INC_DEC(pwgts[to], pwgts[from], vwgt[u]);
                nmoves += 1;
            }

            for k in 0..nnbrs {
                nbrmrks[nbrids[k]] = -1;
                nbrwgts[k] = 0;
            }
        }

        if (ctrl.dbglvl & METIS_DBG_REFINE) {
            print!(
                "     nmoves: {:8}, bal: {:7.4}, cut: {:9}\n",
                nmoves,
                1.0 * imax(nparts, pwgts, 1) * nparts / tvwgt,
                ComputeCut(graph, where_)
            );
        }
        if (nmoves == 0) {
            break;
        }
    }

    /* perform a fixed number of refinement LP iterations */
    if (ctrl.dbglvl & METIS_DBG_REFINE) {
        print!(
            "RLP: nparts: {:}, min-max: [{:}, {:}], bal: {:7.4}, cut: {:9}\n",
            nparts,
            minpwgt,
            maxpwgt,
            1.0 * imax(nparts, pwgts, 1) * nparts / tvwgt,
            ComputeCut(graph, where_)
        );
    }
    for iter in 0..ctrl.niter {
        irandArrayPermute(nvtxs, perm, nvtxs / 8, 1);
        nmoves = 0;

        for ii in 0..nvtxs {
            u = perm[ii];

            from = where_[u];
            if (pwgts[from] - vwgt[u] < minpwgt) {
                continue;
            }
            nnbrs = 0;
            for j in xadj[u]..xadj[u + 1] {
                v = adjncy[j];
                to = where_[v];

                if (to != from && pwgts[to] + vwgt[u] > maxpwgt) {
                    continue; /* skip if 'to' is overweight */
                }
                if ((k = nbrmrks[to]) == -1) {
                    k = nnbrs;
                    nbrmrks[to] = nnbrs;
                    nnbrs += 1;

                    nbrids[k] = to;
                }
                nbrwgts[k] += adjwgt[j];
            }
            if (nnbrs == 0) {
                continue;
            }
            to = nbrids[iargmax(nnbrs, nbrwgts, 1)];
            if (from != to) {
                where_[u] = to;
                INC_DEC(pwgts[to], pwgts[from], vwgt[u]);
                nmoves += 1;
            }

            for k in 0..nnbrs {
                nbrmrks[nbrids[k]] = -1;
                nbrwgts[k] = 0;
            }
        }

        if (ctrl.dbglvl & METIS_DBG_REFINE) {
            print!(
                "     nmoves: {:8}, bal: {:7.4}, cut: {:9}\n",
                nmoves,
                1.0 * imax(nparts, pwgts, 1) * nparts / tvwgt,
                ComputeCut(graph, where_)
            );
        }
        if (nmoves == 0) {
            break;
        }
    }

    WCOREPOP;
}
