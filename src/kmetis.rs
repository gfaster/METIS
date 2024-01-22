/*!
\file
\brief The top-level routines for  multilevel k-way partitioning that minimizes
       the edge cut.

\date   Started 7/28/1997
\author George
\author Copyright 1997-2011, Regents of the University of Minnesota
\version\verbatim $Id: kmetis.c 20398 2016-11-22 17:17:12Z karypis $ \endverbatim
*/

use std::slice;

use crate::*;

// mimics the use of a goto and label in kmetis.c
// I also put it at the end of the function for the return value
// macro_rules! goto_sigthrow {
//     ($sigrval:ident) => {
//         gk_siguntrap();
//         gk_malloc_cleanup(0);
//
//         return util::metis_rcode($sigrval);
//     };
// }

/*************************************************************************/
/* This function is the entry point for MCKMETIS */
/*************************************************************************/
#[allow(non_snake_case)]
pub unsafe extern "C" fn METIS_PartGraphKway(
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
    let sigrval: std::ffi::c_int = 0;
    // let renumber: std::ffi::c_int = 0;

    /* set up malloc cleaning code and signal catchers */
    // seems like calls to malloc init begin a scope where only allocations made within the scope
    // are able to be freed (using gk_malloc/gk_free).
    if gk_malloc_init() == 0 {
        return METIS_ERROR_MEMORY;
    }
    // gk_sigtrap();

    // sigrval = gk_sigcatch();
    // if sigrval != 0 {
    //     goto_sigthrow!(sigrval);
    // }

    /* set up the run parameters */
    let ctrl: *mut ctrl_t = SetupCtrl(METIS_OP_KMETIS, options, *ncon, *nparts, tpwgts, ubvec);
    if ctrl.is_null() {
        // gk_siguntrap();
        return METIS_ERROR_INPUT;
    }
    let ctrl = ctrl.as_mut().unwrap();

    // I don't feel like it
    if ctrl.numflag == 1 {
        panic!("renumbering not supported");
    }

    /* set up the graph */
    let graph: *mut graph_t =
        graph::SetupGraph(ctrl, *nvtxs, *ncon, xadj, adjncy, vwgt, vsize, adjwgt);

    /* set up multipliers for making balance computations easier */
    SetupKWayBalMultipliers(ctrl, graph);

    /* set various run parameters that depend on the graph */
    ctrl.CoarsenTo = ((*nvtxs) / (40 * (*nparts).ilog2() as idx_t)).max(30 * (*nparts));
    ctrl.nIparts = if ctrl.nIparts != -1 {
        ctrl.nIparts
    } else if ctrl.CoarsenTo == 30 * (*nparts) {
        4
    } else {
        5
    };

    /* take care contiguity requests for disconnected graphs */
    if ctrl.contig != 0 && contig::IsConnected(graph, 0) == 0 {
        panic!(
            "METIS Error: A contiguous partition is requested for a non-contiguous input graph."
        );
    }
    /* allocate workspace memory */
    AllocateWorkSpace(ctrl, graph);

    /* start the partitioning */
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.TotalTmr));

    // iset(*nvtxs, 0, part);
    part.write_bytes(0, *nvtxs as usize);
    if ctrl.dbglvl & 512 != 0 {
        *objval = if *nparts == 1 {
            0
        } else {
            BlockKWayPartitioning(ctrl, graph, part)
        };
    } else {
        *objval = if *nparts == 1 {
            0
        } else {
            MlevelKWayPartitioning(ctrl, graph, part)
        };
    }
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.TotalTmr));
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

    /* clean up */
    FreeCtrl(&mut (ctrl as *mut ctrl_t) as *mut *mut ctrl_t);
    gk_malloc_cleanup(0);

    util::metis_rcode(sigrval)
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
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    let objval: idx_t = 0;
    let mut bestobj: idx_t = 0;
    let mut bestbal: real_t = 0.0;
    // let status: i32;

    for i in 0..ctrl.ncuts {
        let cgraph = coarsen::CoarsenGraph(ctrl, graph).as_mut().unwrap();

        // ifset!(
        //     ctrl.dbglvl,
        //     METIS_DBG_TIME,
        //     gk_startcputimer(ctrl.InitPartTmr),
        // );
        kwayrefine::AllocateKWayPartitionMemory(ctrl, cgraph);

        /* Release the work space */
        FreeWorkSpace(ctrl);

        /* Compute the initial partitioning */
        InitKWayPartitioning(ctrl, cgraph);

        /* Re-allocate the work space */
        AllocateWorkSpace(ctrl, graph);
        AllocateRefinementWorkSpace(ctrl, graph.nedges, 2 * cgraph.nedges);

        // ifset!(
        //     ctrl.dbglvl,
        //     METIS_DBG_TIME,
        //     gk_stopcputimer(ctrl.InitPartTmr)
        // );
        ifset!(
            ctrl.dbglvl,
            METIS_DBG_IPART,
            println!("Initial {:}-way partitioning cut: {:}", ctrl.nparts, objval)
        );

        kwayrefine::RefineKWay(ctrl, graph, cgraph);

        let curobj = match ctrl.objtype {
            METIS_OBJTYPE_CUT => graph.mincut,
            METIS_OBJTYPE_VOL => graph.minvol,
            _ => panic!("Unknown objtype: {}", ctrl.objtype),
        };

        let curbal = mcutil::ComputeLoadImbalanceDiff(graph, ctrl.nparts, ctrl.pijbm, ctrl.ubfactors);

        if i == 0
            || (curbal <= 0.0005 && bestobj > curobj)
            || (bestbal > 0.0005 && curbal < bestbal)
        {
            // icopy(graph.nvtxs, graph.where_, part);
            part.copy_from(graph.where_, graph.nvtxs as usize);
            bestobj = curobj;
            bestbal = curbal;
        }

        graph::FreeRData(graph);

        if bestobj == 0 {
            break;
        }
    }

    graph::FreeGraph(&mut (graph as *mut graph_t) as *mut *mut graph_t);

    bestobj
}

/*************************************************************************/
/* This function computes the initial k-way partitioning using PMETIS
*/
/*************************************************************************/
#[metis_func]
pub fn InitKWayPartitioning(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // let ntrials: idx_t;
    let mut options: [idx_t; METIS_NOPTIONS as usize];
    let mut curobj: idx_t = 0;
    // let bestobj: idx_t = 0;
    // let bestwhere_: *mut idx_t = ptr::null_mut();
    let status: std::ffi::c_int;

    options = [-1; METIS_NOPTIONS as usize];
    METIS_SetDefaultOptions(options.as_mut_ptr());
    //options[METIS_OPTION_NITER]     = 10;
    options[METIS_OPTION_NITER as usize] = ctrl.niter;
    options[METIS_OPTION_OBJTYPE as usize] = METIS_OBJTYPE_CUT as idx_t;
    options[METIS_OPTION_NO2HOP as usize] = ctrl.no2hop;
    options[METIS_OPTION_ONDISK as usize] = ctrl.ondisk;
    options[METIS_OPTION_DROPEDGES as usize] = ctrl.dropedges;
    //options[METIS_OPTION_DBGLVL]    = ctrl.dbglvl;

    // ubvec = rmalloc(graph.ncon, "InitKWayPartitioning: ubvec");
    let mut ubvec = vec![0.0 as real_t; graph.ncon as usize];
    mkslice_mut!(ctrl->ubfactors, graph.ncon);
    for i in 0..(graph.ncon as usize) {
        ubvec[i] = (ubfactors[i] as f64).powf(1.0 / (ctrl.nparts as f64).ln()) as real_t;
    }

    match ctrl.objtype {
        METIS_OBJTYPE_CUT | METIS_OBJTYPE_VOL => {
            options[METIS_OPTION_NCUTS as usize] = ctrl.nIparts;
            status = pmetis::METIS_PartGraphRecursive(
                &mut graph.nvtxs,
                &mut graph.ncon,
                graph.xadj,
                graph.adjncy,
                graph.vwgt,
                graph.vsize,
                graph.adjwgt,
                &mut ctrl.nparts,
                ctrl.tpwgts,
                ubvec.as_mut_ptr(),
                options.as_mut_ptr(),
                &mut curobj,
                graph.where_,
            );

            if status != METIS_OK {
                panic!("Failed during initial partitioning");
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
        //         if status != METIS_OK
        //           panic!("Failed during initial partitioning");
        //
        //         curobj = ComputeVolume(graph, graph.where_);
        //
        //         if i == 0 || bestobj > curobj {
        //           bestobj = curobj;
        //           if i < ntrials-1
        //             icopy(graph.nvtxs, graph.where_, bestwhere_);
        //         }
        //
        //         if bestobj == 0
        //           break;
        //       }
        //       if bestobj != curobj
        //         icopy(graph.nvtxs, bestwhere_, graph.where_);
        //
        //       break;
        // #endif
        _ => panic!("Unknown objtype: {}", ctrl.objtype),
    }

    // don't need this - ubvec is RAII, bestwhere_ is commented out
    // gk_free(&mut ubvec as *mut *mut c_void, &mut bestwhere_, LTERM);
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

    let mut queue: pqueue::IPQueue;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices_mut!(graph => vwgt);

    let nparts: idx_t = ctrl.nparts;

    let mynparts = (100 * nparts).min((nvtxs as f32).sqrt() as idx_t);
    mkslice_mut!(part, nvtxs);

    for i in 0..nvtxs {
        part[i] = (i % nparts as usize) as idx_t;
    }
    irandArrayPermute(nvtxs as idx_t, part.as_mut_ptr(), 4 * nvtxs as idx_t, 0);
    println!("Random cut: {}", ComputeCut(graph, part.as_mut_ptr()));

    /* create the initial multi-section */
    let mynparts = GrowMultisection(ctrl, graph, mynparts, part.as_mut_ptr());

    /* balance using label-propagation and refine using a randomized greedy strategy */
    BalanceAndRefineLP(ctrl, graph, mynparts, part.as_mut_ptr());

    /* determine the size of the fine partitions */
    // fpwgts = iset(mynparts, 0, iwspacemalloc(ctrl, mynparts));
    let mut fpwgts = vec![0; mynparts as usize];
    for i in 0..nvtxs {
        fpwgts[part[i] as usize] += vwgt[i];
    }

    /* create and initialize the queue that will determine
    where to put the next one */
    // cpwgts = iset(nparts, 0, iwspacemalloc(ctrl, nparts));
    let mut cpwgts = vec![0; nparts as usize];

    // I am not certain about the behavior of the queue - this will have the most potential for
    // error
    // queue = ipqCreate(nparts);
    queue = pqueue::IPQueue::new(nparts as usize);
    for i in 0..nparts {
        // ipqInsert(queue, i, 0);
        queue.insert(i, 0);
    }

    /* assign the fine partitions into the coarse partitions */
    // fpart = iwspacemalloc(ctrl, mynparts);
    let mut fpart = vec![0; mynparts as usize];
    // perm = iwspacemalloc(ctrl, mynparts);
    let mut perm = vec![0; mynparts as usize];
    irandArrayPermute(mynparts, perm.as_mut_ptr(), mynparts, 1);
    for ii in 0..(mynparts as usize) {
        let i = perm[ii];
        // let j = ipqSeeTopVal(queue);
        let j = queue.peek_val().unwrap_or(-1);
        fpart[i as usize] = j;
        cpwgts[j as usize] += fpwgts[i as usize];
        // ipqUpdate(queue, j, -cpwgts[j]);
        queue.update(j, -cpwgts[j as usize]);
    }
    // ipqDestroy(queue);
    drop(queue);

    for i in 0..(nparts as usize) {
        println!("cpwgts[{}] = {}", i, cpwgts[i]);
    }
    for i in 0..nvtxs {
        part[i] = fpart[part[i] as usize];
    }
    // WCOREPOP;

    ComputeCut(graph, part.as_mut_ptr())
}

/*************************************************************************/
/* This function takes a graph and produces a bisection by using a region
    growing algorithm. The resulting bisection is refined using FM.
    The resulting partition is returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn GrowMultisection(
    _ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    nparts: idx_t,
    where_: *mut idx_t,
) -> idx_t {
    let graph = graph.as_mut().unwrap();
    // let ctrl = ctrl.as_mut().unwrap();
    let mut nparts = nparts;

    // WCOREPUSH;

    let nvtxs: idx_t = graph.nvtxs;
    // xadj = graph.xadj;
    // vwgt = graph.vwgt;
    // adjncy = graph.adjncy;
    get_graph_slices!(graph => xadj vwgt adjncy);

    // get_graph_slices_mut!(graph => where_);
    // using mkslice instead of graph since where_ was passed as argument
    mkslice_mut!(where_, nvtxs);

    // queue = iwspacemalloc(ctrl, nvtxs);
    let mut queue = vec![0; nvtxs as usize];

    /* Select the seeds for the nparts-way BFS */
    let mut nleft = 0;
    for i in 0..(nvtxs as usize) {
        if xadj[i + 1] - xadj[i] > 1
        /* a seed's degree should be > 1 */
        {
            where_[nleft as usize] = i as idx_t;
            nleft += 1;
        }
    }
    nparts = nparts.min(nleft);
    for i in 0..(nparts as usize) {
        let j = irandInRange(nleft) as usize;
        queue[i] = where_[j];
        nleft -= 1;
        where_[j] = nleft;
    }

    let mut pwgts = vec![0; nparts as usize];
    let tvwgt: idx_t = vwgt.iter().sum();
    // tvwgt = isum(nvtxs, vwgt, 1);
    let maxpwgt: idx_t = ((1.5 * tvwgt as f32) / nparts as f32) as idx_t;

    where_.fill(-1);
    for i in 0..(nparts as usize) {
        where_[queue[i] as usize] = i as idx_t;
        pwgts[i] = vwgt[queue[i] as usize];
    }

    let mut first = 0_usize;
    let mut last = nparts as usize;
    nleft = nvtxs - nparts;

    /* Start the BFS from queue to get a partition */
    while first < last {
        let i = queue[first] as usize;
        first += 1;
        let l = where_[i] as usize;
        if pwgts[l] > maxpwgt {
            continue;
        }
        for j in xadj[i]..xadj[i + 1] {
            let j = j as usize;
            let k = adjncy[j] as usize;
            if where_[k] == -1 {
                if pwgts[l] + vwgt[k] > maxpwgt {
                    break;
                }
                pwgts[l] += vwgt[k];
                where_[k] = l as idx_t;
                queue[last] = k as idx_t;
                last += 1;
                nleft -= 1;
            }
        }
    }

    /* Assign the unassigned vertices randomly to the nparts partitions */
    if nleft > 0 {
        for i in 0..(nvtxs as usize) {
            if where_[i] == -1 {
                where_[i] = irandInRange(nparts);
            }
        }
    }

    // WCOREPOP;

    nparts
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
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(ctrl, graph => xadj vwgt adjncy adjwgt);
    let where_ = slice::from_raw_parts_mut(where_, nvtxs);

    let mut pwgts = vec![0; nparts as usize];

    let ubfactor: real_t = util::i2rubfactor(ctrl.ufactor);
    let tvwgt: idx_t = vwgt.iter().sum::<idx_t>();
    let maxpwgt: idx_t = ((ubfactor * tvwgt as real_t) / nparts as real_t) as idx_t;
    let minpwgt: idx_t = ((1.0 * tvwgt as real_t) / (ubfactor * nparts as real_t)) as idx_t;

    for i in 0..nvtxs {
        pwgts[where_[i] as usize] += vwgt[i];
    }
    /* for randomly visiting the vertices */
    // perm = iincset(nvtxs, 0, iwspacemalloc(ctrl, nvtxs));
    let mut perm = Vec::from_iter(0..(nvtxs as idx_t));

    /* for keeping track of adjacent partitions */
    let mut nbrids: Vec<idx_t> = vec![0; nparts as usize];
    let mut nbrwgts: Vec<idx_t> = vec![0; nparts as usize];
    let mut nbrmrks: Vec<idx_t> = vec![-1; nparts as usize];

    /* perform a fixed number of balancing LP iterations */
    if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
        println!(
            "BLP: nparts: {:}, min-max: [{:}, {:}], bal: {:7.4}, cut: {:9}",
            nparts,
            minpwgt,
            maxpwgt,
            *pwgts.iter().max().unwrap_or(&0) as f32 * nparts as f32 / tvwgt as f32,
            // 1.0 * imax(nparts, pwgts, 1) * nparts / tvwgt,
            ComputeCut(graph, where_.as_mut_ptr())
        );
    }
    for _ in 0..ctrl.niter {
        // if imax(nparts, pwgts, 1) * nparts < ubfactor * tvwgt {
        if (*pwgts.iter().max().unwrap_or(&0) as f32) * (nparts as f32) < ubfactor * (tvwgt as f32)
        {
            break;
        }
        // I may choose to recreate this macro impl if it's too slow
        // irandArrayPermute(nvtxs, perm, nvtxs / 8, 1);
        fastrand::shuffle(&mut perm);
        let mut nmoves = 0;

        for ii in 0..nvtxs {
            let u = perm[ii] as usize;

            let from = where_[u] as usize;
            if pwgts[from] - vwgt[u] < minpwgt {
                continue;
            }
            let mut nnbrs = 0;
            for j in xadj[u]..xadj[u + 1] {
                let j = j as usize;
                let v = adjncy[j] as usize;
                let to = where_[v] as usize;

                if pwgts[to] + vwgt[u] > maxpwgt {
                    continue; /* skip if 'to' is overweight */
                }
                let mut k = nbrmrks[to];
                if k == -1 {
                    k = nnbrs;
                    nbrmrks[to] = nnbrs;
                    nnbrs += 1;
                    nbrids[k as usize] = to as idx_t;
                }
                nbrwgts[k as usize] += xadj[v + 1] - xadj[v];
            }
            if nnbrs == 0 {
                continue;
            }
            let to = nbrids[util::iargmax(&nbrwgts, 1)] as usize;
            if from != to {
                where_[u] = to as idx_t;
                inc_dec!(pwgts[to], pwgts[from], vwgt[u]);
                nmoves += 1;
            }

            for k in 0..(nnbrs as usize) {
                nbrmrks[nbrids[k] as usize] = -1;
                nbrwgts[k] = 0;
            }
        }

        if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
            println!(
                "     nmoves: {:8}, bal: {:7.4}, cut: {:9}",
                nmoves,
                *pwgts.iter().max().unwrap_or(&0) as f32 * nparts as f32 / tvwgt as f32,
                ComputeCut(graph, where_.as_mut_ptr())
            );
        }
        if nmoves == 0 {
            break;
        }
    }

    /* perform a fixed number of refinement LP iterations */
    if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
        println!(
            "RLP: nparts: {:}, min-max: [{:}, {:}], bal: {:7.4}, cut: {:9}",
            nparts,
            minpwgt,
            maxpwgt,
            *pwgts.iter().max().unwrap_or(&0) as f32 * nparts as f32 / tvwgt as f32,
            ComputeCut(graph, where_.as_mut_ptr())
        );
    }
    for _ in 0..ctrl.niter {
        irandArrayPermute(nvtxs as idx_t, perm.as_mut_ptr(), nvtxs as idx_t / 8, 1);
        let mut nmoves = 0;

        for ii in 0..nvtxs {
            let u = perm[ii] as usize;

            let from = where_[u] as usize;
            if pwgts[from] - vwgt[u] < minpwgt {
                continue;
            }
            let mut nnbrs = 0;
            for j in xadj[u]..xadj[u + 1] {
                let j = j as usize;
                let v = adjncy[j] as usize;
                let to = where_[v] as usize;

                if to != from && pwgts[to] + vwgt[u] > maxpwgt {
                    continue; /* skip if 'to' is overweight */
                }
                let mut k = nbrmrks[to];
                if k == -1 {
                    k = nnbrs;
                    nbrmrks[to] = nnbrs;
                    nnbrs += 1;

                    nbrids[k as usize] = to as idx_t;
                }
                nbrwgts[k as usize] += adjwgt[j];
            }
            if nnbrs == 0 {
                continue;
            }
            let to = nbrids[util::iargmax(&nbrwgts, 1)] as usize;
            if from != to {
                where_[u] = to as idx_t;
                inc_dec!(pwgts[to], pwgts[from], vwgt[u]);
                nmoves += 1;
            }

            for k in 0..(nnbrs as usize) {
                nbrmrks[nbrids[k] as usize] = -1;
                nbrwgts[k] = 0;
            }
        }

        if ctrl.dbglvl & METIS_DBG_REFINE != 0 {
            println!(
                "     nmoves: {:8}, bal: {:7.4}, cut: {:9}",
                nmoves,
                *pwgts.iter().max().unwrap_or(&0) as f32 * nparts as f32 / tvwgt as f32,
                ComputeCut(graph, where_.as_mut_ptr())
            );
        }
        if nmoves == 0 {
            break;
        }
    }

    // WCOREPOP;
}
