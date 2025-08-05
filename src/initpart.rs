/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * initpart.c
 *
 * This file contains code that performs the initial partition of the
 * coarsest graph
 *
 * Started 7/23/97
 * George
 *
 */

use crate::*;

/*************************************************************************/
/* This function computes the initial bisection of the coarsest graph */
/*************************************************************************/
#[metis_func]
pub fn Init2WayPartition(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niparts: idx_t,
) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    get_graph_slices!(graph => tvwgt);

    assert!(tvwgt[0] >= 0);

    let dbglvl = ctrl.dbglvl;
    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        ctrl.dbglvl -= METIS_DBG_REFINE
    );
    ifset!(
        ctrl.dbglvl,
        METIS_DBG_MOVEINFO,
        ctrl.dbglvl -= METIS_DBG_MOVEINFO
    );

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.InitPartTmr));

    match ctrl.iptype {
        METIS_IPTYPE_RANDOM => {
            if graph.ncon == 1 {
                RandomBisection(ctrl, graph, ntpwgts, niparts);
            } else {
                McRandomBisection(ctrl, graph, ntpwgts, niparts);
            }
        }

        METIS_IPTYPE_GROW => {
            if graph.nedges == 0 {
                if graph.ncon == 1 {
                    RandomBisection(ctrl, graph, ntpwgts, niparts);
                } else {
                    McRandomBisection(ctrl, graph, ntpwgts, niparts);
                }
            } else if graph.ncon == 1 {
                GrowBisection(ctrl, graph, ntpwgts, niparts);
            } else {
                McGrowBisection(ctrl, graph, ntpwgts, niparts);
            }
        }

        _ => panic!("Unknown initial partition type: {}", ctrl.iptype),
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_IPART,
        println!("Initial Cut: {}", graph.mincut)
    );
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.InitPartTmr));
    ctrl.dbglvl = dbglvl;
}

/*************************************************************************/
/* This function computes the initial separator of the coarsest graph */
/*************************************************************************/
#[metis_func]
pub fn InitSeparator(ctrl: *mut ctrl_t, graph: *mut graph_t, niparts: idx_t) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    let mut ntpwgts = [0.5, 0.5];

    let dbglvl: mdbglvl_et = ctrl.dbglvl;
    ifset!(
        ctrl.dbglvl,
        METIS_DBG_REFINE,
        ctrl.dbglvl -= METIS_DBG_REFINE
    );
    ifset!(
        ctrl.dbglvl,
        METIS_DBG_MOVEINFO,
        ctrl.dbglvl -= METIS_DBG_MOVEINFO
    );

    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl.InitPartTmr));

    /* this is required for the cut-based part of the refinement */
    options::Setup2WayBalMultipliers(ctrl, graph, ntpwgts.as_mut_ptr());

    match ctrl.iptype {
        METIS_IPTYPE_EDGE => {
            if graph.nedges == 0 {
                RandomBisection(ctrl, graph, ntpwgts.as_mut_ptr(), niparts);
            } else {
                GrowBisection(ctrl, graph, ntpwgts.as_mut_ptr(), niparts);
            }
            refine::Compute2WayPartitionParams(ctrl, graph);
            separator::ConstructSeparator(ctrl, graph);
        }

        METIS_IPTYPE_NODE => {
            GrowBisectionNode(ctrl, graph, ntpwgts.as_mut_ptr(), niparts);
        }

        _ => panic!("Unknown iptype of {}", ctrl.iptype),
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_IPART,
        println!("Initial Sep: {}", graph.mincut)
    );
    // ifset!(ctrl.dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl.InitPartTmr));

    ctrl.dbglvl = dbglvl;
}

/*************************************************************************/
/* This function computes a bisection of a graph by randomly assigning
    the vertices followed by a bisection refinement.
    The resulting partition is returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn RandomBisection(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niparts: idx_t,
) {
    // idx_t i, ii, j, k, nvtxs, pwgts[2], zeromaxpwgt, from, me,
    //       bestcut=0, icut, mincut, inbfs;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_;
    // idx_t *perm, *bestwhere;
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    let mut bestcut = 0;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;

    refine::Allocate2WayPartitionMemory(ctrl, graph);

    // not used but declared in C function: xadj adjncy adjwgt
    get_graph_slices!(graph => vwgt);
    get_graph_slices_mut!(graph => where_);

    let mut bestwhere = vec![0; nvtxs];
    let mut perm = vec![0; nvtxs];

    let zeromaxpwgt = (*ctrl.ubfactors) * (*graph.tvwgt as real_t) * *ntpwgts;

    for inbfs in 0..niparts {
        // iset(nvtxs, 1, where_);
        where_.fill(1);

        if inbfs > 0 {
            irandArrayPermute(nvtxs as idx_t, perm.as_mut_ptr(), nvtxs as idx_t / 2, 1);
            let mut pwgts = [0, *graph.tvwgt];

            for ii in 0..nvtxs {
                let i = perm[ii] as usize;
                if pwgts[0] + vwgt[i] < zeromaxpwgt as idx_t {
                    where_[i] = 0;
                    pwgts[0] += vwgt[i];
                    pwgts[1] -= vwgt[i];
                    if pwgts[0] as f32 > zeromaxpwgt {
                        break;
                    }
                }
            }
        }

        /* Do some partition refinement  */
        refine::Compute2WayPartitionParams(ctrl, graph);
        /* println!("IPART: %3"PRIDX" [%5"PRIDX" %5"PRIDX"] [%5"PRIDX" %5"PRIDX"] %5"PRIDX"", graph.nvtxs, pwgts[0], pwgts[1], graph.pwgts[0], graph.pwgts[1], graph.mincut); */

        balance::Balance2Way(ctrl, graph, ntpwgts);
        /* println!("BPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"", graph.pwgts[0], graph.pwgts[1], graph.mincut); */

        fm::FM_2WayRefine(ctrl, graph, ntpwgts, 4);
        /* println!("RPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"", graph.pwgts[0], graph.pwgts[1], graph.mincut); */

        if inbfs == 0 || bestcut > graph.mincut {
            bestcut = graph.mincut;
            // icopy(nvtxs, where_, bestwhere);
            bestwhere.copy_from_slice(where_);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.mincut = bestcut;
    // icopy(nvtxs, bestwhere, where_);
    where_.copy_from_slice(&bestwhere);

    // WCOREPOP;
}

/*************************************************************************/
/* This function takes a graph and produces a bisection by using a region
    growing algorithm. The resulting bisection is refined using FM.
    The resulting partition is returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn GrowBisection(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();

    // idx_t i, j, k, nvtxs, drain, nleft, first, last,
    //       pwgts[2], oneminpwgt, onemaxpwgt,
    //       from, me, bestcut=0, icut, mincut, inbfs;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_;
    // idx_t *queue, *touched, *gain, *bestwhere;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    // adjwgt not used in C function
    get_graph_slices!(graph => xadj adjncy vwgt);

    refine::Allocate2WayPartitionMemory(ctrl, graph);
    get_graph_slices_mut!(graph => where_);

    let mut bestwhere = vec![0; nvtxs];
    let mut queue = vec![0; nvtxs];
    let mut touched = vec![0; nvtxs];
    let mut bestcut = 0;

    let onemaxpwgt;
    let oneminpwgt;
    {
        // I'm not sure what the actual dimensions are - it's probably the [2; real_t] declared
        // elsewhere, but I'm not quite sure so I'm isolating this for clarity.
        mkslice!(ntpwgts, 2);
        onemaxpwgt = ((*ctrl.ubfactors) * (*graph.tvwgt as f32) * ntpwgts[1]) as idx_t;
        oneminpwgt = ((1.0 / (*ctrl.ubfactors)) * (*graph.tvwgt as f32) * ntpwgts[1]) as idx_t;
    }

    for inbfs in 0..niparts {
        // iset(nvtxs, 1, where_);
        where_.fill(1);
        // iset(nvtxs, 0, touched);
        touched.fill(0);

        let mut pwgts = [0, *graph.tvwgt];

        let val = irandInRange(nvtxs as idx_t);
        queue[0] = val;
        touched[queue[0] as usize] = 1;
        let mut first = 0;
        let mut last = 1;
        let mut nleft = nvtxs - 1;
        let mut drain = false;

        /* Start the BFS from queue to get a partition */
        loop {
            if first == last {
                /* Empty. Disconnected graph! */
                if nleft == 0 || drain {
                    break;
                }
                let mut k = irandInRange(nleft as idx_t);
                let mut i = 0;
                for ii in 0..nvtxs {
                    i = ii;
                    if touched[i] == 0 {
                        if k == 0 {
                            break;
                        } else {
                            k -= 1;
                        }
                    }
                }

                queue[0] = i as idx_t;
                touched[i] = 1;
                first = 0;
                last = 1;
                nleft -= 1;
            }

            let i = queue[first] as usize;
            first += 1;
            if pwgts[0] > 0 && pwgts[1] - vwgt[i] < oneminpwgt {
                drain = true;
                continue;
            }

            where_[i] = 0;
            inc_dec!(pwgts[0], pwgts[1], vwgt[i]);
            if pwgts[1] <= onemaxpwgt {
                break;
            }
            drain = false;
            for j in xadj[i]..xadj[i + 1] {
                let k = adjncy[j as usize] as usize;
                if touched[k] == 0 {
                    queue[last] = k as idx_t;
                    last += 1;
                    touched[k] = 1;
                    nleft -= 1;
                }
            }
        }

        /* Check to see if we hit any bad limiting cases */
        if pwgts[1] == 0 {
            where_[irandInRange(nvtxs as idx_t) as usize] = 1;
        }
        if pwgts[0] == 0 {
            where_[irandInRange(nvtxs as idx_t) as usize] = 0;
        }
        /*************************************************************
         * Do some partition refinement
         **************************************************************/
        refine::Compute2WayPartitionParams(ctrl, graph);
        /*
        println!("IPART: %3"PRIDX" [%5"PRIDX" %5"PRIDX"] [%5"PRIDX" %5"PRIDX"] %5"PRIDX"",
            graph.nvtxs, pwgts[0], pwgts[1], graph.pwgts[0], graph.pwgts[1], graph.mincut);
        */

        balance::Balance2Way(ctrl, graph, ntpwgts);
        /*
        println!("BPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"", graph.pwgts[0],
            graph.pwgts[1], graph.mincut);
        */

        fm::FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
        /*
        println!("RPART: [%5"PRIDX" %5"PRIDX"] %5"PRIDX"", graph.pwgts[0],
            graph.pwgts[1], graph.mincut);
        */

        if inbfs == 0 || bestcut > graph.mincut {
            bestcut = graph.mincut;
            // icopy(nvtxs, where_, bestwhere);
            bestwhere.copy_from_slice(where_);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.mincut = bestcut;
    // icopy(nvtxs, bestwhere, where_);
    where_.copy_from_slice(&bestwhere);

    // WCOREPOP;
}

/*************************************************************************/
/* This function takes a multi-constraint graph and computes a bisection
    by randomly assigning the vertices and then refining it. The resulting
    partition is returned in graph.where_.
*/
/**************************************************************************/
#[metis_func]
pub fn McRandomBisection(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niparts: idx_t,
) {
    // issues with niparts being corrupted
    let init_niparts = Box::new(niparts);

    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    // idx_t i, ii, j, k, nvtxs, ncon, from, bestcut=0, mincut, inbfs, qnum;
    // idx_t *bestwhere, *where_, *perm, *counts;
    // idx_t *vwgt;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    get_graph_slices!(graph => vwgt);

    refine::Allocate2WayPartitionMemory(ctrl, graph);
    get_graph_slices_mut!(graph => where_);

    let mut bestwhere = vec![0; nvtxs];
    let mut perm = vec![0; nvtxs];
    let mut counts = vec![0; ncon as usize];
    let mut bestcut = 0;

    for inbfs in 0..(2 * niparts) {
        debug_assert_eq!(niparts, *init_niparts);
        irandArrayPermute(nvtxs as idx_t, perm.as_mut_ptr(), nvtxs as idx_t / 2, 1);
        // iset(ncon, 0, counts);
        counts.fill(0);

        /* partition by splitting the queues randomly */
        for ii in 0..nvtxs {
            let i = perm[ii] as usize;
            let idx = i * ncon;

            // qnum = iargmax(ncon, vwgt+i*ncon,1);
            let qnum = util::iargmax(&vwgt[idx..(idx + ncon as usize)], 1) as usize;
            where_[i as usize] = (counts[qnum]) % 2;
            counts[qnum] += 1;
        }

        refine::Compute2WayPartitionParams(ctrl, graph);

        fm::FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
        balance::Balance2Way(ctrl, graph, ntpwgts);
        fm::FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
        balance::Balance2Way(ctrl, graph, ntpwgts);
        fm::FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);

        if inbfs == 0 || bestcut >= graph.mincut {
            bestcut = graph.mincut;
            // icopy(nvtxs, where_, bestwhere);
            bestwhere.copy_from_slice(where_);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.mincut = bestcut;
    // icopy(nvtxs, bestwhere, where_);
    where_.copy_from_slice(&bestwhere);

    // WCOREPOP;
}

/*************************************************************************/
/* This function takes a multi-constraint graph and produces a bisection
    by using a region growing algorithm. The resulting partition is
    returned in graph.where_.
*/
/*************************************************************************/
#[metis_func]
pub fn McGrowBisection(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niparts: idx_t,
) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    // idx_t i, j, k, nvtxs, ncon, from, bestcut=0, mincut, inbfs;
    // idx_t *bestwhere, *where_;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;

    refine::Allocate2WayPartitionMemory(ctrl, graph);
    get_graph_slices_mut!(graph => where_);

    let mut bestwhere = vec![0; nvtxs];
    let mut bestcut = 0;

    for inbfs in 0..2 * niparts {
        // iset(nvtxs, 1, where_);
        where_.fill(1);
        where_[irandInRange(nvtxs as idx_t) as usize] = 0;

        refine::Compute2WayPartitionParams(ctrl, graph);

        balance::Balance2Way(ctrl, graph, ntpwgts);
        fm::FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
        balance::Balance2Way(ctrl, graph, ntpwgts);
        fm::FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);

        if inbfs == 0 || bestcut >= graph.mincut {
            bestcut = graph.mincut;
            // icopy(nvtxs, where_, bestwhere);
            bestwhere.copy_from_slice(where_);
            if bestcut == 0 {
                break;
            }
        }
    }

    graph.mincut = bestcut;
    // icopy(nvtxs, bestwhere, where_);
    where_.copy_from_slice(&bestwhere);

    // WCOREPOP;
}

/*************************************************************************/
/* This function takes a graph and produces a tri-section into left, right,
   and separator using a region growing algorithm. The resulting separator
   is refined using node FM.
   The resulting partition is returned in graph.where_.
*/
/**************************************************************************/
#[metis_func]
pub fn GrowBisectionNode(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    ntpwgts: *mut real_t,
    niparts: idx_t,
) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    // idx_t i, j, k, nvtxs, drain, nleft, first, last, pwgts[2], oneminpwgt,
    //       onemaxpwgt, from, me, bestcut=0, icut, mincut, inbfs;
    // idx_t *xadj, *vwgt, *adjncy, *adjwgt, *where_, *bndind;
    // idx_t *queue, *touched, *gain, *bestwhere;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs;
    // adjwgt not used in C function
    get_graph_slices!(graph => xadj vwgt adjncy);

    let mut bestwhere = vec![0; nvtxs as usize];
    let mut queue = vec![0; nvtxs as usize];
    let mut touched = vec![0; nvtxs as usize];
    let mut bestcut = 0;

    let onemaxpwgt = (*ctrl.ubfactors) * (*graph.tvwgt as f32) * 0.5;
    let oneminpwgt = (1.0 / (*ctrl.ubfactors)) * (*graph.tvwgt as f32) * 0.5;

    /* Allocate refinement memory. Allocate sufficient memory for both edge and node */
    graph.pwgts = imalloc(3, c"GrowBisectionNode: pwgts".as_ptr()) as _;
    graph.where_ = imalloc(nvtxs as usize, c"GrowBisectionNode: where_".as_ptr()) as _;
    graph.bndptr = imalloc(nvtxs as usize, c"GrowBisectionNode: bndptr".as_ptr()) as _;
    graph.bndind = imalloc(nvtxs as usize, c"GrowBisectionNode: bndind".as_ptr()) as _;
    graph.id = imalloc(nvtxs as usize, c"GrowBisectionNode: id".as_ptr()) as _;
    graph.ed = imalloc(nvtxs as usize, c"GrowBisectionNode: ed".as_ptr()) as _;
    graph.nrinfo = gk_malloc(
        nvtxs as usize * std::mem::size_of::<nrinfo_t>(),
        c"GrowBisectionNode: nrinfo".as_ptr(),
    ) as _;

    get_graph_slices_mut!(graph => where_ bndind);

    for inbfs in 0..niparts {
        // iset(nvtxs, 1, where_);
        where_.fill(1);
        // iset(nvtxs, 0, touched);
        touched.fill(0);

        let mut pwgts = [0, *graph.tvwgt];

        queue[0] = irandInRange(nvtxs);
        touched[queue[0] as usize] = 1;
        let mut first = 0;
        let mut last = 1;
        let mut nleft = nvtxs - 1;
        let mut drain = 0;

        /* Start the BFS from queue to get a partition */
        loop {
            if first == last {
                /* Empty. Disconnected graph! */
                if nleft == 0 || drain != 0 {
                    break;
                }
                let mut k = irandInRange(nleft);
                let mut i = 0;
                for ii in 0..nvtxs {
                    i = ii;
                    /* select the kth untouched vertex */
                    if touched[i as usize] == 0 {
                        if k == 0 {
                            break;
                        } else {
                            k -= 1;
                        }
                    }
                }

                queue[0] = i;
                touched[i as usize] = 1;
                first = 0;
                last = 1;
                nleft -= 1;
            }

            let i = queue[first] as usize;
            first += 1;
            if pwgts[1] as f32 - (vwgt[i] as f32) < oneminpwgt {
                drain = 1;
                continue;
            }

            where_[i] = 0;
            inc_dec!(pwgts[0], pwgts[1], vwgt[i]);
            if pwgts[1] as real_t <= onemaxpwgt {
                break;
            }
            drain = 0;
            for j in xadj[i]..xadj[i + 1] {
                let k = adjncy[j as usize] as usize;
                if touched[k] == 0 {
                    queue[last] = k as idx_t;
                    last += 1;
                    touched[k] = 1;
                    nleft -= 1;
                }
            }
        }

        /*************************************************************
         * Do some partition refinement
         **************************************************************/
        refine::Compute2WayPartitionParams(ctrl, graph);
        balance::Balance2Way(ctrl, graph, ntpwgts);
        fm::FM_2WayRefine(ctrl, graph, ntpwgts, 4);

        /* Construct and refine the vertex separator */
        for i in 0..graph.nbnd {
            let j = bndind[i as usize] as usize;
            if xadj[j + 1] - xadj[j] > 0 {
                /* ignore islands */
                where_[j] = 2;
            }
        }

        srefine::Compute2WayNodePartitionParams(ctrl, graph);
        sfm::FM_2WayNodeRefine2Sided(ctrl, graph, 1);
        sfm::FM_2WayNodeRefine1Sided(ctrl, graph, 4);

        /*
        println!("ISep: [{} {} {} {}] {}",
            inbfs, graph.pwgts[0], graph.pwgts[1], graph.pwgts[2], bestcut);
        */

        if inbfs == 0 || bestcut > graph.mincut {
            bestcut = graph.mincut;
            // icopy(nvtxs, where_, bestwhere);
            bestwhere.copy_from_slice(where_);
        }
    }

    graph.mincut = bestcut;
    // icopy(nvtxs, bestwhere, where_);
    where_.copy_from_slice(&bestwhere);

    // WCOREPOP;
}

/*************************************************************************/
/* This function takes a graph and produces a tri-section into left, right,
   and separator using a region growing algorithm. The resulting separator
   is refined using node FM.
   The resulting partition is returned in graph.where_.
*/
/**************************************************************************/
//     #[metis_func]
// pub fn GrowBisectionNode2(ctrl: *mut ctrl_t, graph: *mut graph_t, ntpwgts: *mut real_t, niparts: idx_t)
// {
//   // idx_t i, j, k, nvtxs, bestcut=0, mincut, inbfs;
//   // idx_t *xadj, *where_, *bndind, *bestwhere;
//
//   // WCOREPUSH;
//
//   let nvtxs  = graph.nvtxs;
//   xadj   = graph.xadj;
//
//   /* Allocate refinement memory. Allocate sufficient memory for both edge and node */
//   graph.pwgts  = imalloc(3, "GrowBisectionNode: pwgts");
//   graph.where_  = imalloc(nvtxs, "GrowBisectionNode: where_");
//   graph.bndptr = imalloc(nvtxs, "GrowBisectionNode: bndptr");
//   graph.bndind = imalloc(nvtxs, "GrowBisectionNode: bndind");
//   graph.id     = imalloc(nvtxs, "GrowBisectionNode: id");
//   graph.ed     = imalloc(nvtxs, "GrowBisectionNode: ed");
//   graph.nrinfo = gk_malloc(nvtxs*sizeof(nrinfo_t), "GrowBisectionNode: nrinfo");
//
//   bestwhere = vec![0; nvtxs as usize];
//
//   where_  = graph.where_;
//   bndind = graph.bndind;
//
//   for inbfs in 0..niparts {
//     iset(nvtxs, 1, where_);
//     if inbfs > 0
//         {
//             where_[irandInRange(nvtxs)] = 0;
//         }
//     Compute2WayPartitionParams(ctrl, graph);
//     General2WayBalance(ctrl, graph, ntpwgts);
//     FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);
//
//     /* Construct and refine the vertex separator */
//     for i in 0..graph.nbnd {
//       j = bndind[i];
//       if xadj[j+1]-xadj[j] > 0 /* ignore islands */
//             {
//                 where_[j] = 2;
//             }
//         }
//
//     Compute2WayNodePartitionParams(ctrl, graph);
//     FM_2WayNodeRefine2Sided(ctrl, graph, 4);
//
//     /*
//     println!("ISep: [{} {} {} {}] {}",
//         inbfs, graph.pwgts[0], graph.pwgts[1], graph.pwgts[2], bestcut);
//     */
//
//     if inbfs == 0 || bestcut > graph.mincut {
//       bestcut = graph.mincut;
//       icopy(nvtxs, where_, bestwhere);
//     }
//   }
//
//   graph.mincut = bestcut;
//   icopy(nvtxs, bestwhere, where_);
//
//   // WCOREPOP;
// }


#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use crate::tests::ab_test_partition_test_graphs;

    use super::*;

    #[test]
    fn ab_Init2WayPartition() {
        ab_test_partition_test_graphs("Init2WayPartition:rs", Optype::Pmetis, 8, 2, |mut g| {
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }
    
    #[test]
    fn ab_InitSeparator() {
        ab_test_partition_test_graphs("InitSeparator:rs", Optype::Ometis, 3, 1, |mut g| {
            g.random_vwgt();
            g
        });
    }

    #[test]
    fn ab_RandomBisection() {
        ab_test_partition_test_graphs("RandomBisection:rs", Optype::Pmetis, 8, 1, |mut g| {
            g.set_initial_part_strategy(Iptype::Random);
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    fn ab_McRandomBisection() {
        ab_test_partition_test_graphs("McRandomBisection:rs", Optype::Pmetis, 8, 3, |mut g| {
            g.set_initial_part_strategy(Iptype::Random);
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    fn ab_GrowBisection() {
        ab_test_partition_test_graphs("GrowBisection:rs", Optype::Pmetis, 8, 1, |mut g| {
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    fn ab_McGrowBisection() {
        ab_test_partition_test_graphs("McGrowBisection:rs", Optype::Pmetis, 8, 3, |mut g| {
            g.random_adjwgt();
            g.random_tpwgts();
            g
        });
    }

    #[test]
    fn ab_GrowBisectionNode() {
        ab_test_partition_test_graphs("GrowBisectionNode:rs", Optype::Ometis, 3, 1, |mut g| {
            g.set_initial_part_strategy(Iptype::Node);
            g
        });
    }
}
