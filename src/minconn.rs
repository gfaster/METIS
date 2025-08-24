/*
\file
\brief Functions that deal with prunning the number of adjacent subdomains in kmetis

\date Started 7/15/98
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota
\version $Id: minconn.c 17513 2014-08-05 16:20:50Z dominique $
*/

use std::ptr::NonNull;

use crate::*;

/// The original is a gklib routine. Since it's only used here, I'm not putting it in utils (this
/// is wrong).
///
/// I think I may be able to make this an associated function with DAL
// TODO: move to utils
pub fn iarray2csr(n: usize, range: usize, array: &[idx_t], ptr: &mut [idx_t], ind: &mut [idx_t]) {
    assert_eq!(range + 1, ptr.len(), "assumption about preconditions");
    assert_eq!(n, array.len(), "assumption about preconditions");
    ptr.fill(0);

    for &a in array {
        ptr[a as usize] += 1;
    }

    util::make_csr(range, ptr);
    for i in 0..n {
        ind[ptr[array[i] as usize] as usize] = i as idx_t;
        ptr[array[i] as usize] += 1;
    }
    util::shift_csr(range, ptr);
}

/*************************************************************************/
/* This function computes the subdomain graph storing the result in the
pre-allocated worspace arrays */
/*************************************************************************/
#[metis_func]
pub extern "C" fn ComputeSubDomainGraph(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_ref().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, pid, other, nparts, nvtxs, nnbrs;
    // idx_t *xadj, *adjncy, *adjwgt, *where_;
    // idx_t *pptr, *pind;
    // idx_t nads=0, *vadids, *vadwgts;

    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    // xadj   = graph.xadj;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // where_  = graph.where_;
    get_graph_slices!(graph => where_);

    let nparts = ctrl.nparts as usize;

    mkslice_mut!(vadids: ctrl->pvec1, nparts + 1);
    // vadwgts = iset(nparts, 0, ctrl.pvec2);
    mkslice_mut!(vadwgts: ctrl->pvec2, nparts + 1);
    vadwgts.fill(0);

    let mut pptr = vec![0; nparts + 1 as usize];
    let mut pind = vec![0; nvtxs as usize];
    iarray2csr(nvtxs, nparts, where_, &mut pptr, &mut pind);

    for pid in (0)..(nparts) {
        let mut nads = 0;
        match ctrl.objtype {
            METIS_OBJTYPE_CUT => {
                // ckrinfo_t *rinfo;
                // cnbr_t *nbrs;

                let rinfo = get_graph_slice!(graph => ckrinfo);
                // rinfo = graph.ckrinfo;

                // nads=0;
                for ii in (pptr[pid as usize])..(pptr[(pid + 1) as usize]) {
                    let i = pind[ii as usize] as usize;
                    debug_assert_eq!(pid, where_[i as usize] as usize);

                    if rinfo[i as usize].ed > 0 {
                        let nnbrs = rinfo[i as usize].nnbrs as usize;
                        let nbrs = std::slice::from_raw_parts(
                            ctrl.cnbrpool.add(rinfo[i as usize].inbr as usize),
                            nnbrs,
                        );

                        for j in (0)..(nnbrs) {
                            let other = nbrs[j as usize].pid;
                            if vadwgts[other as usize] == 0 {
                                vadids[(nads) as usize] = other;
                                nads += 1;
                            }
                            vadwgts[other as usize] += nbrs[j as usize].ed;
                        }
                    }
                }
            }

            METIS_OBJTYPE_VOL => {
                // vkrinfo_t *rinfo;
                // vnbr_t *nbrs;

                let rinfo = get_graph_slice!(graph => vkrinfo);
                for ii in (pptr[pid as usize])..(pptr[(pid + 1) as usize]) {
                    let i = pind[ii as usize];
                    debug_assert_eq!(pid as idx_t, where_[i as usize]);

                    if rinfo[i as usize].ned > 0 {
                        let nnbrs = rinfo[i as usize].nnbrs as usize;
                        let nbrs = std::slice::from_raw_parts(
                            ctrl.vnbrpool.add(rinfo[i as usize].inbr as usize),
                            nnbrs,
                        );

                        for j in (0)..(nnbrs) {
                            let other = nbrs[j as usize].pid;
                            if vadwgts[other as usize] == 0 {
                                vadids[(nads) as usize] = other;
                                nads += 1;
                            }
                            vadwgts[other as usize] += nbrs[j as usize].ned;
                        }
                    }
                }
            }

            _ => panic!("Unknown objtype: {}", ctrl.objtype),
        }

        /* See if you have enough memory to store the adjacent info for that subdomain */
        {
            mkslice_mut!(ctrl->maxnads, nparts);
            if maxnads[pid as usize] < nads {
                maxnads[pid as usize] = 2 * nads;

                {
                    mkslice_mut!(ctrl->adids, nparts);
                    let old_adids = NonNull::slice_from_raw_parts(
                        NonNull::new(adids[pid]).unwrap(),
                        nads as usize,
                    );
                    adids[pid as usize] = gk::irealloc(
                        old_adids,
                        maxnads[pid] as usize,
                        c"ComputeSubDomainGraph: adids[pid]",
                    )
                    .as_buf_ptr();
                }

                {
                    mkslice_mut!(ctrl->adwgts, nparts);
                    let old_adwgts = NonNull::slice_from_raw_parts(
                        NonNull::new(adwgts[pid]).unwrap(),
                        nads as usize,
                    );
                    adwgts[pid as usize] = gk::irealloc(
                        old_adwgts,
                        maxnads[pid] as usize,
                        c"ComputeSubDomainGraph: adids[pid]",
                    )
                    .as_buf_ptr();
                }
            }

            // have to remake the slice since it's potentially realloc-d
            mkslice_mut!(cnads: ctrl->nads, nparts);
            mkslice_mut!(ctrl->adids, nparts);
            mkslice_mut!(ctrl->adwgts, nparts);
            cnads[pid] = nads;
            for j in (0)..(nads as usize) {
                *adids[pid].add(j) = vadids[j];
                *adwgts[pid].add(j) = vadwgts[vadids[j] as usize];

                vadwgts[vadids[j] as usize] = 0;
            }
        }
    }

    // WCOREPOP;
}

/*************************************************************************/
/* This function updates the weight of an edge in the subdomain graph by
adding to it the value of ewgt. The update can either increase or
decrease the weight of the subdomain edge based on the value of ewgt.

\param u is the ID of one of the incident subdomains to the edge
\param v is the ID of the other incident subdomains to the edge
\param ewgt is the weight to be added to the subdomain edge
\param nparts is the number of subdomains
\param r_maxndoms is the maximum number of adjacent subdomains and is
updated as necessary. The update is skipped if a std::ptr::null_mut() value is
supplied.
*/
/*************************************************************************/
#[metis_func]
pub extern "C" fn UpdateEdgeSubDomainGraph(
    ctrl: *mut ctrl_t,
    u: idx_t,
    v: idx_t,
    ewgt: idx_t,
    r_maxndoms: *mut idx_t,
) {
    if ewgt == 0 {
        return;
    }

    let ctrl = ctrl.as_mut().unwrap();
    let mut u = u as usize;
    let mut v = v as usize;

    // idx_t i, j, nads;
    let nparts = ctrl.nparts as usize;
    mkslice_mut!(ctrl->adids, nparts);
    mkslice_mut!(ctrl->adwgts, nparts);
    mkslice_mut!(cnads: ctrl->nads, nparts);
    mkslice_mut!(ctrl->maxnads, nparts);

    for _ in (0)..(2) {
        let mut nads = cnads[u] as usize;
        mkslice_mut!(adids_u: adids[u], nads);
        let mut adids_u = adids_u;
        mkslice_mut!(adwgts_u: adwgts[u], nads);
        let mut adwgts_u = adwgts_u;

        /* Find the edge */
        let j = adids_u.iter().position(|&adid| adid == v as idx_t);

        if let Some(j) = j {
            /* See if the updated edge becomes 0 */
            adwgts_u[j] += ewgt;
            debug_assert!(adwgts_u[j] >= 0);
            if adwgts_u[j] == 0 {
                // this "pops", but adids_u and adwgts_u are not read again before being recreated
                adids_u[j] = adids_u[nads - 1];
                adwgts_u[j] = adwgts_u[nads - 1];
                nads -= 1;
                if !r_maxndoms.is_null() && (nads + 1) as idx_t == *r_maxndoms {
                    *r_maxndoms = cnads[(util::iargmax(cnads, 1)) as usize];
                }
            }
        } else {
            /* Deal with the case in which the edge was not found */
            debug_assert!(ewgt > 0);
            if maxnads[u] == nads as idx_t {
                maxnads[u] = 2 * (nads + 1) as idx_t;
                // adids[u]   = irealloc(ctrl.adids[u], ctrl.maxnads[u], "IncreaseEdgeSubDomainGraph: adids[pid as usize]");
                // adwgts[u]  = irealloc(ctrl.adwgts[u], ctrl.maxnads[u], "IncreaseEdgeSubDomainGraph: adids[pid as usize]");

                let old_adids =
                    NonNull::slice_from_raw_parts(NonNull::new(adids[u]).unwrap(), nads);
                let new_adids = gk::irealloc(
                    old_adids,
                    maxnads[u] as usize,
                    c"IncreaseEdgeSubDomainGraph: adids[pid]",
                );
                adids[u] = new_adids.as_buf_ptr();

                let old_adwgts =
                    NonNull::slice_from_raw_parts(NonNull::new(adwgts[u]).unwrap(), nads);
                let new_adwgts = gk::irealloc(
                    old_adwgts,
                    maxnads[u] as usize,
                    c"IncreaseEdgeSubDomainGraph: adwgts[pid]",
                );
                adwgts[u] = new_adwgts.as_buf_ptr();
            }
            adids_u = std::slice::from_raw_parts_mut(adids[u], nads + 1);
            adids_u[nads] = v as idx_t;
            adwgts_u = std::slice::from_raw_parts_mut(adwgts[u], nads + 1);
            adwgts_u[nads] = ewgt;
            nads += 1;
            if !r_maxndoms.is_null() && nads as idx_t > *r_maxndoms {
                println!(
                    "You just increased the maxndoms: {:} {:}",
                    nads, *r_maxndoms
                );
                *r_maxndoms = nads as idx_t;
            }
        }
        cnads[u] = nads as idx_t;

        std::mem::swap(&mut u, &mut v);
        // SWAP(u, v, j);
    }
}

/*************************************************************************/
/* This function computes the subdomain graph */
/*************************************************************************/
#[metis_func]
pub extern "C" fn EliminateSubDomainEdges(ctrl: *mut ctrl_t, graph: *mut graph_t) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, k, ncon, nparts, scheme, pid_from, pid_to, me, other, nvtxs,
    // total, max, avg, totalout, nind=0, ncand=0, ncand2, target, target2,
    // nadd, bestnadd=0;
    // idx_t min, move, *cpwgt;
    // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *pwgts, *where_, *maxpwgt,
    // *mypmat, *otherpmat, *kpmat, *ind;
    // idx_t *nads, **adids, **adwgts;
    // ikv_t *cand, *cand2;
    // ipq_t queue;
    // real_t *tpwgts, badfactor=1.4;
    // idx_t *pptr, *pind;
    // idx_t *vmarker=std::ptr::null_mut(), *pmarker=std::ptr::null_mut(), *modind=std::ptr::null_mut();  /* volume specific work arrays */
    // WCOREPUSH;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;

    let badfactor = 1.4;

    // xadj   = graph.xadj;
    // adjncy = graph.adjncy;
    // vwgt   = graph.vwgt;
    // adjwgt = ( if (ctrl.objtype == METIS_OBJTYPE_VOL)  { (std::ptr::null_mut()) } else { (graph.adjwgt) });
    // where_ = graph.where_;
    // pwgts = graph.pwgts;  /* We assume that this is properly initialized */
    get_graph_slices!(ctrl, graph => xadj adjncy vwgt where_ pwgts tvwgt);
    let adjwgt = (ctrl.objtype != METIS_OBJTYPE_VOL).then(|| get_graph_slice!(graph => adjwgt));

    let nparts = ctrl.nparts as usize;

    let mut cpwgt: Vec<idx_t> = vec![0; ncon];
    let mut maxpwgt: Vec<idx_t> = vec![0; nparts * ncon];
    // let mut ind: Vec<idx_t> = vec![0; nvtxs];
    let mut ind: Vec<idx_t> = Vec::with_capacity(nvtxs);
    let mut otherpmat: Vec<idx_t> = vec![0; nparts];

    // cand  = ikvwspacemalloc(ctrl, nparts);
    // cand2 = ikvwspacemalloc(ctrl, nparts);
    let mut cand: Vec<ikv_t> = Vec::with_capacity(nparts);
    let mut cand2: Vec<ikv_t> = Vec::with_capacity(nparts);

    let mut pptr = vec![0; nparts + 1 as usize];
    let mut pind = vec![0; nvtxs as usize];
    iarray2csr(nvtxs, nparts, where_, &mut pptr, &mut pind);

    let mut modind: Option<Vec<idx_t>> = None;
    let mut vmarker: Option<Vec<idx_t>> = None;
    let mut pmarker: Option<Vec<idx_t>> = None;
    if ctrl.objtype == METIS_OBJTYPE_VOL {
        /* Vol-refinement specific working arrays */
        modind = Some(vec![0; nvtxs as usize]);
        vmarker = Some(vec![0; nvtxs as usize]);
        pmarker = Some(vec![-1; nparts as usize]);
    }

    /* Compute the pmat matrix and ndoms */
    ComputeSubDomainGraph(ctrl, graph);

    // let adids  = ctrl.adids;
    // let adwgts = ctrl.adwgts;
    mkslice!(ctrl->nads, nparts);
    mkslice!(ctrl->adids, nparts);
    mkslice!(ctrl->adwgts, nparts);
    mkslice!(ctrl->tpwgts, nparts);
    mkslice!(ctrl->ubfactors, ncon);

    // mypmat = iset(nparts, 0, ctrl.pvec1);
    // kpmat  = iset(nparts, 0, ctrl.pvec2);
    mkslice_mut!(mypmat: ctrl->pvec1, nparts);
    mkslice_mut!(kpmat: ctrl->pvec2, nparts);
    mypmat.fill(0);
    kpmat.fill(0);

    /* Compute the maximum allowed weight for each domain */
    for i in (0)..(nparts) {
        for j in (0)..(ncon) {
            maxpwgt[(i * ncon + j) as usize] = ((if ncon == 1 { 1.25 } else { 1.025 })
                * tpwgts[i as usize]
                * (tvwgt[j as usize] as real_t)
                * (ubfactors[j as usize])) as idx_t;
        }
    }

    // ipqInit(&queue, nparts);
    let mut queue = pqueue::IPQueue::new(nparts);
    let mut bestnadd = 0;

    /* Get into the loop eliminating subdomain connections */
    loop {
        let total: idx_t = nads.iter().sum();
        let avg = total / (nparts as idx_t);
        // let max = nads[util::iargmax(nads, 1)];
        let max = nads.iter().copied().max().unwrap();

        ifset!(
            ctrl.dbglvl,
            METIS_DBG_CONNINFO,
            println!(
                "Adjacent Subdomain Stats: Total: {:3}, Max: {:3}[{}], Avg: {:3}",
                total,
                max,
                util::iargmax(nads, 1),
                avg
            )
        );

        if (max as real_t) < (badfactor * avg as real_t) {
            break;
        }

        /* Add the subdomains that you will try to reduce their connectivity */
        queue.reset();
        for i in (0)..(nparts) {
            if nads[i] >= avg + (max - avg) / 2 {
                // ipqInsert(&queue, i, nads[i as usize]);
                queue.insert(i as idx_t, nads[i as usize]);
            }
        }

        let mut move_ = 0;
        while let Some(me) = queue.pop() {
            let me = me as usize;
            let adwgts_me = adwgts[me];
            mkslice!(adwgts_me, nads[me]);
            let adids_me = adids[me];
            mkslice!(adids_me, nads[me]);

            // totalout = isum(nads[me], adwgts_me, 1);
            let totalout: idx_t = adwgts_me.iter().sum::<idx_t>();
            // let mut ncand2 = 0;
            cand2.clear();
            for i in 0..(nads[me] as usize) {
                mypmat[adids_me[i as usize] as usize] = adwgts_me[i as usize];

                /* keep track of the weakly connected adjacent subdomains */
                if 2 * nads[me] * adwgts_me[i as usize] < totalout {
                    // cand2[ncand2 as usize].val   = adids_me[i as usize];
                    // cand2[(ncand2) as usize].key = adwgts_me[i as usize];
                    // ncand2+=1;
                    cand2.push(ikv_t {
                        val: adids_me[i],
                        key: adwgts_me[i],
                    })
                }
            }

            ifset!(
                ctrl.dbglvl,
                METIS_DBG_CONNINFO,
                print!(
                    "Me: {:}, Degree: {:4}, TotalOut: {:},\n",
                    me, nads[me], totalout
                )
            );

            /* Sort the connections according to their cut */
            // ikvsorti(ncand2, cand2);
            ikvsorti(&mut cand2);

            /* Two schemes are used for eliminating subdomain edges.
            The first, tries to eliminate subdomain edges by moving remote groups
            of vertices to subdomains that 'me' is already connected to.
            The second, tries to eliminate subdomain edges by moving entire sets of
            my vertices that connect to the 'other' subdomain to a subdomain that
            I'm already connected to.
            These two schemes are applied in sequence. */
            let mut target = usize::MAX;
            let mut target2 = usize::MAX;
            // let mut nind = 0;
            for scheme in [0, 1] {
                for min in 0..cand2.len() {
                    let other = cand2[min].val as usize;

                    /* pid_from is the subdomain from where_ the vertices will be removed.
                    pid_to is the adjacent subdomain to pid_from that defines the
                    (me, other) subdomain edge that needs to be removed */
                    let pid_from;
                    let pid_to;
                    if scheme == 0 {
                        pid_from = other;
                        pid_to = me;
                    } else {
                        pid_from = me;
                        pid_to = other;
                    }

                    /* Go and find the vertices in 'other' that are connected in 'me' */
                    // nind = 0;
                    ind.clear();
                    for ii in pptr[pid_from]..pptr[pid_from + 1] {
                        let i = pind[ii as usize] as usize;
                        debug_assert_eq!(where_[i], pid_from as idx_t);
                        for j in xadj[i]..xadj[i + 1] {
                            if where_[adjncy[j as usize] as usize] == pid_to as idx_t {
                                ind.push(i as idx_t);
                                // ind[nind] = i as idx_t;
                                // nind += 1;
                                break;
                            }
                        }
                    }

                    /* Go and construct the otherpmat to see where these nind vertices are
                    connected to */
                    // iset(ncon, 0, cpwgt);
                    cpwgt.fill(0);
                    cand.clear();
                    // ncand=0;
                    for &i in &ind {
                        let i = i as usize;
                        // blas::iaxpy(ncon, 1, vwgt+i*ncon, 1, cpwgt, 1);
                        blas::iaxpy(ncon, 1, &vwgt[cntrng!(i * ncon, ncon)], 1, &mut cpwgt, 1);

                        for j in xadj[i]..xadj[i + 1] {
                            let k = where_[adjncy[j as usize] as usize] as usize;
                            if k == pid_from {
                                continue;
                            }
                            if otherpmat[k] == 0 {
                                cand.push(ikv_t {
                                    val: k as idx_t,
                                    key: -1,
                                })
                                // cand[(ncand) as usize].val = k;
                                // ncand+=1;
                            }
                            // otherpmat[k as usize] += (adjwgt ? adjwgt[j as usize] : 1);
                            otherpmat[k] += adjwgt.map_or(1, |adjwgt| adjwgt[j as usize]);
                        }
                    }

                    for cand in &mut cand {
                        cand.key = otherpmat[cand.val as usize];
                        debug_assert!(cand.key > 0);
                    }

                    // ikvsortd(ncand, cand);
                    ikvsortd(&mut cand);

                    ifset!(
                        ctrl.dbglvl,
                        METIS_DBG_CONNINFO,
                        println!(
                            "\tMinOut: {:4}, to: {:3}, TtlWgt: {:5}[#:{:}]",
                            mypmat[other as usize],
                            other,
                            cpwgt.iter().sum::<idx_t>(),
                            ind.len()
                        )
                    );

                    /* Go through and select the first domain that is common with 'me', and does
                    not increase the nads[target] higher than nads[me], subject to the maxpwgt
                    constraint. Traversal is done from the mostly connected to the least. */
                    for i in &cand {
                        let k = i.val as usize;

                        if mypmat[k] > 0 {
                            /* Check if balance will go off */
                            // if (!mcutil::ivecaxpylez(ncon, 1, cpwgt, pwgts+k*ncon, maxpwgt+k*ncon))
                            if !util::ivecaxpylez(
                                1,
                                &cpwgt,
                                &pwgts[cntrng!(k * ncon, ncon)],
                                &maxpwgt[cntrng!(k * ncon, ncon)],
                            ) {
                                continue;
                            }

                            /* get a dense vector out of k's connectivity */
                            for j in 0..nads[k] {
                                let j = j as usize;
                                mkslice!(adids_k: adids[k], nads[k]);
                                mkslice!(adwgts_k: adwgts[k], nads[k]);
                                kpmat[adids_k[j] as usize] = adwgts_k[j];
                            }

                            /* Check if the move to domain k will increase the nads of another
                            subdomain j that the set of vertices being moved are connected
                            to but domain k is not connected to. */
                            // for j in (0)..(nparts) {
                            //   if (otherpmat[j] > 0 && kpmat[j] == 0 && nads[j]+1 >= nads[me])  {
                            //     break;
                            //   }
                            // }
                            let bad_effects = (0..nparts).any(|j| {
                                otherpmat[j] > 0 && kpmat[j] == 0 && nads[j] + 1 >= nads[me]
                            });

                            /* There were no bad second level effects. See if you can find a
                            subdomain to move to. */
                            if !bad_effects {
                                let mut nadd = 0;
                                for j in (0)..(nparts) {
                                    if otherpmat[j] > 0 && kpmat[j] == 0 {
                                        nadd += 1;
                                    }
                                }

                                ifset!(
                                    ctrl.dbglvl,
                                    METIS_DBG_CONNINFO,
                                    print!(
                                        "\t\tto={:}, nadd={:}, {:}\n",
                                        k, nadd, nads[k as usize]
                                    )
                                );

                                if nads[k as usize] + nadd < nads[me]
                                    && (target2 == usize::MAX
                                        || nads[target2] + bestnadd > nads[k] + nadd
                                        || (nads[target2] + bestnadd == nads[k] + nadd
                                            && bestnadd > nadd))
                                    {
                                        target2 = k;
                                        bestnadd = nadd;
                                    }

                                if nadd == 0 {
                                    target = k;
                                }
                            }

                            /* reset kpmat for the next iteration */
                            for j in (0)..(nads[k]) {
                                mkslice!(adids_k: adids[k], nads[k]);
                                kpmat[adids_k[j as usize] as usize] = 0;
                            }
                        }

                        if target != usize::MAX {
                            break;
                        }
                    }

                    /* reset the otherpmat for the next iteration */
                    for cand in &cand {
                        otherpmat[cand.val as usize] = 0;
                    }

                    if target == usize::MAX && target2 != usize::MAX {
                        target = target2;
                    }

                    if target != usize::MAX {
                        ifset!(
                            ctrl.dbglvl,
                            METIS_DBG_CONNINFO,
                            print!("\t\tScheme: {:}. Moving to {:}\n", scheme, target)
                        );
                        move_ = 1;
                        break;
                    }
                }

                if target != usize::MAX {
                    break;
                } /* A move was found. No need to try the other scheme */
            }

            /* reset the mypmat for next iteration */
            for i in (0)..(nads[me]) {
                mypmat[adids_me[i as usize] as usize] = 0;
            }

            /* Note that once a target is found the above loops exit right away. So the
            following variables are valid */
            if target != usize::MAX {
                match ctrl.objtype {
                    METIS_OBJTYPE_CUT => MoveGroupMinConnForCut(
                        ctrl,
                        graph,
                        target as idx_t,
                        ind.len() as idx_t,
                        ind.as_mut_ptr(),
                    ),
                    METIS_OBJTYPE_VOL => MoveGroupMinConnForVol(
                        ctrl,
                        graph,
                        target as idx_t,
                        ind.len() as idx_t,
                        ind.as_mut_ptr(),
                        vmarker.as_nullable_ptr_mut(),
                        pmarker.as_nullable_ptr_mut(),
                        modind.as_nullable_ptr_mut(),
                    ),
                    _ => panic!("Unknown objtype of {}\n", ctrl.objtype),
                }

                /* Update the csr representation of the partitioning vector */
                iarray2csr(nvtxs, nparts, where_, &mut pptr, &mut pind);
            }
        }

        if move_ == 0 {
            break;
        }
    }

    // ipqFree(&queue);

    // WCOREPOP;
}

/*************************************************************************/
/* This function moves a collection of vertices and updates their rinfo */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MoveGroupMinConnForCut(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    to: idx_t,
    nind: idx_t,
    ind: *mut idx_t,
) {
    let ctrl = ctrl.as_mut().unwrap();
    let graph = graph.as_mut().unwrap();
    // idx_t i, ii, j, jj, k, l, nvtxs, nbnd, from, me;
    // idx_t *xadj, *adjncy, *adjwgt, *where_, *bndptr, *bndind;
    // ckrinfo_t *myrinfo;
    // cnbr_t *mynbrs;

    // xadj   = graph.xadj;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    //
    // where_  = graph.where_;
    // bndptr = graph.bndptr;
    // bndind = graph.bndind;

    get_graph_slices!(graph => xadj adjncy adjwgt);
    get_graph_slices_mut!(graph => bndptr bndind where_);

    let mut nbnd: usize = graph.nbnd as usize;
    mkslice!(ind, nind);
    for &i in ind.iter().rev() {
        let i = i as usize;
        let from = where_[i] as usize;

        kwayfm::ensure_crinfo_init(ctrl, graph.ckrinfo, i, xadj);
        let (myrinfo, mynbrs) = kwayfm::crinfos_mut(graph.ckrinfo, ctrl.cnbrpool, i);

        /* find the location of 'to' in myrinfo or create it if it is not there */
        let k = if let Some(k) = mynbrs.iter().position(|nbr| nbr.pid == to) {
            k
        } else {
            let k = mynbrs.len();
            myrinfo.nnbrs += 1;
            // we may temporarily need an extra nbr
            debug_assert!((k as idx_t) <= xadj[i + 1] - xadj[i]);

            let (_myrinfo, mynbrs) = kwayfm::crinfos_mut(graph.ckrinfo, ctrl.cnbrpool, i);
            mynbrs[k].pid = to;
            mynbrs[k].ed = 0;
            k
        };
        // re-construct mynbrs since we may have adjusted nnbrs
        let (myrinfo, mynbrs) = kwayfm::crinfos_mut(graph.ckrinfo, ctrl.cnbrpool, i);

        /* Update pwgts */
        {
            let ncon = graph.ncon as usize;
            get_graph_slices!(graph => vwgt);
            get_graph_slices_mut!(ctrl, graph => pwgts);
            blas::iaxpy(
                ncon,
                1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(to as usize * ncon, ncon)],
                1,
            );
            blas::iaxpy(
                ncon,
                -1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(from * ncon, ncon)],
                1,
            );
        }

        /* Update mincut */
        graph.mincut -= mynbrs[k].ed - myrinfo.id;

        /* Update subdomain connectivity graph to reflect the move of 'i' */
        UpdateEdgeSubDomainGraph(
            ctrl,
            from as idx_t,
            to,
            myrinfo.id - mynbrs[k].ed,
            std::ptr::null_mut(),
        );

        /* Update ID/ED and BND related information for the moved vertex */
        kwayfm::UpdateMovedVertexInfoAndBND(
            i,
            from,
            k,
            to as usize,
            myrinfo,
            mynbrs,
            where_,
            &mut nbnd,
            bndptr,
            bndind,
            BNDTYPE_REFINE,
        );

        /* Update the degrees of adjacent vertices */
        for j in xadj[i]..(xadj[i + 1]) {
            let ii = adjncy[j as usize] as usize;
            let me = where_[ii];
            let orinfo = &mut *graph.ckrinfo.add(ii);

            kwayfm::UpdateAdjacentVertexInfoAndBND(
                ctrl,
                ii,
                xadj[ii + 1] - xadj[ii],
                me as usize,
                from,
                to as usize,
                orinfo,
                adjwgt[j as usize],
                &mut nbnd,
                bndptr,
                bndind,
                BNDTYPE_REFINE,
            );

            /* Update subdomain graph to reflect the move of 'i' for domains other
            than 'from' and 'to' */
            if me as usize != from && me != to {
                UpdateEdgeSubDomainGraph(
                    ctrl,
                    from as idx_t,
                    me,
                    -adjwgt[j as usize],
                    std::ptr::null_mut(),
                );
                UpdateEdgeSubDomainGraph(ctrl, to, me, adjwgt[j as usize], std::ptr::null_mut());
            }
        }

        // debug_assert!(debug::CheckRInfoExtended(ctrl, graph, i as idx_t) != 0);
        debug_assert!(myrinfo.nnbrs <= xadj[i + 1] - xadj[i]);
    }

    debug_assert_eq!(debug::ComputeCut(graph, where_.as_ptr()), graph.mincut);

    graph.nbnd = nbnd as idx_t;
}

/*************************************************************************/
/* This function moves a collection of vertices and updates their rinfo */
/*************************************************************************/
#[metis_func]
pub extern "C" fn MoveGroupMinConnForVol(
    ctrl: *mut ctrl_t,
    graph: *mut graph_t,
    to: idx_t,
    nind: idx_t,
    ind: *mut idx_t,
    vmarker: *mut idx_t,
    pmarker: *mut idx_t,
    modind: *mut idx_t,
) {
    let graph = graph.as_mut().unwrap();
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, j, jj, k, l, nvtxs, from, me, other, xgain, ewgt;
    // idx_t *xadj, *vsize, *adjncy, *where_;
    // vkrinfo_t *myrinfo, *orinfo;
    // vnbr_t *mynbrs, *onbrs;

    // xadj   = graph.xadj;
    // vsize  = graph.vsize;
    // adjncy = graph.adjncy;
    // where_  = graph.where_;
    get_graph_slices!(graph => xadj vsize adjncy);
    get_graph_slices_mut!(graph => where_);
    mkslice!(ind, nind);
    for nind in (0..nind).rev() {
        let i = ind[nind as usize] as usize;
        let from = where_[i as usize] as usize;

        let myrinfo = &mut *graph.vkrinfo.add(i);
        if myrinfo.inbr == -1 {
            myrinfo.inbr = wspace::vnbrpoolGetNext(ctrl, xadj[i + 1] - xadj[i]);
            myrinfo.nnbrs = 0;
        }
        let mynbrs = std::slice::from_raw_parts(
            ctrl.vnbrpool.add(myrinfo.inbr as usize),
            myrinfo.nnbrs as usize,
        );

        let mut xgain = if myrinfo.nid == 0 && myrinfo.ned > 0 {
            vsize[i as usize]
        } else {
            0
        };

        //print!("Moving {:} from {:} to {:} [(vsize: {:}) as usize] [(xgain: {:}) as usize]\n",
        //    i, from, to, vsize[i as usize], xgain);

        /* find the location of 'to' in myrinfo or create it if it is not there */
        // for k in (0)..((*myrinfo).nnbrs as usize) {
        //   if (mynbrs[k as usize].pid == to) {
        //     break;
        //   }
        // }

        let ewgt;

        // control flow inverted from original
        if let Some(k) = mynbrs.iter().position(|nbr| nbr.pid == to) {
            // found the location of 'to' in myrinfo
            graph.minvol -= xgain + mynbrs[k].gv;
            graph.mincut -= mynbrs[k].ned - myrinfo.nid;
            ewgt = myrinfo.nid - mynbrs[k].ned;
        } else {
            // no 'to' in myrinfo, so we'll create it
            //print!("Missing neighbor\n");

            // if myrinfo.nid > 0 {
            xgain -= vsize[i as usize];
            // }

            /* determine the volume gain resulting from that move */
            for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
                let ii = adjncy[j as usize] as usize;
                let other = where_[ii as usize] as usize;
                let (_orinfo, onbrs) = kwayfm::vrinfos_mut(graph.vkrinfo, ctrl.vnbrpool, ii);
                debug_assert!(other != to as usize);

                //print!("  %8d %8d %3d\n", (int)ii, (int)vsize[ii as usize], (int)other);

                if from == other {
                    /* Same subdomain vertex: Decrease the gain if 'to' is a new neighbor. */
                    // for l in (0)..(orinfo.nnbrs) {
                    //   if (onbrs[l as usize].pid == to) {
                    //     break;
                    //   }
                    // }
                    // if (l == orinfo.nnbrs)  {
                    //   xgain -= vsize[ii as usize];
                    // }
                    if onbrs.iter().all(|o| o.pid != to) {
                        xgain -= vsize[ii as usize];
                    }
                } else {
                    /* Remote vertex: increase if 'to' is a new subdomain */
                    // for l in (0)..(orinfo.nnbrs) {
                    //   if (onbrs[l as usize].pid == to) {
                    //     break;
                    //   }
                    // }
                    // if (l == orinfo.nnbrs)  {
                    //   xgain -= vsize[ii as usize];
                    // }
                    // TODO: should this really be subtraction?? That contradicts the comment
                    if onbrs.iter().all(|o| o.pid != to) {
                        xgain -= vsize[ii as usize];
                    }

                    /* Remote vertex: decrease if i is the only connection to 'from' */
                    // for l in (0)..(orinfo.nnbrs) {
                    //   if (onbrs[l as usize].pid == from && onbrs[l as usize].ned == 1) {
                    //     xgain += vsize[ii as usize];
                    //     break;
                    //   }
                    // }
                    if onbrs.iter().any(|o| o.pid == from as idx_t && o.ned == 1) {
                        xgain += vsize[ii as usize];
                    }
                }
            }
            graph.minvol -= xgain;
            graph.mincut -= -myrinfo.nid;
            ewgt = myrinfo.nid;
        }

        /* Update where_ and pwgts */
        where_[i as usize] = to;
        {
            get_graph_slices!(graph => vwgt);
            get_graph_slices_mut!(ctrl, graph => pwgts);
            let ncon = graph.ncon as usize;
            // blas::iaxpy(graph.ncon,  1, vwgt+i*graph.ncon, 1, pwgts+to*graph.ncon,   1);
            // blas::iaxpy(graph.ncon, -1, vwgt+i*graph.ncon, 1, pwgts+from*graph.ncon, 1);
            blas::iaxpy(
                ncon,
                1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(to as usize * ncon, ncon)],
                1,
            );
            blas::iaxpy(
                ncon,
                -1,
                &vwgt[cntrng!(i * ncon, ncon)],
                1,
                &mut pwgts[cntrng!(from * ncon, ncon)],
                1,
            );
        }

        /* Update subdomain connectivity graph to reflect the move of 'i' */
        UpdateEdgeSubDomainGraph(ctrl, from as idx_t, to, ewgt, std::ptr::null_mut());

        /* Update the subdomain connectivity of the adjacent vertices */
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            let me = where_[adjncy[j as usize] as usize];
            if me != from as idx_t && me != to {
                UpdateEdgeSubDomainGraph(ctrl, from as idx_t, me, -1, std::ptr::null_mut());
                UpdateEdgeSubDomainGraph(ctrl, to, me, 1, std::ptr::null_mut());
            }
        }

        /* Update the id/ed/gains/bnd of potentially affected nodes */
        mkslice_mut!(vmarker, graph.nvtxs);
        mkslice_mut!(pmarker, ctrl.nparts);
        mkslice_mut!(modind, graph.nvtxs);
        kwayfm::KWayVolUpdate(
            ctrl,
            graph,
            i,
            from,
            to as usize,
            None,
            None,
            None,
            None,
            None,
            BNDTYPE_REFINE,
            vmarker,
            pmarker,
            modind,
        );

        /*CheckKWayVolPartitionParams(ctrl, graph);*/
    }
    debug_assert_eq!(
        debug::ComputeCutUnweighted(graph, where_.as_ptr()),
        graph.mincut
    );
    debug_assert_eq!(debug::ComputeVolume(graph, where_.as_ptr()), graph.minvol);
}

/*************************************************************************/
/* This function computes the subdomain graph. For deubugging purposes. */
/*************************************************************************/
#[metis_func]
pub extern "C" fn PrintSubDomainGraph(graph: *mut graph_t, nparts: idx_t, where_: *mut idx_t) {
    let nparts = nparts as usize;
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, k, me, nvtxs, total, max;
    // idx_t *xadj, *adjncy, *adjwgt, *pmat;

    let nvtxs = graph.nvtxs as usize;
    // xadj   = graph.xadj;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    get_graph_slices!(graph => xadj adjncy adjwgt);
    mkslice!(where_, nvtxs);

    // pmat = ismalloc(nparts*nparts, 0, "ComputeSubDomainGraph: pmat");
    let mut pmat = vec![0; nparts * nparts];

    for i in (0)..(nvtxs) {
        let me = where_[i as usize] as usize;
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            let k = adjncy[j as usize] as usize;
            if where_[k as usize] as usize != me {
                pmat[(me * nparts + where_[k] as usize) as usize] += adjwgt[j as usize];
            }
        }
    }

    /* print!("Subdomain Info\n"); */
    let mut total = 0;
    let mut max = 0;
    for i in (0)..(nparts) {
        let mut k = 0;
        for j in (0)..(nparts) {
            if pmat[(i * nparts + j) as usize] > 0 {
                k += 1;
            }
        }
        total += k;

        if k > max {
            max = k;
        }
        /*
        print!("{:2} . {:2}  ", i, k);
        for j in (0)..(nparts) {
        if (pmat[(i*nparts+j) as usize] > 0) {
        print!("[({:2} {:4}) as usize] ", j, pmat[(i*nparts+j) as usize]);
        }
        }
        print!("\n");
        */
    }
    println!("Total adjacent subdomains: {:}, Max: {:}", total, max);

    // gk_free((void **)&pmat, LTERM);
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use super::*;
    use crate::tests::{TestGraph, ab_test_partition_test_graphs};

    #[test]
    fn ab_ComputeSubDomainGraph_cut() {
        ab_test_partition_test_graphs(
            "ComputeSubDomainGraph:rs",
            Optype::Kmetis,
            20,
            1,
            |mut g| {
                g.set_minconn(true);
                g.random_adjwgt();
                g.random_tpwgts();
                g
            },
        );
    }

    #[test]
    fn ab_ComputeSubDomainGraph_vol() {
        ab_test_partition_test_graphs(
            "ComputeSubDomainGraph:rs",
            Optype::Kmetis,
            20,
            1,
            |mut g| {
                g.set_minconn(true);
                g.random_vsize();
                g.random_tpwgts();
                g
            },
        );
    }

    #[test]
    fn ab_UpdateEdgeSubDomainGraph() {
        ab_test_partition_test_graphs(
            "UpdateEdgeSubDomainGraph:rs",
            Optype::Kmetis,
            20,
            1,
            |mut g| {
                g.set_minconn(true);
                g.random_adjwgt();
                g.random_tpwgts();
                g
            },
        );
    }

    #[test]
    fn ab_EliminateSubDomainEdges_cut() {
        ab_test_partition_test_graphs(
            "EliminateSubDomainEdges:rs",
            Optype::Kmetis,
            20,
            1,
            |mut g| {
                g.set_minconn(true);
                g.set_objective(Objtype::Cut);
                g.random_adjwgt();
                g.random_tpwgts();
                g
            },
        );
    }

    #[test]
    fn ab_EliminateSubDomainEdges_vol_1() {
        ab_test_partition_test_graphs(
            "EliminateSubDomainEdges:rs",
            Optype::Kmetis,
            80,
            1,
            |mut g| {
                g.set_minconn(true);
                g.set_objective(Objtype::Vol);
                g.random_vsize();
                g.random_tpwgts();
                g
            },
        );
    }

    #[test]
    fn ab_EliminateSubDomainEdges_vol_2() {
        // had some different failures at nparts=20
        ab_test_partition_test_graphs(
            "EliminateSubDomainEdges:rs",
            Optype::Kmetis,
            20,
            1,
            |mut g| {
                g.set_minconn(true);
                g.set_objective(Objtype::Vol);
                g.random_vsize();
                g.random_tpwgts();
                // g.enable_dbg(DbgLvl::ConnInfo);
                g
            },
        );
    }

    #[test]
    fn ab_MoveGroupMinConnForCut() {
        ab_test_partition_test_graphs(
            "MoveGroupMinConnForCut:rs",
            Optype::Kmetis,
            80, // need a decently large npart
            1,
            |mut g| {
                g.set_minconn(true);
                g.set_objective(Objtype::Cut);
                g.random_adjwgt();
                g.random_tpwgts();
                g
            },
        );
    }

    #[test]
    fn ab_MoveGroupMinConnForVol() {
        ab_test_partition_test_graphs(
            "MoveGroupMinConnForVol:rs",
            Optype::Kmetis,
            80,
            1,
            |mut g| {
                g.set_minconn(true);
                g.set_objective(Objtype::Vol);
                g.random_vsize();
                g.random_tpwgts();
                g
            },
        );
    }
}
