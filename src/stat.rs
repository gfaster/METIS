/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * stat.c
 *
 * This file computes various statistics
 *
 * Started 7/25/97
 * George
 *
 * $Id: stat.c 17513 2014-08-05 16:20:50Z dominique $
 *
 */

use crate::{params::params_t, *};

/*************************************************************************
* This function computes cuts and balance information
**************************************************************************/
#[metis_func]
pub extern "C" fn ComputePartitionInfoBipartite(
    graph: *const graph_t,
    nparts: idx_t,
    where_: *const idx_t,
) {
    let graph = graph.as_ref().unwrap();
    let nparts = nparts as usize;
    // idx_t i, j, k, nvtxs, ncon, mustfree=0;
    // idx_t *xadj, *adjncy, *vwgt, *vsize, *adjwgt, *kpwgts, *tmpptr;
    // idx_t *padjncy, *padjwgt, *padjcut;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    get_graph_slices!(graph => xadj adjncy vsize);
    get_graph_slices_optional!(graph => vwgt adjwgt);
    mkslice!(where_, nvtxs);

    let vwgt_vec: Vec<idx_t>;
    let vwgt: &[idx_t] = if vwgt.is_none() {
        vwgt_vec = vec![1; nvtxs];
        &vwgt_vec
    } else {
        vwgt.unwrap()
    };
    let adjwgt_vec: Vec<idx_t>;
    let adjwgt: &[idx_t] = if adjwgt.is_none() {
        adjwgt_vec = vec![1; xadj[nvtxs] as usize];
        &adjwgt_vec
    } else {
        adjwgt.unwrap()
    };

    print!(
        "{:}-way Cut: {:5}, Vol: {:5}, ",
        nparts,
        debug::ComputeCut(graph, where_.as_ptr()),
        debug::ComputeVolume(graph, where_.as_ptr())
    );

    /* Compute balance information */
    // let kpwgts = ismalloc(ncon*nparts, 0, "ComputePartitionInfo: kpwgts");
    let mut kpwgts: Vec<idx_t> = vec![0; ncon * nparts];

    for i in (0)..(nvtxs) {
        for j in (0)..(ncon) {
            kpwgts[where_[i as usize] as usize * ncon + j] += vwgt[(i * ncon + j) as usize];
        }
    }

    if ncon == 1 {
        print!(
            "\tBalance: {:5.3} out of {:5.3}\n",
            (nparts as real_t) * (kpwgts[util::iargmax(&kpwgts[..nparts], 1)] as real_t)
                / (kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t),
            (nparts as real_t) * (vwgt[util::iargmax(&vwgt[..nvtxs], 1)] as real_t)
                / (kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t)
        );
    } else {
        print!("\tBalance:");
        for j in (0)..(ncon) {
            print!(
                " ({:5.3} out of {:5.3})",
                (nparts as real_t)
                    * (kpwgts[(ncon * util::iargmax(&kpwgts[j..], ncon) + j) as usize] as real_t)
                    / (kpwgts[cntrng!((j), (nparts))]
                        .iter()
                        .step_by(ncon)
                        .sum::<idx_t>() as real_t),
                (nparts as real_t)
                    * (vwgt[(ncon * util::iargmax(&vwgt[j..], ncon) + j) as usize] as real_t)
                    / (kpwgts[cntrng!((j), (nparts))]
                        .iter()
                        .step_by(ncon)
                        .sum::<idx_t>() as real_t)
            );
        }
        print!("\n");
    }

    /* Compute p-adjncy information */
    // padjncy = ismalloc(nparts*nparts, 0, "ComputePartitionInfo: padjncy");
    // padjwgt = ismalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");
    // padjcut = ismalloc(nparts*nparts, 0, "ComputePartitionInfo: padjwgt");
    let mut padjncy: Vec<idx_t> = vec![0; nparts * nparts];
    let mut padjwgt: Vec<idx_t> = vec![0; nparts * nparts];
    let mut padjcut: Vec<idx_t> = vec![0; nparts * nparts];

    kpwgts[..nparts].fill(0);
    for i in (0)..(nvtxs) {
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            if where_[i as usize] != where_[adjncy[j as usize] as usize] {
                padjncy[where_[i as usize] as usize * nparts
                    + where_[adjncy[j as usize] as usize] as usize] = 1;
                padjcut[where_[i as usize] as usize * nparts
                    + where_[adjncy[j as usize] as usize] as usize] += adjwgt[j as usize];
                if kpwgts[where_[adjncy[j as usize] as usize] as usize] == 0 {
                    padjwgt[where_[i as usize] as usize * nparts
                        + where_[adjncy[j as usize] as usize] as usize] += vsize[i as usize];
                    kpwgts[where_[adjncy[j as usize] as usize] as usize] = 1;
                }
            }
        }
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            kpwgts[where_[adjncy[j as usize] as usize] as usize] = 0;
        }
    }

    for i in (0)..(nparts) {
        kpwgts[i as usize] = padjncy[cntrng!((i * nparts), (nparts))]
            .iter()
            .sum::<idx_t>();
    }
    print!(
        "Min/Max/Avg/Bal # of adjacent     subdomains: {:5} {:5} {:5} {:7.3}\n",
        kpwgts[util::iargmin(&kpwgts[..nparts], 1)],
        kpwgts[util::iargmax(&kpwgts[..nparts], 1)],
        kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() / (nparts as idx_t),
        (nparts as real_t) * (kpwgts[util::iargmax(&kpwgts[..nparts], 1)] as real_t)
            / (kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t)
    );

    for i in (0)..(nparts) {
        kpwgts[i as usize] = padjcut[cntrng!((i * nparts), (nparts))]
            .iter()
            .sum::<idx_t>();
    }
    print!(
        "Min/Max/Avg/Bal # of adjacent subdomain cuts: {:5} {:5} {:5} {:7.3}\n",
        kpwgts[util::iargmin(&kpwgts[..nparts], 1)],
        kpwgts[util::iargmax(&kpwgts[..nparts], 1)],
        kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() / (nparts as idx_t),
        (nparts as real_t) * (kpwgts[util::iargmax(&kpwgts[..nparts], 1)] as real_t)
            / (kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t)
    );

    for i in (0)..(nparts) {
        kpwgts[i as usize] = padjwgt[cntrng!(i * nparts, nparts)].iter().sum::<idx_t>();
    }
    print!(
        "Min/Max/Avg/Bal/Frac # of interface    nodes: {:5} {:5} {:5} {:7.3} {:7.3}\n",
        kpwgts[util::iargmin(&kpwgts[..nparts], 1)],
        kpwgts[util::iargmax(&kpwgts[..nparts], 1)],
        kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() / (nparts as idx_t),
        (nparts as real_t) * (kpwgts[util::iargmax(&kpwgts[..nparts], 1)] as real_t)
            / (kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t),
        kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t / (nvtxs as real_t)
    );
}

/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
#[metis_func]
pub extern "C" fn ComputePartitionBalance(
    graph: *const graph_t,
    nparts: idx_t,
    where_: *const idx_t,
    ubvec: *mut real_t,
) {
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, nvtxs, ncon;
    // idx_t *kpwgts, *vwgt;
    // real_t balance;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    get_graph_slices_optional!(graph => vwgt);
    let nparts = nparts as usize;
    mkslice_mut!(ubvec, ncon);
    mkslice!(where_, nvtxs);

    let mut kpwgts = vec![0; nparts];

    if let Some(vwgt) = vwgt {
        for j in (0)..(ncon) {
            kpwgts[..nparts].fill(0);
            for i in (0)..(nvtxs) {
                kpwgts[where_[i as usize] as usize] += vwgt[(i * ncon + j) as usize];
            }

            ubvec[j as usize] = (nparts as real_t)
                * (kpwgts[util::iargmax(&kpwgts[..nparts], 1)] as real_t)
                / (kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t);
        }
    } else {
        for i in (0)..(nvtxs) {
            kpwgts[where_[i as usize] as usize] += 1;
        }
        ubvec[0 as usize] = (nparts as real_t)
            * (kpwgts[util::iargmax(&kpwgts[..nparts], 1)] as real_t)
            / (nvtxs as real_t);
    }
}

/*************************************************************************
* This function computes the balance of the element partitioning
**************************************************************************/
#[metis_func]
pub extern "C" fn ComputeElementBalance(ne: idx_t, nparts: idx_t, where_: *mut idx_t) -> real_t {
    let ne = ne as usize;
    let nparts = nparts as usize;
    // idx_t i;
    // idx_t *kpwgts;
    // real_t balance;

    let mut kpwgts = vec![0; nparts];
    mkslice!(where_, ne);
    for &i in where_ {
        kpwgts[i as usize] += 1;
    }

    (nparts as real_t) * (kpwgts[util::iargmax(&kpwgts[..nparts], 1)] as real_t)
        / (kpwgts[..nparts].iter().step_by(1).sum::<idx_t>() as real_t)
}

#[expect(non_snake_case)]
pub unsafe fn ComputePartitionInfo(params: &params_t, graph: &graph_t, where_: &[idx_t]) {
    // idx_t i, ii, j, k, nvtxs, ncon, nparts, tvwgt;
    // idx_t *xadj, *adjncy, *vwgt, *adjwgt, *kpwgts;
    // real_t *tpwgts, unbalance;
    // idx_t pid, ndom, maxndom, minndom, tndom, *pptr, *pind, *pdom;
    // idx_t ncmps, nover, *cptr, *cind, *cpwgts;

    let nvtxs = graph.nvtxs as usize;
    let ncon = graph.ncon as usize;
    // TODO: I don't think the C version checked adjwgt is null, and reading in a graph may not
    // allocate adjwgt. what was happening there?
    get_graph_slices!(graph => xadj adjncy vwgt);
    get_graph_slices_optional!(graph => adjwgt);
    let nparts = params.nparts as usize;
    let tpwgts = &params.tpwgts;

    /* Compute objective-related information */
    print!(
        " - Edgecut: {:}, communication volume: {:}.\n\n",
        debug::ComputeCut(graph, where_.as_ptr()),
        debug::ComputeVolume(graph, where_.as_ptr())
    );

    /* Compute constraint-related information */
    // let mut kpwgts = ismalloc(ncon*nparts, 0, "ComputePartitionInfo: kpwgts");
    let mut kpwgts = vec![0; ncon * nparts];

    for i in 0..nvtxs {
        for j in 0..ncon {
            kpwgts[where_[i] as usize * ncon + j] += vwgt[i * ncon + j];
        }
    }

    /* Report on balance */
    print!(" - Balance:\n");
    for j in 0..ncon {
        let tvwgt = kpwgts[cntrng!((j), (nparts))]
            .iter()
            .step_by(ncon)
            .sum::<idx_t>();
        // let mut k = 0;
        let mut unbalance =
            kpwgts[0 * ncon + j] as real_t / (tpwgts[0 * ncon + j] * tvwgt as real_t);
        for i in 1..nparts {
            if unbalance < kpwgts[i * ncon + j] as real_t / (tpwgts[i * ncon + j] * tvwgt as real_t)
            {
                unbalance =
                    kpwgts[i * ncon + j] as real_t / (tpwgts[i * ncon + j] * tvwgt as real_t);
                // k = i;
            }
        }
        print!(
            "     constraint #{:}:  {:5.3} out of {:5.3}\n",
            j,
            unbalance,
            nparts as real_t * vwgt[ncon * util::iargmax(&vwgt[j..], ncon) + j] as real_t
                / (kpwgts[cntrng!((j), (nparts))]
                    .iter()
                    .step_by(ncon)
                    .sum::<idx_t>() as real_t)
        );
    }
    print!("\n");

    if ncon == 1 {
        let tvwgt = kpwgts[..nparts].iter().step_by(1).sum::<idx_t>();
        let mut k = 0;
        let mut unbalance = kpwgts[k] as real_t / (tpwgts[k] * tvwgt as real_t);
        for i in 1..nparts {
            if unbalance < kpwgts[i] as real_t / (tpwgts[i] * tvwgt as real_t) {
                unbalance = kpwgts[i] as real_t / (tpwgts[i] * tvwgt as real_t);
                k = i;
            }
        }

        print!(
            " - Most overweight partition:\n     \
           pid: {:}, actual: {:}, desired: {:}, ratio: {:.2}.\n\n",
            k,
            kpwgts[k],
            (tvwgt as real_t * tpwgts[k]) as idx_t,
            unbalance
        );
    }

    drop(kpwgts);

    /* Compute subdomain adjacency information */
    let mut pptr: Vec<idx_t> = vec![0; nparts + 1];
    let mut pind: Vec<idx_t> = vec![0; nvtxs];
    let mut pdom: Vec<idx_t> = vec![0; nparts];

    minconn::iarray2csr(nvtxs, nparts, where_, &mut pptr, &mut pind);

    let mut maxndom = nparts + 1;
    let mut minndom = 0;
    let mut tndom = 0;
    for pid in 0..nparts {
        pdom[..nparts].fill(0);
        for ii in pptr[pid]..pptr[pid + 1] {
            let i = pind[ii as usize] as usize;
            for j in xadj[i]..xadj[i + 1] {
                pdom[where_[adjncy[j as usize] as usize] as usize] += adjwgt.map_or(1, |a| a[j as usize]);
            }
        }
        pdom[pid] = 0;
        let mut ndom = 0;
        for i in 0..nparts {
            ndom += if pdom[i] > 0 { 1 } else { 0 };
        }
        tndom += ndom;
        if pid == 0 || maxndom < ndom {
            maxndom = ndom;
        }
        if pid == 0 || minndom > ndom {
            minndom = ndom;
        }
    }

    print!(
        " - Subdomain connectivity: max: {:}, min: {:}, avg: {:.2}\n\n",
        maxndom,
        minndom,
        (tndom as f64) / (nparts as f64)
    );

    drop(pptr);
    drop(pind);
    drop(pdom);

    /* Compute subdomain adjacency information */
    // cptr   = imalloc(nvtxs+1, "ComputePartitionInfo: cptr");
    // cind   = imalloc(nvtxs, "ComputePartitionInfo: cind");
    // cpwgts = ismalloc(nparts, 0, "ComputePartitionInfo: cpwgts");
    let mut cptr = vec![0; nvtxs + 1];
    let mut cind = vec![0; nvtxs];
    let mut cpwgts = vec![0; nparts];

    let ncmps = contig::FindPartitionInducedComponents(
        graph,
        where_.as_ptr(),
        cptr.as_mut_ptr(),
        cind.as_mut_ptr(),
    ) as usize;
    if ncmps == nparts {
        print!(" - Each partition is contiguous.\n");
    } else {
        if contig::IsConnected(graph, 0) != 0 {
            let mut nover = 0;
            for i in 0..ncmps {
                cpwgts[where_[cind[cptr[i] as usize] as usize] as usize] += 1;
                if cpwgts[where_[cind[cptr[i] as usize] as usize] as usize] == 2 {
                    nover += 1;
                }
            }
            print!(
                " - There are {:} non-contiguous partitions.\n   \
                    Total components after removing the cut edges: {:},\n   \
                    max components: {:} for pid: {:}.\n",
                nover,
                ncmps,
                cpwgts[..nparts]
                    .iter()
                    .copied()
                    .step_by(1)
                    .max()
                    .unwrap_or_default(),
                util::iargmax(&cpwgts, 1)
            );
        } else {
            print!(
                " - The original graph had {:} connected components and the resulting\n   \
                 partitioning after removing the cut edges has {:} components.",
                contig::FindPartitionInducedComponents(
                    graph,
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    std::ptr::null_mut()
                ),
                ncmps
            );
        }
    }

    // gk_free((void **)&cptr, &cind, &cpwgts, LTERM);
}
