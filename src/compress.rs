/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * compress.c
 *
 * This file contains code for compressing nodes with identical adjacency
 * structure and for prunning dense columns
 *
 * Started 9/17/97
 * George
 */

use std::ptr;

use crate::*;

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct KeyVal {
    key: idx_t,
    val: idx_t,
}

/// This function compresses a graph by merging identical vertices
///    The compression should lead to at least 10% reduction.
///
///    The compressed graph that is generated has its adjwgts set to 1.
///
///    returns 1 if compression was performed, otherwise it returns 0.
#[metis_func]
pub extern "C" fn CompressGraph(
    ctrl: *mut ctrl_t,
    nvtxs: idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    cptr: *mut idx_t,
    cind: *mut idx_t,
) -> *mut graph_t {
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, ii, iii, j, jj, k, l, cnvtxs, cnedges;
    // idx_t *cxadj, *cadjncy, *cvwgt, *mark, *map;
    // ikv_t *keys;
    // graph_t *graph=NULL;

    // mark = ismalloc(nvtxs, -1, "CompressGraph: mark");
    // map = ismalloc(nvtxs, -1, "CompressGraph: map");
    // keys = ikvmalloc(nvtxs, "CompressGraph: keys");
    let nvtxs = nvtxs as usize;
    let mut mark: Vec<idx_t> = vec![-1; nvtxs];
    let mut map: Vec<idx_t> = vec![-1; nvtxs];
    let mut keys = vec![KeyVal::default(); nvtxs];

    mkslice!(xadj, nvtxs + 1);
    mkslice!(adjncy, xadj[nvtxs]);
    mkslice_mut!(cptr, nvtxs);
    mkslice_mut!(cind, nvtxs);

    let mut graph: *mut graph_t = ptr::null_mut();

    /* Compute a key for each adjacency list */
    for i in (0)..(nvtxs) {
        let mut k = 0;
        for j in (xadj[i])..(xadj[i + 1]) {
            k += adjncy[j as usize];
        }
        keys[i].key = k + i as idx_t; /* Add the diagonal entry as well */
        keys[i].val = i as idx_t;
    }

    // ikvsorti(nvtxs, keys);
    keys.sort_unstable_by_key(|kv| kv.key);

    cptr[0] = 0;
    let mut l = 0;
    let mut cnvtxs: usize = 0;
    for i in (0)..(nvtxs) {
        let ii = keys[i].val as usize;
        if map[ii] == -1 {
            mark[ii] = i as idx_t; /* Add the diagonal entry */
            for j in (xadj[ii])..(xadj[ii + 1]) {
                mark[adjncy[j as usize] as usize] = i as idx_t;
            }

            map[ii] = cnvtxs as idx_t;
            cind[l] = ii as idx_t;
            l += 1;

            for j in (i + 1)..(nvtxs) {
                let iii = keys[j].val as usize;

                if keys[i].key != keys[j].key
                    || xadj[ii + 1] - xadj[ii] != xadj[iii + 1] - xadj[iii]
                {
                    break; /* Break if keys or degrees are different */
                }

                if map[iii] == -1 {
                    /* Do a comparison if iii has not been mapped */
                    let mut jj = xadj[iii];
                    for jjj in (xadj[iii])..(xadj[iii + 1]) {
                        jj = jjj;
                        if mark[adjncy[jj as usize] as usize] != i as idx_t {
                            break;
                        }
                    }

                    if jj == xadj[iii + 1] {
                        /* Identical adjacency structure */
                        map[iii] = cnvtxs as idx_t;
                        cind[l] = iii as idx_t;
                        l += 1;
                    }
                }
            }

            cnvtxs += 1;
            cptr[cnvtxs as usize] = l as idx_t;
        }
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_INFO,
        println!(
            "  Compression: reduction in # of vertices: {:}.\n",
            nvtxs - cnvtxs,
        ),
    );

    if (cnvtxs as real_t) < COMPRESSION_FRACTION * nvtxs as real_t {
        /* Sufficient compression is possible, so go ahead and create the
        compressed graph */

        graph = graph::CreateGraph();
        let graph = graph.as_mut().unwrap();

        let mut cnedges = 0;
        for i in (0)..(cnvtxs) {
            let ii = cind[cptr[i as usize] as usize] as usize;
            cnedges += xadj[ii + 1] - xadj[ii];
        }

        /* Allocate memory for the compressed graph */
        graph.xadj = imalloc(cnvtxs as usize + 1, "CompressGraph: xadj\0".as_ptr()) as _;
        graph.vwgt = imalloc(cnvtxs as usize, "CompressGraph: vwgt\0".as_ptr()) as _;
        graph.vwgt.write_bytes(0, cnvtxs as usize); // prolly unnecessary, but matches og
        graph.adjncy = imalloc(cnedges as usize, "CompressGraph: adjncy\0".as_ptr()) as _;
        mkslice_mut!(cxadj: graph->xadj, cnvtxs + 1);
        mkslice_mut!(cvwgt: graph->vwgt, cnvtxs);
        mkslice_mut!(cadjncy: graph->adjncy, cnedges);
        {
            graph.adjwgt = imalloc(cnedges as usize, "CompressGraph: adjwgt\0".as_ptr()) as _;
            graph.adjwgt.write_bytes(1, cnedges as usize);
        }

        /* Now go and compress the graph */
        // iset(nvtxs, -1, mark);
        mark.fill(-1);

        let mut l = 0;
        cxadj[0] = 0;
        for i in (0)..(cnvtxs) {
            mark[i] = i as idx_t; /* Remove any dioganal entries in the compressed graph */
            for j in (cptr[i])..(cptr[i + 1]) {
                let ii = cind[j as usize] as usize;

                /* accumulate the vertex weights of the consistuent vertices */
                cvwgt[i] += if vwgt.is_null() { 1 } else { *vwgt.add(ii) };

                /* generate the combined adjancency list */
                for jj in (xadj[ii])..(xadj[ii + 1]) {
                    let k = map[adjncy[jj as usize] as usize] as usize;
                    if mark[k] != i as idx_t {
                        mark[k] = i as idx_t;
                        cadjncy[l] = k as idx_t;
                        l += 1;
                    }
                }
            }
            cxadj[i + 1] = l as idx_t;
        }

        graph.nvtxs = cnvtxs as idx_t;
        graph.nedges = l as idx_t;
        graph.ncon = 1;

        graph::SetupGraph_tvwgt(graph);
        graph::SetupGraph_label(graph);
    };

    // gk_free((void **)&keys, &map, &mark, LTERM);

    return graph;
}

/// This function prunes all the vertices in a graph with degree greater
///    than `factor * average`.
///
///    returns the number of vertices that were prunned.
#[metis_func]
pub extern "C" fn PruneGraph(
    ctrl: *mut ctrl_t,
    nvtxs: idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    iperm: *mut idx_t,
    factor: real_t,
) -> *mut graph_t {
    let ctrl = ctrl.as_mut().unwrap();
    // idx_t i, j, k, l, nlarge, pnvtxs, pnedges;
    // idx_t *pxadj, *padjncy, *padjwgt, *pvwgt;
    // idx_t *perm;
    // graph_t *graph=NULL;

    let mut graph: *mut graph_t = ptr::null_mut();

    let nvtxs = nvtxs as usize;
    // let perm = imalloc(nvtxs, "PruneGraph: perm\0".as_ptr()) as _;
    let mut perm: Vec<idx_t> = vec![0; nvtxs];

    mkslice!(xadj, nvtxs + 1);
    mkslice!(adjncy, xadj[nvtxs]);
    mkslice_mut!(iperm, nvtxs);

    let factor = factor * (xadj[nvtxs] / nvtxs as idx_t) as real_t;

    let mut pnvtxs: usize = 0;
    let mut pnedges: usize = 0;
    let mut nlarge: usize = 0;
    for i in (0)..(nvtxs) {
        if ((xadj[i + 1] - xadj[i]) as real_t) < factor {
            perm[i] = pnvtxs as idx_t;
            iperm[pnvtxs] = i as idx_t;
            pnvtxs += 1;
            pnedges += xadj[i + 1] as usize - xadj[i] as usize;
        } else {
            nlarge += 1;
            perm[i] = nvtxs as idx_t - nlarge as idx_t;
            iperm[nvtxs - nlarge] = i as idx_t;
        }
    }

    ifset!(
        ctrl.dbglvl,
        METIS_DBG_INFO,
        println!("  Pruned {:} of {:} vertices.\n", nlarge, nvtxs),
    );

    if nlarge > 0 && nlarge < nvtxs {
        /* Prunning is possible, so go ahead and create the prunned graph */
        graph = graph::CreateGraph();
        let graph = graph.as_mut().unwrap();

        /* Allocate memory for the prunned graph*/
        graph.xadj = imalloc(pnvtxs + 1, "PruneGraph: xadj\0".as_ptr()) as _;
        graph.vwgt = imalloc(pnvtxs, "PruneGraph: vwgt\0".as_ptr()) as _;
        graph.adjncy = imalloc(pnedges, "PruneGraph: adjncy\0".as_ptr()) as _;
        // graph.vwgt = ismalloc(pnedges, 1, "PruneGraph: adjwgt") as _;
        graph.adjwgt = imalloc(pnedges, "PruneGraph: adjwgt\0".as_ptr()) as _;
        mkslice_mut!(pxadj: graph->xadj, pnvtxs + 1);
        mkslice_mut!(pvwgt: graph->vwgt, pnvtxs);
        pvwgt.fill(1);
        mkslice_mut!(padjncy: graph->adjncy, pnedges);

        pxadj[0] = 0;
        pnedges = 0;
        let mut l = 0;
        for i in (0)..(nvtxs) {
            if ((xadj[i + 1] - xadj[i]) as real_t) < factor {
                pvwgt[l] = if vwgt.is_null() { 1 } else { *vwgt.add(i) };

                for j in (xadj[i])..(xadj[i + 1]) {
                    let k = perm[adjncy[j as usize] as usize];
                    if k < pnvtxs as idx_t {
                        padjncy[pnedges] = k;
                        pnedges += 1;
                    }
                }
                l += 1;
                pxadj[l as usize] = pnedges as idx_t;
            }
        }

        graph.nvtxs = pnvtxs as idx_t;
        graph.nedges = pnedges as idx_t;
        graph.ncon = 1;

        graph::SetupGraph_tvwgt(graph);
        graph::SetupGraph_label(graph);
    } else if nlarge > 0 && nlarge == nvtxs {
        ifset!(
            ctrl.dbglvl,
            METIS_DBG_INFO,
            println!("  Pruning is ignored as it removes all vertices.\n")
        );
        // set in original, but doesn't appear to do anything
        // nlarge = 0;
    }

    // gk_free((void **)&perm, LTERM);

    return graph;
}
