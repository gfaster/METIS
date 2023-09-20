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

use crate::*;

/*************************************************************************/
/* This function compresses a graph by merging identical vertices
    The compression should lead to at least 10% reduction.

    The compressed graph that is generated has its adjwgts set to 1.

    \returns 1 if compression was performed, otherwise it returns 0.

*/
/**************************************************************************/
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
    // idx_t i, ii, iii, j, jj, k, l, cnvtxs, cnedges;
    // idx_t *cxadj, *cadjncy, *cvwgt, *mark, *map;
    // ikv_t *keys;
    // graph_t *graph=NULL;

    mark = ismalloc(nvtxs, -1, "CompressGraph: mark");
    map = ismalloc(nvtxs, -1, "CompressGraph: map");
    keys = ikvmalloc(nvtxs, "CompressGraph: keys");

    /* Compute a key for each adjacency list */
    for i in (0)..(nvtxs) {
        k = 0;
        for j in (xadj[i])..(xadj[i + 1]) {
            k += adjncy[j];
        }
        keys[i].key = k + i; /* Add the diagonal entry as well */
        keys[i].val = i;
    }

    ikvsorti(nvtxs, keys);

    l = cptr[0] = 0;
    cnvtxs = 0;
    for i in (0)..(nvtxs) {
        ii = keys[i].val;
        if (map[ii] == -1) {
            mark[ii] = i; /* Add the diagonal entry */
            for j in (xadj[ii])..(xadj[ii + 1]) {
                mark[adjncy[j]] = i;
            }

            map[ii] = cnvtxs;
            cind[l] = ii;
            l += 1;

            for j in (i + 1)..(nvtxs) {
                iii = keys[j].val;

                if (keys[i].key != keys[j].key
                    || xadj[ii + 1] - xadj[ii] != xadj[iii + 1] - xadj[iii])
                {
                    break; /* Break if keys or degrees are different */
                }

                if (map[iii] == -1) {
                    /* Do a comparison if iii has not been mapped */
                    for jj in (xadj[iii])..(xadj[iii + 1]) {
                        if (mark[adjncy[jj]] != i) {
                            break;
                        }
                    }

                    if (jj == xadj[iii + 1]) {
                        /* Identical adjacency structure */
                        map[iii] = cnvtxs;
                        cind[l] = iii;
                        l += 1;
                    }
                }
            }

            cnvtxs += 1;
            cptr[cnvtxs] = l;
        }
    }

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_INFO,
        printf(
            "  Compression: reduction in # of vertices: {:}.\n",
            nvtxs - cnvtxs,
        ),
    );

    if (cnvtxs < COMPRESSION_FRACTION * nvtxs) {
        /* Sufficient compression is possible, so go ahead and create the
        compressed graph */

        graph = CreateGraph();

        cnedges = 0;
        for i in (0)..(cnvtxs) {
            ii = cind[cptr[i]];
            cnedges += xadj[ii + 1] - xadj[ii];
        }

        /* Allocate memory for the compressed graph */
        cxadj = graph.xadj = imalloc(cnvtxs + 1, "CompressGraph: xadj");
        cvwgt = graph.vwgt = ismalloc(cnvtxs, 0, "CompressGraph: vwgt");
        cadjncy = graph.adjncy = imalloc(cnedges, "CompressGraph: adjncy");
        graph.adjwgt = ismalloc(cnedges, 1, "CompressGraph: adjwgt");

        /* Now go and compress the graph */
        iset(nvtxs, -1, mark);
        l = cxadj[0] = 0;
        for i in (0)..(cnvtxs) {
            mark[i] = i; /* Remove any dioganal entries in the compressed graph */
            for j in (cptr[i])..(cptr[i + 1]) {
                ii = cind[j];

                /* accumulate the vertex weights of the consistuent vertices */
                cvwgt[i] += (if vwgt == NULL { 1 } else { vwgt[ii] });

                /* generate the combined adjancency list */
                for jj in (xadj[ii])..(xadj[ii + 1]) {
                    k = map[adjncy[jj]];
                    if (mark[k] != i) {
                        mark[k] = i;
                        cadjncy[l] = k;
                        l += 1;
                    }
                }
            }
            cxadj[i + 1] = l;
        }

        graph.nvtxs = cnvtxs;
        graph.nedges = l;
        graph.ncon = 1;

        SetupGraph_tvwgt(graph);
        SetupGraph_label(graph);
    }

    // gk_free((void **)&keys, &map, &mark, LTERM);

    return graph;
}

/*************************************************************************/
/* This function prunes all the vertices in a graph with degree greater
    than factor*average.

    \returns the number of vertices that were prunned.
*/
/*************************************************************************/
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
    // idx_t i, j, k, l, nlarge, pnvtxs, pnedges;
    // idx_t *pxadj, *padjncy, *padjwgt, *pvwgt;
    // idx_t *perm;
    // graph_t *graph=NULL;

    perm = imalloc(nvtxs, "PruneGraph: perm");

    factor = factor * xadj[nvtxs] / nvtxs;

    pnvtxs = pnedges = nlarge = 0;
    for i in (0)..(nvtxs) {
        if (xadj[i + 1] - xadj[i] < factor) {
            perm[i] = pnvtxs;
            iperm[pnvtxs] = i;
            pnvtxs += 1;
            pnedges += xadj[i + 1] - xadj[i];
        } else {
            nlarge += 1;
            perm[i] = nvtxs - nlarge;
            iperm[nvtxs - nlarge] = i;
        }
    }

    IFSET(
        ctrl.dbglvl,
        METIS_DBG_INFO,
        printf("  Pruned {:} of {:} vertices.\n", nlarge, nvtxs),
    );

    if (nlarge > 0 && nlarge < nvtxs) {
        /* Prunning is possible, so go ahead and create the prunned graph */
        graph = CreateGraph();

        /* Allocate memory for the prunned graph*/
        pxadj = graph.xadj = imalloc(pnvtxs + 1, "PruneGraph: xadj");
        pvwgt = graph.vwgt = imalloc(pnvtxs, "PruneGraph: vwgt");
        padjncy = graph.adjncy = imalloc(pnedges, "PruneGraph: adjncy");
        graph.adjwgt = ismalloc(pnedges, 1, "PruneGraph: adjwgt");

        pxadj[0] = pnedges = l = 0;
        for i in (0)..(nvtxs) {
            if (xadj[i + 1] - xadj[i] < factor) {
                pvwgt[l] = (if vwgt == NULL { 1 } else { vwgt[i] });

                for j in (xadj[i])..(xadj[i + 1]) {
                    k = perm[adjncy[j]];
                    if (k < pnvtxs) {
                        padjncy[pnedges] = k;
                        pnedges += 1;
                    }
                }
                l += 1;
                pxadj[l] = pnedges;
            }
        }

        graph.nvtxs = pnvtxs;
        graph.nedges = pnedges;
        graph.ncon = 1;

        SetupGraph_tvwgt(graph);
        SetupGraph_label(graph);
    } else if (nlarge > 0 && nlarge == nvtxs) {
        IFSET(
            ctrl.dbglvl,
            METIS_DBG_INFO,
            printf("  Pruning is ignored as it removes all vertices.\n"),
        );
        nlarge = 0;
    }

    // gk_free((void **)&perm, LTERM);

    return graph;
}
