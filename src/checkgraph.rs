/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * checkgraph.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 */

use crate::*;
use std::ffi::c_int;

/// This function checks if a graph is valid. A valid graph must satisfy
/// the following constraints:
/// - It should contain no self-edges.
/// - It should be undirected; i.e., (u,v) and (v,u) should be present.
/// - The adjacency list should not contain multiple edges to the same
///   other vertex.
/// 
/// \param graph is the graph to be checked, whose numbering starts from 0.
/// \param numflag is 0 if error reporting will be done using 0 as the
///        numbering, or 1 if the reporting should be done using 1.
/// \param verbose is 1 the identified errors will be displayed, or 0, if
///        it should run silently.
#[metis_func]
pub extern "C" fn CheckGraph(graph: *mut graph_t, numflag: c_int, verbose: c_int) -> c_int {
    let graph = graph.as_ref().unwrap();
    let verbose = verbose != 0;
    // idx_t i, j, k, l;
    // idx_t nvtxs, err=0;
    // idx_t minedge, maxedge, minewgt, maxewgt;
    // idx_t *xadj, *adjncy, *adjwgt, *htable;
    let mut minedge = 0;
    let mut maxedge = 0;
    let mut minewgt = 0;
    let mut maxewgt = 0;
    let mut err = 0;

    let numflag: usize = if numflag == 0 { 0 } else { 1 }; /* make sure that numflag is 0 or 1 */

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy);
    get_graph_slices_optional!(graph => adjwgt);

    // htable = ismalloc(nvtxs, 0, "htable");
    let mut htable = vec![0; nvtxs];

    if graph.nedges > 0 {
        minedge = adjncy[0 as usize];
        maxedge = adjncy[0 as usize];
        if let Some(adjwgt) = adjwgt {
            minewgt = adjwgt[0 as usize];
            maxewgt = adjwgt[0 as usize];
        }
    }

    for i in (0)..(nvtxs) {
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            let j = j as usize;
            let k = adjncy[j as usize] as usize;

            // minedge = (k < minedge) ? k : minedge;
            // maxedge = (k > maxedge) ? k : maxedge;
            minedge = minedge.min(k as idx_t);
            maxedge = maxedge.max(k as idx_t);
            if let Some(adjwgt) = adjwgt {
                // minewgt = (adjwgt[j as usize] < minewgt) ? adjwgt[j as usize] : minewgt;
                // maxewgt = (adjwgt[j as usize] > maxewgt) ? adjwgt[j as usize] : maxewgt;
                minewgt = if adjwgt[j as usize] < minewgt {
                    adjwgt[j as usize]
                } else {
                    minewgt
                };
                maxewgt = if adjwgt[j as usize] > maxewgt {
                    adjwgt[j as usize]
                } else {
                    maxewgt
                };
            }

            if i == k {
                if verbose {
                    println!(
                        "Vertex {:} contains a self-loop \
                  (i.e., diagonal entry in the matrix)!",
                        i + numflag
                    );
                }
                err += 1;
            } else {
                let mut found = false;
                for l in (xadj[k as usize])..(xadj[(k + 1) as usize]) {
                    if adjncy[l as usize] == i as idx_t {
                        if let Some(adjwgt) = adjwgt {
                            if adjwgt[l as usize] != adjwgt[j as usize] {
                                if verbose {
                                    println!(
                                        "Edges (u:{:} v:{:} wgt:{:}) and (v:{:} u:{:} wgt:{:}) do not have the same weight!",
                                        i + numflag,
                                        k + numflag,
                                        adjwgt[j as usize],
                                        k + numflag,
                                        i + numflag,
                                        adjwgt[l as usize]
                                    );
                                }
                                err += 1;
                            }
                        }
                        found = true;
                        break;
                    }
                }
                // if (l == xadj[(k+1) as usize])
                if !found {
                    if verbose {
                        println!("Missing edge: ({:} {:})!", k + numflag, i + numflag);
                    }
                    err += 1;
                }
            }

            if htable[k as usize] == 0 {
                htable[k as usize] += 1;
            } else {
                if verbose {
                    println!(
                        "Edge {:} from vertex {:} is repeated {:} times",
                        k + numflag,
                        i + numflag,
                        htable[k as usize]
                    );
                }
                htable[k] += 1;
                err += 1;
            }
        }

        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            htable[adjncy[j as usize] as usize] = 0;
        }
    }

    if err > 0 && verbose {
        println!(
            "A total of {:} errors exist in the input file. \
            Correct them, and run again!",
            err
        );
    }

    // gk_free((void **)&htable, LTERM);

    return if err == 0 { 1 } else { 0 };
}

/// This function performs a quick check of the weights of the graph - it's broken in the original
/// and unimplemented here
#[metis_func]
#[allow(unused_variables)]
pub extern "C" fn CheckInputGraphWeights(
    nvtxs: idx_t,
    ncon: idx_t,
    xadj: *mut idx_t,
    adjncy: *mut idx_t,
    vwgt: *mut idx_t,
    vsize: *mut idx_t,
    adjwgt: *mut idx_t,
) -> c_int {
    // The original implementation is buggy (out of bounds memory accesses) and unused, so it
    // probably isn't worth bothering to port it
    unimplemented!();

    //   // idx_t i;
    //
    //   if (ncon <= 0) {
    //     print!("Input Error: ncon must be >= 1.\n");
    //     return 0;
    //   }
    //
    //   if (vwgt) {
    // let x = ncon * nvtxs;
    //     for (i=x; i>=0; i--) {
    //       if (vwgt[i as usize] < 0) {
    //         print!("Input Error: negative vertex weight(s).\n");
    //         return 0;
    //       }
    //     }
    //   }
    //   if (vsize) {
    //     for (i=nvtxs; i>=0; i--) {
    //       if (vsize[i as usize] < 0) {
    //         print!("Input Error: negative vertex sizes(s).\n");
    //         return 0;
    //       }
    //     }
    //   }
    //   if (adjwgt) {
    //     for (i=xadj[nvtxs as usize]-1; i>=0; i--) {
    //       if (adjwgt[i as usize] < 0) {
    //         print!("Input Error: non-positive edge weight(s).\n");
    //         return 0;
    //       }
    //     }
    //   }
    //
    //   return 1;
}

// The usage of this means that we don't actually need to compare the `w` field for equality or
// ordering, which may allow for a trivial speedup in the future
#[allow(non_camel_case_types)]
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct uvw_t {
    u: idx_t,
    v: idx_t,
    w: idx_t,
}

/// This function creates a graph whose topology is consistent with
/// Metis' requirements that:
/// - There are no self-edges.
/// - It is undirected; i.e., (u,v) and (v,u) should be present and of the
///   same weight.
/// - The adjacency list should not contain multiple edges to the same
///   other vertex.
/// 
/// Any of the above errors are fixed by performing the following operations:
/// - Self-edges are removed.
/// - The undirected graph is formed by the union of edges.
/// - One of the duplicate edges is selected.
/// 
/// The routine does not change the provided vertex weights.
#[metis_func]
pub extern "C" fn FixGraph(graph: *mut graph_t) -> *mut graph_t {
    // FIXME: I'm taking it on faith that this function is correct since it's only used in the
    // graphchk program. It might be fine since this function is so simple

    let graph = graph.as_mut().unwrap();
    // idx_t i, j, k, l, nvtxs, nedges;
    // idx_t *xadj, *adjncy, *adjwgt;
    // idx_t *nxadj, *nadjncy, *nadjwgt;
    // graph_t *ngraph;
    // uvw_t *edges;

    let nvtxs = graph.nvtxs as usize;
    assert!(!graph.adjwgt.is_null());
    get_graph_slices!(graph => xadj adjncy adjwgt vwgt);
    get_graph_slices_optional!(graph => vsize);
    // xadj   = graph.xadj;
    // adjncy = graph.adjncy;
    // adjwgt = graph.adjwgt;
    // ASSERT(adjwgt != NULL);

    let ngraphp = graph::CreateGraph();
    let ngraph = ngraphp.as_mut().unwrap();

    ngraph.nvtxs = graph.nvtxs;

    /* deal with vertex weights/sizes */
    ngraph.ncon = graph.ncon;
    {
        let mut nvwgt = gk::imalloc(nvtxs * graph.ncon as usize, c"FixGraph: vwgt");
        nvwgt.as_mut().copy_from_slice(vwgt);
        ngraph.vwgt = nvwgt.as_buf_ptr();
        // ngraph.vwgt  = icopy(nvtxs*graph.ncon, graph.vwgt,
        //   imalloc(nvtxs*graph.ncon, "FixGraph: vwgt"));
    }

    {
        let mut nvsize = gk::ismalloc(nvtxs, 1, c"FixGraph: vsize");
        if let Some(vsize) = vsize {
            nvsize.as_mut().copy_from_slice(vsize);
        }
        ngraph.vsize = nvsize.as_buf_ptr()
    }

    /* fix graph by sorting the "superset" of edges */
    // edges = gk_malloc(size_of::<uvw_t>()*2*(xadj[nvtxs] as usize), c"FixGraph: edges").cast::<uvw_t>();
    let mut edges = vec![uvw_t::default(); 2 * (xadj[nvtxs] as usize)];

    let mut nedges = 0;
    for i in (0)..(nvtxs as idx_t) {
        for j in (xadj[i as usize])..(xadj[(i + 1) as usize]) {
            /* keep only the upper-trianglular part of the adjacency matrix */
            if i < adjncy[j as usize] {
                edges[nedges as usize].u = i;
                edges[nedges as usize].v = adjncy[j as usize];
                edges[nedges as usize].w = adjwgt[j as usize];
                nedges += 1;
            } else if i > adjncy[j as usize] {
                edges[nedges as usize].u = adjncy[j as usize];
                edges[nedges as usize].v = i;
                edges[nedges as usize].w = adjwgt[j as usize];
                nedges += 1;
            }
        }
    }

    // could maybe make a trivial speedup here if we make it ignore the `w` field for the key
    edges[..nedges].sort_unstable();
    // uvwsorti(nedges, edges);

    // keep the unique subset
    {
        let mut k = 0;
        for i in (1)..(nedges) {
            if edges[k as usize].v != edges[i as usize].v
                || edges[k as usize].u != edges[i as usize].u
            {
                k += 1;
                edges[k as usize] = edges[i as usize];
            }
        }
        nedges = k + 1;
    }

    /* allocate memory for the fixed graph */
    ngraph.xadj = gk::ismalloc(nvtxs + 1, 0, c"FixGraph: nxadj").as_buf_ptr();
    ngraph.adjncy = gk::imalloc(2 * nedges, c"FixGraph: nadjncy").as_buf_ptr();
    ngraph.adjwgt = gk::imalloc(2 * nedges, c"FixGraph: nadjwgt").as_buf_ptr();
    let nxadj = get_graph_slice_mut!(ngraph => xadj);
    let nadjncy = get_graph_slice_mut!(ngraph => adjncy);
    let nadjwgt = get_graph_slice_mut!(ngraph => adjwgt);

    /* create the adjacency list of the fixed graph from the upper-triangular
    part of the adjacency matrix */
    for k in (0)..(nedges) {
        nxadj[edges[k as usize].u as usize] += 1;
        nxadj[edges[k as usize].v as usize] += 1;
    }
    util::make_csr(nvtxs, nxadj);

    for k in (0)..(nedges) {
        nadjncy[nxadj[edges[k as usize].u as usize] as usize] = edges[k as usize].v;
        nadjncy[nxadj[edges[k as usize].v as usize] as usize] = edges[k as usize].u;
        nadjwgt[nxadj[edges[k as usize].u as usize] as usize] = edges[k as usize].w;
        nadjwgt[nxadj[edges[k as usize].v as usize] as usize] = edges[k as usize].w;
        nxadj[edges[k as usize].u as usize] += 1;
        nxadj[edges[k as usize].v as usize] += 1;
    }
    util::shift_csr(nvtxs, nxadj);

    // gk_free((void **)&edges, LTERM);

    return ngraph;
}
