//! partial port of debug.c

use crate::*;

/// This function checks whether or not the boundary information is correct
#[metis_func]
pub extern "C" fn CheckBnd2(graph: *const graph_t) -> idx_t {
    let graph = graph.as_ref().unwrap();
    let nvtxs = graph.nvtxs;
    mkslice!(graph->xadj, nvtxs + 1);
    mkslice!(graph->adjncy, xadj[xadj.len() - 1]);
    let where_ = graph.where_;
    mkslice!(graph->bndptr, nvtxs);
    mkslice!(graph->bndind, nvtxs);
    mkslice!(graph->adjwgt, graph.nedges);
    let mut id;
    let mut ed;

    let mut nbnd = 0;

    for i in 0..(nvtxs as usize) {
        id = 0;
        ed = 0;
        for j in xadj[i]..xadj[i + 1] {
            let j = j as usize;
            if *where_.add(i) != *where_.add(adjncy[j] as usize) {
                ed += adjwgt[j];
            } else {
                id += adjwgt[j];
            }
        }

        if ed - id >= 0 && xadj[i] < xadj[i + 1] {
            nbnd += 1;

            assert_ne!(bndptr[i], -1, "{i} {id} {ed}");
            assert_eq!(bndind[bndptr[i] as usize], i as idx_t);
        }
    }

    assert_eq!(nbnd, graph.nbnd);

    return 1;
}
