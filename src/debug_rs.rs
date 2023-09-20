//! partial port of debug.c

use crate::*;

#[metis_func]
pub extern "C" fn CheckBnd(graph: *const graph_t) -> idx_t {
    let graph = graph.as_ref().unwrap();
    // idx_t i, j, nvtxs, nbnd;
    // idx_t *xadj, *adjncy, *where, *bndptr, *bndind;

    let nvtxs = graph.nvtxs as usize;
    get_graph_slices!(graph => xadj adjncy where_ bndptr bndind);

    let mut nbnd = 0;
    for i in 0..nvtxs {
        if xadj[i + 1] - xadj[i] == 0 {
            nbnd += 1; /* Islands are considered to be boundary vertices */
        }
        for j in xadj[i]..xadj[i + 1] {
            if where_[i] != where_[adjncy[j as usize] as usize] {
                nbnd += 1;
                assert!(bndptr[i] != -1);
                assert!(bndind[bndptr[i] as usize] == i as idx_t);
                break;
            }
        }
    }

    assert!(
        nbnd == graph.nbnd,
        "calculated boundry size vs actual: {} vs {}",
        nbnd,
        graph.nbnd
    );

    1
}

/// This function checks whether or not the boundary information is correct
#[metis_func]
pub extern "C" fn CheckBnd2(graph: *const graph_t) -> idx_t {
    let graph = graph.as_ref().unwrap();
    let nvtxs = graph.nvtxs;
    get_graph_slices!(graph => xadj adjncy where_ bndptr bndind adjwgt);
    let mut id;
    let mut ed;

    let mut nbnd = 0;
    for i in 0..(nvtxs as usize) {
        id = 0;
        ed = 0;
        for j in xadj[i]..xadj[i + 1] {
            let j = j as usize;
            if where_[i] != where_[adjncy[j] as usize] {
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

    1
}
